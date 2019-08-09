/*
 * Cheap maxtrix library for C
 *
 *  Copyright (C) 2019 Hiroshi Kuwagata <kgt9221@gmail.com>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "cmat.h"

#define DEFAULT_ERROR       __LINE__
#define DEFAULT_CUTOFF      1e-10

#define GROW(n)             ((n * 13) / 10)
#define SHRINK(n)           ((n * 10) / 13)
#define SWAP(a,b,t)         do {t c; c = (a); (a) = (b); (b) = c;} while(0)

static int
alloc_object(int rows, int cols, cmat_t* org, cmat_t** dst)
{
  int ret;
  double* tbl;
  double** row;
  cmat_t* obj;
  int i;

  /*
   * initialize
   */
  ret = 0;
  tbl = NULL;
  row = NULL;
  obj = NULL;

  do {
    /*
     * alloc memory
     */
    obj = (cmat_t*)malloc(sizeof(cmat_t));
    if (obj == NULL) {
      ret = CMAT_ERR_NOMEM;
      break;
    }

    if (rows > 0) {
      tbl = (double*)malloc(sizeof(double) * rows * cols);
      if (tbl == NULL) {
        ret = CMAT_ERR_NOMEM;
        break;
      }

      row = (double**)malloc(sizeof(double*) * rows);
      if (row == NULL) {
        ret = CMAT_ERR_NOMEM;
        break;
      }
    }

    /*
     * setup object
     */

    for (i = 0; i < rows; i++) {
      row[i] = tbl + (i * cols);
    }

    obj->tbl  = tbl;
    obj->row  = row;
    obj->rows = rows;
    obj->cols = cols;
    obj->capa = rows;

    if (org) {
      obj->coff = org->coff;
    } else {
      obj->coff = DEFAULT_CUTOFF;
    }

    *dst = obj;
  } while (0);

  /*
   * post process
   */
  if (ret) {
    if (obj) free(obj);
    if (tbl) free(tbl);
    if (row) free(row);
  }

  return ret;
}

static void
free_object(cmat_t* ptr)
{
  if (ptr->tbl) free(ptr->tbl);
  if (ptr->row) free(ptr->row);
  free(ptr);
}

static void
replace_object(cmat_t* ptr, cmat_t** src)
{
  free(ptr->tbl);
  free(ptr->row);
  memcpy(ptr, *src, sizeof(cmat_t));
  free(*src);

  *src = NULL;
}

static int
alloc_table(double** src, int rows, int cols, double** dt, double*** dr)
{
  int ret;
  double* tbl;
  double** row;
  int i;

  ret = 0;
  tbl = NULL;
  row = NULL;

  do {
    tbl = (double*)malloc(sizeof(double) * rows * cols);
    if (tbl == NULL) {
      ret = CMAT_ERR_NOMEM;
      break;
    }

    row = (double**)malloc(sizeof(double) * rows);
    if (row == NULL) {
      ret = CMAT_ERR_NOMEM;
      break;
    }
  } while (0);

  if (!ret) {
    for (i = 0; i < rows; i++) {
      row[i] = tbl + (i * cols);
      if (src) memcpy(row[i], src[i], sizeof(double) * cols);
    }

    *dt = tbl;
    *dr = row;
  }

  if(ret) {
    if (tbl) free(tbl);
    if (row) free(row);
  }

  return ret;
}

static int
format(double val, char* dst, double thr)
{
  int ret;
  int i;

  if (fabs(val) > thr) {
    sprintf(dst, "% f", val);

    for (i = strlen(dst) - 1; i > 0; i--) {
      switch (dst[i]) {
      case '0':
        break;

      case '.':
        dst[i + 0] = '\0';
        ret = i;
        goto loop_out;

      default:
        dst[i + 1] = '\0';
        ret = i + 1;
        goto loop_out;
      }
    }

  } else {
    strcpy(dst, " 0");
    ret = 1;
  }
  loop_out:

  return ret;
}

static inline int
fcmp(double a, double b, double coff)
{
  int ret;

  double v1;
  int e1;

  double v2;
  int e2;

  do {
    ret = !0;

    v1  = frexp(a, &e1);
    v2  = frexp(b, &e2);

    if (e1 != e2) {
      v1 = a;
      v2 = b;
    }

    if (fabs(v1 - v2) > coff) break;

    ret = 0;
  } while (0);

  return ret;
}

/*
 * http://hooktail.org/computer/index.php?LU%CA%AC%B2%F2
 */
static int
lu_decomp(double** row, int sz, double thr, int* piv)
{
  int ret;
  int i;
  int j;
  int k;
  double max;
  double tmp;

  double* pi;
  double* pj;

  ret = 0;

  if (piv) {
    for (i = 0; i < sz; i++) piv[i] = i;
  }

  for (i = 0; i < sz; i++) {
    pi  = row[i];
    max = fabs(pi[i]);
    k   = i;

    /* 注目行以降で最大の値（絶対値）の存在する行を探す */
    for (j = i + 1; j < sz; j++) {
      tmp = fabs(row[j][i]);

      /*
       * 浮動小数点数の丸め誤差の蓄積のため、極小差の場合に大小比較がうまくい
       * かないことがある。これを避けるために差の閾値比較で大小比較を行ってい
       * る。以下のif文は tmp > max を評価している。
       */
      if (tmp - max > thr) {
        max = tmp;
        k   = j;
      }
    }

    /* 注目行と最大値のあった行を入れ替える */
    if (k != i) {
      SWAP(row[i], row[k], double*);
      if (piv) SWAP(piv[i], piv[k], int);

      pi = row[i];
      ret++;
    }

    /* この時点で対角成分が0の場合は注目行に対する分解は終わってると
       考えてよいので次の行に移動する */
    if (pi[i] == 0.0) continue;

    /* forwarding erase */
    for (j = i + 1; j < sz; j++) {
      pj  = row[j];
      tmp = (pj[i] /= pi[i]);

      for (k = i + 1; k < sz; k++) {
        pj[k] -= tmp * pi[k];
      }
    }
  }

  return ret;
}

static double
det(double a, double b, double c, double d)
{
  return (a * d) - (b * c);
}

static double
calc_det_dim2(double* r1, double* r2)
{
  return det(r1[0], r1[1], r2[0], r2[1]);
}

static double
calc_det_dim3(double* r1, double* r2, double* r3)
{
  double ret;

  ret = (r1[0] * det(r2[1], r2[2], r3[1], r3[2])) -
        (r2[0] * det(r1[1], r1[2], r3[1], r3[2])) +
        (r3[0] * det(r1[1], r1[2], r2[1], r2[2]));

  return ret;
}

static int
calc_det(double** row, int sz, double thr, double* dst)
{
  int ret;
  double* wt;   // as "Work Table"
  double** wr;   // as "Work Rows"
  double det;
  int i;
  int j;
  int n;

  do {
    ret = 0;

    /* alloc work buffer */
    ret = alloc_table(row, sz, sz, &wt, &wr);
    if (ret) break;

    /* do LU decomposition */
    n   = lu_decomp(wr, sz, thr, NULL);

    /* calc diagonal multiplier */
    det = (n & 1)? -1.0: 1.0;

    for (i = 0; i < sz; i++) {
      det *= wr[i][i];
    }
  } while (0);

  if (wr) free(wr);
  if (wt) free(wt);
  if (!ret) *dst = det;

  return ret;
}

static void
calc_inverse(double** src, int n, double** dst)
{
  int i;
  int j;
  int k;

  double max;
  double tmp;

  double* si;
  double* sj;
  double* di;
  double* dj;

  /* create identity matrix */
  for (i = 0; i < n; i++) {
    memset(dst[i], 0, sizeof(double) * n);
    dst[i][i] = 1.0;
  }

  /* do row reduction method */
  for (i = 0; i < n; i++) {
    si  = src[i];
    di  = dst[i];

    max = fabs(si[i]);
    k   = i;

    /* ピボット操作 */
    // 注目行以降で最大の値（絶対値）の存在する行を探す
    for (j = i + 1; j < n; j++) {
      tmp = fabs(src[j][i]);

      if (tmp > max) {
        max = tmp;
        k   = j;
      }
    }

    // 注目行と最大値のあった行を入れ替える
    if (i != k) {
      SWAP(src[i], src[k], double*);
      SWAP(dst[i], dst[k], double*);

      si = src[i];
      di = dst[i];
    }

    if (si[i] == 0.0) continue;

    /* ここからガウス・ジョルダン法 */
    tmp = 1.0 / si[i];

    for (j = 0; j < n; j++) {
      si[j] *= tmp;
      di[j] *= tmp;
    }

    for (j = 0; j < n; j++) {
      if (i == j) continue;

      sj  = src[j];
      dj  = dst[j];
      tmp = sj[i];

      for (k = 0; k < n; k++) {
        sj[k] -= si[k] * tmp;
        dj[k] -= di[k] * tmp;
      }
    }
  }
}

static void
sort(int* a, size_t n)
{
  int h;
  int f;
  int i;

  /*
   * sort by ascending order
   */

  h = n;

  do {
    if (h > 1) {
      h = SHRINK(h);
    } else if (!f) {
      break;
    }

    f = 0;

    if (h == 9 || h == 10) h = 11;

    for (i = 0; i < ((int)n - h); i++) {
      if (a[i] > a[i + h]) {
        SWAP(a[i], a[i + h], int);
        f = !0;
      }
    }
  } while (1);
}

/**
 * 行列オブジェクトの生成
 *
 * @param n   列数の指定
 * @param dst   生成したオブジェクトの格納先のポインタ
 *
 * @return エラーコード(0で正常終了)
 */
int
cmat_new(double* src, int rows, int cols, cmat_t** dst)
{
  int ret;
  cmat_t* obj;

  /*
   * initialize
   */
  ret = 0;
  obj = NULL;

  /*
   * argument check
   */
  do {
    if (rows < 0) {
      ret = CMAT_ERR_BSIZE;
      break;
    }

    if (cols <= 0) {
      ret = CMAT_ERR_BSIZE;
      break;
    }

    if (dst == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }
  } while(0);

  /*
   * alloc memory
   */
  if (!ret) {
    ret = alloc_object(rows, cols, NULL, &obj);
  }

  /*
   * set initial values
   */
  if (!ret) {
    if (src) {
      memcpy(obj->tbl, src, sizeof(double) * rows * cols);
    } else {
      memset(obj->tbl, 0, sizeof(double) * rows * cols);
    }
  }

  /*
   * put return parameter
   */
  if (!ret) {
    *dst = obj;
  }

  /*
   * post process
   */
  if (ret) {
    if (obj) free_object(obj);
  }

  return ret;
}

/**
 * 行列オブジェクトの複製
 *
 * @param ptr   複製元になる行列オブジェクト
 * @param dst   生成したオブジェクトの格納先のポインタ
 *
 * @return エラーコード(0で正常終了)
 */
int
cmat_clone(cmat_t* ptr, cmat_t** dst)
{
  int ret;
  cmat_t* obj;
  int i;

  /*
   * initialize
   */
  ret = 0;
  obj = NULL;

  /*
   * argument check
   */
  do {
    if (ptr == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }

    if (dst == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }
  } while(0);

  /*
   * alloc memory
   */
  if (!ret) {
    ret = alloc_object(ptr->rows, ptr->cols, ptr, &obj);
  }

  /*
   * copy values
   */
  if (!ret) {
    memcpy(obj->tbl, ptr->tbl, sizeof(double) * ptr->rows * ptr->cols);

    for (i = 0; i < ptr->rows; i ++) {
      obj->row[i] = obj->tbl + (ptr->row[i] - ptr->tbl);
    }
  }

  /*
   * put return parameter
   */
  if (!ret) {
    *dst = obj;
  }

  return ret;
}

/**
 * 行列オブジェクトの削除
 *
 * @param dst   削除するオブジェクトのポインタ
 *
 * @return エラーコード(0で正常終了)
 */
int
cmat_destroy(cmat_t* ptr)
{
  int ret;

  /*
   * initialize
   */
  ret = 0;

  /*
   * argument check
   */
  if (ptr == NULL) ret = CMAT_ERR_BADDR;

  /*
   * release memory
   */
  if (!ret) {
    free_object(ptr);
  }

  return ret;
}

/**
 * 行列オブジェクトの内容表示
 *
 * @param ptr   対象の行列オブジェクト
 *
 * @return エラーコード(0で正常終了)
 */
int
cmat_print(cmat_t* ptr, char* label)
{
  int ret;
  int r;
  int c;

  double* rp;   // as "Row Pointer"
  char fmt[32];
  char str[32];
  int len;
  int max;
  int i;

  /*
   * initialize
   */
  ret = 0;

  /*
   * argument check
   */
  if (ptr == NULL) ret = CMAT_ERR_BADDR;

  /*
   * make format string
   */
  if (!ret) {
    max = 0;

    for (r = 0; r < ptr->rows; r++) {
      rp = ptr->row[r];

      for (c = 0; c < ptr->cols; c++) {
        len = format(rp[c], str, ptr->coff);
        if (len > max) max = len;
      }
    }

    sprintf(fmt, "%%%ds", max);
  }

  /*
   * show content
   */
  if (!ret) {
    if (label != NULL) printf("%s:\n", label);

    for (r = 0; r < ptr->rows; r++) {
      rp = ptr->row[r];

      if (label != NULL) printf("  ");
      printf("[");

      for (c = 0; c < ptr->cols; c++) {
        format(rp[c], str, ptr->coff);
        printf(fmt, str);
        if (c < (ptr->cols - 1)) printf(" ");
      }

      printf(" ]\n");
    }
  }

  return ret;
}

/**
 * 行の追加
 *
 * @param ptr   追加対象の行列オブジェクト
 * @param src   追加する行のデータ
 *
 * @return エラーコード(0で正常終了)
 */
int
cmat_append(cmat_t* ptr, double* src)
{
  int ret;
  double* tbl;
  double** row;
  int capa;
  int i;

  /*
   * initialize
   */
  ret = 0;
  tbl = NULL;
  row = NULL;

  /*
   * argument check
   */
  do {
    if (ptr == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }

    if (src == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }
  } while (0);

  /*
   * grow table 
   */
  if (!ret) do {
    if (ptr->capa == ptr->rows) {
      capa = (ptr->capa < 10)? 10: GROW(ptr->capa);

      tbl  = (double*)realloc(ptr->tbl, sizeof(double) * capa * ptr->cols);
      if (tbl == NULL) {
        ret = CMAT_ERR_NOMEM;
        break;
      }

      row  = (double**)realloc(ptr->row, sizeof(double*) * capa);
      if (row == NULL) {
        ret = CMAT_ERR_NOMEM;
        break;
      }

      for (i = 0; i < capa; i++) {
        row[i] = tbl + (i * ptr->cols);
      }

      ptr->tbl  = tbl;
      ptr->row  = row;
      ptr->capa = capa;
    }
  } while (0);

  /*
   * update context
   */
  if (!ret) {
    memcpy(ptr->row[ptr->rows++], src, sizeof(double) * ptr->cols);
  }

  /*
   * post process
   */
  if (ret) {
    if (tbl) free(tbl);
    if (row) free(row);
  }

  return ret;
}

/**
 * 行列の和
 *  ptr + op → dst       (dst != NULL)
 *  ptr + op → ptr       (dst == NULL)
 *
 * @param ptr   対象の行列オブジェクト
 * @param op    和行列
 * @param dst   転置結果の格納先
 *
 * @return エラーコード(0で正常終了)
 */
int
cmat_add(cmat_t* ptr, cmat_t* op, cmat_t** dst)
{
  int ret;
  cmat_t* obj;
  int r;
  int c;

  double* d;
  double* s;
  double* o;

  /*
   * initialize
   */
  ret = 0;
  obj = NULL;

  /*
   * check argument
   */
  do {
    if (ptr == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }

    if (op == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }
  } while (0);

  /*
   * check shape
   */
  if (!ret) {
    if (ptr->rows != op->rows || ptr->cols != op->cols) ret = CMAT_ERR_SHAPE;
  }

  /*
   * select target
   */
  if (!ret) {
    if (dst) {
      ret = alloc_object(ptr->rows, ptr->cols, ptr, &obj);
    } else {
      obj = ptr;
    }
  }

  /*
   * do add operation
   */
  if (!ret) {
    for (r = 0; r < ptr->rows; r++) {
      d = obj->row[r];
      s = ptr->row[r];
      o = op->row[r];

      for (c = 0; c < ptr->cols; c++) {
        d[c] = s[c] + o[c];
      }
    }
  }

  /*
   * put return parameter
   */
  if (!ret) {
    if (dst) *dst = obj;
  }

  /*
   * post process
   */
  if (ret) {
    if (dst && obj) free_object(obj);
  }

  return ret;
}

/**
 * 行列の差
 *  ptr - op → dst       (dst != NULL)
 *  ptr - op → ptr       (dst == NULL)
 *
 * @param ptr   対象の行列オブジェクト
 * @param op    和行列
 * @param dst   転置結果の格納先
 *
 * @return エラーコード(0で正常終了)
 */
int
cmat_sub(cmat_t* ptr, cmat_t* op, cmat_t** dst)
{
  int ret;
  cmat_t* obj;
  int r;
  int c;

  double* d;
  double* s;
  double* o;

  /*
   * initialize
   */
  ret = 0;
  obj = NULL;

  /*
   * check argument
   */
  do {
    if (ptr == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }

    if (op == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }
  } while (0);

  /*
   * check shape
   */
  if (!ret) {
    if (ptr->rows != op->rows || ptr->cols != op->cols) ret = CMAT_ERR_SHAPE;
  }

  /*
   * select target
   */
  if (!ret) {
    if (dst) {
      ret = alloc_object(ptr->rows, ptr->cols, ptr, &obj);
    } else {
      obj = ptr;
    }
  }

  /*
   * do add operation
   */
  if (!ret) {
    for (r = 0; r < ptr->rows; r++) {
      d = obj->row[r];
      s = ptr->row[r];
      o = op->row[r];

      for (c = 0; c < ptr->cols; c++) {
        d[c] = s[c] - o[c];
      }
    }
  }

  /*
   * put return parameter
   */
  if (!ret) {
    if (dst) *dst = obj;
  }

  /*
   * post process
   */
  if (ret) {
    if (dst && obj) free_object(obj);
  }

  return ret;
}

/**
 * 行列のスカラー積
 *  ptr * op → dst       (dst != NULL)
 *  ptr * op → ptr       (dst == NULL)
 *
 * @param ptr   対象の行列オブジェクト
 * @param op    スカラー値
 * @param dst   転置結果の格納先
 *
 * @return エラーコード(0で正常終了)
 */
int
cmat_mul(cmat_t* ptr, double op, cmat_t** dst)
{
  int ret;
  cmat_t* obj;
  int r;
  int c;

  double* d;
  double* s;

  /*
   * initialize
   */
  ret = 0;
  obj = NULL;

  /*
   * check argument
   */
  do {
    if (ptr == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }

    if (isnan(op)) {
      ret = CMAT_ERR_INVAL;
      break;
    }
  } while (0);

  /*
   * select target
   */
  if (!ret) {
    if (dst) {
      ret = alloc_object(ptr->rows, ptr->cols, ptr, &obj);
    } else {
      obj = ptr;
    }
  }

  /*
   * do add operation
   */
  if (!ret) {
    for (r = 0; r < ptr->rows; r++) {
      d = obj->row[r];
      s = ptr->row[r];

      for (c = 0; c < ptr->cols; c++) {
        d[c] = s[c] * op;
      }
    }
  }

  /*
   * put return parameter
   */
  if (!ret) {
    if (dst) *dst = obj;
  }

  /*
   * post process
   */
  if (ret) {
    if (dst && obj) free_object(obj);
  }

  return ret;
}

/**
 * 行列の積
 *  ptr * op → dst       (dst != NULL)
 *  ptr * op → ptr       (dst == NULL)
 *
 * @param ptr   転置対象の行列オブジェクト
 * @param op    積行列
 * @param dst   演算結果の格納先
 *
 * @return エラーコード(0で正常終了)
 */
int
cmat_product(cmat_t* ptr, cmat_t* op, cmat_t** dst)
{
  int ret;
  cmat_t* obj;
  int r;
  int c;
  int i;

  double* d;
  double* s;

  /*
   * initialize
   */
  ret = 0;
  obj = NULL;

  /*
   * argument check
   */
  do {
    if (ptr == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }

    if (op == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }
  } while (0);

  /*
   * check op value
   */
  if (!ret) {
    if (ptr->cols != op->rows) ret = CMAT_ERR_SHAPE;
  }

  /*
   * alloc result object
   */
  if (!ret) {
    ret = alloc_object(ptr->rows, op->cols, ptr, &obj);
  }

  /*
   * do multiple operation
   */
  if (!ret) {

    for (r = 0; r < ptr->rows; r++) {
      d = obj->row[r];
      s = ptr->row[r];

      for (c = 0; c < op->cols; c++) {
        d[c] = 0.0;

        for (i = 0; i < ptr->cols; i++) {
          d[c] += s[i] * op->row[i][c];
        }
      }
    }
  }

  /*
   * put return parameter
   */
  if (!ret) {
    if (dst) {
      *dst = obj;
    } else {
      replace_object(ptr, &obj);
    }
  }

  /*
   * post process
   */
  if (ret) {
    if (obj) free_object(obj);
  }

  return ret;
}

/**
 * 行列の転置
 *  transpose(ptr) → dst       (dst != NULL)
 *  transpose(ptr) → ptr       (dst == NULL)
 *
 * @param ptr   転置対象の行列オブジェクト
 * @param dst   転置結果の格納先
 *
 * @return エラーコード(0で正常終了)
 */
int
cmat_transpose(cmat_t* ptr, cmat_t** dst)
{
  int ret;
  cmat_t* obj;
  int r;
  int c;

  double* s;

  /*
   * initialize
   */
  ret = 0;
  obj = NULL;

  /*
   * argument check
   */
  if (ptr == NULL) ret = CMAT_ERR_BADDR;

  /*
   * alloc result object
   */
  if (!ret) {
    ret = alloc_object(ptr->cols, ptr->rows, ptr, &obj);
  }

  /*
   * do transpose operation
   */
  if (!ret) {
    for (r = 0; r < ptr->rows; r++) {
      s = ptr->row[r];

      for (c = 0; c < ptr->cols; c++) {
        obj->row[c][r] = s[c];
      }
    }
  }

  /*
   * put return parameter
   */
  if (!ret) {
    if (dst) {
      *dst = obj;
    } else {
      replace_object(ptr, &obj);
    }
  }

  /*
   * post process
   */
  if (ret) {
    if (obj) cmat_destroy(obj);
  }

  return ret;
}

/**
 * 逆行列の算出
 *  inverse(ptr) → dst       (dst != NULL)
 *  inverse(ptr) → ptr       (dst == NULL)
 *
 * @param ptr   対象の行列オブジェクト
 * @param dst   逆行列の格納先
 *
 * @return エラーコード(0で正常終了)
 *
 * @refer http://thira.plavox.info/blog/2008/06/_c.html
 *        http://www.yamamo10.jp/yamamoto/lecture/2006/5E/
 *                      Linear_eauations/gaussj_html/node2.html
 */
int
cmat_inverse(cmat_t* ptr, cmat_t** dst)
{
  int ret;
  cmat_t* obj;

  double det;
  int i;

  double* st;   // as "Source Table"
  double** sr;  // as "Source Row"

  double* dt;   // as "Destination Table"
  double** dr;  // as "destination Row"

  /*
   * initialize
   */
  ret = 0;
  obj = NULL;
  st  = NULL;
  sr  = NULL;
  dt  = NULL;
  dr  = NULL;

  /*
   * argument check
   */
  if (ptr == NULL) ret = CMAT_ERR_BADDR;
  
  /*
   * check if it's a regular matrix
   */
  if (!ret) {
    ret = cmat_det(ptr, &det);
  }

  if (!ret) {
    if (fabs(det) < ptr->coff) ret = CMAT_ERR_NREGL;
  }

  /*
   * alloc work(or output) memory
   */
  if (!ret) {
    ret = alloc_table(NULL, ptr->capa, ptr->cols, &st, &sr);
  }

  /*
   * alloc result object
   */
  if (!ret) {
    if (dst) {
      ret = alloc_object(ptr->rows, ptr->cols, ptr, &obj);
      if (!ret) {
        for (i = 0; i < ptr->rows; i++) {
          memcpy(sr[i], ptr->row[i], sizeof(double) * ptr->cols);
        }

        dt = obj->tbl;
        dr = obj->row;
      }

    } else {
      dt = st;
      dr = sr;
      st = ptr->tbl;
      sr = ptr->row;
    }
  }

  /*
   * calculate inverse matrix
   */
  if (!ret) {
    calc_inverse(sr, ptr->rows, dr);
  }

  /*
   * put return parameter
   */
  if (!ret) {
    if (dst) {
      *dst = obj;
    } else {
      free(ptr->tbl);
      free(ptr->row);
      ptr->tbl = dt;
      ptr->row = dr;
    }
  }

  /*
   * post process
   */
  if (ret) {
    if (dst) {
      if (obj) free_object(obj);
    } else {
      if (dt) free(dt);
      if (dr) free(dr);
    }
  }

  if (dst) {
    if (st) free(st);
    if (sr) free(sr);
  }

  return ret;
}

/**
 * 行列のLU分解
 *  LU_decomp(ptr) → dst       (dst != NULL)
 *  LU_decomp(ptr) → ptr       (dst == NULL)
 *
 * @param ptr   対象の行列オブジェクト
 * @param dst   分解行列の格納先
 * @param piv   置換数列の格納先（必要ない場合はNULLを指定)
 *
 * @return エラーコード(0で正常終了)
 *
 * @note 上三角行列とした三角行列を合成した状態で出力するので注意
 *       出力行列は以下のようになる。
 *        UUUUUU
 *        LUUUUU
 *        LLUUUU 
 *        LLLUUU
 *        LLLLUU
 *        LLLLLU
 *      ※LU分解時の下三角行列の対角要素はすべて1なので省略している点に注意。
 *
 * @note 置換数列は置換先の数列で返される。置換行列への変換は呼び出し側で
 *       行う必要がある。
 *
 * @refer http://thira.plavox.info/blog/2008/06/_c.html
 */
int
cmat_lu_decomp(cmat_t* ptr, cmat_t** dst, int* piv)
{
  int ret;
  cmat_t* obj;

  int i;

  double** row;  // as "Source Row"

  /*
   * initialize
   */
  ret = 0;
  obj = NULL;
  row = NULL;

  /*
   * argument check
   */
  if (ptr == NULL) ret = CMAT_ERR_BADDR;
 
  /*
   * alloc result object
   */
  if (!ret) {
    if (dst) {
      ret = alloc_object(ptr->rows, ptr->cols, ptr, &obj);
      if (!ret) {
        for (i = 0; i < ptr->rows; i++) {
          memcpy(obj->row[i], ptr->row[i], sizeof(double) * ptr->cols);
        }

        row = obj->row;
      }

    } else {
      row = ptr->row;
    }
  }

  /*
   * do LU decompression
   */
  if (!ret) {
    lu_decomp(row, ptr->rows, ptr->coff, piv);
  }

  /*
   * put return parameter
   */
  if (!ret) {
    if (dst) *dst = obj;
  }

  /*
   * post process
   */
  if (ret) {
    if (dst) {
      if (obj) free_object(obj);
    }
  }

  return ret;
}

/**
 * 行列式の計算
 *  det(ptr) → dst
 *
 * @param ptr   対象の行列オブジェクト
 * @param dst   算出結果の格納先のポインタ
 *
 * @return エラーコード(0で正常終了)
 */
int
cmat_det(cmat_t* ptr, double* dst)
{
  int ret;
  double det;

  /*
   * initialize
   */
  ret = 0;
  det = -1.0;

  /*
   * argument check
   */
  do {
    if (ptr == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }

    if (dst == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }
  } while (0);

  /*
   * check shape
   */
  if (!ret) {
    if (ptr->rows != ptr->cols) ret = CMAT_ERR_SHAPE;
  }

  /*
   * calc determinant
   */
  if (!ret) {
    switch (ptr->rows) {
    case 1:             // when 1x1
      det = ptr->tbl[0];
      break;

    case 2:             // when 2x2
      det = calc_det_dim2(ptr->row[0], ptr->row[1]);
      break;

    case 3:             // when 3x3
      det = calc_det_dim3(ptr->row[0], ptr->row[1], ptr->row[2]);
      break;

    default:            // when nxn
      ret = calc_det(ptr->row, ptr->rows, ptr->coff, &det);
      break;
    }
  }

  /*
   * put return parameter
   */
  if (!ret) {
    *dst = det;
  }

  return ret;
}

/**
 * 行列のドット積の計算
 *  ptr * op → dst
 *
 * @param ptr   対象の行列オブジェクト
 * @param op    オペランド
 * @param dst   算出結果の格納先のポインタ
 *
 * @return エラーコード(0で正常終了)
 */
int
cmat_dot(cmat_t* ptr, cmat_t* op, double* dst)
{
  int ret;
  double dot;
  int i;
  int n;

  int r1;
  int r2;

  double* s;
  double* o;

  /*
   * initialize
   */
  ret = 0;
  dot = 0.0;

  /*
   * argument check
   */
  do {
    if (ptr == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }

    if (op == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }

    if (dst == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }
  } while (0);

  /*
   * check shape
   */
  if (!ret) {
    if ((ptr->rows * ptr->cols) != (op->rows * op->cols)) ret = CMAT_ERR_SHAPE;
  }

  /*
   * calc dot product
   */
  if (!ret) {
    n  = ptr->rows * ptr->cols;
    r1 = 0;
    r2 = 0;

    for (i = 0; i < n; i++) {
      if (i % ptr->cols == 0) s = ptr->row[r1++];
      if (i % op->cols == 0) o = op->row[r2++];

      dot += s[i % ptr->cols] * o[i % op->cols];
    }
  }

  /*
   * put return parameter
   */
  if (!ret) {
    *dst = dot;
  }

  return ret;
}

#ifdef DEBUG
/**
 * 行の置換
 *
 * @param ptr   置換対象の行列オブジェクト
 * @param _piv  置換数列
 *
 * @return エラーコード
 *
 * @note 本館数は呼び出しオブジェクトを書き換える。
 * @note 置換数列にはcmat_lu_decomp()が返す数列をそのまま使用できる。
 */
int
cmat_permute_row(cmat_t* ptr, int* _piv)
{
  int ret;
  int r;
  int* piv;
  int i;

  /*
   * initialize
   */
  ret = 0;
  piv = NULL;

  /*
   * argument check
   */
  do {
    if (ptr == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }

    if (_piv == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }
  } while (0);

  /*
   * alloc pivots array
   */
  if (!ret) {
    piv = (int*)malloc(sizeof(int) * ptr->rows);
    if (piv == NULL) ret = CMAT_ERR_NOMEM;
  }

  /*
   * check pivots
   */
  if (!ret) {
    memcpy(piv, _piv, sizeof(int) * ptr->rows);

    sort(piv, ptr->rows);
    for (i = 0; i < ptr->rows; i++) {
      if (piv[i] != i) {
        ret = CMAT_ERR_INVAL;
        break;
      }
    }
  }

  /*
   * do permutation row
   */
  if (!ret) {
    memcpy(piv, _piv, sizeof(int) * ptr->rows);

    for (r = 0; r < ptr->rows; r++) {
      for (i = r; i < ptr->rows; i++) {
        if (piv[i] == r) break;
      }

      if (r != i) {
        SWAP(ptr->row[r], ptr->row[i], double*);
        SWAP(piv[r], piv[i], int);
      }
    }
  }

  /*
   * post process
   */
  if (piv) free(piv);

  return ret;
}

/**
 * 列の置換
 *
 * @param ptr   置換対象の行列オブジェクト
 * @param piv   置換数列
 *
 * @return エラーコード
 *
 * @note 本館数は呼び出しオブジェクトを書き換える。
 * @note 置換数列にはcmat_lu_decomp()が返す数列をそのまま使用できる。
 */
int
cmat_permute_column(cmat_t* ptr, int* _piv)
{
  int ret;
  int r;
  int* piv;
  double* p;
  int i;
  int j;

  /*
   * initialize
   */
  ret = 0;
  piv = NULL;

  /*
   * argument check
   */
  do {
    if (ptr == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }

    if (_piv == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }
  } while (0);

  /*
   * alloc pivots array
   */
  if (!ret) {
    piv = (int*)malloc(sizeof(int) * ptr->cols);
    if (piv == NULL) ret = CMAT_ERR_NOMEM;
  }

  /*
   * check pivots
   */
  if (!ret) {
    memcpy(piv, _piv, sizeof(int) * ptr->cols);

    sort(piv, ptr->rows);
    for (i = 0; i < ptr->cols; i++) {
      if (piv[i] != i) {
        ret = CMAT_ERR_INVAL;
        break;
      }
    }
  }

  /*
   * do permutation column
   */
  if (!ret) {
    memcpy(piv, _piv, sizeof(int) * ptr->cols);

    for (r = 0; r < ptr->rows; r++) {
      for (i = r; i < ptr->cols; i++) {
        if (piv[i] == r) break;
      }

      if (r != i) {
        for (j = 0; j < ptr->rows; j++) {
          SWAP(ptr->row[j][r], ptr->row[j][i], double);
        }
        SWAP(piv[r], piv[i], int);
      }
    }
  }

  /*
   * post process
   */
  if (piv) free(piv);

  return ret;
}
#endif /* defined(DEBUG) */

/**
 * 行列の比較
 *
 * @param ptr   対象の行列オブジェクト
 * @param op    比較対象の行列オブジェクト
 * @param dst   チェック結果先のポインタ(0で一致)
 *
 * @return エラーコード(0で正常終了)
 */
int
cmat_compare(cmat_t* ptr, cmat_t* op, int* dst)
{
  int ret;
  int res;
  int r;
  int c;
  double* p;
  double* o;

  /*
   * initialize
   */
  ret = 0;
  res = !0;

  /*
   * argument check
   */
  do {
    if (ptr == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }

    if (op == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }

    if (dst == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }
  } while (0);

  /*
   * compare matrixies
   */
  if (!ret) do {
    /* check shape */
    if (ptr->rows != op->rows || ptr->cols != op->cols) break;

    /* check values */
    for (r = 0; r < ptr->rows; r++) {
      p = ptr->row[r];
      o = op->row[r];

      for (c = 0; c < ptr->cols; c++) {
        if (fcmp(p[c], o[c], ptr->coff)) goto loop_out;
      }
    }

    /* mark matched */
    res = 0;

  } while (0);
  loop_out:

  /*
   * put return parameter
   */
  if (!ret) {
    *dst = res;
  }

  return ret;
}

/**
 * 行列内容のチェック
 *
 * @param ptr   転置対象の行列オブジェクト
 * @param val   チェックする行列を一次元展開した配列
 * @param dst   チェック結果先のポインタ(0で一致)
 *
 * @return エラーコード(0で正常終了)
 */
int
cmat_check(cmat_t* ptr, double* val, int* dst)
{
  int ret;
  int res;
  int r;
  int c;

  double* p;

  /*
   * initialize
   */
  ret = 0;
  res = !0;

  /*
   * argument check
   */
  do {
    if (ptr == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }

    if (val == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }

    if (dst == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }
  } while (0);

  /*
   * put return parameter
   */
  if (!ret) {
    /* check values */
    for (r = 0; r < ptr->rows; r++) {
      p = ptr->row[r];

      for (c = 0; c < ptr->cols; c++) {
        if (fcmp(p[c], *val++, ptr->coff)) goto loop_out;
      }
    }

    /* mark matched */
    res = 0;
  }
  loop_out:

  /*
   * put return parameter
   */
  if (!ret) {
    *dst = res;
  }

  return ret;
}

/**
 * 切り捨て処理の閾値の設定
 *
 * @param ptr   対象の行列オブジェクト
 * @param val   閾値の値
 *
 * @return エラーコード(0で正常終了)
 */
int
cmat_set_cutoff_threshold(cmat_t* ptr, double val)
{
  int ret;

  /*
   * initialize
   */
  ret = 0;

  /*
   * check argument
   */
  if (ptr == NULL) ret = CMAT_ERR_BADDR;

  /*
   * update context
   */
  if (!ret) {
    ptr->coff = fabs(val);
  }

  return ret;
} 

