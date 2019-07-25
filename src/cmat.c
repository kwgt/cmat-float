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
#define SWAP(a,b)           do {double t; t = a; a = b; b = t;} while(0)

static void
free_table(double** tbl, int rows)
{
  int i;

  for (i = 0; i < rows; i++) {
    if (tbl[i]) free(tbl[i]);
  }

  free(tbl);
}

int
alloc_table(int rows, int cols, double*** dst)
{
  int ret;
  double **tbl;
  int err;
  int i;

  ret = 0;
  tbl = NULL;

  do {
    tbl = (double**)malloc(sizeof(double*) * rows);
    if (!tbl) {
      ret = CMAT_ERR_NOMEM;
      break;
    }

    memset(tbl, 0, sizeof(double*) * rows);

    for (i = 0; i < rows; i++) {
      tbl[i] = (double*)malloc(sizeof(double) * cols);
      if (!tbl[i]) {
        ret = CMAT_ERR_NOMEM;
        goto loop_out;
      }
    }

    *dst = tbl;
  } while (0);

  loop_out:

  if (ret) {
    if (tbl) free_table(tbl, rows);
  }

  return ret;
}

static void
free_object(cmat_t* obj)
{
  if (obj->tbl) free_table(obj->tbl, obj->rows);
  free(obj);
}

static int
alloc_object(int rows, int cols, cmat_t* org, cmat_t** dst)
{
  int ret;
  double** tbl;
  cmat_t* obj;
  int i;

  /*
   * initialize
   */
  ret = 0;
  obj = NULL;
  tbl = NULL;

  do {
    obj = (cmat_t*)malloc(sizeof(cmat_t));
    if (obj == NULL) {
      ret = CMAT_ERR_NOMEM;
      break;
    }

    if (rows > 0) {
      ret = alloc_table(rows, cols, &tbl);
      if (ret) break;

    } else {
      tbl = NULL;
    }

    obj->tbl  = tbl;
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
    if (tbl) free_table(tbl, rows);
    if (obj) free(obj);
  }

  return ret;
}

static void
replace_object(cmat_t* ptr, cmat_t** src)
{
  free_table(ptr->tbl, ptr->rows);
  memcpy(ptr, *src, sizeof(cmat_t));
  free(*src);

  *src = NULL;
}


static double
cdot(double* op1[], int r, double* op2[], int c, int n)
{
  double ret;
  int i;

  ret = 0.0;

  for (i = 0; i < n; i++) {
    ret += op1[r][i] * op2[i][c];
  }

  return ret;
}

static double
det(double a, double b, double c, double d)
{
  return (a * d) - (b * c);
}

static double
calc_det_dim2(double* t[])
{
  return det(t[0][0], t[0][1], t[1][0], t[1][1]);
}

static double
calc_det_dim3(double* t[])
{
  double ret;

  ret = (t[0][0] * det(t[1][1], t[1][2], t[2][1], t[2][2])) -
        (t[1][0] * det(t[0][1], t[0][2], t[2][1], t[2][2])) +
        (t[2][0] * det(t[0][1], t[0][2], t[1][1], t[1][2]));

  return ret;
}

static double*
alloc_work(int n, double* src[])
{
  double* ret;
  int i;
  int j;

  ret = (double*)malloc(sizeof(double) * n * n);
  if (ret) {
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        *ret++ = src[i][j];
      }
    }

    ret -= (n * n);
  }

  return ret;
}

/*
 * http://hooktail.org/computer/index.php?LU%CA%AC%B2%F2
 */
static int
lu_decomp(double* p, int sz)
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

  for (i = 0, pi = p; i < sz; i++, pi += sz) {
    max = fabs(pi[i]);
    k   = i;

    /* 注目行以降で最大の値（絶対値）の存在する行を探す */
    for (j = i + 1, pj = p + (j * sz); j < sz; j++, pj += sz) {
      tmp = fabs(pj[i]);

      if (tmp > max) {
        max = tmp;
        k   = j;
      }
    }

    /* 注目行と最大値のあった行を入れ替える */
    if (i != k) {
      pj = p + (k * sz);

      for (j = 0; j < sz; j++) {
        SWAP(pj[j], pi[j]);
      }

      ret++;
    }

    /* この時点で対角成分が0の場合は注目行に対する分解は終わってると
       考えてよいので次の行に移動する */
    if (pi[i] == 0.0) continue;

    /* forwarding erase */
    for (j = i + 1, pj = p + (j * sz); j < sz; j++, pj += sz) {
      tmp = (pj[i] /= pi[i]);

      for (k = i + 1; k < sz; k++) {
        pj[k] -= tmp * pi[k];
      }
    }
  }

  return ret;
}

static int
calc_det(double* src[], int sz, double* dst)
{
  int ret;
  double* wrk;
  double det;
  int i;
  int j;
  int tmp;

  do {
    ret = 0;

    /* alloc work buffer */
    wrk = alloc_work(sz, src);
    if (wrk == NULL) {
      ret = CMAT_ERR_NOMEM;
      break;
    }

    /* do LU decomposition */
    tmp = lu_decomp(wrk, sz);

    /* calc diagonal multiplier */
    det = (tmp & 1)? -1.0: 1.0;
    tmp = sz + 1; 

    for (i = 0, j = 0; i < sz; i++, j += tmp) {
      det *= wrk[j];
    }
  } while (0);

  if (wrk) free(wrk);
  if (!ret) *dst = det;

  return ret;
}

/*
 * http://thira.plavox.info/blog/2008/06/_c.html
 */
static int
calc_inverse(double* src[], int n, double* dst[])
{
  int ret;
  double* wrk;
  double tmp;
  int i;
  int j;
  int k;

  double* wi;
  double* wj;
  double* di;
  double* dj;

  do {
    ret = 0;

    /* alloc work buffer */
    wrk = alloc_work(n, src);
    if (wrk == NULL) {
      ret = CMAT_ERR_NOMEM;
      break;
    }

    /* create identity matrix */
    for (i = 0; i < n; i++) {
      di = dst[i];

      for (j = 0; j < n; j++) {
        di[j] = (i == j)? 1.0: 0.0;
      }
    }

    /* do row reduction method */
    for (i = 0, wi = wrk; i < n; i++, wi += n) {
      di  = dst[i];
      tmp = 1.0 / wi[i];

      for (j = 0; j < n; j++) {
        wi[j] *= tmp;
        di[j] *= tmp;
      }

      for (j = 0, wj = wrk; j < n; j++, wj += n) {
        if (i == j) continue;

        dj  = dst[j];
        tmp = wj[i];

        for (k = 0; k < n; k++) {
          wj[k] -= wi[k] * tmp;
          dj[k] -= di[k] * tmp;
        }
      }
    }
  } while (0);

  if (wrk) free(wrk);

  return ret;
}

static void
calc_inverse2(double* src[], int n, double* dst[])
{
  double tmp;
  int i;
  int j;
  int k;

  double* si;
  double* sj;
  double* di;
  double* dj;

  /* create identity matrix */
  for (i = 0; i < n; i++) {
    di = dst[i];

    for (j = 0; j < n; j++) {
      di[j] = (i == j)? 1.0: 0.0;
    }
  }

  /* do row reduction method */
  for (i = 0; i < n; i++) {
    si  = src[i];
    di  = dst[i];
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

/**
 * 行列オブジェクトの生成
 *
 * @param n   列数の指定
 * @param dst   生成したオブジェクトの格納先のポインタ
 *
 * @return エラーコード(0で正常終了)
 */
int
cmat_new(int n, cmat_t** dst)
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
    if (n <= 0) {
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
    ret = alloc_object(0, n, NULL, &obj);
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
 * 行列オブジェクトの生成
 *
 * @param val   行列の初期値を一次元展開した配列
 * @param dst   生成したオブジェクトの格納先のポインタ
 *
 * @return エラーコード(0で正常終了)
 */
int
cmat_new2(double* val, int rows, int cols, cmat_t** dst)
{
  int ret;
  cmat_t* obj;
  int r;
  int c;

  /*
   * initialize
   */
  ret = 0;
  obj = NULL;

  /*
   * argument check
   */
  do {
    if (val == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }

    if (rows < 0) {
      ret = CMAT_ERR_BSIZE;
      break;
    }

    if (cols < 0) {
      ret = CMAT_ERR_BSIZE;
      break;
    }

    if (dst == NULL) {
      ret = CMAT_ERR_BADDR;
      break;
    }
  } while (0);

  /*
   * alloc new object
   */
  if (!ret) {
    ret = alloc_object(rows, cols, NULL, &obj);
  }

  /*
   * put return parameter
   */
  if (!ret) {
    for (r = 0; r < rows; r++) {
      for (c = 0; c < cols; c++) {
        obj->tbl[r][c] = *val++;
      }
    }

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
      for (c = 0; c < ptr->cols; c++) {
        len = format(ptr->tbl[r][c], str, ptr->coff);
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
      if (label != NULL) printf("  ");
      printf("[");

      for (c = 0; c < ptr->cols; c++) {
        format(ptr->tbl[r][c], str, ptr->coff);
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
  double** tbl;
  double* row;
  int capa;

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
   * alloc row
   */
  if (!ret) {
    row = (double*)malloc(sizeof(double) * ptr->cols);
    if (row == NULL) ret = CMAT_ERR_NOMEM;
  }

  /*
   * grow table
   */
  if (!ret) {
    if (ptr->capa == ptr->rows) {
      capa = (ptr->capa < 10)? 10: GROW(ptr->capa);
      tbl  = (double**)realloc(ptr->tbl, sizeof(double*) * capa);
      if (tbl) {
        ptr->tbl  = tbl;
        ptr->capa = capa;

      } else {
        ret = CMAT_ERR_NOMEM;
      }
    }
  }

  /*
   * update context
   */
  if (!ret) {
    memcpy(row, src, sizeof(double) * ptr->cols);
    ptr->tbl[ptr->rows++] = row;
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
      d = obj->tbl[r];
      s = ptr->tbl[r];
      o = op->tbl[r];

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
      d = obj->tbl[r];
      s = ptr->tbl[r];
      o = op->tbl[r];

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
      d = obj->tbl[r];
      s = ptr->tbl[r];

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
      for (c = 0; c < op->cols; c++) {
        obj->tbl[r][c] = cdot(ptr->tbl, r, op->tbl, c, ptr->cols);
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
      for (c = 0; c < ptr->cols; c++) {
        obj->tbl[c][r] = ptr->tbl[r][c];
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
 */
int
cmat_inverse(cmat_t* ptr, cmat_t** dst)
{
  int ret;
  cmat_t* obj;
  double** tbl;
  double det;

  /*
   * initialize
   */
  ret = 0;
  obj = NULL;
  tbl = NULL;

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
    if (fabs(det) < ptr->coff) ret = DEFAULT_ERROR;
  }

  /*
   * alloc result object
   */
  if (!ret) {
    if (dst) {
      ret = alloc_object(ptr->rows, ptr->cols, ptr, &obj);
    } else {
      ret = alloc_table(ptr->rows, ptr->cols, &tbl);
    }
  }

  /*
   * calculate inverse matrix
   */
  if (!ret) {
    if (dst) {
      ret = calc_inverse(ptr->tbl, ptr->rows, obj->tbl);
    } else {
      calc_inverse2(ptr->tbl, ptr->rows, tbl);
    }
  }

  /*
   * put return parameter
   */
  if (!ret) {
    if (dst) {
      *dst = obj;
    } else {
      free_table(ptr->tbl, ptr->rows);
      ptr->tbl = tbl;
    }
  }

  /*
   * post process
   */
  if (ret) {
    if (obj) free_object(obj);
    if (tbl) free_table(tbl, ptr->rows);
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
  double** tbl;
  double det;
  double tmp;
  int i;
  int j;
  int k;
  int n;

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
    tbl = ptr->tbl;

    switch (ptr->rows) {
    case 1:
      det = tbl[0][0];
      break;

    case 2:
      det = calc_det_dim2(ptr->tbl);
      break;

    case 3:
      det = calc_det_dim3(ptr->tbl);
      break;

    default:
      ret = calc_det(ptr->tbl, ptr->rows, &det);
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
    n = ptr->rows * ptr->cols;

    for (i = 0; i < n; i++) {
      if (i % ptr->cols == 0) s = ptr->tbl[i / ptr->cols];
      if (i % op->cols == 0) o = op->tbl[i / op->cols];

      dot += *s++ * *o++;
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

  double* p1;
  double* p2;

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
      p1 = ptr->tbl[r];
      p2 = op->tbl[r];

      for (c = 0; c < ptr->cols; c++) {
        if (fabs(p1[c] - p2[c]) > ptr->coff) goto loop_out;
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
 * @param dst   チェック結果先のポインタ(非零で一致)
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

  double* p1;

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
      p1 = ptr->tbl[r];

      for (c = 0; c < ptr->cols; c++) {
        if (fabs(p1[c] - *val++) > ptr->coff)  goto loop_out;
      }
    }

    /* mark mached */
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

