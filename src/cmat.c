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

static int
alloc_object(int rows, int cols, cmat_t* org, cmat_t** dst)
{
  int ret;
  double* tbl;
  cmat_t* obj;
  int i;

  /*
   * initialize
   */
  ret = 0;
  tbl = NULL;
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

    } else {
      tbl = NULL;
    }

    /*
     * setup object
     */
    obj->tbl = tbl;
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
    if (tbl) free(tbl);
    if (obj) free(obj);
  }

  return ret;
}

static void
free_object(cmat_t* ptr)
{
  if (ptr->tbl) free(ptr->tbl);
  free(ptr);
}

static void
replace_object(cmat_t* ptr, cmat_t** src)
{
  free(ptr->tbl);
  memcpy(ptr, *src, sizeof(cmat_t));
  free(*src);

  *src = NULL;
}

static double*
alloc_work(int n, double* src)
{
  double* ret;
  int i;
  int j;

  n *= n;

  ret = (double*)malloc(sizeof(double) * n);
  if (ret) memcpy(ret, src, sizeof(double) * n);

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

static double
det(double a, double b, double c, double d)
{
  return (a * d) - (b * c);
}

static double
calc_det_dim2(double* t)
{
  return det(t[0], t[1], t[2], t[3]);
}

static double
calc_det_dim3(double* t)
{
  double ret;

  ret = (t[0] * det(t[4], t[5], t[7], t[8])) -
        (t[3] * det(t[1], t[2], t[7], t[8])) +
        (t[6] * det(t[1], t[2], t[4], t[5]));

  return ret;
}

static int
calc_det(double* src, int sz, double* dst)
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

static void
calc_inverse(double* src, int n, double* dst)
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
  memset(dst, 0, sizeof(double) * n * n);
  for (i = 0, di = dst; i < n; i++, di += n) {
    di[i] = 1.0;
  }

  /* do row reduction method */
  for (i = 0, si = src, di = dst; i < n; i++, si += n, di += n) {
    max = fabs(si[i]);
    k   = i;

    /* ピボット操作 */
    // 注目行以降で最大の値（絶対値）の存在する行を探す
    for (j = i + 1, sj = src + (j * n); j < n; j++, sj += n) {
      tmp = fabs(sj[i]);

      if (tmp > max) {
        max = tmp;
        k   = j;
      }
    }

    // 注目行と最大値のあった行を入れ替える
    if (i != k) {
      sj = src + (k * n);
      dj = dst + (k * n);

      for (j = 0; j < n; j++) {
        SWAP(sj[j], si[j]);
        SWAP(dj[j], di[j]);
      }
    }

    if (si[i] == 0.0) continue;

    /* ここからガウス・ジョルダン法 */
    tmp = 1.0 / si[i];

    for (j = 0; j < n; j++) {
      si[j] *= tmp;
      di[j] *= tmp;
    }

    for (j = 0, sj = src, dj = dst; j < n; j++, sj += n, dj += n) {
      if (i == j) continue;

      tmp = sj[i];

      for (k = 0; k < n; k++) {
        sj[k] -= si[k] * tmp;
        dj[k] -= di[k] * tmp;
      }
    }
  }
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

    for (r = 0, rp = ptr->tbl; r < ptr->rows; r++, rp += ptr->cols) {
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

    for (r = 0, rp = ptr->tbl; r < ptr->rows; r++, rp += ptr->cols) {
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
   * grow table 
   */
  if (!ret) {
    if (ptr->capa == ptr->rows) {
      capa = (ptr->capa < 10)? 10: GROW(ptr->capa);
      tbl  = (double*)realloc(ptr->tbl, sizeof(double*) * ptr->cols * capa);
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
    memcpy(ptr->tbl + (ptr->rows++ * ptr->cols),
           src,
           sizeof(double) * ptr->cols);
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
  int i;
  int n;

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
    n = ptr->rows * ptr->cols;
    d = obj->tbl;
    s = ptr->tbl;
    o = op->tbl;

    for (i = 0; i < n; i++) {
      d[i] = s[i] + o[i];
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
  int i;
  int n;

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
    n = ptr->rows * ptr->cols;
    d = obj->tbl;
    s = ptr->tbl;
    o = op->tbl;

    for (i = 0; i < n; i++) {
      d[i] = s[i] - o[i];
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
  int i;
  int n;

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
    n = ptr->rows * ptr->cols;
    d = obj->tbl;
    s = ptr->tbl;

    for (i = 0; i < n; i++) {
      d[i] = s[i] * op;
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
  double* o;

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
    d = obj->tbl;
    s = ptr->tbl;

    for (r = 0; r < ptr->rows; r++) {
      for (c = 0; c < op->cols; c++) {

        d[c] = 0.0;
        o = op->tbl + c;

        for (i = 0; i < ptr->cols; i++) {
          d[c] += s[i] * o[0];
          o += op->cols;
        }
      }

      d += obj->cols;
      s += ptr->cols;
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
    for (r = 0, s = ptr->tbl; r < ptr->rows; r++, s += ptr->cols) {
      for (c = 0, d = obj->tbl + r; c < ptr->cols; c++, d += obj->cols) {
        *d = s[c];
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
  double* sp;
  double* dp;
  double det;

  /*
   * initialize
   */
  ret = 0;
  obj = NULL;
  sp  = NULL;
  dp  = NULL;

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
    sp = (double*)malloc(sizeof(double) * ptr->capa * ptr->cols);
    if (sp == NULL) ret = CMAT_ERR_NOMEM;
  }

  if (!ret) {
    if (dst) {
      ret = alloc_object(ptr->rows, ptr->cols, ptr, &obj);
      if (!ret) {
        memcpy(sp, ptr->tbl, sizeof(double) * ptr->rows * ptr->cols);
        dp = obj->tbl;
      }

    } else {
      dp = sp;
      sp = ptr->tbl;
    }
  }

  /*
   * calculate inverse matrix
   */
  if (!ret) {
    calc_inverse(sp, ptr->cols, dp);
  }

  /*
   * put return parameter
   */
  if (!ret) {
    if (dst) {
      *dst = obj;
    } else {
      free(ptr->tbl);
      ptr->tbl = dp;
    }
  }

  /*
   * post process
   */
  if (ret) {
    if (dst) {
      if (obj) free_object(obj);
    } else {
      if (dp) free(dp);
    }
  }

  if (dst) {
    if (sp) free(sp);
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
      det = calc_det_dim2(ptr->tbl);
      break;

    case 3:             // when 3x3
      det = calc_det_dim3(ptr->tbl);
      break;

    default:            // when nxn
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
    s = ptr->tbl;
    o = op->tbl;

    for (i = 0; i < n; i++) {
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
  int i;
  int n;
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
    n = ptr->rows * ptr->cols;
    p = ptr->tbl;
    o = op->tbl;

    for (i = 0; i < n; i++) {
      if (fcmp(p[i], o[i], ptr->coff)) goto loop_out;
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
  int i;
  int n;

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
    n = ptr->rows * ptr->cols;
    p = ptr->tbl;

    for (i = 0; i < n; i++) {
      if (fcmp(p[i], val[i], ptr->coff)) goto loop_out;
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

