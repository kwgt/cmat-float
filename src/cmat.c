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
      ret = DEFAULT_ERROR;
      break;
    }

    if (rows > 0) {
      tbl = (double**)malloc(sizeof(double*) * rows);
      if (tbl == NULL) {
        ret = DEFAULT_ERROR;
        break;
      }

      memset(tbl, 0, sizeof(double*) * rows);

      for (i = 0; i < rows; i++) {
        tbl[i] = (double*)malloc(sizeof(double) * cols);
        if (tbl[i] == NULL) {
          ret = DEFAULT_ERROR;
          break;
        }
      }

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
    if (obj) free(obj);

    if (tbl) {
      for (i = 0; i < rows; i++) {
        if (tbl[i]) {
          free(tbl[i]);
        } else {
          break;
        }
      }

      free(tbl);
    }
  }

  return ret;
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
det_dim2(double* t[])
{
  return det(t[0][0], t[0][1], t[1][0], t[1][1]);
}

static double
det_dim3(double* t[])
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
 * http://thira.plavox.info/blog/2008/06/_c.html
 */
static int
det_matrix(double* src[], int n, double* dst)
{
  int ret;
  double* wrk;
  double det;
  double tmp;
  int i;
  int j;
  int k;

  double *wi;
  double *wj;

  do {
    ret = 0;

    /* alloc work buffer */
    wrk = alloc_work(n, src);
    if (wrk == NULL) {
      ret = DEFAULT_ERROR;
      break;
    }

    /* calc triagular matrix */
    for (i = 0, wi = wrk; i < n; i++, wi += n) {
      for (j = 0, wj = wrk; j < n; j++, wj += n) {
        if (i == j) continue;

        tmp = wj[i] / wi[i];
        for (k = 0; k < n; k++) {
          wj[k] -= wi[k] * tmp;
        }
      }
    }

    /* calc diagonal multiplier */
    det = 1.0;

    for (i = 0, wi = wrk; i < n; i++, wi += (n + 1)) {
      det *= *wi;
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
inverse_matrix(double* src[], int n, double* dst[])
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
      ret = DEFAULT_ERROR;
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
    strcpy(dst, "0");
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
      ret = DEFAULT_ERROR;
      break;
    }

    if (dst == NULL) {
      ret = DEFAULT_ERROR;
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
    if (obj) free(obj);
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
      ret = DEFAULT_ERROR;
      break;
    }

    if (rows < 0) {
      ret = DEFAULT_ERROR;
      break;
    }

    if (cols < 0) {
      ret = DEFAULT_ERROR;
      break;
    }

    if (dst == NULL) {
      ret = DEFAULT_ERROR;
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
    if (obj) free(obj);
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
  if (ptr == NULL) ret = DEFAULT_ERROR;

  /*
   * release memory
   */
  if (!ret) {
    if (ptr->tbl) {
      for (i = 0; i < ptr->rows; i++) {
        free(ptr->tbl[i]);
      }

      free(ptr->tbl);
    }

    free(ptr);
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
  if (ptr == NULL) ret = DEFAULT_ERROR;

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
 * 列の追加
 *
 * @param ptr   追加対象の行列オブジェクト
 * @param row   追加する行のデータ
 *
 * @return エラーコード(0で正常終了)
 */
int
cmat_add_row(cmat_t* ptr, double* src)
{
  int ret;
  double** tbl;
  double* dst;
  int capa;

  /*
   * initialize
   */
  ret = 0;
  tbl = NULL;
  dst = NULL;

  /*
   * argument check
   */
  do {
    if (ptr == NULL) {
      ret = DEFAULT_ERROR;
      break;
    }

    if (src == NULL) {
      ret = DEFAULT_ERROR;
      break;
    }
  } while (0);

  /*
   * alloc row
   */
  if (!ret) {
    dst = (double*)malloc(sizeof(double) * ptr->cols);
    if (dst == NULL) ret = DEFAULT_ERROR;
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
        ret = DEFAULT_ERROR;
      }
    }
  }

  /*
   * update context
   */
  if (!ret) {
    memcpy(dst, src, sizeof(double) * ptr->cols);
    ptr->tbl[ptr->rows++] = dst;
  }

  /*
   * post process
   */
  if (ret) {
    if (tbl) free(tbl);
    if (dst) free(dst);
  }

  return ret;
}

/**
 * 行列オブジェクトの転置
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

  /*
   * argument check
   */
  do {
    if (ptr == NULL) {
      ret = DEFAULT_ERROR;
      break;
    }

    if (dst == NULL) {
      ret = DEFAULT_ERROR;
      break;
    }
  } while (0);

  /*
   * alloc result object
   */
  if (!ret) {
    ret = alloc_object(ptr->cols, ptr->rows, ptr, &obj);
  }

  /*
   * put return parameter
   */
  if (!ret) {
    for (r = 0; r < ptr->rows; r++) {
      for (c = 0; c < ptr->cols; c++) {
        obj->tbl[c][r] = ptr->tbl[r][c];
      }
    }

    *dst = obj;
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
  double det;

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
      ret = DEFAULT_ERROR;
      break;
    }

    if (dst == NULL) {
      ret = DEFAULT_ERROR;
      break;
    }
  } while (0);
  
  /*
   * check if it's a regular matrix
   */
  if (!ret) {
    ret = cmat_det(ptr, &det);
  }

  if (!ret) {
    if (fabs(det) < 1.0e-6) ret = DEFAULT_ERROR;
  }

  /*
   * alloc result object
   */
  if (!ret) {
    ret = alloc_object(ptr->rows, ptr->cols, ptr, &obj);
  }

  /*
   * calculate inverse matrix
   */
  if (!ret) {
    ret = inverse_matrix(ptr->tbl, ptr->rows, obj->tbl);
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
    if (obj) cmat_destroy(obj);
  }

  return ret;
}

/**
 * 行列オブジェクトの積
 *
 * @param ptr   転置対象の行列オブジェクト
 * @param op    積行列
 * @param dst   演算結果の格納先
 *
 * @return エラーコード(0で正常終了)
 */
int
cmat_mul(cmat_t* ptr, cmat_t* op, cmat_t** dst)
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
      ret = DEFAULT_ERROR;
      break;
    }

    if (op == NULL) {
      ret = DEFAULT_ERROR;
      break;
    }

    if (dst == NULL) {
      ret = DEFAULT_ERROR;
      break;
    }
  } while (0);

  /*
   * check op value
   */
  if (!ret) {
    if (ptr->cols != op->rows) ret = DEFAULT_ERROR;
  }

  /*
   * alloc result object
   */
  if (!ret) {
    ret = alloc_object(ptr->rows, op->cols, ptr, &obj);
  }

  /*
   * set return parameter
   */
  if (!ret) {
    for (r = 0; r < ptr->rows; r++) {
      for (c = 0; c < op->cols; c++) {
        obj->tbl[r][c] = cdot(ptr->tbl, r, op->tbl, c, ptr->cols);
      }
    }

    *dst = obj;
  }

  return ret;
}

/**
 * 行列式の計算
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
      ret = DEFAULT_ERROR;
      break;
    }

    if (dst == NULL) {
      ret = DEFAULT_ERROR;
      break;
    }
  } while (0);

  /*
   * check shape
   */
  if (!ret) {
    if (ptr->rows != ptr->cols) ret = DEFAULT_ERROR;
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
      det = det_dim2(ptr->tbl);
      break;

    case 3:
      det = det_dim3(ptr->tbl);
      break;

    default:
      ret = det_matrix(ptr->tbl, ptr->rows, &det);
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

  /*
   * initialize
   */
  ret = 0;
  res = 0;

  /*
   * argument check
   */
  do {
    if (ptr == NULL) {
      ret = DEFAULT_ERROR;
      break;
    }

    if (val == NULL) {
      ret = DEFAULT_ERROR;
      break;
    }

    if (dst == NULL) {
      ret = DEFAULT_ERROR;
      break;
    }
  } while (0);

  /*
   * put return parameter
   */
  if (!ret) {
    for (r = 0; r < ptr->rows; r++) {
      for (c = 0; c < ptr->cols; c++) {
        if (fabs(ptr->tbl[r][c] - *val++) > ptr->coff) {
          res = !0;
          goto loop_out;
        }
      }
    }
    loop_out:

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
  if (ptr == NULL) ret = DEFAULT_ERROR;

  /*
   * update context
   */
  if (!ret) {
    ptr->coff = fabs(val);
  }

  return ret;
} 

