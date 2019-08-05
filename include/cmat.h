/*
 * Cheap maxtrix library for C
 *
 *  Copyright (C) 2019 Hiroshi Kuwagata <kgt9221@gmail.com>
 */

#ifndef __CHEAP_MATRIX_H__
#define __CHEAP_MATRIX_H__

typedef struct {
  double* tbl;
  int rows;
  int cols;

  int capa;
  double coff;  // as cutoff
} cmat_t;

#define CMAT_ERR_NOMEM      -1
#define CMAT_ERR_BADDR      -2
#define CMAT_ERR_BSIZE      -3
#define CMAT_ERR_INVAL      -4
#define CMAT_ERR_SHAPE      -5

int cmat_new(double* src, int rows, int cols, cmat_t** dst);
int cmat_destroy(cmat_t* dst);
int cmat_append(cmat_t* ptr, double* r);

int cmat_add(cmat_t* ptr, cmat_t* op, cmat_t** dst);
int cmat_sub(cmat_t* ptr, cmat_t* op, cmat_t** dst);
int cmat_product(cmat_t* ptr, cmat_t* op, cmat_t** dst);
int cmat_mul(cmat_t* ptr, double op, cmat_t** dst);
int cmat_transpose(cmat_t* ptr, cmat_t** dst);
int cmat_det(cmat_t* ptr, double* dst);
int cmat_dot(cmat_t* ptr, cmat_t* op, double* dst);
int cmat_inverse(cmat_t* ptr, cmat_t** dst);

int cmat_print(cmat_t* ptr, char* label);
int cmat_compare(cmat_t* ptr, cmat_t* op, int* dst);
int cmat_check(cmat_t* ptr, double* val, int* dst);
int cmat_set_cutoff_threshold(cmat_t* ptr, double val);

#endif /* !defined(__CHEAP_MATRIX_H__) */
