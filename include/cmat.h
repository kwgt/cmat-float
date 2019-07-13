/*
 * Cheap maxtrix library for C
 *
 *  Copyright (C) 2019 Hiroshi Kuwagata <kgt9221@gmail.com>
 */

#ifndef __CHEAP_MATRIX_H__
#define __CHEAP_MATRIX_H__

typedef struct {
  double** tbl;
  int rows;
  int cols;

  int capa;
  double coff;  // as cutoff
} cmat_t;

int cmat_new(int n, cmat_t** dst);
int cmat_new2(double* src, int rows, int cols, cmat_t** dst);
int cmat_destroy(cmat_t* dst);
int cmat_add_row(cmat_t* ptr, double* r);
int cmat_transpose(cmat_t* ptr, cmat_t** dst);
int cmat_inverse(cmat_t* ptr, cmat_t** dst);
int cmat_mul(cmat_t* ptr, cmat_t* op, cmat_t** dst);
int cmat_det(cmat_t* ptr, double* dst);


int cmat_print(cmat_t* ptr, char* label);
int cmat_is_equal(cmat_t* ptr, cmat_t* op);
int cmat_check(cmat_t* ptr, double* val, int* dst);
int cmat_set_cutoff_threshold(cmat_t* ptr, double val);

#endif /* !defined(__CHEAP_MATRIX_H__) */
