/*
 * Cheap maxtrix library for C
 *
 *  Copyright (C) 2019 Hiroshi Kuwagata <kgt9221@gmail.com>
 */

#ifndef __CHEAP_MATRIX_H__
#define __CHEAP_MATRIX_H__

typedef struct {
  float* tbl;
  float** row;

  int rows;
  int cols;

  int capa;
  float coff;  // as cutoff
} cmat_t;

#define CMAT_ERR_NOMEM      -1    // NO MEMORY
#define CMAT_ERR_BADDR      -2    // BAD ADDRESS
#define CMAT_ERR_BSIZE      -3    // BAD SIZE
#define CMAT_ERR_INVAL      -4    // INVALID VALUE
#define CMAT_ERR_SHAPE      -5    // MATRIX SHAPE ERROR
#define CMAT_ERR_NREGL      -6    // NOT REGULAR MATRIX

#define CMAT_ROW(p,i)       ((p)->row[(i)])

int cmat_new(float* src, int rows, int cols, cmat_t** dst);
int cmat_clone(cmat_t* src, cmat_t** dst);
int cmat_destroy(cmat_t* ptr);
int cmat_append(cmat_t* ptr, float* r);

int cmat_add(cmat_t* ptr, cmat_t* op, cmat_t** dst);
int cmat_sub(cmat_t* ptr, cmat_t* op, cmat_t** dst);
int cmat_product(cmat_t* ptr, cmat_t* op, cmat_t** dst);
int cmat_mul(cmat_t* ptr, float op, cmat_t** dst);
int cmat_transpose(cmat_t* ptr, cmat_t** dst);
int cmat_det(cmat_t* ptr, float* dst);
int cmat_dot(cmat_t* ptr, cmat_t* op, float* dst);
int cmat_inverse(cmat_t* ptr, cmat_t** dst);
int cmat_lu_decomp(cmat_t* ptr, cmat_t** dst, int* piv);

int cmat_abs_max(cmat_t*ptr, float* dst);
int cmat_abs_min(cmat_t*ptr, float* dst);
int cmat_permute_row(cmat_t* ptr, int* piv);
int cmat_permute_column(cmat_t* ptr, int* piv);

int cmat_print(cmat_t* ptr, char* label);
int cmat_compare(cmat_t* ptr, cmat_t* op, int* dst);
int cmat_check(cmat_t* ptr, float* val, int* dst);
int cmat_set_cutoff_threshold(cmat_t* ptr, float val);

#endif /* !defined(__CHEAP_MATRIX_H__) */
