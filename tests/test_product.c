﻿#include <CUnit/CUnit.h>

#include <stdio.h>
#include "cmat.h"
#include "test_product.h"

#define N(x)        (sizeof(x) / sizeof(*x))

static int
create_matrix(const matrix_info_t* info, cmat_t** dst)
{
  return cmat_new(info->val, info->rows, info->cols, dst);
}

static void
check_matrix(cmat_t* ptr, const matrix_info_t* info)
{
  int res;

  cmat_check(ptr, info->val, &res);
  CU_ASSERT(res == 0);
  CU_ASSERT(ptr->rows == info->rows);
  CU_ASSERT(ptr->cols == info->cols);
}

static void
test_normal_1(void)
{
  cmat_t* m1;
  cmat_t* m2;
  cmat_t* m3;
  int i;

  for (i = 0; i < N(data); i++) {
    create_matrix(&data[i].op1, &m1);
    create_matrix(&data[i].op2, &m2);

    cmat_product(m1, m2, &m3);

#if 0
    printf("\n");
    cmat_print(m1, NULL);
    printf("\n");
    cmat_print(m2, NULL);
    printf("\n");
    cmat_print(m3, NULL);
#endif

    check_matrix(m3, &data[i].ans);

    cmat_destroy(m1);
    cmat_destroy(m2);
    cmat_destroy(m3);
  }
}

static void
test_normal_2(void)
{
  cmat_t* m1;
  cmat_t* m2;
  int i;

  for (i = 0; i < N(data); i++) {
    create_matrix(&data[i].op1, &m1);
    create_matrix(&data[i].op2, &m2);

#if 0
    printf("\n");
    cmat_print(m1, NULL);
    printf("\n");
    cmat_print(m2, NULL);
#endif

    cmat_product(m1, m2, NULL);

#if 0
    printf("\n");
    cmat_print(m1, NULL);
#endif

    check_matrix(m1, &data[i].ans);

    cmat_destroy(m1);
    cmat_destroy(m2);
  }
}

static void
test_error_1(void)
{
  cmat_t* m;
  int err;

  create_matrix(&data[0].op1, &m);

  err = cmat_product(NULL, m, NULL);
  cmat_destroy(m);

  CU_ASSERT(err == CMAT_ERR_BADDR);
}

static void
test_error_2(void)
{
  cmat_t* m;
  int err;

  create_matrix(&data[0].op1, &m);

  err = cmat_product(m, NULL, NULL);
  cmat_destroy(m);

  CU_ASSERT(err == CMAT_ERR_BADDR);
}

static void
test_error_3(void)
{
  float v[] = {
    1, 2, 3,
    4, 5, 6,
  };

  cmat_t* m1;
  cmat_t* m2;
  cmat_t* m3;
  int err;

  cmat_new(v, 2, 3, &m1);
  cmat_new(v, 2, 3, &m2);

  err = cmat_product(m1, m2, &m3);
  cmat_destroy(m1);
  cmat_destroy(m2);

  CU_ASSERT(err == CMAT_ERR_SHAPE);
}

void
init_test_product()
{
  CU_pSuite suite;

  suite = CU_add_suite("product", NULL, NULL);
  CU_add_test(suite, "product#1", test_normal_1);
  CU_add_test(suite, "product#2", test_normal_2);
  CU_add_test(suite, "product#E1", test_error_1);
  CU_add_test(suite, "product#E2", test_error_2);
  CU_add_test(suite, "product#E3", test_error_3);
}
