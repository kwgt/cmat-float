﻿#include <CUnit/CUnit.h>

#include <stdio.h>
#include "cmat.h"
#include "test_sub.h"

#define N(x)        (sizeof(x) / sizeof(*x))

static int
create_matrix(const matrix_info_t* info, cmat_t** dst)
{
  return cmat_new2(info->val, info->rows, info->cols, dst);
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

    cmat_sub(m1, m2, &m3);

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

    cmat_sub(m1, m2, NULL);

#if 0
    printf("\n");
    cmat_print(m1, NULL);
#endif

    check_matrix(m1, &data[i].ans);

    cmat_destroy(m1);
    cmat_destroy(m2);
  }
}



void
init_test_sub()
{
  CU_pSuite suite;

  suite = CU_add_suite("sub", NULL, NULL);
  CU_add_test(suite, "sub#1", test_normal_1);
  CU_add_test(suite, "sub#2", test_normal_2);
}