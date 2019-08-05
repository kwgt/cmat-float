#include <CUnit/CUnit.h>

#include <stdio.h>
#include "cmat.h"
#include "test_transpose.h"

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
  int i;

  for (i = 0; i < N(data); i++) {
    create_matrix(&data[i].op, &m1);

    cmat_transpose(m1, &m2);

#if 0
    printf("\n");
    cmat_print(m1, NULL);
#endif

    check_matrix(m2, &data[i].ans);

#if 0
    printf("\n");
    cmat_print(m2, NULL);
#endif

    cmat_destroy(m1);
    cmat_destroy(m2);
  }
}

static void
test_normal_2(void)
{
  cmat_t* m1;
  int i;

  for (i = 0; i < N(data); i++) {
    create_matrix(&data[i].op, &m1);

#if 0
    printf("\n");
    cmat_print(m1, NULL);
#endif

    cmat_transpose(m1, NULL);

#if 0
    printf("\n");
    cmat_print(m1, NULL);
#endif

    check_matrix(m1, &data[i].ans);

    cmat_destroy(m1);
  }
}

static void
test_error_1(void)
{
  int err;

  err = cmat_transpose(NULL, NULL);
  CU_ASSERT(err == CMAT_ERR_BADDR);
}

static void
test_error_2(void)
{
  cmat_t* m;
  int err;

  create_matrix(&data[0].op, &m);

#if 0
  printf("\n");
  cmat_print(m, NULL);
#endif

  err = cmat_transpose(NULL, &m);
  cmat_destroy(m);

  CU_ASSERT(err == CMAT_ERR_BADDR);
}

void
init_test_transpose()
{
  CU_pSuite suite;

  suite = CU_add_suite("transpose", NULL, NULL);
  CU_add_test(suite, "transpose#1", test_normal_1);
  CU_add_test(suite, "transpose#2", test_normal_2);
  CU_add_test(suite, "transpose#E1", test_error_1);
  CU_add_test(suite, "transpose#E2", test_error_2);
}
