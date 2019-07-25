#include <CUnit/CUnit.h>

#include <stdio.h>
#include "cmat.h"
#include "test_mul.h"

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
  int i;

  for (i = 0; i < N(data); i++) {
    create_matrix(&data[i].op1, &m1);

    cmat_mul(m1, data[i].op2, &m2);

#if 0
    printf("\n");
    cmat_print(m1, NULL);
    printf("\n");
    cmat_print(m2, NULL);
    printf("\n");
    cmat_print(m3, NULL);
#endif

    check_matrix(m2, &data[i].ans);

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
    create_matrix(&data[i].op1, &m1);

#if 0
    printf("\n");
    cmat_print(m1, NULL);
    printf("\n");
    cmat_print(m2, NULL);
#endif

    cmat_mul(m1, data[i].op2, NULL);

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
  cmat_t* m;

  err = cmat_mul(NULL, 0.0, &m);

  CU_ASSERT(err == CMAT_ERR_BADDR);

  err = cmat_mul(NULL, 0.0, NULL);

  CU_ASSERT(err == CMAT_ERR_BADDR);
}

static void
test_error_2(void)
{
  int err;
  cmat_t* m;

  create_matrix(&data[0].op1, &m);

  err = cmat_mul(m, NAN, NULL);
  cmat_destroy(m);

  CU_ASSERT(err == CMAT_ERR_INVAL);
}

void
init_test_mul()
{
  CU_pSuite suite;

  suite = CU_add_suite("mul", NULL, NULL);
  CU_add_test(suite, "mul#1", test_normal_1);
  CU_add_test(suite, "mul#2", test_normal_2);
  CU_add_test(suite, "mul#E1", test_error_1);
  CU_add_test(suite, "mul#E2", test_error_2);
}
