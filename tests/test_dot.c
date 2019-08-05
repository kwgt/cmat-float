#include <CUnit/CUnit.h>

#include <stdio.h>
#include "cmat.h"
#include "test_dot.h"

#define N(x)        (sizeof(x) / sizeof(*x))

static int
create_matrix(const matrix_info_t* info, cmat_t** dst)
{
  return cmat_new(info->val, info->rows, info->cols, dst);
}

static void
test_normal_1(void)
{
  int err;
  cmat_t* m1;
  cmat_t* m2;
  double dot;
  int i;

  for (i = 0; i < N(data); i++) {
    create_matrix(&data[i].op1, &m1);
    create_matrix(&data[i].op2, &m2);

    err = cmat_dot(m1, m2, &dot);

#if 0
    printf("\n");
    cmat_print(m1, NULL);
    printf("\n");
    cmat_print(m2, NULL);
    printf("\n");
    printf("%f\n", dot);
#endif

    CU_ASSERT(err == 0);
    CU_ASSERT(dot == data[i].ans);

    cmat_destroy(m1);
    cmat_destroy(m2);
  }
}

static void
test_error_1(void)
{
  int err;
  cmat_t* m;
  double dot;
  int i;

  create_matrix(&data[0].op1, &m);

  err = cmat_dot(NULL, m, &dot);
  cmat_destroy(m);

  CU_ASSERT(err == CMAT_ERR_BADDR);
}

static void
test_error_2(void)
{
  cmat_t* m1;
  cmat_t* m2;
  double dot;
  int err;

  create_matrix(&data[0].op1, &m1);
  create_matrix(&data[0].op2, &m2);

  err = cmat_dot(m1, m2, NULL);
  cmat_destroy(m1);
  cmat_destroy(m2);

  CU_ASSERT(err == CMAT_ERR_BADDR);
}

static void
test_error_3(void)
{
  double v1[] = {
    1, 2,
    3, 4
  };

  cmat_t* m1;

  double v2[] = {
    1, 2, 3,
    4, 5, 6
  };

  cmat_t* m2;
  double dot;
  int err;

  cmat_new(v1, 2, 2, &m1);
  cmat_new(v1, 2, 3, &m2);

  err = cmat_dot(m1, m2, &dot);
  cmat_destroy(m1);
  cmat_destroy(m2);

  CU_ASSERT(err == CMAT_ERR_SHAPE);
}

void
init_test_dot()
{
  CU_pSuite suite;

  suite = CU_add_suite("dot product", NULL, NULL);
  CU_add_test(suite, "dot#1", test_normal_1);
  CU_add_test(suite, "dot#E1", test_error_1);
  CU_add_test(suite, "dot#E2", test_error_2);
  CU_add_test(suite, "dot#E3", test_error_3);
}
