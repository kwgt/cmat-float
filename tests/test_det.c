#include <CUnit/CUnit.h>

#include <stdio.h>
#include <math.h>
#include "cmat.h"
#include "test_det.h"

#define N(x)        (sizeof(x) / sizeof(*x))

static int
create_matrix(const matrix_info_t* info, cmat_t** dst)
{
  return cmat_new(info->val, info->rows, info->cols, dst);
}

static int
check_det(float a, float b)
{
  int s1;
  float v1;
  int e1;

  int s2;
  float v2;
  int e2;

  s1 = signbit(a);
  v1 = fabs(frexp(a, &e1));

  s2 = signbit(b);
  v2 = fabs(frexp(b, &e2));

  // 行列式の計算結果が非常に大きな値になる場合があるので
  // 数値そのものの一致では無く、有効桁数で確認を行っている
  return s1 == s2 && e1 == e2 && fabs(v1 - v2) < 1e-5;
}

static void
test_normal_1(void)
{
  cmat_t* m1;
  float det;
  int i;

  for (i = 0; i < N(data); i++) {
    create_matrix(&data[i].op, &m1);

    cmat_det(m1, &det);

    CU_ASSERT(check_det(det, data[i].ans));

#if 0
    if (!check_det(det, data[i].ans)) {
      printf("\n");
      cmat_print(m1, NULL);
      printf("%f %f\n", det, data[i].ans);
      exit(0);
    }
#endif

    cmat_destroy(m1);
  }
}

static void
test_error_1(void)
{
  float det;
  int err;

  det = NAN;

  err = cmat_det(NULL, &det);
  CU_ASSERT(err == CMAT_ERR_BADDR);
  CU_ASSERT(isnan(det));
}

static void
test_error_2(void)
{
  cmat_t* m;
  int err;

  create_matrix(&data[0].op, &m);

  err = cmat_det(m, NULL);
  cmat_destroy(m);

  CU_ASSERT(err == CMAT_ERR_BADDR);
}


void
init_test_det()
{
  CU_pSuite suite;

  suite = CU_add_suite("det", NULL, NULL);
  CU_add_test(suite, "det#1", test_normal_1);
  CU_add_test(suite, "det#E1", test_error_1);
  CU_add_test(suite, "det#E2", test_error_2);
}
