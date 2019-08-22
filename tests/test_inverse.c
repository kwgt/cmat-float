#include <CUnit/CUnit.h>

#include <stdio.h>
#include <stdlib.h>
#include "cmat.h"
#include "test_inverse.h"

#define N(x)        (sizeof(x) / sizeof(*x))

static int
create_matrix(const matrix_info_t* info, cmat_t** dst)
{
  return cmat_new(info->val, info->size, info->size, dst);
}

static void
check_matrix(cmat_t* org, cmat_t* inv, const matrix_info_t* info)
{
  int res1;
  int res2;
  int size;
  cmat_t* im;
  cmat_t* cm;
  int i;
  int n;
  float* row;

  n = info->size;

  /* check matrix values */
  cmat_check(inv, info->val, &res1);

#if 0
  if (res1) {
    cmat_t* m;
    cmat_t* d;
    float max;

    cmat_new(info->val, inv->rows, inv->cols, &m);
    cmat_sub(inv, m, &d);
    cmat_abs_max(d, &max);
    printf("%.20f\n", max);

    cmat_destroy(m);
    cmat_destroy(d);
  }
#endif

  /* check matrix product result */
  cmat_new(NULL, 0, n, &im);
  row = (float*)malloc(sizeof(float) * n);

  for (i = 0; i < n; i++) {
    memset(row, 0, sizeof(float) * n);
    row[i] = 1.0;

    cmat_append(im, row);
  }

  free(row);

  cmat_product(org, inv, &cm);
  cmat_compare(im, cm, &res2);

  cmat_destroy(im);
  cmat_destroy(cm);

  CU_ASSERT(res1 == 0);
  CU_ASSERT(inv->rows == n);
  CU_ASSERT(inv->cols == n);
  CU_ASSERT(res2 == 0);
}

static void
test_normal_1(void)
{
  int err;
  cmat_t* m1;
  cmat_t* m2;
  int i;

  for (i = 0; i < N(data); i++) {
    create_matrix(&data[i].op, &m1);

    err = cmat_inverse(m1, &m2);
    CU_ASSERT(err == 0);

#if 0
    printf("#%d\n", i);
    cmat_print(m1, NULL);

    printf("\n");
    cmat_print(m2, NULL);
    printf("\n");
#endif

    check_matrix(m1, m2, &data[i].ans);

    cmat_destroy(m1);
    if (err == 0) cmat_destroy(m2);
  }
}

static void
test_normal_2(void)
{
  int err;
  cmat_t* m1;
  cmat_t* m2;
  int i;

  for (i = 0; i < N(data); i++) {
    create_matrix(&data[i].op, &m1);
    create_matrix(&data[i].op, &m2);

#if 0
    printf("\n");
    cmat_print(m1, NULL);
#endif

    err = cmat_inverse(m1, NULL);
    CU_ASSERT(err == 0);

#if 0
    printf("\n");
    cmat_print(m1, NULL);
#endif

    check_matrix(m2, m1, &data[i].ans);

    cmat_destroy(m1);
    cmat_destroy(m2);
  }
}

static void
test_error_1(void)
{
  int err;
  cmat_t* m;

  m   = NULL;
  err = cmat_inverse(NULL, &m);

  CU_ASSERT(err == CMAT_ERR_BADDR);
  CU_ASSERT(m == NULL);
}

static void
test_error_2(void)
{
  int err;
  float v[] = {
    1, 2, 3,
    4, 5, 6
  };

  cmat_t* m;

  cmat_new(v, 2, 3, &m);

  err = cmat_inverse(m, NULL);
  cmat_destroy(m);

  CU_ASSERT(err == CMAT_ERR_SHAPE);
}

void
init_test_inverse()
{
  CU_pSuite suite;

  suite = CU_add_suite("inverse", NULL, NULL);
  CU_add_test(suite, "inverse#1", test_normal_1);
  CU_add_test(suite, "inverse#2", test_normal_2);
  CU_add_test(suite, "inverse#E1", test_error_1);
  CU_add_test(suite, "inverse#E2", test_error_2);
}
