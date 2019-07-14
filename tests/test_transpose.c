#include <CUnit/CUnit.h>

#include <stdio.h>
#include "cmat.h"

void
test_transpose_1(void)
{
  cmat_t* m1;
  cmat_t* m2;

  double r1[] = {1, 2};
  double r2[] = {3, 4};
  double r3[] = {5, 6};
  double v[] = {1, 3, 5, 2, 4, 6};
  int res;

  cmat_new(2, &m1);
  cmat_append(m1, r1);
  cmat_append(m1, r2);
  cmat_append(m1, r3);
  cmat_transpose(m1, &m2);

#if 0
  printf("\n");
  cmat_print(m1, NULL);
  printf("\n");
  cmat_print(m2, NULL);
#endif

  cmat_check(m2, v, &res);
  CU_ASSERT(res == 0);

  cmat_destroy(m1);
  cmat_destroy(m2);
}

void
test_transpose_2(void)
{
  cmat_t* m1;
  cmat_t* m2;

  double v1[] = {
    1, 2, 3,
    4, 5, 6,
    7, 8, 9,
  };

  double v2[] = {
    1, 4, 7,
    2, 5, 8,
    3, 6, 9,
  };

  int res;

  cmat_new2(v1, 3, 3, &m1);
  cmat_transpose(m1, &m2);

#if 0
  printf("\n");
  cmat_print(m1, NULL);
  printf("\n");
  cmat_print(m2, NULL);
#endif

  cmat_check(m2, v2, &res);
  CU_ASSERT(res == 0);

  cmat_destroy(m1);
  cmat_destroy(m2);
}

void
test_transpose_3(void)
{
  cmat_t* m1;

  double src[] = {
    1, 2, 3,
    4, 5, 6
  };

  double mst[] = {
    1, 4,
    2, 5,
    3, 6
  };

  int res;

  cmat_new2(src, 2, 3, &m1);

#if 0
  printf("\n");
  cmat_print(m1, "original");
#endif

  cmat_transpose(m1, NULL);

#if 0
  printf("\n");
  cmat_print(m1, "transposed");
#endif

  cmat_check(m1, mst, &res);
  CU_ASSERT(res == 0);

  cmat_destroy(m1);
}


void
init_test_transpose()
{
  CU_pSuite suite;

  suite = CU_add_suite("transpose", NULL, NULL);
  CU_add_test(suite, "transpose#1", test_transpose_1);
  CU_add_test(suite, "transpose#2", test_transpose_2);
  CU_add_test(suite, "transpose#3", test_transpose_3);
}
