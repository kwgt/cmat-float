#include <CUnit/CUnit.h>

#include <stdio.h>
#include "cmat.h"

void
test_sub_1(void)
{
  cmat_t* m1;
  double v1[] = {1, 2, 3, 4};

  cmat_t* m2;
  double v2[] = {1, 2, 3, 4};

  cmat_t* m3;

  double mst[] = {0, 0, 0, 0};
  int res;

  cmat_new2(v1, 2, 2, &m1);
  cmat_new2(v2, 2, 2, &m2);
  cmat_sub(m1, m2, &m3);

#if 0
  printf("\n");
  cmat_print(m1, NULL);
  printf("\n");
  cmat_print(m2, NULL);
  printf("\n");
  cmat_print(m3, NULL);
#endif

  cmat_check(m3, mst, &res);
  CU_ASSERT(res == 0);

  cmat_destroy(m1);
  cmat_destroy(m2);
  cmat_destroy(m3);
}

void
test_sub_2(void)
{
  cmat_t* m1;
  double v1[] = {1, 2, 3, 4};

  cmat_t* m2;
  double v2[] = {1, 2, 3, 4};

  double mst[] = {0, 0, 0, 0};
  int res;

  cmat_new2(v1, 2, 2, &m1);
  cmat_new2(v2, 2, 2, &m2);
  cmat_sub(m1, m2, NULL);

#if 0
  printf("\n");
  cmat_print(m1, NULL);
  printf("\n");
  cmat_print(m2, NULL);
  printf("\n");
  cmat_print(m3, NULL);
#endif

  cmat_check(m1, mst, &res);
  CU_ASSERT(res == 0);

  cmat_destroy(m1);
  cmat_destroy(m2);
}

void
init_test_sub()
{
  CU_pSuite suite;

  suite = CU_add_suite("sub", NULL, NULL);
  CU_add_test(suite, "sub#1", test_sub_1);
  CU_add_test(suite, "sub#2", test_sub_2);
}
