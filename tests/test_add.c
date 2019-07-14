#include <CUnit/CUnit.h>

#include <stdio.h>
#include "cmat.h"

void
test_add_1(void)
{
  cmat_t* m1;
  double v1[] = {1, 2, 3, 4};

  cmat_t* m2;
  double v2[] = {1, 2, 3, 4};

  cmat_t* m3;

  double mst[] = {2, 4, 6, 8};
  int res;

  cmat_new2(v1, 2, 2, &m1);
  cmat_new2(v2, 2, 2, &m2);
  cmat_add(m1, m2, &m3);

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
test_add_2(void)
{
  cmat_t* m1;
  double v1[] = {1, 2, 3, 4};

  cmat_t* m2;
  double v2[] = {1, 2, 3, 4};

  double mst[] = {2, 4, 6, 8};
  int res;

  cmat_new2(v1, 2, 2, &m1);
  cmat_new2(v2, 2, 2, &m2);
  cmat_add(m1, m2, NULL);

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
init_test_add()
{
  CU_pSuite suite;

  suite = CU_add_suite("add", NULL, NULL);
  CU_add_test(suite, "add#1", test_add_1);
  CU_add_test(suite, "add#2", test_add_2);
}
