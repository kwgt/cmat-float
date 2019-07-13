#include <CUnit/CUnit.h>

#include <stdio.h>
#include "cmat.h"

void
test_mul_1(void)
{
  cmat_t* m1;
  double v1[] = {1, 2, 3, 4};

  cmat_t* m2;
  double v2[] = {5, 6, 7, 8};

  cmat_t* m3;
  double v3[] = {19, 22, 43, 50};
  int res;

  cmat_new2(v1, 2, 2, &m1);
  cmat_new2(v2, 2, 2, &m2);
  cmat_mul(m1, m2, &m3);

#if 0
  printf("\n");
  cmat_print(m1, NULL);
  printf("\n");
  cmat_print(m2, NULL);
  printf("\n");
  cmat_print(m3, NULL);
#endif

  cmat_check(m3, v3, &res);
  CU_ASSERT(res == 0);

  cmat_destroy(m1);
  cmat_destroy(m2);
  cmat_destroy(m3);
}
