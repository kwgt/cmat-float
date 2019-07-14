#include <CUnit/CUnit.h>

#include <stdio.h>
#include "cmat.h"

void
test_inverse_1(void)
{
  cmat_t* m1;
  cmat_t* m2;
  cmat_t* m3;
  int res;

  {
    double v[] = {
       1,  2,  0, -1,
      -1,  1,  2, -0,
       2,  0,  1,  1,
       1, -2, -1,  1
    };

    double mst[] = {
      1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 1,
    };

    cmat_new2(v, 4, 4, &m1);
    cmat_inverse(m1, &m2);
    cmat_mul(m1, m2, &m3);
    cmat_check(m3, mst, &res);

#if 1
    printf("\n");
    cmat_print(m2, "inverse");
    printf("\n");
    cmat_print(m3, "check");
#endif

    cmat_destroy(m1);
    cmat_destroy(m2);
    cmat_destroy(m3);

    CU_ASSERT(res == 0);
  }
}

void
init_test_inverse()
{
  CU_pSuite suite;

  suite = CU_add_suite("inverse", NULL, NULL);
  CU_add_test(suite, "inverse#1", test_inverse_1);
}
