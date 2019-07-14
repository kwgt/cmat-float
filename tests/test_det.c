#include <CUnit/CUnit.h>

#include <stdio.h>
#include "cmat.h"

void
test_det_1(void)
{
  cmat_t* m;
  double v[] = {1};
  double det;

  cmat_new2(v, 1, 1, &m);
  cmat_det(m, &det);

  CU_ASSERT(det == 1.0);

  cmat_destroy(m);
}

void
test_det_2(void)
{
  cmat_t* m;
  double v[] = {2, -6, -1, 3};
  double det;

  cmat_new2(v, 2, 2, &m);
  cmat_det(m, &det);

  CU_ASSERT(det == 0.0);

  cmat_destroy(m);
}

void
test_det_3(void)
{
  cmat_t* m;
  double v[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  double det;

  cmat_new2(v, 3, 3, &m);
  cmat_det(m, &det);

  CU_ASSERT(det == 0.0);

  cmat_destroy(m);
}

void
test_det_4(void)
{
  cmat_t* m;
  double det;

  {
    double v[] = {
       2, -2,  4,  2,
       2, -1,  6,  3,
       3, -2, 12, 12,
      -1,  3, -4,  4
    };

    cmat_new2(v, 4, 4, &m);
    cmat_det(m, &det);
    cmat_destroy(m);

    CU_ASSERT(det == 120.0);
  }

  {
    double v[] = {
       2,  1,  3,  2,
       1,  0,  4,  0,
       3,  2,  0,  2,
      -2, -1,  1,  3,
    };

    cmat_new2(v, 4, 4, &m);
    cmat_det(m, &det);
    cmat_destroy(m);

    CU_ASSERT(det == 2);
  }

  {
    double v[] = {
       3, 1, 2, 5,
      -1, 1, 3, 6,
       4, 0, 2, 1,
       5, 1, 0, 4
    };

    cmat_new2(v, 4, 4, &m);
    cmat_det(m, &det);
    cmat_destroy(m);

    // printf("%f\n", det);

    CU_ASSERT(fabs(det- 6.0) < 0.00001);
  }
}

void
init_test_det()
{
  CU_pSuite suite;

  suite = CU_add_suite("det", NULL, NULL);

  CU_add_test(suite, "det#1", test_det_1);
  CU_add_test(suite, "det#2", test_det_2);
  CU_add_test(suite, "det#3", test_det_3);
  CU_add_test(suite, "det#4", test_det_4);
}
