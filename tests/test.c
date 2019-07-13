#include <CUnit/CUnit.h>
#include <CUnit/Console.h>

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
  cmat_add_row(m1, r1);
  cmat_add_row(m1, r2);
  cmat_add_row(m1, r3);
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

int
main(int argc, char* argv[])
{
  CU_pSuite suite;

  CU_initialize_registry();

  suite = CU_add_suite("transpose", NULL, NULL);
  CU_add_test(suite, "transpose#1", test_transpose_1);
  CU_add_test(suite, "transpose#2", test_transpose_2);

  suite = CU_add_suite("mul", NULL, NULL);
  CU_add_test(suite, "mul#1", test_mul_1);

  suite = CU_add_suite("det", NULL, NULL);
  CU_add_test(suite, "det#1", test_det_1);
  CU_add_test(suite, "det#2", test_det_2);
  CU_add_test(suite, "det#3", test_det_3);
  CU_add_test(suite, "det#4", test_det_4);

  suite = CU_add_suite("inverse", NULL, NULL);
  CU_add_test(suite, "inverse#1", test_inverse_1);

  CU_console_run_tests();
  CU_cleanup_registry();

  return 0;
}
