#include <CUnit/CUnit.h>

#include <stdio.h>
#include <stdlib.h>
#include "cmat.h"

#define N(x)          (sizeof(x) / sizeof(*x))
#define SWAP(a,b,t)   do {t c; c = (a); (a) = (b); (b) = c;} while(0)

static void
test_normal_1(void)
{
  int err;
  cmat_t* m;
  int res;

  int piv[4] = {3, 2, 1, 0};

  {
    double v[] = {
      1, 1, 1, 1,
      2, 2, 2, 2,
      3, 3, 3, 3,
      4, 4, 4, 4,
    };

    cmat_new(v, 4, 4, &m);
  }

#if 0
  printf("\n");
  cmat_print(m, "orig");
#endif

  err = cmat_permute_row(m, piv);
  CU_ASSERT(err == 0);

#if 0
  printf("\n");
  cmat_print(m, "permuted");
#endif

  {
    double v[] = {
      4, 4, 4, 4,
      3, 3, 3, 3,
      2, 2, 2, 2,
      1, 1, 1, 1,
    };

    cmat_check(m, v, &res);
    CU_ASSERT(res == 0);
  }

  CU_ASSERT(piv[0] == 3);
  CU_ASSERT(piv[1] == 2);
  CU_ASSERT(piv[2] == 1);
  CU_ASSERT(piv[3] == 0);

  cmat_destroy(m);
}

static void
test_normal_2(void)
{
  int err;
  cmat_t* m;
  int res;

  int piv[4] = {1, 3, 2, 0};

  {
    double v[] = {
      1, 1, 1, 1,
      2, 2, 2, 2,
      3, 3, 3, 3,
      4, 4, 4, 4,
    };

    cmat_new(v, 4, 4, &m);
  }

#if 0
  printf("\n");
  cmat_print(m, "orig");
#endif

  err = cmat_permute_row(m, piv);
  CU_ASSERT(err == 0);

#if 0
  printf("\n");
  cmat_print(m, "permuted");
#endif

  {
    double v[] = {
      4, 4, 4, 4,
      1, 1, 1, 1,
      3, 3, 3, 3,
      2, 2, 2, 2,
    };

    cmat_check(m, v, &res);
    CU_ASSERT(res == 0);
  }

  CU_ASSERT(piv[0] == 1);
  CU_ASSERT(piv[1] == 3);
  CU_ASSERT(piv[2] == 2);
  CU_ASSERT(piv[3] == 0);

  cmat_destroy(m);
}

static void
test_error_1(void)
{
  int err;

  int piv[4] = {1, 3, 2, 0};

  err = cmat_permute_row(NULL, piv);
  CU_ASSERT(err == CMAT_ERR_BADDR);
}

static void
test_error_2(void)
{
  int err;
  cmat_t* m;
  int res;

  {
    double v[] = {
      1, 1, 1, 1,
      2, 2, 2, 2,
      3, 3, 3, 3,
      4, 4, 4, 4,
    };

    cmat_new(v, 4, 4, &m);
  }

  err = cmat_permute_row(m, NULL);
  CU_ASSERT(err == CMAT_ERR_BADDR);

  cmat_destroy(m);
}

static void
test_error_3(void)
{
  int err;
  cmat_t* m;
  int res;

  {
    double v[] = {
      1, 1, 1, 1,
      2, 2, 2, 2,
      3, 3, 3, 3,
      4, 4, 4, 4,
    };

    cmat_new(v, 4, 4, &m);
  }

  {
    int piv[4] = {0, 0, 0, 0};

    err = cmat_permute_row(m, piv);
    CU_ASSERT(err == CMAT_ERR_INVAL);
  }

  {
    int piv[4] = {1, 2, 3, 4};

    err = cmat_permute_row(m, piv);
    CU_ASSERT(err == CMAT_ERR_INVAL);
  }

  {
    int piv[4] = {0, 1, 2, 0};

    err = cmat_permute_row(m, piv);
    CU_ASSERT(err == CMAT_ERR_INVAL);
  }

  {
    int piv[4] = {1, 1, 2, 2};

    err = cmat_permute_row(m, piv);
    CU_ASSERT(err == CMAT_ERR_INVAL);
  }

  cmat_destroy(m);
}


void
init_test_permute_row()
{
  CU_pSuite suite;

  suite = CU_add_suite("Permute row", NULL, NULL);
  CU_add_test(suite, "Permute row#1", test_normal_1);
  CU_add_test(suite, "Permute row#2", test_normal_2);
  CU_add_test(suite, "Permute row#E1", test_error_1);
  CU_add_test(suite, "Permute row#E2", test_error_2);
  CU_add_test(suite, "Permute row#E3", test_error_3);
}
