#include <CUnit/CUnit.h>

#include <stdio.h>
#include "cmat.h"

static void
test_normal_1(void)
{
  int err;
  cmat_t* m;

  float r1[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  float r2[] = {10, 9, 8, 7, 6, 5, 4, 3, 2, 1};

  m = NULL;

  cmat_new(NULL, 0, 10, &m);

  err = cmat_append(m, r1);
  CU_ASSERT(err == 0);
  CU_ASSERT(m->rows == 1);

#if 0
  cmat_print(m, NULL);
#endif

  err = cmat_append(m, r2);
  CU_ASSERT(err == 0);
  CU_ASSERT(m->rows == 2);

#if 0
  cmat_print(m, NULL);
#endif

  cmat_destroy(m);
}

static void
test_normal_2(void)
{
  int err;
  cmat_t* m;
  float v[] = {
    1, 2, 3,
    4, 5, 6,
  };
  float r1[] = { 7,  8,  9};
  float r2[] = {10, 11, 12};

  m = NULL;

  cmat_new(v, 2, 3, &m);

  err = cmat_append(m, r1);

#if 0
  cmat_print(m, NULL);
#endif

  CU_ASSERT(err == 0);
  CU_ASSERT(m->rows == 3);

  err = cmat_append(m, r2);

#if 0
  cmat_print(m, NULL);
#endif

  CU_ASSERT(err == 0);
  CU_ASSERT(m->rows == 4);

  cmat_destroy(m);
}

static void
test_error_1(void)
{
  int err;
  float r[] = {1, 2, 3};

  err = cmat_append(NULL, r);
  CU_ASSERT(err == CMAT_ERR_BADDR);
}

static void
test_error_2(void)
{
  int err;
  cmat_t* m;

  cmat_new(NULL, 0, 10, &m);

  err = cmat_append(m, NULL);
  CU_ASSERT(err == CMAT_ERR_BADDR);
}

void
init_test_append()
{
  CU_pSuite suite;

  suite = CU_add_suite("append row", NULL, NULL);
  CU_add_test(suite, "append#1", test_normal_1);
  CU_add_test(suite, "append#2", test_normal_2);
  CU_add_test(suite, "append#E1", test_error_1);
  CU_add_test(suite, "append#E2", test_error_2);
}
