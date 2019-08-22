#include <CUnit/CUnit.h>

#include <stdio.h>
#include "cmat.h"
#include "test_add.h"

static void
test_normal_1(void)
{
  int err;

  cmat_t* m1;
  cmat_t* m2;

  m1  = NULL;
  m2  = NULL;

  cmat_new(NULL, 0, 10, &m1);
  err = cmat_clone(m1, &m2);

  CU_ASSERT(err == 0);
  CU_ASSERT(m2 != NULL);

  if (m2 != NULL) {
    CU_ASSERT(m2->rows == m1->rows);
    CU_ASSERT(m2->cols == m1->cols);
  }

  cmat_destroy(m1);
  cmat_destroy(m2);
}

static void
test_normal_2(void)
{
  int err;
  cmat_t* m1;
  cmat_t* m2;
  float v[] = {
    1, 2, 3,
    4, 5, 6,
  };

  m1  = NULL;
  m2  = NULL;

  cmat_new(v, 2, 3, &m1);
  err = cmat_clone(m1, &m2);

  CU_ASSERT(err == 0);
  CU_ASSERT(m2 != NULL);

  if (m2 != NULL) {
    CU_ASSERT(m2->rows == m1->rows);
    CU_ASSERT(m2->cols == m1->cols);
    CU_ASSERT(!memcmp(m1->tbl, m2->tbl, sizeof(float) * m2->rows * m2->cols));
  }

  cmat_destroy(m1);
  cmat_destroy(m2);
}

static void
test_error_e1(void)
{
  int err;
  cmat_t* m;

  err = cmat_clone(NULL, &m);
  CU_ASSERT(err == CMAT_ERR_BADDR);
}

static void
test_error_e2(void)
{
  int err;
  cmat_t* m;
  float v[] = {
    1, 2, 3,
    4, 5, 6,
  };

  cmat_new(v, 2, 3, &m);

  err = cmat_clone(m, NULL);
  CU_ASSERT(err == CMAT_ERR_BADDR);

  cmat_destroy(m);
}

void
init_test_clone()
{
  CU_pSuite suite;

  suite = CU_add_suite("clone object", NULL, NULL);
  CU_add_test(suite, "clone#1", test_normal_1);
  CU_add_test(suite, "clone#2", test_normal_2);
  CU_add_test(suite, "new#E1", test_error_e1);
  CU_add_test(suite, "new#E2", test_error_e2);
}
