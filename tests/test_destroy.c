#include <CUnit/CUnit.h>

#include <stdio.h>
#include "cmat.h"
#include "test_add.h"

static void
test_normal_1(void)
{
  int err;

  cmat_t* m;

  cmat_new(10, &m);

  err = cmat_destroy(m);

  CU_ASSERT(err == 0);
}

static void
test_normal_2(void)
{
  int err;
  cmat_t* m;
  double v[] = {
    1, 2, 3,
    4, 5, 6,
  };

  cmat_new2(v, 2, 3, &m);

  err = cmat_destroy(m);

  CU_ASSERT(err == 0);
}

static void
test_error_1(void)
{
  int err;

  err = cmat_destroy(NULL);
  CU_ASSERT(err == CMAT_ERR_BADDR);
}

void
init_test_destroy()
{
  CU_pSuite suite;

  suite = CU_add_suite("destroy object", NULL, NULL);
  CU_add_test(suite, "destroy#1", test_normal_1);
  CU_add_test(suite, "destroy#1", test_normal_2);
  CU_add_test(suite, "destroy#E1", test_error_1);
}
