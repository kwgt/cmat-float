#include <CUnit/CUnit.h>
#include <CUnit/Console.h>

#include "test.h"

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
