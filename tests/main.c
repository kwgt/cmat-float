#include <CUnit/CUnit.h>
#include <CUnit/Console.h>

extern void init_test_add();
extern void init_test_sub();
extern void init_test_mul();
extern void init_test_transpose();
extern void init_test_det();
extern void init_test_inverse();

int
main(int argc, char* argv[])
{
  CU_initialize_registry();

  init_test_add();
  init_test_sub();
  init_test_mul();
  init_test_transpose();
  init_test_det();
  init_test_inverse();

  CU_console_run_tests();
  CU_cleanup_registry();

  return 0;
}
