#include <CUnit/CUnit.h>
#include <CUnit/Console.h>

extern void init_test_new();
extern void init_test_destroy();
extern void init_test_append();
extern void init_test_add();
extern void init_test_sub();
extern void init_test_mul();
extern void init_test_product();
extern void init_test_transpose();
extern void init_test_det();
extern void init_test_dot();
extern void init_test_inverse();

int
main(int argc, char* argv[])
{
  CU_initialize_registry();

  init_test_new();
  init_test_destroy();
  init_test_append();
  init_test_add();
  init_test_sub();
  init_test_mul();
  init_test_product();
  init_test_transpose();
  init_test_det();
  init_test_dot();
  init_test_inverse();

  CU_console_run_tests();
  CU_cleanup_registry();

  return 0;
}
