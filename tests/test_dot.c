#include <CUnit/CUnit.h>

#include <stdio.h>
#include "cmat.h"
#include "test_dot.h"

#define N(x)        (sizeof(x) / sizeof(*x))

static int
create_matrix(const matrix_info_t* info, cmat_t** dst)
{
  return cmat_new2(info->val, info->rows, info->cols, dst);
}

static void
test_normal_1(void)
{
  cmat_t* m1;
  cmat_t* m2;
  double dot;
  int i;

  for (i = 0; i < N(data); i++) {
    create_matrix(&data[i].op1, &m1);
    create_matrix(&data[i].op2, &m2);

    cmat_dot(m1, m2, &dot);

#if 0
    printf("\n");
    cmat_print(m1, NULL);
    printf("\n");
    cmat_print(m2, NULL);
    printf("\n");
    printf("%f\n", dot);
#endif

    CU_ASSERT(dot == data[i].ans);

    cmat_destroy(m1);
    cmat_destroy(m2);
  }
}

void
init_test_dot()
{
  CU_pSuite suite;

  suite = CU_add_suite("dot product", NULL, NULL);
  CU_add_test(suite, "dot#1", test_normal_1);
}
