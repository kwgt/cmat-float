#include <CUnit/CUnit.h>

#include <stdio.h>
#include <stdlib.h>
#include "cmat.h"
#include "test_lu_decomp.h"

#define N(x)          (sizeof(x) / sizeof(*x))
#define SWAP(a,b,t)   do {t c; c = (a); (a) = (b); (b) = c;} while(0)

static void
check_rule1(cmat_t* p, cmat_t* a, cmat_t* l, cmat_t* u)
{
  cmat_t* pa;
  cmat_t* lu;
  int res;

  /*
   * PA = LUをチェック
   */

  cmat_product(p, a, &pa);
  cmat_product(l, u, &lu);

  cmat_compare(pa, lu, &res);

  cmat_destroy(pa);
  cmat_destroy(lu);

  CU_ASSERT(res == 0);
}

static void
check_rule2(cmat_t* a, cmat_t* u, cmat_t* l, int* piv)
{
  cmat_t* ia;
  cmat_t* iu;
  cmat_t* il;
  cmat_t* out;
  int res;

  ia  = NULL;
  iu  = NULL;
  il  = NULL;
  out = NULL;

  /*
   * inv(A) = inv(L) * inv(U)をチェック
   */
  cmat_inverse(a, &ia);
  cmat_inverse(u, &iu);
  cmat_inverse(l, &il);

  if (ia && iu && il) {
    cmat_product(iu, il, &out);
    cmat_permute_column(out, piv);

    cmat_compare(ia, out, &res);
    CU_ASSERT(res == 0);

#if 0
    if (res) {
      cmat_t* d;

      cmat_sub(ia, out, &d);
      cmat_print(d, "diff");
      cmat_destroy(d);
    }
#endif
  }

  if (ia) cmat_destroy(ia);
  if (iu) cmat_destroy(iu);
  if (il) cmat_destroy(il);
  if (out) cmat_destroy(out);
}

static void
test_normal_1(void)
{
  int err;
  cmat_t* m1;
  cmat_t* m2;
  int res;
  int i;

  for (i = 0; i < N(data); i++) {
    cmat_new(data[i].op.val, data[i].op.size, data[i].op.size, &m1);

    err = cmat_lu_decomp(m1, &m2, NULL);
    CU_ASSERT(err == 0);

#if 0
    printf("#%d\n", i);
    cmat_print(m1, NULL);

    printf("\n");
    cmat_print(m2, NULL);
    printf("\n");
#endif

    cmat_check(m2, data[i].ans.val, &res);
    CU_ASSERT(res  == 0);

    cmat_destroy(m1);
    if (err == 0) cmat_destroy(m2);
  }
}

static void
test_normal_2(void)
{
  int err;
  cmat_t* a;
  cmat_t* o;
  cmat_t* l;
  cmat_t* u;
  cmat_t* p;

  int i;
  int j;
  int k;
  int sz;
  double* val;
  int res;
  int piv[20];

  for (i = 0; i < N(data); i++) {
    a   = NULL;
    o   = NULL;
    l   = NULL;
    u   = NULL;
    p   = NULL;

    val = data[i].op.val; 
    sz  = data[i].op.size;

    cmat_new(val, sz, sz, &a);

    err = cmat_lu_decomp(a, &o, piv);
    CU_ASSERT(err == 0);

    /*
     * LとUの分離を行う。
     * また同時に置換行列を作成する。
     */

    cmat_new(NULL, sz, sz, &l);
    cmat_new(NULL, sz, sz, &u);
    cmat_new(NULL, sz, sz, &p);

    for (j = 0; j < sz; j++) {
      for (k = 0; k < sz; k++) {
        if (k >= j) CMAT_ROW(u, j)[k] = CMAT_ROW(o, j)[k];
        if (k == j) CMAT_ROW(l, j)[k] = 1.0;
        if (k < j)  CMAT_ROW(l, j)[k] = CMAT_ROW(o, j)[k];
      }

      CMAT_ROW(p, j)[piv[j]] = 1.0;
    }

    check_rule1(p, a, l, u);
    check_rule2(a, u, l, piv);

    if (a) cmat_destroy(a);
    if (o) cmat_destroy(o);
    if (l) cmat_destroy(l);
    if (u) cmat_destroy(u);
    if (p) cmat_destroy(p);
  }
}

static void
test_normal_3(void)
{
  int err;
  cmat_t* m;
  int res;
  int i;

  for (i = 0; i < N(data); i++) {
    cmat_new(data[i].op.val, data[i].op.size, data[i].op.size, &m);

    err = cmat_lu_decomp(m, NULL, NULL);
    CU_ASSERT(err == 0);

#if 0
    printf("#%d\n", i);
    cmat_print(m1, NULL);

    printf("\n");
    cmat_print(m2, NULL);
    printf("\n");
#endif

    cmat_check(m, data[i].ans.val, &res);
    CU_ASSERT(res  == 0);

    cmat_destroy(m);
  }
}

static void
test_error_1(void)
{
  int err;
  cmat_t* m;
  int piv[5];

  err = cmat_lu_decomp(NULL, NULL, NULL);
  CU_ASSERT(err == CMAT_ERR_BADDR);

  err = cmat_lu_decomp(NULL, &m, NULL);
  CU_ASSERT(err == CMAT_ERR_BADDR);

  err = cmat_lu_decomp(NULL, &m, piv);
  CU_ASSERT(err == CMAT_ERR_BADDR);
}

#if 0
static void
test_error_2(void)
{
  int err;
  double v[] = {
    1, 2, 3,
    4, 5, 6
  };

  cmat_t* m;

  cmat_new(v, 2, 3, &m);

  err = cmat_inverse(m, NULL);
  cmat_destroy(m);

  CU_ASSERT(err == CMAT_ERR_SHAPE);
}
#endif

void
init_test_lu_decomp()
{
  CU_pSuite suite;

  suite = CU_add_suite("LU decomp", NULL, NULL);
  CU_add_test(suite, "LU decomp#1", test_normal_1);
  CU_add_test(suite, "LU decomp#2", test_normal_2);
  CU_add_test(suite, "LU decomp#3", test_normal_3);
  CU_add_test(suite, "LU decomp#E1", test_error_1);
  //CU_add_test(suite, "LU decomp#E2", test_error_2);
}
