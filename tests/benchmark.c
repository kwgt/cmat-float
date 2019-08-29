#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "cmat.h"

typedef struct {
  long long tm;
  struct timespec t0;
} tmmes_t;

void
start_measure(tmmes_t* tm)
{
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &tm->t0);
}

void
stop_measure(tmmes_t* tm)
{
  struct timespec t1;

  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t1);

  tm->tm = (((t1.tv_sec - tm->t0.tv_sec) * 1000000000LL) +
                              (t1.tv_nsec - tm->t0.tv_nsec));
}

void
bench_add(tmmes_t* tm)
{
  cmat_t* m1;
  cmat_t* m2;
  int i;

  cmat_new(NULL, 1000, 1000, &m1);
  cmat_new(NULL, 1000, 1000, &m2);

  start_measure(tm);

  for (i = 0; i < 200; i++) {
    cmat_add(m1, m2, NULL);
  }

  stop_measure(tm);

  cmat_destroy(m1);
  cmat_destroy(m2);
}

void
bench_sub(tmmes_t* tm)
{
  cmat_t* m1;
  cmat_t* m2;
  int i;

  cmat_new(NULL, 1000, 1000, &m1);
  cmat_new(NULL, 1000, 1000, &m2);

  start_measure(tm);

  for (i = 0; i < 200; i++) {
    cmat_sub(m1, m2, NULL);
  }

  stop_measure(tm);

  cmat_destroy(m1);
  cmat_destroy(m2);
}

void
bench_mul(tmmes_t* tm)
{
  cmat_t* m1;
  int i;

  cmat_new(NULL, 1000, 1000, &m1);

  start_measure(tm);

  for (i = 0; i < 200; i++) {
    cmat_mul(m1, 2.0, NULL);
  }

  stop_measure(tm);

  cmat_destroy(m1);
}

void
bench_product(tmmes_t* tm)
{
  cmat_t* m1;
  cmat_t* m2;
  int i;

  cmat_new(NULL, 100, 100, &m1);
  cmat_new(NULL, 100, 100, &m2);

  start_measure(tm);

  for (i = 0; i < 200; i++) {
    cmat_product(m1, m2, NULL);
  }

  stop_measure(tm);

  cmat_destroy(m1);
  cmat_destroy(m2);
}

void
bench_det(tmmes_t* tm)
{
  cmat_t* m;
  float det;
  int i;
  int j;
  float* row;

  cmat_new(NULL, 100, 100, &m);

  srand(0);
  for (i = 0; i < 100; i++) {
    row = CMAT_ROW(m, i);

    for (j = 0; j < 100; j++) {
      row[j] = (float)rand() / RAND_MAX;
    }
  }

  start_measure(tm);

  for (i = 0; i < 2000; i++) {
    cmat_det(m, &det);
  }

  stop_measure(tm);

  cmat_destroy(m);
}

void
bench_inverse(tmmes_t* tm)
{
  cmat_t* m1;
  cmat_t* m2;
  int i;
  int j;
  float* row;

  cmat_new(NULL, 100, 100, &m1);

  srand(0);
  for (i = 0; i < 100; i++) {
    row = CMAT_ROW(m1, i);

    for (j = 0; j < 100; j++) {
      row[j] = (float)rand() / RAND_MAX;
    }
  }

  start_measure(tm);

  for (i = 0; i < 2000; i++) {
    cmat_inverse(m1, &m2);
    cmat_destroy(m2);
  }

  stop_measure(tm);

  cmat_destroy(m1);
}



int
main(int argc, char* argv[])
{
  tmmes_t tm;

  bench_add(&tm);
  printf("add     %10fmsec\n", tm.tm / 1000000.0);

  bench_sub(&tm);
  printf("sub     %10fmsec\n", tm.tm / 1000000.0);

  bench_mul(&tm);
  printf("mul     %10fmsec\n", tm.tm / 1000000.0);

  bench_product(&tm);
  printf("product %10fmsec\n", tm.tm / 1000000.0);

  bench_det(&tm);
  printf("det     %10fmsec\n", tm.tm / 1000000.0);

  bench_inverse(&tm);
  printf("inverse %10fmsec\n", tm.tm / 1000000.0);

  return 0;
}
