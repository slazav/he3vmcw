extern "C" {
#include "grad.h"
#include "check.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

const char *xyz[3]={"x","y","z"};


// check if absolute value is less then E
void checkz(const char *name, double V1, double E, bool nl){
  printf(" [%s", name);
  if (fabs(V1) > E){
    printf("\nFAIL: V1=%e; E=%e\n", V1, E);
    exit(1);
  }
  printf("]");
  if (nl) printf("\n");
}

// check if relative difference between two (non-zero) values is less then E
void check(const char *name, double V1, double V2, double E, bool nl){
  printf(" [%s", name);
  if (fabs(V1-V2)/(fabs(V1)+fabs(V2)) > E){
    printf("\n FAIL: V1=%e; V2=%e, E=%e\n", V1, V2, E);
    exit(1);
  }
  printf("]");
  if (nl) printf("\n");
}

void check3z(const char *name, double V1[DIM], double E, bool nl){
  int i;
  printf(" [%s: ", name);
  FOR(i) checkz(xyz[i], V1[i], E, false);
  printf("]");
  if (nl) printf("\n");
}

void check3(const char *name, double V1[DIM], double V2[DIM], double E, bool nl){
  int i;
  printf(" [%s: ", name);
  FOR(i) check(xyz[i], V1[i], V2[i], E, false);
  printf("]");
  if (nl) printf("\n");
}
}