#ifndef CHECK_H
#define CHECK_H
extern "C" {

// check if absolute value is less then E
void checkz(const char *name, double V1, double E, bool nl);

// check if relative difference between two (non-zero) values is less then E
void check(const char *name, double V1, double V2, double E, bool nl);

// check if absolute value of all vector elements is less then E
void check3z(const char *name, double V1[DIM], double E, bool nl);

// check if relative difference of all vector elements is less then E
void check3(const char *name, double V1[DIM], double V2[DIM], double E, bool nl);

}
#endif