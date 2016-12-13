#ifndef UTIL_H
#define UTIL_H

#include <stdlib.h>

#define Alloc2D(Variable,Type,n1,n2) { Variable=(Type **)malloc(sizeof(Type *)*n1); Variable[0]=(Type *)malloc(sizeof(Type)*n1*n2); for(unsigned int A2D=1;A2D<n1;A2D++) Variable[A2D]=Variable[A2D-1]+n2; }
#define Free2D(Variable) { delete Variable[0]; delete Variable;}

#ifndef MINFLOAT
#define MINFLOAT        ((float)1.17549435e-38)
#endif

#ifndef MAXFLOAT
#define MAXFLOAT        ((float)3.40282347e+38)
#endif

#define MINDOUBLE       2.2250738585072014e-308
#define MAXDOUBLE       1.7976931348623157e+308
#define M_PI            3.14159265358979323846



#endif
