#ifndef UTIL_H
#define UTIL_H

#include <stdlib.h>

#define Alloc2D(Variable,Type,n1,n2) { Variable=(Type **)malloc(sizeof(Type *)*n1); Variable[0]=(Type *)malloc(sizeof(Type)*n1*n2); for(unsigned int A2D=1;A2D<n1;A2D++) Variable[A2D]=Variable[A2D-1]+n2; }

#define Free2D(Variable) { delete Variable[0]; delete Variable;}

#endif
