#ifndef find_critical_times_h
#define find_critical_times_h

#include <stdlib.h>
#include <QVector>

QVector<int> find_critical_times(int N,const double *X,int radius);
QVector<int> find_critical_times_brute_force(int N,const double *X,int radius);

#endif