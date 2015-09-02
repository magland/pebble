#ifndef find_critical_times_h
#define find_critical_times_h

#include <stdlib.h>
#include <QVector>
#include "mda.h"

QVector<int> find_critical_times(Mda &X,int radius);
QVector<int> find_critical_times_brute_force(Mda &X,int radius);

#endif
