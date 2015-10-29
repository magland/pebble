#ifndef PEBBLE_HELPERS_H
#define PEBBLE_HELPERS_H

#include <QString>
#include "mda.h"
#include <QVector>

void split_into_channels(const QString &inpath,const QString &outpath);
Mda compute_adjacency_matrix(const Mda &locations,double radius);
QVector<int> find_patch_indices(const Mda &adjacency_matrix,int channel_num);
Mda extract_clips(Mda &X,int clip_size,QVector<int> clip_times);
int find_max(const QVector<int> &X);
Mda compute_features_from_clips(Mda &clips,int num_features);
QVector<int> remove_two_closest_to_zero(Mda &features,const QVector<int> labels);

#endif // PEBBLE_HELPERS_H

