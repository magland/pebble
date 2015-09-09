#include "pebble_helpers.h"
#include <stdio.h>
#include <math.h>
#include "do_pca.h"
#include <QTime>

void split_into_channels(const QString &inpath,const QString &outpath) {
	QTime timer; timer.start();
	printf("split_into_channels... ");
	Mda X;
	X.read(inpath.toLatin1().data());
	double *Xptr=X.dataPtr();
	int M=X.N1();
	int N=X.N2();
	for (int ch=0; ch<M; ch++) {
		printf("channel %d/%d\n",ch,M);
		Mda AA; AA.allocate(1,N);
		double *AAptr=AA.dataPtr();
		//for (int j=0; j<N; j++) AA.setValue1(X.value(ch,j),j);
		int k=ch;
		for (int j=0; j<N; j++) {
			AAptr[j]=Xptr[k];
			k+=M;
		}
		QString path0=QString("%1/%2.mda").arg(outpath).arg(ch);
		AA.write(path0.toLatin1().data());
	}
	int elapsed=timer.elapsed();
	printf("Elapsed (ms): %d; %g MB/sec\n",elapsed,M*N*sizeof(double)*1e-6/(elapsed*1e-3));
}

Mda compute_adjacency_matrix(const Mda &locations) {
	Mda ret;
	int M=locations.N1();
	ret.allocate(M,M);
	for (int j=0; j<M; j++)
		for (int i=0; i<M; i++) {
			double dx=locations.value(i,0)-locations.value(j,0);
			double dy=locations.value(i,1)-locations.value(j,1);
			double dist=sqrt(dx*dx+dy*dy);
			if (dist<100) ret.setValue(1,i,j);
		}
	return ret;
}

QVector<int> find_patch_indices(const Mda &adjacency_matrix,int channel_num) {
	QVector<int> ret;
	int M=adjacency_matrix.N1();
	for (int i=0; i<M; i++) {
		if (adjacency_matrix.value(channel_num,i))
			ret << i;
	}
	return ret;
}

Mda extract_clips(Mda &X,int clip_size,QVector<int> clip_times) {
	double *Xptr=X.dataPtr();
	int num_clips=clip_times.count();
	int num_channels=X.N1();
	Mda clips; clips.allocate(num_channels,clip_size,num_clips);
	double *Cptr=clips.dataPtr();
	int jj=0;
	for (int i=0; i<num_clips; i++) {
		int kk=(clip_times[i]+(-clip_size/2))*num_channels;
		for (int t=0; t<clip_size; t++) {
			for (int ch=0; ch<num_channels; ch++) {
				Cptr[jj]=Xptr[kk];
				kk++;
				jj++;
			}
		}
	}
	return clips;
}

int find_max(const QVector<int> &X) {
	int ret=X.value(0);
	for (int i=0; i<X.count(); i++) {
		ret=qMax(ret,X.value(i));
	}
	return ret;
}

Mda compute_features_from_clips(Mda &clips,int num_features) {
	int M=clips.N1();
	int T=clips.N2();
	int N=clips.N3();
	int MT=M*T;
	Mda features; features.allocate(num_features,N);
	do_pca(MT,N,num_features,features.dataPtr(),clips.dataPtr());
	return features;
}

QVector<int> remove_two_closest_to_zero(Mda &features,const QVector<int> labels) {
	int M=features.N1();
	int N=features.N2();
	int num_labels=find_max(labels);
	QVector<double> minvals(num_labels+1);
	for (int n=0; n<N; n++) {
		double dist=0;
		for (int m=0; m<M; m++) {
			dist+=features.value(m,n)*features.value(m,n);
		}
		int label=labels[n];
		if ((minvals[label]==0)||(dist<minvals[label])) {
			minvals[label]=dist;
		}
	}
	int i1=-1; double v1=0;
	int i2=-1; double v2=0;
	for (int j=1; j<=num_labels; j++) {
		double val=minvals[j];
		if (i1<0) {i1=j; v1=val;}
		else if (i2<0) {i2=j; v2=val;}
		else {
			if (val<v1) {
				if (v1<v2) {
					i2=i1; v2=v1;
				}
				i1=j; v1=val;
			}
			else if (val<v2) {
				i2=j; v2=val;
			}
		}
	}
	int mapping[num_labels+1];
	int k=1;
	for (int j=1; j<=num_labels; j++) {
		if ((i1==j)||(i2==j)) {
			mapping[j]=0;
		}
		else {
			mapping[j]=k; k++;
		}
	}
	QVector<int> ret;
	for (int i=0; i<N; i++) ret << mapping[labels[i]];
	return ret;
}
