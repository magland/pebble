//#include <QCoreApplication>
#include "mda.h"
#include <QString>
#include <QDebug>
#include <math.h>
#include <QTime>

void split_into_channels(const QString &inpath,const QString &outpath) {
	Mda X;
	X.read(inpath.toLatin1().data());
	qDebug() << X.N1() << X.N2() << X.N3();
	int M=X.N1();
	int N=X.N2();
	for (int ch=0; ch<M; ch++) {
		printf("channel %d/%d\n",ch,M);
		Mda AA; AA.allocate(1,N);
		for (int j=0; j<N; j++) AA.setValue1(X.value(ch,j),j);
		QString path0=QString("%1/%2.mda").arg(outpath).arg(ch);
		AA.write(path0.toLatin1().data());
	}
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

QList<int> find_critical_times(int N,double *X,int radius) {
	QList<int> ret;
	double *localmax=(double *)malloc(sizeof(double)*N);
	for (int t=0; t<N; t++) {
		double maxval=0;
		for (int dt=-radius; dt<=radius; dt++) {
			int t0=t+dt;
			if ((0<=t0)&&(t0<N)) {
				if (qAbs(X[t0])>maxval) maxval=qAbs(X[t0]);
			}
		}
		localmax[t]=maxval;
	}
	for (int t=0; t<N; t++) {
		if (qAbs(X[t])==localmax[t]) ret << t;
	}
	free(localmax);
	return ret;
}

int main(int argc, char *argv[])
{
	QString channels_path="/home/magland/data/EJ/channels";
	if (0) {
		split_into_channels(
					"/home/magland/data/EJ/Spikes_all_channels_filtered.mda",
					channels_path
					);
	}

	QString locations_path="/home/magland/dev/pebble/testdata/locations.mda";
	Mda locations; locations.read(locations_path.toLatin1().data());
	Mda AM=compute_adjacency_matrix(locations);
	AM.write("testdata/AM.mda");

	printf("Reading channel 0\n");
	QString path0=QString("%1/%2.mda").arg(channels_path).arg(0);
	Mda DD; DD.read(path0.toLatin1().data());

	{
		printf("Finding critical times...\n");
		QTime timer;
		timer.start();
		int N=DD.totalSize();
		QList<int> critical_times=find_critical_times(N,DD.dataPtr(),50);
		qDebug() << QString("Found %1 of %2 (%3\%)").arg(critical_times.count()).arg(N).arg(critical_times.count()*100.0/N);
		qDebug() << "Elapsed: " << timer.elapsed();
	}

	return 0;
}


