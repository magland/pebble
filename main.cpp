//#include <QCoreApplication>
#include "mda.h"
#include <QString>
#include <QDebug>
#include <math.h>
#include <QTime>
#include <QDir>

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
	for (int t=0; t<N; t++) {
		double maxval=0;
		int maxind=t;
		for (int dt=-radius; dt<=radius; dt++) {
			int t0=t+dt;
			if ((0<=t0)&&(t0<N)) {
				if (qAbs(X[t0])>maxval) {
					maxval=qAbs(X[t0]);
					maxind=t0;
				}
			}
		}
		if (maxind==t) ret << t;

	}
	return ret;
}

QList<int> find_critical_times_fast(int N,double *X,int radius) {
	QList<int> ret;
	if (N<=0) return ret;
	int last_best_ind=0; //the index of last maximum within the radius
	double last_best_val=qAbs(X[0]); //the value of the last maximum within the radius
	bool is_pseudo=false; //it true, this means that the above is the last maximum within the radius, but it cannot be counted as a critical point
	for (int i=1; i<N; i++) {
		if (qAbs(X[i])>last_best_val) {
			//we have a new maximum!
			if ((!is_pseudo)&&(i-last_best_ind>radius)) { //if the old one is not a pseudo, and it is now outside the range, add that one to our list
				ret << last_best_ind;
			}
			last_best_ind=i;
			last_best_val=qAbs(X[i]);
			is_pseudo=false; //not pseudo, because it was bigger than the previous biggest within the radius
		}
		else {
			//not a new maximum
			if (i-last_best_ind>radius) { //the old maximum is now outside the radius
				if (!is_pseudo) { //add it to our list if it wasn't a pseudo
					ret << last_best_ind;
				}
				double maxval=0; //we need to find the new best
				int maxind=i;
				for (int dt=-radius; dt<=0; dt++) {
					int t0=i+dt;
					if (0<=t0) {
						if (qAbs(X[t0])>maxval) {
							maxval=qAbs(X[t0]);
							maxind=t0;
						}
					}
				}
				last_best_ind=maxind;
				last_best_val=maxval;
				if (last_best_ind==i) is_pseudo=false; //it's not a pseudo if it is the current index
				else is_pseudo=true; //otherwise it is a pseudo because it was within the radius of the previous max
			}
		}
	}
	if (!is_pseudo) { //add the final index if not a pseudo
		ret << last_best_ind;
	}
	return ret;
}

int main(int argc, char *argv[])
{
	QString channels_path="/home/magland/data/EJ/channels";
	if (0) {
		if (!QFile::exists(channels_path)) QDir(QFileInfo(channels_path).path()).mkdir(QFileInfo(channels_path).fileName());
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

		int N=DD.totalSize();
		{
			QTime timer;
			timer.start();
			QList<int> critical_times=find_critical_times(N,DD.dataPtr(),50);
			qDebug() << critical_times.mid(0,20);
			qDebug() << critical_times[critical_times.count()-1];
			qDebug() << QString("Found %1 of %2 (%3\%)").arg(critical_times.count()).arg(N).arg(critical_times.count()*100.0/N);
			qDebug() << "Elapsed: " << timer.elapsed();
		}

		{
			QTime timer;
			timer.start();
			QList<int> critical_times=find_critical_times_fast(N,DD.dataPtr(),50);
			qDebug() << critical_times.mid(0,20);
			qDebug() << critical_times[critical_times.count()-1];
			qDebug() << QString("Found %1 of %2 (%3\%)").arg(critical_times.count()).arg(N).arg(critical_times.count()*100.0/N);
			qDebug() << "Elapsed: " << timer.elapsed();
		}
	}

	return 0;
}


