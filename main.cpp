//#include <QCoreApplication>
#include "mda.h"
#include <QString>
#include <QDebug>
#include <math.h>
#include <QTime>
#include <QDir>
#include "find_critical_times.h"

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

/*
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
*/

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
		int kk=clip_times[i]+(-clip_size/2);
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

int main(int argc, char *argv[])
{
	QString channels_path="/home/magland/data/EJ/channels";
	QTime timer;

	// Split into channels
	if (0) {
		if (!QFile::exists(channels_path)) QDir(QFileInfo(channels_path).path()).mkdir(QFileInfo(channels_path).fileName());
		split_into_channels(
					"/home/magland/data/EJ/Spikes_all_channels_filtered.mda",
					channels_path
					);
	}

	int clip_size=50;
	int channel_num=0;
	int find_critical_times_radius=150;

	// compute adjacency matrix
	QString locations_path="/home/magland/dev/pebble/testdata/locations.mda";
	Mda locations; locations.read(locations_path.toLatin1().data());
	Mda AM=compute_adjacency_matrix(locations);
	AM.write("testdata/AM.mda");

	// Read patch data
	printf("Reading patch data... "); timer.start();
	QVector<int> patch_indices=find_patch_indices(AM,channel_num);
	int patch_size=patch_indices.count();
	Mda patch_data; int N=0;
	Mda channel0_data;
	for (int p=0; p<patch_size; p++) {
		QString path0=QString("%1/%2.mda").arg(channels_path).arg(patch_indices[p]);
		Mda channel_data; channel_data.read(path0.toLatin1().data()); //1xN
		if (p==0) {
			N=channel_data.N2();
			patch_data.allocate(patch_size,N);
			channel0_data=channel_data;
		}
		for (int n=0; n<N; n++) patch_data.setValue(channel_data.value(0,n),p,n);
	}
	printf("Elapsed (ms): %d\n",timer.elapsed());

	//Find critical times
	printf("Finding critical_times... "); timer.start();
	QVector<int> critical_times=find_critical_times(channel0_data,find_critical_times_radius);
	printf("Elapsed (ms): %d\n",timer.elapsed());

	//Extracting clips
	printf("Extracting clips... "); timer.start();
	Mda clips=extract_clips(patch_data,clip_size,critical_times); //patch_size x clip_size
	printf("Elapsed (ms): %d\n",timer.elapsed());

	//here is where we would upsample

	//Clustering
	printf("Clustering... "); timer.start();
	printf("Elapsed (ms): %d\n",timer.elapsed());

	return 0;
}


