#include <QCoreApplication>
#include "mda.h"
#include <QString>
#include <QDebug>
#include <math.h>
#include <QTime>
#include <QDir>
#include "find_critical_times.h"
#include "isosplit.h"
#include "do_pca.h"

void split_into_channels(const QString &inpath,const QString &outpath) {
	Mda X;
	X.read(inpath.toLatin1().data());
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

int main(int argc, char *argv[])
{
    QCoreApplication app(argc,argv);

    qsrand(time(NULL));

    printf("starting pebble...\n");

	QString channels_path="/home/magland/data/EJ/channels";
    QString testdata_path="/home/magland/dev/pebble/testdata";
	QTime timer;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Split into channels
    if (0) {
		if (!QFile::exists(channels_path)) QDir(QFileInfo(channels_path).path()).mkdir(QFileInfo(channels_path).fileName());
		split_into_channels(
					"/home/magland/data/EJ/Spikes_all_channels_filtered.mda",
					channels_path
					);
	}

	int clip_size=50;

	int find_critical_times_radius=150;
    int num_features=3;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
	// compute adjacency matrix
    printf("computing adjacency matrix...\n");
    QString locations_path=testdata_path+"/locations.mda";
	Mda locations; locations.read(locations_path.toLatin1().data());
	Mda AM=compute_adjacency_matrix(locations);
    AM.write(testdata_path+"/AM.mda");


    bool write_intermediate=false;

    for (int channel_num=0; channel_num<20; channel_num++) {
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

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////
        //Find critical times
        printf("Finding critical_times... "); timer.start();
        QVector<int> critical_times=find_critical_times(channel0_data,find_critical_times_radius);
        int num_clips=critical_times.count();
        printf("found %d critical times.... ",num_clips);
        printf("Elapsed (ms): %d\n",timer.elapsed());

        //Extracting clips
        printf("Extracting clips... "); timer.start();
        Mda clips=extract_clips(patch_data,clip_size,critical_times); //patch_size x clip_size x num_clips
        if (write_intermediate) clips.write(testdata_path+"/clips.mda");
        printf("Elapsed (ms): %d\n",timer.elapsed());


        ///////////////////////////////////////////////////////////////////////////////////////////////////////////
        //here is where we would upsample and time align?

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////
        //Features
        printf("Computing features... "); timer.start();
        Mda features=compute_features_from_clips(clips,num_features); //num_features x num_clips
        if (write_intermediate) features.write(testdata_path+"/features.mda");
        printf("Elapsed (ms): %d\n",timer.elapsed());

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////
        //Clustering
        printf("Clustering... "); timer.start();
        QVector<int> labels=isosplit(features);
        Mda labels_mda; labels_mda.allocate(1,labels.count()); for (int ii=0; ii<labels.count(); ii++) labels_mda.setValue(labels[ii],0,ii);
        if (write_intermediate) labels_mda.write(testdata_path+"/labels.mda");
        printf("Elapsed (ms): %d\n",timer.elapsed());

        printf("Found %d clusters.\n",find_max(labels));
    }

	return 0;
}
