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
#include "parse_command_line_params.h"

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

int main(int argc, char *argv[])
{
    QCoreApplication app(argc,argv);

    qsrand(time(NULL));

    QStringList required_params;
    QStringList optional_params; optional_params << "channel";
    CLParams CLP=parse_command_line_params(argc,argv,required_params,optional_params);
    QString command=CLP.unnamed_parameters.value(0);

    printf("starting pebble...\n");

    //QString channels_path="/home/magland/data/EJ/channels";
    QString channels_path="/mnt/xfs1/home/magland/data/EJ/channels";
    //QString channels_path="/dev/shm/channels";
    //QString testdata_path="/home/magland/dev/pebble/testdata";
    QString testdata_path="/mnt/xfs1/home/magland/dev/pebble/testdata";
	QTime timer;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Split into channels
    if (command=="split_into_channels") {
        printf("Splitting into channels...\n");
		if (!QFile::exists(channels_path)) QDir(QFileInfo(channels_path).path()).mkdir(QFileInfo(channels_path).fileName());
		split_into_channels(
                    "/mnt/xfs1/home/magland/data/EJ/Spikes_all_channels_filtered.mda",
					channels_path
					);
        return 0;
	}

    int clip_size=80;
	int find_critical_times_radius=150;
    int num_features=3;

    // compute adjacency matrix
    printf("computing adjacency matrix...\n");
    QString locations_path=testdata_path+"/locations.mda";
    Mda locations; locations.read(locations_path.toLatin1().data());
    Mda AM=compute_adjacency_matrix(locations);
    AM.write(testdata_path+"/AM.mda");

    if (command=="stage3") {
        printf(":::::::::: STAGE3\n");
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///// Combine channels
        Mda XX;
        {
            printf("Combining channels... "); timer.start();
            QString spath=QString("%1/shapes").arg(testdata_path);
            QStringList list1=QDir(spath).entryList(QStringList("shapes-*.mda"),QDir::Files,QDir::Name);
            int num_channels=list1.count();
            int num_neurons=0;
            for (int ch=0; ch<num_channels; ch++) {
                QString str=QString("shapes-%1.mda").arg(ch);
                Mda A; A.read(spath+"/"+str);
                if (ch==0) {
                    num_neurons=A.N3();
                    XX.allocate(num_channels,clip_size,num_neurons);
                }
                for (int nn=0; nn<num_neurons; nn++) {
                    for (int tt=0; tt<clip_size; tt++) {
                        XX.setValue(A.value(0,tt,nn),ch,tt,nn);
                    }
                }
            }
            printf("Elapsed (ms): %d\n",timer.elapsed());
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///// Writing combined file
        {
            printf("Writing combined file... "); timer.start();
            XX.write(QString("%1/shapes.mda").arg(testdata_path));
            printf("Elapsed (ms): %d\n",timer.elapsed());
        }

        return 0;
    }

    if (command=="stage2") {
        //Take all of the times and convert them into spike shapes by averaging over the data for this channel
        int channel=CLP.named_parameters["channel"].toInt();

        printf(":::::::::: STAGE2 for CHANNEL %d\n",channel);
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///// Read channel data
        Mda channel_data;
        {
            printf("Reading channel data... "); timer.start();
            QString path0=QString("%1/%2.mda").arg(channels_path).arg(channel);
            channel_data.read(path0); //1xN
            int N=channel_data.N2();
            printf("Elapsed (ms): %d\n",timer.elapsed());
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///// Creating spike shapes
        Mda spike_shapes;
        Mda channel_ids_mda;
        {
            printf("Creating spike shapes... "); timer.start();

            QStringList paths;
            QList<int> channel_ids;
            QString tpath=QString("%1/times").arg(testdata_path);
            QStringList list1=QDir(tpath).entryList(QStringList("*"),QDir::Dirs|QDir::NoDotAndDotDot,QDir::Name);
            foreach (QString str1,list1) {
                QStringList list2=QDir(tpath+"/"+str1).entryList(QStringList("*.mda"),QDir::Files,QDir::Name);
                foreach (QString str2,list2) {
                    paths << QString("%1/times/%2/%3").arg(testdata_path).arg(str1).arg(str2);
                    channel_ids << str1.toInt();
                }
            }
            spike_shapes.allocate(1,clip_size,paths.count());
            channel_ids_mda.allocate(1,paths.count());
            for (int ii=0; ii<paths.count(); ii++) {
                channel_ids_mda.setValue(channel_ids[ii],0,ii);
                Mda A; A.read(paths[ii]);
                QVector<int> clip_times; for (int i=0; i<A.N2(); i++) clip_times << (int)A.value(0,i);
                Mda clips=extract_clips(channel_data,clip_size,clip_times);
                QVector<double> avg(clip_size); for (int jj=0; jj<clip_size; jj++) avg[jj]=0;
                for (int jj=0; jj<clip_times.count(); jj++) {
                    for (int tt=0; tt<clip_size; tt++) {
                        avg[tt]+=clips.value(0,tt,jj);
                    }
                }
                for (int tt=0; tt<clip_size; tt++) avg[tt]/=clip_size;
                for (int tt=0; tt<clip_size; tt++) {
                    spike_shapes.setValue(avg[tt],0,tt,ii);
                }
            }
            printf("Created array (%d x %d x %d)... ",spike_shapes.N1(),spike_shapes.N2(),spike_shapes.N3());
            printf("Elapsed (ms): %d\n",timer.elapsed());
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///// Writing spike shapes
        {
            printf("Writing spike shapes... "); timer.start();
            if (!QFile::exists(testdata_path+"/shapes")) QDir(testdata_path).mkdir("shapes");
            spike_shapes.write(QString("%1/shapes/shapes-%2.mda").arg(testdata_path).arg(channel));
            if (channel==0) channel_ids_mda.write(QString("%1/channel_ids.mda").arg(testdata_path));
            printf("Elapsed (ms): %d\n",timer.elapsed());
        }

        return 0;
    }

    if (command=="stage1") {

        int channel=CLP.named_parameters["channel"].toInt();
        for (int channel_num=channel; channel_num<=channel; channel_num++) {
            printf(":::::::::: STAGE 1 for CHANNEL %d\n",channel_num);
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///// Read patch data
            printf("Reading patch data... "); timer.start();
            QVector<int> patch_indices=find_patch_indices(AM,channel_num);
            int patch_size=patch_indices.count();
            Mda patch_data; int N=0;
            Mda channel0_data;
            for (int p=0; p<patch_size; p++) {
                QString path0=QString("%1/%2.mda").arg(channels_path).arg(patch_indices[p]);
                Mda channel_data; channel_data.read(path0); //1xN
                if (p==0) {
                    N=channel_data.N2();
                    patch_data.allocate(patch_size,N);
                    channel0_data=channel_data;
                    printf("allocated %g MB... ",patch_size*N*8*1.0/1000000);
                }
                double *patch_data_ptr=patch_data.dataPtr();
                double *channel_data_ptr=channel_data.dataPtr();
                int iii=p;
                for (int n=0; n<N; n++) {
                    patch_data_ptr[iii]=channel_data_ptr[n];
                    iii+=patch_size;
                }
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
            printf("Elapsed (ms): %d\n",timer.elapsed());

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////
            //here is where we would upsample and time align?

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////
            //Features
            printf("Computing features... "); timer.start();
            Mda features=compute_features_from_clips(clips,num_features); //num_features x num_clips
            printf("Elapsed (ms): %d\n",timer.elapsed());

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////
            //Clustering
            printf("Clustering... "); timer.start();
            QVector<int> labels=isosplit(features);
            labels=remove_two_closest_to_zero(features,labels);
            int num_clusters=find_max(labels);
            printf("Found %d clusters. ",num_clusters);
            printf("Elapsed (ms): %d\n",timer.elapsed());

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////
            //Spike times
            printf("Writing spike times... "); timer.start();
            for (int cc=1; cc<=num_clusters; cc++) {
                QVector<int> indices; for (int i=0; i<labels.count(); i++) if (labels[i]==cc) indices << i;
                Mda times_mda; times_mda.allocate(1,indices.count());
                for (int i=0; i<indices.count(); i++) times_mda.setValue(critical_times[indices[i]],0,i);
                if (!QFile::exists(testdata_path+"/times")) QDir(testdata_path).mkdir("times");
                if (!QFile::exists(testdata_path+"/times/"+QString("%1").arg(channel_num))) QDir(testdata_path+"/times").mkdir(QString("%1").arg(channel_num));
                times_mda.write(QString("%1/times/%2/times-%2-%3.mda").arg(testdata_path).arg(channel_num).arg(cc));
            }
            printf("Elapsed (ms): %d\n",timer.elapsed());
        }
        return 0;
    }

    qWarning() << "Unknown command: "+command;

	return 0;
}
