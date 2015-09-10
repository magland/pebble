#include <QCoreApplication>
#include "mda.h"
#include <QString>
#include <QDebug>
#include <math.h>
#include <QTime>
#include <QDir>
#include "parse_command_line_params.h"
#include "pebble_helpers.h"
#include "find_critical_times.h"
#include "isosplit.h"

void multi_channel_convolution(int M,int N,int T,double *out,double *in,double *kernel) {
	for (int ii=0; ii<N; ii++) out[ii]=0;
	int ii_out=0;
	int ii_in=0;
	for (int n=0; n<N-T; n++) {
		int ii_kernel=0;
		for (int dt=0; dt<T; dt++) {
			for (int m=0; m<M; m++) {
				out[ii_out]+=in[ii_in]*kernel[ii_kernel];
				ii_in++;
				ii_kernel++;
			}
		}
		ii_in += -M*T + M;
		ii_out++;
	}
}

bool lock_file(const QString &path,int timeout) {
	QString path0=path+".lock";
	QTime timer; timer.start();
	while (timer.elapsed()<timeout) {
		if (!QFile::exists(path0)) {
			FILE *inf=fopen(path0.toLatin1(),"w");
			if (inf) {
				fclose(inf);
				return true;
			}
		}
	}
	return false;
}
bool unlock_file(const QString &path) {
	QString path0=path+".lock";
	return QFile::remove(path0);
}

int main(int argc, char *argv[])
{
    QCoreApplication app(argc,argv);

    qsrand(time(NULL));

    QStringList required_params;
	QStringList optional_params; optional_params << "channel" << "neuron";
    CLParams CLP=parse_command_line_params(argc,argv,required_params,optional_params);
    QString command=CLP.unnamed_parameters.value(0);

    printf("starting pebble...\n");

	int clip_size=80;
	int find_critical_times_radius=150;
	int num_features=3;

    //QString channels_path="/home/magland/data/EJ/channels";
	//QString channels_path="/home/magland/data/EJ/channels";
	//QString channels_path="/dev/shm/channels";
	QString channels_path="/mnt/xfs1/home/magland/data/EJ/channels";
    //QString channels_path="/dev/shm/channels";
    //QString testdata_path="/home/magland/dev/pebble/testdata";
	//QString testdata_path="/home/magland/dev/pebble/testdata";
	//QString testdata_path="/dev/shm/testdata";
	QString testdata_path="/mnt/xfs1/home/magland/dev/pebble/testdata";
	QString locations_path=app.applicationDirPath()+"/../testdata/locations.mda";
	QTime timer;

	if (!QFile::exists(testdata_path)) QDir(QFileInfo(testdata_path).path()).mkdir(QFileInfo(testdata_path).fileName());

	printf("Reading locations file...\n");
	Mda locations; locations.read(locations_path.toLatin1().data());
	qDebug() << locations.N1() << locations.N2();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Split into channels
    if (command=="split_into_channels") {
        printf("Splitting into channels...\n");
		if (!QFile::exists(channels_path)) QDir(QFileInfo(channels_path).path()).mkdir(QFileInfo(channels_path).fileName());
		split_into_channels(
					"/home/magland/data/EJ/Spikes_all_channels_filtered.mda",
					channels_path
					);
		// compute adjacency matrix
		printf("computing adjacency matrix...\n");
		Mda AM0=compute_adjacency_matrix(locations);
		qDebug() << QString("%1 x %2").arg(AM0.N1()).arg(AM0.N2());
		AM0.write(testdata_path+"/AM.mda");

        return 0;
	}

	Mda AM; AM.read(testdata_path+"/AM.mda");

    if (command=="stage1") {
		if (!CLP.named_parameters.contains("channel")) {
			printf("Missing parameter: channel\n");
			return -1;
		}
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
				/*if (!lock_file(path0,1000)) { //not sure if necessary
					qCritical() << "Unable to lock file: " << path0;
					return -1;
				}*/
                Mda channel_data; channel_data.read(path0); //1xN
				/*if (!unlock_file(path0)) { //not sure if necessary
					qCritical() << "Unable to unlock file: " << path0;
					return -1;
				}*/
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
			qDebug() << QString("%1 x %2 x %3").arg(patch_data.N1()).arg(patch_data.N2()).arg(patch_data.N3());
			int elapsed=timer.elapsed();
			printf("Elapsed (ms): %d; %g MB/sec\n",elapsed,N*patch_size*1e-6/(elapsed*1e-3));

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
				QString fname0=QString("%1/times/%2/times-%2-%3.mda").arg(testdata_path).arg(channel_num).arg(cc);
				printf("Writing %s\n",fname0.toLatin1().data());
				times_mda.write(fname0);
            }
            printf("Elapsed (ms): %d\n",timer.elapsed());
        }
        return 0;
    }

	if (command=="stage2") {
		//Take all of the times and convert them into spike shapes by averaging over the data for this channel
		if (!CLP.named_parameters.contains("channel")) {
			printf("Missing parameter: channel\n");
			return -1;
		}

		int channel=CLP.named_parameters["channel"].toInt();

		printf(":::::::::: STAGE2 for CHANNEL %d\n",channel);
		///////////////////////////////////////////////////////////////////////////////////////////////////////////
		///// Read channel data
		Mda channel_data;
		{
			printf("Reading channel data... "); timer.start();
			QString path0=QString("%1/%2.mda").arg(channels_path).arg(channel);
			channel_data.read(path0); //1xN
			printf("Elapsed (ms): %d\n",timer.elapsed());
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////
		///// Creating spike shapes
		Mda spike_shapes;
		Mda channel_ids_mda;
		QStringList paths;
		{
			printf("Creating spike shapes... "); timer.start();

			QList<int> channel_ids;
			QString tpath=QString("%1/times").arg(testdata_path);
			QStringList list1=QDir(tpath).entryList(QStringList("*"),QDir::Dirs|QDir::NoDotAndDotDot,QDir::Name);
			foreach (QString str1,list1) {
				QStringList list2=QDir(tpath+"/"+str1).entryList(QStringList("*.mda"),QDir::Files,QDir::Name);
				foreach (QString str2,list2) {
					paths << QString("%1/times/%2/%3").arg(testdata_path).arg(str1).arg(str2);
					int channel0=str1.toInt();
					channel_ids << channel0;
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
                for (int tt=0; tt<clip_size; tt++) avg[tt]/=clip_times.count();
				for (int tt=0; tt<clip_size; tt++) {
					spike_shapes.setValue(avg[tt],0,tt,ii);
				}
			}

			printf("Created array (%d x %d x %d)... ",spike_shapes.N1(),spike_shapes.N2(),spike_shapes.N3());
			printf("Elapsed (ms): %d\n",timer.elapsed());
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////
		///// Writing spike shapes
		///// also write channel_ids.mda if channel==0
		{
			printf("Writing spike shapes... "); timer.start();
			if (!QFile::exists(testdata_path+"/shapes")) QDir(testdata_path).mkdir("shapes");
			spike_shapes.write(QString("%1/shapes/shapes-%2.mda").arg(testdata_path).arg(channel));
			if (channel==0) channel_ids_mda.write(QString("%1/channel_ids.mda").arg(testdata_path));
			printf("Elapsed (ms): %d\n",timer.elapsed());
		}

		return 0;
	}

	if (command=="stage3") {
		// Combine the channels into a single combined file: shapes.mda
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

	if (command=="compute_scores") {
		//Compute the scores for a particular neuron
		if (!CLP.named_parameters.contains("neuron")) {
			printf("Missing parameter: neuron\n");
			return -1;
		}

		int neuron=CLP.named_parameters["neuron"].toInt();
		Mda AM; AM.read(testdata_path+"/AM.mda");
		Mda channel_ids_mda; channel_ids_mda.read(QString("%1/channel_ids.mda").arg(testdata_path));
		int channel_num=(int)channel_ids_mda.value(0,neuron);
		QVector<int> patch_indices=find_patch_indices(AM,channel_num);
		printf(":::::::::: COMPUTE SCORES for neuron %d at channel %d\n",neuron,channel_num);

		///////////////////////////////////////////////////////////////////////////////////////////////////////////
		//// Reading shapes file
		Mda shape;
		int T=0,M=0;
		{
			printf("Reading shapes file... "); timer.start();
			Mda shapes;
			shapes.read(QString("%1/shapes.mda").arg(testdata_path));
			M=shapes.N1();
			T=shapes.N2();
			shape.allocate(M,T);
			for (int t=0; t<T; t++) {
				for (int m=0; m<M; m++) {
                    double val=shapes.value(m,t,neuron);
					shape.setValue(val,m,t);
				}
			}
			printf("Elapsed (ms): %d\n",timer.elapsed());
		}

        shape.write(testdata_path+"/shape.mda"); //debug

		///////////////////////////////////////////////////////////////////////////////////////////////////////////
		//// Computing patch scores
		Mda patch_scores;
		{
			printf("Computing patch scores... "); timer.start();

			Mda patch_data;
            Mda patch_shape;
			int N=0;
			int MM=patch_indices.count();
			QTime read_timer; read_timer.start();
            for (int mm=0; mm<MM; mm++) {
                int m=patch_indices[mm];
				Mda X; X.read(QString("%1/%2.mda").arg(channels_path).arg(m));
                if (mm==0) {
					N=X.N2();
					patch_data.allocate(MM,N);
                    patch_shape.allocate(MM,T);
				}
				for (int n=0; n<N; n++) {
                    patch_data.setValue(X.value(0,n),mm,n);
				}
                for (int t=0; t<T; t++) {
                    patch_shape.setValue(shape.value(m,t),mm,t);
                }
			}
			printf("Time for reading patch data (ms): %d\n",read_timer.elapsed());

            patch_data.write(testdata_path+"/patch_data.mda"); //debug
            patch_shape.write(testdata_path+"/patch_shape.mda"); //debug

			double *kernel=(double *)malloc(sizeof(double)*MM*T);
			double *data=(double *)malloc(sizeof(double)*MM*N);
            Mda convolution; convolution.allocate(1,N);
            double *convolutionPtr=convolution.dataPtr();
			int ct=0;
			for (int t=0; t<T; t++) {
				for (int mm=0; mm<MM; mm++) {
                    kernel[ct]=patch_shape.value(mm,t); ct++;
				}
			}
			ct=0;

			QTime timerA; timerA.start();
			for (int n=0; n<N; n++) {
				for (int mm=0; mm<MM; mm++) {
					data[ct]=patch_data.value(mm,n); ct++;
				}
			}
			printf("Time for preparing data (ms): %d\n",timerA.elapsed());

			QTime MCC_timer; MCC_timer.start();
            multi_channel_convolution(MM,N,T,convolutionPtr,data,kernel);
			printf("Time for convolution (ms): %d\n",MCC_timer.elapsed());

			QTime timerB; timerB.start();
			double kernel_sumsqr=0;
			ct=0;
			for (int t=0; t<T; t++) {
				for (int mm=0; mm<MM; mm++) {
					kernel_sumsqr+=kernel[ct]*kernel[ct];
					ct++;
				}
			}
			//sum(a_i)^2 - sum(a_i-k_i)^2 = Saa - Saa - Skk + 2Sak = 2Sak - Skk

			patch_scores.allocate(1,N);
			double *ptr=patch_scores.dataPtr();
			for (int n=0; n<N-T; n++) {
                ptr[n]=2*convolutionPtr[n]-kernel_sumsqr;
			}
			printf("Time for setting scores (ms): %d\n",timerB.elapsed());

			free(kernel);
			free(data);
			printf("Elapsed (ms): %d\n",timer.elapsed());

            convolution.write(testdata_path+"/convolution.mda"); //debug
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////////
		//// Computing full scores
		Mda full_scores;
		{
			printf("Computing full scores... "); timer.start();

			Mda full_data;
			int N=0;
			QList<int> full_indices; for (int m=0; m<M; m++) full_indices << m;
			int MM=full_indices.count();
			QTime read_timer; read_timer.start();
			for (int ii=0; ii<MM; ii++) {
				int m=full_indices[ii];
				Mda X; X.read(QString("%1/%2.mda").arg(channels_path).arg(m));
				if (ii==0) {
					N=X.N2();
					full_data.allocate(MM,N);
				}
				double *full_data_ptr=full_data.dataPtr();
				double *Xptr=X.dataPtr();
				int jj=ii;
				int kk=0;
				for (int n=0; n<N; n++) {
					full_data_ptr[jj]=Xptr[kk];
					jj+=MM;
					kk++;
				}
			}
			printf("Time for reading full data (ms): %d\n",read_timer.elapsed());

			double *data=(double *)malloc(sizeof(double)*MM*N);
			double *kernel=(double *)malloc(sizeof(double)*MM*T);
			int ct=0;
			for (int t=0; t<T; t++) {
				for (int mm=0; mm<MM; mm++) {
					int m=full_indices[mm];
					kernel[ct]=shape.value(m,t); ct++;
				}
			}

			QTime timerA; timerA.start();
			ct=0;
			for (int n=0; n<N; n++) {
				for (int mm=0; mm<MM; mm++) {
//                    data[ct]=full_data_ptr[ct]; ct++; //finish!!!
				}
			}
			printf("Time for preparing full data (ms): %d\n",timerA.elapsed());

			QTime timerB; timerB.start();
			double kernel_sumsqr=0;
			ct=0;
			for (int t=0; t<T; t++) {
				for (int mm=0; mm<MM; mm++) {
					kernel_sumsqr+=kernel[ct]*kernel[ct];
					ct++;
				}
			}

			full_scores.allocate(1,N);
			double *ptr=full_scores.dataPtr();
			double *ptr2=patch_scores.dataPtr();
			int ii_data=0;
			int num_positive=0;
			for (int n=0; n<N-T; n++) {
				if (ptr2[n]>0) {
					num_positive++;
					double val=0;
					int jj=0;
					for (int t=0; t<T; t++) {
						for (int mm=0; mm<MM; mm++) {
							val+=data[ii_data+jj]*kernel[jj]; jj++;
						}
					}
					ptr[n]=2*val-kernel_sumsqr;
				}
				else {
                    //ptr[n]=0;
                    ptr[n]=ptr2[n];
				}
				ii_data+=MM;
			}
			printf("Time for setting full scores (ms): %d, # positive: %d\n",timerB.elapsed(),num_positive);

            full_scores.write(testdata_path+"/full_scores.mda");

			free(kernel);
			free(data);
			printf("Elapsed (ms): %d\n",timer.elapsed());
		}


		return 0;
	}



    qWarning() << "Unknown command: "+command;

	return 0;
}
