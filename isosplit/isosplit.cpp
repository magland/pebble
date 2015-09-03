#include "isosplit.h"
#include <QSet>
#include <QVector>
#include <QDebug>
#include <math.h>
#include "isosplit1d.h"

//choose K distinct (sorted) integers between 0 and N-1. If K>N then it will repeat the last integer a suitable number of times
QList<int> choose_random_indices(int N,int K) {;
	QList<int> ret;
	if (K>=N) {
		for (int i=0; i<N; i++) ret << i;
		while (ret.count()<K) ret << N-1;
		return ret;
	}
	QSet<int> theset;
	while (theset.count()<K) {
		int ind=(qrand()%N);
		theset.insert(ind);
	}
	ret=theset.toList();
	qSort(ret);
	return ret;
}

//do k-means with K clusters -- X is MxN representing N points in M-dimensional space. Returns a labels vector of size N.
QVector<int> do_kmeans(Mda &X,int K) {
	int M=X.N1();
	int N=X.N2();
	double *Xptr=X.dataPtr();
	Mda centroids_mda; centroids_mda.allocate(M,K); double *centroids=centroids_mda.dataPtr();
	QVector<int> labels; for (int i=0; i<N; i++) labels << -1;
	int *counts=(int *)malloc(sizeof(int)*K);

	//initialize the centroids
	QList<int> initial=choose_random_indices(N,K);
	for (int j=0; j<K; j++) {
		int ind=initial[j];
        int jj=ind*M;
        int ii=j*M;
		for (int m=0; m<M; m++) {
            centroids[m+ii]=Xptr[m+jj];
		}
	}

	bool something_changed=true;
	while (something_changed) {
		something_changed=false;
		//Assign the labels
		for (int n=0; n<N; n++) {
			int jj=n*M;
			double best_distsqr=0;
			int best_k=0;
			for (int k=0; k<K; k++) {
				int ii=k*M;
				double tmp=0;
				for (int m=0; m<M; m++) {
					tmp+=(centroids[m+ii]-Xptr[m+jj])*(centroids[m+ii]-Xptr[m+jj]);
				}
				if ((k==0)||(tmp<best_distsqr)) {
					best_distsqr=tmp;
					best_k=k;
				}
			}
			if (labels[n]!=best_k) {
				labels[n]=best_k;
				something_changed=true;
			}
		}

		if (something_changed) {
			//Compute the centroids
			for (int k=0; k<K; k++) {
				int ii=k*M;
				for (int m=0; m<M; m++) {
					centroids[m+ii]=0;
				}
				counts[k]=0;
			}
			for (int n=0; n<N; n++) {
				int jj=n*M;
				int k=labels[n];
				int ii=k*M;
				for (int m=0; m<M; m++) {
					centroids[m+ii]+=Xptr[m+jj];
				}
				counts[k]++;
			}
			for (int k=0; k<K; k++) {
				int ii=k*M;
				if (counts[k]) {
					for (int m=0; m<M; m++)
						centroids[m+ii]/=counts[k];
				}
			}
		}
	}

	free(counts);

	return labels;
}

//compute the centroids: X is MxN, labels are between 0 and K-1, and the return array is MxK
Mda compute_centroids(Mda &X,const QVector<int> &labels,int K) {
    double *Xptr=X.dataPtr();
    int M=X.N1();
    int N=X.N2();
    Mda ret; ret.allocate(M,K); double *retptr=ret.dataPtr();
    if (labels.count()!=N) {
        qWarning() << "Unexpected problem in compute_centroids at line:" << __LINE__;
        return ret;
    }

    int *counts=(int *)malloc(sizeof(int)*K);
	for (int k=0; k<K; k++) counts[k]=0;
	


	for (int n=0; n<N; n++) {
		int ii=n*M;
		int k=labels[n];
		int jj=k*M;
		for (int m=0; m<M; m++) {
			retptr[m+jj]+=Xptr[m+ii];
		}
		counts[k]++;
	}
	for (int k=0; k<K; k++) {
		if (counts[k]) {
			int jj=k*M;
			for (int m=0; m<M; m++) {
				 retptr[m+jj]/=counts[k];
			}
		}
	}

	free(counts);

	return ret;
}

//Compute the square matrix of distances between K vectors in M-dimensional space. Centroids is MxK and the output is KxK
Mda compute_distances(Mda &centroids) {
	int M=centroids.N1();
	int K=centroids.N2();
	double *Cptr=centroids.dataPtr();
	Mda ret; ret.allocate(K,K);
	for (int k2=0; k2<K; k2++) {
		int ii2=k2*M;
		for (int k1=k2; k1<K; k1++) {
			int ii1=k1*M;
			double tmp=0;
			for (int m=0; m<M; m++) {
				tmp+=(Cptr[m+ii1]-Cptr[m+ii2])*(Cptr[m+ii1]-Cptr[m+ii2]);
			}
			double dist=sqrt(tmp);
			ret.setValue(dist,k1,k2);
			ret.setValue(dist,k2,k1);
		}
	}
	return ret;
}

//find the best pair of labels (label1,label2) for the next iteration of isosplit
//Chooses them based on minimizing the distance (distances: KxK) for the active labels not yet attempted
//If none found, returns label1=label2=-1
void find_best_pair(int &label1,int &label2,bool *active_labels,Mda &distances,QSet<double> &attempted_redistributions) {
	double best_dist=-1;
	int K=distances.N1();
	double *Dptr=distances.dataPtr();
	label1=-1;
	label2=-1;
	for (int k2=0; k2<K; k2++) {
		if (active_labels[k2]) {
			int ii2=k2*K;
			for (int k1=k2+1; k1<K; k1++) {
				if (active_labels[k1]) {
					double dist=Dptr[k1+ii2];
					if (dist>=0) {
						if (!attempted_redistributions.contains(dist)) {
							if ((best_dist<0)||(dist<best_dist)) {
								best_dist=dist;
								label1=k1;
								label2=k2;
							}
						}
					}
				}
			}
		}
	}
}

//Return the vector of all the indices where labels[i]==label
QVector<int> find_inds_for_label(const QVector<int> &labels,int label) {
	QVector<int> ret;
	for (int i=0; i<labels.count(); i++)
		if (labels[i]==label) ret << i;
	return ret;
}

void attempt_to_redistribute_two_clusters(QVector<int> &ii1,QVector<int> &ii2,bool &redistributed,  Mda &X,const QVector<int> &inds1,const QVector<int> &inds2,double *C1,double *C2,double split_threshold) {
	ii1.clear();
	ii2.clear();
	redistributed=false;

	int M=X.N1();
	//int N=X.N2();
	double *Xptr=X.dataPtr();
	QVector<int> inds12=inds1; inds12.append(inds2);
	//X1=X(:,inds1);
	//X2=X(:,inds2);
	//the vector from one centroid to another
	QVector<double> V(M);
	for (int i=0; i<M; i++) V[i]=C1[i]-C2[i];
	double sumsqrV=0; for (int i=0; i<M; i++) sumsqrV+=V[i]*V[i];
	if (!sumsqrV) {
		qWarning() << "ISOSPLIT: vector V is null";
		return;
	}
	for (int i=0; i<M; i++) V[i]/=sqrt(sumsqrV);

	//project onto the line connecting the centroids
	int NN=inds12.count();
	double *XX=(double *)malloc(sizeof(double)*NN);
	int *labels=(int *)malloc(sizeof(int)*NN);
	for (int ii=0; ii<NN; ii++) {
		double tmp=0;
		int aa=inds12[ii]*M;
		for (int m=0; m<M; m++) {
			tmp+=Xptr[m+aa]*V[m];
		}
		XX[ii]=tmp;
	}

	//This is the core procedure -- split based on isotonic regression
	double pp;
	if (!isosplit1d(NN,labels,pp,XX)) {
		qWarning() << "Unexpected problem in isosplit1d: NN =" << NN;
	}
	if (pp>split_threshold) {
        //It was a statistically significant split -- so let's redistribute!
		for (int ii=0; ii<NN; ii++) {
			if (labels[ii]==1) {
				ii1 << inds12[ii];
			}
			else {
				ii2 << inds12[ii];
			}
		}
	}
	else {
		//Otherwise, merge into the first
		for (int ii=0; ii<NN; ii++) {
			ii1 << inds12[ii];
		}
	}

	//We have redistributed, unless everything is still the same.
	redistributed=true;
	if ((ii1.count()==inds1.count())&&(ii2.count()==inds2.count())) {
		QSet<int> s1=QSet<int>::fromList(ii1.toList());
		bool same1=true;
		for (int ii=0; ii<inds1.count(); ii++) {
			if (!s1.contains(inds1[ii])) same1=false;
		}
		if (same1) {
			QSet<int> s2=QSet<int>::fromList(ii2.toList());
			bool same2=true;
			for (int ii=0; ii<inds2.count(); ii++) {
				if (!s2.contains(inds2[ii])) same2=false;
			}
			if (same2) redistributed=false;
		}
	}

	free(XX);
	free(labels);
}

//compute the distance between two centrods. centroids is MxK, k1,k2 are between 0 and K-1
double compute_centroid_distance(Mda &centroids,int k1,int k2) {
	int M=centroids.N1();
	//int K=centroids.N2();
	double *Cptr=centroids.dataPtr();
	int ii1=k1*M;
	int ii2=k2*M;
	double ret=0;
	for (int m=0; m<M; m++) {
		ret+=(Cptr[m+ii1]-Cptr[m+ii2])*(Cptr[m+ii1]-Cptr[m+ii2]);
	}
	return sqrt(ret);
}

QVector<int> isosplit(Mda &X) {
    int K_initial=25;
	double split_threshold=0.9;

	int M=X.N1();
    int N=X.N2();
	QVector<int> labels=do_kmeans(X,K_initial);

	bool active_labels[K_initial];
	for (int ii=0; ii<K_initial; ii++) active_labels[ii]=true;
	Mda centroids=compute_centroids(X,labels,K_initial); //M x K_initial
	Mda distances=compute_distances(centroids); //K_initial x K_initial

	double *Cptr=centroids.dataPtr();

	//Here is a set of codes for the attempted cluster splits/redistributions -- we don't
	//ever want to repeat any of these.
	QSet<double> attempted_redistributions;

	int num_iterations=0;
	while (true) {
		num_iterations++;
        //if (num_iterations>35) return labels;
		QVector<int> old_labels=labels;
		//find the closest two clusters to check for redistribution/merging
		//we want the centroids to be as close as possible, but we want to
		//exclude the pairs we have tried previously
		int label1,label2;
		find_best_pair(label1,label2,active_labels,distances,attempted_redistributions);
		if (label1<0) break; //this means we've tried everything... so we are done!
		attempted_redistributions.insert(distances.value(label1,label2)); //we assume the distance between centroids uniquely defines this attempt!
		QVector<int> inds1=find_inds_for_label(labels,label1);
		QVector<int> inds2=find_inds_for_label(labels,label2);
		QVector<int> ii1,ii2;
		bool redistributed;
        attempt_to_redistribute_two_clusters(ii1,ii2,redistributed, X,inds1,inds2,&Cptr[label1*M],&Cptr[label2*M],split_threshold);
		if (redistributed) {
			//okay, we've changed something. Now let's updated the labels
			for (int jj=0; jj<ii1.count(); jj++) {
				labels[ii1[jj]]=label1;
			}
			for (int jj=0; jj<ii2.count(); jj++) {
				labels[ii2[jj]]=label2;
			}
			//recompute the centroids
			{
				double tmp1[M];
                for (int jj=0; jj<M; jj++) tmp1[jj]=0;
				for (int jj=0; jj<ii1.count(); jj++) {
					for (int mm=0; mm<M; mm++) {
						tmp1[mm]+=X.value(mm,ii1[jj]);
					}
				}
				for (int mm=0; mm<M; mm++) {
                    if (ii1.count())
                        centroids.setValue(tmp1[mm]/ii1.count(),mm,label1);
				}
			}
			if (ii2.count()>0) {
				double tmp2[M];
                for (int jj=0; jj<M; jj++) tmp2[jj]=0;
				for (int jj=0; jj<ii2.count(); jj++) {
					for (int mm=0; mm<M; mm++) {
						tmp2[mm]+=X.value(mm,ii2[jj]);
					}
				}
				for (int mm=0; mm<M; mm++) {
                    if (ii2.count())
                        centroids.setValue(tmp2[mm]/ii2.count(),mm,label2);
				}
			}
			else {
				for (int mm=0; mm<M; mm++) {
					centroids.setValue(0,mm,label2);
				}
				active_labels[label2]=false;
			}
			//recompute the distances
			for (int ii=0; ii<K_initial; ii++) {
				if (active_labels[ii]) {
					double val1=compute_centroid_distance(centroids,ii,label1);
					distances.setValue(val1,ii,label1);
					distances.setValue(val1,label1,ii);
					double val2=compute_centroid_distance(centroids,ii,label2);
					distances.setValue(val2,ii,label2);
					distances.setValue(val2,label2,ii);
				}
				else {
					distances.setValue(-1,ii,label1);
					distances.setValue(-1,ii,label2);
				}
			}
		}
		else {
			distances.setValue(-1,label1,label2);
		}
	}
    //remap the labels
    int labels_map[K_initial];
    int kk=1;
    for (int k=0; k<K_initial; k++) {
        if (active_labels[k]) {
            labels_map[k]=kk;
            kk++;
        }
        else labels_map[k]=0;
    }
    QVector<int> labels2;
    for (int n=0; n<N; n++) labels2 << labels_map[labels[n]];
    return labels2;
}
