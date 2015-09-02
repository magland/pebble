#include "isosplit.h"
#include <QSet>

QList<int> do_kmeans(Mda &X,K) {
	int M=X.N1();
	int N=X.N2();
	double *Xptr=X.dataPtr();
	Mda centroids_mda; centroids_mda.allocate(M,K); double *centroids=centroids_mda.dataPtr();

	//initialize the centroids
	QList<int> initial=choose_random_indices(N,K);
	int ct1=0;
	for (int j=0; j<K; j++) {
		int ind=initial[j];
		int ct2=Xptr[ind*M];
		for (int m=0; m<M; m++) {
			centroids[ct1]=Xptr[ct2];
			ct1++; ct2++;
		}
	}
}

QList<int> isosplit(Mda &X) {
	int K_initial=10;
	double split_threshold=0.9;

	int M=X.N1();
	int N=X.N2();
	QList<int> labels;

	labels=do_kmeans(X,K_initial);

	bool active_labels[K_initial];
	for (int ii=0; ii<K_initial; ii++) active_labels[ii]=true;
	Mda centroids=compute_centroids(X,labels,K_initial); //M x K_initial
	Mda distances=compute_distances(centroids); //K_initial x K_initial

	//Here is a set of codes for the attempted cluster splits/redistributions -- we don't
	//ever want to repeat any of these.
	QSet<double> attempted_redistributions;

	int num_iterations=0;
	while (true) {
		num_iterations++;
		QList<int> old_labels=labels;
		//find the closest two clusters to check for redistribution/merging
		//we want the centroids to be as close as possible, but we want to
		//exclude the pairs we have tried previously
		int label1,label2;
		find_best_pair(label1,label2,active_labels,distances,attempted_redistributions);
		if (label1==0) break; //this means we've tried everything... so we are done!
		attempted_redistributions.insert(distances.value(label1,label2)); //we assume the distance between centroids uniquely defines this attempt!
		QList<int> inds1=find_inds_for_label(labels,label1);
		QList<int> inds2=find_inds_for_label(labels,label2);
		QList<int> ii1,ii2;
		bool redistributed;
		attempt_to_redistribute_two_clusters(ii1,ii2,redistributed, X,inds1,inds2,centroid1,centroid2,split_threshold);
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
				for (int jj=0; jj<ii1.count(); jj++) {
					for (int mm=0; mm<M; mm++) {
						tmp1[mm]+=X.value(mm,ii1[jj]);
					}
				}
				for (int mm=0; mm<M; mm++) {
					centroids.setValue(tmp1[mm]/M,mm,label1);
				}
			}
			if (ii2.count()>0) {
				double tmp2[M];
				for (int jj=0; jj<ii2.count(); jj++) {
					for (int mm=0; mm<M; mm++) {
						tmp2[mm]+=X.value(mm,ii2[jj]);
					}
				}
				for (int mm=0; mm<M; mm++) {
					centroids.setValue(tmp2[mm]/M,mm,label2);
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
}

/*

	if (redistributed)
		%okay, we've changed something. Now let's update the labels
		if (opts.verbose)
			fprintf('  Redistributed (%d,%d)->(%d,%d)\n',length(inds1),length(inds2),length(ii1),length(ii2));
		end;
		labels(ii1)=label1;
		labels(ii2)=label2;
		centroids(:,label1)=mean(X(:,ii1),2);
		tmp0=sum((centroids-repmat(centroids(:,label1),1,size(centroids,2))).^2,1);
		tmp0(ismember(tmp0,attempted_redistributions))=inf;
		distances(:,label1)=tmp0;
		distances(label1,:)=distances(:,label1); distances(label1,label1)=inf;
		if (length(ii2>0))
			centroids(:,label2)=mean(X(:,ii2),2);
			tmp0=sum((centroids-repmat(centroids(:,label2),1,size(centroids,2))).^2,1);
			tmp0(ismember(tmp0,attempted_redistributions))=inf;
			distances(:,label2)=tmp0;
			distances(label2,:)=distances(:,label2); distances(label2,label2)=inf;
		else
			centroids(:,label2)=0;
		end;
		[labels,centroids,distances]=normalize_labels(labels,centroids,distances); % we may have eliminated a label, so let's shift the labelings down
	else
		distances(label1,label2)=inf;
		distances(label2,label1)=inf;
	end;

	% The following might be a bad idea.
	if (num_iterations_with_same_number_of_clusters>opts.max_iterations_per_number_clusters)
		warning(sprintf('%d iterations with same number of clusters.... stopping',num_iterations_with_same_number_of_clusters));
		return;
	end;
	%toc

*/

