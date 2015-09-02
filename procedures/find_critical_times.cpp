#include "find_critical_times.h"

QVector<int> find_critical_times_brute_force(Mda &X,int radius) {
	int N=X.N2();
	QVector<int> ret;
	double *Xptr=X.dataPtr();
	for (int t=0; t<N; t++) {
		double maxval=0;
		int maxind=t;
		for (int dt=-radius; dt<=radius; dt++) {
			int t0=t+dt;
			if ((0<=t0)&&(t0<N)) {
				if (qAbs(Xptr[t0])>maxval) {
					maxval=qAbs(Xptr[t0]);
					maxind=t0;
				}
			}
		}
		if (maxind==t) {
			if ((radius<=maxind)&&(maxind+radius<N))
				ret << t; //add this index if it is the maximum within the radius
		}

	}
	return ret;
}

QVector<int> find_critical_times(Mda &X,int radius) {
	int N=X.N2();
	QVector<int> ret;
	double *Xptr=X.dataPtr();
	if (N<=0) return ret;
	int last_best_ind=0; //the index of last maximum within the radius
	double last_best_val=qAbs(Xptr[0]); //the value of the last maximum within the radius
	bool is_pseudo=false; //it true, this means that the above is the last maximum within the radius, but it cannot be counted as a critical point
	for (int i=1; i<N; i++) {
		if (qAbs(Xptr[i])>last_best_val) {
			//we have a new maximum!
			if ((!is_pseudo)&&(i-last_best_ind>radius)) { //if the old one is not a pseudo, and it is now outside the range, add that one to our list
				if ((radius<=last_best_ind)&&(last_best_ind+radius<N))
					ret << last_best_ind;
			}
			last_best_ind=i;
			last_best_val=qAbs(Xptr[i]);
			is_pseudo=false; //not pseudo, because it was bigger than the previous biggest within the radius
		}
		else {
			//not a new maximum
			if (i-last_best_ind>radius) { //the old maximum is now outside the radius
				if (!is_pseudo) { //add it to our list if it wasn't a pseudo
					if ((radius<=last_best_ind)&&(last_best_ind+radius<N))
						ret << last_best_ind;
				}
				double maxval=0; //we need to find the new best
				int maxind=i;
				for (int dt=-radius; dt<=0; dt++) {
					int t0=i+dt;
					if (0<=t0) {
						if (qAbs(Xptr[t0])>maxval) {
							maxval=qAbs(Xptr[t0]);
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
		if ((radius<=last_best_ind)&&(last_best_ind+radius<N))
			ret << last_best_ind;
	}
	return ret;
}
