#include <QCoreApplication>
#include "isosplit1d.h"
#include <QVector>
#include "mda.h"
#include <math.h>
#include "jisotonic.h"
#include <QDir>
#include <QMap>
#include <stdio.h>
#include <QStringList>
#include <QDebug>
#include <QTime>

void jisotonic_updown(int N,double *out,double *in,double *weights) {
	double *B1=(double *)malloc(sizeof(double)*N);
	double *MSE1=(double *)malloc(sizeof(double)*N);
	double *B2=(double *)malloc(sizeof(double)*N);
	double *MSE2=(double *)malloc(sizeof(double)*N);
	double *in_reversed=(double *)malloc(sizeof(double)*N);
    double *weights_reversed=0;

    for (int j=0; j<N; j++) {in_reversed[j]=in[N-1-j];}
	if (weights) {
        weights_reversed=(double *)malloc(sizeof(double)*N);
        for (int j=0; j<N; j++) {weights_reversed[j]=weights[N-1-j];}
	}
	jisotonic(N,B1,MSE1,in,weights);
	jisotonic(N,B2,MSE2,in_reversed,weights_reversed);
	for (int j=0; j<N; j++) MSE1[j]+=MSE2[N-1-j];
	double bestval=MSE1[0];
	int best_ind=0;
	for (int j=0; j<N; j++) {
		if (MSE1[j]<bestval) {
			bestval=MSE1[j];
			best_ind=j;
		}
	}
	jisotonic(best_ind+1,B1,MSE1,in,weights);
	jisotonic(N-best_ind,B2,MSE2,in_reversed,weights_reversed);
	for (int j=0; j<=best_ind; j++) out[j]=B1[j];
	for (int j=0; j<N-best_ind-1; j++) out[N-1-j]=B2[j];

	free(B1);
	free(MSE1);
	free(B2);
	free(MSE2);
	free(in_reversed);
	if (weights_reversed) free(weights_reversed);
}

void jisotonic_downup(int N,double *out,double *in,double *weights) {
	double *in_neg=(double *)malloc(sizeof(double)*N);

	for (int j=0; j<N; j++) in_neg[j]=-in[j];
	jisotonic_updown(N,out,in_neg,weights);
	for (int j=0; j<N; j++) out[j]=-out[j];

	free(in_neg);
}

void quick_sort (double *a, int n) {
	int i, j;
	double p,t;
	if (n < 2)
		return;
	p = a[n / 2];
	for (i = 0, j = n - 1;; i++, j--) {
		while (a[i] < p)
			i++;
		while (p < a[j])
			j--;
		if (i >= j)
			break;
		t = a[i];
		a[i] = a[j];
		a[j] = t;
	}
	quick_sort(a, i);
	quick_sort(&a[i], n - i);
}

void sort(int N,double *out,double *in) {
	for (int i=0; i<N; i++) out[i]=in[i];
	quick_sort(out,N);
//	QVector<double> in0(N);
//	for (int j=0; j<N; j++) in0[j]=in[j];
//	qSort(in0);
//	for (int j=0; j<N; j++) out[j]=in0[j];
}

bool get_isosplit_curve(int curve_len,int N,double *curve,double &cutpoint,double *X_in,int minsize) {
	QTime timer; timer.start();
	if (N<minsize*2) {
		for (int jj=0; jj<curve_len; jj++) curve[jj]=0;
		cutpoint=0;
		return true;
	}
    //qDebug() << "ELAPSED" << __FILE__ << __LINE__ << timer.elapsed(); timer.start();
    Mda X_mda; X_mda.allocate(1,N); double *X=X_mda.dataPtr();
    sort(N,X,X_in);
    //qDebug() << "ELAPSED" << __FILE__ << __LINE__ << timer.elapsed(); timer.start();
	Mda spacings_mda; spacings_mda.allocate(1,N); double *spacings=spacings_mda.dataPtr();
	spacings[0]=X[1]-X[0];
	spacings[N-1]=X[N-1]-X[N-2];
	for (int j=1; j<N-1; j++) {
		spacings[j]=(X[j+1]-X[j-1])/2;
		if (spacings[j]<=0) {
            printf ("ERROR: in get_isosplit_curve, spacings must be strictly positive -- did you forget to sort, or are the datapoints not distinct? spacings[%d]=%g\n",j,spacings[j]);
			return false;
		}
	}
    //qDebug() << "ELAPSED" << __FILE__ << __LINE__ << timer.elapsed(); timer.start();
	Mda logdensity_mda; logdensity_mda.allocate(1,N); double *logdensity=logdensity_mda.dataPtr();
	Mda sqrt_spacings_mda; sqrt_spacings_mda.allocate(1,N); double *sqrt_spacings=sqrt_spacings_mda.dataPtr();
    //qDebug() << "ELAPSED" << __FILE__ << __LINE__ << timer.elapsed(); timer.start();
	for (int j=0; j<N; j++) {
		logdensity[j]=log(1/spacings[j]);
		sqrt_spacings[j]=sqrt(spacings[j]);
	}
    //qDebug() << "ELAPSED" << __FILE__ << __LINE__ << timer.elapsed(); timer.start();
	Mda fit1_mda; fit1_mda.allocate(1,N); double *fit1=fit1_mda.dataPtr();
	Mda resid1_mda; resid1_mda.allocate(1,N); double *resid1=resid1_mda.dataPtr();
	Mda fit2_mda; fit2_mda.allocate(1,N); double *fit2=fit2_mda.dataPtr();
	Mda fit2_sorted_mda; fit2_sorted_mda.allocate(1,N); double *fit2_sorted=fit2_sorted_mda.dataPtr();
    //qDebug() << "ELAPSED" << __FILE__ << __LINE__ << timer.elapsed(); timer.start();

	jisotonic_updown(N,fit1,logdensity,0);

    //qDebug() << "ELAPSED" << __FILE__ << __LINE__ << timer.elapsed(); timer.start();
	for (int j=0; j<N; j++) resid1[j]=logdensity[j]-fit1[j];
	jisotonic_downup(N,fit2,resid1,sqrt_spacings);
    //qDebug() << "ELAPSED" << __FILE__ << __LINE__ << timer.elapsed(); timer.start();
	double mean_fit2=0;
	for (int j=0; j<N; j++) {
		mean_fit2+=fit2[j];
	}
	mean_fit2/=N;
	for (int j=0; j<N; j++) {
		fit2[j]-=mean_fit2;
	}
	fit2=&fit2[minsize-1];
	N=N-2*minsize+2;
    //qDebug() << "ELAPSED" << __FILE__ << __LINE__ << timer.elapsed(); timer.start();
	sort(N,fit2_sorted,fit2);
    //qDebug() << "ELAPSED" << __FILE__ << __LINE__ << timer.elapsed(); timer.start();
	double sum0=0;
	for (int k=0; k<curve_len; k++) {
		if (k<curve_len) {
			sum0+=fit2_sorted[k];
			curve[k]=sum0/(k+1);
		}
		else curve[k]=0;
	}
	double minval=fit2[0];
	int minind=0;
	for (int j=0; j<N; j++) {
		if (fit2[j]<minval) {
			minval=fit2[j];
			minind=j;
		}
	}
    //qDebug() << "ELAPSED" << __FILE__ << __LINE__ << timer.elapsed(); timer.start();
	minind+=minsize-1;
	if (minind==0) {cutpoint=X[0]; return true;}
	if (minind==N-1) {cutpoint=X[N-1]; return true;}
	if (spacings[minind-1]>spacings[minind+1]) {
		cutpoint=(X[minind-1]+X[minind])/2;
		return true;
	}
	else {
		cutpoint=(X[minind]+X[minind+1])/2;
		return true;
	}
}

struct calibration_file {
	int curve_len;
	int num_trials;
	QVector<double> avg;
	QVector<double> stdev;
	QVector<double> scores;
};
static QMap<int,calibration_file> s_calibration_files;
static QString s_calibration_dir;
static QList<int> s_calibration_Ns;


QString find_directory_path(QString base_path,QString dirname) {
	QStringList folders=QDir(base_path).entryList(QStringList("*"),QDir::Dirs|QDir::NoDotAndDotDot,QDir::Name);
	if (folders.contains(dirname)) return base_path+"/"+dirname;
	for (int i=0; i<folders.count(); i++) {
		QString tmp=find_directory_path(base_path+"/"+folders[i],dirname);
		if (!tmp.isEmpty()) return tmp;
	}
	return "";
}

void set_calibration_dir() {
	QString base_path=qApp->applicationDirPath()+"/..";
	QString path=find_directory_path(base_path,"isosplit_calibration");
	if (path.isEmpty()) return;
	s_calibration_dir=path;
}

void set_calibration_Ns() {
	s_calibration_Ns.clear();
	if (s_calibration_dir.isEmpty()) return;
	QStringList list=QDir(s_calibration_dir).entryList(QStringList("*"),QDir::Files,QDir::Name);
	for (int i=0; i<list.count(); i++) {
		QString fname=list[i];
		int ind1=fname.lastIndexOf("_");
		int ind2=fname.lastIndexOf(".");
		if ((ind1>=0)&&(ind2>=ind1+2)) {
			int N0=fname.mid(ind1+1,ind2-ind1-1).toInt();
			if (N0>0) s_calibration_Ns << N0;
		}
	}
	qSort(s_calibration_Ns);
}

void load_calibration_file(int N) {
	QString path=QString("%1/isosplit_calibration_%2.txt").arg(s_calibration_dir).arg(N);
	QFile FF(path);
	if (FF.open(QFile::ReadOnly|QFile::Text)) {
		QString txt=QString(FF.readAll());
		FF.close();
		QList<QString> lines=txt.split("\n");
		if (lines.count()>10) {
			calibration_file CF;
			CF.curve_len=lines[0].toInt();
			CF.num_trials=lines[1].toInt();
			for (int j=0; j<CF.curve_len; j++) {
				QString line=lines.value(2+j);
				QList<QString> tmp=line.split(",");
				if (tmp.count()!=2) return;
				CF.avg << tmp.value(0).toDouble();
				CF.stdev << tmp.value(1).toDouble();
			}
			for (int j=0; j<CF.num_trials; j++) {
				QString line=lines.value(2+CF.curve_len+j);
				CF.scores << line.toDouble();
			}
			s_calibration_files[N]=CF;
		}
		else return;
	}

}

bool read_isosplit_calibration(int N,QVector<double> &avg,QVector<double> &stdev,QVector<double> &scores) {
	if (s_calibration_dir.isEmpty()) set_calibration_dir();
    if (s_calibration_dir.isEmpty()) {printf ("Error: unable to find isosplit calibration directory (isosplit_calibration)\n"); return false;}
	if (s_calibration_Ns.isEmpty()) set_calibration_Ns();
	int ii=-1;
	while ((ii+2<s_calibration_Ns.count())&&(s_calibration_Ns[ii+1]<=N)) ii++;
    if ((ii<0)||(ii+1>=s_calibration_Ns.count())) {printf ("Error: calibration out of range (N=%d)\n",N); return false;}
	int N1=s_calibration_Ns[ii];
	int N2=s_calibration_Ns[ii+1];
	if (N1==N) N2=N1;
	if (!s_calibration_files.contains(N1)) {
		load_calibration_file(N1);
        if (!s_calibration_files.contains(N1)) {printf ("Error: unable to find calibration file (N=%d)\n",N1); return false;}
	}
	if (!s_calibration_files.contains(N2)) {
		load_calibration_file(N2);
        if (!s_calibration_files.contains(N2)) {printf ("Error: unable to find calibration file (N=%d)\n",N2); return false;}
	}
	if (N1==N2) {
		calibration_file *CF=&s_calibration_files[N1];
		avg=CF->avg; stdev=CF->stdev; scores=CF->scores;
	}
	else {
		double pct=(N-N1)*1.0/(N2-N1);
		calibration_file *CF1=&s_calibration_files[N1];
		calibration_file *CF2=&s_calibration_files[N2];
		avg=CF1->avg; stdev=CF1->stdev; scores=CF1->scores;
		for (int i=0; i<avg.count(); i++) {
			avg[i]=avg[i]+pct*(CF2->avg[i]-avg[i]);
			stdev[i]=stdev[i]+pct*(CF2->stdev[i]-stdev[i]);
		}
		for (int i=0; i<scores.count(); i++) {
			scores[i]=scores[i]+pct*(CF2->scores[i]-scores[i]);
		}
	}
	return true;
}

bool isosplit1d(int N,int *labels,double &pp,double *X) {
	QTime timer; timer.start();
	int minsize=4; //hard coded to be consistent with calibration

    //qDebug() << "ELAPSED" << __FILE__ << __LINE__ << timer.elapsed(); timer.start();
	QVector<double> avg(N),stdev(N),scores(N);
	if (!read_isosplit_calibration(N,avg,stdev,scores)) {
        printf ("ERROR: problem in read_isosplit_calibration.\n");
		return false;
	}
	int curve_len=avg.count();
	if (curve_len==0) {
        printf ("ERROR: calibration not performed for N=%d\n",N);
		return false;
	}
	int num_trials=scores.count();
    //qDebug() << "ELAPSED" << __FILE__ << __LINE__ << timer.elapsed(); timer.start();

	//compute the curve and cutpoint
	Mda curve_mda; curve_mda.allocate(1,curve_len); double *curve=curve_mda.dataPtr();
	double cutpoint;
    //qDebug() << "ELAPSED" << __FILE__ << __LINE__ << timer.elapsed(); timer.start();
	if (!get_isosplit_curve(curve_len,N,curve,cutpoint,X,minsize)) {
        printf ("ERROR in get_isosplit_curve\n");
		return false;
	}
    //qDebug() << "ELAPSED" << __FILE__ << __LINE__ << timer.elapsed(); timer.start();

	//compute the score by comparing to calibration avg and stdev and then compute the pp
	double score0=0;
	for (int ii=0; ii<curve_len; ii++) {
		if (stdev[ii]) {
			double diff0=avg[ii]-curve[ii];
			//we need a change of at least 0.1 to count it
			//this is important in the case where stdev is extremely close to 0
			if (qAbs(diff0)<0.1) diff0=0;
			double tmp_score=diff0/stdev[ii];
			if (tmp_score>score0) score0=tmp_score;
		}
	}
	if (score0!=0) {
		double tmp_ct=0;
		for (int jj=0; jj<scores.count(); jj++) {
			if (scores[jj]<score0) tmp_ct++;
		}
		pp=tmp_ct/num_trials;
	}
	else {
		pp=0;
	}
    //qDebug() << "ELAPSED" << __FILE__ << __LINE__ << timer.elapsed(); timer.start();

	//set the labels
	for (int kk=0; kk<N; kk++) {
		if (X[kk]<cutpoint) labels[kk]=1;
		else labels[kk]=2;
	}
    //qDebug() << "ELAPSED" << __FILE__ << __LINE__ << timer.elapsed(); timer.start();

	return true;
}
