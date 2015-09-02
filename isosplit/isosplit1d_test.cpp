#include <stdio.h>
#include "isosplit1d.h"
#include "mda.h"
#include <QDebug>
#include <QCoreApplication>
#include <QTime>

int main(int argc,char *argv[]) {
	QCoreApplication app(argc,argv);
	int N=100000;
	int *labels=(int *)malloc(sizeof(int)*N);
	double *X=(double *)malloc(sizeof(double)*N);
	double pp;
	for (int i=0; i<N; i++) {
		double tmp=i*(1.0/N); //uniform between 0 and 1
		if (tmp<0.5) X[i]=tmp*tmp;
		else X[i]=(tmp*tmp+0.001);
	}
	QTime timer; timer.start();
	isosplit1d(N,labels,pp,X);
	qDebug() << pp << timer.elapsed();
	free(labels);
	free(X);
	return 0;
}

