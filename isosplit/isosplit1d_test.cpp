#include <stdio.h>
#include "isosplit1d.h"
#include "mda.h"
#include <QDebug>
#include <QCoreApplication>
#include <QTime>
#include "isosplit.h"
#include <time.h>

int main(int argc,char *argv[]) {
	QCoreApplication app(argc,argv);

    srand (time(NULL));

    if (0) {
        int N=100000;
        int *labels=(int *)malloc(sizeof(int)*N);
        double *X=(double *)malloc(sizeof(double)*N);
        double pp;
        for (int i=0; i<N; i++) {
            double tmp=i*(1.0/N); //uniform between 0 and 1
            if (tmp<0.5) X[i]=tmp*tmp;
            else X[i]=(tmp*tmp+2);
        }
        QTime timer; timer.start();
        isosplit1d(N,labels,pp,X);
        qDebug() << pp << timer.elapsed();
        free(labels);
        free(X);
    }

    if (1) {
        int M=4;
        int N=500;
        Mda A; A.allocate(M,N);
        for (int n=0; n<N; n++) {
            for (int m=0; m<M; m++) {
                double val=(rand()*1.0/RAND_MAX);
                if ((n<100)&&(m==0)) val=val+3;
                A.setValue(val,m,n);
            }
        }
        QVector<int> labels=isosplit(A);
        qDebug() << labels;
    }

	return 0;
}

