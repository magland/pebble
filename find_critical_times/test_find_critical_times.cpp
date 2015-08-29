#include "find_critical_times.h"
#include <QVector>
#include <QDebug>
#include <QTime>

int main(int argc,char *argv[]) {
	{
		int N=25;
		const double X[]={1,2,3,2,1,2,3,4,5,4,3,2,1,2,3,4,5,6,7,6,5,4,3,2,1};
		qDebug() << find_critical_times(N,X,3);
		qDebug() << find_critical_times_brute_force(N,X,3);
	}

	{
		int NN=1e6;
		double A[NN];
		for (int i=0; i<NN; i++) A[i]=rand()*1.0/rand();
		QTime timer;
		
		timer.start();
		for (int jj=0; jj<15; jj++) {
			find_critical_times(NN,A,50);
		}
		qDebug() << "Elapsed time for find_critical_times (ms): " << timer.elapsed();

		timer.start();
		for (int jj=0; jj<15; jj++) {
			find_critical_times_brute_force(NN,A,50);
		}
		qDebug() << "Elapsed time for find_critical_times_brute_force (ms): " << timer.elapsed();
	}
	
	return 0;
}
