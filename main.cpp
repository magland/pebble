//#include <QCoreApplication>
#include "mda.h"
#include <QString>
#include <QDebug>

void split_into_channels(const QString &inpath,const QString &outpath) {
	Mda X;
	X.read(inpath.toLatin1().data());
	qDebug() << X.N1() << X.N2() << X.N3();
}

int main(int argc, char *argv[])
{
	//QCoreApplication a(argc, argv);

	split_into_channels(
				"/home/magland/data/EJ/Spikes_all_channels_filtered.mda",
				"/home/magland/data/EJ/channels"
				);

	//return a.exec();

	return 0;
}
