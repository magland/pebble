#include "mda.h"

#include <stdlib.h>
#include <stdio.h>
#define UNUSED(expr) do { (void)(expr); } while (0);
#include <QDebug>
#include "usagetracking.h"

class MdaPrivate {
public:
	Mda *q;

	int *m_size;
	double *m_data_real;
	int m_data_type;	
	
	void construct() {
		m_size=(int *)jmalloc(sizeof(int)*MDA_MAX_DIMS);
		for (int i=0; i<MDA_MAX_DIMS; i++) m_size[i]=1;
		
		m_data_real=(double *)jmalloc(sizeof(float)*1);
		m_data_real[0]=0;
	
		m_data_type=MDA_TYPE_REAL;
	}
	
	bool do_read(FILE *inf);
	int read_int(FILE *inf);
	float read_float(FILE *inf);
	short read_short(FILE *inf);	
	unsigned short read_unsigned_short(FILE *inf);
	unsigned char read_unsigned_char(FILE *inf);
	bool do_write(FILE *outf);
	void write_int(FILE *outf,int val);
	void write_float(FILE *outf,float val);
	void write_short(FILE *outf,short val);	
	void write_unsigned_short(FILE *outf,unsigned short val);
	void write_unsigned_char(FILE *outf,unsigned char val);
};

Mda::Mda() 
{
	d=new MdaPrivate;
	d->q=this;
	
	d->construct();
}

Mda::Mda(const Mda &X) {
	d=new MdaPrivate;
	d->q=this;
	
	d->construct();
	
	setDataType(X.dataType());
	allocate(X.dimCount(),X.d->m_size);
	setValues(X.d->m_data_real);
}

Mda::~Mda()
{
	jfree(d->m_size);
	jfree(d->m_data_real);
	delete d;
}

void Mda::operator=(const Mda &X) {
	setDataType(X.dataType());
	allocate(X.dimCount(),X.d->m_size);
	setValues(X.d->m_data_real);
}

void Mda::setDataType(int dt) {
	d->m_data_type = dt;
	allocate(1, 1);
}

int Mda::dataType() const {
	return d->m_data_type;
}
void Mda::allocate(int N1,int N2,int N3,int N4,int N5,int N6) {
	int tmp[6];
	tmp[0]=N1; tmp[1]=N2; tmp[2]=N3; tmp[3]=N4; tmp[4]=N5; tmp[5]=N6;
	allocate(6,tmp);
}
void Mda::allocate(int num_dims,int *size) {
	for (int i=0; i<MDA_MAX_DIMS; i++) d->m_size[i]=1;
	
	int NN=1;
	for (int i = 0; i < num_dims; i++) {
		if (size[i] <= 0) {
			size[i] = 1;
		}
		d->m_size[i] = size[i];
		NN *= size[i];
	}
	if (NN >MDA_MAX_SIZE) {
		printf ("Unable to allocate mda. Size is too large: %d.\n", NN);
		allocate(1, 1);
		return;
	}
	if (NN > 0) {
		jfree(d->m_data_real);
		d->m_data_real = (double *)jmalloc(sizeof(double)*NN);
		for (int i=0; i<NN; i++) d->m_data_real[i]=0;
	}
}
int Mda::size(int dim) const {
	if (dim >= MDA_MAX_DIMS) {
		return 1;
	}
	if (dim < 0) {
		return 0;
	}
	return d->m_size[dim];
}
int Mda::N1() const {
	return size(0);
}
int Mda::N2() const {
	return size(1);
}
int Mda::N3() const {
	return size(2);
}
int Mda::N4() const {
	return size(3);
}
int Mda::N5() const {
	return size(4);
}
int Mda::N6() const {
	return size(5);
}
int Mda::dimCount() const {
	int ret = 2;
	for (int i = 2; i < MDA_MAX_DIMS; i++) {
		if (d->m_size[i] > 1) {
			ret = i + 1;
		}
	}
	return ret;
}
int Mda::totalSize() const {
	int ret = 1;
	for (int j = 0; j < MDA_MAX_DIMS; j++) {
		ret *= d->m_size[j];
	}
	return ret;
}

void Mda::reshape(int N1, int N2, int N3, int N4, int N5, int N6)
{
	if (N1*N2*N3*N4*N5*N6!=totalSize()) {
		qWarning() << "Unable to reshape. Inconsistent size." << N1 << N2 << N3 << N4 << N5 << N6 << totalSize();
		return;
	}
	d->m_size[0]=N1;
	d->m_size[1]=N2;
	d->m_size[2]=N3;
	d->m_size[3]=N4;
	d->m_size[4]=N5;
	d->m_size[5]=N6;
	for (int j=6; j<MDA_MAX_DIMS; j++) d->m_size[j]=1;
}
double Mda::value1(int i) const {
	return d->m_data_real[i];
}
double Mda::value(int i1,int i2,int i3,int i4) const {
	if ((i1<0)||(i1>=size(0))) return 0;
	if ((i2<0)||(i2>=size(1))) return 0;
	if ((i3<0)||(i3>=size(2))) return 0;
	if ((i4<0)||(i4>=size(3))) return 0;
	return d->m_data_real[i1+size(0)*i2+size(0)*size(1)*i3+size(0)*size(1)*size(2)*i4];
}
double Mda::value(int num_dims,int *ind) const {
	int tmp=0;
	int factor=1;
	for (int i=0; i<num_dims; i++) {
		if ((ind[i]<0)||(ind[i]>=size(i))) return 0;
		tmp+=ind[i]*factor;
		factor*=size(i);
	}
	return d->m_data_real[tmp];
}
void Mda::setValue1(double val,int i) {
	d->m_data_real[i]=val;
}
void Mda::setValue(double val,int i1,int i2,int i3,int i4) {
	if ((i1<0)||(i1>=size(0))) return;
	if ((i2<0)||(i2>=size(1))) return;
	if ((i3<0)||(i3>=size(2))) return;
	if ((i4<0)||(i4>=size(3))) return;
	d->m_data_real[i1+size(0)*i2+size(0)*size(1)*i3+size(0)*size(1)*size(2)*i4]=val;
}
void Mda::setValue(double val,int num_dims,int *ind) {
	int tmp=0;
	int factor=1;
	for (int i=0; i<num_dims; i++) {
		if ((ind[i]<0)||(ind[i]>=size(i))) return;
		tmp+=ind[i]*factor;
		factor*=size(i);
	}
	d->m_data_real[tmp]=val;
}
void Mda::setValues(double *vals) {
	int ts=totalSize();
	for (int i=0; i<ts; i++) {
		d->m_data_real[i]=vals[i];
	}
}
void Mda::setValues(int *vals) {
	int ts=totalSize();
	for (int i=0; i<ts; i++) {
		d->m_data_real[i]=(float)vals[i];
	}
}
void Mda::setValues(short *vals) {
	int ts=totalSize();
	for (int i=0; i<ts; i++) {
		d->m_data_real[i]=(float)vals[i];
	}
}
void Mda::setValues(unsigned char *vals) {
	int ts=totalSize();
	for (int i=0; i<ts; i++) {
		d->m_data_real[i]=(float)vals[i];
	}
}

Mda Mda::getDataXY(int num_inds,int *inds) const {
	Mda ret;
	ret.allocate(N1(), N2());
	int *inds0=(int *)jmalloc(sizeof(int)*(num_inds+2));
	for (int j = 0; j < num_inds; j++) {
		inds0[j + 2] = inds[j];
	}
	int n1=size(0);
	int n2=size(1);	
	for (int y = 0; y < n2; y++) {
		inds0[1] = y;
		for (int x = 0; x < n1; x++) {
			inds0[0] = x;
			ret.setValue(value((num_inds+2),inds0), x, y);
		}
	}
	jfree(inds0);
	return ret;
}
Mda Mda::getDataXZ(int num_inds,int *inds) const {
	Mda ret;
	ret.allocate(N1(), N3());
	int *inds0=(int *)jmalloc(sizeof(int)*(num_inds+2));
	inds0[1]=inds[0];
	for (int j = 1; j < num_inds; j++) {
		inds0[j + 2] = inds[j];
	}
	int n1=size(0);
	int n2=size(2);	
	for (int z = 0; z < n2; z++) {
		inds0[2] = z;
		for (int x = 0; x < n1; x++) {
			inds0[0] = x;
			ret.setValue(value((num_inds+2),inds0), x, z);
		}
	}
	jfree(inds0);
	return ret;
}
Mda Mda::getDataYZ(int num_inds,int *inds) const {
	Mda ret;
	ret.allocate(N2(), N3());
	int *inds0=(int *)jmalloc(sizeof(int)*(num_inds+2));
	inds0[0]=inds[0];
	for (int j = 1; j < num_inds; j++) {
		inds0[j + 2] = inds[j];
	}
	int n1=size(1);
	int n2=size(2);	
	for (int z = 0; z < n2; z++) {
		inds0[2] = z;
		for (int y = 0; y < n1; y++) {
			inds0[1] = y;
			ret.setValue(value((num_inds+2),inds0), y, z);
		}
	}
	jfree(inds0);
	return ret;
}
Mda Mda::transpose() const {
	Mda ret;
	int *size=(int *)jmalloc(sizeof(int)*MDA_MAX_DIMS);
	for (int i=0; i<MDA_MAX_DIMS; i++) size[i]=d->m_size[i];
	int s1=size[0];
	int s2=size[1];
	size[0] = s2;
	size[1] = s1;
	ret.allocate(MDA_MAX_DIMS,size);
	int tot_size = totalSize();
	int num_planes = tot_size / (size[0] * size[1]);
	for (int i = 0; i < num_planes; i++) {
		int offset = i * size[0] * size[1];
		for (int y = 0; y < d->m_size[1]; y++) {
			for (int x = 0; x < d->m_size[0]; x++) {
				ret.setValue1(value1(offset + x + y * d->m_size[0]), offset + y + x * d->m_size[1]);
			}
		}
	}
	jfree(size);
	return ret;
}
bool Mda::read(char *path) {
	allocate(1, 1);

	FILE *inf=jfopen(path,"rb");
	if (!inf) return false;
	
	bool ret=d->do_read(inf);
	
	jfclose(inf);
	
	return ret;
	
}
bool MdaPrivate::do_read(FILE *inf) {
	int hold_num_dims;
	int hold_dims[MDA_MAX_DIMS];
	for (int i=0; i<MDA_MAX_DIMS; i++) hold_dims[i]=1;
	
	hold_num_dims = read_int(inf);
	
	int data_type;
	if (hold_num_dims < 0) {
		data_type = hold_num_dims;
		int num_bytes = read_int(inf);
		UNUSED(num_bytes)
		hold_num_dims = read_int(inf);
	} else {
		data_type = MDA_TYPE_COMPLEX;
	}
	if (hold_num_dims > MDA_MAX_DIMS) {
		printf ("number of dimensions exceeds maximum: %d\n", hold_num_dims);
		return false;
	}
	if (hold_num_dims <= 0) {
		printf ("unexpected number of dimensions: %d\n", hold_num_dims);
		return false;
	}
	for (int j = 0; j < hold_num_dims; j++) {
		int holdval = read_int(inf);
		hold_dims[j] = holdval;
	}
	{
		q->setDataType(data_type);
		q->allocate(hold_num_dims,hold_dims);
		int N = q->totalSize();
		if (data_type == MDA_TYPE_COMPLEX) {
			for (int ii = 0; ii < N; ii++) {
				float re0 = read_float(inf);
				float im0 = read_float(inf);
				UNUSED(im0);
				m_data_real[ii] = re0;
			}
		} else if (data_type == MDA_TYPE_REAL) {
			for (int ii = 0; ii < N; ii++) {
				float re0 = read_float(inf);
				m_data_real[ii] = re0;
			}
		} else if (data_type == MDA_TYPE_SHORT) {
			for (int ii = 0; ii < N; ii++) {
				float re0 = read_short(inf);
				m_data_real[ii] = re0;
			}
		} else if (data_type == MDA_TYPE_UINT16) {
			for (int ii = 0; ii < N; ii++) {
				int re0 = read_unsigned_short(inf);
				m_data_real[ii] = re0;
			}
		} else if (data_type == MDA_TYPE_INT32) {
			for (int ii = 0; ii < N; ii++) {
				int re0 = read_int(inf);
				m_data_real[ii] = re0;
			}
		} else if (data_type == MDA_TYPE_BYTE) {
			for (int ii = 0; ii < N; ii++) {
				unsigned char re0 = read_unsigned_char(inf);
				m_data_real[ii] = re0;
			}
		} else {
			printf ("Unrecognized data type %d\n", data_type);
			return false;
		}
	}
	
	return true;
}
int MdaPrivate::read_int(FILE *inf) {
	int ret=0;
	int b0=jfread(&ret,4,1,inf);
	UNUSED(b0)
	return ret;
}
float MdaPrivate::read_float(FILE *inf) {
	float ret=0;
	int b0=jfread(&ret,4,1,inf);
	UNUSED(b0)
	return ret;
}
short MdaPrivate::read_short(FILE *inf) {
	short ret=0;
	int b0=jfread(&ret,2,1,inf);
	UNUSED(b0)
	return ret;
}
unsigned short MdaPrivate::read_unsigned_short(FILE *inf) {
	unsigned short ret=0;
	int b0=jfread(&ret,2,1,inf);
	UNUSED(b0)
	return ret;
}
unsigned char MdaPrivate::read_unsigned_char(FILE *inf) {
	unsigned char ret=0;
	int b0=jfread(&ret,1,1,inf);
	UNUSED(b0)
	return ret;
}

bool Mda::write(char *path) {

	FILE *outf=jfopen(path,"wb");
	if (!outf) {
		printf ("Unable to write mda file: %s\n",path);
		return false;
	}
	
	bool ret=d->do_write(outf);
	
	jfclose(outf);
	
	return ret;
}

bool MdaPrivate::do_write(FILE *outf) {
	write_int(outf, m_data_type);
	int num_bytes = 4;
	if (m_data_type == MDA_TYPE_COMPLEX) {
		num_bytes = 8;
	} else if (m_data_type == MDA_TYPE_BYTE) {
		num_bytes = 1;
	} else if (m_data_type == MDA_TYPE_SHORT) {
		num_bytes = 2;
	} else if (m_data_type == MDA_TYPE_UINT16) {
		num_bytes = 2;
	}
	write_int(outf, num_bytes);
	int num_dims = 2;
	for (int i = 2; i < MDA_MAX_DIMS; i++) {
		if (m_size[i] > 1) {
			num_dims = i + 1;
		}
	}
	write_int(outf, num_dims);
	for (int ii = 0; ii < num_dims; ii++) {
		write_int(outf, m_size[ii]);
	}
	int N = q->totalSize();
	if (m_data_type == MDA_TYPE_COMPLEX) {
		for (int i = 0; i < N; i++) {
			float re0 = (float) m_data_real[i];
			write_float(outf, re0);
			write_float(outf, 0);
		}
	} else if (m_data_type == MDA_TYPE_REAL) {
		for (int i = 0; i < N; i++) {
			float re0 = (float) m_data_real[i];
			write_float(outf, re0);
		}
	} else if (m_data_type == MDA_TYPE_BYTE) {
		for (int i = 0; i < N; i++) {
			unsigned char re0 = (unsigned char) m_data_real[i];
			write_unsigned_char(outf,re0);
		}
	} else if (m_data_type == MDA_TYPE_SHORT) {
		for (int i = 0; i < N; i++) {
			short re0 = (short) m_data_real[i];
			write_short(outf, (short) re0);
		}
	} else if (m_data_type == MDA_TYPE_UINT16) {
		for (int i = 0; i < N; i++) {
			unsigned short re0 = (unsigned short) m_data_real[i];
			write_unsigned_short(outf, (unsigned short) re0);
		}
	} else if (m_data_type == MDA_TYPE_INT32) {
		for (int i = 0; i < N; i++) {
			int re0 = (int) m_data_real[i];
			write_int(outf, re0);
		}
	}
	else {
		printf ("Problem in do_write... unexpected data type: %d\n",m_data_type);
		return false;
	}
	return true;
}
void MdaPrivate::write_int(FILE *outf,int val) {
	fwrite(&val,4,1,outf);
}
void MdaPrivate::write_float(FILE *outf,float val) {
	fwrite(&val,4,1,outf);
}
void MdaPrivate::write_short(FILE *outf,short val) {
	fwrite(&val,2,1,outf);
}
void MdaPrivate::write_unsigned_short(FILE *outf,unsigned short val) {
	fwrite(&val,2,1,outf);
}
void MdaPrivate::write_unsigned_char(FILE *outf,unsigned char val) {
	fwrite(&val,1,1,outf);
}
