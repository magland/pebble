#include "do_pca.h"
#include <math.h>
#include <stdlib.h>

void define_random_vector(long N,double *v) {
    for (long j=0; j<N; j++) {
        v[j]=(rand()*1.0/RAND_MAX)*2-1;
    }
}
double vector_norm(long N,double *v) {
    double norm=0;
    for (long j=0; j<N; j++) norm+=v[j]*v[j];
    norm=sqrt(norm);
    return norm;
}
void vector_normalize(long N,double *v) {
    double norm=vector_norm(N,v);
    if (norm==0) return;
    for (long j=0; j<N; j++) v[j]/=norm;
}
double inner_product(long N,double *v,double *w) {
    double ret=0;
    for (long j=0; j<N; j++)
        ret+=v[j]*w[j];
    return ret;
}
void subtract_component(long N,double *v,double *comp) {
    double ip=inner_product(N,v,comp);
    for (long j=0; j<N; j++) {
        v[j]-=ip*comp[j];
    }
}


void do_pca(int M,int N,int C,double *out,double *in) {
    int num_iterations=10; //not sure how to set

    int MN=M*N;
    int MC=M*C;
    int CN=C*N;


    //Define the working vectors
    double *working_vectors=(double *)malloc(sizeof(double)*M*N);
    for (int ii=0; ii<MN; ii++) working_vectors[ii]=in[ii];

    //Define the working components
    double *working_components=(double *)malloc(sizeof(double)*MC);
    for (int ii=0; ii<MC; ii++) working_components[ii]=0;

    //allocate the coefficients and energies
    double *coeffs=(double *)malloc(sizeof(double)*CN);
    double *energies=(double *)malloc(sizeof(double)*C);

    for (int cc=0; cc<C; cc++) {
        double *component_vector=&working_components[cc*M];
        double component_norm=vector_norm(M,component_vector);

        if (component_norm<0.1) {
            define_random_vector(M,component_vector);
            vector_normalize(M,component_vector);
        }

        for (long it=0; it<num_iterations; it++) {
            double *hold=(double *)malloc(sizeof(double)*M);
            for (int j=0; j<M; j++) hold[j]=0;
            for (int n=0; n<N; n++) {
                double *ptr=&working_vectors[n*M];
                double ip=inner_product(M,ptr,component_vector);
                for (int j=0; j<M; j++) hold[j]+=ip*ptr[j];
            }
            vector_normalize(M,hold);
            for (int j=0; j<M; j++) component_vector[j]=hold[j];
            free(hold);
        }

        //Compute coefficients
        for (long n=0; n<N; n++) {
            double ip0=inner_product(M,&working_vectors[n*M],&working_components[cc*M]);
            coeffs[cc+C*n]=ip0;
        }
        //Compute energy (lambda)
        double val=0;
        for (int n=0; n<N; n++) {
            val+=coeffs[cc+C*n]*coeffs[cc+C*n];
        }
        energies[cc]=val;
        //Subtract this component from the working vectors
        for (int n=0; n<N; n++) {
            subtract_component(M,&working_vectors[n*M],&working_components[cc*M]);
        }
    }

    for (int ii=0; ii<CN; ii++) out[ii]=coeffs[ii];

    //free working components, working vectors, coefficients, and energies
    free(working_components);
    free(working_vectors);
    free(coeffs);
    free(energies);
}
