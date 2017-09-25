#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include <cblas.h>
#include <sys/time.h>

double getETime(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + (double)tv.tv_usec*1e-6;
}

int main (int argc, char * argv[]){
    FILE * fpr;
    FILE * ofpr;
    char * fname = argv[1];
    long i, j, size, N, ret;
    double _s, s, _e, e;

    printf("DataAll\n");

    ofpr = fopen(argv[2],"w");
    _s = getETime();
    s = getETime();

    fpr = fopen(fname, "rb" );
    int m = 0;
    ret = fread( &m, sizeof(int), 1, fpr );
    printf("%d\n", m);
    N = (long)m;

    int * D;
    D = malloc((long)sizeof(int)*(N*N));
    ret = fread( D, sizeof(int), N*N, fpr );

    fclose( fpr );
    e = getETime();

    fprintf(ofpr, "Load Time: %lf\n", (e-s));

    double * D2;
    D2 = malloc((long)sizeof(double)*(N)*(N));

    s = getETime();

    int edge = 0;

    for(i = 0; i < N; i++){
        D2[i+i*N] = 0.0;
        for(j = i+1; j < N; j++){
            D2[i+j*N] = D2[j+i*N] = (double)(D[i+j*N]*D[i+j*N]);
            if(D[i+j*N] == 1) edge++;
        }
    }

    printf("%d\n", edge);

    for(i = 0; i < N; i++){
        for(j = i+1; j < N; j++){
            if(D[i+j*N] == 1) printf("%d %d\n", i, j);
        }
    }

    e = getETime();
    fprintf(ofpr, "Edge Output Time: %lf\n", (e-s));
    
    free(D);
    
    s = getETime();
    double * b; double * H;
    b = malloc((long)sizeof(double)*(N)*(N));
    H = malloc((long)sizeof(double)*(N)*(N));

    double c = -1.0/(double)N;
    for(i = 0; i < N; i++){
        b[i+i*N] =  0.0;
        H[i+i*N] = 1.0 + c;
        for(j = i+1; j < N; j++){
            b[i+j*N] = b[j+i*N] = 0.0; 
            H[i+j*N] = H[j+i*N] = c;
        }
    }
    double alpha = 1.0, beta = 0;
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
            N, N, N, alpha, H, N, D2, N, beta, b, N);
    free(D2);
    double * B;
    B = malloc((long)sizeof(double)*(N)*(N));
    for(i = 0; i < N; i++){
        B[i+i*N] = 0.0;
        for(j = i+1; j < N; j++){
            B[i+j*N] = B[j+i*N] = 0.0; 
        }
    }
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
            N, N, N, alpha, b, N, H, N, beta, B, N);
    free(b); free(H);

    double * A;
    double * w;
    A = malloc((long)sizeof(double) * (N) * (N));
    w = malloc((long)sizeof(double) * (N));

    for(i = 0; i < N; i++){
        A[i+i*N] = -0.5*B[i+i*N];
        for(j = i+1; j < N; j++){
            A[i+j*N] = -0.5*B[i+j*N];
            A[j+i*N] = -0.5*B[j+i*N];
        }
    }
    free(B);
    
    e = getETime();
    fprintf(ofpr, "Set Matrix Time: %lf\n", (e-s));

    s = getETime();
    lapack_int info, n, lda;
    n = N; lda = n;

    for(i = 0; i < n; i++){
        w[i] = 0;
    }

    info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', n, A, n, w);

    e = getETime();

    fprintf(ofpr, "Lapack Time: %lf\n", (e-s));

    s = getETime();

    int dim = 0;
    double eps = 0.01;

    for(i = N-1; i >= 0; i--){
        if(w[i] > eps) dim++;
        else break;
    }

    printf("%d\n", dim);

    double * lambdas;
    lambdas = malloc((long)sizeof(double)*dim);

    for(i = 0; i < dim; i++){
        lambdas[i] = w[N-1-i];
    }

    free(w); 

    double * U;
    U = malloc((long)sizeof(double)*N*N);

    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            U[i+j*N] = A[j+(N-1-i)*N];
        }
    }

    free(A);

    double * P; double * L; 
    P = malloc((long)sizeof(double)*N*dim);
    L = malloc((long)sizeof(double)*N*dim);

    for(i = 0; i < N; i++){
        for(j = 0; j < dim; j++){
            if(i == j) L[i+j*N] = (double)sqrt(lambdas[j]);
            else L[i+j*N] = 0.0;
            P[i+j*N] = 0.0;
        }
    }

    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
            N, dim, N, alpha, U, N, L, N, beta, P, N);

    e = getETime();
    fprintf(ofpr, "Projection Time: %lf\n", (e-s));

    s = getETime();
    for(i = 0; i < dim; i++){
        printf("%f\n", lambdas[i]);
    }

    for(i = 0; i < N; i++){
        for(j = 0; j < dim-1; j++){
            printf("%f ", P[i+N*j]);
        }
        printf("%f\n", P[i+N*(dim-1)]);
    }

    e = getETime();
    fprintf(ofpr, "OUTPUT Time: %lf\n", (e-s));
    
    _e = getETime();
    fprintf(ofpr, "Total Time: %lf\n", (_e-_s));

    fclose(ofpr);
    free(lambdas); free(U); free(P); free(L);
    return 0;
}
