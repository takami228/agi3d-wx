#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda.h>
#include <cublas.h>
#include <cuda_runtime_api.h>
#include "magma.h"
#include "magma_lapack.h"
#include "testings.h"
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main(int argc, char * argv[]){
    TESTING_CUDA_INIT();

    cout << "Data" << endl;

    ifstream file(argv[1]);
    int n; file >> n;
    magma_int_t N = n;
    int ** D; double **D2; double **A;
    D = new int*[N]; D2 = new double*[N]; A = new double*[N];

    for(int i = 0; i < N; i++){
        D[i] = new int[N]; D2[i] = new double[N]; A[i] = new double[N];
        for(int j = 0; j < N; j++){
            file >> D[i][j];
        }
    }

    int E = 0;
    vector<int> from, to;

    for(int i = 0; i < N; i++){
        for(int j = i+1; j < N; j++){
            if(D[i][j]==1){
                from.push_back(i);
                to.push_back(j);
                E++;
            }
        }
    }

    cout << N << endl;
    cout << E << endl;

    for(int i = 0; i < E; i++){
        cout << from[i] << " " << to[i] << endl;
    }

    for(int i = 0; i < N; i++){
        D2[i][i] = 0.0f;
        for(int j = i+1; j < N; j++){
            D2[i][j] = (double)(D[i][j]*D[i][j]);
            D2[j][i] = D2[i][j];
        }
    }

    double sum_all = 0.0f;
    double * sum_colums = new double[N];
    double * sum_rows = new double[N];
    for(int i = 0; i < N; i++){
        sum_colums[i] = 0.0f;
        sum_rows[i] = 0.0f;
    }
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            sum_colums[i] += D2[i][j];
            sum_rows[j] += D2[i][j];
            sum_all += D2[i][j];
        }
    }

    double k = 1.0f/(double)n;
    double c = sum_all*k*k;
    for(int i = 0; i < n; i++){
        for(int j = i; j < n; j++){
            double a_ij = (0.5)*(k*sum_colums[i] + k*sum_rows[j] - c - D2[i][j]);
            A[i][j] = a_ij; A[j][i] = a_ij;
        }
    }

    double *h_R, *h_work, *w;
    magma_int_t *iwork;
    int lda = N;
    magma_int_t info;
    const char *uplo = MagmaUpperStr; const char *jobz = MagmaVectorsStr;

    double aux_work[1]; magma_int_t aux_iwork[1];
    magma_dsyevd(jobz[0], uplo[0], N, h_R, N, w, aux_work, -1, aux_iwork, -1, &info);
    magma_int_t lwork, liwork;
    lwork  = (magma_int_t) aux_work[0]; liwork = aux_iwork[0];

    TESTING_MALLOC( w, double, N);
    TESTING_HOSTALLOC( h_R, double, N*N);
    TESTING_HOSTALLOC( h_work, double, lwork);
    TESTING_MALLOC( iwork, magma_int_t, liwork);

    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            h_R[i+j*lda] = A[i][j];
        }
    }
    //lapackf77_ssyevd(jobz, uplo, &N, h_R, &N, w, h_work, &lwork,iwork, &liwork,&info);
    magma_dsyevd( jobz[0], uplo[0], N, h_R, N, w, h_work, lwork, iwork, liwork, &info );

   int dim = 0; double eps = 0.01;

    for(int i = N-1; i >= 0; i--){
        if(w[i] > eps) dim++;
        else break;
    }

    //if(dim > 50) dim = 50;

    double * lambdas = new double[dim];

    cout << dim << endl;

    for(int i = 0; i < dim; i++){
        lambdas[i] = w[N-i-1];
        cout << lambdas[i] << endl;
    }

    //multiply matrix
    magma_int_t _M = N, _N = dim, _K = N;
    magma_int_t Am, An, Bm, Bn;
    magma_int_t ldb, ldc, ldda, lddb, lddc;

    double *h_U, *h_L, *h_P, *h_P2;
    double *d_A, *d_B, *d_C;
    double alpha = 1.0f; double beta  = 0.0f;
    char transA = MagmaTrans; char transB = MagmaNoTrans;

    Am = _M; An = _K; Bm = _K; Bn = _N; lda = _M; ldb = Bm; ldc = _M;

    ldda = ((lda+31)/32)*32; lddb = ((ldb+31)/32)*32; lddc = ((ldc+31)/32)*32;

    TESTING_MALLOC( h_U,  double, lda*_K );
    TESTING_MALLOC( h_L,  double, ldb*Bn );
    TESTING_MALLOC( h_P,  double, ldc*_N );
    TESTING_MALLOC( h_P2, double, ldc*_N );
    TESTING_DEVALLOC( d_A, double, ldda*_K );
    TESTING_DEVALLOC( d_B, double, lddb*Bn );
    TESTING_DEVALLOC( d_C, double, lddc*_N );

    /*
    printf("U=[");
    for(int i = 0; i < N; i++){
        cout << "[";
        for(int j = 0; j < N; j++){
            printf("%f,", h_R[i+j*N]);
        }
        printf("%f];\n", h_R[i+(N-1)*N]);
    }
    cout << "]" << endl;
    */

    /* Initialize the matrices */
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            h_U[i+lda*j] = h_R[j+(N-1-i)*N];
        }
    }


    for(int i = 0; i < N; i++){
        for(int j = 0; j < dim; j++){
            if(i == j) h_L[i+i*ldb] = (double)sqrt(lambdas[j]);
            else h_L[i+j*ldb] = 0.0f;
            h_P[i+ldc*j] = 0.0f;
            h_P2[i+ldc*j] = 0.0f;
        }
    }

    magma_dsetmatrix( Am, An, h_U, lda, d_A, ldda );
    magma_dsetmatrix( Bm, Bn, h_L, ldb, d_B, lddb );
    magma_dsetmatrix( _M, _N, h_P, ldc, d_C, lddc );

    magmablas_dgemm( transA, transB, _M, _N, _K, alpha, d_A, ldda, d_B, lddb, beta, d_C, lddc );
    magma_dgetmatrix( _M, _N, d_C, lddc, h_P2, ldc );

    /*
    printf("A=[");
    for(int i = 0; i < N; i++){
        cout << "[";
        for(int j = 0; j < dim-1; j++){
            printf("%f,", h_P2[i+j*ldc]);
        }
        printf("%f]\n", h_P2[i+(dim-1)*ldc]);
    }
    cout << "]" << endl;
    */


    for(int i = 0; i < N; i++){
        for(int j = 0; j < dim-1; j++){
            printf("%f ", h_P2[i+j*ldc]);
        }
        printf("%f\n", h_P2[i+(dim-1)*ldc]);
    }

    /* Memory clean up */
    TESTING_FREE( w );
    TESTING_HOSTFREE( h_R );
    TESTING_FREE( iwork );
    TESTING_HOSTFREE( h_work );
    TESTING_FREE( h_U );
    TESTING_FREE( h_L );
    TESTING_FREE( h_P );
    TESTING_FREE( h_P2 );
    TESTING_DEVFREE( d_A );
    TESTING_DEVFREE( d_B );
    TESTING_DEVFREE( d_C );
    TESTING_CUDA_FINALIZE();
}
