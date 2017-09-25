#include <stdio.h>
#include <stdlib.h>

int main(int argc, char * args[]){
    FILE * fpr;
    char * fname = args[1];
    long i, j, size, size2, N;

    fpr = fopen(fname, "rb" );
    int n = 0;
    size = fread( &n, sizeof(int), 1, fpr );
    printf("%d\n", n);
    N = (long)n;

    int * D;
    D = malloc((long)sizeof(int)*(N*N));
    
    size2 = fread( D, sizeof(int), N*N, fpr );
    for(i = 0; i < N; i++){
        for(j = 0; j < N-1; j++){
            printf("%d ", D[i+j*N]);
        }
        printf("%d\n", D[i+(N-1)*N]);
    }
    fclose( fpr );

    return 0;
}
