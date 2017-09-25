#include <stdlib.h>
#include <stdio.h>

int main(int argc, char * args[]){
    FILE *fp;
    FILE *fpw;
    char *fname = args[1];
    char *fname_w = args[2];
    long i, j, size, N; int ret;
     
    fp = fopen(fname, "r");

    ret = fscanf(fp, "%ld", &N);

    long size_N = (long)N;

    int * D;
    size = (long)sizeof(int)*((size_N)*(size_N)+1);
    D = malloc(size);
    
    for(i = 0; i < N*N+1; i++){
        D[i] = 0;
    }
    
    D[0] = N;

    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            ret = fscanf(fp, "%d", &D[i+N*j+1]);
        }
    }

    fpw = fopen( fname_w, "wb" );

    fwrite( D, sizeof( int ), size, fpw );

    fclose( fp );
    fclose( fpw );

    return 0;
}
