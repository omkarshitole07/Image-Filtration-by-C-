ARCHITECTURE HYPERCUBE(7); /*128 processors*/
#include <mpi.h>

#define n 1280 /* dimension of image */
#define m 10 /* number of rows in each partition */

int inrows[m+2][n+2], outrows[m][n];
int filter[3][3] = {1,1,1,1,1,1,1,1,1};

MPI_Status status;

int myrank, totalproc, blocksize;

void input_image( ) { /*input and partition the image*/
    int image[n+2][n+2];
    int k, dest, tag;
    /* input pixel values for image */
    /* send partition to each process */
    for (k = 0; k < totalproc; k++) {
        dest = k; tag = 1;
        MPI_Send(&image[k*m][0], blocksize, MPI_INT, dest, tag, MPI_COMM_WORLD);
    }
}

main( ) {
    MPI_Init( );
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &totalproc);
    int source, tag, dest, i, j, k;
    blocksize = (m+2)*(n+2);
    if (myrank == 0) 
        input_image( );
    /* receive my block of image*/
    source = 0; 
    tag = 1;
    MPI_Recv(&inrows[0][0], blocksize, MPI_INT, source, tag, MPI_COMM_WORLD, &status);

    /* apply filter mask to my partition*/
    for (i = 1; i <= m; i++)
        for (j = 1; j <= n; j++)
            outrows[i-1,j-1] =
                (int) ((filter[0,0]*inrows[i-1,j-1]
                + filter[0,1]*inrows[i-1,j]
                + filter[0,2]*inrows[i-1,j+1]
                + filter[1,0]*inrows[i,j-1]
                + filter[1,1]*inrows[i,j]
                + filter[1,2]*inrows[i,j+1]
                + filter[2,0]*inrows[i+1,j-1]
                + filter[2,1]*inrows[i+1,j]
                + filter[2,2]*inrows[i+1,j+1]) / 9);
    /* send filtered partition back to process 0 */
    blocksize = m*n;
    dest = 0; 
    tag = 2;
    MPI_Send(&outrows[0][0], blocksize, MPI_INT, dest, tag, MPI_COMM_WORLD);
    /*assemble partitions into filtered image as output*/
    if (myrank == 0) {
        for (k = 0; k < totalproc; k++) {
            source = k; tag = 2;
            MPI_Recv(&outrows[0][0], blocksize, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
            for (i = 0; i < m; i++) {
                for (j = 0; j <= n; j++)
                    cout << outrows[i][j];
                cout << endl;
            }
        }
    }
    MPI_Finalize( );
}