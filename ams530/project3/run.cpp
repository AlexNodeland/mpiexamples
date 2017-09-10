#include<cstdlib>
#include<cmath>
#include<cassert>
#include<cfloat>
#include<cstring>
#include<ctime>
#include<map>
#include<random>
#include<iomanip>
#include<iostream>
using namespace std;

#include<mpi.h>

int main(int argc, char* argv[]){

int i, j, k, new_nrow, new_ncol, p, my_rank;

double total_time = 0;

int nrow = atoi(argv[1]);
int ncol = atoi(argv[2]);
int iterations = atoi(argv[3]);


//Initializing MPI//

MPI::Init(argc, argv);
MPI_Status status;
MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
MPI_Comm_size(MPI_COMM_WORLD, &p);

if(nrow%p == 0){
        new_nrow = nrow;
        new_ncol = ncol;
}else{
        new_nrow = nrow + p - nrow%p;
        new_ncol = ncol + p - ncol%p;
}

//Create matrices dynamically//
float** matrix_A = new float*[new_nrow];
float** matrix_B = new float*[new_nrow];

for(i=0; i<new_nrow; i++)
        matrix_A[i] = new float[new_ncol];

for(i=0; i<new_nrow; i++)
        matrix_B[i] = new float[new_ncol];

float* row_matrix = new float[new_nrow*new_ncol];
float* col_matrix = new float[new_nrow*new_ncol];
float* result_mat = new float[new_nrow*new_ncol];
float* recv_buf_A = new float[new_nrow*new_ncol/p];
float* recv_buf_B = new float[new_nrow*new_ncol/p];
float* recv_buf_C = new float[new_nrow*new_ncol/p];
float* recv_buf_t = new float[new_nrow*new_ncol/p];

// Random number generation //
 if (my_rank ==0){
 random_device rd;
 mt19937 mt(rd());
 uniform_real_distribution<float> dist(-0.9999999,1);


//populate  matrices//
for(i=0; i<new_nrow; i++)
        for(j=0; j<new_ncol; j++){
                if(i >= nrow || j >= ncol){
                        matrix_A[i][j] = 0;
                        row_matrix[i*new_nrow+j] = matrix_A[i][j];
        		matrix_B[i][j] = 0;
                        col_matrix[j*new_nrow+i] = matrix_B[i][j]; 
	       }else{
                        matrix_A[i][j] = dist(mt);
                        row_matrix[i*new_nrow+j] = matrix_A[i][j];
			matrix_B[i][j] = dist(mt);
                        col_matrix[j*new_nrow+i] = matrix_B[i][j];
                }
        }
}
for(i=0; i<iterations; i++){

total_time -= MPI_Wtime();

MPI_Scatter(row_matrix, new_nrow*new_ncol/p, MPI::FLOAT, recv_buf_A, new_nrow*new_ncol/p, MPI::FLOAT, 0, MPI_COMM_WORLD);

MPI_Scatter(col_matrix, new_nrow*new_ncol/p, MPI::FLOAT, recv_buf_B, new_nrow*new_ncol/p, MPI::FLOAT, 0, MPI_COMM_WORLD);

for(int roll =0; roll<p; roll++){


 float* c = new float[new_nrow*new_ncol/p];
 for(i=0; i<p; i++){
 	for(j=0; j<p; j++){
        	double temp = 0;
                for(k=0; k<new_nrow; k++){
			temp = temp+recv_buf_A[i*new_nrow+k]*recv_buf_B[j*new_ncol+k];
                }
                        c[i*new_nrow+j] = temp;			
        }

MPI_Gather(c, (new_nrow/p)*(new_ncol/p), MPI::FLOAT, recv_buf_C, new_nrow*new_ncol/p, MPI::FLOAT, 0, MPI_COMM_WORLD);

int c_part;
for(c_part = roll*(new_nrow*new_ncol/p); c_part<(roll+1)*new_nrow*new_ncol/p; c_part++)
        result_mat[c_part] = recv_buf_C[c_part];


if(my_rank%2 ==0){
	MPI_Send(recv_buf_A, (new_nrow/p)*new_ncol/p, MPI::FLOAT, (my_rank+1)%p, 0, MPI_COMM_WORLD);
        MPI_Recv(recv_buf_t, (new_nrow/p)*new_ncol/p, MPI::FLOAT, (my_rank+p-1)%p, 1, MPI_COMM_WORLD, &status);
}else{
        MPI_Recv(recv_buf_t, (new_nrow/p)*new_ncol/p, MPI::FLOAT, (my_rank+p-1)%p, 0, MPI_COMM_WORLD, &status);
        MPI_Send(recv_buf_A, (new_nrow/p)*new_ncol/p, MPI::FLOAT, (my_rank+1)%p, 1, MPI_COMM_WORLD);
}



recv_buf_A = recv_buf_t;

}

total_time += MPI_Wtime();
}
}
if(my_rank == 0)
cout<<"Total time"<<total_time/iterations<<endl;


//Delete Matrices//
for(i=0; i<new_nrow; i++){
        delete [] matrix_A[i];
        delete [] matrix_B[i];
}      
delete [] matrix_A;
delete [] matrix_B;
delete [] row_matrix;
delete [] col_matrix;

MPI::Finalize();
return 0;

}








