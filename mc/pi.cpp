#include<iostream>
#include<cstdlib>
#include<cmath>
#include<cassert>
#include<cfloat>
#include<cstring>
#include<ctime>
#include<map>
#include<random>
#include<iomanip>

using namespace std;

#include"mpi.h"


int main(int argc, char* argv[])
{
int p, my_rank;


MPI::Init(argc, argv);
MPI_Status status;
MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
MPI_Comm_size(MPI_COMM_WORLD, &p);

double total_time = 0;
int min_dots = atoi(argv[1]);
int iterations = atoi(argv[2]);

int num_dots = min_dots+p+min_dots%p;
int num_coords = 2*num_dots;

float* my_in = new float[1];
float* in = new float[1];
my_in[0]=0;
in[0]=0;

float my_pi;
float* my_coords = new float[num_coords/p];

random_device rd;
mt19937 mt(rd());
uniform_real_distribution<float> dist(0,1);

for(int i=0; i<num_coords/p; i++)
        my_coords[i] = dist(mt);

total_time -= MPI_Wtime();

for(int i=0;i<num_coords/p;i+=2)
        if(sqrt(pow(my_coords[i],2)+pow(my_coords[i+1],2))<=1)
		my_in[0]++;

MPI_Reduce(my_in,in,1,MPI::FLOAT,MPI_SUM,0,MPI_COMM_WORLD);

total_time += MPI_Wtime();		


if(my_rank==0)
{
        my_pi = 4*in[0]/(num_dots*iterations);
        cout<<"Total time= "<<total_time<<endl;
        cout<<"pi= "<<my_pi<<endl;
	cout<<"# processors = "<<p<<endl;
}

delete [] my_coords;
MPI::Finalize();
return 0;
}

