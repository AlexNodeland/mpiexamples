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
system("hostname");
cout<<"Hello from processor "<<my_rank<<" of "<<p<<" processes."<<endl;

MPI::Finalize();
return 0;
}

