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

double total_time = 0;

int mindots = 3*atoi(argv[1]);
int num_dots;
double ratz, my_r, vol, ratio;
int* my_in = new int;
int* in = new int;

int R = 6;
int r = 3;
int h = 5;

MPI::Init(argc, argv);
MPI_Status status;
MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
MPI_Comm_size(MPI_COMM_WORLD, &p);

if(mindots%p == 0)
{
        num_dots = mindots;
}
else
{
        num_dots = mindots + p - mindots%p;
}

float* dots = new float[num_dots];

float* my_dots = new float[num_dots/p];


if(my_rank==0)
{

        random_device rd;
        mt19937 mt(rd());
        uniform_real_distribution<float> dist(0,R);

        for(int i=0; i<num_dots; i++)
                dots[i] = dist(mt);

}

total_time -= MPI_Wtime();

MPI_Scatter(dots, num_dots/p, MPI::FLOAT, my_dots, num_dots/p, MPI::FLOAT, 0, MPI_COMM_WORLD);

for(int i = 0; i < num_dots/p; i+=3)
{
        if (my_dots[i] > 5)
                ratz = 0;
        else
        {
                ratz = R - ((R-r)/h)*my_dots[i];
                my_r = sqrt(pow(my_dots[i+1],2)+pow(my_dots[i+2],2));
        }
        if(my_r<ratz)
                my_in++;
}
cout<<*my_in<<endl;
MPI_Reduce(&my_in,&in,1,MPI::INT,MPI_SUM,0,MPI_COMM_WORLD);

total_time += MPI_Wtime();
if(my_rank==0)
{
cout<<num_dots<<endl;
        ratio = 3*(*in)/num_dots;
        cout<<"found in"<<in<<"ratio"<<ratio<<"rank"<<my_rank<<endl;
        vol = pow(R,3)*ratio;
        cout<<"found vol"<<endl;

        cout<<"Total time: "<<total_time<<endl;
        cout<<"Volume: "<<vol<<endl;
}
cout<<"did if"<<my_rank<<endl;
delete [] dots;
delete [] my_dots;
MPI::Finalize();
return 0;


}

