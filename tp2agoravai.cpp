#include "mpi.h"
#include <stdio.h>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;



int *escolhido;
struct cidades{
	int id;
	float x;
	float y;

};

float calculaDist(cidades c1, cidades c2)
{
	float aux;
	aux = sqrt(pow(c1.x-c2.x,2)+pow(c1.y-c2.y,2));
	return aux;
}

void swap (int *x, int *y)
{
    int temp;
    temp = *x;
    *x = *y;
    *y = temp;
}
float copy_array(int *a, int n, float **matDist,float cost)
{
    int i, sum = 0;
    //for(int i=0;i<=n;i++)printf("%d\n",a[i]);
    for(i = 0; i <= n; i++)
    {
        sum += matDist[a[i % (n+1)]][a[(i + 1) % (n+1)]];
    }
    if (cost > sum)
    {
        cost = sum;   
    }
    //printf("\n\n\n");
    return cost;
}  

float permute(int *a, int i, int n,float **mD,float cost) 
{
   int j, k; 
   if (i == n)  
   {
       cost = copy_array(a, n, mD,cost);
   }
   else
   {
        for (j = i; j <= n; j++)
        { 
            swap((a + i), (a + j));
            cost = permute(a, i + 1, n,mD,cost);
            swap((a + i), (a + j));
        }
    }
    return cost;
} 

int main (int argc,char* argv[])
{
    int numtask, rank, dest, source, rc, count, tag = 1;
    MPI_Status Stat;
    MPI_Request request;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtask);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank ==0){
    	string nav;
    	int tam_city;
    	cidades *pc;
    	int *IndCity;
        float **matDist;
        int aux;
        float cost = 1000000,cost2=1000000,*comparaCost;
        int fatia,fatiaAux;
    	ifstream ifs;

    	ifs.open(argv[1]);
    	if (ifs.good()){
    		while(nav!="DIMENSION:"){
    			ifs >> nav;
    		}
    		ifs>>tam_city;
    		while (nav != "NODE_COORD_SECTION") {
                //printf("%s\n", nav.c_str());
                ifs >> nav;
            }
            fatia = tam_city/numtask;
            escolhido = new int [tam_city];
            pc = new cidades[tam_city];
            matDist = new float*[tam_city];
            IndCity = new int[tam_city];
            comparaCost = new float[numtask];
            for (int i =0;i<tam_city;i++){
            	ifs >> pc[i].id;
            	ifs >> pc[i].x;
            	ifs >> pc[i].y;
            	matDist[i] = new float[tam_city];
            	IndCity[i] = i;
            }
            for (int i=0;i<tam_city;i++)
            {
            	for(int j =0;j<tam_city;j++)
            	{
            		matDist[i][j] = calculaDist(pc[i],pc[j]);
            		printf("Cidade %d para Cidade %d: %f\n",i+1,j+1,matDist[i][j]);
            	}
            }
    	}
        for (int i =1;i<numtask;i++){
            rc = MPI_Send(&tam_city,1,MPI_INT,i,i,MPI_COMM_WORLD);
            for(int j=0;j<tam_city;j++)
                rc = MPI_Send(matDist[j],tam_city,MPI_INT,i,i+1,MPI_COMM_WORLD);
        }
        for(int i=0;i<fatia;i++)
        {
            aux = (rank+1)+(i*numtask);
            if(aux<tam_city){
                swap(IndCity[1],IndCity[aux]);
                cost2 = permute(IndCity,2,tam_city-1,matDist,cost);
                if(cost2<cost)cost=cost2;
                swap(IndCity[aux],IndCity[1]);    
            }
        }

        comparaCost[0]=cost;
        for(int i =1;i<numtask;i++)
        {
           rc = MPI_Recv(&comparaCost[i],1,MPI_FLOAT,i,i+3,MPI_COMM_WORLD,&Stat);
        }
        for(int i =1;i<numtask;i++)
        {
            if(comparaCost[i]<cost)
                cost = comparaCost[i];
        }
    	printf("Custo minimo %f\n",cost);
    }else{
        int tamC,*cities,fatia,auxInd;
        float **distCities,costRV=1000000,costAux=1000000;
        rc = MPI_Recv(&tamC,1,MPI_INT,0,rank,MPI_COMM_WORLD, &Stat);
        cities = new int[tamC];
        distCities = new float *[tamC];
        for (int i=0;i<tamC;i++){
            distCities[i] = new float[tamC];
            cities[i] = i;
        }
        for(int i=0;i<tamC;i++)
            rc = MPI_Recv(distCities[i],tamC,MPI_INT,0,rank+1,MPI_COMM_WORLD,&Stat);

        for(int i=0;i<tamC;i++)
        {
            for(int j=0;j<tamC;j++)
            {
                printf("Cidade %d para a cidade %d é igual a: %f\n",i,i+1,distCities[i][j]);
            }
        }
        fatia = tamC-1/numtask;
        for(int i=0;i<fatia;i++)
        {   
            auxInd = (rank+1)+(i*numtask);
            if(auxInd<tamC){
                swap(cities[1],cities[auxInd]);
                costAux = permute(cities,2,tamC-1,distCities,costRV);
                if(costAux<costRV)costRV=costAux;
                swap(cities[auxInd],cities[1]);    
            }
        }
        printf("O menor do rank %d eh %f\n",rank,costRV);
        rc = MPI_Send(&costRV,1,MPI_FLOAT,0,rank+3,MPI_COMM_WORLD);
    }


    MPI_Finalize();
}

