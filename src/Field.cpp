/*
 * generate.cppx
 *
 *  Created on: Dec 12, 2013
 *      Author: navid
 */

#include <iostream>
#include <cstdio>
#include "string.h"
#include <stdlib.h>
#include <cmath>
#include <fstream>

#include <sstream>

using namespace std;


double distance_squared(double x1, double y1, double x2, double y2) {
	return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
}


double distance2d(double x1, double y1, double x2, double y2) {
    return sqrt(distance_squared(x1,y1,x2,y2));
}



double Green(double x1, double y1, double x2, double y2){
    double r2 = distance_squared(x1,y1,x2,y2);
    return exp(-r2*10);
}


/*instead of caluculating dx dy each time in weight2d,
we first calculate the matrix of all distances once and then access it to construct all the weights.

*/
void dists(double* X, double* Y, double* Dis2, int N) {
  /*calculate all distances once, to reference later
  */
  for (int i = 0; i < N - 1; i++) {
		for (int j = 0; j < N -1; j++) {
		  Dis2[i+j*N]=distance_squared(X[i],Y[i] , X[j], Y[j]);
		}
  }
}

void Green1(double *x, double *y, double *Dist2, double *Grs , int N){
    /*
     Calculates Green's function weights and returns them in Grs
     */
    double Dt;
    Dt=.01; //diffusion * time
    for (int i=0; i< N*N; i++){
      Grs[i]=exp(-Dist2[i]/(4*Dt));
    }
    //return exp(-r2*10);
}


double weight_hub(double x, double y){
    double r2 = distance_squared(x,y,0.0,0.0);
    
    // return r^(-4)
    return 1.0/r2/r2;
}




double score2d(double xi, double yi, double xj,double yj, double *X_hub, double *Y_hub, int N_hub){
    /*
    input:
        xi,yi,xj,yj: coordinates of two nodes
        X_hub, Y_hub: pointer to the array of hub coordinates
        N_hub: number of hubs
    */

    double x,y;
    double score = 0;
    // Loop over all hubs
    for (int counter_hub = 0; counter_hub < N_hub; counter_hub ++ ){
        x = X_hub[counter_hub];
        y = Y_hub[counter_hub];
        score += Green(xi,yi,x,y) * Green(xj,yj,x,y) * weight_hub(x,y);
    }
    return score;
}
//double score2d1(double *x, double *y, double *X_hub, double *Y_hub, int N_hub, int N){
    /*
    input:
        x,y: coordinates of nodes
        X_hub, Y_hub: pointer to the array of hub coordinates
        N_hub: number of hubs
    */
/*
    double x1,y1;
    double score = 0;
    // Loop over all hubs
    for (int counter_hub = 0; counter_hub < N_hub; counter_hub ++ ){
        x1 = X_hub[counter_hub];
        y1 = Y_hub[counter_hub];
        score += Green(xi,yi,x,y) * Green(xj,yj,x,y) * weight_hub(x,y);
    }
    return score;
}
*/



void get_point(double &x, double &y) {
	double r = 1.0 * rand() / RAND_MAX;
	double q = 1.0 * rand() / RAND_MAX;

	x = r  - 0.5;
	y = q  - 0.5;

}


void generate_hubs(double* X, double* Y, int N_hub){
    for (int i = 0 ; i < N_hub; i++)
        get_point(X[i],Y[i]);
}

//OLD
double weight2d(double x1, double y1, double x2, double y2) {
	double d = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
	double w = exp(-10.0 * d);
	return w;

}

void simulate2d(int id, int N, int N_hub) {
  /*
  Input:
      id: the integer id of the current batch
      N: number of nodes
      alpha: blah blah
      N_hub: number of hubs
  */


  // X coordinates of the nodes
  double X[N];

  // y coordinates of the nodes
  double Y[N];

  double degree[N];

  
  // X coordinates of the hubs
  double X_hub[N_hub];

  // y coordinates of the hubs
  double Y_hub[N_hub];

  
  // generate hubs
  generate_hubs(X_hub, Y_hub, N_hub);


	ofstream f;


	stringstream filename;
	filename << "degrees-" <<  N << "-" << N_hub << ".txt";


//	const char* filename = "degrees.txt";
	cout << filename.str().c_str() << "\n";
	if (id == 0)
		f.open(filename.str().c_str());
	else
		f.open(filename.str().c_str(), std::fstream::app);

	for (int i = 0; i < N; i++)
		degree[i] = 0;

	for (int i = 0; i < N; i++) {
		get_point(X[i], Y[i]);
	}
	printf("Done getting points.\n");
	double s;
	double *dist2;
	dist2=(double *)malloc(N*N*sizeof(double));
	double *grs;
	grs=(double *)malloc(N*N*sizeof(double));
	//for(int j=1;j<N;j++)
    //dist2[j]=dist2[j-1]+N;
	printf("calculating distance^2 matrix...\n");
	dists(X,Y,dist2, N);
	printf("distances computed.\n");
	Green1(X,Y,dist2,grs,N);
	printf("Greens computed computed.\n");
	
	for (int i = 0; i < N - 1; i++) {
		for (int j = i + 1; j < N; j++) {
			s = score2d(X[i],Y[i],X[j],Y[j], X_hub, Y_hub, N_hub);
			// doing j>i makes sure each pai is counted only once. each side gets this score added to their degree.
			degree[i] = degree[i] + s ;
			degree[j] = degree[j] + s ;
		}
	}

	for (int i = 0; i < N; i++) {
		if (X[i] * X[i] + Y[i] * Y[i] < 0.1)
			f << degree[i] << " " << X[i] << " " << Y[i] << '\n';
	}

	f.close();
}
//





void runBatch2D(int num_repeats, int N, int N_hub) {
	for (int i = 0; i < num_repeats; i++){
        cout<< "Simulating batch no "<< i << "\n";
		    simulate2d(i, N,  N_hub);
    }
        

}

int main(int argc, char *argv[]) {
	/* default values */
	int N = 1200;
    int N_hub = 100;

	/* process commandline arguments */
	if (argc >= 3) {
		N = atoi(argv[1]);
		N_hub = atoi(argv[2]);
		
	}
	int batch=10;
	if (argc==4){
	  batch= atoi(argv[3]);
	}

	runBatch2D(batch, N, N_hub);

}
