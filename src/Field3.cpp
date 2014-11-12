/*
 * //generate.cppx
 *  compile from root directory with
 *  g++ -o ./bin/field.exe ./src/Field3.cpp
 *  Created on: Dec 12, 2013
 *      Author: Navid Dianati
 *
 *  Modified on: Nov 1, 2014
 *           by: Nima Dehmamy
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
void dists_to_hubs(double* X, double* Y,double* Xh, double* Yh, double* Dis2, int N, int N_hub) {
  /*calculate all distances once, to reference later
  */
  for (int i = 0; i < N ; i++) {
		for (int j = 0; j < N_hub ; j++) {
		  Dis2[i*N_hub+j]=distance_squared(X[i],Y[i] , Xh[j], Yh[j]);
		}
  }
}

void Green1(double *Dist2, double *Grs , int N, int N_hub, double Dt){
    /*
     Calculates Green's function weights and returns them in Grs
     */
    //double Dt;
    //Dt=.01; //diffusion * time
    for (int i=0; i< N*N_hub; i++){
      Grs[i]=exp(-Dist2[i]/(4*Dt));
    }
    
}


double weight_hub(double x, double y){
    double r2 = distance_squared(x,y,0.0,0.0);
    
    // return r^(-4)
    return 1.0/r2/r2;
    //double r2t= r2 + t;
    //return - 4*t*t*(r2-t)/ pow(r2t,3) - t*t/pow(r2t,2) +2*t/ r2t;
}

void weight_hub1(double *Xh, double *Yh, double *weights, int N_hub, double t){
  /* weights is an array of length N_hub which stores the weight of each hub.
  */
  for (int i=0; i<N_hub; i++){
    double r2 = distance_squared(Xh[i],Yh[i],0.0,0.0);
    double r2t= r2 + t;
    weights[i] = - 4*t*t*(r2-t)/ pow(r2t,3) - t*t/pow(r2t,2) +2*t/ r2t; //=1.0/r2/r2;
  }
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



void simulate2d(int id, int N, int N_hub, double Dt_max, double dt) {
  /*
  Input:
      id: the integer id of the current batch
      N: number of nodes
      N_hub: number of hubs
      Dt_max: Diffusion * time_max
      dt: steps for incrementing Dt to Dt_max
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
  
	for (int i = 0; i < N; i++)
		degree[i] = 0;
		
	for (int i = 0; i < N; i++) {
		get_point(X[i], Y[i]);
	}
	printf("Done getting points.\n");
	double s;
	//****************************
	
	double *dist2;
	dist2=(double *)malloc(N*N_hub*sizeof(double));// for distances to hubs
	double *grs;
	grs=(double *)malloc(N*N_hub*sizeof(double)); // for Green's functions to hubs
	double *wgh;
	wgh=(double *)malloc(N_hub*sizeof(double)); // for weights of hubs
	double *wgh_Gr;
	wgh_Gr=(double *)malloc(N_hub*sizeof(double)); // for weights multiplied by sum of Green's f's to get what needs to be multiplied into Gr to get degree
	
	/* We want to calculate the weights for each time and form the network for each time
	each gives us a degree contribution. The total degree comes from summing over time
	
	The  max time times diffusion constant "Dt" must be << size of space squared.
	Here we discard points with r^2 > 0.1 Therefore max Dt <<.1
	*/
	
	
	printf("calculating distance^2 matrix...\n");
	dists_to_hubs(X,Y, X_hub,Y_hub ,dist2, N, N_hub);
	printf("distances to hubs computed.\n");
	//weight_hub1(X_hub,Y_hub, wgh, N_hub); // weights of hubs put in wgh
	/* Now that we have Gr's and weights we can use them to calculate the degrees
	*/
	
	printf("Weights computed. Multiplying Greens...\n");
	//double Dt=0.0001;
	//double inc=0.001;
	//double Dt_max=0.04;
	
	FILE *ff;
	char fnam[100];
	int fnam_len = sprintf( fnam, "output/degrees-N=%d,N_hub=%d,Dt_max=%.2g,dt=%.2g.txt",N,N_hub,Dt_max,dt);
	
	ff=fopen(fnam,"a") ;
  //fprintf(ff,"%.8g %.8g %.8g\n",degree, X, Y);
  //fflush(ff);
  //fclose(ff);
	
	/*
	ofstream f;
	stringstream filename;
	filename << "output/degrees-" <<  N << "-" << N_hub << ".txt";
  if (id == 0)
		f.open(filename.str().c_str());
	else
		f.open(filename.str().c_str(), std::fstream::app);
  cout << filename.str().c_str() << "\n";
  */
  cout << "output file: " << fnam << "\n";
  char outs[10];
  cout << "Greens for Dt=\n";
	for (double Dt=1.0e-7+dt; Dt< Dt_max; Dt+=dt){
	    weight_hub1(X_hub,Y_hub, wgh, N_hub, Dt); // weights of hubs put in wgh
    	Green1(dist2,grs,N, N_hub, Dt); // takes distances to hubs and puts Green's f's in grs
    	//printf("Greens for Dt=%.4g computed.\n", Dt );
    	int ndt=sprintf(outs,"%.5g, ", Dt);
    	cout << outs << "\n" ;
    	for (int j = 0; j < N_hub; j++) {
    		wgh_Gr[j]=0;
    		for (int i = 0; i < N; i++) {
    		  wgh_Gr[j]+=wgh[j]*grs[i*N_hub+j]; // weights *Greens
    		  // This is the " \lambda * \int Green " which is the same for all degree calculations.
    		}
    	}
    	/* The degrees are then the product of this and another Green's f summed over hubs
    	*/
    	//printf("Greens multiplied. Calculating degrees...\n");
    	for (int i = 0; i < N ; i++) {
    	  for (int j = 0; j < N_hub; j++) {
    	    degree[i]+= grs[i*N_hub+j]*wgh_Gr[j]; // Gr[i --> j] * \lambda[j] * \int Gr[all N --> j]
    	  }
    	}
  }// Dt sum
  cout << "\nDone!\n";

	for (int i = 0; i < N; i++) {
		if (X[i] * X[i] + Y[i] * Y[i] < 0.1)
			//f << degree[i] << " " << X[i] << " " << Y[i] << '\n';
			fprintf(ff,"%.8g %.8g %.8g\n",degree[i], X[i], Y[i]);
	}
  
  fflush(ff);
  fclose(ff);
	//f.close();
	
}
//





void runBatch2D(int num_repeats, int N, int N_hub, double Dt_max, double dt) {
	for (int i = 0; i < num_repeats; i++){
        cout<< "Simulating batch no "<< i << "\n";
		    simulate2d(i, N,  N_hub, Dt_max, dt);
    }
        

}

int main(int argc, char *argv[]) {
	srand (time(NULL)); // seed rand with time
	/* default values */
	
	int N = 1200;
  int N_hub = 100;
  double Dt_max = 5.e-3;
  double dt = 5.e-4;

	/* process commandline arguments */
	if (argc >= 3) {
		N = atoi(argv[1]);
		N_hub = atoi(argv[2]);
		
	}
	int batch=10;
	if (argc >= 4){
	  batch = atoi(argv[3]);
	}
	if (argc >= 5){
	  Dt_max = atof(argv[4]); //!!! float, thus atof()
	}
  if (argc >= 6){
	  dt = atof(argv[5]);
	}
	runBatch2D(batch, N, N_hub, Dt_max, dt);

}
