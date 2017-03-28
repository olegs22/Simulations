#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


double *corr_2p(int len_DD, double *pos_x, double *pos_y, double *pos_z, int limit, int bin){

  int i,j,k,l,r;
  double RR[len_DD][3];
  double DD[len_DD][3];
  double r_data[3];
  double DD_r, RR_r, DR_r;
  double coords[3];
  //double range = RAND_MAX/limit;
  double max_d = sqrt(3.0)*limit;
  double bin_width = (double)(max_d/bin);
  double hist0[bin];
  double hist1[bin];
  double hist2[bin];
  static double corr[10];

  srand(12345);

  for (i=0; i<len_DD; i++){
    /*
    coords[0] = rand() / range;
    coords[1] = rand() / range;
    coords[2] = rand() / range;
    */
    coords[0] = drand48()*limit;
    coords[1] = drand48()*limit;
    coords[2] = drand48()*limit;

    r_data[0] = pos_x[i];
    r_data[1] = pos_y[i];
    r_data[2] = pos_z[i];

    for (j=0; j<3; j++){
      RR[i][j] = coords[j];
      DD[i][j] = r_data[j];
    }
  }


  for (i=0; i<len_DD; i++){
    for (k=i+1; k<len_DD; k++){

      DD_r = sqrt(pow(DD[k][0]-DD[i][0],2) + pow(DD[k][1]-DD[i][1],2) + pow(DD[k][2]-DD[i][2],2));
      RR_r = sqrt(pow(RR[k][0]-RR[i][0],2) + pow(RR[k][1]-RR[i][1],2) + pow(RR[k][2]-RR[i][2],2));
      DR_r = sqrt(pow(DD[k][0]-RR[i][0],2) + pow(DD[k][1]-RR[i][1],2) + pow(DD[k][2]-RR[i][2],2));
    }

    for (l=0; l<bin; l++){
      if(DD_r>=bin_width*l && DD_r<bin_width*(l+1)){
        hist0[l] += 1;
      }
      if(RR_r>=bin_width*l && RR_r<bin_width*(l+1)){
        hist1[l] += 1;
      }
      if(DR_r>=bin_width*l && DR_r<bin_width*(l+1)){
        hist2[l] += 1;
      }
    }


  }


  for (r=0; r<bin; r++){
    if (hist1[r] == 0.0){
      corr[r] = 0.0;
    }
    else{
      corr[r] = (hist0[r] - 2.0*hist2[r] + hist1[r]) / hist1[r];
    }

  }

  return corr;

}
