#include <stdio.h>
#include <math.h>
#include <stdlib.h>


double * corr_2p(int len_DD, double *pos_x, double *pos_y, double *pos_z, int limit, int bin){

  int i,j,k,l,r;
  double RR[len_DD][3];
  double DD[len_DD][3];
  double r_data[3];
  double DD_r, RR_r, DR_r;
  double coords[3];
  double range = RAND_MAX/limit;
  double max_d = sqrt(2.0*pow(limit,2.0));
  double bin_width = (double)(max_d/bin);
  double hist0[bin];
  double hist1[bin];
  double hist2[bin];
  static double corr[10];

  for (i=0; i<len_DD; i++){
    coords[0] = rand() / range;
    coords[1] = rand() / range;
    coords[2] = rand() / range;

    r_data[0] = pos_x[i];
    r_data[1] = pos_y[i];
    r_data[2] = pos_z[i];

    for (j=0; j<3; j++){
      RR[i][j] = coords[j];
      DD[i][j] = r_data[j];
    }
  }


  for (i=0; i<len_DD; i++){
    for (k=1; k<len_DD; k++){

      DD_r = sqrt(pow(DD[i][0]-DD[k][0],2.0) + pow(DD[i][1]-DD[k][1],2.0) + pow(DD[i][2]-DD[k][2],2.0));
      RR_r = sqrt(pow(RR[i][0]-RR[k][0],2.0) + pow(RR[i][1]-RR[k][1],2.0) + pow(RR[i][2]-RR[k][2],2.0));
      DR_r = sqrt(pow(DD[i][0]-RR[k][0],2.0) + pow(DD[i][1]-RR[k][1],2.0) + pow(DD[i][2]-RR[k][2],2.0));
    }

    for (l=0; l<bin; l++){
      if(DD_r>=(double)(l*bin_width) && DD_r<(double)((l+1)*bin_width)){
        hist0[l] += 1;
      }
      if(RR_r>=(double)(l*bin_width) && RR_r<(double)((l+1)*bin_width)){
        hist1[l] += 1;
      }
      if(DR_r>=(double)(l*bin_width) && DR_r<(double)((l+1)*bin_width)){
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

int main(){
  FILE *file;
  //int len = 3102;
  int len = 10000;
  double box_size = 50.0;
  int bining = 10;
  double *hist;
  double x[len];
  double y[len];
  double z[len];
  double pos1[len], pos2[len], pos3[len];
  int i,k,j;

  //file = fopen("/Users/Oleg/Documents/Simulations/positions.txt","r");
  file = fopen("/Users/Oleg/Documents/Simulations/random_data.txt","r");


  //while (!feof(file)) {
  for(k=0; k<len; k++){
    /* code */
    fscanf(file,"%lf\t%lf\t%lf",&x[k],&y[k],&z[k]);
  }

  hist = corr_2p(len+1, x, y, z, box_size, bining);
  for (i=0; i<bining; i++){
    printf("*(hist + %d) : %lf\n", i, *(hist+i));
  }
  fclose(file);
  return 0;
}
