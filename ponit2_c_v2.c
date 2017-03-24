#include <stdio.h>
#include <math.h>
#include <stdlib.h>


double * corr_2p(double pos_x[2999], double pos_y[2999], double pos_z[2999], int len_DD, int limit, int bin){

  int i,j,k,l,r;
  double RR[len_DD][3];
  double DD[len_DD][3];
  double r_data[3];
  double DD_r[len_DD], RR_r[len_DD], DR_r[len_DD];
  double x,y,z;
  double coords[3];
  double bin_width;
  double range = RAND_MAX/limit;

  bin_width = limit/bin;
  double hist0[bin];
  double hist1[bin];
  double hist2[bin];
  static double corr[10];

  for (i=0; i<len_DD; i++){
    x = rand() / range;
    y = rand() / range;
    z = rand() / range;
    coords[0] = x;
    coords[1] = y;
    coords[2] = z;

    for (j=0; j<3; j++){
      RR[i][j] = coords[j];
    }
  }

  for (i=0; i<len_DD; i++){
    r_data[0] = pos_x[i];
    r_data[1] = pos_y[i];
    r_data[2] = pos_z[i];

    for (j=0; j<3; j++){
      DD[i][j] = r_data[j];
    }
  }

  for (i=0; i<len_DD; i++){
    for (k=1; k<len_DD; k++){

      DD_r[i] = sqrt(pow(DD[i][0]-DD[k][0],2.0) + pow(DD[i][1]-DD[k][1],2.0) + pow(DD[i][2]-DD[k][2],2.0));
      RR_r[i] = sqrt(pow(RR[i][0]-RR[k][0],2.0) + pow(RR[i][1]-RR[k][1],2.0) + pow(RR[i][2]-RR[k][2],2.0));
      DR_r[i] = sqrt(pow(DD[i][0]-RR[k][0],2.0) + pow(DD[i][1]-RR[k][1],2.0) + pow(DD[i][2]-RR[k][2],2.0));


    }
  }
  for (i=0; i<len_DD; i++){
    for (l=0; l<bin; l++){
      if(DD_r[l]>=(double)(l*bin_width) && DD_r[i]<(double)((l+1)*bin_width)){
        hist0[l] += 1;
      }
      if(RR_r[l]>=(double)(l*bin_width) && RR_r[i]<(double)((l+1)*bin_width)){
        hist1[l] += 1;
      }
      if(DR_r[l]>=(double)(l*bin_width) && DR_r[i]<(double)((l+1)*bin_width)){
        hist2[l] += 1;
      }

    }

  }


  for (r=0; r<bin; r++){
    printf("%lf\n",hist2[r] );
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
  int len = 3000-1;
  double box_size = 70.71068;
  int bining = 10;
  double *hist;
  double x[len];
  double y[len];
  double z[len];
  double pos1[len-1], pos2[len-1], pos3[len-1];
  int i,j;

  //file = fopen("/Users/Oleg/Documents/Simulations/positions.txt","r");
  file = fopen("/Users/Oleg/Documents/Simulations/random_data.txt","r");
  i = 0;

  while (!feof(file)) {
    /* code */
    fscanf(file,"%lf\t%lf\t%lf",&x[i],&y[i],&z[i]);
    i++;
  }
  for (j=0; j<len; j++){
      pos1[j]=x[j];
      pos2[j]=y[j];
      pos3[j]=z[j];
    }

  hist = corr_2p(pos1, pos2, pos3, len-1, box_size, bining);
  for (i=0; i<bining;i++){
    printf("*(hist + %d) : %lf\n", i, *(hist+i));
  }
  fclose(file);
  return 0;
}
