#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


double *corr_ls(int len_DD, double *pos_x, double *pos_y, double *pos_z, int limit, int bin, int model){

  int i,k,r,s;
  double DD_r, RR_r, DR_r;
  int n_dd, n_rr, n_dr;
  double dbin = (double)bin;
  double max_d = (double)(sqrt(3.0)*limit);
  double bin_width = (double)(max_d/dbin);
  double *rand_x = (double*)malloc(len_DD *sizeof(double));
  double *rand_y = (double*)malloc(len_DD *sizeof(double));
  double *rand_z = (double*)malloc(len_DD *sizeof(double));
  int *hist0 = (int*)malloc(bin *sizeof(int));
  int *hist1 = (int*)malloc(bin *sizeof(int));
  int *hist2 = (int*)malloc(bin *sizeof(int));
  double *corr = (double*)malloc(dbin *sizeof(double));

  // initialization of all arrays
  for (i=0; i<len_DD; i++){
    rand_x[i] = 0.0;
    rand_y[i] = 0.0;
    rand_z[i] = 0.0;
  }

  for (r=0; r<bin; r++){
    hist0[r] = 0.0;
    hist1[r] = 0.0;
    hist2[r] = 0.0;
    corr[r] = 0.0;
  }
  srand(time(NULL));

  for (i=0; i<len_DD; i++){

    rand_x[i] = drand48()*limit;
    rand_y[i] = drand48()*limit;
    rand_z[i] = drand48()*limit;
  }

  for (i=0; i<len_DD; i++){
    for (k=i+1; k<len_DD; k++){

      DD_r = sqrt(pow(pos_x[i]-pos_x[k],2) + pow(pos_y[i]-pos_y[k],2) + pow(pos_z[i]-pos_z[k],2));
      RR_r = sqrt(pow(rand_x[i]-rand_x[k],2) + pow(rand_y[i]-rand_y[k],2) + pow(rand_z[i]-rand_z[k],2));

      if(DD_r <= 200.0){
      n_dd = (int)(DD_r / bin_width);
      hist0[n_dd] += 1;
      }

      if (RR_r <= 200.0){
        n_rr = (int)(RR_r / bin_width);
        hist1[n_rr] += 1;
      }
    }
    for (s=0; s<len_DD;s++){
      DR_r = sqrt(pow(pos_x[i]-rand_x[s],2) + pow(pos_y[i]-rand_y[s],2) + pow(pos_z[i]-rand_z[s],2));

      if (DR_r <= 200.0){
      n_dr = (int)(DR_r / bin_width);
      hist2[n_dr] += 1;
    }
    }

  }
  free(rand_x);
  free(rand_y);
  free(rand_z);

  FILE *histdr;
  histdr = fopen("hist2.txt","w");
  for (r=0; r<bin; r++){
    fprintf(histdr, "%d\t%d\t%d\n", hist0[r], hist1[r], hist2[r] );
    if (hist1[r] == 0.0){
      corr[r] = 0.0;
    }
    else{
      if (model == 0){
        corr[r] = (double)(hist0[r] - hist2[r] + hist1[r]) /  hist1[r];
      }
      if (model == 1){
        corr[r] = (double)(hist0[r] / hist1[r]) - 1.0;
      }

    }

  }
  fclose(histdr);
  free(hist0);
  free(hist1);
  free(hist2);
  return corr;

}
