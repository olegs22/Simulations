#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


double *corr_2p(int len_DD, double *pos_x, double *pos_y, double *pos_z, int limit, int bin, int model){

  int i,j,k,l,r,s;
  double **RR;
  double **DD;
  double r_data[3];
  double DD_r, RR_r, DR_r;
  int n_dd, n_rr, n_dr;
  double coords[3];
  //double range = RAND_MAX/limit;
  double max_d = sqrt(3.0)*limit;
  double bin_width = (double)(max_d/bin);
  double *hist0 = (double*)malloc(bin *sizeof(double));
  double *hist1 = (double*)malloc(bin *sizeof(double));
  double *hist2 = (double*)malloc(bin *sizeof(double));
  double *corr = (double*)malloc(bin *sizeof(double));

  RR = (double **) malloc(len_DD*sizeof(double *));
  DD = (double **) malloc(len_DD*sizeof(double *));
  for (i=0; i<len_DD; i++){

    RR[i] = (double *) malloc(3*sizeof(double));
    DD[i] = (double *) malloc(3*sizeof(double));
  }
  srand(time(NULL));

  for (i=0; i<len_DD; i++){

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

      n_dd = (int)(DD_r / bin_width);
      n_rr = (int)(RR_r / bin_width);
      hist0[n_dd] += 1.0;
      hist1[n_rr] += 1.0;
    }
    for (s=0; s<len_DD;s++){
      DR_r = sqrt(pow(DD[s][0]-RR[i][0],2) + pow(DD[s][1]-RR[i][1],2) + pow(DD[s][2]-RR[i][2],2));

      n_dr = (int)(DR_r / bin_width);
      hist2[n_dr] += 1.0;
    }

  }

  free(DD);
  free(RR);

  for (r=0; r<bin; r++){
    if (hist1[r] == 0.0){
      corr[r] = 0.0;
    }
    else{
      if (model == 0){
        corr[r] = (hist0[r] - hist2[r] + hist1[r]) /  hist1[r];
      }
      if (model == 1){
        corr[r] = (hist0[r] / hist1[r]) - 1.0;
      }

    }

  }
  free(hist0);
  free(hist1);
  free(hist2);
  return corr;

}


double *corr_2p_norm(int len_d, int len_r, double *pos_x, double *pos_y, double * pos_z, double *r_x, double *r_y, double *r_z, int limit, int bin, int model){

  int i,j,k,l,m,r;
  double DD, RR, DR;
  double **data;
  double **random;
  int ndd, nrr, ndr;
  double r_data[3];
  double r_random[3];
  double *hist_dd = (double*)malloc(bin *sizeof(double));
  double *hist_rr = (double*)malloc(bin *sizeof(double));
  double *hist_dr = (double*)malloc(bin *sizeof(double));
  double max_d = sqrt(3.0)*limit;
  double bin_width = (double)(max_d/bin);
  double factor0 = (double) (len_r/len_d) * ((len_r - 1.0)/(len_d - 1.0));
  double factor1 = (double) (len_r - 1.0) / len_d;
  double *corr = (double*)malloc(bin *sizeof(double));

  data = (double **) malloc(len_d*sizeof(double *));
  random = (double **) malloc(len_r*sizeof(double *));
  for (i=0; i<len_d; i++){

    data[i] = (double *) malloc(3*sizeof(double));
  }
  for(l=0;l<len_r;l++){

    random[l] = (double *) malloc(3*sizeof(double));
  }

  for (i=0; i<len_d; i++){

    r_data[0] = pos_x[i];
    r_data[1] = pos_y[i];
    r_data[2] = pos_z[i];

    for (j=0; j<3; j++){
      data[i][j] = r_data[j];
    }
  }
  for(l=0;l<len_r;l++){

    r_random[0] = r_x[l];
    r_random[1] = r_y[l];
    r_random[2] = r_z[l];

    for(j=0;j<3;j++){
      random[l][j] = r_random[j];
    }
  }


  for (i=0; i<len_d; i++){
    for (k=i+1; k<len_d; k++){
      DD = sqrt(pow(data[k][0] - data[i][0],2) + pow(data[k][1]-data[i][1],2) + pow(data[k][2]-data[i][2],2));
      ndd = (int)(DD / bin_width);
      hist_dd[ndd] += 1.0;
    }
  }
  for (l=0; l<len_r;l++){
    for (m=j+1;m<len_r;m++){
      RR = sqrt(pow(random[m][0]-random[l][0],2) + pow(random[m][1]-random[l][1],2) + pow(random[m][2]-random[l][2],2));
      nrr = (int)(RR / bin_width);
      hist_rr[nrr] += 1.0;
    }
  }
  for(i=0;i<len_d;i++){
    for (l=0; l<len_r;l++){
      DR = sqrt(pow(data[i][0]-random[l][0],2) + pow(data[i][1]-random[l][1],2) + pow(data[i][2]-random[l][2],2));
      ndr = (int)(DR / bin_width);
      hist_dr[ndr] += 1.0;
    }

  }

  free(data);
  free(random);

  for (r=0; r<bin; r++){
    if (hist_rr[r] == 0.0){
      corr[r] = 0.0;
    }
    else{
      if (model == 0){
        corr[r] = ((factor0*hist_dd[r]) - (factor1*hist_dr[r]) + hist_rr[r]) / hist_rr[r];
      }
      if (model == 1){
        corr[r] = factor0 * (hist_dd[r] / hist_rr[r]) - 1.0;
      }

    }

  }

  free(hist_dd);
  free(hist_rr);
  free(hist_dr);
  return corr;


}
