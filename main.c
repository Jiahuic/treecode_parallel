/* Jiahui Chen
   Advisor Dr. Weihua Geng
   2/5/2015 */

/* Inclusions */
#include <stdlib.h> /* malloc(), free() */
#include <stdio.h>  /* printf() */
#include <time.h>   /* time */
#include <math.h>   /* pow() */
#include <stdint.h> /* INT#@_MAX */

#include "mpi.h"

/* Writen Inclusions */
#include "treecode.h" /* writen in molecule */


int main(int argc, char *argv[]) {

  /* local variables */
  int i,j,k,err,idx[3],ileverl,istep;
  int ii,jj,kk;
  double t1,abserr,relerr,absinf_err,relinf_err;
  double ***f_inferr, ***f_relinferr,t;
  double tnorm;
  int maxint;
  double sttime,ettime,sdtime,edtime;
  double tttime,tdtime;

  char fname[16],density[16];

  /* functions */
  int readin(char fname[16],char density[16]);
  int compute_direct(MPI_Comm comm);
  int treecode3d_yukawa(MPI_Comm comm);

  /* MPI variables */
  int numprocs, myid, is, ie;
  int *numparsend, *numparrecv;

  /* initialize MPI */
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    fprintf(stderr,"Error in MPI_Init\n");
    return 1;
  }
  if (MPI_Comm_size(MPI_COMM_WORLD, &numprocs) != MPI_SUCCESS) {
    fprintf(stderr,"Error in MPI_Comm_size\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if (MPI_Comm_rank(MPI_COMM_WORLD, &myid) != MPI_SUCCESS) {
    fprintf(stderr,"Error in MPIC_Comm_rank\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  /* set constant */
  kappa=0.0;             // screening coefficient
  maxparnode=500;
  order = 3;
  theta = 0.5;

  /* by readin to get the surface of molecule */
  sprintf(fname,"1a63");
  sprintf(density,"3");

  if (myid == 0){
    err = readin(fname,density);
    if (err != 0){
      fprintf(stderr,"Error in readin.c\n");
      return 1;
    }
  }

  MPI_Bcast(&nface,1,MPI_INT,0,MPI_COMM_WORLD);
  numpars=nface;         // number of partucakes = number of faces

  x = (double*)calloc(numpars, sizeof(double));
  y = (double*)calloc(numpars, sizeof(double));
  z = (double*)calloc(numpars, sizeof(double));
  q = (double*)calloc(numpars, sizeof(double));
  if (x==NULL){
    fprintf(stderr, "setup error in main.c: x empty data array");
    return 1;
  }
  if (y==NULL){
    fprintf(stderr, "setup error in main.c: y empty data array");
    return 1;
  }
  if (z==NULL){
    fprintf(stderr, "setup error in main.c: z empty data array");
    return 1;
  }
  if (q==NULL){
    fprintf(stderr, "setup error in main.c: q empty data array");
    return 1;
  }

  if (myid == 0){

    printf("finished readin\n");

    for (i=0;i<numpars;i++){
      for (j=0;j<3;j++){
        idx[j] = nvert[j][i];
        r0[j] = 0.0;
        v0[j] = 0.0;
      }
      for (j=0;j<3;j++){
        for (k=0;k<3;k++){
          r0[k] += 1.0/3.0*sptpos[k][idx[j]-1];
          v0[k] += 1.0/3.0*sptnrm[k][idx[j]-1];
        }
      }
  
      /* normlaize */
      double dot_product = 0.0;
      for (j=0;j<3;j++)
        dot_product += v0[j]*v0[j];
      dot_product = sqrt(dot_product);
      for (j=0;j<3;j++)
        v0[j] = v0[j]/dot_product;
  
      x[i] = r0[0];
      y[i] = r0[1];
      z[i] = r0[2];
    }

  }

  MPI_Bcast(x,numpars,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(y,numpars,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(z,numpars,MPI_DOUBLE,0,MPI_COMM_WORLD);

  /* two tepy charge */
  for (i=0;i<numpars;i++)
    q[i]=1.0;
//    q[i] = random()/(pow(2.0,31.0)-1.0)-0.5;
//  printf("%f\n",q[i]);

  tpoten = (double*)calloc(numpars, sizeof(double));
  dpoten = (double*)calloc(numpars, sizeof(double));
  if (tpoten==NULL){
    fprintf(stderr, "setup error in main.c: tpoten empty data array");
    return 1;
  }
  if (dpoten==NULL){
    fprintf(stderr, "setup error in main.c: dpoten empty data array");
    return 1;
  }

  if (myid == 0){
    printf("  \n");
    printf("Run parameteres:  ");
    printf("                   numpars    = %d\n ",numpars);
    printf("                   kappa      = %f\n ",kappa);
    printf("                   theta      = %f\n ",theta);
    printf("                   order      = %d\n ",order);
    printf("                   maxparnode = %d\n",maxparnode);
  }


  /* compute potential by treecode */
  for (i=0;i<numpars;i++) tpoten[i]=0.0;

  sttime = MPI_Wtime();
  treecode3d_yukawa(MPI_COMM_WORLD);
  ettime = MPI_Wtime();
  tttime = ettime-sttime;

  if (myid == 0){
    printf("  \n");
    printf("Runtime for treecode is %f\n",tttime);

  /* compute potential directly */
    printf("  \n");
    printf("Computing potential - directly\n");
  }

  sdtime = MPI_Wtime();
  compute_direct(MPI_COMM_WORLD);
  edtime = MPI_Wtime();
  tdtime = edtime-sdtime;

  if (myid == 0){
    printf("  \n");
    printf("Runtime for treecode is %f\n",tdtime);
  
    printf("  \n");
    printf("Computing potential error\n");
    printf("  \n");
    
  
    abserr=0.0;
    relerr=0.0;
    relinf_err=0.0;
    absinf_err=0.0;
    for (i=0;i<numpars;i++){
      tnorm = fabs(tpoten[i]-dpoten[i]);
      relerr += tnorm*tnorm;
      abserr+=dpoten[i]*dpoten[i];
      if (tnorm>relinf_err){
        relinf_err = tnorm;
      }
      tnorm = fabs(dpoten[i]);
      if (tnorm>absinf_err){
        absinf_err = tnorm;
        maxint = i;
      }
    }
  
    relerr = sqrt(relerr/abserr);
    relinf_err = relinf_err/absinf_err;
    printf("get max %f @ %d\n",absinf_err,maxint);
    printf("Relative L2 and Inf error: %e,%e\n",relerr,relinf_err);
    printf("  \n");
  }

  if (myid == 0){
    for (i=0;i<3;i++){
      free(extr_v[i]);
      free(sptpos[i]);
      free(sptnrm[i]);
      free(atmpos[i]);
      free(nvert[i]);
    }
    for (i=0;i<2;i++){
      free(extr_f[i]);
    }
    free(extr_v);
    free(sptpos);
    free(sptnrm);
    free(extr_f);
    free(atmpos);
    free(atmrad);
    free(nvert);
  }

  free(x);

  free(y);

  free(z);

  free(q);

  free(tpoten);

  free(dpoten);

  MPI_Finalize();

  return 0;
} // end main

int compute_direct(MPI_Comm comm){
  int i,j,k;
  double tx,ty,tz,fx,fy,fz,teng,dist,t1;
  double dpeng,temp,peng,pi,tempx,tempq;

  /* MPI variables */
  int myid, numprocs,is,ie;
  int *recvcounts,*recv_disp;
  double *dpoten_sub;

  /* get MPI id */
  if (MPI_Comm_rank(comm, &myid) != MPI_SUCCESS) {
    fprintf(stderr,"MPI_Comm_rank failure\n");
    return 1;
  }

  /* get the total number of processes */
  if (MPI_Comm_size(comm, &numprocs) != MPI_SUCCESS) {
    fprintf(stderr,"MPI_Comm_size failure\n");
    return 1;
  }

  /* decompose the iteration space */
  is = myid*numpars/numprocs;
  ie = (myid+1)*numpars/numprocs;

  /* setting the position of send & recieve buff */
  dpoten_sub = (double*)calloc(ie-is,sizeof(double));
  recvcounts = (int*)calloc(numprocs,sizeof(int));
  recv_disp = (int*)calloc(numprocs,sizeof(int));
  if (dpoten_sub==NULL){
    fprintf(stderr, "setup error in direct: dpoten_sub empty data array");
    return 1;
  }
  if (recvcounts==NULL){
    fprintf(stderr, "setup error in direct: recvcounts empty data array");
    return 1;
  }
  if (recv_disp==NULL){
    fprintf(stderr, "setup error in direct: recv_disp empty data array");
    return 1;
  }
  for (i=0;i<numprocs;i++){
    recvcounts[i] = (i+1)*numpars/numprocs-i*numpars/numprocs;
    recv_disp[i] = i*numpars/numprocs;
  }

  pi = 3.141592653589793238462643;/* 24 digits of point */

  for (i=0;i<ie-is;i++){
    peng = 0.0;
    tempx=x[is+i];
    tempq=q[is+i];
    q[is+i]=0.0;
    x[is+i]=100.0;

    for (j=0;j<numpars;j++){
      tx = x[j]-tempx;
      ty = y[j]-y[is+i];
      tz = z[j]-z[is+i];
      dist = sqrt(tx*tx+ty*ty+tz*tz);
      temp = exp(-kappa*dist)/dist/4/pi;
      peng += q[j]*temp;
    }
    dpoten_sub[i]=tempq*peng;
    q[is+i]=tempq;
    x[is+i]=tempx;
  }

  MPI_Allgatherv(dpoten_sub,ie-is,MPI_DOUBLE,dpoten,recvcounts,
                 recv_disp,MPI_DOUBLE,comm);

  free(dpoten_sub);
  free(recvcounts);
  free(recv_disp);

  return 0;
}


