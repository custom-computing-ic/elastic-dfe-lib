#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "Maxfiles.h"
#include "MaxSLiCInterface.h"
#include "OptionPricing.hpp"

#define real double
#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

double optionPricing(double strike, double sigma, double timestep, int numMaturity,
                     int paraNode, int numPathGroup, double T,
                     double **out_rand1, double **out_rand2
                     ) {

  //need to change the value in the java file
  int numPE = 4;

  real *maturity = malloc(sizeof(real) * numPE * numMaturity);
  real *maturity_diff = malloc(sizeof(real) * numPE * numMaturity);
  real *f = malloc(sizeof(real) * numPE * numMaturity);
  real *fout = malloc(sizeof(real) * 2048);
  unsigned int *seed1 = malloc(sizeof(unsigned int) * 64*100);
  unsigned int *seed2 = malloc(sizeof(unsigned int) * 64*100);
  int i, j, k;
  for(i=0; i<64; i++){
    seed1[i] = i;
    seed2[i] = 7-i;
    /* printf("%d %d, ", i, 7-i); */
  }

  for(i=0;i<16;i++){
    for(j=0;j<numPE;j++){
      seed1[i*4*numPE+j*4] = i*4;
      seed1[i*4*numPE+j*4+1] = i*4+1;
      seed1[i*4*numPE+j*4+2] = i*4+2;
      seed1[i*4*numPE+j*4+3] = i*4+3;
      /* printf("%d %d %d %d\n", seed1[i*4*numPE+j*4], seed1[i*4*numPE+j*4+1], seed1[i*4*numPE+j*4+2], seed1[i*4*numPE+j*4+3]); */
      seed2[i*4*numPE+j*4]   = 7- i*4;
      seed2[i*4*numPE+j*4+1] = 7-(i*4+1);
      seed2[i*4*numPE+j*4+2] = 7-(i*4+2);
      seed2[i*4*numPE+j*4+3] = 7-(i*4+3);
      /* printf("%d %d %d %d\n", seed2[i*4*numPE+j*4], seed2[i*4*numPE+j*4+1], seed2[i*4*numPE+j*4+2], seed2[i*4*numPE+j*4+3]); */
    }
  }

  /* for(i=0;i<16*4*4; i++){ */
  /*   printf("%d ", seed1[i]); */
  /* } */


  for(j=0;j<numPE;j++){
    maturity[0*numPE+j] = 0;
    maturity[1*numPE+j] = 1;
    maturity[2*numPE+j] = 3;
    maturity[3*numPE+j] = 5;
    maturity[4*numPE+j] = 10;
    maturity[5*numPE+j] = 30;

    f[0*numPE+j] = 0.0695;
    f[1*numPE+j] = 0.075;
    f[2*numPE+j] = 0.078;
    f[3*numPE+j] = 0.0797;
    f[4*numPE+j] = 0.082;
    f[5*numPE+j] = 0.086;
    f[6*numPE+j] = 0.086;
    f[7*numPE+j] = 0.086;
  }

  for(j = 0; j<numMaturity-1; j++){
    for(k=0;k<numPE;k++){
      maturity_diff[j*numPE+k] = 0;
      if(j>0)
        maturity_diff[j*numPE+k] = maturity[j*numPE+k]- maturity[(j-1)*numPE+k];
    }
  }

  int outputRand = 1;
  int initMax = 8;
  real* rand1 = (real*)malloc(numPE*(numPathGroup)*paraNode*sizeof(real));
  real* rand2 = (real*)malloc(numPE*(numPathGroup)*paraNode*sizeof(real));
  printf("strike=%f, run for %d times\n", strike, numPathGroup*paraNode);

  if (out_rand1 != NULL) {
    *out_rand1 = rand1;
  }

  if (out_rand2 != NULL) {
    *out_rand2 = rand2;
  }

  OptionPricing(
                numPE*initMax,//initsize
                numPE*paraNode,//nodesize
                numPE*(numPathGroup)*paraNode,//pathsize
                numPE*64,//seedsize
                initMax + (numMaturity-1)*(numPathGroup)*paraNode,//ticks
                T,//T
                exp(-f[0]*T),//discount
                numMaturity, // numMaturity
                numPathGroup,//numPath
                outputRand,
                sigma,
                sqrt(T),
                strike,
                f,//fin
                maturity,
                maturity_diff,
                seed1,
                seed2,
                rand1,
                rand2,
                fout//result
                );

  printf("done!\n");
  //  for(i = 0; i< numMaturity; i++){
  //    printf("fin[%d] = %lf\n", i, f[i*numPE]);
  //  }

  //add them together
  real sum = 0;
  for(i = 0; i< paraNode; i++){
    //for(j=0;j<numPE;j++){
    //  if(j==0){
    //          printf("fout[%d] = %lf\n", i, fout[i*numPE+j]);
    //          sum += fout[i*numPE+j];
    //  }
    //}
    sum += fout[i*numPE];
  }
  /*for(i = 0; i< (numPathGroup)*paraNode; i++){
    printf("rand1[%d] = %lf\n", i, rand1[i]);
    }*/

  double dfeResult = sum/(numPathGroup)/paraNode;
  printf("result = %lf\n", dfeResult);


  return dfeResult;
}


double cpuOptionPricing(
                        double strike, double sigma, double timestep, int numMaturity,
                        int paraNode, int numPathGroup, double T,
                        double *rand1, double *rand2)
{
  int numPE = 4;
  real *maturity = malloc(sizeof(real) * numPE * numMaturity);
  real *maturity_diff = malloc(sizeof(real) * numPE * numMaturity);
  real *f = malloc(sizeof(real) * numPE * numMaturity);

  double sum = 0.0;
  int N = (numPathGroup)*paraNode;
  maturity[0] = 0;
  maturity[1] = 1;
  maturity[2] = 3;
  maturity[3] = 5;
  maturity[4] = 10;
  maturity[5] = 30;

  f[0] = 0.0695;
  f[1] = 0.075;
  f[2] = 0.078;
  f[3] = 0.0797;
  f[4] = 0.082;
  f[5] = 0.086;
  f[6] = 0.086;
  f[7] = 0.086;

  for(int j = 0; j<numMaturity-1; j++){
    maturity_diff[j] = 0;
    if(j>0)
      maturity_diff[j] = maturity[j]- maturity[(j-1)];
  }

  for(int k = 0; k< N;k++){
    /*-------------one MC iteration----------------------*/
    real discount =1;
    real swapPrice = 0;
    //at timeline[0], f[0 to numMaturity-1] is observed
    f[0] = 0.0695;
    f[1] = 0.075;
    f[2] = 0.078;
    f[3] = 0.0797;
    f[4] = 0.082;
    f[5] = 0.086;
    //just single jump, no time steps involved
    //for(int i = 1; i< numTimestep+1; i++)
    //{
    real A1 = 0;
    real A2 = 0;
    real Bprev = 0;
    //discount based on short rate up until t-1, constant across all paths
    //discount *= exp(-f[0] * T);
    discount *= exp(-f[0] * T);
    //printf("%f\n", f[(i-1)*numMaturity]*timeline[i]);
    real Z1 = rand1[k*numPE];
    real Z2 = rand2[k*numPE];

    real zeroRate = 0;
    real zeroPrice_i_j = 0;
    real annuity=0;
    for(int j = 0; j<numMaturity-1; j++)
      {
        /*-------------------miu-----------------*/
        //time between to maturity date, i.e. viewed from timeline[i]
        //double timespan = 0;
        //if(j>0)
        //  timespan = maturity[j]- maturity[j-1];
        //double vol1 = fvol1(maturity[j]);
        //double vol2 = fvol2(maturity[j]);

        //use textbook miu for testing, old value overwritten
        //miu[j]= fmiu(maturity[j]);
        real realT = maturity[j]+T;
        real miu= 0.5*sigma*sigma*(realT*realT- maturity[j]*maturity[j]);
        real vol1 = sigma;
        real vol2 = sigma;
        /*------------------------------------*/
        //f_new = f_old+miu*dt+vol*sqrt(dt)*dW
        f[j] = f[j+1]+
          miu + sqrt(T)*(vol1*Z1 + vol2*Z2);
        //printf("f1 = %f, miu = %f, sqrt_T=%f, N1 = %f, N2 = %f, fnew = %f, ftime=%f\n", f[j+1],  miu, sqrt(T), Z1, Z2, f[j], f[j]*maturity_diff[j]);

        zeroRate += f[j]*maturity_diff[j];
        zeroPrice_i_j = exp(-zeroRate);

        //printf("j = %d, Local zero rate is %f,\n", j,zeroRate);
        annuity += zeroPrice_i_j;
        //printf("timespan=%f, miu=%f\n", timespan, miu[j]);
      }
    real swapRate= (1-zeroPrice_i_j)/annuity;
    //printf("Local annuity is %f, swap rate is %f, zeroPrice is %f\n", annuity, swapRate, zeroPrice_i_j);
    real payoff_i = (swapRate-strike)>0? (swapRate-strike)*annuity : 0;
    swapPrice = payoff_i*discount;
    //printf("Local swap rate is %f, pay off is %f, discount is %f, price is %f\n", swapRate, payoff_i, discount, swapPrice);
    //}
    //printf("*******Local swap price[%d] is %f\n", k, swapPrice);
    sum += swapPrice;
  }

  free(rand1);
  free(rand2);

  return sum/N;
}

/* int main(int argc, char *argv[]) { */

/*   // -- Parameters -- */
/*   real strike = 0.01; */
/*   real sigma = 0.02; */
/*   real timestep = 0.05; */
/*   int numTimeStep = (int)(10/0.05); */
/*   int numMaturity = 2000000; */
/*   int paraNode = 50; */
/*   int numPathGroup = 20; */
/*   int numPE = 4; */
/*   double *rand1, *rand2; */
/*   real T = 10; */

/*   double dfeResult = optionPricing(strike, sigma, timestep, numMaturity, */
/*                                    paraNode, numPathGroup, T, */
/*                                    &rand1, &rand2); */

/*   double cpuResult = cpuOptionPricing(strike, sigma, timestep, numMaturity, */
/*                                       paraNode, numPathGroup, T, */
/*                                       rand1, rand2); */

/*   if (fabs(cpuResult - dfeResult) > 1E-6) { */
/*     printf("Error! Expected: %lf Got: %lf\n", cpuResult, dfeResult); */
/*     return 1; */
/*   } */

/*   return 0; */
/* } */
