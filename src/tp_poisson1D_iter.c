/******************************************/
/* tp2_poisson1D_iter.c                   */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"

int main(int argc,char *argv[])
/* ** argc: Number of arguments */
/* ** argv: Values of arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, lab, kv;
  int *ipiv;
  int info;
  int NRHS;
  int IMPLEM = 0;
  double T0, T1;
  double *RHS, *SOL, *EX_SOL, *X;
  double *AB;

  if (argc == 2) {
    IMPLEM = atoi(argv[1]);
  } else if (argc > 2) {
    perror("Application takes at most one argument");
    exit(1);
  }

  /* Size of the problem */
  NRHS=1;
  nbpoints=12;
  la=nbpoints-2;

  /* Dirichlet Boundary conditions */
  T0=5.0;
  T1=20.0;

  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  SOL=(double *) calloc(la, sizeof(double)); 
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);

  /* Setup the Poisson 1D problem */
  /* General Band Storage */
  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  kv=0;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;
  
  AB = (double *) malloc(sizeof(double)*lab*la);
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  
  /* uncomment the following to check matrix A */
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

  /* Write solution */
  write_vec(SOL, &la, "SOL.dat");

  free(RHS);
  free(SOL);
  free(EX_SOL);
  free(X);
  free(AB);
  printf("\n\n--------- End -----------\n");
}
