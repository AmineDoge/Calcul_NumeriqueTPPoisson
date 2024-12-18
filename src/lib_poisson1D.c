/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  int i, j;

  for (j = 0; j < *lab; j ++) 
  {
    for (i = 0; i < *la; i ++) 
    {
      AB[i * (*lab) + j] = 0.0;
    }
  }
  for (i = 0; i < *la; i ++) 
  {
    AB[i * (*lab) + (*kv) + 1] = 2.0; 
  }
  for (i = 0; i < *la - 1; i ++) 
  {
    AB[(i + 1) * (*lab) + (*kv)] = -1.0; 
    AB[i * (*lab) + (*kv) + 2] = -1.0;   
  }
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  int i, j;
  for (j = 0; j < *la; j++) {
      for (i = 0; i < *lab; i++) {
            AB[i + j * (*lab)] = 0.0;
        }
  }

  for (j = 0; j < *la; j++) {
      AB[1 + j * (*lab)] = 1.0;
  }
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  RHS[0] = *BC0;
  for (int i = 1; i < *la-1; i ++) 
  {
    RHS[i] = 0.0;
  }
  RHS[*la - 1] = *BC1;
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  double T0 = *BC0;
  double T1 = *BC1;
  for (int i = 0; i < *la; i++) 
  {
    EX_SOL[i] = X[i] * (T1 - T0) + T0;
  }
}  

void set_grid_points_1D(double* x, int* la){
  for (int i = 0; i < *la; i++) 
  {
    x[i] = (i + 1) * 1.0 / (*la + 1);
  }
}

double relative_forward_error(double* x, double* y, int* la){
  double num = 0.0;
  double denum = 0.0;

  for (int i = 0; i < *la; i++) {
    num += (x[i] - y[i]) * (x[i] - y[i]);   
    denum += y[i] * y[i];
  }

  return sqrt(num / denum);
}

int indexABCol(int i, int j, int *lab){
  return (1 + i - j) + j * (*lab);
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
    int i;

    *info = 0;

    for (i = 0; i < (*n - 1); i++) {
      if (AB[1 + i * (*lab)] == 0.0) {
          *info = i + 1;
          return *info;
      }

        
      AB[0 + (i + 1) * (*lab)] /= AB[1 + i * (*lab)];

        
      AB[1 + (i + 1) * (*lab)] -= AB[0 + (i + 1) * (*lab)] * AB[2 + i * (*lab)];

        
      ipiv[i] = i + 1;
    }
    ipiv[*n - 1] = *n; 
    return *info;
}

