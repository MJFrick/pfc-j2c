// Define dependencies
#include <stdlib.h>
#include <math.h>
#include <fftw3-mpi.h>
#include <float.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include "pfc_header.h"
#define pi 2.0*acos(0.0)
#define third 1.0 / 3.0
 void Calc_halfReals(state* s)
{
double Lx;
double Ly;
double kx2;
double ky2;
double k2;
double M;
double dt;
double Lambda;
int index;
    Lx = s->Nx * s->dx;
    Ly = s->Ny * s->dy;
for(int i = 0; i < s->local_n1; i++)
   {
        kx2 = i < ((s->Ny>>1) + 1) ? pow((2*pi * i / Lx),2) : pow((2*pi * (i-s->Nx) / Lx),2);
for(int j = 0; j < s->Nx; j++)
       {
index = i * s->Nx + j;
            ky2 = j < ((s->Ny>>1) + 1) ? pow((2*pi * j / Ly),2) : pow((2*pi * (j-s->Ny) / Ly),2);
            s->k2[index] = (kx2 + ky2);
           }
       }
for(int i = 0; i < s->local_n1; i++)
   {
for(int j = 0; j < s->Nx; j++)
       {
index = i * s->Nx + j;
    s->C2 [index] = s->Bx*(2.0*s->k2[index] - s->k2[index]*s->k2[index]);
       }
   }
for(int i = 0; i < s->local_n1; i++)
   {
for(int j = 0; j < s->Nx; j++)
       {
index = i * s->Nx + j;
            k2 = s->k2[index];
            M = s->M;
            dt = s->dt;
            Lambda = 1.0 - s->C2[index];
            s->P[index] = -M * dt * k2 * Lambda / (1.0 + M * dt* k2 * Lambda);
            s->Q[index] = -M * dt * k2 / (1.0 + M * dt * k2 * Lambda);
           }
       }
    return ;
}
 void Calc_NL(state* s)
{
double norm;
int index;
    norm = 1.0 / (s->Nx * s->Ny);
for(int i = 0; i < s->local_n0; i++)
   {
for(int j = 0; j < s->Ny; j++)
       {
index = 2 * i * ((s->Ny>>1)+1) + j;
    s->NL [index] = -0.5 * s->t * s->n[index] * s->n[index] + third * s->v * s->n[index] * s->n[index] * s->n[index];
       }
   }
    fftw_mpi_execute_dft_r2c(s->forward, s->NL, s->NL_k );
for(int i = 0; i < s->local_n1; i++)
   {
for(int j = 0; j < s->Nx; j++)
       {
index = i * s->Nx + j;
    s->NL_k [index][0] *=       norm;
    s->NL_k [index][1] *=       norm;
       }
   }
    return ;
}
 void step(state* s)
{
double norm;
int index;
    norm = 1.0 / (s->Nx * s->Ny);
    fftw_mpi_execute_dft_r2c(s->forward, s->n, s->nk );
for(int i = 0; i < s->local_n1; i++)
   {
for(int j = 0; j < s->Nx; j++)
       {
index = i * s->Nx + j;
    s->nk [index][0] *=       norm;
    s->nk [index][1] *=       norm;
       }
   }
    Calc_NL(s);
for(int i = 0; i < s->local_n1; i++)
   {
for(int j = 0; j < s->Nx; j++)
       {
index = i * s->Nx + j;
    s->nk [index][0] += s->P[index] * s->nk[index][0] + s->Q[index] * s->NL_k[index][0];
    s->nk [index][1] += s->P[index] * s->nk[index][1] + s->Q[index] * s->NL_k[index][1];
       }
   }
    fftw_mpi_execute_dft_c2r(s->back, s->nk, s->n );
    s->step += 1;
    s->tm += s->dt;
    return ;
}
 void seed(state* s, double q, double A)
{
double x;
double y;
int index;
for(int i = 0; i < s->local_n0; i++)
   {
for(int j = 0; j < s->Ny; j++)
       {
index = 2 * i * ((s->Ny>>1)+1) + j;
            x = i * s->dx;
            y = j * s->dy;
            s->n[index] += exp(-(pow((i - (s->Nx >>1)),2)+pow((j - (s->Ny >> 1)),2)) / 300) * 2.0 * A * (cos(q * x) + 2.0 * cos(0.5*(q * x)) * cos(0.5 * sqrt(3) * q * y));
           }
       }
    return ;
}
