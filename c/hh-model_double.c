 /*
 *==============================================================================
 *
 * Copyright (C) 2018 Joshua Vega (@jsvcycling)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 *
 *------------------------------------------------------------------------------
 *
 * An implementation of the Hodkin-Huxley model in C using doubles.
 *
 *==============================================================================
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/* Using preprocessor statements to reduce memory calls. */
#define C_m 1.f

#define g_K   36.0
#define g_Na  120.0
#define g_L   0.3

#define V_K   -12.0
#define V_Na  115.0
#define V_L   10.6

#define t_max   10000.0
#define dt      0.01

#define I_start_time  1000.0
#define I_end_time    5000.0
#define I             12.0

double *t;
double *V;
double *N;
double *M;
double *H;

double alpha_n(double v) {
  return 0.01f*(10.0 - v) / (exp((10.0 - v)/10.0) - 1.0);
}

double alpha_m(double v) {
  return 0.1f*(25.0 - v) / (exp((25.0 - v)/10.0) - 1.0);
}

double alpha_h(double v) {
  return 0.07 * expf(-v / 20.);
}

double beta_n(double v) {
  return 0.125 * exp(-v / 80.0);
}

double beta_m(double v) {
  return 4.0 * exp(-v / 18.0);
}

double beta_h(double v) {
  return 1.0 / (exp((30.0 - v)/10.0) + 1.0);
}

double heaviside(double x) {
  if (x >= 1.0) return 1.0;
  else return 0.0;
}

int main(int argc, char **argv) {
  int num_ts = (int)ceil(t_max / dt);

  t = (double *)malloc(sizeof(double) * num_ts);
  V = (double *)malloc(sizeof(double) * num_ts);
  N = (double *)malloc(sizeof(double) * num_ts);
  M = (double *)malloc(sizeof(double) * num_ts);
  H = (double *)malloc(sizeof(double) * num_ts);

  t[0] = 0.0;
  
  for (int i = 1; i < num_ts; i++) {
    t[i] = t[i-1] + dt;
  }

  V[0] = V_L;
  N[0] = alpha_n(V[0]);
  M[0] = alpha_m(V[0]);
  H[0] = alpha_h(V[0]);

  for (int i = 0; i < num_ts - 1; i++) {
    double Iapp = I*heaviside(t[i] - I_start_time)*heaviside(I_end_time - t[i]);

    double I_K1 = g_K*pow(N[i], 4)*(V[i] - V_K);
    double I_Na1 = g_Na*pow(M[i], 3)*H[i]*(V[i] - V_Na);
    double I_L1 = g_L*(V[i] - V_L);

    double V_1 = (Iapp - I_K1 - I_Na1 - I_L1)/C_m;
    double N_1 = alpha_n(V[i])*(1.0 - N[i]) - beta_n(V[i])*N[i];
    double M_1 = alpha_m(V[i])*(1.0 - M[i]) - beta_m(V[i])*M[i];
    double H_1 = alpha_h(V[i])*(1.0 - H[i]) - beta_h(V[i])*H[i];

    double aV = V[i] + V_1*dt;
    double aN = N[i] + N_1*dt;
    double aM = M[i] + M_1*dt;
    double aH = H[i] + H_1*dt;

    double I_K2 = g_K*pow(aN, 4)*(aV - V_K);
    double I_Na2 = g_Na*pow(aM, 3)*aH*(aV - V_Na);
    double I_L2 = g_L*(aV - V_L);

    double V_2 = (Iapp - I_K2 - I_Na2 - I_L2)/C_m;
    double N_2 = alpha_n(aV)*(1.0 - aN) - beta_n(aV)*aN;
    double M_2 = alpha_m(aV)*(1.0 - aM) - beta_m(aV)*aM;
    double H_2 = alpha_h(aV)*(1.0 - aH) - beta_h(aV)*aH;

    V[i+1] = V[i] + (V_1 + V_2) * dt / 2.0;
    N[i+1] = N[i] + (N_1 + N_2) * dt / 2.0;
    M[i+1] = M[i] + (M_1 + M_2) * dt / 2.0;
    H[i+1] = H[i] + (H_1 + H_2) * dt / 2.0;
  }

  printf("%f\n", V[num_ts - 1]);

  /*
   * Optional plotting functionality. Uncomment to enable (requires gnuplot).
   */
  /*
  FILE *gnuplot = popen("gnuplot -persistent", "w");
  fprintf(gnuplot, "set term png large size 1280,720\n");
  fprintf(gnuplot, "set output 'output_double.png'\n");
  fprintf(gnuplot, "set xrange [%f:%f]\n", 1000.0, 1050.0);
  fprintf(gnuplot, "set yrange [%f:%f]\n", -20.0, 120.0);
  fprintf(gnuplot, "set datafile separator ','\n");
  fprintf(gnuplot, "plot '-' with lines\n");

  for (int i = 0; i < num_ts; i++) {
    fprintf(gnuplot, "%f,%f\n", t[i], V[i]);
    fflush(gnuplot);
  }

  fprintf(gnuplot, "e\n");
  fflush(gnuplot);
  
  pclose(gnuplot);
  */

  free(t);
  free(V);
  free(N);
  free(M);
  free(H);

  exit(EXIT_SUCCESS);
}
