#===============================================================================
#
# Copyright (C) 2018 Joshua Vega (@jsvcycling)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#
#-------------------------------------------------------------------------------
#
# An implementation of the Hodkin-Huxley model in Python 3 using numpy.
#
#===============================================================================

import numpy as np
import matplotlib.pyplot as plt


# alpha_x functions
alpha_n = lambda v: 0.01*(10 - v)/(np.exp((10-v)/10) - 1)
alpha_m = lambda v: 0.1*(25 - v)/(np.exp((25-v)/10) - 1)
alpha_h = lambda v: 0.07*np.exp(-v/20)

# beta_x functions
beta_n = lambda v: 0.125*np.exp(-v/80)
beta_m = lambda v: 4*np.exp(-v/18)
beta_h = lambda v: 1/(np.exp((30 - v)/10) + 1)

# heaviside function
heaviside = lambda x: 1 if x >= 1 else 0

C_m = 1

g_K = 36
g_Na = 120
g_L = 0.3

V_K = -12
V_Na = 115
V_L = 10.6

t_max = 10000
dt = 0.01
t = np.arange(0, t_max, dt)

I_start_time = 1000
I_end_time = 5000
I = 12

V = np.zeros(t.shape)
V[0] = V_L

N = np.zeros(t.shape[0])
N[0] = alpha_n(V[0])

M = np.zeros(t.shape[0])
M[0] = alpha_m(V[0])

H = np.zeros(t.shape[0])
H[0] = alpha_h(V[0])

for i in range(0, t.shape[0] - 1):
  Iapp = I*heaviside(t[i] - I_start_time)*heaviside(I_end_time - t[i])

  I_K1 = g_K*np.power(N[i], 4)*(V[i] - V_K)
  I_Na1 = g_Na*np.power(M[i], 3)*H[i]*(V[i] - V_Na)
  I_L1 = g_L*(V[i] - V_L)

  V_1 = (Iapp - I_K1 - I_Na1 - I_L1)/C_m
  N_1 = alpha_n(V[i])*(1 - N[i]) - beta_n(V[i])*N[i]
  M_1 = alpha_m(V[i])*(1 - M[i]) - beta_m(V[i])*M[i]
  H_1 = alpha_h(V[i])*(1 - H[i]) - beta_h(V[i])*H[i]

  aV = V[i] + V_1*dt
  aN = N[i] + N_1*dt
  aM = M[i] + M_1*dt
  aH = H[i] + H_1*dt

  I_K2 = g_K*np.power(aN, 4)*(aV - V_K)
  I_Na2 = g_Na*np.power(aM, 3)*aH*(aV - V_Na)
  I_L2 = g_L*(aV - V_L)

  V_2 = (Iapp - I_K2 - I_Na2 - I_L2)/C_m
  N_2 = alpha_n(aV)*(1 - aN) - beta_n(aV)*aN
  M_2 = alpha_m(aV)*(1 - aM) - beta_m(aV)*aM
  H_2 = alpha_h(aV)*(1 - aH) - beta_h(aV)*aH

  V[i+1] = V[i] + (V_1 + V_2)*dt/2
  N[i+1] = N[i] + (N_1 + N_2)*dt/2
  M[i+1] = M[i] + (M_1 + M_2)*dt/2
  H[i+1] = H[i] + (H_1 + H_2)*dt/2

print(V[-1])

# Optional plotting functionality. Uncomment to enable (requires matplotlib).
# plt.plot(t, V)
# plt.show()