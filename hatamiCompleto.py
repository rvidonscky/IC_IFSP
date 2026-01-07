import CoolProp.CoolProp as CP
import numpy as np
from devito import *

import matplotlib.pyplot       as plt
import matplotlib.ticker       as mticker
from   mpl_toolkits.axes_grid1 import make_axes_locatable
from   matplotlib              import cm
from   matplotlib              import ticker

# Coeficientes
i = 32 * np.pi/180  #Latitude
S = 37 * np.pi/180  #Angulo do coletor
n = 200   #Dia do ano
k_g = 5.9 #Condutividade Térmica do Vidro
C_g = 753 #Calor específico do vidro
ro_g = 2466  #Densidade do Vidro
delta_g = 4e-3  #Espessura do vidro
l_p = 2 #Largura do Vidro
ro_a = 1.09 #Densidade do ar
k_a = 0.025 #Condutividade do ar
C_a = 1005  #Calor específico do ar
alpha_g = 0.1 #Coeficiente de absorvição do ar
alpha_p = 0.96  #Coeficiente de absorvição da placa
tau_g = 0.8 #Coeficiente transiente do vidro
C_p = 900 #Calor específico da placa
ro_w = 1000 #Densidade da água
L = 2 #Largura do coletor
b = 1 #Comprimento do coletor
ro_d = 960 #Densidade do material seco
a = 0.1 #Altura do duto coletor
k_p = 250 #Condutividade térmica da placa
M_0 = 3.5 #Conteúdo úmido inicial
sigma = 5.57e-8
A_C = b * L #Área do coletor
epslon_p = 1 #Emitancia, como é um corpo negro, é 1
ro_p = 2700 #Densidade da placa (alumínio)
P_atm = 101325 #Pressão atmosférica (Pa)

VelocidadeAr = 1 #Pode ser alterada

#import ay_eq_radiacao
# Radiação solar total fixa pra testar inicialmente:
It = 800.0

# time grid
from examples.seismic.source import Receiver, TimeAxis
t0 = 0
tn = 86400
dt = 0.01                       #Devido ao número de Courant, para o sistema entregar valores estável
nt = int((tn-t0)/dt)            # o precisa ser dt < 0.1. Menor, melhor

time_range = TimeAxis(start=t0,stop=tn,num=nt+1)

# Grid
grid = Grid(shape=(10), extent=(1.))

# Temperaturas
TemperaturaVidro = TimeFunction(name='T_g', grid=grid, space_order=2)
TemperaturaPlaca = TimeFunction(name='T_p', grid=grid, space_order=2)
TemperaturaAr = TimeFunction(name='T_a', grid=grid, space_order=2)
TemperaturaAmbiente = Function(name='T_am', grid=grid, space_order=2)

#Fatores Transferencias de calor dependentes das temperaturas
h_rpg = TimeFunction(name='h_rpg', grid=grid, space_order=1)
h_rpsky = TimeFunction(name='h_rpsky', grid=grid)

# Termos para trabalhar nas equações de contorno
x = grid.dimensions[0]
t = grid.stepping_dim
dz = grid.spacing[0]

#Condição inicial (eq. 16) considerando T0 = 35°C
TemperaturaVidro.data[0,:] = 308.15
TemperaturaPlaca.data[0,:] = 308.15
TemperaturaAmbiente.data[:] = 308.15          #Por hora considerando que a temperatura ambiente não varia com o tempo

# Termos dependentes das temperaturas
TemperaturaCeu = 0.0552*(TemperaturaAmbiente**1.5)
T_ref = 300.0
viscosidadeDinamica = CP.PropsSI("V", "T", T_ref, "P", P_atm, "Air")
calorEspecifico = CP.PropsSI("CPMASS", "T", T_ref, "P", P_atm, "Air")
condutividadeTermica = CP.PropsSI("CONDUCTIVITY", "T", T_ref, "P", P_atm, "Air")

Re = ro_a * VelocidadeAr * L / viscosidadeDinamica
Pr = calorEspecifico * viscosidadeDinamica / condutividadeTermica

if Re * Pr * a/L > 70:
    Nu = 7.6
else:
    Nu = 1.85*pow(Re * Pr * a / L, 1/3)

h_ga = Nu * k_a / L
h_pa = h_ga

Nu_gam = 0.86 * pow(Re, 0.5) * pow(Pr, 1/3)
h_gam = Nu_gam * k_a / L

# Equações Diferenciais
eq11 = Eq(TemperaturaVidro.dt,
          (k_g/(C_g*ro_g))*TemperaturaVidro.dx2 + It*alpha_g/(ro_g*C_g*delta_g))

eq12 = Eq(TemperaturaPlaca.dt,
          It*tau_g*alpha_p/(ro_p*C_p*l_p)
          + (h_pa/(ro_p*C_p*l_p))*(TemperaturaAr - TemperaturaPlaca)
          - (h_rpsky/(ro_p*C_p*l_p)) * (TemperaturaPlaca - TemperaturaCeu)
          - (h_rpg/(ro_p*C_p*l_p)) * (TemperaturaPlaca - TemperaturaVidro))

eq13 = Eq(TemperaturaAr.dt,
          (-1) * VelocidadeAr * TemperaturaAr.dxl
          + (h_ga*b/(ro_a*C_a*A_C))*(TemperaturaVidro - TemperaturaAr)
          + (h_pa*b/(ro_p*C_p*A_C))*(TemperaturaPlaca - TemperaturaAr))

eq_1 = Eq(TemperaturaVidro.forward, solve(eq11, TemperaturaVidro.forward))
eq_2 = Eq(TemperaturaPlaca.forward, solve(eq12, TemperaturaPlaca.forward))
eq_3 = Eq(TemperaturaAr.forward, solve(eq13, TemperaturaAr.forward))

#Equações dos coeficientes de transf. de calor
eq_rpsky = Eq(h_rpsky.forward,
              sigma * (TemperaturaPlaca**2 + TemperaturaCeu**2)
              * (TemperaturaPlaca + TemperaturaCeu)
              / (1/epslon_p + 1/tau_g - 1))

eq_rpg = Eq(h_rpg.forward,
            sigma * (TemperaturaPlaca**2 + TemperaturaVidro**2)
            * (TemperaturaPlaca + TemperaturaVidro)
            / (1/epslon_p + 1/tau_g - 1))

#Condições de contorno
eq_entradaAr = Eq(TemperaturaAr.forward.subs(x, 0),                             #A função subs diz onde eu quero
                TemperaturaAmbiente.subs(x, 0))                                 #que a igualdade seja válida

eq_17 = Eq(TemperaturaVidro.forward.subs(x, 0),                                 #Isolando TempVidro(z=0)
           (k_g/dz * TemperaturaVidro.forward.subs(x, x+dz) + h_ga * TemperaturaAr.forward.subs(x, 0))
           / (h_ga + k_g/dz))

eq_18 = Eq(TemperaturaVidro.forward.subs(x, x.symbolic_max),                    #Isolando TempVidro(z=z_max)
           (k_g/dz * TemperaturaVidro.forward.subs(x, x.symbolic_max - dz) + h_ga * TemperaturaAr.forward.subs(x, x.symbolic_max))
           / (h_ga + k_g/dz))

#Recebedor
nrec = 1
rec_T = Receiver(name="rec_T", grid=grid, npoint=nrec, time_range=time_range)
rec_T.coordinates.data[:, 0] = 0.5
rec_term = rec_T.interpolate(expr=TemperaturaVidro.forward)

#Operador
op = Operator([eq_1, eq_2, eq_3,
                eq_entradaAr, eq_17, eq_18, 
                eq_rpsky, eq_rpg,] + rec_term)

op(dt=dt, time=nt)

# Plot Configuration
#==============================================================================
plt.rc('text' , usetex=False)
plt.rc('font' , family='serif')
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)

plt.plot(time_range.time_values, rec_T.data, label=f"Receiver")
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.title("Receiver Data")
plt.legend()
plt.grid()
plt.show()
