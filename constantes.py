import CoolProp.CoolProp as CP
import numpy as np
from devito import *

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
alpha_g = 0.1 #Coeficiente de absorção do ar
alpha_p = 0.96  #Coeficiente de absorção da placa
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

VelocidadeAr = 0.1468 #Pode ser alterada

# time grid
from examples.seismic.source import Receiver, TimeAxis
t0 = 0 + 60 * 60 * 9
tn = 60*60*8 + 60 * 60 * 9
dt = 0.05                     #Devido ao número de Courant, para o sistema entregar valores estável o dt deve estar em um intervalo ótimo
nt = int((tn-t0)/dt)

time_range = TimeAxis(start=t0,stop=tn,num=nt+1)

# Discretização da malha
nx = 50  # Número de pontos ao longo de x
nz = 5  # Número de pontos ao longo de z (resulta em passos de 1mm: 0, 1, 2, 3, 4)

# Termos dependentes das temperaturas
T_ref = 308.15
viscosidadeDinamica = CP.PropsSI("V", "T", T_ref, "P", P_atm, "Air")
calorEspecifico = CP.PropsSI("CPMASS", "T", T_ref, "P", P_atm, "Air")
condutividadeTermica = CP.PropsSI("CONDUCTIVITY", "T", T_ref, "P", P_atm, "Air")

#Cálculo dos termos de Tranfêrencia de Calor
Re = ro_a * VelocidadeAr * L / viscosidadeDinamica
Pr = calorEspecifico * viscosidadeDinamica / condutividadeTermica
if Re * Pr * a/L > 70:
    Nu = 7.6
else:
    Nu = 1.85 * pow(Re * Pr * a / L, 1/3)

h_ga = Nu * k_a / L
h_pa = h_ga


Nu_gam = 0.86 * pow(Re, 0.5) * pow(Pr, 1/3)
h_gam = Nu_gam * k_a / L
