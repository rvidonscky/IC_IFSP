# RADIAÇÃO DIRETA ---------------------------------------------------------------------------------------------------------------------

import math

# conversões úteis:
def deg2rad(x):
    return x * math.pi / 180.0

def rad2deg(x):
    return x * 180.0 / math.pi

def SolarDeclination(day_of_year):
    deg = 23.45 * math.sin(deg2rad(360.0 * (284 + day_of_year) / 365.0))
    return deg2rad(deg)

def hour_angle(hora_decimal):                       #A alteração da hora de pico do Sol e da alteração do angulo com a variação da hora
    return deg2rad((hora_decimal - 12.5) * 17.5)    #Serviram para ajustar o gráfico plotado

def solar_altitude(latitude_deg, decl_rad, hour_ang_rad):
    i = deg2rad(latitude_deg)
    return math.asin(math.sin(i) * math.sin(decl_rad) + math.cos(i) * math.cos(decl_rad) * math.cos(hour_ang_rad))

# radiação direta
def direct_radiation(Ic, B, b):
    return Ic * math.exp(-B / math.sin(b))

# RADIAÇÃO DIFUSA ---------------------------------------------------------------------------------------------------------------------

def FSS(tilt_deg):
    S = deg2rad(tilt_deg)
    return (1.0 + math.cos(S)) / 2.0

def diffuse_radiation(IDN, tilt_deg):
    C = 0.135
    fss = FSS(tilt_deg)
    Idq = C * IDN * fss
    return Idq

# RADIAÇÃO REFLETIDA ----------------------------------------------------------------------------------------------------------------

def Fsg(tilt_deg):
    S = deg2rad(tilt_deg)
    return (1.0 - math.cos(S)) / 2.0

def reflected_radiation(IDN, Idq, tilt_deg):
    krq = 0.2
    fsg = Fsg(tilt_deg)
    Irq = (IDN + Idq) * krq * fsg
    return Irq

# RADIAÇÃO TOTAL ----------------------------------------------------------------------------------------------------------------------

import matplotlib.pyplot as plt

# Função para calcular o cosseno do ângulo de incidência via vetores
# O vetor do Sol (s) e o vetor normal do painel (n) são definidos em coordenadas locais
def cos_incidence_vector(latitude_deg, decl_rad, tilt_deg, hr, b):

    i = deg2rad(latitude_deg)
    S = deg2rad(tilt_deg)
    
    Azimute = math.asin(math.cos(decl_rad) * math.sin(hr) / math.cos(b))
    cosq = math.cos(b) * math.cos(Azimute) * math.sin(S) + math.sin(b) * math.cos(S)
    return cosq
    

# Radiação total (Eq. 8)
def total_radiation(IDN, Idq, Irq, cosq):
    return IDN * cosq + Idq + Irq

# ----------------------------------------------------------------------------------------------
# PARÂMETROS FIXOS DO ARTIGO (caso de verao, agosto)
Ic = 1357          
latitude_deg = 32  #deg
tilt_deg = 37      #deg
B = 0.21           
day_of_year = 200  
# ----------------------------------------------------------------------------------------------

# CÁLCULOS HORA A HORA (com passo menor pra suavizar)
horas = [h for h in [x * 0.25 for x in range(9 * 4, 15 * 4 + 1)]]  # 9h até 15h, a cada 15 min
It_values = []
horas_solares = []  # só para as horas em que há sol (sinb > 0)

decl_rad = SolarDeclination(day_of_year)

for hora in horas:

    Hr = hour_angle(hora)
    b = solar_altitude(latitude_deg, decl_rad, Hr)
    
    IDN = direct_radiation(Ic, B, b)          
    Idq = diffuse_radiation(IDN, tilt_deg)    
    Irq = reflected_radiation(IDN, Idq, tilt_deg) 
    cosq = cos_incidence_vector(latitude_deg, decl_rad, tilt_deg, Hr, b)
    It = total_radiation(IDN, Idq, Irq, cosq)
    
    It_values.append(It)
    horas_solares.append(hora)

# ----------------------------------------------------------------------------------------------
# PLOT DO GRÁFICO FINAL (como o do artigo)
plt.figure(figsize=(8,5))
plt.plot(horas_solares, It_values, color='orange', linewidth=2)
plt.grid(True)
plt.tight_layout()
plt.show()
