"""
        Nomeclaturas
    i = Latidude
    Hr = Angulo da Hora (Definido como a diferença entre a hora observada e o meio dia, vezes 15 )
    delta = Declinação Solar
    Ic = Constante Solar (Assumida em 1357 W*m^(-2) )
    B = Coeficiente de extinção atmosférica ( O quanto de luz é absorvida e espalhada pela atmosfera)
    Beta = Angulo da altitude solar
    IdTeta = Radiação difusa
    C = the ratio of diffusion radiation to the direct ones on the horizontal surface
    Fss = the ratio of the diffusion radiation on the tilted surface of the collector to the horizontal ones
    S = Angulo entre o coletor e a horizontal
    IrTeta = Radiação Refletida
    krTeta = Coeficiente de reflexão do meio
    Fsg = 'related to the collector angle'
    It = Radiação Total
    phi = Angulo de Azimute (Não sei explicar isso)
"""

# O artigo usa todos os angulos em graus

import numpy as np

i = 1
Hr = 1
Ic = 1
B = 1
C = 1
S = 1
krTeta = 1
DiaDoAno = 1

def eq_delta():
    return 23.45*np.sin(np.deg2rad(360 * (284 + DiaDoAno)/365))
delta = eq_delta()

def eq_Beta():
    return np.arcsin(np.cos(i) * np.cos(Hr) * np.cos(delta) + np.sin(i) * np.sin(delta))
beta = eq_Beta()

def eq_Idn():
    return Ic*np.exp((-1)*B/np.sin(beta))
Idn = eq_Idn()

def eq_Fss():
    return (1 + np.cos(S))/2
Fss = eq_Fss()

def eq_IdTeta():
    return C * Idn * Fss
IdTeta = eq_IdTeta()

def eq_Fsg():
    return (1 - np.cos(S))/2
Fsg = eq_Fsg()

def eq_IrTeta():
    return (Idn + IdTeta) * krTeta * Fsg
IrTeta = eq_IrTeta()
    
def eq_phi():
    return np.arcsin(np.cos(delta) * np.sin(Hr) / np.cos(beta))
phi = eq_phi()

def eq_Teta():
    return np.cos(beta) * np.cos(phi) * np.sin(S) + np.sin(beta) * np.cos(S)
teta = eq_Teta()
    
def eq_It():
    return Idn*np.cos(teta) + IrTeta + IdTeta
It = eq_It()