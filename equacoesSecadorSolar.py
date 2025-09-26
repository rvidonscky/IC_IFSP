def SignificadoDasIncógnitas():
    """
    
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
    return

# O artigo usa todos os angulos em graus

import numpy as np

def eq_Beta(*, i, Hr, delta):
    return np.arcsin(np.cos(i) * np.cos(Hr) * np.cos(delta) + np.sin(i) * np.sin(delta))

def eq_Idn(*, Ic, B, Beta):
    return Ic*np.exp((-1)*B/np.sin(Beta))
    
def eq_delta(*, DiaDoAno):
    return 23.45*np.sin(np.deg2rad(360 * (284 + n)/365))
    
def eq_IdTeta(*, C, Idn, Fss):
    return C * Idn * Fss
    
def eq_Fss(*, S):
    return (1 + np.cos(S))/2
    
def eq_IrTeta(*, Idn, IdTeta, krTeta, Fsg):
    return (Idn + IdTeta) * krTeta * Fsg
    
def eq_Fsg(*, S):
    return (1 - np.cos(S))/2
    
def eq_It(*, Idn, IdTeta, IrTeta, teta):
    return Idn*np.cos(teta) + IrTeta + IdTeta
    
def eq_Teta(*, beta, phi, S):
    return np.cos(beta) * np.cos(phi) * np.sin(S) + np.sin(beta) * np.cos(s)
    
def eq_phi(*, delta, Hr, Beta):
    return np.arcsin(np.cos(delta) * np.sin(Hr) / cos(Beta))
    
