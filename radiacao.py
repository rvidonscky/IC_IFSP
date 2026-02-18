# radiacao.py - completo (salvar para importação nas eq 11, 12 e 13)

import math
import numpy as np

def deg2rad(x): return x * math.pi / 180.0
def rad2deg(x): return x * 180.0 / math.pi

def SolarDeclination(day_of_year):
    deg = 23.45 * math.sin(deg2rad(360.0 * (284 + day_of_year) / 365.0))
    return deg2rad(deg)

def hour_angle(hora_decimal):
    return deg2rad((hora_decimal - 12.0) * 15.0)

def sin_solar_altitude(latitude_deg, decl_rad, hour_ang_rad):
    i = deg2rad(latitude_deg)
    return math.cos(i) * math.cos(hour_ang_rad) * math.cos(decl_rad) + math.sin(i) * math.sin(decl_rad)

def direct_radiation(Ic, B, sinb):
    if sinb <= 0.0: return 0.0
    return Ic * math.exp(-B / sinb)

def FSS(tilt_deg):
    S = deg2rad(tilt_deg)
    return (1.0 + math.cos(S)) / 2.0

def diffuse_radiation(IDN, tilt_deg):
    C = 0.135
    return C * IDN * FSS(tilt_deg)

def Fsg(tilt_deg):
    S = deg2rad(tilt_deg)
    return (1.0 - math.cos(S)) / 2.0

def reflected_radiation(IDN, Idq, tilt_deg):
    krq = 0.2
    return (IDN + Idq) * krq * Fsg(tilt_deg)

def cos_solar_azimuth(latitude_deg, decl_rad, sinb):
    i = deg2rad(latitude_deg)
    numerador = (sinb - math.sin(i) * math.sin(decl_rad))
    denominador = math.cos(i) * math.cos(decl_rad)
    if abs(denominador) < 1e-8: return 1.0
    cosf = numerador / denominador
    return max(-1.0, min(1.0, cosf))

def cos_incidence_full(latitude_deg, decl_rad, tilt_deg, hour_ang_rad, cosf):
    i = deg2rad(latitude_deg)
    S = deg2rad(tilt_deg)
    sinf = math.sqrt(max(0.0, 1.0 - cosf**2))
    term1 = math.sin(decl_rad) * math.sin(i) * math.cos(S)
    term2 = -math.sin(decl_rad) * math.cos(i) * math.sin(S) * cosf
    term3 = math.cos(decl_rad) * math.cos(i) * math.cos(S) * math.cos(hour_ang_rad)
    term4 = math.cos(decl_rad) * math.sin(i) * math.sin(S) * cosf * math.cos(hour_ang_rad)
    term5 = math.cos(decl_rad) * math.sin(S) * sinf * math.sin(hour_ang_rad)
    return term1 + term2 + term3 + term4 + term5

def total_radiation(IDN, Idq, Irq, cosq):
    if cosq < 0: cosq = 0.0
    return IDN * cosq + Idq + Irq

def compute_It_series(time_seconds,
                      latitude_deg=32,
                      tilt_deg=37,
                      day_of_year=200,
                      Ic=1357.0,
                      B=0.21):
    
    times = np.asarray(time_seconds)
    # converte os segundos para hora decimal (0-24) por exemplo, 43200 segundos vira 12 
    hours = (times / 3600.0)  # 0..24
    decl = SolarDeclination(day_of_year)
    It_series = np.zeros_like(hours)

    for idx, h in enumerate(hours):
        # se estiver fora do intervalo do dia o sinb dará <=0 e já retorna 0
        Hr = hour_angle(h)
        sinb = sin_solar_altitude(latitude_deg, decl, Hr)
        if sinb <= 0:
            It_series[idx] = 0.0
            continue
        IDN = direct_radiation(Ic, B, sinb)
        Idq = diffuse_radiation(IDN, tilt_deg)
        Irq = reflected_radiation(IDN, Idq, tilt_deg)
        cosf = cos_solar_azimuth(latitude_deg, decl, sinb)
        cosq = cos_incidence_full(latitude_deg, decl, tilt_deg, Hr, cosf)
        It = total_radiation(IDN, Idq, Irq, cosq)
        It_series[idx] = It

    return It_series
