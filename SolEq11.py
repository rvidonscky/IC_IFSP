#Entender e resolver as condições de borda e de contorno

#As letras e nomes das funções são baseadas no livro ""Numerical Partial Differential Equations: Finite Difference Methods""
#Paginas 5 ate 17

# A letra k se refere a quantidade de passos em 'x' e a letra n se refere à quantidade de passos em 't'
#Deu confusão

import EqAlSemParam as easp

#CONSTANTES
M = 10
X_MINIMO = 0
X_MAXIMO = 1000
kg = 1
rog = 1
Cg = 1
alfag = 1
deltag = 1
hgam = 1
C = kg / (rog * Cg)
TgDelta = 35 #Fig. 4

deltaX = 1/M
deltaT = 0.02

def F(x, t):
    return easp.It * alfag / (rog * Cg * deltag)

def f(x):
    return # Tg = Tp = Tam
    
def a(t):
    return hgam * (Tg0 - Ta)
    
def b(t):
    return hgam * (Tam - TgDelta)

def u(k, n + 1):
    if n == -1:
        return f(k*deltaX)
    elif k == 0:
        return a((n+1)*deltaT)
    elif k == M:
        return b((n+1)*deltaT)
    else:
        return u(k, n) + C * (u(k+1, n) - 2 * u(k, n) + u(k - 1, n)) * deltaT / (deltaX ** 2) + deltaT * F(k*deltaX, n * deltaT)