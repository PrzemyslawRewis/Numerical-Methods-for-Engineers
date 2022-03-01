# IMN projekt 10 Przemyslaw_Rewis 23.11.2021
# Rownanie falowe - symulacja drgan struny metoda Verleta.

import numpy as np
import wykres

# zadane parametry
nx = 150
nt = 1000
delta = 0.1
delta_t = 0.05
xA = 7.5
sigma = 0.5
t_max = delta_t * nt  # punkt 3.4
XF = 2.5  # wzor 3 punkt 3.4


# implementacja delty Kroneckera
def licz_delte_kroneckera(x):
    if (abs(x - XF) < 1e-8):  # potencjalny punkt blednego outputu 1e-3 moze byc zamiast 1e-3
        return 1
    else:
        return 0


# implementacja wymuszenia z punktu 1 wzory (2) (3)
def wymuszenie(x, t):
    delta_Kroneckera = licz_delte_kroneckera(x)
    return np.cos(50 * t / t_max) * delta_Kroneckera


# implementacja funkcji liczacej energie zgromadzona w strunie wzor 17
def licz_energie_w_strunie(u, v):
    E = 0.0
    prawyskladniksumy = 0.0
    lewyskladniksumy = delta / 4.0 * (((u[1] - u[0]) / delta) ** 2 + ((u[nx] - u[nx - 1]) / delta) ** 2)
    for i in np.arange(1, nx):
        prawyskladniksumy += v[i] ** 2 + ((u[i + 1] - u[i - 1]) / (2 * delta)) ** 2
    prawyskladniksumy *= delta / 2.0
    E = lewyskladniksumy + prawyskladniksumy
    return E


# funkcja zapisujaca wynik do pliku
def mapazapiszdopilku(u, t, plik):
    for i in np.arange(nx + 1):
        x = delta * i
        plik.write(f'{t},{x},{u[i]}\n')


# funkcja liczaca a wektor przyspieszenia wzor 11
def przyspieszenie(a, u, u0, t, alpha, betha):
    for i in np.arange(1, nx):
        x = delta * i
        a[i] = (u[i + 1] - 2 * u[i] + u[i - 1]) / (delta ** 2) - betha * (u[i] - u0[i]) / delta_t + alpha * wymuszenie(
            x, t)


# punkt 2.6 implementacja algorytmu
def varlette_schema(alfa, beta):
    plikE = open(f'Energia_alfa={alfa}_beta={beta}.csv', 'w')
    plikU = open(f'Mapau_alfa={alfa}_beta={beta}.csv', 'w')
    # tworzenie tablic 1D
    u0 = np.zeros(nx + 1)
    u = np.zeros(nx + 1)
    v = np.zeros(nx + 1)
    vp = np.zeros(nx + 1)
    a = np.zeros(nx + 1)
    # warunek poczatkowy punkt 2.4
    if alfa != 1.0:
        for i in np.arange(1, nx):
            x = delta * i
            u[i] = np.exp(-(x - xA) ** 2 / (2 * sigma ** 2))  # wzor 14
    u0 = u.copy()  # zachowanie poprzedniego wyniku
    # dalsza czesc algorytmu
    for n in np.arange(nt + 1):
        t = delta_t * n
        vp = v + delta_t / 2.0 * a
        u0 = u.copy()
        u = u + delta_t * vp
        przyspieszenie(a, u, u0, t, alfa, beta)
        v = vp + delta_t / 2.0 * a
        E = licz_energie_w_strunie(u, v)
        mapazapiszdopilku(u, t, plikU)
        plikE.write(f'{t},{E}\n')
    plikE.close()
    plikU.close()
    wykres.rysujmape(alfa, beta, f'Mapau_alfa={alfa}_beta={beta}.csv')


if __name__ == '__main__':
    varlette_schema(0.0, 0.0)  # punkt 3.3
    varlette_schema(0.0, 0.1)  # punkt 3.3
    varlette_schema(0.0, 1.0)  # punkt 3.3
    varlette_schema(1.0, 1.0)  # punkt 3.4
    wykres.rysujwykres()
