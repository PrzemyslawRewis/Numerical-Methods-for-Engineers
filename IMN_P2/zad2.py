# IMN projekt 2 Przemyslaw_Rewis 2.11.2021
# zadanie 2.  Metoda trapezow z iteracja Newtona (wzór 13)
import numpy as np

# zadane parametry
Beta = 0.001
N = 500.0
y = 0.1
t_lewy_kraniec = 0.0
t_prawy_kraniec = 100.0
dt = 0.1
mi0 = 1.0
TOL = 0.000001  # tolerancja
Maksymalna_liczba_itreacji = 20
alfa = Beta * N - y


# pomocnicza funkcja liczaca liczbe krokow param dt-krok czasowy

def ile_krokow(dt):
    argument = float(dt)
    return int((t_prawy_kraniec - t_lewy_kraniec) / argument + 1)


# dana funkcja f(t,u) wzor 2
def f(t, u):
    argument = float(t)
    argument2 = float(u)
    return (Beta * N - y) * argument2 - Beta * argument2 ** 2


# pomocnicza funkcja do sprawdzania warunku stopu (wynik<tolerancja) zwraca bool

def stop(obecny, poprzedni):
    return np.abs(obecny - poprzedni) < TOL


# funkcja implementujaca wzor 10
def newton(x, y):
    u_n = float(x)
    u_nu = float(y)
    return u_nu - (u_nu - u_n - dt / 2.0 * ((alfa * u_n - Beta * u_n ** 2) + (alfa * u_nu - Beta * u_nu ** 2))) / (
                1 - dt / 2.0 * (alfa - 2.0 * Beta * u_nu))


# funkcja do iteracyjnego poprawiania rozwiazania

def popraw(x):
    u_nu = float(x)
    counter = 0
    while True:
        counter = counter + 1
        u_n = u_nu
        u_nu = newton(x, u_n)
        if stop(u_nu, u_n) or counter > Maksymalna_liczba_itreacji:
            break
    return u_nu


# main
if __name__ == '__main__':
    n = ile_krokow(dt)
    plik_wynikowy = open(f'./newton.csv', 'w+')
    z = np.empty(n, dtype=float)
    u = np.empty(n, dtype=float)
    u[0] = mi0
    z[0] = N - u[0]
    plik_wynikowy.write(f'{0:.9f},{u[0]:.9f},{z[0]:.9f}\n')
    for i in np.arange(1, n):
        u_nu = popraw(u[i - 1])
        u[i] = u[i - 1] + dt / 2.0 * (f((i - 1) * dt, u[i - 1]) + f(i * dt, u_nu))  # 9
        z[i] = N - u[i]
        plik_wynikowy.write(f'{i * dt:.9f},{u[i]:.9f},{z[i]:.9f}\n')
