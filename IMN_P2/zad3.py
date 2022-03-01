# IMN projekt 2 Przemyslaw_Rewis 2.11.2021
# zadanie 3.  Niejawna metoda RK2 wzór korektora (17)
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

# tablica Butchera

a11 = 0.25
a12 = 0.25 - np.sqrt(3) / 6.0
a21 = 0.25 + np.sqrt(3) / 6.0
a22 = 0.25
b1 = 0.5
b2 = 0.5
c1 = 0.5 - np.sqrt(3) / 6.0
c2 = 0.5 + np.sqrt(3) / 6.0


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


# pomocnicze funkcje do zapisania równań predyktora jako układ równań nieliniowych

def f1(x, y, z):
    U1 = float(x)
    U2 = float(y)
    u_n = float(z)
    return U1 - u_n - dt * (a11 * (alfa * U1 - Beta * U1 ** 2) + a12 * (alfa * U2 - Beta * U2 ** 2))


def f2(x, y, z):
    U1 = float(x)
    U2 = float(y)
    u_n = float(z)
    return U2 - u_n - dt * (a21 * (alfa * U1 - Beta * U1 ** 2) + a22 * (alfa * U2 - Beta * U2 ** 2))


# main
if __name__ == '__main__':
    n = ile_krokow(dt)
    plik_wynikowy = open(f'./RK2.csv', 'w+')
    z = np.empty(n, dtype=float)
    u = np.empty(n, dtype=float)
    u[0] = mi0
    z[0] = N - u[0]
    plik_wynikowy.write(f'{0:.9f},{u[0]:.9f},{z[0]:.9f}\n')
    for i in np.arange(1, n):
        U1 = u[i - 1]
        U2 = u[i - 1]
        dU1 = 0.0
        dU2 = 0.0
        U1_u = U1
        U2_u = U2
        counter = 0
        while True:
            counter = counter + 1
            m11 = 1 - dt * a11*(alfa - 2 * Beta * U1)
            m12 = -dt * a12*(alfa - 2 * Beta * U2)
            m21 = -dt * a21*(alfa - 2 * Beta * U1)
            m22 = 1 - dt * a22*(alfa - 2 * Beta * U2)
            dU1 = (f2(U1_u, U2_u, u[i - 1]) * m12 - f1(U1_u, U2_u, u[i - 1]) * m22) / (m11 * m22 - m12 * m21)
            dU2 = (f1(U1_u, U2_u, u[i - 1]) * m21 - f2(U1_u, U2_u, u[i - 1]) * m11) / (m11 * m22 - m12 * m21)
            U1 = U1_u + dU1
            U2 = U2_u + dU2
            U1_u = U1
            U2_u = U2
            if stop(U1, U1_u) or counter > Maksymalna_liczba_itreacji:
                break
        u[i] = u[i - 1] + dt * (b1 * f((i - 1) * dt + c1 * dt, U1) + b2 * f((i - 1) * dt + c2 * dt, U2))  # wzor 17
        z[i] = N - u[i]

        plik_wynikowy.write(f'{i * dt:.9f},{u[i]:.9f},{z[i]:.9f}\n')
