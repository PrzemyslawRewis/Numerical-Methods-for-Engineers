# IMN projekt 1 Przemyslaw_Rewis 22.10.2021
# zadanie 3. Metoda jawna RK4
import numpy as np

# zadane parametry

lamb = -1.0
t_lewy_kraniec = 0.0
t_prawy_kraniec = 5.0


# dana funkcja y(t)
def y(t):
    argument = float(t)
    return np.exp(lamb * argument)


# pomocnicza funkcja liczaca liczbe krokow param delt-krok czasowy

def ile_krokow(dt):
    argument = float(dt)
    return int((t_prawy_kraniec - t_lewy_kraniec) / argument + 1)


# main
if __name__ == '__main__':
    kroki_czasowe = np.array([1.0, 0.1, 0.01])
    for delt in kroki_czasowe:
        n = ile_krokow(delt)
        krok = np.empty(n, dtype=float)  # indeks dolny przy y tablica
        y_arg = np.empty(n, dtype=float)  # wartosc rozwiazania tablica
        krok[0] = 0.0  # wartosc poczatkowa
        y_arg[0] = 1.0  # wartosc poczatkowa
        numeryczne = open(f'./z3rozwiazanie_numeryczne{delt}.csv', 'w+')  # w+ tryb odczyt zapis
        analityczne = open('./z3rozwiazanie_analityczne.csv', 'w+')
        numeryczne.write(f'{krok[0]:.9f},{y_arg[0]:.9f},{(y_arg[0] - y(0.0)):.9f}\n')
        # :.9f fromatowanie pola 9 pol po przecinku pola zapisane do pliku w kolumnach przyklad t0_01, y0_01, blad0_01
        analityczne.write(f'{krok[0]:.9f},{y(0.0):.9f}\n')
        for i in np.arange(1, n):

            k_1=lamb*y_arg[i-1]
            k_2=lamb*(y_arg[i-1]+delt/2*k_1)
            k_3=lamb*(y_arg[i-1]+delt/2*k_2)
            k_4=lamb*(y_arg[i-1]+delt*k_3)
            y_arg[i] = y_arg[i - 1] + delt/6 *(k_1+2*k_2+2*k_3+k_4)
            krok[i] = i * delt
            analityczne.write(f'{krok[i]:.9f},{y(i * delt):.9f}\n')
            numeryczne.write(f'{krok[i]:.9f},{y_arg[i]:.9f},{(y_arg[i] - y(i * delt)):.9f}\n')
