# IMN projekt 1 Przemyslaw_Rewis 22.10.2021
# zadanie 4. RRZ 2 rzędu
import numpy as np

# dane parametry

rezystancja = 100  # R
indukcyjnosc = 0.1  # L
pojemnosc = 0.001  # C
omega_0 = 1.0 / np.sqrt(indukcyjnosc * pojemnosc)
T_0 = 2 * np.pi / omega_0
t_lewy_kraniec = 0.0
t_prawy_kraniec = 4 * T_0
delt = 0.0001
n = int((t_prawy_kraniec - t_lewy_kraniec) / delt + 1)


def potencjal_zrodla_napiecia(t, omega_v):
    krok_czasowy = float(t)
    czestosc_zrodla = float(omega_v)
    return 10 * np.sin(krok_czasowy * czestosc_zrodla)


# warunki poczatkowe

Q_0 = 0
I_0 = 0

if __name__ == '__main__':
    czestosci_zrodla = np.array([0.5, 0.8, 1.0, 1.2])
    for omega in czestosci_zrodla:
        numeryczne = open(f'./z4rozwiazanie_numeryczne{omega}.csv', 'w+')
        krok = np.empty(n, dtype=float)
        # tablice przechowujace wartosci
        wynik_Q = np.empty(n, dtype=float)
        wynik_I = np.empty(n, dtype=float)
        wynik_V = np.empty(n, dtype=float)
        aktualna_czestosc = omega * omega_0
        krok[0] = 0.0
        wynik_V[0] = potencjal_zrodla_napiecia(0.0, aktualna_czestosc)
        wynik_I[0] = I_0
        wynik_Q[0] = Q_0
        numeryczne.write(f'{krok[0]:.9f},{wynik_Q[0]:.9f},{wynik_I[0]:.9f}\n')
        for i in np.arange(1, n):
            krok[i] = i * delt
            k_Q1 = wynik_Q[i - 1]
            k_I1 = (wynik_V[i - 1] / indukcyjnosc) - (1.0 / (indukcyjnosc * pojemnosc) * wynik_Q[i - 1]) - (
                        rezystancja / indukcyjnosc * wynik_I[i - 1])
            k_Q2 = wynik_I[i - 1] + delt / 2.0 * k_I1
            k_I2 = (wynik_V[i - 1] / indukcyjnosc) - (
                        1.0 / (indukcyjnosc * pojemnosc) * (wynik_Q[i - 1] + delt / 2.0 * k_Q1)) - (
                           rezystancja / indukcyjnosc * (wynik_I[i - 1] + delt / 2.0 * k_I1))
            k_Q3 = wynik_I[i - 1] + delt / 2.0 * k_I2
            k_I3 = (wynik_V[i - 1] / indukcyjnosc) - (
                        1.0 / (indukcyjnosc * pojemnosc) * (wynik_Q[i - 1] + delt / 2.0 * k_Q2)) - (
                           rezystancja / indukcyjnosc * (wynik_I[i - 1] + delt / 2.0 * k_I2))
            k_Q4 = wynik_I[i - 1] + delt * k_I3
            k_I4 = (wynik_V[i - 1] / indukcyjnosc) - (
                        1.0 / (indukcyjnosc * pojemnosc) * (wynik_Q[i - 1] + delt * k_Q3)) - (
                           rezystancja / indukcyjnosc * (wynik_I[i - 1] + delt * k_I3))

            wynik_Q[i] = wynik_Q[i - 1] + delt / 6.0 * (k_Q1 + 2 * k_Q2 + 2 * k_Q3 + k_Q4)
            wynik_I[i] = wynik_I[i - 1] + delt / 6.0 * (k_I1 + 2 * k_I2 + 2 * k_I3 + k_I4)
            wynik_V[i] = potencjal_zrodla_napiecia(delt * i, aktualna_czestosc)

            numeryczne.write(f'{krok[i]:.9f},{wynik_Q[i]:.9f},{wynik_I[i]:.9f}\n')
