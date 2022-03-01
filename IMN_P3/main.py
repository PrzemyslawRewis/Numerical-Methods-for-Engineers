# IMN projekt 3 Przemyslaw_Rewis 13.11.2021
# RRZ - kontrola kroku czasowego w problemach sztywnych.


# implementacja wzoru 2
def f(v):
    arg = float(v)
    return arg


# implementacja wzoru 3
def g(x, v, alfa):
    xarg = float(x)
    varg = float(v)
    alfaarg = float(alfa)
    return alfaarg * (1.0 - x ** 2) * varg - xarg


# implementacja wzoru 11
def F(xn, vn, xn_poprzedni, vn_poprzedni, delt, alfa):
    wynik = float(xn - xn_poprzedni - delt / 2.0 * (f(vn_poprzedni) + f(vn)))
    return wynik


# implementacja wzoru 12
def G(xn, vn, xn_poprzedni, vn_poprzedni, delt, alfa):
    wynik = float(vn - vn_poprzedni - delt / 2.0 * (g(xn_poprzedni, vn_poprzedni, alfa) + g(xn, vn, alfa)))
    return wynik


# implementacja metody trapezow punkt 2.1
def mtrapezow(xn_poprzedni, vn_poprzedni, delt, alfa, delta=1e-10):
    xn = xn_poprzedni
    vn = vn_poprzedni
    while True:
        a11 = 1.0
        a12 = -delt / 2.0
        a21 = -delt / 2.0 * (-2.0 * alfa * xn * vn - 1.0)
        a22 = 1.0 - delt / 2.0 * alfa * (1.0 - xn ** 2)
        dx = (-F(xn, vn, xn_poprzedni, vn_poprzedni, delt, alfa) * a22 + G(xn, vn, xn_poprzedni, vn_poprzedni, delt,
                                                                           alfa) * a12) / (a11 * a22 - a12 * a21)
        dv = (-G(xn, vn, xn_poprzedni, vn_poprzedni, delt, alfa) * a11 + F(xn, vn, xn_poprzedni, vn_poprzedni, delt,
                                                                           alfa) * a21) / (a11 * a22 - a12 * a21)
        xn += dx
        vn += dv
        if abs(dx) < delta and abs(dv) < delta:
            break

    return xn, vn


# implementacja metody rk2 punkt2.2
def rk2(xn, vn, delt, alfa):
    k1x = vn
    k1v = alfa * (1.0 - xn ** 2) * vn - xn
    k2x = vn + delt * k1v
    k2v = alfa * (1.0 - (xn + delt * k1x) ** 2) * (vn + delt * k1v) - (xn + delt * k1x)
    xn_nastepny = float(xn + delt / 2.0 * (k1x + k2x))
    vn_nastepny = float(vn + delt / 2.0 * (k1v + k2v))
    return xn_nastepny, vn_nastepny


# implementacja ogolnego algorytmu numerycznego rozwiazywania rownania rozniczkowego z doborem kroku czasowego
def krok_czasowy(schematnumeryczny, TOL, nazwaplikuwynikowego):
    xn = 0.01
    vn = 0
    delt = 1
    S = 0.75
    p = 2
    t_max = 40
    alfa = 5
    t = 0
    plik = open(nazwaplikuwynikowego, "w+")

    while True:
        xn1, vn1 = schematnumeryczny(xn, vn, delt, alfa)
        xn2_2, vn2_2 = schematnumeryczny(xn1, vn1, delt, alfa)
        xn2_1, vn2_1 = schematnumeryczny(xn, vn, 2 * delt, alfa)
        Ex = (xn2_2 - xn2_1) / (2 ** p - 1)
        Ev = (vn2_2 - vn2_1) / (2 ** p - 1)
        if max(abs(Ex), abs(Ev)) < TOL:
            t += 2 * delt
            xn = xn2_2
            vn = vn2_2
            plik.write(f'{t},{delt},{xn},{vn}\n')
        delt = ((S * TOL) / (max(abs(Ex), abs(Ev)))) ** (1 / (p + 1)) * delt
        if t > t_max:
            break
    plik.close()


# main
if __name__ == '__main__':
    krok_czasowy(mtrapezow, 1e-2, 'metoda_trapezow_1.csv')
    krok_czasowy(mtrapezow, 1e-5, 'metoda_trapezow_2.csv')
    krok_czasowy(rk2, 1e-2, 'rk2_1.csv')
    krok_czasowy(rk2, 1e-5, 'rk2_2.csv')
