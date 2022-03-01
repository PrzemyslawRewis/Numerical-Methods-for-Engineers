// IMN projekt 8
// Wymagany minimalny standard jezyka -std=c++11
// Przemyslaw_Rewis 14.12.2021
// Rownanie adwekcji-dyfuzji – symulacja transportu masy metoda Cranka-Nicolson.

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <array>
#include <iostream>

//parametry punkt 3.1
constexpr short int nx=400;
constexpr short int ny=90;
constexpr short int i1=200;
constexpr short int i2=210;
constexpr short int j_1=50;
constexpr double delta=0.01;
constexpr double sigma=10*delta;
constexpr double xa=0.45;
constexpr double ya=0.45;
constexpr short int ITMAX=10000;

void solver(double D, FILE *plik1, FILE *plik2);

int main()
{
    FILE *plik1,*plik2,*plik3,*plik4;
    plik1 = fopen("c_xsr.txt","w");
    plik2 = fopen("mapa.txt","w");
    plik3 = fopen("c_xsr2.txt","w");
    plik4 = fopen("mapa2.txt","w");
    solver(0.0, plik1, plik2);
    solver(0.1, plik3, plik4);
    return 0;
}


void solver(double D, FILE *plik1, FILE *plik2)
{
    // tworzenie tablic
    std::array<std::array<double,ny+1>,nx+1> u0;
    std::array<std::array<double,ny+1>,nx+1> u1;
    std::array<std::array<double,ny+1>,nx+1> vx;
    std::array<std::array<double,ny+1>,nx+1> vy;
    std::array<std::array<double,ny+1>,nx+1> psi;

    FILE *p1,*p2;
    p1 = fopen("vx.txt","w");
    p2 = fopen("vy.txt","w");


    // wczytanie funkcji strumienia punkt 3.2
    FILE *funkcja;
    funkcja = fopen("psi.dat","r");
    int i, j;
    double ps;
    while(fscanf(funkcja, "%d%d%lf", &i, &j, &ps) == 3) 
    {
        psi[i][j]=ps;
    }
    //Pole predkosci punkt 3.3 i 2.4
    for(int i = 1; i<=nx-1; i++)
    {
        for(int j = 1; j<=ny-1; j++)
        {
            vx[i][j] = (psi[i][j+1] - psi[i][j-1])/(2.0*delta);//wzor 11
            vy[i][j] = - (psi[i+1][j] - psi[i-1][j])/(2.0*delta);//wzor 12
        }
    }
    //wzor 13 zastawka
    for(int i = i1; i <= i2; i++)
    {
        for(int j = 0; j <= j_1; j++)
        {
            vx[i][j] = 0.0;
            vy[i][j] = 0.0;
        }
    }
    //wzor 14 dolny i gorny brzeg
    for(int i = 1; i <= nx-1; i++)
    {
        vx[i][0] = 0.0;
        vy[i][ny] = 0.0;
    }
    // lewy i prawy brzeg przypisanie wartosci z poprzednich wezlow
    for(int j = 0; j <=ny; j++)
    {
        vx[0][j] = vx[1][j]; //wzor 15
        vx[nx][j] = vx[nx-1][j];//wzor 16
    }
    //zapis do pliku wynikow
    for(int i = 0; i<=nx; i++)
    {
        for(int j = 0; j<=ny; j++)
        {
            fprintf(p1,"%g\t",vx[i][j]);
            fprintf(p2,"%g\t",vy[i][j]);
        }
         fprintf(p1,"\n");
         fprintf(p2,"\n");
    }
    //v max szukamy punkt i,j dla ktorego moduł predkosci jest najwiekszy
    double v_max = 0.0;
    for(int i = 0; i <= nx; i++)
    {
        for(int j = 0; j <= ny; j++)
        {
            if(sqrt(pow(vx[i][j], 2) + pow(vy[i][j], 2)) > v_max)
            {
                v_max = sqrt(pow(vx[i][j], 2) + pow(vy[i][j], 2));//ptk 3.3
            }
                
        }
    }

    double delta_t = delta/(4.0*v_max);//wyznaczenie kroku czasowego wzor 18
    std::cout<<"krok czasowy"<<delta_t;

    //punkt 2.5 warunek poczatkowy inicjalizacja gestosci 
    for(int i = 0; i <= nx; i++)
    {
        for(int j = 0; j <= ny; j++)
        {
            u0[i][j] = 1.0/(2.0*M_PI*pow(sigma, 2)) * exp( -(pow((delta*i) - xa, 2) + pow((delta*j) - ya, 2))/(2.0*pow(sigma, 2)));//wzor 17
        }
    }
    // implementacja algorytmu punkt 2.6
    for(int it = 0; it <= ITMAX; it++)
    {
        for(int i = 0; i <= nx; i++)
        {
            for(int j = 0; j <= ny; j++)
            {
                u1[i][j] = u0[i][j];//inicjalizacja kolejnego kroku
            }
        }
        for(int k = 1; k <= 20; k++)
        {
            for(int i = 0; i <= nx; i++)
            {
                for(int j = 1; j <= ny-1; j++)
                {
                    if(i<i1 || i>i2 || j>j_1)//rysunek (i,j) E zastawka
                    {
                        if(i==0) 
                        {
                            u1[i][j] = (1.0/(1.0+((2.0*D*delta_t) / pow(delta, 2)))) * (u0[i][j] - (delta_t/2.0) * vx[i][j] *
                            (((u0[i+1][j] - u0[nx][j])/(2.0*delta)) + (u1[i+1][j] - u1[nx][j])/(2.0*delta)) - (delta_t / 2.0) * vy[i][j] *
                             ((u0[i][j+1] - u0[i][j-1])/(2.0*delta) + (u1[i][j+1] - u1[i][j-1])/(2.0*delta)) + (delta_t/2.0) * D *
                             ((u0[i+1][j] + u0[nx][j] + u0[i][j+1] + u0[i][j-1] - 4.0*u0[i][j])/pow(delta,2) + (u1[i+1][j] + u1[nx][j] + u1[i][j+1] + u1[i][j-1] )/pow(delta,2)));
                        }
                        else if(i==nx)
                        {
                            u1[i][j] = (1.0/(1.0+( (2.0*D*delta_t) / pow(delta, 2)))) * (u0[i][j] - (delta_t/2.0) * vx[i][j] *
                            (((u0[0][j] - u0[i-1][j])/(2.0*delta)) + (u1[0][j] - u1[i-1][j])/(2.0*delta)) - (delta_t / 2.0) * vy[i][j] *
                             ((u0[i][j+1] - u0[i][j-1])/(2.0*delta) + (u1[i][j+1] - u1[i][j-1])/(2.0*delta)) + (delta_t/2.0) * D *
                             ((u0[0][j] + u0[i-1][j] + u0[i][j+1] + u0[i][j-1] - 4.0*u0[i][j])/pow(delta,2) + (u1[0][j] + u1[i-1][j] + u1[i][j+1] + u1[i][j-1])/pow(delta,2)));
                        }
                        else
                        {
                            u1[i][j] = (1.0/(1.0+((2.0*D*delta_t )/ pow(delta, 2)))) * (u0[i][j] - (delta_t/2.0) * vx[i][j] *
                            (((u0[i+1][j] - u0[i-1][j])/(2.0*delta)) + (u1[i+1][j] - u1[i-1][j])/(2.0*delta)) - (delta_t / 2.0) * vy[i][j] *
                            ((u0[i][j+1] - u0[i][j-1] )/(2.0*delta) + (u1[i][j+1] - u1[i][j-1])/(2.0*delta)) + (delta_t/2.0) * D *
                            ((u0[i+1][j] + u0[i-1][j] + u0[i][j+1] + u0[i][j-1] - 4.0*u0[i][j])/pow(delta,2) + (u1[i+1][j] + u1[i-1][j] + u1[i][j+1] + u1[i][j-1])/pow(delta,2)));
                        }
                    }
                }
            }

        }

        for(int i = 0; i <= nx; i++)
        {
            for(int j = 0; j <= ny; j++)
            {
                u0[i][j] = u1[i][j];//zachowanie rozwiazania do nastepnego wywolania
            }
        }

        double c=0.0;
        double x_sr=0.0;
        for(int i = 0; i <= nx; i++)
        {
            for(int j = 0; j <=ny; j++)
            {
                c += u0[i][j];//wzor 19 lewy czynnik
                x_sr += (i*delta) * u0[i][j];//wzor 20 lewy czynnik
            }
        }
        c *= pow(delta, 2);//cd wzor 19 
		x_sr *= pow(delta, 2);//cd wzor 20 
        fprintf(plik1, "%g\t%g\t%g\n",delta_t*it, c, x_sr);//zapis wyniku
        int k=1;

        if(it%2000==0 && it<=10000)
        {
            for (int l = 0; l <= nx; l++) 
            {
                for (int m = 0;m <= ny; m++) 
                {
                    fprintf(plik2,"%g\t",u1[l][m]);
                }
                fprintf(plik2,"\n");
            }
            fprintf(plik2,"\n\n\n");
            k++;    
        }
        
    }
    fclose(p1);
    fclose(p2);
    fclose(funkcja);
}