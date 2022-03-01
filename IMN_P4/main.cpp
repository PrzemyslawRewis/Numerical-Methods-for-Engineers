// IMN projekt 4
//Przemyslaw_Rewis 16.11.2021
//Równanie Poissona - relaksacja globalna i lokalna.
#include <cstdio>
#include <cstdlib>
#include <cmath>

//definicje stalych jako dyrektywy preprocesora w celu zwiększenia wydajności
#define epsilon 1.0
#define delta 0.1
#define nx 150
#define ny 100
#define V1 10.0
#define V2 0.0
#define xmax delta * nx
#define ymax delta * ny
#define sigmax 0.1 * xmax
#define sigmay 0.1 * ymax
#define TOL 1e-8

inline double Ro(double i, double j);
void metodaRelaksacjiGlobalnej(double wg, FILE *plik1,FILE *plik2,FILE *plik3);
void metodaRelaksacjiLokalnej(double wl, FILE *plik);


int main()
{
    FILE *globalna1;
    FILE *globalna2;
    FILE *globalna3;
    FILE *globalna4;
    FILE *globalna5;
    FILE *globalna6;
    FILE *lokalna1;
    FILE *lokalna2;
    FILE *lokalna3;
    FILE *lokalna4;
    globalna1 = fopen("globalna1_blad.txt", "w");
	globalna2 = fopen("globalna2_blad.txt", "w");
	globalna3 = fopen("globalna1.txt", "w");
	globalna4 = fopen("globalna2.txt", "w");
    globalna5 = fopen("globalna1_potencjal.txt", "w");
	globalna6 = fopen("globalna2_potencjal.txt", "w");
    metodaRelaksacjiGlobalnej(0.6, globalna3,globalna1,globalna5);
    metodaRelaksacjiGlobalnej(1.0, globalna4,globalna2,globalna6);
    lokalna1 = fopen("lokalna1.txt", "w");
    lokalna2 = fopen("lokalna2.txt", "w");
    lokalna3 = fopen("lokalna3.txt", "w");
    lokalna4 = fopen("lokalna4.txt", "w");
    metodaRelaksacjiLokalnej(1.0, lokalna1);
    metodaRelaksacjiLokalnej(1.4, lokalna2);
    metodaRelaksacjiLokalnej(1.8, lokalna3);
    metodaRelaksacjiLokalnej(1.9, lokalna4); 
    return 0;
}
//potegowanie wykonane jako operacja mnozenia wzory 19 i 20
inline double Ro(double i, double j)
{
    double p1 = exp( - (((i * delta - 0.35 * xmax)*(i * delta - 0.35 * xmax))/(sigmax*sigmax)) - (((j * delta - 0.5 * ymax)*(j * delta - 0.5 * ymax))/ (sigmay*sigmay)));
    double p2 = -exp( - (((i * delta - 0.65 * xmax)*(i * delta - 0.65 * xmax))/(sigmax*sigmax)) - (((j * delta - 0.5 * ymax)*(j * delta - 0.5 * ymax))/(sigmay*sigmay)));
    return p1 + p2;
}

void metodaRelaksacjiGlobalnej(double wg, FILE *plik1,FILE *plik2,FILE *plik3)
{
    int iterator =0;
    double Vn[nx+1][ny+1];
    double Vs[nx+1][ny+1];
    double p=0.0;
    double S=0.0;
    double S1=0.0;
    //wypelnianie macierzy zerami
    for(int i = 0; i <= nx; i++)
    {
        for(int j = 0; j <= ny; j++) 
        {
            Vn[i][j] = 0.0;
            Vs[i][j] = 0.0;
        }
    }
    // wypelnienie wierszy V1 V2 rys1
    for(int i = 0; i <= nx; i++)
    {
            Vn[i][0] = V1;
            Vs[i][0] = V1;
            Vn[i][ny] = V2;
            Vs[i][ny] = V2;
    }
    do
    {
        for(int i = 1; i < nx; i++)
        {
            for(int j = 1; j < ny; j++) 
            {
                p = Ro(i,j);
                Vn[i][j] = (1.0/4.0) * (Vs[i+1][j] + Vs[i-1][j] + Vs[i][j+1] + Vs[i][j-1] + ((delta*delta)/epsilon)* p);//wzor 9
            }
        }
        for(int j = 1; j < ny; j++) 
        {
            Vn[0][j] =  Vn[1][j];// wzor 10
            Vn[nx][j] = Vn[nx-1][j];//wzor 11
        }

        for(int i = 0; i <= nx; i++)
        {
            for(int j = 1; j < ny; j++) 
            {
                Vs[i][j] = (1.0 - wg) * Vs[i][j] + wg * Vn[i][j];//wzor 12
            }
        }

        S1 = S;
        S = 0.0;
        p = 0.0;
        for(int i =0; i < nx; i++)
        {
            for(int j = 0; j < ny; j++)
            {
                p = Ro(i,j);
                S += pow(delta,2.0) * ( (1.0/2.0) * pow(((Vs[i+1][j] - Vs[i][j])/delta), 2.0) + (1.0/2.0) * pow(((Vs[i][j+1] - Vs[i][j])/delta), 2.0) - p * Vs[i][j]);//wzor 17
            }
        }
        fprintf(plik1, "%d,%g\n",iterator,S);
        iterator++;
    }while(fabs((S - S1)/S1) > TOL);
    double blad[nx+1][ny+1];
    for(int i = 0; i <= nx; i++)
    {
        for(int j = 0; j <= ny; j++) 
        {
            blad[i][j] = 0.0;
        }
    }
    for(int i=1; i<nx; i++)
    {
        for(int j=1; j<ny; j++)
        {
            p = Ro(i,j);
            blad[i][j]= ((Vn[i+1][j] - 2*Vn[i][j] + Vn[i-1][j])/(delta*delta) + (Vn[i][j+1] - 2*Vn[i][j] + Vn[i][j-1])/(delta*delta)) + (p/epsilon);
            fprintf(plik2, "%g\t",blad[i][j]);
            fprintf(plik3,"%g\t",Vn[i][j]);
        }
         fprintf(plik2, "\n");
         fprintf(plik3,"\n");
    }
}

void metodaRelaksacjiLokalnej(double wl, FILE *plik)
{
    int iterator = 0;
    double V[nx+1][ny+1];
    double p = 0.0;
    double S = 0.0;
    double S1 = 0.0;
    //wypelnianie macierzy zerami
    for(int i = 0; i <= nx; i++)
    {
        for(int j = 0; j <= ny; j++) 
        {
            V[i][j] = 0.0;
        }
    }
    // wypelnienie wierszy V1 V2 rys1
    for(int i = 0; i <= nx; i++)
    {
        V[i][0] = V1;
        V[i][ny] = V2;
    }
    do
    {
        for(int i = 1; i < nx; i++)
        {
            for(int j = 1; j < ny; j++)
            {
                p = Ro(i,j);
                V[i][j] = (1 - wl) * V[i][j] + (wl/4.0) * (V[i+1][j] + V[i-1][j] + V[i][j+1] + V[i][j-1] + ((delta*delta)/epsilon) * p);//wzor 13
            }
        }
            for(int j = 1; j < ny; j++) 
            {
                V[0][j] =  V[1][j];
                V[nx][j] = V[nx-1][j];
            }
        S1 = S;
        S = 0.0;
        p = 0.0;
        for(int i =0; i < nx; i++)
        {
            for(int j = 0; j < ny; j++)
            {
                p = Ro(i,j);
                S += (delta*delta) * ( (1.0/2.0) * pow(((V[i+1][j] - V[i][j])/delta), 2.0) + (1.0/2.0) * pow(((V[i][j+1] - V[i][j])/delta), 2.0) - p * V[i][j] );//wzor 17
            }
        }
        fprintf(plik, "%d,%g\n",iterator,S);
        iterator++;
    }while(fabs((S - S1)/S1) > TOL);
}