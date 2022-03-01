// IMN projekt 5
// Przemyslaw_Rewis 17.11.2021
// Równanie Poissona - relaksacja wielosiatkowa.
#include <cstdio>
#include <cstdlib>
#include <cmath>

//zadane parametry
const double delta = 0.2;
const int nx = 128;
const int ny = 128;
const double xmax = delta*nx;
const double ymax = delta*ny;
const double TOL = 1e-8;

void relaksacjaWielosiatkowa(FILE *plik,FILE *plik2,FILE *plik3,FILE *plik4,FILE *plik5,FILE *plik6);

int main()
{
    FILE *plik,*plik2,*plik3,*plik4,*plik5,*plik6;
    plik  = fopen("relaksacjaWielosiatkowa.txt", "w");
    plik2 = fopen("relaksacjaWielosiatkowa_k16.txt", "w");
    plik3 = fopen("relaksacjaWielosiatkowa_k8.txt", "w");
    plik4 = fopen("relaksacjaWielosiatkowa_k4.txt", "w");
    plik5 = fopen("relaksacjaWielosiatkowa_k2.txt", "w");
    plik6 = fopen("relaksacjaWielosiatkowa_k1.txt", "w");
    relaksacjaWielosiatkowa(plik, plik2, plik3, plik4, plik5, plik6);
    return 0;
}

void relaksacjaWielosiatkowa(FILE *plik,FILE *plik2,FILE *plik3,FILE *plik4,FILE *plik5,FILE *plik6)
{
    double S = 0.0; 
    double Spoprzednie = 0.0;
    int iterator = 0, k;
    double V[nx+1][ny+1];

    // warunki brzegowe Dirichleta
    for(int i = 0; i <= nx; i++)
    {     
        V[i][ny] = -sin(2.0* M_PI * (i * delta / xmax));//wzor 18 + rys1
        V[i][0] = sin(2.0* M_PI * (i * delta / xmax));//wzor 20 + rys1
    }
    for(int i = 0; i <= ny; i++)
    {
        V[0][i] = sin(M_PI * (i * delta / ymax));//wzor 17 + rys1
        V[nx][i] = sin(M_PI * (i * delta / ymax));//wzor 19 + rys1
    }
    //wypelniamy macierz zerami
    for( int i=1; i<nx; i++ )
    {
        for( int j=1; j<ny; j++ )
        {
            V[i][j] = 0.0;
        }
    }
    // petla zmieniajaca k w sposob k=16,k=8,k=4,k=2,k=1
    for(k = 16; k>=1; k/=2)
    {
        do
        {
            for(int i = k; i <= nx-k; i+=k)
            {
                for(int j = k; j<= ny-k; j+=k)
                {
                    V[i][j] = (1.0/4.0) * (V[i+k][j] + V[i-k][j] + V[i][j+k] + V[i][j-k]);// wzor 7 zmodyfikowany podstawowy przepis relaksacji
                }
            }
            Spoprzednie = S;
            S=0.0;
            //implementacja wzoru 15
            for(int i =0; i <= nx - k; i+=k)
            {
                for(int j = 0; j <= ny - k; j+=k)
                {
                    S += (pow(k*delta,2.0)/2.0) * (pow(((V[i+k][j]-V[i][j])/(2.0*k*delta) + (V[i+k][j+k] -V[i][j+k])/(2.0*k*delta)) , 2.0) + pow(((V[i][j+k]-V[i][j])/(2.0*k*delta) + (V[i+k][j+k] -V[i+k][j])/(2.0*k*delta)) , 2.0) );
                }
            }
            fprintf(plik, "%d\t%g\n",iterator,S);
            iterator++;
        }while(fabs((S - Spoprzednie)/Spoprzednie) > TOL);//sprawdzanie warunku stopu wzor 16
        fprintf(plik, "\n");
        // Proces zmiany siatki i następującej po niej relaksacji powtarzamy aż do uzyskania rozwiązania na najgęstszej siatce (k = 1). rys 1 czerwone
        if( k != 1 )
        {
            for(int i = 0; i <= nx-k; i+=k)
            {
                for(int j = 0; j <= ny-k; j+=k)
                {
                    V[i+k/2][j+k/2] = 1.0/4.0 * (V[i][j]   + V[i+k][j] + V[i][j+k] + V[i+k][j+k]);// wzor 8
                    if( i+k != nx )
                        V[i+k][j+k/2] = 1.0/2.0 * (V[i+k][j] + V[i+k][j+k]);// wzor 9
                    if( j+k != ny )
                        V[i+k/2][j+k]   = 1.0/2.0 * (V[i][j+k] + V[i+k][j+k]);// wzor 10
                    if (j != 0)
                        V[i+k/2][j] = 1.0/2.0  * (V[i][j] + V[i+k][j]);// wzor 11
                    if (i != 0)
                        V[i][j+k/2] = 1.0/2.0  * (V[i][j] + V[i][j+k]);// wzor 12
                }
            }
        }
        //zapis wynikow
        for(int i=0; i <=nx; i+=k)
        {
            for(int j=0; j<=ny; j+=k)
            {
                switch (k)
                {
                    case 16:
                        fprintf(plik2, "%g\t",V[i][j]);
                        break;
                    case 8:
                        fprintf(plik3, "%g\t",V[i][j]);
                        break;
                    case 4:
                        fprintf(plik4, "%g\t",V[i][j]);
                        break;  
                    case 2:
                        fprintf(plik5, "%g\t",V[i][j]);
                        break; 
                    case 1:
                        fprintf(plik6, "%g\t",V[i][j]);
                        break;             
                    default:
                        break;
                }
            }
            switch (k)
            {
                case 16:
                    fprintf(plik2, "\n");
                    break;
                case 8:
                    fprintf(plik3, "\n");
                    break;
                case 4:
                    fprintf(plik4, "\n");
                    break;  
                case 2:
                    fprintf(plik5, "\n");
                    break; 
                case 1:
                    fprintf(plik6, "\n");
                    break;             
                default:
                    break;
            }
        }
    }
}