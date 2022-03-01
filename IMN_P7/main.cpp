// IMN projekt 7
// Przemyslaw_Rewis 23.11.2021
// Rownania Naviera-Stokesa - symulacja przeplywu cieczy lepkiej.
#include <cstdio>
#include <cstdlib>
#include <cmath>

//zadane parametry punkt 3.1 
const double delta = 0.01;
const double ro = 1.0;
const double mi = 1.0; 
const short int nx = 200;
const short int ny = 90;
const short int i_1 = 50;
const short int j_1 = 55;
const short int IT_MAX = 20000;

inline double y(int j){ return j * delta; } // wzor 5
void Wb_Psi(double Psi[nx+1][ny+1], double Qwy, double Qwe);
void Wb_Dzeta(double Dzeta[nx+1][ny+1], double Psi[nx+1][ny+1], double Qwy, double Qwe);
void RelaksacjaNS(double Qwe, FILE * plik1, FILE * plik2, FILE * plik3, FILE * plik4);


int main()
{
    FILE *plik1,*plik2,*plik3,*plik4,*plik5,*plik6,*plik7,*plik8,*plik9,*plik10,*plik11,*plik12;
    plik1 = fopen("Relaksacja_psi_Q-1000.txt", "w");
    plik2 = fopen("Relaksacja_dzeta_Q-1000.txt", "w");
    plik3 = fopen("Relaksacja_u_Q-1000.txt", "w");
    plik4 = fopen("Relaksacja_v_Q-1000.txt", "w");
    plik5 = fopen("Relaksacja_psi_Q-4000.txt", "w");
    plik6 = fopen("Relaksacja_dzeta_Q-4000.txt", "w");
    plik7 = fopen("Relaksacja_u_Q-4000.txt", "w");
    plik8 = fopen("Relaksacja_v_Q-4000.txt", "w");
    plik9 = fopen("Relaksacja_psi_Q4000.txt", "w");
    plik10 = fopen("Relaksacja_dzeta_Q4000.txt", "w");
    plik11 = fopen("Relaksacja_u_Q4000.txt", "w");
    plik12 = fopen("Relaksacja_v_Q4000.txt", "w");
    RelaksacjaNS(-1000.0, plik1, plik2, plik3, plik4);// punkt 3.4
    RelaksacjaNS(-4000.0, plik5, plik6, plik7, plik8);// punkt 3.5
    RelaksacjaNS(4000.0, plik9, plik10, plik11, plik12);// punkt 3.6
    return 0;
}

// punkt 1.2.1
void Wb_Psi(double Psi[nx+1][ny+1], double Qwy, double Qwe)
{
    //brzeg A (wejscie)
    for(int j=j_1; j<=ny; j++)
    {
        Psi[0][j] = (Qwe/(2.0*mi))*((pow(y(j),3)/3.0)-(pow(y(j),2)/2.0)*(y(j_1)+y(ny)) + y(j)*y(j_1)*y(ny));//wzor 13
    }
    //brzeg C (wyjscie)
    for(int j=0; j<=ny; j++)
    {
        Psi[nx][j] = (Qwy/(2.0*mi))*((pow(y(j), 3)/3.0) - (pow(y(j), 2)/2.0)* y(ny))+(Qwe*pow(y(j_1), 2)*(-y(j_1)+3.0*y(ny)) )/(12.0 * mi);//wzor 14
    }
    //brzeg B
    for(int i=1; i<=nx-1; i++)
    {
        Psi[i][ny]=Psi[0][ny];//wzor 15
    }
    //brzeg D
    for(int i=i_1; i<=nx-1; i++)
    {
        Psi[i][0]=Psi[0][j_1];//wzor 16
    }
    //brzeg E
    for(int j=1; j<=j_1; j++)
    {
        Psi[i_1][j]=Psi[0][j_1];//wzor 17
    }
    //brzeg F
    for(int i=1; i<=i_1; i++)
    {
        Psi[i][j_1]=Psi[0][j_1];//wzor 18
    }
}
// punkt 1.2.2
void Wb_Dzeta(double Dzeta[nx+1][ny+1], double Psi[nx+1][ny+1], double Qwy, double Qwe)
{
    // brzeg A (wejscie)
    for(int j=j_1; j<=ny; j++)
    {
        Dzeta[0][j]=(Qwe/(2.0*mi))*(2.0*y(j) - y(j_1) -y(ny));//wzor 19
    }
    // brzeg C (wyjscie)
    for(int j=0; j<=ny; j++)
    {
        Dzeta[nx][j]=(Qwy/(2.0*mi))*(2.0*y(j) - y(ny));//wzor 20
    }
    // brzeg B
    for(int i=1; i<=nx-1; i++)
    {
        Dzeta[i][ny] = (2.0/(delta*delta))*(Psi[i][ny-1] - Psi[i][ny]);//wzor 21
    }
    //brzeg D
    for(int i = i_1+1; i<=nx-1; i++)
    {
        Dzeta[i][0] = (2.0/(delta*delta))*(Psi[i][1] - Psi[i][0]);//wzor 22
    }
    //brzeg E
    for(int j=1; j<=j_1-1; j++)
    {
        Dzeta[i_1][j] = (2.0/(delta*delta))*(Psi[i_1+1][j] - Psi[i_1][j]);//wzor 23
    }
    //brzeg F
    for(int i=1; i<=i_1; i++)
    {
        Dzeta[i][j_1] = (2.0/(delta*delta))*(Psi[i][j_1+1] - Psi[i][j_1]);//wzor 24
    }
    //wierzcholek E/F
    Dzeta[i_1][j_1] = (0.5)*(Dzeta[i_1-1][j_1]+Dzeta[i_1][j_1-1]);//wzor 25
}

void RelaksacjaNS(double Qwe, FILE * plik1, FILE * plik2, FILE * plik3, FILE * plik4)
{
    double Psi[nx+1][ny+1];
    double Dzeta[nx+1][ny+1];
    double u[nx+1][ny+1]; //[V=(u,v)] punkt 1.2
    double v[nx+1][ny+1];
    //inicjalizacja tablic zerami
    for(int i=0; i<=nx; i++)
    {
        for(int j=0; j<=ny; j++)
        {
            Dzeta[i][j]=0.0;
            Psi[i][j]=0.0;
            u[i][j]=0.0;
            v[i][j]=0.0;
        }
    }
    // punkt 1.2
    double Qwy = Qwe*((pow(y(ny), 3) - pow(y(j_1), 3) - 3.0*y(j_1)*pow(y(ny),2) + 3.0* pow(y(j_1),2)*y(ny))/(pow(y(ny),3))); //wzor 12
    double omega=0.0;
    // implementacja algorytmu relaksacji rownan NS
    Wb_Psi(Psi,Qwy,Qwe);//ustalamy WB

    for(int it=1; it<IT_MAX; it++)
    {
        omega = (it < 2000) ? 0.0 : 1.0;
        for(int i =1; i<=nx-1; i++)
        {
            for(int j=1; j<=ny-1; j++)
            {
                if( i>i_1 || j>j_1)//(i,j)!=BRZEG
                {
                    Psi[i][j] = (0.25) * (Psi[i+1][j] + Psi[i-1][j] + Psi[i][j+1] + Psi[i][j-1]- (delta*delta)*Dzeta[i][j]);// wzor 8
                    Dzeta[i][j] = (0.25) * (Dzeta[i+1][j] + Dzeta[i-1][j] + Dzeta[i][j+1]+Dzeta[i][j-1]) - omega* (ro/(16.0*mi))*((Psi[i][j+1] - Psi[i][j-1])*(Dzeta[i+1][j] -Dzeta[i-1][j]) - (Psi[i+1][j] - Psi[i-1][j])*(Dzeta[i][j+1]-Dzeta[i][j-1]));
                    u[i][j] =(Psi[i][j+1] - Psi[i][j-1])/(2.0*delta); // punkt 1.1 zastapienie pochodnej ilorazem roznicowym
                    v[i][j] = -(Psi[i+1][j] - Psi[i-1][j])/(2.0*delta); // punkt 1.1 zastapienie pochodnej ilorazem roznicowym
                }
                else
                {
                    u[i][j] = 0.0;
                    v[i][j] = 0.0;
                }
            }
        }
        Wb_Dzeta(Dzeta,Psi,Qwy,Qwe);//modyfikacja WB
        //kontrola bledu gamma
        double gama =0.0;
        int j2 =j_1+2;
        for(int i =1; i<=nx-1; i++)
        {
            gama += (Psi[i+1][j2] + Psi[i-1][j2] + Psi[i][j2+1] + Psi[i][j2-1] - 4.0*Psi[i][j2] - (delta*delta)*Dzeta[i][j2]);//wzor 27
        }
    }
    //zapis wynikow
    for(int i=0; i<=nx; i++)
    {
        for(int j=0; j<=ny; j++)
        {
            fprintf(plik1,"%g\t", Psi[i][j]);
            fprintf(plik2,"%g\t", Dzeta[i][j]);
            fprintf(plik3,"%g\t", u[i][j]);
            fprintf(plik4,"%g\t", v[i][j]);
        }
        fprintf(plik1,"\n");
        fprintf(plik2,"\n");
        fprintf(plik3,"\n");
        fprintf(plik4,"\n");
    }
}