// IMN projekt 6
// Przemyslaw_Rewis 19.11.2021
// Równanie Poissona - rozwiązanie metodą algebraiczną.

#include "mgmres.c"

// dane stale
const double delta = 0.1;
const short int itr_max = 500;
const short int mr = 500;
const double tol_abs = 1e-8;
const double tol_rel = 1e-8;
static int control=1;

void poissonMetodaAlgebraiczna(int nx, int ny, double epsilon1, double epsilon2, double V1, double V2, double V3, double V4, FILE* plikmatrixA, FILE* plikvectorB, FILE* plikvectorV);

int main()
{
    FILE *plikMaczierzA, *plikWektorb, *plikWektorRozwiazanV, *plikWektorRozwiazanV_nx50, *plikWektorRozwiazanV_nx100, *plikWektorRozwiazanV_nx200, *plikWektorRozwiazanV_epsilony11, *plikWektorRozwiazanV_epsilony12, *plikWektorRozwiazanV_epsilony110;
    plikMaczierzA = fopen("macierzA.txt", "w");//punkt 2.3
    plikWektorb = fopen("vektorb.txt", "w");//punkt 2.3
    //punkt 2.1 
    plikWektorRozwiazanV = fopen("wektorRozwiazanV_nx4.txt", "w");
    poissonMetodaAlgebraiczna(4,4,1.0,1.0,10,-10,10,-10, plikMaczierzA, plikWektorb, plikWektorRozwiazanV);
    //punkt 2.5
    plikWektorRozwiazanV_nx50= fopen("wektorRozwiazanV_nx50.txt", "w");
    poissonMetodaAlgebraiczna(50,50,1.0,1.0,10,-10,10,-10, plikMaczierzA, plikWektorb, plikWektorRozwiazanV_nx50);
    plikWektorRozwiazanV_nx100= fopen("wektorRozwiazanV_nx100.txt", "w");
    poissonMetodaAlgebraiczna(100,100,1.0,1.0,10,-10,10,-10, plikMaczierzA, plikWektorb, plikWektorRozwiazanV_nx100);
    plikWektorRozwiazanV_nx200= fopen("wektorRozwiazanV_nx200.txt", "w");
    poissonMetodaAlgebraiczna(200,200,1.0,1.0,10,-10,10,-10, plikMaczierzA, plikWektorb, plikWektorRozwiazanV_nx200);
    //punkt 2.6
    plikWektorRozwiazanV_epsilony11= fopen("wektorRozwiazanV_epsilony11.txt", "w");
    poissonMetodaAlgebraiczna(100,100,1.0,1.0,0.0,0.0,0.0,0.0, plikMaczierzA, plikWektorb, plikWektorRozwiazanV_epsilony11);
    plikWektorRozwiazanV_epsilony12= fopen("wektorRozwiazanV_epsilony12.txt", "w");
    poissonMetodaAlgebraiczna(100,100,1.0,2.0,0,0,0,0, plikMaczierzA, plikWektorb, plikWektorRozwiazanV_epsilony12);
    plikWektorRozwiazanV_epsilony110= fopen("wektorRozwiazanV_epsilony110.txt", "w");
    poissonMetodaAlgebraiczna(100,100,1.0,10.0,0,0,0,0, plikMaczierzA, plikWektorb, plikWektorRozwiazanV_epsilony110);
    return 0;
}

void poissonMetodaAlgebraiczna(int nx, int ny, double epsilon1, double epsilon2, double V1, double V2, double V3, double V4, FILE* plikmatrixA,FILE* plikvectorB,FILE* plikvectorV)
{
    //parametry punkt 2.1 punkt 2.6
    const int N =(nx + 1) * (ny + 1);
    int i = 0;
    int j = 0;
    int nz_num = 0; //ilosc niezerowych elementow
    double sigma = (delta*nx)/10;
    double* epsilonL= new double[N];//(wzor 21)
    double ro = 0.0; //wartosc ro dla innych punktow niz 6
    //punkt 2.2
    double* a = new double[5*N];
    int* ja = new int[5*N];
    int* ia = new int [N+1];
    double* b= new double [N];
    double* V= new double [N];
    //wypelnienie zerami wektorow b i V
    for(int s=0; s<N; s++)
    {
        V[s] = 0.0;
        b[s] = 0.0;
    }
    //wypelnienie -1 wektoru ia
    for(int s = 0; s < N+1; s++)
    {
        ia[s] = -1;
    }
    int k = -1; // numeruje niezerowe elementy A
    //implementacja wzoru 21
    for(int l = 0; l < N; l++)
    {
        j = static_cast<int>(floor(l/(nx+1)));//wzor 12
        i = l-j*(nx+1);//wzor 13
        epsilonL[l] = (i <= nx / 2) ? epsilon1 : epsilon2; //zastapienie if else operatorem (warunek) ? vT : vF;    
    }
    // punkt 3 implementacja algorytmu wypelniania macierzy rzadkiej w formacie CSR + WB Dirichleta
    for(int l = 0; l < N; l++)
    {
        j = static_cast<int>(floor(l/(nx+1)));
        i = l-j*(nx+1);
        int brzeg = 0; //wskaznik położenia: 0 - srodek obszaru; 1 - brzeg
        double vb = 0; //potencjal na brzegu
        if(i == 0)     //lewy brzeg
        {
            brzeg = 1;
            vb = V1;
        }
        if(j == ny)   //gorny brzeg
        {
            brzeg = 1;
            vb = V2;
        }
        if(i == nx)   //prawy brzeg
        {
            brzeg = 1;
            vb = V3;
        }
        if(j == 0)    //dolny brzeg
        {
            brzeg = 1;
            vb = V4;
        }
        if(V1 == 0.0) //obliczenie wartosci ro dla punktu 6 wzory: (4),(5),(25),(26)
        {
          ro = exp(-(pow(delta*i - 0.25*delta*nx, 2)/pow(sigma,2))- (pow(delta*j - 0.5*delta*ny, 2)/pow(sigma,2)))+
          (-exp(-(pow(delta*i - 0.75*delta*nx, 2)/pow(sigma,2))- (pow(delta*j - 0.5*delta*ny, 2)/pow(sigma,2))));
        }
         
        //wypelniamy wektor wyrazow wolnych
        b[l]= -ro; //jesli w srodku jest gestosc

        if(brzeg == 1)
        {
            b[l] =  vb; //wymuszamy wartosc potencjalu na brzegu
        }
        //wypelniamy elementy macierzy A
        ia[l] = -1;//wskanik dla pierwszego elementu w wierszu
        //lewa skrajna przekatna
        if(l - nx - 1 >=0 && brzeg ==0)
        {
            k++;
            if(ia[l]<0)
            { 
                ia[l] =k;
            }    
            a[k]= epsilonL[l]/(delta*delta);// wzor 16
            ja[k] = l - nx - 1;
        }
        //poddiagonala
        if(l-1 >= 0 && brzeg == 0)
        {
            k++;
            if(ia[l]<0)
            {
                ia[l]=k;
            }     
            a[k]= epsilonL[l]/(delta*delta);//wzor 17
            ja[k] = l - 1;
        }
        //diagonala
        k++;
        if(ia[l]<0)
        {
            ia[l]=k;   
        } 
        if(brzeg ==0)
        {
            a[k] = -(2*epsilonL[l] + epsilonL[l+1]+ epsilonL[l+nx+1])/(delta*delta);//wzor 18
        }
        else
        {
            a[k]=1;
        }

        ja[k] = l;
        //naddiagonala
        if(l < N && brzeg == 0)
        {
            k++;
            a[k] = epsilonL[l+1]/(delta*delta);//wzor 19
            ja[k] = l+1;
        }
        // prawa skrajna przekatna
        if(l< N -nx -1 && brzeg ==0)
        {
            k++;
            a[k] = epsilonL[l+nx+1]/(delta*delta);//wzor 20
            ja[k] = l + nx +1;
        }

         //zapis wynikow
        if(control==1)
            fprintf(plikvectorB, "%d\t%d\t%d\t%f\n",l,i,j,b[l]); // zapis wynikow

    }
    for(int l = 0; l < 5*N; l++)
    {
        j = static_cast<int>(floor(l/(nx+1)));
        i = l-j*(nx+1);
        if(control==1)
            fprintf(plikmatrixA, "%d\t%f\n",l,a[l]);
    }

    control++;
    
    nz_num = k+1;
    ia[N] = nz_num;

    pmgmres_ilu_cr(N, nz_num, ia, ja, a, V, b, itr_max, mr, tol_abs, tol_rel);//2.4
    //zapis przerabiamy wynik z wektora 1D na 2D macierz do wykresu matlab
    for(int h = 0; h < N; h++)
    {
        fprintf(plikvectorV, "%g\t", V[h]);
        if((h+1)%(nx+1) == 0)
        {
            fprintf(plikvectorV, "\n");
        }
    }
    delete [] a;
    delete [] ja;
    delete [] ia;
    delete [] b;
    delete [] V;
    delete [] epsilonL;


}