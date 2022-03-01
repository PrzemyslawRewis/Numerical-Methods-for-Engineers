// IMN projekt 9
// Przemyslaw_Rewis 29.12.2021
// Wymagany minimalny standard jezyka -std=c++11
// Dyfuzja ciepla - metoda Cranck-Nicloson.

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>


//parametry punkt 3.1

constexpr short int nx=40;
constexpr short int ny=40;
constexpr short int N=(nx+1)*(ny+1);
constexpr double delta=1.0;
constexpr double delta_t=1.0;
constexpr short int TA=40;
constexpr short int TB=0;
constexpr short int TC=30;
constexpr short int TD=0;
constexpr double kB=0.1;
constexpr double kD=0.6;
constexpr short int IT_MAX=2000;


void solver();

int main()
{
    solver();
    return 0;
}

void solver()
{
    int l = 0;
    int signum = 0;
    //tworzenie macierzy A, B i wektora c,d oraz T 3.2 3.3
    gsl_matrix *A = gsl_matrix_calloc(N, N);
    gsl_matrix *B = gsl_matrix_calloc(N, N);
    gsl_vector *c = gsl_vector_calloc(N);
    gsl_vector *d = gsl_vector_calloc(N);
    gsl_vector *T = gsl_vector_calloc(N);
    

    gsl_permutation *p = gsl_permutation_alloc(N);
    

    // wypelnienie macierzy A i B 1.1.1
    //szary obszar na rysunku
    for(int i = 1; i <= nx - 1; i++)
    {
        for(int j = 1; j <=ny - 1; j++)
        {
            l = i + j*(nx+1);//wzor 9
            //wzor 17
            gsl_matrix_set(A, l, l-nx-1, delta_t / (2*delta*delta));
            gsl_matrix_set(A, l, l-1,    delta_t / (2*delta*delta));
            gsl_matrix_set(A, l, l+1,    delta_t / (2*delta*delta));
            gsl_matrix_set(A, l, l+nx+1, delta_t / (2*delta*delta));
            //wzor 18
            gsl_matrix_set(A, l, l,      -(2*delta_t)/(delta*delta) - 1);
            //wzor 19
            gsl_matrix_set(B, l, l-nx-1, -delta_t / (2*delta*delta));
            gsl_matrix_set(B, l, l-1,    -delta_t / (2*delta*delta));
            gsl_matrix_set(B, l, l+1,    -delta_t / (2*delta*delta));
            gsl_matrix_set(B, l, l+nx+1, -delta_t / (2*delta*delta));
            //wzor 20
            gsl_matrix_set(B, l, l,      (2*delta_t)/(delta*delta) - 1);      
        }
    }
    //lewy i prawy brzeg
    for(int i=0; i<=nx; i+=nx)
    {
        for(int j=0; j<=ny; j++)
        {
            l = i + j*(nx+1);//wzor 9
            gsl_matrix_set(A, l, l, 1);//wzor 23
            gsl_matrix_set(B, l, l, 1);//wzor24
            gsl_vector_set(c, l, 0);//wzor 25
        }
    }
    // WB von Neumanna na gornym brzegu dla chwili n + 1
    for(int i=1; i<=nx-1; i++)
    {
        l = i + ny*(nx+1);//wzor 9, j=ny wzor 31
        gsl_matrix_set(A, l, l-nx-1, - 1.0/(kB*delta)); //wzor 32
        gsl_matrix_set(A, l, l, 1 + 1.0/(kB*delta));//wzor 33 
        gsl_vector_set(c, l, TB);//wzor 34
        for( int k=0; k<N; k++ )
        {
            gsl_matrix_set(B, l, k, 0);//wzor 35
        }
    }
    // WB von Neumanna na dolnym brzegu dla chwili n + 1
    for(int i=1; i<=nx-1; i++)
    {
        l = i + 0*(nx+1);//wzor 9, j=0 wzor 41
        gsl_matrix_set(A, l, l, 1 + 1.0/(kD*delta)); //wzor 42
        gsl_matrix_set(A, l, l+nx+1, - 1.0/(kD*delta)); //wzor 43
        gsl_vector_set(c, l, TD);//wzor 44
        for( int k=0; k<N; k++ )
        {
            gsl_matrix_set(B, l, k, 0);//wzor 45
        }

    }


    //punkt 3.3 narzucenie warunkow poczatkowych
    for(int j=0; j<=ny; j++)
    {
        l = 0 + j*(nx+1);
        gsl_vector_set(T, l, TA);//wzor 46
    }
    for(int j=0; j<=ny; j++)
    {
        l = nx + j*(nx+1);
        gsl_vector_set(T, l, TC);//wzor 47
    }
    for(int i=1; i<=nx-1; i++)
    {
        for(int j=0; j<=ny; j++)
        {
            l = i + j*(nx+1);
            gsl_vector_set(T, l, 0);//wzor 48
        }
    }


    FILE* plik1,*plik2;
    plik1 = fopen("T.txt", "w");
    plik2 = fopen("d2T.txt", "w");


    gsl_linalg_LU_decomp(A, p, &signum);//punkt 3.4
    
    // punkt3.5
    for(int iteracja=0; iteracja<=IT_MAX; iteracja++)//punkt 3.6
    {
        gsl_blas_dgemv( CblasNoTrans, 1, B, T, 0, d );//d=B*T
        gsl_blas_daxpy( 1, c, d );// d = c + d 
        gsl_linalg_LU_solve( A, p, d, T );//rozwiazanie ukladu rownan uzywajac LU
        if( iteracja == 100 || iteracja == 200 || iteracja == 500 || iteracja == 1000 || iteracja == 2000 )//3.7
        { 
            for( int i=0; i<N; i++ )
            {
                fprintf(plik1, "%g\t", gsl_vector_get(T, i));
                if((i+1)%(nx+1)==0)
                {
                    fprintf(plik1, "\n");
                }
            }
            for(int i=1; i<=nx-1; i++)
            {
                for(int j=1; j<=ny-1; j++)
                {
                    l = i + j*(nx+1);
                    fprintf(plik2, "%g\t",(( gsl_vector_get(T, l+1) - 2 * gsl_vector_get(T, l) + gsl_vector_get(T, l-1))/(delta*delta))
                                      + (( gsl_vector_get(T, l+nx+1) - 2*gsl_vector_get(T, l) + gsl_vector_get(T, l-nx-1))/(delta*delta))); 
                }
                fprintf(plik2, "\n");
            }
            fprintf(plik2, "\n");  
        }
    }
    gsl_matrix_free(A);
    gsl_matrix_free(B);
    gsl_vector_free(c);
    gsl_vector_free(T);
    gsl_permutation_free(p);   
}

