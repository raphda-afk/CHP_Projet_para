#include "operation.h"
#include "fonction.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace std;


//produit matrice vecteur dans le cas explicite
void matrice_vecteur(int space_scheme, std::vector<double>& U, std::vector<double>& result, int Nx, int Ny, double dx, double dy,double xmin, double ymin, double t, int cas, double dt,int iBeg, int iEnd)
{

    for (int k = iBeg ; k <= iEnd ; k++)
    {
        int I = iBeg + k;
        int j = I%Nx;
        int i = I - Nx*j;

        double xi = xmin + i * dx;
        double yj = ymin + j * dy;

        double vx_ij = vx(cas, xi, yj, t);
        double vy_ij = vy(cas, xi, yj, t);

        
        int I_est = (j)*Nx + (i+1) % Nx - (iBeg - Nx);
        int I_ouest = (j)*Nx + (i + (Nx-1))%Nx - (iBeg - Nx);
        int I_nord = ((j+1)%Ny)*Nx + i - (iBeg - Nx);
        int I_sud = ((j + (Ny-1))%Ny)*Nx + i - (iBeg - Nx);

        double dU_dx, dU_dy;



        if ( space_scheme == 1 )
        {
            //cout << "Le schéma choisi est centré" << endl;

            dU_dx = - ( dt / (2.0 * dx) ) * ( U[I_est] - U[I_ouest] );
            dU_dy = - ( dt / (2.0 * dy) ) * ( U[I_nord] - U[I_sud] );

        }
        if ( space_scheme == 2 )
        {
            //cout << "Le schéma choisi est l'upwind" << endl;
            
            if (vx_ij >= 0)
            {
                dU_dx = -( dt / dx ) *( U[I] - U[I_ouest] );
            }
            else
            {
                dU_dx = -( dt / dx ) * ( U[I_est] - U[I] );
            }

            if (vy_ij >= 0)
            {
                dU_dy = -( dt / dy ) * ( U[I] - U[I_sud] );
            }
            else
            {
                dU_dy = -( dt / dy ) * ( U[I_nord] - U[I] );
            }

        }
        result[I] = U[I] + vx_ij * dU_dx + vy_ij * dU_dy;
    }
}


void matrice_vecteur_imp(int space_scheme, std::vector<double>& U, std::vector<double>& result, int Nx, int Ny, double dx, double dy,double xmin, double ymin, double t, int cas, double dt)
{
    for (int i = 1 ; i <= Nx ; i++)
    {
        for (int j = 1 ; j <= Ny ; j++)
        {
            int I = (j-1)*Nx + i-1;

            double xi = xmin + i * dx;
            double yj = ymin + j * dy;

            double vx_ij = vx(cas, xi, yj, t);
            double vy_ij = vy(cas, xi, yj, t);

            int I_est = (j-1)*Nx + i % Nx;
            int I_ouest = (j-1)*Nx + (i-1 + (Nx-1))%Nx;
            int I_nord = (j%Ny)*Nx + i-1;
            int I_sud = ((j-1 + (Ny-1))%Ny)*Nx + i-1;
            

            if ( space_scheme == 1 )
            {
                //cout << "Le schéma choisi est centré" << endl;

                double cx = ( vx_ij * dt ) / ( 2.0 * dx );
                double cy = ( vy_ij * dt ) / ( 2.0 * dy );

                result[I] = U[I] + cx*U[I_est] - cx*U[I_ouest] + cy*U[I_nord] - cy*U[I_sud];

            }
            if ( space_scheme == 2 )
            {
                //cout << "Le schéma choisi est l'upwind" << endl;

                double cx = ( vx_ij * dt ) / ( dx );
                double cy = ( vy_ij * dt ) / ( dy );
                double a = 1.0 + cx + cy;
                
                result[I] = a*U[I] - cx*U[I_ouest] - cy*U[I_sud];

            }
        }
    }
}    


//somme matrice
vector<double> mat_somm(const vector<double>& A, const vector<double>& B, int N)
{
    vector<double> result(N);

    for (int i=0; i < N ; i++)
    {
        result[i] = A[i] + B[i];
    }

    return result;
}

//produit const matrice
vector<double> mat_const(const vector<double>& A, double b, int N)
{
    vector<double> result(N);

    for (int i=0; i < N ; i++)
    {
        result[i] = b*A[i];
    }

    return result;
}

//soustraction de matrice A - B
vector<double> mat_sous(const vector<double>& A, const vector<double>& B, int N)
{
    vector<double> result(N);

    for (int i=0; i < N ; i++)
    {
        result[i] = A[i] - B[i];
    }

    return result;
}

//produit scalaire
double mat_scal(const vector<double>& A, const vector<double>& B, int N)
{
    double result = 0.0;

    for (int i = 0; i < N; i++)
    {
        result = result + A[i]*B[i];
    }
    return result;
}