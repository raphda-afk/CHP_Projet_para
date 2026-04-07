#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <string>
#include <mpi.h>
#include <stdio.h>
#include "fonction.h"
#include "param.h"
#include "operation.h"
#include "BICGstab.h"


using namespace std;


void ecrire_fichier(int ite, double Ny , double Nx , double xmin, double ymin, const std::vector<double>& U, double dx, double dy)
{
    // Construire le nom du fichier
    std::string nom_fichier = "sol." + std::to_string(ite) + ".dat";

    std::ofstream fichier(nom_fichier);

    if (!fichier) {
        std::cerr << "Erreur ouverture fichier !" << std::endl;
        return;
    }

    for (int j = 1; j <= Ny; j++) 
    {
        for (int i = 1; i <= Nx; i++) 
        {
            double xi = xmin + i * dx;
            double yj = ymin + j * dy;
            fichier << xi << " " << yj << " " <<  U[(j-1)*Nx + (i-1)] << endl;
        }
    }

    fichier.close();
}



void charge_a(int me, int n,int np, int *iBeg, int *iEnd)
{
    int k = n / np;
    int reste = n % np;

    if (me < reste)
    {
        *iBeg = me * (k+1);
        *iEnd = (me+1) * (k+1)- 1;
    }
    else
    {
        *iBeg = me * k + reste ;
        *iEnd = (me+1) * k- 1 + reste;
    }
}

int main(int argc, char** argv)
{   
    // pour compil:
    // g++ main.cpp fonction.cpp param.cpp operation.cpp -o main
    // convert -delay 10 -loop 0 sol.*.png animation.gif
    
    //----------------------------Lecture du fichier----------------------------
    cout << "------------------------------------" << endl;
    cout << "----------Lecture du fichier----------"<<endl;
    // ifstream fichier("param.dat");

    int cas;
    double xmin, xmax, ymin, ymax;
    double Tf;
    int Nx, Ny;
    double CFL;
    int space_scheme, time_scheme;

    // Lire les param dans le fichier des param
    lireParametres(cas, xmin, xmax, ymin, ymax, Tf, Nx, Ny, CFL, space_scheme, time_scheme);

    //----------------------------boucle de temps----------------------------
 
    // définition des param restant
    int N = Nx * Ny;
    double dx = (xmax - xmin) / Nx;
    double dy = (ymax - ymin) / Ny;

    double dt;
    dt = cfl(cas, dx, dy, CFL); 
    int Nt = int(Tf / dt);

    vector<double> U_final(N);



    int tag(200);
    MPI_Status status;
    int rank;
    int Np;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &Np);
    MPI_Barrier(MPI_COMM_WORLD);

    int iBeg;
    int iEnd;
    charge_a(rank, N ,Np, &iBeg,&iEnd);

    // vector<double> U_rank(iEnd - iBeg + 1); //partie de matrice pour la para
    vector<double> U_loc( iEnd-iBeg + 1 + 2*Nx );
    int proc_prec = (rank - 1 + Np)%Np; 
    int proc_sui = (rank + 1)%Np;



    vector<double> U(iEnd-iBeg + 1 );       
    vector<double> LU(iEnd-iBeg + 1 );
    vector<double> U0(iEnd-iBeg + 1 );

    // Initialisation
    for (int k = iBeg ; k <= iEnd ; k++)
    {
        int I = iBeg + k;
        int j = I%Nx;
        int i = I - Nx*j;

        double xi = xmin + i * dx;
        double yj = ymin + j * dy;
        U[k-iBeg] = v_init(cas, xi, yj);
    }

    //ecrire_fichier(0, Ny , Nx , xmin, ymin, U, dx, dy);
    int nb_sol = int( Nt / 49);
    int compteur = 1;






    double t = 0.0;
    for (int n = 0; n < Nt ; n++)
    {   

        MPI_Send(&U[0], Nx, MPI_DOUBLE, proc_prec , tag, MPI_COMM_WORLD);
        MPI_Send(&U[iEnd-iBeg - Nx], Nx, MPI_DOUBLE, proc_sui , tag, MPI_COMM_WORLD);
        
        //construction de U_loc
        MPI_Recv(&U_loc[0], Nx, MPI_DOUBLE, proc_prec , tag, MPI_COMM_WORLD, &status);
        MPI_Recv(&U_loc[iEnd-iBeg + 1 + Nx ], Nx , MPI_DOUBLE, proc_sui , tag, MPI_COMM_WORLD, &status);
        for (int k = Nx; k <= iEnd-iBeg + Nx; k++ )
        {
            U_loc[k] = U[k-Nx];
        }

        U0 = U;

        // cas Explicite
        if ( time_scheme == 1 )
        {
            //produit matrice vecteur U^n+1=AU^n 
            matrice_vecteur(space_scheme, U_loc, LU, Nx, Ny, dx, dy, xmin, ymin, t, cas, dt, iBeg, iEnd);
        }

        //cas Implicite
        if ( time_scheme == 2)
        {
            //produit matrice vecteur U^n+1=AU^n 
            BiCG( space_scheme, U, LU, Nx, Ny, dx, dy, xmin, ymin, t, cas, dt, Nt);
        }
        U = LU;

        t = t + dt;
       
        // //ecriture des fichier 
        // if ( n%nb_sol == 0)
        // {
        //     ecrire_fichier(compteur, Ny , Nx , xmin, ymin, U, dx, dy);
        //     compteur = compteur + 1;
        // }
    }
    //ecrire_fichier( compteur+1 , Ny , Nx , xmin, ymin, U, dx, dy);

    // re écrire le vecteur U en partageant les valeurs des différentes parties

    // Taille locale
    int local_size = iEnd - iBeg + 1;

    // Tableaux pour Gatherv
    vector<int> recvcounts(Np);
    vector<int> displs(Np);

    // Chaque proc envoie sa taille au proc 0
    MPI_Gather(&local_size, 1, MPI_INT, recvcounts.data(), 1, MPI_INT,0, MPI_COMM_WORLD);

    // Construire les déplacements (uniquement sur proc 0)
    if (rank == 0)
    {
        displs[0] = 0;
        for (int i = 1; i < Np; i++)
        {
            displs[i] = displs[i-1] + recvcounts[i-1];
        }
    }

    // Gather des données
    MPI_Gatherv(U.data(), local_size, MPI_DOUBLE,U_final.data(), recvcounts.data(), displs.data(), MPI_DOUBLE,0, MPI_COMM_WORLD);



    MPI_Finalize();
    return 0;
}