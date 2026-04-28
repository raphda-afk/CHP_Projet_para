#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <string>
#include <mpi.h>
#include <stdio.h>
#include <chrono>
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



// Rassemble le vecteur distribué U sur le rang 0 et écrit le doc solvoid 
void regroupe_et_ecrit(int ite, const vector<double>& U, int local_size,int Nx, int Ny, double xmin, double ymin, double dx, double dy, int rank, int Np)
{
    vector<double> U_global(Nx * Ny);
    vector<int> recvcounts(Np), displs(Np);
 
    MPI_Gather(&local_size, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
 
    if (rank == 0)
    {
        displs[0] = 0;
        for (int i = 1; i < Np; i++)
        {
            displs[i] = displs[i-1] + recvcounts[i-1];
        }
        
    }
 
    MPI_Gatherv(U.data(), local_size, MPI_DOUBLE, U_global.data(), recvcounts.data(), displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
 
    if (rank == 0)
    {
        ecrire_fichier(ite, Ny, Nx, xmin, ymin, U_global, dx, dy);
    }
}



int main(int argc, char** argv)
{   
    // pour compil:
    // g++ main.cpp fonction.cpp param.cpp operation.cpp -o main
    // convert -delay 10 -loop 0 sol.*.png animation.gif
    // plot 'sol.0.dat' u 2:3, 'sol.1.dat' u 2:3, exp(-((x-0.2-0.5*0.2)/0.05)**2)
    
    //----------------------------ouverture de l'espace //----------------------------

    int tag(200);
    MPI_Status status;
    int rank;
    int Np;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &Np);
    MPI_Barrier(MPI_COMM_WORLD);

    //----------------------------Lecture du fichier----------------------------
    // ifstream fichier("param.dat");

    int cas;
    double xmin, xmax, ymin, ymax;
    double Tf;
    int Nx, Ny;
    double CFL;
    int space_scheme, time_scheme;

    // Lire les param dans le fichier des param
    lireParametres(cas, xmin, xmax, ymin, ymax, Tf, Nx, Ny, CFL, space_scheme, time_scheme);

    if ( rank == 0)
    {
        cout << "---------------------------------------" << endl;
        cout << "----------Lecture du fichier-----------"<<endl;
        printf("\nLe cas choisi est le n° %d\n", cas);
        cout << "Le domaine est [" << xmin << " ; " << xmax << "]x[" << ymin << " ; " << ymax << "].\n" << endl;
        cout << "Le temps final est " << Tf << "\n" << endl;
        cout << "La CFL est " << CFL << "\n" << endl;
        if ( space_scheme == 1 )
        {
            cout << "Le schéma choisi est : centré" << endl;
        }
        if ( space_scheme == 2 )
        {
            cout << "Le schéma choisi est : upwind" << endl;
        }
        if ( time_scheme == 1 )
        {
            cout << "Le schéma choisi est : explicite\n" << endl;
        }
        if ( time_scheme == 2 )
        {
            cout << "Le schéma choisi est : implicite\n" << endl;
        }
    }
    //a améliorer pour le tps d'exec (faire lire 1 seul + Bcast pour distrib)

    //----------------------------def des vecteur et constante restante----------------------------
    
    // définition des param utile
    int N = Nx * Ny;
    double dx = (xmax - xmin) / Nx;
    double dy = (ymax - ymin) / Ny;
    
    // def et calcul de la CFL
    double dt;
    dt = cfl(cas, dx, dy, CFL); 
    int Nt = int(Tf / dt);
    
    // repartition des charges
    int iBeg;
    int iEnd;
    charge_a(rank, N ,Np, &iBeg,&iEnd);
    
    // vector<double> U_rank(iEnd - iBeg + 1); //partie de matrice pour la para
    // vecteur pour la résolution 
    vector<double> U_final(N);                      //vecteur pour le résultat final
    vector<double> U_loc( iEnd-iBeg + 1 + 2*Nx );   // vecteur pour la // 
    vector<double> U(iEnd-iBeg + 1 );               // vecteur local du proc 
    vector<double> LU(iEnd-iBeg + 1 );              // vecteur local du proc utilisé pour le t^n+1
    vector<double> U0(iEnd-iBeg + 1 );              //vecteur local du proc pour la condition init

    // calcul des proc précédent et suivant pour les MPI_SEND et MPI_RECV
    int proc_prec = (rank - 1 + Np)%Np;     
    int proc_sui = (rank + 1)%Np;
    
    //----------------------------initialisation de U----------------------------
    if ( rank == 0)
    {
        cout << "---------------------------------------" << endl;
        cout << "----------initialisation de U----------"<<endl;
    }

    auto start = std::chrono::high_resolution_clock::now();

    for (int k = iBeg ; k <= iEnd ; k++)
    {
        int I = k;
        int i = I%Nx + 1;
        int j = I/Nx + 1;

        double xi = xmin + i * dx;
        double yj = ymin + j * dy;
        U[k-iBeg] = v_init(cas, xi, yj);
    }

    if ( rank == 0)
    {
        cout << "\n l'initialisation est effectuée \n" <<endl;
    }


    // Écriture de la solution initiale (t = 0)
    regroupe_et_ecrit(0, U, iEnd - iBeg + 1 , Nx, Ny, xmin, ymin, dx, dy, rank, Np);
    
    if ( rank == 0)
    {
        cout << "\n la solution pour t=0 est enregistré : sol.O.dat \n" <<endl;
    }

    //----------------------------boucle en temps----------------------------
    if ( rank == 0)
    {
        cout << "---------------------------------------" << endl;
        cout << "------------boucle en temps------------"<<endl;
    }

    //ecrire_fichier(0, Ny , Nx , xmin, ymin, U, dx, dy);

    //choix du nombre de solution enregistré
    int ns = 49;
    int nb_sol = int( Nt / ns);
    int compteur = 1;

    if ( rank == 0)
    {
        cout << "\n " << ns << " solutions seront enregistrées" <<endl;
    }



    double t = 0.0;
    for (int n = 0; n < Nt ; n++)
    {   
        // remplissage de la matrice U_loc
        //envoie et reception (cas paire et impaire diff pour eviter blockage)
        if (rank % 2 == 0)
        {
            MPI_Send(&U[iEnd-iBeg-Nx+1], Nx, MPI_DOUBLE, proc_sui,  tag, MPI_COMM_WORLD);
            MPI_Recv(&U_loc[0], Nx, MPI_DOUBLE, proc_prec, tag, MPI_COMM_WORLD, &status);
        }
        else
        {
            MPI_Recv(&U_loc[0], Nx, MPI_DOUBLE, proc_prec, tag, MPI_COMM_WORLD, &status);
            MPI_Send(&U[iEnd-iBeg-Nx+1], Nx, MPI_DOUBLE, proc_sui,  tag, MPI_COMM_WORLD);
        }

        if (rank % 2 == 0)
        {
            MPI_Send(&U[0], Nx, MPI_DOUBLE, proc_prec, tag, MPI_COMM_WORLD);
            MPI_Recv(&U_loc[iEnd-iBeg+1+Nx], Nx, MPI_DOUBLE, proc_sui,  tag, MPI_COMM_WORLD, &status);
        }
        else
        {
            MPI_Recv(&U_loc[iEnd-iBeg+1+Nx], Nx, MPI_DOUBLE, proc_sui,  tag, MPI_COMM_WORLD, &status);
            MPI_Send(&U[0], Nx, MPI_DOUBLE, proc_prec, tag, MPI_COMM_WORLD);
        }

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
            BiCG( space_scheme, U_loc, LU, Nx, Ny, dx, dy, xmin, ymin, t, cas, dt, Nt, iBeg, iEnd);
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

    // Écriture de la solution finale (t = Tf)
    regroupe_et_ecrit(1, U, iEnd - iBeg + 1 , Nx, Ny, xmin, ymin, dx, dy, rank, Np);
    
    if ( rank == 0)
    {
        cout << "\n la solution pour t=Tf est enregistré : sol.1.dat \n" <<endl;
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> float_ms = end - start;

    //MPI_WTIME
    if (rank == 0)
    {
        printf("==============================\n");
        cout << "Np = " << Np << endl;
        cout << "Temps d'execution: " << float_ms.count() << " milliseconds" << endl;
        printf("==============================\n");
    }

    MPI_Finalize();
    return 0;
}