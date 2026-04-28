#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <string>
#include <cmath>
#include <mpi.h>
#include "fonction.h"
#include "param.h"
#include "operation.h"
#include "BICGstab.h"


void BiCG(int space_scheme, std::vector<double>& U0, std::vector<double>& U1, int Nx, int Ny, double dx, double dy,double xmin, double ymin, double t, int cas, double dt, int Nt,int iBeg,int iEnd) 
{
    int N = iEnd - iBeg + 1 ;
    std::vector<double> x = U0;
    std::vector<double> r(N), p(N), r_bis(N), Ax(N), Atx(N), mu(N), h(N), s(N), t2(N), U(N);
    std::vector<double> p_plus( N + 2*Nx ),s_plus( N + 2*Nx );
    double rho, rho_n, norme_s, norme_r, w;
    double epsilon = 10.0e-15;

    int tag(100);
    MPI_Status status;
    int rank;
    int Np;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &Np);

    // calcul des proc précédent et suivant pour les MPI_SEND et MPI_RECV
    int proc_prec = (rank - 1 + Np)%Np;     
    int proc_sui = (rank + 1)%Np;

    norme_s = 10.0;
    norme_r = 10.0;

    matrice_vecteur_imp( space_scheme, x, Ax,  Nx,  Ny,  dx,  dy, xmin,  ymin,  t,  cas,  dt, iBeg, iEnd);
    for (int k = 0 ; k < iEnd - iBeg + 1; k++)
    {
        U[k] = U0[k + Nx];
        
    }

    std::cout << "t=" << t << std::endl;
    if ( t == 0.1)
    {
        std::cout << "==============================================\n" << std::endl;
        for (int k = 0 ; k < iEnd - iBeg + 1; k++)
        {
            std::cout << "U[" << iBeg + k << "]= " << U[k] << std::endl;
        }

        std::cout << "\n==============================================\n" << std::endl;
    }

    r = mat_sous(U, Ax, N);
    r_bis = r;
    rho_n = mat_scal(r_bis, r, N);
    MPI_Allreduce(&rho_n, &rho ,1,  MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);
    // std::cout << "U[0]= " << rho << std::endl;

    p = r;
    //transmettre les val de p qui sont dans les Nx avant et Nx apres comme fait dans le main 
    //!! norme s et norme r ? 

    // Boucle while
    while ( (norme_s > epsilon) and (norme_r > epsilon) )
    {
        double alpha, beta, rho_next ;  
        double alpha_n, w_n_haut, w_n_bas, w_n1, w_n2;  
        double rho_n, norme_s_n, norme_r_n ;

        // remplissage de la matrice p_plus
        //envoie et reception (cas paire et impaire diff pour eviter blockage)
        if (rank % 2 == 0)
        {
            MPI_Send(&p[iEnd-iBeg-Nx+1], Nx, MPI_DOUBLE, proc_sui,  tag, MPI_COMM_WORLD);
            MPI_Recv(&p_plus[0], Nx, MPI_DOUBLE, proc_prec, tag, MPI_COMM_WORLD, &status);
        }
        else
        {
            MPI_Recv(&p[0], Nx, MPI_DOUBLE, proc_prec, tag, MPI_COMM_WORLD, &status);
            MPI_Send(&p_plus[iEnd-iBeg-Nx+1], Nx, MPI_DOUBLE, proc_sui,  tag, MPI_COMM_WORLD);
        }

        if (rank % 2 == 0)
        {
            MPI_Send(&p[0], Nx, MPI_DOUBLE, proc_prec, tag, MPI_COMM_WORLD);
            MPI_Recv(&p_plus[iEnd-iBeg+1+Nx], Nx, MPI_DOUBLE, proc_sui,  tag, MPI_COMM_WORLD, &status);
        }
        else
        {
            MPI_Recv(&p_plus[iEnd-iBeg+1+Nx], Nx, MPI_DOUBLE, proc_sui,  tag, MPI_COMM_WORLD, &status);
            MPI_Send(&p[0], Nx, MPI_DOUBLE, proc_prec, tag, MPI_COMM_WORLD);
        }

        for (int k = Nx; k <= iEnd-iBeg + Nx; k++ )
        {
            p_plus[k] = p[k-Nx];
        }


        matrice_vecteur_imp( space_scheme, p_plus , mu, Nx, Ny, dx, dy, xmin, ymin, t, cas, dt, iBeg, iEnd);
        
        alpha_n = mat_scal(r_bis, mu, N);
        //pour que a = somme de alpha_proc
        MPI_Allreduce(&alpha_n, &alpha ,1,  MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);
        alpha = rho / alpha;
        // std::cout << "U[0]= " << alpha << std::endl;


        h = mat_somm(x, mat_const(p, alpha, N), N);
        s = mat_sous(r, mat_const(mu, alpha, N), N);

        // remplissage de la matrice p_plus
        //envoie et reception (cas paire et impaire diff pour eviter blockage)
        if (rank % 2 == 0)
        {
            MPI_Send(&s[iEnd-iBeg-Nx+1], Nx, MPI_DOUBLE, proc_sui,  tag, MPI_COMM_WORLD);
            MPI_Recv(&s_plus[0], Nx, MPI_DOUBLE, proc_prec, tag, MPI_COMM_WORLD, &status);
        }
        else
        {
            MPI_Recv(&s[0], Nx, MPI_DOUBLE, proc_prec, tag, MPI_COMM_WORLD, &status);
            MPI_Send(&s_plus[iEnd-iBeg-Nx+1], Nx, MPI_DOUBLE, proc_sui,  tag, MPI_COMM_WORLD);
        }

        if (rank % 2 == 0)
        {
            MPI_Send(&s[0], Nx, MPI_DOUBLE, proc_prec, tag, MPI_COMM_WORLD);
            MPI_Recv(&s_plus[iEnd-iBeg+1+Nx], Nx, MPI_DOUBLE, proc_sui,  tag, MPI_COMM_WORLD, &status);
        }
        else
        {
            MPI_Recv(&s_plus[iEnd-iBeg+1+Nx], Nx, MPI_DOUBLE, proc_sui,  tag, MPI_COMM_WORLD, &status);
            MPI_Send(&s[0], Nx, MPI_DOUBLE, proc_prec, tag, MPI_COMM_WORLD);
        }

        for (int k = Nx; k <= iEnd-iBeg + Nx; k++ )
        {
            s_plus[k] = s[k-Nx];
        }

        matrice_vecteur_imp( space_scheme,  s_plus,  t2,  Nx,  Ny, dx, dy, xmin, ymin, t, cas, dt, iBeg, iEnd);

        w_n_bas =mat_scal(t2, t2, N);
        MPI_Allreduce(&w_n_bas, &w_n1 ,1,  MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);
        w_n_haut = mat_scal(t2, s, N);
        MPI_Allreduce(&w_n_haut, &w_n2 ,1,  MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);
        w = w_n2 / w_n1;
        

        x = mat_somm(h , mat_const(s, w, N), N);

        r = mat_sous(s, mat_const(t2, w, N), N);

        rho_n = mat_scal(r_bis, r, N);
        MPI_Allreduce(&rho_n, &rho_next ,1,  MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);

        beta = (rho_next / rho) * (alpha / w);

        p = mat_somm(r, mat_const(mat_sous(p, mat_const(mu, w, N), N),beta , N), N);



        norme_s_n = mat_scal(s, s, N);
        MPI_Allreduce(&norme_s_n, &norme_s ,1,  MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);
        norme_s = sqrt(norme_s);


        if ( norme_s <= epsilon)
        {
            x = h;
            break;
        }

        norme_r_n = mat_scal(r, r, N);
        MPI_Allreduce(&norme_r_n, &norme_r ,1,  MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);
        norme_r = sqrt(norme_r);

        rho = rho_next;
        //cout << "alpha= " << alpha << endl;
    } 
    U1 = x;
}