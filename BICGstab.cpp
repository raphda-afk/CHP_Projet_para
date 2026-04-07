#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <string>
#include <cmath>
#include "fonction.h"
#include "param.h"
#include "operation.h"
#include "BICGstab.h"


void BiCG(int space_scheme, std::vector<double>& U0, std::vector<double>& U1, int Nx, int Ny, double dx, double dy,double xmin, double ymin, double t, int cas, double dt, int Nt) 
{
    int N = Nx * Ny;
    std::vector<double> x = U0;
    std::vector<double> r(N), p(N), r_bis(N), Ax(N), Atx(N), mu(N), h(N), s(N), t2(N);
    double rho, norme_s, norme_r, w;
    double epsilon = 10.0e-15;

    norme_s = 10.0;
    norme_r = 10.0;
    matrice_vecteur_imp( space_scheme, x, Ax,  Nx,  Ny,  dx,  dy, xmin,  ymin,  t,  cas,  dt);
    // Initialisation
    
    r = mat_sous(U0, Ax, N);
    r_bis = r;
    rho = mat_scal(r_bis, r, N);
    p = r;

    // Boucle while
    while ( (norme_s > epsilon) and (norme_r > epsilon) )
    {
        double alpha, beta, rho_next ;                  

        matrice_vecteur_imp( space_scheme, p, mu, Nx, Ny, dx, dy, xmin, ymin, t, cas, dt);
        
        alpha = rho / mat_scal(r_bis, mu, N);
        // pour scal il faut faire reduction (r_bis pas complet)
        h = mat_somm(x, mat_const(p, alpha, N), N);
        s = mat_sous(r, mat_const(mu, alpha, N), N);

        matrice_vecteur_imp( space_scheme,  s,  t2,  Nx,  Ny, dx, dy, xmin, ymin, t, cas, dt);

        w = mat_scal(t2, s, N) / mat_scal(t2, t2, N);
        // pour scal il faut faire reduction (pas complet)

        x = mat_somm(h , mat_const(s, w, N), N);

        r = mat_sous(s, mat_const(t2, w, N), N);

        rho_next = mat_scal(r_bis, r, N);
        // pour scal il faut faire reduction (pas complet)

        beta = (rho_next / rho) * (alpha / w);

        p = mat_somm(r, mat_const(mat_sous(p, mat_const(mu, w, N), N),beta , N), N);

        norme_s = sqrt(mat_scal(s, s, N));
        // j'ai pas compris mais attention

        if ( norme_s <= epsilon)
        {
            x = h;
            break;
        }
        norme_r = sqrt(mat_scal(r, r, N));
        // pareil
        rho = rho_next;
        
    }

    U1 = x;
}