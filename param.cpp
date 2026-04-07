#include "param.h"
#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

void lireParametres(int& cas, double& xmin, double& xmax, double& ymin, double& ymax, double& Tf, int& Nx, int& Ny, double& CFL, int& space_scheme, int& time_scheme)
{   

    ifstream fichier("param.dat");

    if (!fichier)
    {
        cout << "Erreur : document non lisible."<<endl;
    }

    //int cas;
    fichier >> cas;
    printf("Le cas choisi est le n° %d\n", cas);

    //double xmin , xmax , ymin , ymax;
    fichier >> xmin >> xmax >> ymin >> ymax;
    cout << "Le domaine est [" << xmin << " ; " << xmax << "]x[" << ymin << " ; " << ymax << "].\n" << endl;

    //double Tf;
    fichier >> Tf;
    cout << "Le temps final est " << Tf << "\n" << endl;

    //int Nx , Ny;
    fichier >> Nx >> Ny;

    //double CFL;
    fichier >> CFL;
    cout << "La CFL est " << CFL << "\n" << endl;

    //int space_scheme , time_scheme;
    fichier >> space_scheme >> time_scheme;
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
        cout << "Le schéma choisi est : explicite" << endl;
    }
    if ( time_scheme == 2 )
    {
        cout << "Le schéma choisi est : implicite" << endl;
    }
}