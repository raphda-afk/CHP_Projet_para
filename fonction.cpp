#include "fonction.h"
#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

double vx(int cas, double x, double y, double t)
{
    //cout << "vx est lancé" << endl;

    // Validation
    if (cas == 0)
    {
        return 0;
    }
    // Translation d’une gaussienne
    if (cas == 1)
    {
        return 1.;
    }
    // Translation d’un cylindre
    if (cas == 2)
    {
        return 0.5;
    }
    // Rotation d’un cylindre
    if (cas == 3)
    {
        return -y;
    }
    else
    {
        return 0;
    }
}


double vy(int cas, double x, double y, double t)
{   
    //cout << "vy est lancé" << endl;

    // Validation
    if (cas == 0)
    {
        return 0.5;
    }
    // Translation d’une gaussienne
    if (cas == 1)
    {
        return 0;
    }
    // Translation d’un cylindre
    if (cas == 2)
    {
        return 0.5;
    }
    // Rotation d’un cylindre
    if (cas == 3)
    {
        return x;
    }
    else
    {
        return 0;
    }
}

double v_init(int cas, double x, double y)
{   
    //cout << "v_init est lancé" << endl;

    // Validation
    if (cas == 0)
    {
        double x0 = 0.2;
        return exp( (-(y - x0)*(y - x0)/((0.05)*0.05) ));
    }
    // Translation d’une gaussienne
    if (cas == 1)
    {
        return exp( (-(x*x)/(0.0075) ) + (-(y*y)/0.0075) );
    }
    // Translation d’un cylindre
    if (cas == 2)
    {
        if ( (sqrt( x*x + y*y)) <= 0.4  )
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
    // Rotation d’un cylindre
    if (cas == 3)
    {
        double x0 = 0.5;
        double y0 = 0.3;

        if (sqrt((x - x0)*(x - x0) + (y - y0)*(y - y0)) <= 0.4)
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
    else
    {
        return 0;
    }
}

double v_exacte(int cas, double x, double y, double t)
{
    //cout << "v_exacte est lancé" << endl;
    if (cas == 0)
    {
        double x0 = 0.2;
        double y_exact = y - 0.5 * t;
        return exp(-(y_exact - x0)*(y_exact - x0) / (0.0025));
    }
    return 0;
}

double cfl(int cas,double dx,double dy, double CFL)
{
    double vxmax;
    double vymax;

    // Validation
    if (cas == 0)
    {
        vxmax = 0;
        vymax = 0.5;
    }
    // Translation d’une gaussienne
    else if  (cas == 1)
    {
        vxmax = 1;
        vymax = 0;
    }
    // Translation d’un cylindre
    else if  (cas == 2)
    {
        vxmax = 0.5;
        vymax = 0.5;
    }
    // Rotation d’un cylindre
    else if  (cas == 3)
    {
        vxmax = 1;
        vymax = 1;
    }
    else 
    {
        cout << "\n!!!!!!!!!!\nCAS NON VALIDE\n!!!!!!!!!!\n" << endl;
        vxmax = 1;
        vymax = 1;
    }
    
    return CFL *  1 / ( (vxmax / dx)  + (vymax / dy) );
}