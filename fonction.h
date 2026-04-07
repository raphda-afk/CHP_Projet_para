#ifndef _FONCTION_H
#define _FONCTION_H

#include <fstream>
#include <vector>

double vx(int cas, double x, double y, double t);
double vy(int cas, double x, double y, double t);
double v_init(int cas, double x, double y);
double v_exacte(int cas,double x, double y, double t);
double cfl(int cas,double dx,double dy, double CFL);



#endif