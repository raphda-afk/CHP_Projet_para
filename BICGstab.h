#ifndef _BICGstab_H
#define _BICGstab_H

#include <fstream>
#include <vector>

using namespace std;

void BiCG(int space_scheme, std::vector<double>& U0, std::vector<double>& U1, int Nx, int Ny, double dx, double dy,double xmin, double ymin, double t, int cas, double dt, int Nt);


#endif