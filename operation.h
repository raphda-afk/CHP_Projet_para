#ifndef _OPERATION_H
#define _OPERATION_H

#include <fstream>
#include <vector>

using namespace std;

void matrice_vecteur(int space_scheme, std::vector<double>& U, std::vector<double>& result, int Nx, int Ny, double dx, double dy,double xmin, double ymin, double t, int cas, double dt,int iBeg, int iEnd);
void matrice_vecteur_imp(int space_scheme, std::vector<double>& U, std::vector<double>& result, int Nx, int Ny, double dx, double dy,double xmin, double ymin, double t, int cas, double dt,int iBeg, int iEnd);
vector<double> mat_somm(const vector<double>& A, const vector<double>& B, int N);
vector<double> mat_const(const vector<double>& A, double b, int N);
vector<double> mat_sous(const vector<double>& A, const vector<double>& B, int N);
double mat_scal(const vector<double>& A, const vector<double>& B, int N);

#endif
