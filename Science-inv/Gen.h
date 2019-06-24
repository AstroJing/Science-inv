#ifndef HEAD
#define HEAD
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<direct.h>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include<ctime>
#define pi 3.1415926535
#define hbar    197.3286
std::pair<double, double> GetTrho(double L, double Ksym, double Jsym, double J0, double K0);
double solver(double U, double L, double Ksym, double Jsym, double J0, double K0, int c, double* point);
double PSolver(double rho, double delta, double Esym, double dE0, double dEsym,double* erhol);
#endif
