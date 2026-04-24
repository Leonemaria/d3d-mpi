#ifndef UTILITIES_H
#define UTILITIES_H
//utilities.h
#include <iostream>
#include <fstream>
#include <cmath>
#include "matrix.h"
#include "vectors.h"
//
// this is a collection of usefull functions
//typedef double (*iniPar)(vector2D, int);
void chk(std::ifstream &fin, std::string str);
void skipLine(std::ifstream &fin, int n);
double sq(double x);
int sq(int n);
#endif

