//utilities.cpp
#include "utilities.h"
void chk(std::ifstream &fin, std::string str)
{
    if (fin.fail()) {std::cout << str+" file not found" << std::endl; std::exit(0);}
}
void skipLine(std::ifstream &fin, int n)
{
    for (int i=0; i<n; i++)
    {
        fin.ignore(256,'\n');
    }
}
double sq(double x)
{
    return x*x;
}
int sq(int n)
{
    return n*n;
}

