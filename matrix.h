#ifndef MATRIX_H
#define MATRIX_H
// matrix.h
#include<cmath>
#include <cstdlib>
#include <iostream>
#include "utilities.h"

// this class defines the matrix (n x m)
//
class matrix
{
    public:
// constructors/destructor declarations
        matrix();
        matrix(int n, int m);
        matrix(int n);
        matrix(const matrix& other);
        matrix& operator=(const matrix& other);
        ~matrix();
// public methods declarations
        matrix eByE(const matrix& m) const; //element by element product
        matrix cByC(const matrix& m, int c) const; //column by column product
        matrix col(int j) const;
        void combineRows(int i1, double c1, int i2, double c2, int i3);
        double complement(int i, int j) const;
        double det() const;
        void dim(int n, int m);
        void dim(int n);
        double get(int i,int j) const;
        double get(int i) const;
        void identity();
        matrix inv() const;
        double mean() const;
        matrix minor(int i, int j) const;
        int nR() const;
        int nC() const;
        matrix part(int i, int j, int h, int l) const; 
        void pivot(int i);
        matrix pw(int n) const;
        void print() const;
        void print(const std::string s) const;
        matrix row(int i) const;
        void set(int i, int j, double aij);
        void set(int i, double ai);
        void add(int i, int j, double aij);
        void add(int i, double ai);
        void add(int i, int j, matrix m);
        void set(int i, int j, matrix m);
        int size() const;
        void switchRows(int i1, int i2);
        matrix T() const;
        double trace() const;
        void zero();
        void operator+=(const matrix& m);
        void operator-=(const matrix& m);
        matrix operator*(const matrix& m) const;
        matrix operator/(const matrix& m) const;
        matrix operator+(const matrix& m) const;
        matrix operator-(const matrix& m) const;
        matrix operator*(const double c) const;
        matrix operator/(const double c) const;
        matrix operator*(const int c) const;
        matrix operator/(const int c) const;
        matrix operator-() const;
        friend matrix operator*(double c, const matrix& m);
        friend matrix operator*(int c, const matrix& m);
    private:
// private variable declarations
        int nr=0, nc=0;
        double* a=nullptr;       
};
//
// this class defines the matrix (n x m) of integer numbers
//
class intMatrix
{
    public:
// constructors/destructor declarations
        intMatrix();
        intMatrix(int n, int m);
        intMatrix(int n);
        intMatrix(const intMatrix& other);
        intMatrix& operator=(const intMatrix& other);
        ~intMatrix();
// public methods declarations
        void dim(int n, int m);
        void dim(int n);
        int get(int i,int j) const;
        int get(int i) const;
        int nR() const;
        int nC() const;
        void print() const;
        void set(int i, int j, int aij);
        void set(int i, int ai);
        int size() const;
        void zero();
    private:
// private variable declarations
        int nr=0, nc=0;
        int* a=nullptr;       
};
#endif
