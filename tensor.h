#ifndef TENSOR_H
#define TENSOR_H
// matrix.h
#include<cmath>
#include <cstdlib>
#include <iostream>
#include "utilities.h"
//
// this class defines the symmetrical tensor (3 x 3)
//
class tensor; //pre-declaration to avoid incongruences
//
class symTensor
{
    public:
// constructors/destructor declarations
        symTensor();
        symTensor(double a00, double a10, double a11, double a20, double a21, double a22);
        symTensor(const symTensor& other);
        symTensor& operator=(const symTensor& other);
        // public methods declarations
        void add(int i, int j, double aij);
        symTensor eByE(const symTensor& t) const; //element by element product
        double autoContr() const;
        double doubleContr(const symTensor& t) const;
        double doubleContr(const tensor& t) const;
        double Frobenius() const;
        double get(int i,int j) const;
        symTensor noTrace() const;
        void print() const;
        void set(int i, int j, double aij);
        void set(double a00, double a10, double a11, double a20, double a21, double a22);
        double trace() const;
        void zero();
        void operator+=(const symTensor& t);
        void operator-=(const symTensor& t);
        symTensor operator*(const double c) const;
        symTensor operator/(const double c) const;
        symTensor operator+(const symTensor& t) const;
        symTensor operator-(const symTensor& t) const;
        symTensor operator*(const int c) const;
        symTensor operator/(const int c) const;
        symTensor operator-() const;
        symTensor operator*(const symTensor& t) const;
        tensor operator*(const tensor& t) const;
        friend symTensor operator*(double c, const symTensor& t);
        friend symTensor operator*(int c, const symTensor& t);
    private:
// private methods declarations
        double get(int i) const;
// private variable declarations
        double a[6];       
};
//
// this class defines the tensor (3 x 3)
//
class tensor
{
    public:
// constructors/destructor declarations
        tensor();
        tensor(double a00, double a01, double a02, double a10, double a11, double a12, double a20, double a21, double a22);
        tensor(const tensor& other);
        tensor& operator=(const tensor& other);
        void operator=(const symTensor& other);
// public methods declarations
        void add(int i, int j, double aij);
        tensor eByE(const tensor& t) const; //element by element product
        double autoContr() const;
        double doubleContr(const tensor& t) const;
        double doubleContr(const symTensor& t) const;
        double Frobenius() const;
        double get(int i,int j) const;
        void print() const;
        void set(int i, int j, double aij);
        void set(double a00, double a01, double a02, double a10, double a11, double a12, double a20, double a21, double a22);
        symTensor symmetrize() const;
        double trace() const;
        void zero();
        void operator+=(const tensor& t);
        void operator-=(const tensor& t);
        tensor operator*(const double c) const;
        tensor operator/(const double c) const;
        tensor operator+(const tensor& t) const;
        tensor operator-(const tensor& t) const;
        tensor operator*(const int c) const;
        tensor operator/(const int c) const;
        tensor operator-() const;
        tensor operator*(const tensor& t) const; // double contraction
        tensor operator*(const symTensor& t) const; // double contraction
        friend tensor operator*(double c, const tensor& t);
        friend tensor operator*(int c, const tensor& t);
    private:
// private methods declarations
        double get(int i) const;
// private variable declarations
        double a[9];       
};
#endif
