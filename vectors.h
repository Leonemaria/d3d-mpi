#ifndef VECTORS_H
#define VECTORS_H
// vectors.h
#include <iostream>
#include <cstdlib>
#include <fstream>
#include<cmath>

class vector2D
{
    public:
// constructors
        vector2D();
        vector2D(double c1,double c2);
        vector2D(const vector2D& other);
        vector2D& operator=(const vector2D& other);
// methods
        double norm() const;
        void normalize();
        vector2D normalized() const;
        void set(double c1, double c2);
        void zero();
        double dot(const vector2D& v) const;
        double operator[](int i) const;
        double& operator[](int i);
        vector2D operator-();
        void print() const;
    private:
        double a[2];
};
class vector3D
{
    public:
// constructors
        vector3D();
        vector3D(double c1, double c2, double c3);
        vector3D(const vector3D& other);
        vector3D& operator=(const vector3D& other);
// methods
        vector3D cross(const vector3D& v);
        double norm() const;
        void normalize();
        vector3D normalized() const;
        void set(double c1, double c2, double c3);
        void zero();
        double dot(const vector3D& v) const;
        double operator[](int i) const;
        double& operator[](int i);
        vector3D operator-();
        void print() const;
    private:
        double a[3];
};
vector2D just2D(vector3D v);
vector3D as3D(vector2D v);
//
vector2D operator*(const double c, const vector2D& v);
vector2D operator*(const vector2D& v, const double c);
vector2D operator/(const vector2D& v, const double c);
vector2D operator+(const vector2D& v1, const vector2D& v2);
vector2D operator-(const vector2D& v1, const vector2D& v2);
//
vector3D operator*(const double c, const vector3D& v);
vector3D operator*(const vector3D& v, const double c);
vector3D operator/(const vector3D& v, const double c);
vector3D operator+(const vector3D& v1, const vector3D& v2);
vector3D operator-(const vector3D& v1, const vector3D& v2);
//
#endif
