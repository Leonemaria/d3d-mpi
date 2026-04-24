//matricesAndVectors.cpp
#include "vectors.h"

vector2D::vector2D()
{
    zero();
}
vector2D::vector2D(double c1,double c2)
{
    a[0]=c1; a[1]=c2;
}
// copy constructors
vector2D::vector2D(const vector2D& other)
{
    a[0]=other[0];
    a[1]=other[1];
}
vector2D& vector2D::operator=(const vector2D& other)
{
    a[0]=other[0];
    a[1]=other[1];
    return *this;
}
// methods
double vector2D::norm() const
{
    double n;
    n=sqrt(a[0]*a[0]+a[1]*a[1]);
    return n;
}
void vector2D::normalize()
{
    a[0]=a[0]/norm();
    a[1]=a[1]/norm();
}
vector2D vector2D::normalized() const
{
    vector2D vec;
    vec[0]=a[0]/norm();
    vec[1]=a[1]/norm();
    return vec;
};
void vector2D::set(double c1, double c2)
{
    a[0]=c1; a[1]=c2;
}
void vector2D::zero()
{
    a[0]=0.; a[1]=0.;
}
double vector2D::dot(const vector2D& v) const
{
    return a[0]*v[0]+a[1]*v[1];
}
double vector2D::operator[](int i) const
{
    return a[i];
}
double& vector2D::operator[](int i)
{
    return a[i];
}
vector2D vector2D::operator-()
{
    return (*this)*(-1);
}
void vector2D::print() const
{
    std::cout << a[0] << "  " << a[1] << std::endl;
}
//
vector3D::vector3D()
{
    zero();
}
vector3D::vector3D(double c1, double c2, double c3)
{
    a[0]=c1; a[1]=c2; a[2]=c3;
}
// copy constructors
vector3D::vector3D(const vector3D& other)
{
    a[0]=other[0];
    a[1]=other[1];
    a[2]=other[2];
}
vector3D& vector3D::operator=(const vector3D& other)
{
    a[0]=other[0];
    a[1]=other[1];
    a[2]=other[2];
    return *this;
}
// methods
vector3D vector3D::cross(const vector3D& v)
{
    vector3D vc;
    vc[0]=a[1]*v[2]-a[2]*v[1];
    vc[1]=a[2]*v[0]-a[0]*v[2];
    vc[2]=a[0]*v[1]-a[1]*v[0];
    return vc;
}


double vector3D::norm() const
{
    double n;
    n=sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
    return n;
}
void vector3D::normalize()
{
    a[0]=a[0]/norm();
    a[1]=a[1]/norm();
    a[2]=a[2]/norm();
}
vector3D vector3D::normalized() const
{
    vector3D vec;
    vec[0]=a[0]/norm();
    vec[1]=a[1]/norm();
    vec[2]=a[2]/norm();
    return vec;
};
void vector3D::set(double c1, double c2, double c3)
{
    a[0]=c1; a[1]=c2; a[2]=c3;
}
void vector3D::zero()
{
    a[0]=0.; a[1]=0.; a[2]=0.;
}
double vector3D::dot(const vector3D& v) const
{
    return a[0]*v[0]+a[1]*v[1]+a[2]*v[2];
}
double vector3D::operator[](int i) const
{
    return a[i];
}
double& vector3D::operator[](int i)
{
    return a[i];
}
vector3D vector3D::operator-()
{
    return (*this)*(-1);
}
void vector3D::print() const
{
    std::cout << a[0] << "  " << a[1] << "  " << a[2] << std::endl;
}
// non element functions
vector3D as3D(vector2D v)
{
    vector3D va;
    va[0]=v[0]; va[1]=v[1]; va[2]=0.;
    return va;
}
vector2D just2D(vector3D v)
{
    vector2D vj;
    vj[0]=v[0]; vj[1]=v[1];
    return vj;
}
// non element operators
vector2D operator*(const double c, const vector2D& v)
{
    vector2D vec;
    vec[0]=c*v[0];
    vec[1]=c*v[1];
    return vec;
}
vector2D operator*(const vector2D& v, const double c)
{
    vector2D vec;
    vec[0]=c*v[0];
    vec[1]=c*v[1];
    return vec;
};
vector2D operator/(const vector2D& v, const double c)
{
    vector2D vec;
    vec[0]=v[0]/c;
    vec[1]=v[1]/c;
    return vec;
};
vector2D operator+(const vector2D& v1, const vector2D& v2)
{
    vector2D vec;
    vec[0]=v1[0]+v2[0];
    vec[1]=v1[1]+v2[1];
    return vec;
}
vector2D operator-(const vector2D& v1, const vector2D& v2)
{
    vector2D vec;
    vec[0]=v1[0]-v2[0];
    vec[1]=v1[1]-v2[1];
    return vec;
}
//
vector3D operator*(const double c, const vector3D& v)
{
    vector3D vec;
    vec[0]=c*v[0];
    vec[1]=c*v[1];
    vec[2]=c*v[2];
    return vec;
}
vector3D operator*(const vector3D& v, const double c)
{
    vector3D vec;
    vec[0]=c*v[0];
    vec[1]=c*v[1];
    vec[2]=c*v[2];
    return vec;
};
vector3D operator/(const vector3D& v, const double c)
{
    vector3D vec;
    vec[0]=v[0]/c;
    vec[1]=v[1]/c;
    vec[2]=v[2]/c;
    return vec;
};
vector3D operator+(const vector3D& v1, const vector3D& v2)
{
    vector3D vec;
    vec[0]=v1[0]+v2[0];
    vec[1]=v1[1]+v2[1];
    vec[2]=v1[2]+v2[2];
    return vec;
}
vector3D operator-(const vector3D& v1, const vector3D& v2)
{
    vector3D vec;
    vec[0]=v1[0]-v2[0];
    vec[1]=v1[1]-v2[1];
    vec[2]=v1[2]-v2[2];
    return vec;
}



