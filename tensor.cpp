#include "tensor.h"
// constructors of symTensor
symTensor::symTensor()
{
    zero();
}
symTensor::symTensor(double a00, double a10, double a11, double a20, double a21, double a22)
{
    a[0]=a00; a[1]=a10; a[2]=a11; a[3]=a20; a[4]=a21; a[5]=a22;
}
// copy constructors
symTensor::symTensor(const symTensor& other)
{
    a[0]=other.a[0]; a[1]=other.a[1]; a[2]=other.a[2];
    a[3]=other.a[3]; a[4]=other.a[4]; a[5]=other.a[5];
}
symTensor& symTensor::operator=(const symTensor& other)
{
    if (this != &other)
    {
        a[0]=other.a[0]; a[1]=other.a[1]; a[2]=other.a[2];
        a[3]=other.a[3]; a[4]=other.a[4]; a[5]=other.a[5];
    }
    return *this;
}
// methods of symTensor
void symTensor::add(int i, int j, double aij)
{
//    a[i+j+max(i,j)/2]+=aij;
    int k=i; if (j>i) {k=j;}
    a[i+j+k/2]+=aij;
}
symTensor symTensor::eByE(const symTensor& t) const
{
    symTensor tens;
    for(int i=0; i<6; i++)
    {
        tens.a[i]=a[i]*t.a[i];
    }
    return tens;
}
double symTensor::autoContr() const
{
    return sq(a[0])+sq(a[2])+sq(a[5])+2.*(sq(a[1])+sq(a[3])+sq(a[4]));
}
double symTensor::Frobenius() const
{
    return sqrt(autoContr());
}
double symTensor::doubleContr(const symTensor& t) const
{
    return a[0]*t.a[0]+a[2]*t.a[2]+a[5]*t.a[5]+2.*(a[1]*t.a[1]+a[3]*t.a[3]+a[4]*t.a[4]);
}
double symTensor::doubleContr(const tensor& t) const
{
    return a[0]*t.get(0,0)+a[1]*t.get(1,0)+a[2]*t.get(2,0)+a[3]*t.get(0,1)+a[4]*t.get(1,1)+a[5]*t.get(2,1)+a[6]*t.get(0,2)+a[7]*t.get(1,2)+a[8]*t.get(2,2);
}
double symTensor::get(int i,int j) const
{
    int k=i; if (j>i) {k=j;}
//    return a[i+j+max(i,j)/2];
    return a[i+j+k/2];
}
double symTensor::get(int i) const
{
    return a[i];
}
symTensor symTensor::noTrace() const
{
    symTensor t=(*this);
    double tr=-trace()/3.;
    t.a[0]+=tr; t.a[2]+=tr; t.a[5]+=tr;
    return t;
}
void symTensor::print() const
{
    std::cout << "==================" << std::endl;
    std::cout << 3 << " " << 3 << std::endl;
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            if (j<2) {std::cout << get(i,j) << " ";}
            else {std::cout << get(i,j) << std::endl;}
        }
    }
    std::cout << "==================" << std::endl;
}
void symTensor::set(int i, int j, double aij)
{
    int k=i; if (j>i) {k=j;}
//    a[i+j+max(i,j)/2]=aij;
    a[i+j+k/2]=aij;
}
void symTensor::set(double a00, double a10, double a11, double a20, double a21, double a22)
{
    a[0]=a00; a[1]=a10; a[2]=a11; a[3]=a20; a[4]=a21; a[5]=a22;
}
double symTensor::trace() const
{
    double t=a[0]+a[2]+a[5];
    return t;
}
void symTensor::zero()
{
    a[0]=0.; a[1]=0.; a[2]=0.; a[3]=0.; a[4]=0.; a[5]=0.;
}

// operators
void symTensor::operator+=(const symTensor& t)
{
    for (int i=0; i<6; i++)
    {
        a[i]+=t.a[i];
    }    
}
void symTensor::operator-=(const symTensor& t)
{
    for (int i=0; i<6; i++)
    {
        a[i]-=t.a[i];
    }    
}
symTensor symTensor::operator*(const double c) const
{
    symTensor t;
    for(int i=0; i<6; i++)
    {
        t.a[i]=a[i]*c;
    }
    return t;
}
symTensor symTensor::operator/(const double c) const
{
    symTensor t;
    for(int i=0; i<6; i++)
    {
        t.a[i]=a[i]/c;
    }
    return t;
}
symTensor symTensor::operator+(const symTensor& t) const
{
    symTensor ten;
    for(int i=0; i<6; i++)
    {
        ten.a[i]=a[i]+t.a[i];
    }
    return ten;
}
symTensor symTensor::operator-(const symTensor& t) const
{
    symTensor ten;
    for(int i=0; i<6; i++)
    {
        ten.a[i]=a[i]-t.a[i];
    }
    return ten;
}
symTensor symTensor::operator*(const int c) const
{
    symTensor t;
    for(int i=0; i<6; i++)
    {
        t.a[i]=a[i]*c;
    }
    return t;
}
symTensor symTensor::operator/(const int c) const
{
    symTensor t;
    for(int i=0; i<6; i++)
    {
        t.a[i]=a[i]/c;
    }
    return t;
}
symTensor symTensor::operator-() const
{
    return (*this)*(-1);
}
symTensor symTensor::operator*(const symTensor& t) const
{
    symTensor tens;
    tens.a[0]=a[0]*t.a[0]+a[1]*t.a[1]+a[3]*t.a[3];
    tens.a[1]=a[1]*t.a[0]+a[2]*t.a[1]+a[4]*t.a[3];
    tens.a[2]=a[1]*t.a[1]+a[2]*t.a[2]+a[4]*t.a[4];
    tens.a[3]=a[3]*t.a[0]+a[4]*t.a[1]+a[5]*t.a[3];
    tens.a[4]=a[3]*t.a[1]+a[4]*t.a[2]+a[5]*t.a[4];
    tens.a[5]=a[3]*t.a[3]+a[4]*t.a[4]+a[5]*t.a[5];
    return tens;
}
symTensor operator*(double c, const symTensor& t)
{
    symTensor tens;
    for(int i=0; i<6; i++)
    {
        tens.a[i]=t.a[i]*c;
    }
    return tens;
}
symTensor operator*(int c, const symTensor& t)
{
    symTensor tens;
    for(int i=0; i<6; i++)
    {
        tens.a[i]=t.a[i]*c;
    }
    return tens;
}
//
//
// constructors of tensor
tensor::tensor()
{
    zero();
}
tensor::tensor(double a00, double a01, double a02, double a10, double a11, double a12, double a20, double a21, double a22)
{
    a[0]=a00; a[1]=a01; a[2]=a02; a[3]=a10; a[4]=a11; a[5]=a12; a[6]=a20; a[7]=a21; a[8]=a22;
}
// copy constructors
tensor::tensor(const tensor& other)
{
    a[0]=other.a[0]; a[1]=other.a[1]; a[2]=other.a[2];
    a[3]=other.a[3]; a[4]=other.a[4]; a[5]=other.a[5];
    a[6]=other.a[6]; a[7]=other.a[7]; a[8]=other.a[8];
}
tensor& tensor::operator=(const tensor& other)
{
    if (this != &other)
    {
        a[0]=other.a[0]; a[1]=other.a[1]; a[2]=other.a[2];
        a[3]=other.a[3]; a[4]=other.a[4]; a[5]=other.a[5];
        a[6]=other.a[6]; a[7]=other.a[7]; a[8]=other.a[8];
    }
    return *this;
}
void tensor::operator=(const symTensor& other)
{
    {
        a[0]=other.get(0,0); a[1]=other.get(0,1); a[2]=other.get(0,2);
        a[3]=other.get(1,0); a[4]=other.get(1,1); a[5]=other.get(1,2);
        a[6]=other.get(2,0); a[7]=other.get(2,1); a[8]=other.get(2,2);
    }
}
// methods of tensor
void tensor::add(int i, int j, double aij)
{
    a[3*i+j]+=aij;
}
tensor tensor::eByE(const tensor& t) const
{
    tensor tens;
    for(int i=0; i<9; i++)
    {
        tens.a[i]=a[i]*t.a[i];
    }
    return tens;
}
double tensor::autoContr() const
{
    return sq(a[0])+sq(a[1])+sq(a[2])+sq(a[3])+sq(a[4])+sq(a[5])+sq(a[6])+sq(a[7])+sq(a[8]);
}
double tensor::doubleContr(const tensor& t) const
{
    return a[0]*t.a[0]+a[1]*t.a[3]+a[2]*t.a[6]+a[3]*t.a[1]+a[4]*t.a[4]+a[5]*t.a[7]+a[6]*t.a[2]+a[7]*t.a[5]+a[8]*t.a[8];
}
double tensor::doubleContr(const symTensor& t) const
{
    return a[0]*t.get(0,0)+a[1]*t.get(1,0)+a[2]*t.get(2,0)+a[3]*t.get(0,1)+a[4]*t.get(1,1)+a[5]*t.get(2,1)+a[6]*t.get(0,2)+a[7]*t.get(1,2)+a[8]*t.get(2,2);
}
double tensor::Frobenius() const
{
    return sqrt(autoContr());
}
double tensor::get(int i,int j) const
{
    return a[3*i+j];
}
double tensor::get(int i) const
{
    return a[i];
}
void tensor::print() const
{
    std::cout << "==================" << std::endl;
    std::cout << 3 << " " << 3 << std::endl;
    std::cout << a[0] << " " << a[1] << " " << a[2] << std::endl; 
    std::cout << a[3] << " " << a[4] << " " << a[5] << std::endl; 
    std::cout << a[6] << " " << a[7] << " " << a[8] << std::endl; 
    std::cout << "==================" << std::endl;
}
void tensor::set(int i, int j, double aij)
{
    a[3*i+j]=aij;
}
void tensor::set(double a00, double a01, double a02, double a10, double a11, double a12, double a20, double a21, double a22)
{
    a[0]=a00; a[1]=a01; a[2]=a02; a[3]=a10; a[4]=a11; a[5]=a12; a[6]=a20; a[7]=a21; a[8]=a22;
}
symTensor tensor::symmetrize() const
{
    symTensor t(a[0],(a[3]+a[1])/2.,a[4],(a[6]+a[2])/2.,(a[7]+a[5])/2.,a[8]);
    return t;    
}
double tensor::trace() const
{
    double t=a[0]+a[4]+a[8];
    return t;
}
void tensor::zero()
{
    a[0]=0.; a[1]=0.; a[2]=0.;
    a[3]=0.; a[4]=0.; a[5]=0.;
    a[6]=0.; a[7]=0.; a[8]=0.;
}

// operators
void tensor::operator+=(const tensor& t)
{
    for (int i=0; i<9; i++)
    {
        a[i]+=t.a[i];
    }    
}
void tensor::operator-=(const tensor& t)
{
    for (int i=0; i<9; i++)
    {
        a[i]-=t.a[i];
    }    
}
tensor tensor::operator*(const double c) const
{
    tensor t;
    for(int i=0; i<9; i++)
    {
        t.a[i]=a[i]*c;
    }
    return t;
}
tensor tensor::operator/(const double c) const
{
    tensor t;
    for(int i=0; i<9; i++)
    {
        t.a[i]=a[i]/c;
    }
    return t;
}
tensor tensor::operator+(const tensor& t) const
{
    tensor ten;
    for(int i=0; i<9; i++)
    {
        ten.a[i]=a[i]+t.a[i];
    }
    return ten;
}
tensor tensor::operator-(const tensor& t) const
{
    tensor ten;
    for(int i=0; i<9; i++)
    {
        ten.a[i]=a[i]-t.a[i];
    }
    return ten;
}
tensor tensor::operator*(const int c) const
{
    tensor t;
    for(int i=0; i<9; i++)
    {
        t.a[i]=a[i]*c;
    }
    return t;
}
tensor tensor::operator/(const int c) const
{
    tensor t;
    for(int i=0; i<9; i++)
    {
        t.a[i]=a[i]/c;
    }
    return t;
}
tensor tensor::operator-() const
{
    return (*this)*(-1);
}
tensor tensor::operator*(const tensor& t) const
{
    tensor tens;
    tens.a[0]=a[0]*t.a[0]+a[1]*t.a[3]+a[2]*t.a[6];
    tens.a[1]=a[0]*t.a[1]+a[1]*t.a[4]+a[2]*t.a[7];
    tens.a[2]=a[0]*t.a[2]+a[1]*t.a[5]+a[2]*t.a[8];
    tens.a[3]=a[3]*t.a[0]+a[4]*t.a[3]+a[5]*t.a[6];
    tens.a[4]=a[3]*t.a[1]+a[4]*t.a[4]+a[5]*t.a[7];
    tens.a[5]=a[3]*t.a[2]+a[4]*t.a[5]+a[5]*t.a[8];
    tens.a[6]=a[6]*t.a[0]+a[7]*t.a[3]+a[8]*t.a[6];
    tens.a[7]=a[6]*t.a[1]+a[7]*t.a[4]+a[8]*t.a[7];
    tens.a[8]=a[6]*t.a[2]+a[7]*t.a[5]+a[8]*t.a[8];    
    return tens;
}
tensor operator*(double c, const tensor& t)
{
    tensor tens;
    for(int i=0; i<9; i++)
    {
        tens.a[i]=t.a[i]*c;
    }
    return tens;
}
tensor operator*(int c, const tensor& t)
{
    tensor tens;
    for(int i=0; i<9; i++)
    {
        tens.a[i]=t.a[i]*c;
    }
    return tens;
}


