//matrix.cpp
#include "matrix.h"

// constructors of matrix
matrix::matrix()
{
    a=nullptr;
}
matrix::matrix(int n, int m)
{
    dim(n,m);
}
matrix::matrix(int n)
{
    dim(n);
    zero();
}
// copy constructors
matrix::matrix(const matrix& other)
{
    dim(other.nR(),other.nC());
    for (int i=0; i<size(); i++)
    {
        a[i]=other.a[i];
    }
}
matrix& matrix::operator=(const matrix& other)
{
    if (this != &other)
    {
        dim(other.nR(),other.nC());
        for (int i=0; i<size(); i++)
        {
            a[i]=other.a[i];
        }
    }
    return *this;
}
// destructor of matrix
matrix::~matrix()
{
    delete[] a; a=nullptr;
}
// methods of matrix
matrix matrix::eByE(const matrix& m) const
{
    matrix mat(nr,nc);
    for(int i=0; i<(nr*nc); i++)
    {
        mat.a[i]=a[i]*m.a[i];
    }
    return mat;
}
matrix matrix::cByC(const matrix& m, int c) const
{
    matrix mat(nr,nc);
    for (int i=0; i<nr; i++)
    {
        for (int j=0; j<nc; j++)
        {
            mat.set(i,j,a[nc*i+j]*(m.col(c)).a[i]);
        }
    }
    return mat;
}
matrix matrix::col(int j) const
{
    matrix c(nr,1);
    for (int i=0; i<nr; i++)
    {
        c.set(i,0,a[nc*i+j]);
    }
    return c;
}
void matrix::combineRows(int i1, double c1, int i2, double c2, int i3)
{
    double b;
    for (int j=0; j<nc; j++)
    {
        b=get(i1,j)*c1+get(i2,j)*c2;
        set(i3,j,b);
    }
}
double matrix::complement(int i, int j) const
{
    return pow(-1,i+j)*minor(i,j).det();
}
double matrix::det() const
{
    double d=0.;
    if (nr==1) {d=get(0,0);}
    else
    {
        if (nr==2) {d=get(0,0)*get(1,1)-get(0,1)*get(1,0);}
        else
        {
            int i=0;
            for(int j=0; j<nc; j++)
            {
                d+=get(i,j)*complement(i,j);
            }
        }
    }
    return d;
}        
void matrix::dim(int n, int m)
{
    if (size()>0) {delete[] a; a=nullptr;}
    nr=n; nc=m;
    a=new double[nc*nr];
    zero();
}
void matrix::dim(int n)
{
    if (size()>0) {delete[] a; a=nullptr;}
    nr=n; nc=1;
    a=new double[nr];
    zero();
}
double matrix::get(int i,int j) const
{
    return a[nc*i+j];
}
double matrix::get(int i) const
{
    return a[i];
}
void matrix::identity()
{
    zero();
    for (int i=0; i<nr; i++)
    {
        set(i,i,1.);    
    }
}
matrix matrix::inv() const
{
    matrix M(nr,nc*2); matrix m(nr,nc);
    M.zero(); double c;
    for (int i=0; i<nr; i++)
    {
        M.set(i,nc+i,1.);
        for (int j=0; j<nc; j++)
        {
            M.set(i,j,get(i,j));
        }
    }
    M.pivot(0); c=M.get(0,0);
    for (int j=0; j<nc*2; j++)
    {
        M.set(0,j,M.get(0,j)/c);
    }
    for (int i=1; i<nr; i++)
    {
        for (int ii=i; ii<nr; ii++)
        {
            M.combineRows(i-1,M.get(ii,i-1),ii,-M.get(i-1,i-1),ii);
        }
        M.pivot(i); c=M.get(i,i);
        for (int j=0; j<nc*2; j++)
        {
            M.set(i,j,M.get(i,j)/c);
        }
    }
    for (int i=nr-2; i>-1; i--)
    {
        for (int ii=i; ii>-1; ii--)
        {
            M.combineRows(i+1,M.get(ii,i+1),ii,-M.get(i+1,i+1),ii);
        }
        c=M.get(i,i);
        for (int j=0; j<nc*2; j++)
        {
            M.set(i,j,M.get(i,j)/c);
        }
    }
    for (int i=0; i<nr; i++)
    {
        for (int j=0; j<nc; j++)
        {
            m.set(i,j,M.get(i,nc+j));
        }
    }
    return m;  
}
double matrix::mean() const
{
    double m=0.;
    for(int i=0; i<size(); i++)
    {
        m+=a[i];
    }
    return m/size();
}
matrix matrix::minor(int i, int j) const
{
    matrix m(nr-1,nc-1);
    int di=0, dj;
    for (int ii=0; ii<nr-1; ii++)
    {
        if (ii>=i) {di=1;} else {di=0;}
        dj=0;
        for (int jj=0; jj<nc-1; jj++)
        {
            if (jj>=j) {dj=1;} else {dj=0;}
            m.set(ii,jj,get(ii+di,jj+dj));   
        }
    }
    return m;
}
int matrix::nC() const
{
    return nc;
}
int matrix::nR() const
{
    return nr;
}
matrix matrix::part(int i, int j, int h, int l) const
{
    matrix m(h,l);
    for (int ii=0; ii<h; ii++)
    {
        for (int jj=0; jj<l; jj++)
        {
            m.set(ii,jj,get(i+ii,j+jj));
        }
    }
    return m;
}
void matrix::pivot(int i)
{
    double mx=abs(get(i,i)); int imx=i;
    for (int ii=i+1; ii<nr; ii++)
    {
        if (abs(get(ii,i))>mx) {mx=abs(get(ii,i)); imx=ii;}
    }
    if (imx!=i) {switchRows(i,imx);}
}
matrix matrix::pw(int n) const
{
    matrix mat(nr,nc);
    for(int i=0; i<(nr*nc); i++)
    {
        mat.a[i]=pow(a[i],n);
    }
    return mat;
}
void matrix::print() const
{
    std::cout << "==================" << std::endl;
    std::cout << nr << " " << nc << std::endl;
    for (int i=0; i<nr; i++)
    {
        for (int j=0; j<nc; j++)
        {
            if (j<nc-1) {std::cout << get(i,j) << " ";}
            else {std::cout << get(i,j) << std::endl;}
        }
    }
    std::cout << "==================" << std::endl;
}
void matrix::print(const std::string s) const
{
    std::cout << "===== " << s << " ======" << std::endl;
    std::cout << nr << " " << nc << std::endl;
    for (int i=0; i<nr; i++)
    {
        for (int j=0; j<nc; j++)
        {
            if (j<nc-1) {std::cout << get(i,j) << " ";}
            else {std::cout << get(i,j) << std::endl;}
        }
    }
    std::cout << "==================" << std::endl;
}
matrix matrix::row(int i) const
{
    matrix r(1,nc);
    for (int j=0; j<nc; j++)
    {
        r.set(0,j,a[nc*i+j]);
    }
    return r;
}
void matrix::set(int i, int j, double aij)
{
    a[nc*i+j]=aij;
}
void matrix::set(int i, double ai)
{
    a[i]=ai;
}
void matrix::add(int i, int j, double aij)
{
    a[nc*i+j]+=aij;
}
void matrix::add(int i, double ai)
{
    a[i]+=ai;
}
void matrix::add(int i0, int j0, matrix m)
{
    for (int i=0; i<m.nR(); i++)
    {
        for (int j=0; j<m.nC(); j++)
        {
            a[nc*(i0+i)+j0+j]+=m.a[m.nc*i+j];
        }
    }
}
void matrix::set(int i0, int j0, matrix m)
{
    for (int i=0; i<m.nR(); i++)
    {
        for (int j=0; j<m.nC(); j++)
        {
            a[nc*(i0+i)+j0+j]=m.a[m.nC()*i+j];
        }
    }
}

int matrix::size() const
{
    return nr*nc;
}
matrix matrix::T() const
{
    matrix mat(nc,nr);
    for (int i=0; i<nc; i++)
    {
        for (int j=0; j<nr; j++)
        {
            mat.a[nr*i+j]=a[nc*j+i];
        }
    }
    return mat;
}
void matrix::switchRows(int i1, int i2)
{
    double b;
    for (int j=0; j<nc; j++)
    {
        b=a[nc*i1+j];
        a[nc*i1+j]=a[nc*i2+j];
        a[nc*i2+j]=b;
    }
}
double matrix::trace() const
{
    double t=0.; int ncp=nc+1;
    for (int i=0; i<nc; i++)
    {
        t+=a[ncp*i];
    }
    return t;
}
void matrix::zero()
{
    for (int i=0; i<size(); i++)
    {
        a[i]=0.;
    }
}

// operators
void matrix::operator+=(const matrix& m)
{
    for (int i=0; i<nr*nc; i++)
    {
        a[i]+=m.a[i];
    }    
}
void matrix::operator-=(const matrix& m)
{
    for (int i=0; i<nr*nc; i++)
    {
        a[i]-=m.a[i];
    }    
}
matrix matrix::operator*(const matrix& m) const
{
    double e;
    matrix mat(nr,m.nC());
    int ij=0;
    int ik0=0, ik;
    for(int i=0; i<nr; i++)
    {
        for(int j=0; j<m.nC(); j++)
        {
             ik=ik0;
             e=0.;
             for(int k=0; k<nc; k++)
             {
                 e+=a[ik]*m.a[m.nC()*k+j];
                 ik++;
             }
             mat.a[ij]=e;
             ij++;
        }
        ik0+=nc;
    }
    return mat;
}
matrix matrix::operator/(const matrix& m) const
{
    return m.inv()*(*this);
}
matrix matrix::operator+(const matrix& m) const
{
    matrix mat(nr,nc);
    for(int i=0; i<(nr*nc); i++)
    {
        mat.a[i]=a[i]+m.a[i];
    }
    return mat;
}
matrix matrix::operator-(const matrix& m) const
{
    matrix mat(nr,nc);
    for(int i=0; i<(nr*nc); i++)
    {
        mat.a[i]=a[i]-m.a[i];
    }
    return mat;
}
matrix matrix::operator*(const double c) const
{
    matrix mat(nr,nc);
    for(int i=0; i<(nr*nc); i++)
    {
        mat.a[i]=a[i]*c;
    }
    return mat;
}
matrix matrix::operator/(const double c) const
{
    matrix mat(nr,nc);
    for(int i=0; i<(nr*nc); i++)
    {
        mat.a[i]=a[i]/c;
    }
    return mat;
}
matrix matrix::operator*(const int c) const
{
    matrix mat(nr,nc);
    for(int i=0; i<(nr*nc); i++)
    {
        mat.a[i]=a[i]*c;
    }
    return mat;
}
matrix matrix::operator/(const int c) const
{
    matrix mat(nr,nc);
    for(int i=0; i<(nr*nc); i++)
    {
        mat.a[i]=a[i]/c;
    }
    return mat;
}
matrix matrix::operator-() const
{
    return (*this)*(-1);
}
matrix operator*(double c, const matrix& m)
{
    matrix mat(m.nR(),m.nC());
    for(int i=0; i<m.size(); i++)
    {
        mat.a[i]=m.a[i]*c;
    }
    return mat;
}
matrix operator*(int c, const matrix& m)
{
    matrix mat(m.nR(),m.nC());
    for(int i=0; i<m.size(); i++)
    {
        mat.a[i]=m.a[i]*c;
    }
    return mat;
}





// constructors of intMatrix
intMatrix::intMatrix()
{
}
intMatrix::intMatrix(int n, int m)
{
    dim(n,m);
}
intMatrix::intMatrix(int n)
{
    dim(n);
}
// copy constructors
intMatrix::intMatrix(const intMatrix& other)
{
    dim(other.nR(),other.nC());
    for (int i=0; i<size(); i++)
    {
        a[i]=other.a[i];
    }
}
intMatrix& intMatrix::operator=(const intMatrix& other)
{
    if (this != &other)
    {
        dim(other.nR(),other.nC());
        for (int i=0; i<size(); i++)
        {
            a[i]=other.a[i];
        }
    }
    return *this;
}
// destructor of intMatrix
intMatrix::~intMatrix()
{
    delete[] a; a=nullptr;
}
// methods of intMatrix
void intMatrix::dim(int n, int m)
{
    if (size()>0) {delete[] a; a=nullptr;}
    nr=n; nc=m;
    a=new int[nc*nr];
}
void intMatrix::dim(int n)
{
    if (size()>0) {delete[] a; a=nullptr;}
    nr=n; nc=1;
    a=new int[nr];
}
int intMatrix::get(int i,int j) const
{
    return a[nc*i+j];
}
int intMatrix::get(int i) const
{
    return a[i];
}
int intMatrix::nC() const
{
    return nc;
}
int intMatrix::nR() const
{
    return nr;
}
void intMatrix::print() const
{
    std::cout << "==================" << std::endl;
    std::cout << nr << " " << nc << std::endl;
    for (int i=0; i<nr; i++)
    {
        for (int j=0; j<nc; j++)
        {
            if (j<nc-1) {std::cout << get(i,j) << " ";}
            else {std::cout << get(i,j) << std::endl;}
        }
    }
    std::cout << "==================" << std::endl;
}
void intMatrix::set(int i, int j, int aij)
{
    a[nc*i+j]=aij;
}
void intMatrix::set(int i, int ai)
{
    a[i]=ai;
}
int intMatrix::size() const
{
    return nr*nc;
}
void intMatrix::zero()
{
    for (int i=0; i<size(); i++)
    {
        a[i]=0;
    }
}
//
