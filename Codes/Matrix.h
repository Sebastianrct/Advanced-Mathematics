
// this is a c++ matrix class for further perpuses

// this is the class matrix


#include<vector>
#include<complex>
#include<iostream>
#include<iomanip>
#include<initializer_list>
#include<algorithm>

#pragma once

template<typename T>
class Matrix
{
 // the private attributes in the class given by the Matrix properties
 //using value_type=T;
 private:
 size_t mrows;
 size_t mcolumns;
 std::vector<std::vector<T>> data; 
// The public attributes of the class

 public:
T& operator[](size_t i)
 {
  return data[i];
 }
 // start with the big five construtors
 
 // default constructor

 Matrix()=default;

 Matrix(size_t rows,size_t columns) : data(rows),mcolumns(columns),mrows(rows)
 {
   for(auto i=0;i<mrows;i++)
    data[i].resize(mcolumns);
 }
 

 
 /*T operator()(size_t i,size_t j) const
 {
  return data[i][j];
 }*/
 std::vector<T>& operator[](int r)
 {
   return data[r];
 }

 T& operator()(int r,int c)
 {
  return data[r][c];
 }

 // resize a Matrix
 void resize(size_t Rows,size_t Columns);

// initialize Matrix with list
Matrix( std::initializer_list<std::initializer_list<T>> il) : mrows(il.size()), mcolumns(  il.size() ? il.begin()->size():0),data(il.begin(),il.end())
 {
   for(auto i=0;i<mrows;i++)
   {
    std::vector<T> v(mcolumns);
    for(auto j : v)
    {
     data[i].push_back(j);
    }
   }
 }
// Matri

 size_t cols()const
 {
  return mcolumns;
 }
 size_t row() const
 {
  return mrows;
 }

 // compare two matrices

/* bool operator==(const Matrix& m)
 {
  return is_Matrix;
 }*/

 // copy and move constructor

 Matrix(Matrix && A) : mrows(A.mrows) ,mcolumns(A.mcolumns),data(std::move(A.data))
 {
 }

 Matrix<T>& operator=(Matrix<T>&& A)
 {
  mcolumns=A.mcolumns;
  mrows=A.mrows;
  data=std::move(A.data);
  return *this;
 }

 // copy constructor
 Matrix(const Matrix& A) : mrows(A.mrows),mcolumns(A.mcolumns)
 {
  data=std::copy(A.data);
 }

 Matrix& operator=(const Matrix<T>&A)
 {
   mrows=A.mrows;
   mcolumns=A.mcolumns;
   data=std::move(A.data);
   return *this;

 }

  
 //Matrix operator( Matrix&)=default;
 // destructor of class

 ~Matrix()=default;
  
 // some useful operations on data
 //template<typename U>
 //Matrix& operator*=(const std::vector<T>);
 // assignment operator
 /*Matrix<T>& operator=(const Matrix<T>& A)
 {
   
 }*/
 Matrix<T>& operator*=(double x);
 std::vector<T> operator*(const std::vector<T> &v)
 {
   std::vector<T> result(v.size());
 for(auto i=0;i<mrows;i++)
 {
  for(auto j=0;j<mcolumns;j++)
  {
    result[i]=this->data[i][j]*v[j];
  }
 }

 return result;
 }
 // dot product for matrices
 Matrix<T> operator*(const Matrix<T>& A);
 Matrix<T> operator/(T x);
 // *= operator for matrices
 Matrix<T>& operator*=(const Matrix& A);
 Matrix<T>& operator+=(double x);
// Matrix& operator*=(const std::vector<T>);
 const Matrix<T> operator+(  Matrix<T>& r)
 {
   Matrix<T> result;
   for(auto i=0;i<mrows;i++)
   {
     for(auto j=0;j<mcolumns;j++)
     {
       result(i,j)=this->data[i][j]+r(i,j);
     }
   }
   return result;
 }
 Matrix<T> operator-(const Matrix<T>& r);
 Matrix& operator-=(double scalar);
 Matrix<T>& operator-=(const std::vector<T>);
 Matrix<T>& operator-=(const Matrix<T>& r);
 // the input and output operator
template<typename U>
friend std::ostream& operator<<(std::ostream &os,const Matrix<U>& A);

friend std::istream& operator>>(std::istream &is, Matrix<T> A)
{
 int rows,columns;
 is>>rows>>columns;
 if(!is)
 {
  throw "error in istream no values were inputed";
  return is;
 }
  Matrix tmp(rows,columns);
  std::cout<<"Enter values"<<std::endl;
 for(auto i=0;i<A.row();i++)
  {
   for(auto j=0;j<A.cols();j++)
   {
    is>>tmp.data[i][j];
   }
  }
  A=tmp;
  return is;
}
 
 
 size_t index(size_t r,size_t c)
 {
  return r*cols()+c;
 }
 //std::istream &operator>>(std::istream &is,const Matrix<T> &A);

 // iterators for matrix class written in own class
 class iterator
 {
  friend class Matrix;

  public:
   iterator rowIterBegin(size_t mrows)
   {
    return data[mrows].begin(); 
   }
   auto rowIterEnd(size_t mrows) {return data[mrows].end();}
   
   iterator colIterBegin(size_t mcolumns)
   {
    return data[mcolumns].begin();
   }
  
   iterator colIterEnd(size_t mcolumns) {return data[mcolumns].end();}
   
 };
 
  
 // Matrix product;
 //Matrix operator*(Matrix& A,Matrix& B);

 // another useful functions

 // the transpose,hermitian conjugated,inverse,...

  //void inverse(Matrix<T> & A);
  Matrix transpose()
  {
    Matrix result;
    for(auto i=0;i<mrows;i++)
    {
      for(auto j=0;j<mcolumns;j++)
      {
         result(i,j)=this->data[j][i];
      }
    }
    return result;
  }

  Matrix<std::complex<T>> conjugate()
  {
    Matrix result=this->transpose();
    for(auto i=0;i<mrows;i++)
    {
      for(auto j=0;j<mcolumns;j++)
      {
        std::imag(result(i,j))=-std::imag(result(i,j));
      }
    }
    return result;
  }
  T Trace()
  {
    if(mrows!=mcolumns)
    {
      throw std::runtime_error("Matrix is not quadratic");
    }
    T result;
    for(auto i=0;i<mrows;i++)
    {
      for(auto j=0;j<mcolumns;j++)
      {
        if(i==j)
        {
        result+=this->data[i][i];
        }
      }
    }
    return result;
  }

  // get submatrix in c++
  Matrix submatrix(size_t nrows,size_t ncolumns);
  // ordinary determinant
  T det()
  {
    if(mrows!=mcolumns)
 {
  throw std::invalid_argument("inhomogeneous matrix");
 }
 size_t Dim=mrows;
 
 T result;
 Matrix<T> TempMatrix(Dim,Dim); 
 // default argument with matrix size 2x2

  if(data.size()==1)
  {
   return data[0][0];
  }

  else if(data.size()==2)
  {
   result=this->data[0][0]*data[1][1]-data[0][1]*data[1][0];
  } 

  /*else
  {
   for(auto p=0;p<Dim;p++)
   {
    for(auto i=1;i<Dim;i++)
    {
     auto j1=0;
     for(auto j=0;j<Dim;j++)
     {
      if(j!=p)
      {
       TempMatrix[p]=this->data[i][j];
      }
     }
    }
    result-=(this->data[0][p])*TempMatrix[p].det();
   
     return result;
   }
  }*/
 return result; 
  }
  // inverse of A
  // cofactor

  Matrix<T> Cofactor()
  {
    // result

 Matrix<T> result(mrows,mcolumns);
 // checking whether Matrix was quadratic

 if(mrows!=mcolumns)
 {
  throw std::runtime_error("Matrix is not quadratic");
 }



 if(mrows==2)
 {
  result(0,0)=this->data[1][1];
  result(0,1)=-this->data[0][1];
  result(1,0)=-this->data[1][0];
  result(1,1)=this->data[0][0];
  
 }

 else if(mrows>2)
 {
  throw std::runtime_error("Matrix is not 2x2");
 }

 return result;
  }  
  
  // Inverse
  // kronecker product
  Matrix<T> inverse()
  {
    Matrix<T> inv(mrows,mcolumns);

 // define det
  Matrix<T> cofactor;

  T determinant=det();
  cofactor=this->Cofactor();
  // components of the adjoint matrix

  if(determinant==0.0)
  {
   throw std::runtime_error ("Matrix has poles");
  }
  // set inverse
if(data.size()==1)
{
   inv(0,0)=1.0/data[0][0];
}

  for(auto i=0;i<mrows;i++)
  {
   for(auto j=0;j<mcolumns;j++)
   {
    inv(i,j)=cofactor(i,j)/determinant;
   }
  }

  return inv;
   
  }
  Matrix tensor(const Matrix<T>& B);
  Matrix Solve();
  Matrix LU_dec();
  

};

 // the corresonding output operator
 template<typename U>
std::ostream& operator<<(std::ostream& os, Matrix<U>& A)
 {
  for(auto i=0;i<A.row();i++)
  {
   for(auto j=0;j<A.cols();j++)
   {
    os<<A(i,j)<<" ";
   }
   os<<std::endl;
  }
  return os;
 }



// the input operator
