#include<iostream>
#include<vector>
#include<complex>
#include<functional>
#include<tuple>
#include<experimental/optional>
#include<cmath>
#include "vector_oper.hpp"
#include<boost/numeric/odeint.hpp>
#include<numeric>
#include<type_traits>
#include<fstream>

// This is a dimension-th order Rungekutta class for solving nonlinear systems of differential equations

// closed in namespace Math
namespace Math
{
    namespace Algorithm
    {

   // building class initial condition
  /* template<class T,typename std::enable_if<std::is_default_constructible<T>::value>::type>
   class Initial_cond
   {
       protected:
       T y_0;

       public:
       Initial_cond(y0): y_0{y0}{}
   };*/

   // begin of template class Runge_Kutta 

   template<class Order,class Type,class Vec_Type,class Function>//,typename std::enable_if<std::is_floating_point<Order>::value>::type>  // Order=nth order of method, Type=Type of function
   class Runge_Kutta
   {
       private:
     // attributes of a class specified in a c++ way
     
     // nuber of Rungekutta order
     std::vector<Vec_Type> Initial_cond; // vector of equations going in
     Type x_0; // initial cond
     Type x_max; // maximal condition
     Type step_size; // step size
     std::size_t dimension; //number of dimensions (e.g 2 )
     // function 
      std::function<std::vector<Vec_Type>(Type,std::vector<Vec_Type>)> F;
    // private constructor in class 
    std::size_t N; // number of steps

     
     public:
     // constructor for class
     //Runge_Kutta(std::size_t n_,std::vector<Type> y_,Initial_cond<std::vector<Type>> I_(std::vector<Type> y0),Initial_cond<std::vector<Type>> B_(std::vector<Type> ymax),size_t n_ ) : N(N_),Initial_cond(y_),I(I_),B(B_),dimension(n_){};
     // a Solve method for rungekutta class


     // public constructor with error checks
     Runge_Kutta( std::size_t dimension_,Type x0,Type x_max_,Type step_size_,std::vector<Vec_Type> Initial_cond_,std::function<std::vector<Vec_Type>(Type,std::vector<Vec_Type>)>f) : dimension(dimension_),x_0(x0),step_size(step_size_),Initial_cond(Initial_cond_),F(f),x_max(x_max_)
     {
         //std::cout<<Initial_cond[0]<<" "<<Initial_cond_[0]<<std::endl;
         Initial_cond=Initial_cond_;

         //std::cout<<Initial_cond[0]<<" "<<Initial_cond_[0]<<std::endl;
        try
        {
         if(dimension_>=0 && step_size_>=0 &&Initial_cond_.size()==dimension_)
         {
            
         }
         else
         {
            throw std::runtime_error("Something wrong with constructor");
         }
         
        }
        catch(std::bad_alloc&e)
        {
            std::cerr<<" failure in constrcutor"<<std::endl;
        }
         
     }

     // Final component important for studies for resonances or boundary studies:
     std::vector<Vec_Type> getfinalYComponent();

     std::vector<std::vector<Vec_Type>> Step(Type x,std::vector<Vec_Type> Initial_cond,Type step_size);
     std::vector<Type> Step(Order N);
     std::vector<Vec_Type> Step(Type x,std::vector<Vec_Type> Initial_cond);
 
     std::vector<Type> Solve();
    void Runge_Kutta_4th(std::vector<std::vector<Vec_Type>>& Matrix,std::vector<Type>& x) noexcept;

     std::vector<Vec_Type>  Rungekutta_Fehl();

     int size()
     {
         std::vector<double> Initial_cond;
         Initial_cond=getfinalYComponent();
         return Initial_cond.size();
         
     }

     Vec_Type Index(std::size_t i,std::size_t j);



     
   };

   // implementation of solve function 
   /*template<typename Order,class Type>//,typename std::enable_if<std::is_floating_point<Order>::value>::type>
   std::vector<Type> Runge_Kutta<Order,Type>::Step(Order N)
   {
      
        // here calculation of partial sum matrix

        // coefficients b and c in Rungekutta formula

        std::vector<Type> b;
        std::vector<Type> c;
        std::vector<Type> epsilon;
        std::vector<Type> y_new;


        if(std::partial_sum(1,N,b)!=1)
        {
           throw std::runtime_error("Condition for Taylor expansion not fullfilled");
        }

        // calculate the coefficient matrix

        Matrix<Type,Type> a(std::ranges::size(c),std::ranges::size(b));

        for(auto i=std::ranges::begin(a),i<std::ranges::end(a);i++)
        {
            if(std::partial_sum(1,i-1,a)!=c[i])
            {
                throw std::runtime_error("Condition for Taylor expansion not fullfilled");
            }

        }

    
        // the k-factors considering the Runge_kutta

        std::vector<Type> k(std::ranges::size(k));

        for(auto i : k)
        {
          if(std::ranges::empty(k))
          {
             throw std::runtime_error("range is empty");
          }
          i=Initial_cond+step_size*std::partial_sum(1,N,a[i]*f(k,x+c[i]*step_size));
        }
         
  
          
      



        // coefficients b

        if(N==1)
        {
           b{1,2,3,4};
           a{{0,1},{0,0}};
        }

        else if(N==2)
        {

        }

        if(N==3)
        {

        }

        if(N==4)
        {

            b{1/6,1/3,1/3,1/6};
            c{0,1/2,1/2,1};
            a{{0,0,0,0},{1/2,0,0,0},{0,1/2,0,0},{1/6,1/3,1/3,1/6}};
           
        }

        if(N==5)
        {

        }

        if(N==6)
        {

        }

        else
        {
            throw std::runtime_error("Order Not defined");
        }

        return Initial_cond;
        
 
   }*/

   template<typename Order,class Type,class Vec_Type,class Function>
   std::vector<Vec_Type> Runge_Kutta<Order,Type,Vec_Type,Function>::Rungekutta_Fehl()
   {
       std::vector<Vec_Type> Ret(dimension);
      // this is an algorithm specified for the Runge-Kutta Fehlberg Method in order to get higher precision
       Type hc=step_size;
      // adaptive step-size
      std::vector<Vec_Type> x(N);
       x[0]=x_0;
      std::vector<std::vector<Vec_Type>> yv(2);
      // setting maximal alowed error to step_size;
      Type epsilon_max=std::pow(step_size,4);

     for(auto i=1;i<N;i++)
      {
          yv=Step(x[i],Initial_cond,hc);
           x=x+hc;
          Type epsilon;
          epsilon=Norm(yv[0]-yv[1]);
          // adapt if error is ok
          if(epsilon<epsilon_max)
          {
              hc=hc*0.84*std::pow(epsilon_max/epsilon,0.25);
          }
          // re do step if error is too large
          while(epsilon>epsilon_max)
          {  
              hc=hc*0.84*std::pow(epsilon_max/epsilon,0.25);
              yv=Step(x,Initial_cond,hc);
              epsilon=Norm(yv[0]-yv[1]);
          }

          // adapt step x;
         
        Ret=yv[1];
        
        
         
      }
      return Ret;
  }

   template<typename Order, class Type,class Vec_Type,class Function>
   std::vector<std::vector<Vec_Type>> Runge_Kutta<Order,Type,Vec_Type,Function>::Step(Type x,std::vector<Vec_Type> Initial_cond,Type step_size)
   {
     // define a vector to store the data in:
     std::vector<std::vector<Vec_Type>> Ret(2);
     std::vector<Vec_Type> y_4(dimension);
     std::vector<Vec_Type> y_5(dimension);
     x=x_0;
     std::vector<Vec_Type> y1(dimension),y2(dimension),y3(dimension),y4(dimension),y5(dimension);

     auto k1=F(x,Initial_cond);
     for(auto i=0;i<dimension;i++)
     {
        y1[i]=Initial_cond[i]+k1[i]/4.0;
     }

     auto k2=F(x+step_size/4.0,y1);
     for(auto i=0;i<dimension;i++)
     {
         y2[i]=Initial_cond[i]+(3.0*k1[i]/32.0)+(9.0*k2[i]/32.0);
     }
     
    auto k3=F(x+(3.0*step_size/8.0),y2);
    for(auto i=0;i<dimension;i++)
    {
        y3[i]=Initial_cond[i]+(1932.0*k1[i]/2197.0)-(7200.0*k2[i]/2197.0)+(7296.0*k3[i]/2197.0);
    }
    auto k4=F(x+(12.0*step_size/13.0),y3);
    for(auto i=0;i<dimension;i++)
    {
        y4[i]=Initial_cond[i]+(439.0*k1[i]/216.0)-(8.0)*k2[i]+(3680.0*k3[i]/513.0)-(845.0*k4[i]/4104.0);
    }

    auto k5=F(x+step_size,y4);
    for(auto i=0;i<dimension;i++)
    {
      y5[i]=Initial_cond[i]-(8.0/27.0)*k1[i]+(2*k2[i])-(3544.0*k3[i]/2565.0)+(1859.0*k4[i]/4104.0)-(11.0/40.0)*k5[i];
    }
    auto k6=F(x+step_size/2,y5);
     
     for(auto i=0;i<dimension;i++)
     {
         y_4[i]=Initial_cond[i]+step_size*(25*k1[i]/216.0)+(1408.0*k3[i]/2565.0)+(2197.0*k4[i]/4101.0)-(k5[i]/5.0);
         y_5[i]=Initial_cond[i]+step_size*(16*k1[i]/135.0)+(6656.0*k3[i]/12825.0)+(28561.0*k4[i]/56430.0)-(9.0*k5[i]/50.0)+(2.0*k6[i]/55.0);
     }
    

     //store y_4 and y_5

     Ret[0]=y_4;
     Ret[1]=y_5;

     return Ret;
   }


   template<class Order, class Type,class Vec_Type,class Function>
   std::vector<Vec_Type> Runge_Kutta<Order,Type,Vec_Type,Function>::Step(Type x,std::vector<Vec_Type> Initial_cond)
   {
         
      std::vector<Vec_Type> y1(dimension),y2(dimension),y3(dimension),y4(dimension);

      auto k1=F(x,Initial_cond);
     // std::cout<<k1[0]<<" "<<k1[1]<<std::endl;
    
      for(auto i=0;i<dimension;i++)
      {
          y1[i]=Initial_cond[i]+step_size*k1[i]/2.0;
      }
    //  std::cout<<y1[0]<<" "<<y1[1]<<std::endl;
      auto k2=F(x+step_size/2.0,y1);
     // std::cout<<k2[0]<<" "<<k2[1]<<std::endl;
      for(auto i=0;i<dimension;i++)
      {
          y2[i]=Initial_cond[i]+step_size*k2[i]/2.0;
      }
     // std::cout<<y2[0]<<" "<<y2[1]<<std::endl;
      auto k3=F(x+step_size/2.0,y2);
     //std::cout<<k3[0]<<" "<<k3[1]<<std::endl;
      for(auto i=0;i<dimension;i++)
      {
          y3[i]=Initial_cond[i]+step_size*k3[i];
      }
      //std::cout<<y3[0]<<" "<<y3[1]<<std::endl;
      

      auto k4=F(x+step_size,y3);
      //std::cout<<k4[0]<<" "<<k4[1]<<std::endl;
      for(auto i=0;i<dimension;i++)
      {
          Initial_cond[i]=Initial_cond[i]+step_size/6.0*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i]);
      }


      //std::cout<<Initial_cond[0]<<" "<<Initial_cond[1]<<std::endl;
      return Initial_cond;
   }

   template<class Order, class Type,class Vec_Type,class Function>
  void Runge_Kutta<Order,Type,Vec_Type,Function>::Runge_Kutta_4th(std::vector<std::vector<Vec_Type>>& Matrix,std::vector<Type>& x) noexcept
   {
       unsigned int N=(int)(x_max/step_size);

     x[0]=x_0;
     
      // this is an algorithm specified for the Runge-Kutta Fehlberg Method in order to get higher precision
       //std::cout<<hc<<std::endl;
       Matrix[0]=Initial_cond;
       //std::cout<<Initial_cond<<std::endl;
       //std::cout<<Matrix[0]<<std::endl;

      
      // adaptive step-size
    //iterate over number of Runge-Kutta steps
   for(auto i=0;i<N;i++)
    {
       
        
        Matrix[i+1]=Step(x[i],Matrix[i]);
        x[i+1]=x[i]+step_size;
        
    }
    

    
   }
 
   template<class Order, class Type,class Vec_Type,class Function>
  std::vector<Vec_Type> Runge_Kutta<Order,Type,Vec_Type,Function>::getfinalYComponent()
 {
      // return Vector as the last component of Ret
      std::vector<Type> x(N+1);
      std::vector<std::vector<Vec_Type>> y(N+1,std::vector<Vec_Type>(dimension));
      
      // look for errors
      Runge_Kutta_4th(y,x);
      
     
      
      return y[N];
 }

 template<class Order,class Type,class Vec_Type,class Function>
  Vec_Type  Runge_Kutta<Order,Type,Vec_Type,Function>::Index(std::size_t i,std::size_t j)
  {
       Order N=(x_max-x_0)/step_size;
      std::vector<Type> x(N+1);
      std::vector<std::vector<Vec_Type>> y(N+1,std::vector<Vec_Type>(dimension));
      
      // look for errors
      Runge_Kutta_4th(y,x);

     return y[i][j]; 

  }




























}

}
