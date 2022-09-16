// this is a c++ library for numerical integration

// include statements to make it work

#include<iostream>
#include<vector> // needed for std::vector
#include<functional> // function class
#include<numeric>
#include<fstream>
#include<filesystem>
#include<tuple>
#include<random>
// stored in namespace Math

namespace Math
{
    // first integration method is the Riemann integration
    template<class Type,typename N>
    class Riemann_int
    {
        protected: 
        // attributes are the end point, start point and step size:
        Type a;  // starting point 
        Type b;  // end-point 
        N n; // number of steps;
        std::function<Type(Type)> F;

       


        public:

        // constructor for the Riemann_int;

        Riemann_int(Type a_,Type b_,N n_,std::function<Type(Type)> f): a(a_),b(b_),n(n_),F(f){}

        // function calculating the integral
        Type Solve();


        // destructor

        ~Riemann_int()=default;



    };


    // implementation of function solve:
    template<class Type,typename N>
    Type Riemann_int<Type,N>::Solve()
    {
        // define a step size:
        auto h=(b-a)/n;
        Type x=a;
        Type ret;

        for(auto i=0;i<n;i++)
        {
            ret+=F(x);
            x+=h;
        }

        return h*ret;


    }

    // the next integration procedure is the trapezoidal rule
    template<class Type>
    class Trapeziodal
    {
           
      protected:
      
        Type a;  // starting point 
        Type b;  // end-point 
        std::size_t N; // step-size
        std::function<Type(Type)> F;

     public:

      // constructor for the Riemann_int;

       Trapeziodal(Type a_,Type b_,std::size_t N_,std::function<Type(Type)> f): a(a_),b(b_),N(N_),F(f){}

        // function calculating the integral
        Type Solve();


        // destructor

       // ~Trapezoidal()=default;


    };

    template<class Type>
    Type Trapeziodal<Type>::Solve()
    {
        Type ret;
        Type h=(b-a)/N;
        for(auto x=a;x<b;x+=h)
        {
          ret+=F(x)+F(x+h);
        }
        ret*=0.5*h;
        return ret;
    }

    // simpsons rule

    template<class Type>
    class Simpson
    {
           
      protected:
      
        Type a;  // starting point 
        Type b;  // end-point 
        std::size_t N; // step-size
        std::function<Type(Type)> F;

      public:

      // constructor for the Riemann_int;

       Simpson(Type a_,Type b_,std::size_t N_,std::function<Type(Type)> f): a(a_),b(b_),N(N_),F(f){}

        // function calculating the integral
        Type Solve();


        // destructor

       // ~Trapezoidal()=default;


    };

    template<class Type>
    Type Simpson<Type>::Solve()
    {
        Type ret;
        Type h=(b-a)/N;
        for(auto x=a;x<b;x+=h)
        {
          ret+=F(x)+4*F(x+h/2.0)+F(x+h);
        }
        ret*=h/6.0;
        return ret;
    }




    template<class T,class U,std::size_t N>
    class MC_1D_integration
    {
         
          protected:

          
          T error;

          T average;
        
          std::size_t num_sweeps;

          T x_lower;

          T x_upper;

          std::function<U(T)> f;

          public:

          MC_1D_integration(T x_low, T x_up, std::size_t sweep,std::function<U(T)> F): x_lower(x_low), x_upper(x_up),f(F), num_sweeps(sweep) {}


          auto getlowerboundary()
          {
            return x_lower;
          }


          auto getupperBoundary()
          {
             return x_upper;
          }


          auto getError()
          {
            return error;
          }

          auto getAverage()
          {
            return average;
          }

          std::pair<T,T> MC_Integration();

          void printmessage();

          void writefunctioninfile(int number_of_points);







    };

    template<class T,class U,std::size_t N>
    void MC_1D_integration<T,U,N>::writefunctioninfile(int number_of_points)
    {
      //std::filesystem::create_directory("Results");

      //std::filesystem::current_path("Results");

      std::ofstream file;
      file.open("function.asc");

      T h=(x_upper-x_lower)/number_of_points;

      file<<"x:"<<" "<<"f(x)"<<std::endl;
      file<<std::endl;
      for(auto x=x_lower;x<=x_upper;x+=h)
      {
       
        file<<x<<"    "<<f(x)<<std::endl;
      }
    }

    template<class T, class U, std::size_t N>
    void MC_1D_integration<T,U,N>::printmessage()
    {
         writefunctioninfile(1000);
         char option;


        std::cout<<"This is a c++ class which will calculate the 1 dimensional integral of a function based on the Monte Carlo Method"<<std::endl;

        std::cout<<"For more information enter one of the following options:"<<std::endl;
        std::cout<<"basic algorithm description [B] "<<std::endl;
        std::cout<<"description on usage of class [D]"<<std::endl;

        std::cin>>option;
        int input;
      

        switch(option)
        {

          case 'B': 
          std::cout<<" The Monte Carlo Method was initially developed for the generation of pseudo random numbers in the Computer"<<std::endl;

          case 'D':

          std::cout<<"Quantum theory is cool"<<std::endl;

          input++;


        }

    }

    template<class T,class U, std::size_t N>
    std::pair<T,T> MC_1D_integration<T,U,N>::MC_Integration()
    {

       //printmessage();
         std::random_device rnd_device;
         std::mt19937 generator(rnd_device());
         std::uniform_real_distribution<double> distribution(x_lower,x_upper);
      
        T variance;
        // define array with points

        std::vector<T> x(num_sweeps-1);

        for(auto& elem : x)
        {
          elem=distribution(generator);
          //std::cout<<elem<<std::endl;
          average+=f(elem);
          variance+=1.0/num_sweeps*f(elem)*f(elem);
        }

        average*=1.0/num_sweeps*(x_upper-x_lower);
        error=(x_upper-x_lower)*std::sqrt((variance-std::pow(average,2.0))/num_sweeps);

        //std::cout<<" actual result:"<<M_PI/4.0<<" "<<" average result: "<<average<<" "<<"Error: "<<error<<std::endl;
      




      return std::make_pair(average,error);



    }


    template<class T,class U,std::size_t N>
    class MC_3D_integration
    {
         
          protected:

          
          T error;

          T average;
        
          std::size_t num_sweeps;

          std::tuple<T,T,T> lower_limit;

          std::tuple<T,T,T> upper_limit;

          std::function<U(T,T,T)> f;

          public:

          MC_3D_integration(std::tuple<T,T,T> low, std::tuple<T,T,T> up, std::size_t sweep,std::function<U(T,T,T)> F): lower_limit(low), upper_limit(up),f(F), num_sweeps(sweep) {}


          auto getlowerboundary()
          {
            return lower_limit;
          }


          auto getupperBoundary()
          {
             return upper_limit;
          }


          auto getError()
          {
            return error;
          }

          auto getAverage()
          {
            return average;
          }

          std::pair<T,T> MC_Integration();

          void printmessage();

          void writefunctioninfile(int number_of_points);







    };

    template<class T,class U,std::size_t N>
    void MC_3D_integration<T,U,N>::writefunctioninfile(int number_of_points)
    {
      /*std::filesystem::create_directory("Results");

      std::filesystem::current_path("Results");

      std::ofstream file;
      file.open("function.asc");

      T h=(x_upper-x_lower)/number_of_points;

      file<<"x:"<<" "<<"f(x)"<<std::endl;
      file<<std::endl;
      for(auto x=x_lower;x<=x_upper;x+=h)
      {
       
        file<<x<<"    "<<f(x,0)<<std::endl;
      }*/
    }

    template<class T, class U, std::size_t N>
    void MC_3D_integration<T,U,N>::printmessage()
    {
         writefunctioninfile(1000);
         char option;


        std::cout<<"This is a c++ class which will calculate the 1 dimensional integral of a function based on the Monte Carlo Method"<<std::endl;

        std::cout<<"For more information enter one of the following options:"<<std::endl;
        std::cout<<"basic algorithm description [B] "<<std::endl;
        std::cout<<"description on usage of class [D]"<<std::endl;

        std::cin>>option;
        int input;
      

        switch(option)
        {

          case 'B': 
          std::cout<<" The Monte Carlo Method was initially developed for the generation of pseudo random numbers in the Computer"<<std::endl;

          case 'D':

          std::cout<<"Quantum theory is cool"<<std::endl;

          input++;


        }

    }

    template<class T,class U, std::size_t N>
    std::pair<T,T> MC_3D_integration<T,U,N>::MC_Integration()
    {

       //printmessage();
         T x_lower=std::get<0>(lower_limit);
         T y_lower=std::get<1>(lower_limit);
         T z_lower=std::get<2>(lower_limit);
         T x_upper=std::get<0>(upper_limit);
         T y_upper=std::get<1>(upper_limit);
         T z_upper=std::get<2>(upper_limit);

         std::random_device rnd_device;
         std::mt19937 generator(rnd_device());
         std::uniform_real_distribution<double> distribution_1(x_lower,x_upper);
         std::uniform_real_distribution<double> distribution_2(y_lower,y_upper);
         std::uniform_real_distribution<double> distribution_3(z_lower,z_upper);
      
        T variance;
        // define array with points

        std::vector<T> x(num_sweeps-1);
        std::vector<T> y(num_sweeps-1);
        std::vector<T> z(num_sweeps-1);

        for(auto i=0;i<num_sweeps-1;i++)
        {
          x[i]=distribution_1(generator);
          y[i]=distribution_2(generator);
          z[i]=distribution_3(generator);
          //std::cout<<elem<<std::endl;
          average+=f(x[i],y[i],z[i]);
          variance+=f(x[i],y[i],z[i])*f(x[i],y[i],z[i]);
        }

        average*=1.0/num_sweeps*(x_upper-x_lower)*(y_upper-y_lower)*(z_upper-z_lower);
        error=(x_upper-x_lower)*(y_upper-y_lower)*(z_upper-z_lower)*std::sqrt((1.0/num_sweeps*variance-std::pow(average,2.0))/num_sweeps);
        //error=variance-std::pow(average,2.0);
        //std::cout<<" actual result:"<<M_PI/4.0<<" "<<" average result: "<<average<<" "<<"Error: "<<error<<std::endl;
      




      return std::make_pair(average,error);



    }
    



}