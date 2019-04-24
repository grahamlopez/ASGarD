#include "apply_fmwt.hpp"

#include "connectivity.hpp"
#include "matlab_utilities.hpp"
#include "quadrature.hpp"
#include "tensors.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>

template<typename P>
fk::matrix<P> apply_fmwt(fk::matrix<P> fmwt,fk::matrix<P> matrix2, int const kdeg, int const lev, int const isLeft, int const isTrans)
{
  int const n = kdeg*pow(2,lev);
  int ip = 0;
  int ipend = 2*kdeg-1;
  int col1 = 0;
  int col2 = n-1;
  //std::cout << "n ip ip end cols " << n << " " << ip << " " << ipend << " " << col1 << " " << col2 << std::endl;
  fk::matrix<P> product(n,n);
  if(isLeft)
  {
    if(isTrans)
    {
      fk::matrix<P> fmwt_sub1 = fmwt.extract_submatrix(ip,col1,ipend-ip+1,col2-col1+1);
      fk::matrix<P> fmwt_sub1t = fk::matrix<P>(fmwt_sub1).transpose();
      //fk::matrix<P> fmwt_sub11 = fmwt_sub1.extract_submatrix(0,0,4,4);
      //fk::matrix<P> fmwt_sub1t1 = fmwt_sub1t.extract_submatrix(0,0,4,4);
      //fmwt_sub11.print("sub NO transp");
      //fmwt_sub1t1.print("sub transp");
      //fk::matrix<P> mat2sub1 = matrix2.extract_submatrix(ip,0,ipend-ip+1,n);
      //fk::matrix<P> mat2sub11 = mat2sub1.extract_submatrix(0,0,4,4);
      //mat2sub11.print("sub mat2");
      fk::matrix<P> check1 = fmwt_sub1t*matrix2.extract_submatrix(ip,0,ipend-ip+1,n);
      product.set_submatrix(col1,0,check1);
    }
    else
    {
      fk::matrix<P> check1 = fmwt.extract_submatrix(ip,col1,ipend-ip+1,col2-col1+1)*matrix2.extract_submatrix(col1,0,col2-col1+1,n);
      product.set_submatrix(ip,0,check1);
    }
  }  
  else
  {
    if(isTrans)
    {
      fk::matrix<P> fmwt_sub1 = fmwt.extract_submatrix(ip,col1,ipend-ip+1,col2-col1+1);
      fk::matrix<P> check1 = matrix2.extract_submatrix(0,col1,n,col2-col1+1)*fmwt_sub1.transpose();
      product.set_submatrix(ip,0,check1);
    }
    else
    {
      fk::matrix<P> check1 = matrix2.extract_submatrix(0,ip,n,ipend-ip+1)*fmwt.extract_submatrix(ip,col1,ipend-ip+1,col2-col1+1);
      product.set_submatrix(ip,0,check1);
    }
  }
  //check1.print("check1");
  //fk::matrix<P> product = fmwt*matrix2;

  ip = 2*kdeg;
  for(int iLev=1;iLev<lev;iLev++)
  {
    int ncells = pow(2,iLev);
    int isize = n/ncells;

    for(int icell=0; icell< ncells;icell++)
    {
      //std::cout << "iLev icell " << iLev << " " << icell  << std::endl;
      ipend = ip + kdeg -1;
      col1 = icell*isize;
      col2 = col1 + isize;
      //std::cout << "n ip ip end cols " << n << " " << ip << " " << ipend << " " << col1 << " " << col2 << std::endl;
      if(isLeft)
      {
        if(isTrans)
        {
          fk::matrix<P> fmwt_sub1  = fmwt.extract_submatrix(ip,col1,ipend-ip+1,col2-col1);
          fk::matrix<P> fmwt_sub1t = fk::matrix<P>(fmwt_sub1).transpose();
	  product.set_submatrix(col1,0,product.extract_submatrix(col1,0,col2-col1,n) +
                                 fmwt_sub1t * 
          				 matrix2.extract_submatrix(ip,0,ipend-ip+1,n));
	}
	else
	{
          product.set_submatrix(ip,0,product.extract_submatrix(ip,0,ipend-ip+1,n) +
                                 fmwt.extract_submatrix(ip,col1,ipend-ip+1,col2-col1) * 
				 matrix2.extract_submatrix(col1,0,col2-col1,n));
        }
      }	
      else
      {
        if(isTrans)
	{
          fk::matrix<P> fmwt_sub1  = fmwt.extract_submatrix(ip,col1,ipend-ip+1,col2-col1);
	  product.set_submatrix(0,ip,product.extract_submatrix(0,ip,n,ipend-ip+1) + 
				 matrix2.extract_submatrix(0,col1,n,col2-col1) * fmwt_sub1.transpose());
	}
	else
	{
          product.set_submatrix(0,col1,product.extract_submatrix(0,col1,n,col2-col1) +
	        		 matrix2.extract_submatrix(0,ip,n,ipend-ip+1) *
                                 fmwt.extract_submatrix(ip,col1,ipend-ip+1,col2-col1));
	}
      }

      ip = ipend+1;					     
    }
  }
  return product;
}

template fk::matrix<double>
apply_fmwt(fk::matrix<double> fmwt,fk::matrix<double> matrix2, int const kdeg, int const lev, int const isLeft, int const isTrans);
template fk::matrix<float>
apply_fmwt(fk::matrix<float> fmwt,fk::matrix<float> matrix2, int const kdeg, int const lev, int const isLeft, int const isTrans);

