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
fk::matrix<P> reshape(fk::matrix<P> mat, int const dim1, int const dim2)
{
	  fk::vector<P> X(mat);
	  fk::matrix<P> Xreshape(dim1,dim2);
	  int count = 0;
	  for(int i=0;i<dim2;i++)
	  {
	    for(int j=0;j<dim1;j++)
	    {
	      count = i*dim1+j;
	      Xreshape(j,i) = X(count);
	    }
	  }
  return Xreshape;	  
}

template fk::matrix<double>
reshape(fk::matrix<double> mat,int const dim1, int const dim2);
template fk::matrix<float>
reshape(fk::matrix<float> mat,int const dim1, int const dim2);

template<typename P>
fk::matrix<P> apply_fmwt(fk::matrix<P> fmwt,fk::matrix<P> matrix2, int const kdeg, int const lev, int const isLeft, int const isTrans,int const method)
{
  int const n = kdeg*pow(2,lev);
  fk::matrix<P> product(n,n);
  if(method == 2)
  {
    int ip = 0;
    int ipend = 2*kdeg-1;
    int col1 = 0;
    int col2 = n-1;
    //std::cout << "n ip ip end cols " << n << " " << ip << " " << ipend << " " << col1 << " " << col2 << std::endl;
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
  }
  else if(method == 3)
  {
    int ip = 0;
    int ipend = 2*kdeg-1;
    int col1 = 0;
    int col2 = n-1;
    //std::cout << "n ip ip end cols " << n << " " << ip << " " << ipend << " " << col1 << " " << col2 << std::endl;
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
        check1.print("check1");
      }
      else
      {
        fk::matrix<P> check11 = fmwt.extract_submatrix(ip,col1,ipend-ip+1,col2-col1+1);
        fk::matrix<P> check12 = matrix2.extract_submatrix(col1,0,col2-col1+1,n);
        fk::matrix<P> check1 = fmwt.extract_submatrix(ip,col1,ipend-ip+1,col2-col1+1)*matrix2.extract_submatrix(col1,0,col2-col1+1,n);
        check11.print("check11");
        check12.print("check12");
        check1.print("check1");
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
        check1.print("check1");
      }
      else
      {
        fk::matrix<P> check1 = matrix2.extract_submatrix(0,ip,n,ipend-ip+1)*fmwt.extract_submatrix(ip,col1,ipend-ip+1,col2-col1+1);
        product.set_submatrix(0,col1,check1);
        check1.print("check1");
      }
    }
    //fk::matrix<P> product = fmwt*matrix2;

    ip = 2*kdeg;
    std::cout << "ip lev " << ip << " " << lev << std::endl;
    for(int iLev=1;iLev<=lev-1;iLev++)
    {
      std::cout << "ip iLev " << ip << " " << iLev << std::endl;
      int ncells = pow(2,iLev);
      int isize = n/ncells;
      
      int icell = 1;
      int ipend = ip + kdeg-1;
      col1 = 1 + (icell-1)*isize;
      col2 = col1 + isize-1;

      std::cout << "col1 and 2 " << col1 << " " << col2 << std::endl;
      fk::matrix<P> Fmat = fmwt.extract_submatrix(ip,col1-1,ipend-ip+1,col2-col1+1);
	  Fmat.print("Fmat mat");
      int ip2 = ip + (ncells * isize)-1;
      int ncolX = n;
      int nrowX = n;
      int ncol = (n*ncolX/isize);
      int nrows = (ncells*kdeg);
      ipend = ip + (ncells*kdeg)-1;
      std::cout << "ip and ipend " << ip << " " << ipend << std::endl;
    
      if(isLeft)
      {
        if(isTrans)
        {
	  matrix2.print("matrix2 mat");
	  fk::matrix<P> Xsub = matrix2.extract_submatrix(ip,0,ipend-ip+1,ncolX);
	  Xsub.print("Xsub mat");
          std::cout << "reshape dims " << kdeg << " " << ((ipend-ip+1)*ncolX/kdeg) << std::endl;
	  fk::matrix<P> XsubReshape = reshape(Xsub,kdeg,((ipend-ip+1)*ncolX/kdeg));
	  XsubReshape.print("XsubReshape mat");
          fk::matrix<P> FmatSub = Fmat.extract_submatrix(0,0,kdeg,isize);
          fk::matrix<P> FmatSub_transpose = fk::matrix<P>(FmatSub).transpose();
	  FmatSub.print("FmatSub mat");
          fk::matrix<P> FmatSubX = FmatSub_transpose*XsubReshape;
	  fk::matrix<P> FmatSubXreshape = reshape(FmatSubX,nrowX,ncolX);
	  FmatSubXreshape.print("FmatSubXreshape mat");
	  product.print("product mat");
	  product = product + FmatSubXreshape;
	  //fk::matrix<P> product2 = product + FmatSubXreshape;
	  //product.set_submatrix(0,0,product2);
	  //product.print("product mat");
        }
	else
	{
	  fk::vector<P> X(matrix2);
	  fk::matrix<P> Xreshape(isize,ncol);
	  int count = 0;
	  for(int i=0;i<ncol;i++)
	  {
	    for(int j=0;j<isize;j++)
	    {
	      count = i*isize+j;
	      Xreshape(j,i) = X(count);
	    }
	  }
	  X.print("X vector");
	  Xreshape.print("Xreshape matrix");
          fk::matrix<P> FmatSub = Fmat.extract_submatrix(0,0,kdeg,isize);
	  FmatSub.print("FmatSub mat");
          fk::matrix<P> FmatSubX = FmatSub*Xreshape;
	  FmatSubX.print("FmatSubX mat");
          fk::vector<P> FSXvec(FmatSubX);
	  std::cout << " fsxvec set" << std::endl;
	  fk::matrix<P> FSXreshape(nrows,(kdeg*ncol)/nrows);
	  std::cout << " fsxvec set" << std::endl;
	  count = 0;
	  std::cout << " irange "<< (kdeg*ncol)/nrows << std::endl;
	  std::cout << " jrange "<< nrows << std::endl;
	  for(int i=0;i<(kdeg*ncol)/nrows;i++)
	  {
	    for(int j=0;j<nrows;j++)
	    {
	      count = i*nrows+j;
	      FSXreshape(j,i) = FSXvec(count);
	    }
	  }
	  FSXreshape.print("FSXreshape vector");
          product.set_submatrix(ip,0,FSXreshape);


	}
      }
      else
      {
        if(isTrans)
        {
	   for(int icell=1;icell<=ncells;icell++)
	   {
             int j1 = ip + (icell-1)*kdeg;
             int j2 = j1 + kdeg-1;

             int jx1 = (icell-1)*isize;
             int jx2 = jx1 + isize-1;
	     std::cout << "j1 j2 jx1 jx2 " << j1 << " " << j2 << " " << jx1 << " " << jx2 << " " << std::endl;
             fk::matrix<P> FmatSubR = Fmat.extract_submatrix(0,0,kdeg,isize);
             fk::matrix<P> XsubR = matrix2.extract_submatrix(0,jx1,nrowX,jx2-jx1+1);
             fk::matrix<P> FmatX = XsubR*FmatSubR.transpose();
	     FmatX.print("FmatX");
	     product.set_submatrix(0,j1,FmatX);
             //Y(1:nrowX, j1:j2) = Y(1:nrowX, j1:j2) + ...
                  //X(1:nrowX,  jx1:jx2 ) * Fmat( 1:kdeg, 1:isize );
	   }	  
	}
	else
	{
	   for(int icell=1;icell<=ncells;icell++)
	   {
             int j1 = (icell-1)*isize;
             int j2 = j1 + isize-1;

             int jx1 = ip + (icell-1)*kdeg;
             int jx2 = jx1 + kdeg-1;
	     std::cout << "j1 j2 jx1 jx2 " << j1 << " " << j2 << " " << jx1 << " " << jx2 << " " << std::endl;
             fk::matrix<P> FmatSubR = Fmat.extract_submatrix(0,0,kdeg,isize);
             fk::matrix<P> XsubR = matrix2.extract_submatrix(0,jx1,nrowX,jx2-jx1+1);
             fk::matrix<P> FmatX = XsubR*FmatSubR;
	     FmatX.print("FmatX");
             fk::matrix<P> productPlus = product.extract_submatrix(0,j1,nrowX,j2-j1+1) + FmatX;
	     productPlus.print("productPlus");
	     product.set_submatrix(0,j1,productPlus);
             //Y(1:nrowX, j1:j2) = Y(1:nrowX, j1:j2) + ...
                  //X(1:nrowX,  jx1:jx2 ) * Fmat( 1:kdeg, 1:isize );
	   }	  
	}
      }
      
      ip = ip + (ncells * kdeg);

    }
  }  
  return product;
}

template fk::matrix<double>
apply_fmwt(fk::matrix<double> fmwt,fk::matrix<double> matrix2, int const kdeg, int const lev, int const isLeft, int const isTrans, int const method);
template fk::matrix<float>
apply_fmwt(fk::matrix<float> fmwt,fk::matrix<float> matrix2, int const kdeg, int const lev, int const isLeft, int const isTrans, int const method);

