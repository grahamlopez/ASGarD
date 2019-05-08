#include "matlab_utilities.hpp"
#include "pde.hpp"
#include "tests_general.hpp"
#include "apply_fmwt.hpp"
#include "transformations.hpp"
#include <numeric>

TEMPLATE_TEST_CASE("apply_fmwt", "[apply_fmwt]", double, float)
{
  auto const relaxed_comparison = [](auto const first, auto const second) {
    auto first_it = first.begin();
    std::for_each(second.begin(), second.end(), [&first_it](auto &second_elem) {
      auto const f1 = *first_it++; 
      std::cout.precision(17);
      std::cout << "first and second element " << std::fixed << f1
                   << " " << std::fixed << second_elem << std::endl;
      std::cout << "diff " << std::fixed << std::abs(f1 - second_elem) << std::endl;
      std::cout << "epsilon  "<< std::numeric_limits<TestType>::epsilon() << std::endl;
      std::cout << "max  "<< std::max(std::abs(f1),std::abs(second_elem)) << std::endl;
      std::cout << "compare " << std::numeric_limits<TestType>::epsilon()*1.0e3*std::max(std::abs(f1),std::abs(second_elem)) << std::endl;
      auto const diff = std::abs(f1 - second_elem);///std::max(std::abs(f1),std::abs(second_elem));
      auto const abs_diff = std::abs(f1 - second_elem);
      REQUIRE(diff < std::numeric_limits<TestType>::epsilon()*1.0e3*std::max(std::abs(f1),std::abs(second_elem)));
      //REQUIRE(Approx(f1)
      //            .margin(std::numeric_limits<TestType>::epsilon()) ==
      //        second_elem);
      //REQUIRE(std::abs((*first_it++) - second_elem) < 1.0e-5*std::max(std::abs((*first_it++)),std::abs((second_elem))));
    });
  };
  
  SECTION("Apply fmwt test 1")
  {
    std::cout << " type epsilon " << std::numeric_limits<TestType>::epsilon() << std::endl;
    int const kdeg = 2;
    int const lev = 2;
    std::string out_base =
        "../testing/generated-inputs/apply_fmwt/";

    std::string mat1_string    = out_base + "mat1_k1_lev2.dat";
    //std::string fmwt_string    = out_base + "fmwt_k1_lev2.dat";
    
    dimension const dim =
        make_PDE<TestType>(PDE_opts::continuity_1, lev, kdeg)
            ->get_dimensions()[0];

    fk::matrix<TestType> const fmwt = operator_two_scale<TestType>(dim);

    fk::matrix<TestType> mat1 =
        fk::matrix<TestType>(read_matrix_from_txt_file(mat1_string));
    //fk::matrix<TestType> fmwt =
    //    fk::matrix<TestType>(read_matrix_from_txt_file(fmwt_string));
    
    //mat1.print("mat1");	
    //fmwt.print("fmwt");	
    int isLeft = 1;
    int isTrans = 0;
    int method = 2;
    fk::matrix<TestType> productLeftFull = fmwt*mat1;
    auto const productLeft = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans,method);
    productLeftFull.print("full");	
    productLeft.print("apply left");	
    
    SECTION("degree = 4, lev 5 fmwt apply left") { relaxed_comparison(productLeftFull, productLeft); }
  
  //
    isLeft = 1;
    isTrans = 1;
    fk::matrix<TestType> fmwt_transpose = fk::matrix<TestType>(fmwt).transpose();
    fk::matrix<TestType> productLeftTransFull = fmwt_transpose*mat1;
    auto const productLeftTrans = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans,method);
    //productLeftFull.print("productLeftTransFull");	
    //productLeft.print("productLeftTrans");	
   
    
    //fk::matrix<TestType> sub1 = productLeftTransFull.extract_submatrix(0,0,4,4);
    //fk::matrix<TestType> sub2 = productLeftTrans.extract_submatrix(0,0,4,4);
    //sub1.print("sub1");
    //sub2.print("sub2");
    SECTION("degree = 4, lev 5 fmwt apply left transpose") { relaxed_comparison(productLeftTrans, productLeftTransFull); }
    //SECTION("degree = 4, lev 5 fmwt apply left transpose") { relaxed_comparison(sub1,sub2); }
    
    isLeft = 0;
    isTrans = 0;
    fk::matrix<TestType> productRightFull = mat1*fmwt;
    auto const productRight = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans,method);
    
    SECTION("degree = 4, lev 5 fmwt apply right") { relaxed_comparison(productRightFull, productRight); }
    
    isLeft = 0;
    isTrans = 1;
    fk::matrix<TestType> productRightTransFull = mat1*fmwt_transpose;
    auto const productRightTrans = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans,method);
    //productLeftFull.print("productLeftTransFull");	
    //productLeft.print("productLeftTrans");	
   
    
    //fk::matrix<TestType> sub1 = productLeftTransFull.extract_submatrix(0,0,4,4);
    //fk::matrix<TestType> sub2 = productLeftTrans.extract_submatrix(0,0,4,4);
    //sub1.print("sub1");
    //sub2.print("sub2");
    SECTION("degree = 4, lev 5 fmwt apply right transpose") { relaxed_comparison(productRightTrans, productRightTransFull); }
  
  
  }
  SECTION("Apply fmwt test 2")
  {
    int const kdeg = 2;
    int const lev = 3;
    std::string out_base =
        "../testing/generated-inputs/apply_fmwt/";

    std::string mat1_string    = out_base + "mat1.dat";
    //std::string fmwt_string    = out_base + "fmwt_k4_lev5.dat";

    fk::matrix<TestType> mat1 =
        fk::matrix<TestType>(read_matrix_from_txt_file(mat1_string));
    dimension const dim =
        make_PDE<TestType>(PDE_opts::continuity_1, lev, kdeg)
            ->get_dimensions()[0];

    fk::matrix<TestType> const fmwt = operator_two_scale<TestType>(dim);
    int const n = pow(2,lev);
    fk::matrix<TestType> shape(kdeg*n,kdeg*n);
    fk::matrix<TestType> one(1,1);
    one(0,0) = 1.0;
    one.print("one");
    for(int i=0;i<kdeg;i++)
    {
      for(int j=0;j<kdeg*n;j++)
      {
        shape.set_submatrix(i,j,one);
      }
    }
    
    int sum=1;
    for(int i=0;i<lev;i++)
    {
      std::cout << "i " << i << std::endl;
      for(int j=0;j<pow(2,i);j++)
      {
      std::cout << "j " << j << std::endl;
	int row_start = sum*kdeg;
	int row_end = (sum+1)*kdeg;
	int col_start = j*kdeg*n/pow(2,i);
	int col_end = (j+1)*kdeg*n/pow(2,i);
        sum = sum+1;
    
        for(int k=row_start;k<row_end;k++)
        {
      std::cout << "k " << k << std::endl;
          for(int l=col_start;l<col_end;l++)
          {
      std::cout << "l " << l << std::endl;
            shape.set_submatrix(k,l,one);
          }
        }


      }
    }

    shape.print("shape");
    
      std::cout.precision(17);
    for(int i=0;i<kdeg*n;i++)
    {
      for(int j=0;j<kdeg*n;j++)
      {
        if(shape(i,j) < 1.0)
	{
	  std::cout << "fmwt i j value " << i << " " << j << " " << fmwt(i,j) << std::endl;
	}
      }
    }

   std::cout << "finished printing " << std::endl; 
    //fk::matrix<TestType> fmwt =
    //    fk::matrix<TestType>(read_matrix_from_txt_file(fmwt_string));
    
    //mat1.print("mat1");	
    //fmwt.print("fmwt");	
    int isLeft = 1;
    int isTrans = 0;
    int method = 2;
    fk::matrix<TestType> productLeftFull = fmwt*mat1;
    auto const productLeft = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans,method);
    
    SECTION("degree = 4, lev 5 fmwt apply left") { relaxed_comparison(productLeftFull, productLeft); }
    //
    isLeft = 1;
    isTrans = 1;
    fk::matrix<TestType> fmwt_transpose = fk::matrix<TestType>(fmwt).transpose();
    fk::matrix<TestType> productLeftTransFull = fmwt_transpose*mat1;
    auto const productLeftTrans = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans,method);
    //productLeftFull.print("productLeftTransFull");	
    //productLeft.print("productLeftTrans");	
    fk::matrix<TestType> difference = productLeftTrans - productLeftTransFull;
    TestType maxdiff = *std::max_element(difference.begin(),difference.end());
    //difference.print("difference");
    TestType x1;
    std::cout << "print " << typeid(x1).name() << std::endl; 
    std::cout << "max diff " << maxdiff << std::endl; 
    //fk::matrix<TestType> sub1 = productLeftTransFull.extract_submatrix(0,0,4,4);
    //fk::matrix<TestType> sub2 = productLeftTrans.extract_submatrix(0,0,4,4);
    //sub1.print("sub1");
    //sub2.print("sub2");
    SECTION("degree = 4, lev 5 fmwt apply left transpose") { relaxed_comparison(productLeftTrans, productLeftTransFull); }
    //SECTION("degree = 4, lev 5 fmwt apply left transpose") { relaxed_comparison(sub1,sub2); }
    
    
    isLeft = 0;
    isTrans = 0;
    fk::matrix<TestType> productRightFull = mat1*fmwt;
    auto const productRight = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans,method);
    
    SECTION("degree = 4, lev 5 fmwt apply right") { relaxed_comparison(productRightFull, productRight); }
    
    isLeft = 0;
    isTrans = 1;
    fk::matrix<TestType> productRightTransFull = mat1*fmwt_transpose;
    auto const productRightTrans = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans,method);
    //productLeftFull.print("productLeftTransFull");	
    //productLeft.print("productLeftTrans");	
   
    
    //fk::matrix<TestType> sub1 = productLeftTransFull.extract_submatrix(0,0,4,4);
    //fk::matrix<TestType> sub2 = productLeftTrans.extract_submatrix(0,0,4,4);
    //sub1.print("sub1");
    //sub2.print("sub2");
    SECTION("degree = 4, lev 5 fmwt apply right transpose") { relaxed_comparison(productRightTrans, productRightTransFull); }
  
  }

}
