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
      REQUIRE(Approx(*first_it++)
                  .margin(std::numeric_limits<TestType>::epsilon()*0.1) ==
              second_elem);
    });
  };
  /*
  SECTION("Apply fmwt test 1")
  {
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
    fk::matrix<TestType> productLeftFull = fmwt*mat1;
    auto const productLeft = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans);
    
    SECTION("degree = 4, lev 5 fmwt apply left") { relaxed_comparison(productLeftFull, productLeft); }
    //
    isLeft = 1;
    isTrans = 1;
    fk::matrix<TestType> fmwt_transpose = fk::matrix<TestType>(fmwt).transpose();
    fk::matrix<TestType> productLeftTransFull = fmwt_transpose*mat1;
    auto const productLeftTrans = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans);
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
    auto const productRight = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans);
    
    SECTION("degree = 4, lev 5 fmwt apply right") { relaxed_comparison(productRightFull, productRight); }
    
    isLeft = 0;
    isTrans = 1;
    fk::matrix<TestType> productRightTransFull = mat1*fmwt_transpose;
    auto const productRightTrans = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans);
    //productLeftFull.print("productLeftTransFull");	
    //productLeft.print("productLeftTrans");	
   
    
    //fk::matrix<TestType> sub1 = productLeftTransFull.extract_submatrix(0,0,4,4);
    //fk::matrix<TestType> sub2 = productLeftTrans.extract_submatrix(0,0,4,4);
    //sub1.print("sub1");
    //sub2.print("sub2");
    SECTION("degree = 4, lev 5 fmwt apply right transpose") { relaxed_comparison(productRightTrans, productRightTransFull); }
  }
  */
  SECTION("Apply fmwt test 2")
  {
    int const kdeg = 4;
    int const lev = 5;
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
    //fk::matrix<TestType> fmwt =
    //    fk::matrix<TestType>(read_matrix_from_txt_file(fmwt_string));
    
    //mat1.print("mat1");	
    //fmwt.print("fmwt");	
    int isLeft = 1;
    int isTrans = 0;
    fk::matrix<TestType> productLeftFull = fmwt*mat1;
    auto const productLeft = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans);
    
    SECTION("degree = 4, lev 5 fmwt apply left") { relaxed_comparison(productLeftFull, productLeft); }
    //
    isLeft = 1;
    isTrans = 1;
    fk::matrix<TestType> fmwt_transpose = fk::matrix<TestType>(fmwt).transpose();
    fk::matrix<TestType> productLeftTransFull = fmwt_transpose*mat1;
    auto const productLeftTrans = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans);
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
    /*
    isLeft = 0;
    isTrans = 0;
    fk::matrix<TestType> productRightFull = mat1*fmwt;
    auto const productRight = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans);
    
    SECTION("degree = 4, lev 5 fmwt apply right") { relaxed_comparison(productRightFull, productRight); }
    
    isLeft = 0;
    isTrans = 1;
    fk::matrix<TestType> productRightTransFull = mat1*fmwt_transpose;
    auto const productRightTrans = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans);
    //productLeftFull.print("productLeftTransFull");	
    //productLeft.print("productLeftTrans");	
   
    
    //fk::matrix<TestType> sub1 = productLeftTransFull.extract_submatrix(0,0,4,4);
    //fk::matrix<TestType> sub2 = productLeftTrans.extract_submatrix(0,0,4,4);
    //sub1.print("sub1");
    //sub2.print("sub2");
    SECTION("degree = 4, lev 5 fmwt apply right transpose") { relaxed_comparison(productRightTrans, productRightTransFull); }
  */
  }

}
