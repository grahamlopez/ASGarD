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
      //std::cout.precision(17);
      //std::cout << "first and second element " << std::fixed << f1
      //             << " " << std::fixed << second_elem << std::endl;
      //std::cout << "diff " << std::fixed << std::abs(f1 - second_elem) << std::endl;
      //std::cout << "epsilon  "<< std::numeric_limits<TestType>::epsilon() << std::endl;
      //std::cout << "max  "<< std::max(std::abs(f1),std::abs(second_elem)) << std::endl;
      //std::cout << "compare " << std::numeric_limits<TestType>::epsilon()*1.0e4*std::max(std::abs(f1),std::abs(second_elem)) << std::endl;
      auto const diff = std::abs(f1 - second_elem);///std::max(std::abs(f1),std::abs(second_elem));
      auto const abs_diff = std::abs(f1 - second_elem);
      REQUIRE(diff < std::numeric_limits<TestType>::epsilon()*1.0e5*std::max(std::abs(f1),std::abs(second_elem)));
      //REQUIRE(Approx(f1)
      //            .margin(std::numeric_limits<TestType>::epsilon()) ==
      //        second_elem);
      //REQUIRE(std::abs((*first_it++) - second_elem) < 1.0e-5*std::max(std::abs((*first_it++)),std::abs((second_elem))));
    });
  };
  
  SECTION("Apply fmwt test 1")
  {
    std::cout << " type epsilon " << std::numeric_limits<TestType>::epsilon() << std::endl;
    int kdeg = 2;
    int lev = 2;
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
    
  
  
    isLeft = 1;
    isTrans = 1;
    method=2;
    fk::matrix<TestType> fmwt_transpose = fk::matrix<TestType>(fmwt).transpose();
    fk::matrix<TestType> productLeftTransFull = fmwt_transpose*mat1;
    auto const productLeftTrans = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans,method);
    SECTION("degree = 2, lev 2 fmwt apply left transpose method 2") { relaxed_comparison(productLeftTrans, productLeftTransFull); }
    
    isLeft = 0;
    isTrans = 0;
    fk::matrix<TestType> productRightFull = mat1*fmwt;
    auto const productRight = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans,method);
    
    SECTION("degree = 2, lev 2 fmwt apply right method 2") { relaxed_comparison(productRightFull, productRight); }
    
    isLeft = 0;
    isTrans = 1;
    fk::matrix<TestType> productRightTransFull = mat1*fmwt_transpose;
    auto const productRightTrans = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans,method);
    SECTION("degree = 2, lev 2 fmwt apply right transpose method 2") { relaxed_comparison(productRightTrans, productRightTransFull); }
    
    
    method=3;
    isLeft = 1;
    isTrans = 0;
    auto const productLeft3 = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans,method);
    productLeftFull.print("full");	
    productLeft3.print("apply left 3");	
    
    SECTION("degree = 2, lev 2 fmwt apply left method 3") { relaxed_comparison(productLeftFull, productLeft3); }
    
    isLeft = 1;
    isTrans = 1;
    auto const productLeftTrans3 = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans,method);
    productLeftTransFull.print("full");	
    productLeftTrans3.print("apply left 3");	
    
    SECTION("degree = 2, lev 2 fmwt apply left method 3") { relaxed_comparison(productLeftTransFull, productLeftTrans3); }
    
    isLeft = 0;
    isTrans = 0;
    auto const productRight3 = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans,method);
    productRightFull.print("full");	
    productRight3.print("apply left 3");	
    
    SECTION("degree = 2, lev 2 fmwt apply right method 3") { relaxed_comparison(productRightFull, productRight3); }
    
    isLeft = 0;
    isTrans = 1;
    auto const productRightTrans3 = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans,method);
    productRightTransFull.print("full");	
    productRightTrans3.print("apply right trans 3");	
    
    SECTION("degree = 2, lev 2 fmwt apply right trans method 3") { relaxed_comparison(productRightTransFull, productRightTrans3); }
  }
  
  SECTION("Apply fmwt test 2")
  {
    int kdeg = 4;
    int lev = 5;
    std::string out_base =
        "../testing/generated-inputs/apply_fmwt/";

    std::string mat1_string    = out_base + "mat1.dat";
    
    dimension const dim =
        make_PDE<TestType>(PDE_opts::continuity_1, lev, kdeg)
            ->get_dimensions()[0];

    fk::matrix<TestType> const fmwt = operator_two_scale<TestType>(dim);

    fk::matrix<TestType> mat1 =
        fk::matrix<TestType>(read_matrix_from_txt_file(mat1_string));
    
    int isLeft = 1;
    int isTrans = 0;
    int method = 2;
    fk::matrix<TestType> productLeftFull = fmwt*mat1;
    auto const productLeft = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans,method);
    productLeftFull.print("full");	
    productLeft.print("apply left");	
    
    SECTION("degree = 4, lev 5 fmwt apply left") { relaxed_comparison(productLeftFull, productLeft); }
    
  
     
    isLeft = 1;
    isTrans = 1;
    method=2;
    fk::matrix<TestType> fmwt_transpose = fk::matrix<TestType>(fmwt).transpose();
    fk::matrix<TestType> productLeftTransFull = fmwt_transpose*mat1;
    auto const productLeftTrans = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans,method);
    SECTION("degree = 2, lev 2 fmwt apply left transpose method 2") { relaxed_comparison(productLeftTrans, productLeftTransFull); }
  
    isLeft = 0;
    isTrans = 0;
    fk::matrix<TestType> productRightFull = mat1*fmwt;
    auto const productRight = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans,method);
    
    SECTION("degree = 2, lev 2 fmwt apply right method 2") { relaxed_comparison(productRightFull, productRight); }
    
    isLeft = 0;
    isTrans = 1;
    fk::matrix<TestType> productRightTransFull = mat1*fmwt_transpose;
    auto const productRightTrans = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans,method);
    SECTION("degree = 2, lev 2 fmwt apply right transpose method 2") { relaxed_comparison(productRightTrans, productRightTransFull); }
    
   
    method=3;
    isLeft = 1;
    isTrans = 0;
    auto const productLeft3 = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans,method);
    productLeftFull.print("full");	
    productLeft3.print("apply left 3");	
    
    SECTION("degree = 2, lev 2 fmwt apply left method 3") { relaxed_comparison(productLeftFull, productLeft3); }
    
    isLeft = 1;
    isTrans = 1;
    auto const productLeftTrans3 = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans,method);
    productLeftTransFull.print("full");	
    productLeftTrans3.print("apply left 3");	
    
    SECTION("degree = 2, lev 2 fmwt apply left method 3") { relaxed_comparison(productLeftTransFull, productLeftTrans3); }
    
    isLeft = 0;
    isTrans = 0;
    auto const productRight3 = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans,method);
    productRightFull.print("full");	
    productRight3.print("apply left 3");	
    
    SECTION("degree = 2, lev 2 fmwt apply right method 3") { relaxed_comparison(productRightFull, productRight3); }
    
    isLeft = 0;
    isTrans = 1;
    auto const productRightTrans3 = apply_fmwt<TestType>(fmwt,mat1,kdeg,lev,isLeft,isTrans,method);
    productRightTransFull.print("full");	
    productRightTrans3.print("apply right trans 3");	
    
    SECTION("degree = 2, lev 2 fmwt apply right trans method 3") { relaxed_comparison(productRightTransFull, productRightTrans3); }
  
  }

}
