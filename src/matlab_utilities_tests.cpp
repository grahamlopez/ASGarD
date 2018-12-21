
#include "matlab_utilities.hpp"

#include "tests_general.hpp"
#include <vector>

TEMPLATE_TEST_CASE("linspace() matches matlab implementation", "[matlab]",
                   float, double)
{
  SECTION("linspace(0,1) returns 100 elements")
  {
    fk::vector<TestType> test = linspace<TestType>(0, 1);
    REQUIRE(test.size() == 100);
  }
  SECTION("linspace(-1,1,9)")
  {
    fk::vector<TestType> gold = readVectorFromTxtFile(
        "../testing/generated-inputs/linspace_neg1_1_9.dat");
    REQUIRE(gold.size() == 9);
    fk::vector<TestType> test = linspace<TestType>(-1, 1, 9);
    REQUIRE(test == gold);
  }
  SECTION("linspace(1,-1,9)")
  {
    fk::vector<TestType> gold = readVectorFromTxtFile(
        "../testing/generated-inputs/linspace_1_neg1_9.dat");
    REQUIRE(gold.size() == 9);
    fk::vector<TestType> test = linspace<TestType>(1, -1, 9);
    REQUIRE(test == gold);
  }
  SECTION("linspace(-1,1,8)")
  {
    fk::vector<TestType> gold = readVectorFromTxtFile(
        "../testing/generated-inputs/linspace_neg1_1_8.dat");
    REQUIRE(gold.size() == 8);
    fk::vector<TestType> test = linspace<TestType>(-1, 1, 8);
    REQUIRE(test == gold);
  }
}

// using widening conversions for golden data in order to test integers
// FIXME look for another way
TEMPLATE_TEST_CASE("eye() matches matlab implementation", "[matlab]", float,
                   double, int)
{
  SECTION("eye()")
  {
    fk::matrix<TestType> gold{{1}};
    fk::matrix<TestType> test = eye<TestType>();
    REQUIRE(test == gold);
  }
  SECTION("eye(5)")
  {
    // clang-format off
    fk::matrix<TestType> gold{
      {1, 0, 0, 0, 0},
      {0, 1, 0, 0, 0},
      {0, 0, 1, 0, 0},
      {0, 0, 0, 1, 0},
      {0, 0, 0, 0, 1},
    }; // clang-format on
    fk::matrix<TestType> test = eye<TestType>(5);
    REQUIRE(test == gold);
  }
  SECTION("eye(5,5)")
  {
    // clang-format off
    fk::matrix<TestType> gold{
      {1, 0, 0, 0, 0},
      {0, 1, 0, 0, 0},
      {0, 0, 1, 0, 0},
      {0, 0, 0, 1, 0},
      {0, 0, 0, 0, 1},
    }; // clang-format on
    fk::matrix<TestType> test = eye<TestType>(5, 5);
    REQUIRE(test == gold);
  }
  SECTION("eye(5,3)")
  {
    // clang-format off
    fk::matrix<TestType> gold{
      {1, 0, 0},
      {0, 1, 0},
      {0, 0, 1},
      {0, 0, 0},
      {0, 0, 0},
    }; // clang-format on
    fk::matrix<TestType> test = eye<TestType>(5, 3);
    REQUIRE(test == gold);
  }
  SECTION("eye(3,5)")
  {
    // clang-format off
    fk::matrix<TestType> gold{
      {1, 0, 0, 0, 0},
      {0, 1, 0, 0, 0},
      {0, 0, 1, 0, 0},
    }; // clang-format on
    fk::matrix<TestType> test = eye<TestType>(3, 5);
    REQUIRE(test == gold);
  }
}
TEST_CASE("readVectorFromBinFile returns expected vector", "[matlab]")
{
  SECTION("readVectorFromBinFile gets 100-element row vector")
  {
    fk::vector<double> gold = linspace<double>(-1, 1);
    fk::vector<double> test = readVectorFromBinFile(
        "../testing/generated-inputs/readVectorBin_neg1_1_100.dat");
    REQUIRE(test == gold);
  }
  SECTION("readVectorFromBinFile gets 100-element column vector")
  {
    fk::vector<double> gold = linspace<double>(-1, 1);
    fk::vector<double> test = readVectorFromBinFile(
        "../testing/generated-inputs/readVectorBin_neg1_1_100T.dat");
    REQUIRE(test == gold);
  }
  SECTION("readVectorFromBinFile fails on non-existent path")
  {
    fk::vector<double> test = readVectorFromBinFile("this/path/does/not/exist");
    REQUIRE(test.size() == 0);
  }
}

TEST_CASE("readVectorFromTxtFile returns expected vector", "[matlab]")
{
  SECTION("readVectorFromTxtFile gets 100-element row vector")
  {
    fk::vector<double> gold = linspace<double>(-1, 1);
    fk::vector<double> test = readVectorFromTxtFile(
        "../testing/generated-inputs/readVectorTxt_neg1_1_100.dat");
    REQUIRE(test == gold);
  }
  SECTION("readVectorFromTxtFile gets 100-element column vector")
  {
    fk::vector<double> gold = linspace<double>(-1, 1);
    fk::vector<double> test = readVectorFromTxtFile(
        "../testing/generated-inputs/readVectorTxt_neg1_1_100T.dat");
    REQUIRE(test == gold);
  }
  SECTION("readVectorFromTxtFile fails on non-existent path")
  {
    fk::vector<double> test = readVectorFromTxtFile("this/path/does/not/exist");
    REQUIRE(test.size() == 0);
  }
}

TEST_CASE("readMatrixFromTxtFile returns expected vector", "[matlab]")
{
  SECTION("readMatrixFromTxtFile gets 5,5 matrix")
  {
    auto gold = fk::matrix<double>(5, 5);
    // generate the golden matrix
    for (int i = 0; i < 5; i++)
      for (int j = 0; j < 5; j++)
        gold(i, j) = 17.0 / (i + 1 + j);

    fk::matrix<double> test = readMatrixFromTxtFile(
        "../testing/generated-inputs/readMatrixTxt_5x5.dat");
    REQUIRE(test == gold);
  }
  SECTION("readMatrixFromTxtFile fails on non-existent path")
  {
    fk::matrix<double> test = readMatrixFromTxtFile("this/path/does/not/exist");
    REQUIRE(test.size() == 0);
  }
}
