#pragma once

//-----------------------------------------------------------------------------
//
// a collection of utility functions to make working with matlab-like
// code a little easier
//
// - matlab functions ported to c++
//    - linspace (scalar inputs only)
//
// - matlab/octave file IO
//    - readVectorFromTxtFile (tested for octave)
//    - readVectorFromBinFile (tested for octave and matlab)
//
//-----------------------------------------------------------------------------

#include "tensors.hpp"
#include <string>
#include <vector>

// matlab's "linspace(start, end, N)" function
//-----------------------------------------------------------------------------
//
// c++ implementation of matlab (a subset of) linspace() function
// initial c++ implementation by Tyler McDaniel
//
// -- linspace (START, END)
// -- linspace (START, END, N)
//     Return a row vector with N linearly spaced elements between START
//     and END.
//
//     If the number of elements is greater than one, then the endpoints
//     START and END are always included in the range.  If START is
//     greater than END, the elements are stored in decreasing order.  If
//     the number of points is not specified, a value of 100 is used.
//
//     The 'linspace' function returns a row vector when both START and
//     END are scalars.
//
//  (unsupported)
//     If one, or both, inputs are vectors, then
//     'linspace' transforms them to column vectors and returns a matrix
//     where each row is an independent sequence between
//     'START(ROW_N), END(ROW_N)'.
//
//     For compatibility with MATLAB, return the second argument (END)
//     when only a single value (N = 1) is requested.
//
//-----------------------------------------------------------------------------
template<typename P>
std::enable_if_t<std::is_floating_point<P>::value, fk::vector<P>>
linspace(P const start, P const end, unsigned int const num_elems = 100)
{
  assert(num_elems > 1); // must have at least 2 elements

  // create output vector
  fk::vector<P> points(num_elems);

  // find interval size
  P const interval_size = (end - start) / (num_elems - 1);

  // insert first and last elements
  points(0)             = start;
  points(num_elems - 1) = end;

  // fill in the middle
  for (unsigned int i = 1; i < num_elems - 1; ++i)
  {
    points(i) = start + i * interval_size;
  }

  return points;
}
//-----------------------------------------------------------------------------
//
// c++ implementation of (a subset of) eye() function
// The following are not supported here:
// - providing a third "CLASS" argument
// - providing a vector argument for the dimensions
//
// -- eye (N)
// -- eye (M, N)
// -- eye ([M N])
//     Return an identity matrix.
//
//     If invoked with a single scalar argument N, return a square NxN
//     identity matrix.
//
//     If supplied two scalar arguments (M, N), 'eye' takes them to be the
//     number of rows and columns.  If given a vector with two elements,
//     'eye' uses the values of the elements as the number of rows and
//     columns, respectively.  For example:
//
//          eye (3)
//           =>  1  0  0
//               0  1  0
//               0  0  1
//
//     The following expressions all produce the same result:
//
//          eye (2)
//          ==
//          eye (2, 2)
//          ==
//          eye (size ([1, 2; 3, 4]))
//
//     Calling 'eye' with no arguments is equivalent to calling it with an
//     argument of 1.  Any negative dimensions are treated as zero.  These
//     odd definitions are for compatibility with MATLAB.
//
//-----------------------------------------------------------------------------
template<typename P>
fk::matrix<P> eye(int const M = 1)
{
  fk::matrix<P> id(M, M);
  for (auto i = 0; i < M; ++i)
    id(i, i) = 1.0;
  return id;
}
template<typename P>
fk::matrix<P> eye(int const M, int const N)
{
  fk::matrix<P> id(M, N);
  for (auto i = 0; i < (M < N ? M : N); ++i)
    id(i, i) = 1.0;
  return id;
}

// read a matlab vector from binary file into a std::vector
// note that fk::vector has a copy assignment overload from std::vector
fk::vector<double> readVectorFromBinFile(std::string const &path);

// read an octave vector from text file into a std::vector
// note that fk::vector has a copy assignment overload from std::vector
fk::vector<double> readVectorFromTxtFile(std::string const &path);

// read an octave matrix from text file into a fk::matrix
fk::matrix<double> readMatrixFromTxtFile(std::string const &path);
