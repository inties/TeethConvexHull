// Copyright (c) 1997-2000
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL: https://github.com/CGAL/cgal/blob/v5.6.2/Kernel_d/include/CGAL/Linear_algebraHd.h $
// $Id: Linear_algebraHd.h 45478184de2 2022-11-15T13:39:40+01:00 albert-github
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>

//---------------------------------------------------------------------
// file generated by notangle from Linear_algebra.lw
// please debug or modify noweb file
// based on LEDA architecture by S. Naeher, C. Uhrig
// coding: K. Mehlhorn, M. Seel
// debugging and templatization: M. Seel
//---------------------------------------------------------------------

#ifndef CGAL_LINEAR_ALGEBRAHD_H
#define CGAL_LINEAR_ALGEBRAHD_H

#include <CGAL/Kernel_d/Vector__.h>
#include <CGAL/Kernel_d/Matrix__.h>

// #define CGAL_LA_SELFTEST
namespace CGAL {

/*{\Moptions outfile=Linear_algebra.man}*/
/*{\Manpage {Linear_algebraHd}{RT}{Linear Algebra on RT}{LA}}*/

template <class RT_, class AL_ = CGAL_ALLOCATOR(RT_) >
class Linear_algebraHd
{
/*{\Mdefinition
The data type |\Mname| encapsulates two classes |Matrix|, |Vector|
and many functions of basic linear algebra. It is parametrized by a
number type |RT|. An instance of data type |Matrix| is a matrix of
variables of type |RT|, the so called ring type. Accordingly,
|Vector| implements vectors of variables of type |RT|. The arithmetic
type |RT| is required to behave like integers in the mathematical
sense. The manual pages of |Vector| and |Matrix| follow below.

All functions compute the exact result, i.e., there is no rounding
error.  Most functions of linear algebra are \emph{checkable}, i.e.,
the programs can be asked for a proof that their output is
correct. For example, if the linear system solver declares a linear
system $A x = b$ unsolvable it also returns a vector $c$ such that
$c^T A = 0$ and $c^T b \neq 0$.  All internal correctness checks can
be switched on by the flag [[CGAL_LA_SELFTEST]].}*/

public:

/*{\Mtypes 5.5}*/

typedef RT_ RT;
/*{\Mtypemember the ring type of the components.}*/

typedef Linear_Algebra::Vector_<RT_,AL_> Vector;
/*{\Mtypemember the vector type.}*/

typedef Linear_Algebra::Matrix_<RT_,AL_> Matrix;
/*{\Mtypemember the matrix type.}*/

typedef AL_ allocator_type;
/*{\Mtypemember the allocator used for memory management. |\Mname| is
an abbreviation for |Linear_algebraHd<RT, ALLOC = allocator<RT,LA> >|. Thus
|allocator_type| defaults to the standard allocator offered by the STL.}*/

/*{\Moperations 2 1}*/

static Matrix  transpose(const Matrix& M);
/*{\Mstatic  returns  $M^T$ ($m\times n$ - matrix). }*/

static bool inverse(const Matrix& M, Matrix& I, RT& D, Vector& c);
/*{\Mstatic determines whether |M| has an inverse. It also computes
either the inverse as $(1/D) \cdot |I|$ or when no inverse
exists, a vector $c$ such that $c^T \cdot M = 0 $.  }*/

static Matrix  inverse(const Matrix& M, RT& D)
/*{\Mstatic returns the inverse matrix of |M|. More precisely, $1/D$
            times the matrix returned is the inverse of |M|.\\
            \precond  |determinant(M) != 0|. }*/
{
  Matrix result;
  Vector c;
  if (!inverse(M,result,D,c))
    CGAL_error_msg("inverse(): matrix is singular.");
  return result;
}

static RT  determinant (const Matrix& M, Matrix& L, Matrix& U,
                        std::vector<int>& q, Vector& c);
/*{\Mstatic returns the determinant $D$ of |M| and sufficient information
            to verify that the value of the determinant is correct. If
            the determinant is zero then $c$ is a vector such that
            $c^T \cdot M = 0$. If the determinant is non-zero then $L$
            and $U$ are lower and upper diagonal matrices respectively
            and $q$ encodes a permutation matrix $Q$ with $Q(i,j) = 1$
            iff $i = q(j)$ such that $L \cdot M \cdot Q = U$,
            $L(0,0) = 1$, $L(i,i) = U(i - 1,i - 1)$ for all $i$,
            $1 \le i < n$, and $D = s \cdot U(n - 1,n - 1)$ where $s$ is
            the determinant of $Q$. \precond  |M| is square. }*/

static bool verify_determinant (const Matrix& M, RT D, Matrix& L, Matrix& U,
                                const std::vector<int>& q, Vector& c);
/*{\Mstatic verifies the conditions stated above. }*/

static RT determinant (const Matrix& M);
/*{\Mstatic  returns the determinant of |M|.
         \precond  |M| is square. }*/

static int sign_of_determinant (const Matrix& M);
/*{\Mstatic returns the sign of the determinant of |M|.
        \precond  |M| is square. }*/

static bool linear_solver(const Matrix& M, const Vector& b,
                          Vector& x, RT& D,
                          Matrix& spanning_vectors,
                          Vector& c);
/*{\Mstatic determines the complete solution space of the linear system
            $M\cdot x = b$. If the system is unsolvable then
            $c^T \cdot M = 0$ and $c^T \cdot b \not= 0$.
            If the system is solvable then $(1/D) x$ is a solution, and
            the columns of |spanning_vectors| are a maximal set of linearly
            independent solutions to the corresponding homogeneous system.
            \precond |M.row_dimension() = b.dimension()|. }*/

static bool linear_solver(const Matrix& M, const Vector& b,
                          Vector& x, RT& D,
                          Vector& c)
/*{\Mstatic determines whether the linear system $M\cdot x = b$ is
           solvable. If yes, then $(1/D) x$ is a solution, if not then
           $c^T \cdot M = 0$ and $c^T \cdot b \not= 0$.
           \precond |M.row_dimension() = b.dimension()|. }*/
{
  Matrix spanning_vectors;
  return linear_solver(M,b,x,D,spanning_vectors,c);
}

static bool linear_solver(const Matrix& M, const Vector& b,
                          Vector& x, RT& D)
/*{\Mstatic as above, but without the witness $c$
           \precond |M.row_dimension() = b.dimension()|. }*/
{
  Matrix spanning_vectors; Vector c;
  return linear_solver(M,b,x,D,spanning_vectors,c);
}

static bool is_solvable(const Matrix& M, const Vector& b)
/*{\Mstatic determines whether the system $M \cdot x = b$ is solvable \\
        \precond |M.row_dimension() = b.dimension()|. }*/
{
  Vector x; RT D; Matrix spanning_vectors; Vector c;
  return linear_solver(M,b,x,D,spanning_vectors,c);
}

static bool homogeneous_linear_solver (const Matrix& M, Vector& x);
/*{\Mstatic determines whether the homogeneous linear system
        $M\cdot x = 0$ has a non - trivial solution. If
        yes, then $x$ is such a solution. }*/

static int homogeneous_linear_solver (const Matrix& M, Matrix& spanning_vecs);
/*{\Mstatic determines the solution space of the homogeneous linear
system $M\cdot x = 0$. It returns the dimension of the solution space.
Moreover the columns of |spanning_vecs| span the solution space. }*/

static int independent_columns (const Matrix& M, std::vector<int>& columns);
/*{\Mstatic returns the indices of a maximal subset of independent
columns of |M|.}*/

static int rank (const Matrix & M);
/*{\Mstatic returns the rank of matrix |M| }*/

/*{\Mimplementation The datatype |\Mname| is a wrapper class for the
linear algebra functionality on matrices and vectors.  Operations
|determinant|, |inverse|, |linear_solver|, and |rank| take time
$O(n^3)$, and all other operations take time $O(nm)$. These time
bounds ignore the cost for multiprecision arithmetic operations.

All functions on integer matrices compute the exact result, i.e.,
there is no rounding error. The implementation follows a proposal of
J. Edmonds (J. Edmonds, Systems of distinct representatives and linear
algebra, Journal of Research of the Bureau of National Standards, (B),
71, 241 - 245). Most functions of linear algebra are { \em checkable
}, i.e., the programs can be asked for a proof that their output is
correct. For example, if the linear system solver declares a linear
system $A x = b$ unsolvable it also returns a vector $c$ such that
$c^T A = 0$ and $c^T b \not= 0$.}*/

}; // Linear_algebraHd


} //namespace CGAL

#include <CGAL/Kernel_d/Linear_algebraHd_impl.h>

#endif // CGAL_LINALG_ALGEBRAHD_H
