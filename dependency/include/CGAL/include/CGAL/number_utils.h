// Copyright (c) 1999
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL: https://github.com/CGAL/cgal/blob/v5.6.2/Algebraic_foundations/include/CGAL/number_utils.h $
// $Id: number_utils.h 3fa436459d4 2022-06-10T08:41:19+02:00 Sébastien Loriot
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stefan Schirra

#ifndef CGAL_NUMBER_UTILS_H
#define CGAL_NUMBER_UTILS_H

#include <CGAL/number_type_config.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Real_embeddable_traits.h>
#include <CGAL/Kernel/Same_uncertainty.h>

namespace CGAL {
CGAL_NTS_BEGIN_NAMESPACE


// AST-Functor adapting functions UNARY
template< class AS >
inline
void
simplify( AS& x ) {
    typename Algebraic_structure_traits< AS >::Simplify simplify;
    simplify( x );
}

template< class AS >
inline
typename Algebraic_structure_traits< AS >::Unit_part::result_type
unit_part( const AS& x ) {
    typename Algebraic_structure_traits< AS >::Unit_part unit_part;
    return unit_part( x );
}


template< class AS >
inline
typename Algebraic_structure_traits< AS >::Is_square::result_type
is_square( const AS& x,
           typename Algebraic_structure_traits< AS >::Is_square::second_argument_type y )
{
    typename Algebraic_structure_traits< AS >::Is_square is_square;
    return is_square( x, y );
}

template< class AS >
inline
typename Algebraic_structure_traits< AS >::Is_square::result_type
is_square( const AS& x){
    typename Algebraic_structure_traits< AS >::Is_square is_square;
    return is_square( x );
}


template< class AS >
inline
typename Algebraic_structure_traits< AS >::Square::result_type
square( const AS& x ) {
    typename Algebraic_structure_traits< AS >::Square square;
    return square( x );
}


template< class AS >
inline
typename Algebraic_structure_traits< AS >::Inverse::result_type
inverse( const AS& x ) {
    typename Algebraic_structure_traits< AS >::Inverse inverse;
    return inverse( x );
}

template< class AS >
inline
typename Algebraic_structure_traits<AS>::Is_one::result_type
is_one( const AS& x ) {
    typename Algebraic_structure_traits< AS >::Is_one is_one;
    return is_one( x );
}

template< class AS >
inline
typename Algebraic_structure_traits< AS >::Sqrt::result_type
sqrt( const AS& x ) {
    typename Algebraic_structure_traits< AS >::Sqrt sqrt;
    return sqrt( x );
}

// AST-Functor adapting functions BINARY

template< class A, class B >
inline
typename Algebraic_structure_traits< typename Coercion_traits<A,B>::Type>
::Integral_division::result_type
integral_division( const A& x, const B& y ) {
    typedef typename Coercion_traits<A,B>::Type Type;
    typename Algebraic_structure_traits< Type >::Integral_division
        integral_division;
    return integral_division( x, y );
}

template< class A, class B >
inline
typename Algebraic_structure_traits< typename Coercion_traits<A,B>::Type>
::Divides::result_type
divides( const A& x, const B& y ) {
  typedef typename Coercion_traits<A,B>::Type Type;
  typename Algebraic_structure_traits< Type >::Divides  divides;
  return divides( x, y );
}

template< class Type >
inline
typename Algebraic_structure_traits<Type>::Divides::result_type
divides( const Type& x, const Type& y, Type& q ) {
  typename Algebraic_structure_traits< Type >::Divides  divides;
  return divides( x, y, q);
}

template< class A, class B >
inline
typename Algebraic_structure_traits< typename Coercion_traits<A,B>::Type >
::Gcd::result_type
gcd( const A& x, const B& y ) {
    typedef typename Coercion_traits<A,B>::Type      Type;
    typename Algebraic_structure_traits< Type >::Gcd gcd;
    return gcd( x, y );
}


template< class A, class B >
inline
typename Algebraic_structure_traits< typename Coercion_traits<A,B>::Type >
::Mod::result_type
mod( const A& x, const B& y ) {
    typedef typename Coercion_traits<A,B>::Type Type;
    typename Algebraic_structure_traits<Type >::Mod mod;
    return mod( x, y );
}

template< class A, class B >
inline
typename Algebraic_structure_traits< typename Coercion_traits<A,B>::Type>::Div::result_type
div( const A& x, const B& y ) {
    typedef typename Coercion_traits<A,B>::Type Type;
    typename Algebraic_structure_traits<Type >::Div div;
    return div( x, y );
}

template< class A, class B >
inline
void
div_mod(
        const A& x,
        const B& y,
        typename Coercion_traits<A,B>::Type& q,
        typename Coercion_traits<A,B>::Type& r ) {
    typedef typename Coercion_traits<A,B>::Type Type;
    typename Algebraic_structure_traits< Type >::Div_mod div_mod;
    div_mod( x, y, q, r );
}

// others
template< class AS >
inline
typename Algebraic_structure_traits< AS >::Kth_root::result_type
kth_root( int k, const AS& x ) {
    typename Algebraic_structure_traits< AS >::Kth_root
        kth_root;
    return kth_root( k, x );
}


template< class Input_iterator >
inline
typename Algebraic_structure_traits< typename std::iterator_traits<Input_iterator>::value_type >
::Root_of::result_type
root_of( int k, Input_iterator begin, Input_iterator end ) {
    typedef typename std::iterator_traits<Input_iterator>::value_type AS;
    return typename Algebraic_structure_traits<AS>::Root_of()( k, begin, end );
}

// AST- and RET-functor adapting function
template< class Number_type >
inline
// select a Is_zero functor
typename boost::mpl::if_c<
 ::std::is_same< typename Algebraic_structure_traits< Number_type >::Is_zero,
 Null_functor  >::value ,
  typename Real_embeddable_traits< Number_type >::Is_zero,
  typename Algebraic_structure_traits< Number_type >::Is_zero
>::type::result_type
is_zero( const Number_type& x ) {
    // We take the Algebraic_structure_traits<>::Is_zero functor by default. If it
    //  is not available, we take the Real_embeddable_traits functor
    typename ::boost::mpl::if_c<
        ::std::is_same<
             typename Algebraic_structure_traits< Number_type >::Is_zero,
             Null_functor >::value ,
       typename Real_embeddable_traits< Number_type >::Is_zero,
       typename Algebraic_structure_traits< Number_type >::Is_zero >::type
       is_zero;
return is_zero( x );
}


template <class A, class B>
inline
typename Real_embeddable_traits< typename Coercion_traits<A,B>::Type >
::Compare::result_type
compare(const A& a, const B& b)
{
    typedef typename Coercion_traits<A,B>::Type Type;
    typename Real_embeddable_traits<Type>::Compare compare;
    return compare (a,b);
    // return (a < b) ? SMALLER : (b < a) ? LARGER : EQUAL;
}


// RET-Functor adapting functions
template< class Real_embeddable >
inline
//Real_embeddable
typename Real_embeddable_traits< Real_embeddable >::Abs::result_type
abs( const Real_embeddable& x ) {
    typename Real_embeddable_traits< Real_embeddable >::Abs abs;
    return abs( x );
}

template< class Real_embeddable >
inline
//::Sign
typename Real_embeddable_traits< Real_embeddable >::Sgn::result_type
sign( const Real_embeddable& x ) {
    typename Real_embeddable_traits< Real_embeddable >::Sgn sgn;
    return sgn( x );
}

template< class Real_embeddable >
inline
//bool
typename Real_embeddable_traits< Real_embeddable >::Is_finite::result_type
is_finite( const Real_embeddable& x ) {
    return typename Real_embeddable_traits< Real_embeddable >::Is_finite()( x );
}

template< class Real_embeddable >
inline
typename Real_embeddable_traits< Real_embeddable >::Is_positive::result_type
is_positive( const Real_embeddable& x ) {
    typename Real_embeddable_traits< Real_embeddable >::Is_positive
        is_positive;
    return is_positive( x );
}

template< class Real_embeddable >
inline
typename Real_embeddable_traits< Real_embeddable >::Is_negative::result_type
is_negative( const Real_embeddable& x ) {
    typename Real_embeddable_traits< Real_embeddable >::Is_negative
        is_negative;
    return is_negative( x );
}

/*
template< class Real_embeddable >
inline
typename Real_embeddable_traits< Real_embeddable >::Compare::result_type
//Comparison_result
compare( const Real_embeddable& x, const Real_embeddable& y ) {
    typename Real_embeddable_traits< Real_embeddable >::Compare compare;
    return compare( x, y );
}
*/

template< class Real_embeddable >
inline
typename Real_embeddable_traits< Real_embeddable >::To_double::result_type
//double
to_double( const Real_embeddable& x ) {
    typename Real_embeddable_traits< Real_embeddable >::To_double to_double;
    return to_double( x );
}

template< class Real_embeddable >
inline
typename Real_embeddable_traits< Real_embeddable >::To_interval::result_type
//std::pair< double, double >
to_interval( const Real_embeddable& x) {
    typename Real_embeddable_traits< Real_embeddable >::To_interval
        to_interval;
    return to_interval( x );
}

template <typename NT>
typename Coercion_traits<double, NT>::Type
approximate_sqrt(const NT& x, CGAL::Null_functor)
{
  return sqrt(CGAL::to_double(x));
}

template <typename NT, typename Sqrt>
typename Sqrt::result_type
approximate_sqrt(const NT& nt, Sqrt sqrt)
{
  return sqrt(nt);
}

template <typename NT>
decltype(auto) approximate_sqrt(const NT& nt)
{
  // the initial version of this function was using Algebraic_category
  // for the dispatch but some ring type (like Gmpz) provides a Sqrt
  // functor even if not being Field_with_sqrt.
  typedef CGAL::Algebraic_structure_traits<NT> AST;
  typedef typename AST::Sqrt Sqrt;
  return approximate_sqrt(nt, Sqrt());
}

template <class NT>
typename Same_uncertainty_nt<Comparison_result, NT>::type
compare_quotients(const NT& xnum, const NT& xden,
                  const NT& ynum, const NT& yden)
{
  // No assumptions on the sign of  den  are made

  // code assumes that SMALLER == - 1;
  CGAL_precondition( SMALLER == static_cast<Comparison_result>(-1) );

  int xsign =  sign(xnum) *  sign(xden) ;
  int ysign =  sign(ynum) *  sign(yden) ;
  if (xsign == 0) return static_cast<Comparison_result>(-ysign);
  if (ysign == 0) return static_cast<Comparison_result>(xsign);
  // now (x != 0) && (y != 0)
  int diff = xsign - ysign;
  if (diff == 0)
    {
      int msign =  sign(xden) *  sign(yden);
      NT leftop  = NT(xnum * yden * msign);
      NT rightop = NT(ynum * xden * msign);
      return  CGAL::compare(leftop, rightop);
    }
  else
    {
      return (xsign < ysign) ? SMALLER : LARGER;
    }
}


CGAL_NTS_END_NAMESPACE
} //namespace CGAL

#endif // CGAL_NUMBER_UTILS_H
