// Copyright (c) 1997-2001
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL: https://github.com/CGAL/cgal/blob/v5.6.2/Optimisation_basic/include/CGAL/Optimisation/Access_coordinates_begin_2.h $
// $Id: Access_coordinates_begin_2.h 07793738355 2020-03-26T13:31:46+01:00 Sébastien Loriot
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Schoenherr <sven@inf.ethz.ch>

#ifndef CGAL_OPTIMISATION_ACCESS_COORDINATES_BEGIN_2_H
#define CGAL_OPTIMISATION_ACCESS_COORDINATES_BEGIN_2_H

#include <cstddef>
#include <iterator>

namespace CGAL {

// Class declarations
// ==================
template < class R_ >
class Access_coordinates_begin_2;

template < class R_ >
class Point_2_coordinate_iterator;

// Class interfaces
// ================
template < class R_ >
class Point_2_coordinate_iterator {
  public:
    // self
    typedef  R_                         R;
    typedef  Point_2_coordinate_iterator<R>
                                        Self;

    // types
    typedef  typename R::Point_2        Point;

    // iterator types
    typedef  typename R::RT             value_type;
    typedef  std::ptrdiff_t             difference_type;
    typedef  value_type*                pointer;
    typedef  value_type&                reference;
    typedef  std::random_access_iterator_tag
                                        iterator_category;

    // forward operations
    Point_2_coordinate_iterator( const Point&  point = Point(),
                                 int           index = 0)
        : p( point), i( index) { }

    bool        operator == ( const Self& it) const { return ( i == it.i);}
    bool        operator != ( const Self& it) const { return ( i != it.i);}

    value_type  operator *  ( ) const { return p.homogeneous( i); }

    Self&       operator ++ (    ) {                   ++i; return *this; }
    Self        operator ++ ( int) { Self tmp = *this; ++i; return tmp;   }

    // bidirectional operations
    Self&       operator -- (    ) {                   --i; return *this; }
    Self        operator -- ( int) { Self tmp = *this; --i; return tmp;   }

    // random access operations
    Self&       operator += ( int n) { i += n; return *this; }
    Self&       operator -= ( int n) { i -= n; return *this; }

    Self        operator +  ( int n) const
                                     { Self tmp = *this; return tmp += n; }
    Self        operator -  ( int n) const
                                     { Self tmp = *this; return tmp -= n; }

    difference_type
                operator -  ( const Self& it) const { return i - it.i; }

    value_type  operator [] ( int n) const { return p.homogeneous( i+n); }

    bool   operator <  ( const Self& it) const { return ( i <  it.i); }
    bool   operator >  ( const Self& it) const { return ( i >  it.i); }
    bool   operator <= ( const Self& it) const { return ( i <= it.i); }
    bool   operator >= ( const Self& it) const { return ( i >= it.i); }

private:
    Point         p;
    int           i;
};

template < class R_ >
class Access_coordinates_begin_2 {
  public:
    // self
    typedef  R_                         R;
    typedef  Access_coordinates_begin_2<R>
                                        Self;

    // types
    typedef  typename R::Point_2        Point;
    typedef  Point_2_coordinate_iterator<R>
                                        Coordinate_iterator;

    // unary function class types
    typedef  Coordinate_iterator        result_type;
    typedef  Point                      argument_type;

    // creation
    Access_coordinates_begin_2( ) { }

    // operations
    Coordinate_iterator
    operator() ( const Point& p) const { return Coordinate_iterator( p); }
};

} //namespace CGAL

#endif // CGAL_OPTIMISATION_ACCESS_COORDINATES_BEGIN_2_H

// ===== EOF ==================================================================
