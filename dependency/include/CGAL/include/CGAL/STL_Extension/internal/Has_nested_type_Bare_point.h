// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// Copyright (c) 2016 GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL: https://github.com/CGAL/cgal/blob/v5.6.2/STL_Extension/include/CGAL/STL_Extension/internal/Has_nested_type_Bare_point.h $
// $Id: Has_nested_type_Bare_point.h 98e471849bd 2021-08-26T11:33:39+02:00 Sébastien Loriot
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Jane Tournois

#ifndef CGAL_HAS_NESTED_TYPE_BARE_POINT_H
#define CGAL_HAS_NESTED_TYPE_BARE_POINT_H

#include <boost/mpl/has_xxx.hpp>

namespace CGAL {

  namespace internal {

    BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_type_Bare_point, Bare_point, false)

    template<typename Gt>
    struct Bare_point_type
    {
    public:
      typedef typename Gt::Bare_point type;
    };

  } // end namespace internal
} // end namespace CGAL

#endif // CGAL_HAS_NESTED_TYPE_BARE_POINT_H
