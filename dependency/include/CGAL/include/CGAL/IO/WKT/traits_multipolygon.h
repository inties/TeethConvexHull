// Copyright (c) 2018  GeometryFactory Sarl (France).
// All rights reserved.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL: https://github.com/CGAL/cgal/blob/v5.6.2/Stream_support/include/CGAL/IO/WKT/traits_multipolygon.h $
// $Id: traits_multipolygon.h 98e471849bd 2021-08-26T11:33:39+02:00 Sébastien Loriot
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Maxime Gimeno

#ifndef CGAL_IO_WKT_TRAITS_MULTIPOLYGON_H
#define CGAL_IO_WKT_TRAITS_MULTIPOLYGON_H

#include <CGAL/Stream_support/internal/Geometry_container.h>

#include <boost/geometry/io/wkt/write.hpp>
#include <boost/geometry/io/wkt/read.hpp>

namespace boost {
namespace geometry {
namespace traits {

// WKT traits for MultiPolygon
template< typename R >
struct tag<CGAL::internal::Geometry_container<R, multi_polygon_tag> >
{
  typedef multi_polygon_tag type;
};

} // namespace traits
} // namespace geometry
} // namespace boost

#endif // CGAL_IO_WKT_TRAITS_MULTIPOLYGON_H

