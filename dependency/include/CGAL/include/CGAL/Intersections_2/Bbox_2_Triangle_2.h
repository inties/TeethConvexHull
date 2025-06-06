// Copyright (c) 2019
// GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL: https://github.com/CGAL/cgal/blob/v5.6.2/Intersections_2/include/CGAL/Intersections_2/Bbox_2_Triangle_2.h $
// $Id: Bbox_2_Triangle_2.h 386c6a3ac26 2022-11-22T18:42:13+01:00 Mael
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Maxime Gimeno


#ifndef CGAL_INTERSECTIONS_BBOX_2_TRIANGLE_2_H
#define CGAL_INTERSECTIONS_BBOX_2_TRIANGLE_2_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Intersections_2/Iso_rectangle_2_Triangle_2.h>

namespace CGAL {

template <class K>
inline
typename K::Boolean
do_intersect(const Triangle_2<K>& tr,
             const Bbox_2& box)
{
  typename K::Iso_rectangle_2 rec(box.xmin(), box.ymin(), box.xmax(), box.ymax());
  return do_intersect(rec, tr);
}

template <class K>
inline
typename K::Boolean
do_intersect(const Bbox_2& box,
             const Triangle_2<K>& tr)
{
  return do_intersect(tr, box);
}

template<typename K>
typename Intersection_traits<K, typename K::Triangle_2, Bbox_2>::result_type
intersection(const Bbox_2& box,
             const Triangle_2<K>& tr)
{
  typename K::Iso_rectangle_2 rec(box.xmin(), box.ymin(), box.xmax(), box.ymax());
  return intersection(rec, tr);
}

template<typename K>
typename Intersection_traits<K, typename K::Triangle_2, Bbox_2>::result_type
intersection(const Triangle_2<K>& tr,
             const Bbox_2& box)
{
  return intersection(box, tr);
}

} // namespace CGAL

#endif // CGAL_INTERSECTIONS_BBOX_2_TRIANGLE_2_H
