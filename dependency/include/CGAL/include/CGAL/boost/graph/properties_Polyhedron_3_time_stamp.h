// Copyright (c) 2017  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL: https://github.com/CGAL/cgal/blob/v5.6.2/Polyhedron/include/CGAL/boost/graph/properties_Polyhedron_3_time_stamp.h $
// $Id: properties_Polyhedron_3_time_stamp.h 968e9cf5a96 2022-06-16T13:29:20+02:00 Sébastien Loriot
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_PROPERTIES_POLYHEDRON_3_TIME_STAMP_H
#define CGAL_PROPERTIES_POLYHEDRON_3_TIME_STAMP_H

#include <CGAL/Polyhedron_3.h>

#define CGAL_HDS_PARAM_ template < class Traits, class Items, class Alloc> class HDS

namespace CGAL {

struct Polyhedron_face_time_stamp_pmap
{
  typedef void                               key_type;
  typedef std::size_t                        value_type;
  typedef std::size_t                        reference;
  typedef boost::read_write_property_map_tag category;
};

template <typename Handle_type>
std::size_t get(Polyhedron_face_time_stamp_pmap, Handle_type h)
{
  return h->time_stamp();
}

template <typename Handle_type>
void put(Polyhedron_face_time_stamp_pmap, Handle_type h, std::size_t ts)
{
  h->set_time_stamp(ts);
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
struct HDS_property_map<CGAL::Polyhedron_3<Gt, I, HDS, A>, CGAL::vertex_time_stamp_t>
{
  struct bind_
  {
    typedef Polyhedron_face_time_stamp_pmap type;
    typedef type const_type;
  };
};

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
struct HDS_property_map<CGAL::Polyhedron_3<Gt, I, HDS, A>, CGAL::halfedge_time_stamp_t>
  : public HDS_property_map<CGAL::Polyhedron_3<Gt, I, HDS, A>, CGAL::vertex_time_stamp_t>
{};

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
struct HDS_property_map<CGAL::Polyhedron_3<Gt, I, HDS, A>, CGAL::face_time_stamp_t>
  : public HDS_property_map<CGAL::Polyhedron_3<Gt, I, HDS, A>, CGAL::vertex_time_stamp_t>
{};





} // end namespace CGAL

#undef CGAL_HDS_PARAM_

#endif // CGAL_PROPERTIES_POLYHEDRON_3_TIME_STAMP_H
