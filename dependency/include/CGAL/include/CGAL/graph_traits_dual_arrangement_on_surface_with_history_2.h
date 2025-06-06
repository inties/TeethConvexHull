// Copyright (c) 2018 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL: https://github.com/CGAL/cgal/blob/v5.6.2/Arrangement_on_surface_2/include/CGAL/graph_traits_dual_arrangement_on_surface_with_history_2.h $
// $Id: graph_traits_dual_arrangement_on_surface_with_history_2.h 6d3176e0619 2022-01-07T14:42:25+01:00 Sébastien Loriot
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Ron Wein         <wein@post.tau.ac.il>
//             Ophir Setter     <ophirset@post.tau.ac.il>
//             Sebastien Loriot <sebastien.loriot@cgal.org>
//             Efi Fogel        <efifogel@gmail.com>

#ifndef CGAL_GRAPH_TRAITS_DUAL_ARRANGEMENT_ON_SURFACE_WITH_HISTORY_2_H
#define CGAL_GRAPH_TRAITS_DUAL_ARRANGEMENT_ON_SURFACE_WITH_HISTORY_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 * Definition of:
 * 1. the specialized Dual<Arrangement_on_surface_with_history_2> class,
 * 2. the specialized
 *      boost::graph_traits<Dual<Arrangement_on_surface_with_history_2> >class,
 * 3. The free functions required by the various graph concepts.
 */

// include this to avoid a VC15 warning
#include <CGAL/Named_function_parameters.h>

#include <CGAL/Arrangement_on_surface_with_history_2.h>
#include <CGAL/Arrangement_2/graph_traits_dual.h>
#include <CGAL/disable_warnings.h>

namespace CGAL {

// The specialized Dual<Arrangement_on_surface_with_history_2... class template.
template <typename GeomTraits_2, typename TopolTraits>
class Dual<Arrangement_on_surface_with_history_2<GeomTraits_2, TopolTraits> > :
    public Dual_arrangement_on_surface<Arrangement_on_surface_with_history_2
                                       <GeomTraits_2, TopolTraits> >
{
public:
  typedef Arrangement_on_surface_with_history_2<GeomTraits_2, TopolTraits>                Arrangement;
  typedef typename Arrangement::Geometry_traits_2         Geometry_traits_2;
  typedef typename Arrangement::Topology_traits           Topology_traits;

private:
  typedef Dual_arrangement_on_surface<Arrangement>        Base;

public:
  /*! Default constructor. */
  Dual() : Base() {}

  /*! Constructor from an arrangement. */
  Dual(const Arrangement& arr) : Base(arr) {}
};

}

namespace boost {

// The specialized
// graph_traits<CGAL::Dual<CGAL::Arrangement_on_surface_with_history_2... class
// template.
template <typename GeomTraits_2, typename TopolTraits>
class graph_traits<CGAL::Dual<CGAL::Arrangement_on_surface_with_history_2
                              <GeomTraits_2, TopolTraits> > > :
    public CGAL::Graph_traits_dual_arr_on_surface_impl
             <CGAL::Arrangement_on_surface_with_history_2<GeomTraits_2,
                                                          TopolTraits> >
{};

}

namespace CGAL {

// Templates of free functions that handle
//   graph_traits<Dual<Arrangement_on_surface_with_history_2...
// class template.
CGAL_DUAL_ARRANGEMENT_2_OUT_DEGREE(Arrangement_on_surface_with_history_2)
CGAL_DUAL_ARRANGEMENT_2_OUT_EDGES(Arrangement_on_surface_with_history_2)
CGAL_DUAL_ARRANGEMENT_2_SOURCE(Arrangement_on_surface_with_history_2)
CGAL_DUAL_ARRANGEMENT_2_TARGET(Arrangement_on_surface_with_history_2)
CGAL_DUAL_ARRANGEMENT_2_IN_DEGREE(Arrangement_on_surface_with_history_2)
CGAL_DUAL_ARRANGEMENT_2_IN_EDGES(Arrangement_on_surface_with_history_2)
CGAL_DUAL_ARRANGEMENT_2_DEGREE(Arrangement_on_surface_with_history_2)
CGAL_DUAL_ARRANGEMENT_2_NUM_VERTICES(Arrangement_on_surface_with_history_2)
CGAL_DUAL_ARRANGEMENT_2_VERTICES(Arrangement_on_surface_with_history_2)
CGAL_DUAL_ARRANGEMENT_2_NUM_EDGES(Arrangement_on_surface_with_history_2)
CGAL_DUAL_ARRANGEMENT_2_EDGES(Arrangement_on_surface_with_history_2)

}

#endif
