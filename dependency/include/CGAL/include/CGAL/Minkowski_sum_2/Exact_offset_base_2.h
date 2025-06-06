// Copyright (c) 2006  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL: https://github.com/CGAL/cgal/blob/v5.6.2/Minkowski_sum_2/include/CGAL/Minkowski_sum_2/Exact_offset_base_2.h $
// $Id: Exact_offset_base_2.h 488ba8c4e48 2022-08-10T23:32:42+03:00 Efi Fogel
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Ron Wein   <wein_r@yahoo.com>

#ifndef CGAL_EXACT_OFFSET_BASE_H
#define CGAL_EXACT_OFFSET_BASE_H

#include <CGAL/license/Minkowski_sum_2.h>


#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Gps_traits_2.h>
#include <CGAL/Minkowski_sum_2/Labels.h>
#include <CGAL/Minkowski_sum_2/Arr_labeled_traits_2.h>
#include <CGAL/use.h>

namespace CGAL {

/*! \class
 * A base class for computing the offset of a given polygon by a given
 * radius in an exact manner.
 */
template <typename Traits_, typename Container_>
class Exact_offset_base_2 {
private:
  typedef Traits_                                        Traits_2;

  // Rational kernel types:
  typedef typename Traits_2::Rat_kernel                  Rat_kernel;
  typedef typename Rat_kernel::FT                        Rational;
  typedef typename Rat_kernel::Point_2                   Rat_point_2;
  typedef typename Rat_kernel::Line_2                    Rat_line_2;
  typedef typename Rat_kernel::Circle_2                  Rat_circle_2;

protected:

  typedef Rat_kernel                                     Basic_kernel;
  typedef Rational                                       Basic_NT;

private:
  // Algebraic kernel types:
  typedef typename Traits_2::Alg_kernel                  Alg_kernel;
  typedef typename Alg_kernel::FT                        Algebraic;
  typedef typename Alg_kernel::Point_2                   Alg_point_2;

  typedef typename Traits_2::Nt_traits                   Nt_traits;

  // Traits-class types:
  typedef typename Traits_2::Curve_2                     Curve_2;
  typedef typename Traits_2::X_monotone_curve_2          X_monotone_curve_2;

  typedef CGAL::Gps_traits_2<Traits_2>                   Gps_traits_2;

protected:
  typedef CGAL::Polygon_2<Rat_kernel, Container_>        Polygon_2;
  typedef CGAL::Polygon_with_holes_2<Rat_kernel,
                                     Container_>         Polygon_with_holes_2;
  typedef typename Gps_traits_2::Polygon_2               Offset_polygon_2;

private:

  // Polygon-related types:
  typedef typename Polygon_2::Vertex_circulator          Vertex_circulator;

protected:
  typedef Arr_labeled_traits_2<Traits_2>                 Labeled_traits_2;
  typedef typename Labeled_traits_2::X_monotone_curve_2  Labeled_curve_2;

public:
  /*! Default constructor. */
  Exact_offset_base_2() {}

protected:
  /*! Compute the curves that constitute the offset of a simple polygon by a
   * given radius.
   * \param pgn The polygon.
   * \param orient The orientation to traverse the vertices.
   * \param r The offset radius.
   * \param cycle_id The index of the cycle.
   * \param oi An output iterator for the offset curves.
   * \pre The value type of the output iterator is Labeled_curve_2.
   * \return A past-the-end iterator for the holes container.
   */
  template <typename OutputIterator>
  OutputIterator _offset_polygon(const Polygon_2& pgn,
                                 CGAL::Orientation orient,
                                 const Rational& r,
                                 unsigned int cycle_id,
                                 OutputIterator oi) const
  {
    // Prepare circulators over the polygon vertices.
    const bool forward = (pgn.orientation() == orient);
    Vertex_circulator first, curr, next;

    first = pgn.vertices_circulator();
    curr = first;
    next = first;

    // Traverse the polygon vertices and edges and construct the arcs that
    // constitute the single convolution cycle.
    const Rational sqr_r = CGAL::square (r);
    Rational x1, y1;                    // The source of the current edge.
    Rational x2, y2;                    // The target of the current edge.
    Rational delta_x, delta_y;          // (x2 - x1) and (y2 - y1), resp.
    Algebraic len;                      // The length of the current edge.
    Algebraic trans_x, trans_y;         // The translation vector.
    Alg_point_2 op1, op2;               // The edge points of the offset edge.
    Alg_point_2 first_op;               // The first offset point.
    Algebraic a, b, c;

    unsigned int curve_index(0);
    std::list<Object> xobjs;

    Traits_2 traits;
    auto nt_traits = traits.nt_traits();
    const Algebraic alg_r = nt_traits->convert(r);
    auto f_make_x_monotone = traits.make_x_monotone_2_object();

    auto alg_ker = traits.alg_kernel();
    auto f_equal = alg_ker->equal_2_object();

    bool assign_success;

    do {
      // Get a circulator for the next vertex (in the proper orientation).
      if (forward) ++next;
      else --next;

      // Compute the vector v = (delta_x, delta_y) of the current edge,
      // and compute the edge length ||v||.
      x1 = curr->x();
      y1 = curr->y();
      x2 = next->x();
      y2 = next->y();

      delta_x = x2 - x1;
      delta_y = y2 - y1;
      len = nt_traits->sqrt(nt_traits->convert(CGAL::square(delta_x) +
                                               CGAL::square(delta_y)));

      // The angle theta between the vector v and the x-axis is given by:
      //
      //                 y2 - y1                        x2 - x1
      //   sin(alpha) = ---------         cos(alpha) = ---------
      //                  ||v||                          ||v||
      //
      // To offset the endpoints of the current edge we compute a vector
      // (trans_x, trans_y) perpendicular to v. Since we traverse the polygon
      // in a counterclockwise manner, the angle this vector forms with the
      // x-axis is (alpha - PI/2), and we have:
      //
      //   trans_x = r*cos(alpha - PI/2) = r*sin(alpha)
      //   trans_y = r*sin(alpha - PI/2) = -r*cos(alpha)
      trans_x = nt_traits->convert(r * delta_y) / len;
      trans_y = nt_traits->convert(-r * delta_x) / len;

      // Construct the first offset vertex, which corresponds to the
      // source vertex of the current polygon edge.
      op1 = Alg_point_2(nt_traits->convert(x1) + trans_x,
                        nt_traits->convert(y1) + trans_y);

      if (curr == first) {
        // This is the first edge we visit -- store op1 for future use.
        first_op = op1;
      }
      else {
        if (! f_equal (op2, op1)) {
          // Connect op2 (from the previous iteration) and op1 with a circular
          // arc, whose supporting circle is (x1, x2) with radius r.
          auto ctr_cv = traits.construct_curve_2_object();
          Curve_2 arc = ctr_cv(Rat_circle_2 (*curr, sqr_r),
                               CGAL::COUNTERCLOCKWISE, op2, op1);

          // Subdivide the arc into x-monotone subarcs and append them to the
          // convolution cycle.
          xobjs.clear();
          f_make_x_monotone(arc, std::back_inserter(xobjs));

          for (auto xobj_it = xobjs.begin(); xobj_it != xobjs.end(); ++xobj_it) {
            X_monotone_curve_2 xarc;
            assign_success = CGAL::assign(xarc, *xobj_it);
            CGAL_assertion (assign_success);
            CGAL_USE(assign_success);

            *oi++ = Labeled_curve_2(xarc, X_curve_label(xarc.is_directed_right(),
                                                        cycle_id, curve_index));
            curve_index++;
          }
        }
      }

      // Construct the second offset vertex, which corresponds to the
      // target vertex of the current polygon edge.
      op2 = Alg_point_2(nt_traits->convert(x2) + trans_x,
                        nt_traits->convert(y2) + trans_y);

      // The equation of the line connecting op1 and op2 is given by:
      //
      //   (y1 - y2)*x + (x2 - x1)*y + (r*len - y1*x2 - x1*y2) = 0
      //
      a = nt_traits->convert(-delta_y);
      b = nt_traits->convert(delta_x);
      c = alg_r*len - nt_traits->convert(y1*x2 - x1*y2);

      auto ctr_xcv = traits.construct_x_monotone_curve_2_object();
      X_monotone_curve_2 xarc = ctr_xcv(a, b, c, op1, op2);
      *oi++ = Labeled_curve_2(xarc, X_curve_label(xarc.is_directed_right(),
                                                  cycle_id, curve_index));
      curve_index++;

      // Proceed to the next polygon vertex.
      curr = next;

    } while (curr != first);

    if (! f_equal (op2, first_op)) {
      // Close the convolution cycle by creating the final circular arc,
      // centered at the first vertex.
      auto ctr_cv = traits.construct_curve_2_object();
      Curve_2 arc = ctr_cv(Rat_circle_2 (*first, sqr_r),
                           CGAL::COUNTERCLOCKWISE, op2, first_op);

      // Subdivide the arc into x-monotone subarcs and append them to the
      // convolution cycle.
      bool is_last;

      xobjs.clear();
      f_make_x_monotone(arc, std::back_inserter(xobjs));

      auto xobj_it = xobjs.begin();
      while (xobj_it != xobjs.end()) {
        X_monotone_curve_2 xarc;
        assign_success = CGAL::assign(xarc, *xobj_it);
        CGAL_assertion (assign_success);
        CGAL_USE(assign_success);

        ++xobj_it;
        is_last = (xobj_it == xobjs.end());

        *oi++ = Labeled_curve_2(xarc, X_curve_label(xarc.is_directed_right(),
                                                    cycle_id, curve_index,
                                                    is_last));
        curve_index++;
      }
    }

    return oi;
  }

};


} //namespace CGAL

#endif
