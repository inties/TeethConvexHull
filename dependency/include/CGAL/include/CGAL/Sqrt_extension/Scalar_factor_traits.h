// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL: https://github.com/CGAL/cgal/blob/v5.6.2/Number_types/include/CGAL/Sqrt_extension/Scalar_factor_traits.h $
// $Id: Scalar_factor_traits.h 52164b1fbae 2019-10-19T15:34:59+02:00 Sébastien Loriot
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hemmer   <hemmer@mpi-inf.mpg.de>


#ifndef CGAL_SQRT_EXTENSION_SCALAR_FACTOR_TRAITS_H
#define CGAL_SQRT_EXTENSION_SCALAR_FACTOR_TRAITS_H

#include <CGAL/basic.h>

namespace CGAL {

// This is the specialization for Sqrt_extension
template <class COEFF, class ROOT, class ACDE_TAG,class FP_TAG>
class Scalar_factor_traits< Sqrt_extension<COEFF, ROOT, ACDE_TAG,FP_TAG> > {
public:

    //! the number type for which this instance has been instantiated
    typedef Sqrt_extension<COEFF, ROOT,ACDE_TAG,FP_TAG> NT;
      //! the number type of scalars that can be extracted from NT
    typedef typename Scalar_factor_traits<COEFF>::Scalar Scalar;

    class Scalar_factor
    {
    public:
        //! argument type
        typedef NT argument_type;
        //! first argument type
        typedef NT first_argument_type;
        //! second argument type
        typedef Scalar second_argument_type;
        //! result type
        typedef Scalar result_type;

        Scalar
        operator () (const NT& x, const Scalar& d_ = Scalar(0) ) {
            typename Scalar_factor_traits<COEFF>::Scalar_factor sfac;

            Scalar d(d_);
            Scalar unity(1);
            if(d==unity) return d;
            d=sfac(x.a0(),d);
            if(d==unity) return d;
            if(x.is_extended())
                d=sfac(x.a1(),d);
            return d;
        }
    };

    class Scalar_div
    {
    public:
        //! first_argument_type
        typedef NT first_argument_type;
        //! second_argument_type
        typedef Scalar second_argument_type;
        //! divides an extension \c a by a scalar factor \c b
        void operator () (NT& a, const Scalar& b) {
            CGAL_precondition(b != Scalar(0));
            typename Scalar_factor_traits<COEFF>::Scalar_div sdiv;
            sdiv(a.a0(), b); sdiv(a.a1(), b); // perform division in place
        }
    };
};

} //namespace CGAL

#endif
