// Copyright (c) 2016  GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL: https://github.com/CGAL/cgal/blob/v5.6.2/Installation/include/CGAL/license/Surface_mesh_approximation.h $
// $Id: Surface_mesh_approximation.h c1afb483f58 2022-07-19T09:04:19+02:00 Sébastien Loriot
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Andreas Fabri
//
// Warning: this file is generated, see include/CGAL/license/README.md

#ifndef CGAL_LICENSE_SURFACE_MESH_APPROXIMATION_H
#define CGAL_LICENSE_SURFACE_MESH_APPROXIMATION_H

#include <CGAL/config.h>
#include <CGAL/license.h>

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_COMMERCIAL_LICENSE

#  if CGAL_SURFACE_MESH_APPROXIMATION_COMMERCIAL_LICENSE < CGAL_RELEASE_DATE

#    if defined(CGAL_LICENSE_WARNING)

       CGAL_pragma_warning("Your commercial license for CGAL does not cover "
                           "this release of the Triangulated Surface Mesh Approximation package.")
#    endif

#    ifdef CGAL_LICENSE_ERROR
#      error "Your commercial license for CGAL does not cover this release \
              of the Triangulated Surface Mesh Approximation package. \
              You get this error, as you defined CGAL_LICENSE_ERROR."
#    endif // CGAL_LICENSE_ERROR

#  endif // CGAL_SURFACE_MESH_APPROXIMATION_COMMERCIAL_LICENSE < CGAL_RELEASE_DATE

#else // no CGAL_SURFACE_MESH_APPROXIMATION_COMMERCIAL_LICENSE

#  if defined(CGAL_LICENSE_WARNING)
     CGAL_pragma_warning("\nThe macro CGAL_SURFACE_MESH_APPROXIMATION_COMMERCIAL_LICENSE is not defined."
                          "\nYou use the CGAL Triangulated Surface Mesh Approximation package under "
                          "the terms of the GPLv3+.")
#  endif // CGAL_LICENSE_WARNING

#  ifdef CGAL_LICENSE_ERROR
#    error "The macro CGAL_SURFACE_MESH_APPROXIMATION_COMMERCIAL_LICENSE is not defined.\
            You use the CGAL Triangulated Surface Mesh Approximation package under the terms of \
            the GPLv3+. You get this error, as you defined CGAL_LICENSE_ERROR."
#  endif // CGAL_LICENSE_ERROR

#endif // no CGAL_SURFACE_MESH_APPROXIMATION_COMMERCIAL_LICENSE

#endif // CGAL_LICENSE_SURFACE_MESH_APPROXIMATION_H
