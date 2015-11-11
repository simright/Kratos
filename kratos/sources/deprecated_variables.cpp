/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: JMCarbonell $
//   Date:                $Date:          2014 $
//   Revision:            $Revision:      1.20 $
//   Note:  IT HAS TO BE CLEANED AND SIMPLIFIED
//


// This define must be HERE
#define DKRATOS_EXPORT_INTERFACE_2 1

// System includes
#include <string>
#include <iostream>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/deprecated_variables.h"
#include "includes/kernel.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/constitutive_law.h"
#include "includes/geometrical_object.h"

#include "geometries/line_2d.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"
#include "geometries/point.h"
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"
#include "geometries/sphere_3d_1.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_9.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/prism_3d_6.h"
#include "geometries/prism_3d_15.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"

#include "python/add_deprecated_variables_to_python.h"

#include "includes/convection_diffusion_settings.h"
#include "includes/radiation_settings.h"

#include "includes/kratos_flags.h"

namespace Kratos
{

    //Create Variables by type:

    //bools


    //for Structural application:
    KRATOS_CREATE_VARIABLE( bool, IS_INACTIVE )

    //for Level Set application:
    KRATOS_CREATE_VARIABLE( bool, IS_DUPLICATED )
    KRATOS_CREATE_VARIABLE( bool, SPLIT_ELEMENT )
    KRATOS_CREATE_VARIABLE( bool, SPLIT_NODAL )

    KRATOS_CREATE_VARIABLE( int, IS_CONTACT_MASTER )
    KRATOS_CREATE_VARIABLE( int, IS_CONTACT_SLAVE )

    //for PFEM fluids application:
    KRATOS_CREATE_VARIABLE( int, IS_JACK_LINK )
    KRATOS_CREATE_VARIABLE( int, IMPOSED_PRESSURE )
    KRATOS_CREATE_VARIABLE( int, IMPOSED_VELOCITY_X )
    KRATOS_CREATE_VARIABLE( int, IMPOSED_VELOCITY_Y )
    KRATOS_CREATE_VARIABLE( int, IMPOSED_VELOCITY_Z )
    KRATOS_CREATE_VARIABLE( int, IMPOSED_ANGULAR_VELOCITY_X )
    KRATOS_CREATE_VARIABLE( int, IMPOSED_ANGULAR_VELOCITY_Y )
    KRATOS_CREATE_VARIABLE( int, IMPOSED_ANGULAR_VELOCITY_Z )
    
    //For the DEM Application:
    KRATOS_CREATE_VARIABLE(double, IMPOSED_VELOCITY_X_VALUE)
    KRATOS_CREATE_VARIABLE(double, IMPOSED_VELOCITY_Y_VALUE)
    KRATOS_CREATE_VARIABLE(double, IMPOSED_VELOCITY_Z_VALUE)
    KRATOS_CREATE_VARIABLE(double, IMPOSED_ANGULAR_VELOCITY_X_VALUE)
    KRATOS_CREATE_VARIABLE(double, IMPOSED_ANGULAR_VELOCITY_Y_VALUE)
    KRATOS_CREATE_VARIABLE(double, IMPOSED_ANGULAR_VELOCITY_Z_VALUE)

    KRATOS_CREATE_VARIABLE( double, IS_INLET )
    KRATOS_CREATE_VARIABLE( double, IS_INTERFACE )
    KRATOS_CREATE_VARIABLE( double, IS_VISITED )
    KRATOS_CREATE_VARIABLE( double, IS_EROSIONABLE )

    KRATOS_CREATE_VARIABLE( double, IS_STRUCTURE )
    KRATOS_CREATE_VARIABLE( double, IS_POROUS )
    KRATOS_CREATE_VARIABLE( double, IS_WATER )
    KRATOS_CREATE_VARIABLE( double, IS_FLUID )
    KRATOS_CREATE_VARIABLE( double, IS_BOUNDARY )
    KRATOS_CREATE_VARIABLE( double, IS_FREE_SURFACE )
    KRATOS_CREATE_VARIABLE( double, IS_AIR_EXIT )
    KRATOS_CREATE_VARIABLE( double, IS_LAGRANGIAN_INLET )
    KRATOS_CREATE_VARIABLE( double, IS_WATER_ELEMENT )


    KRATOS_CREATE_VARIABLE( double, IS_BURN )
    KRATOS_CREATE_VARIABLE( double, IS_DRIPPING )
    KRATOS_CREATE_VARIABLE( double, IS_PERMANENT )
    KRATOS_CREATE_VARIABLE( double, IS_WALL )

    KRATOS_CREATE_VARIABLE( double, Ypr ) //var name does not follow standard
    KRATOS_CREATE_VARIABLE( double, Yox )
    KRATOS_CREATE_VARIABLE( double, Yfuel )
    KRATOS_CREATE_VARIABLE( double, Hfuel )
    KRATOS_CREATE_VARIABLE( double, Hpr )
    KRATOS_CREATE_VARIABLE( double, Hpr1 )
    KRATOS_CREATE_VARIABLE( double, Hox )

    KRATOS_CREATE_VARIABLE( double, RADIATIVE_INTENSITY )
    KRATOS_CREATE_VARIABLE( double, RADIATIVE_INTENSITY_1 )
    KRATOS_CREATE_VARIABLE( double, RADIATIVE_INTENSITY_2 )
    KRATOS_CREATE_VARIABLE( double, RADIATIVE_INTENSITY_3 )
    KRATOS_CREATE_VARIABLE( double, RADIATIVE_INTENSITY_4 )
    KRATOS_CREATE_VARIABLE( double, RADIATIVE_INTENSITY_5 )
    KRATOS_CREATE_VARIABLE( double, RADIATIVE_INTENSITY_6 )
    KRATOS_CREATE_VARIABLE( double, RADIATIVE_INTENSITY_7 )
    KRATOS_CREATE_VARIABLE( double, RADIATIVE_INTENSITY_8 )
    KRATOS_CREATE_VARIABLE( double, RADIATIVE_INTENSITY_9 )
    KRATOS_CREATE_VARIABLE( double, RADIATIVE_INTENSITY_10 )
    KRATOS_CREATE_VARIABLE( double, RADIATIVE_INTENSITY_11 )
    KRATOS_CREATE_VARIABLE( double, RADIATIVE_INTENSITY_12 )
    KRATOS_CREATE_VARIABLE( double, RADIATIVE_INTENSITY_13 )
    KRATOS_CREATE_VARIABLE( double, RADIATIVE_INTENSITY_14 )
    KRATOS_CREATE_VARIABLE( double, RADIATIVE_INTENSITY_15 )
    KRATOS_CREATE_VARIABLE( double, RADIATIVE_INTENSITY_16 )
    KRATOS_CREATE_VARIABLE( double, RADIATIVE_INTENSITY_17 )
    KRATOS_CREATE_VARIABLE( double, RADIATIVE_INTENSITY_18 )
    KRATOS_CREATE_VARIABLE( double, RADIATIVE_INTENSITY_19 )
    KRATOS_CREATE_VARIABLE( double, RADIATIVE_INTENSITY_20 )
    KRATOS_CREATE_VARIABLE( double, RADIATIVE_INTENSITY_21 )
    KRATOS_CREATE_VARIABLE( double, RADIATIVE_INTENSITY_22 )
    KRATOS_CREATE_VARIABLE( double, RADIATIVE_INTENSITY_23 )
    KRATOS_CREATE_VARIABLE( double, RADIATIVE_INTENSITY_24 )

    KRATOS_CREATE_VARIABLE( double, rhoD )
    KRATOS_CREATE_VARIABLE( double, xi )
    KRATOS_CREATE_VARIABLE( double, a )
    KRATOS_CREATE_VARIABLE( double, b )


    KRATOS_CREATE_VARIABLE( double, IS_SLIP )

    //for Level Set application:
    KRATOS_CREATE_VARIABLE( double, IS_DIVIDED )

    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( xi_c )

  

  void KratosApplication::RegisterDeprecatedVariables()
  {
      
    KRATOS_REGISTER_VARIABLE(  IS_INACTIVE )

    //for Level Set application:
    KRATOS_REGISTER_VARIABLE(  IS_DUPLICATED )
    KRATOS_REGISTER_VARIABLE(  SPLIT_ELEMENT )
    KRATOS_REGISTER_VARIABLE(  SPLIT_NODAL )

    KRATOS_REGISTER_VARIABLE( IS_CONTACT_MASTER )
    KRATOS_REGISTER_VARIABLE( IS_CONTACT_SLAVE )
    //for PFEM fluids application:
    KRATOS_REGISTER_VARIABLE( IS_JACK_LINK )
    KRATOS_REGISTER_VARIABLE( IMPOSED_PRESSURE )
    KRATOS_REGISTER_VARIABLE( IMPOSED_VELOCITY_X )
    KRATOS_REGISTER_VARIABLE( IMPOSED_VELOCITY_Y )
    KRATOS_REGISTER_VARIABLE( IMPOSED_VELOCITY_Z )
    KRATOS_REGISTER_VARIABLE( IMPOSED_ANGULAR_VELOCITY_X )
    KRATOS_REGISTER_VARIABLE( IMPOSED_ANGULAR_VELOCITY_Y )
    KRATOS_REGISTER_VARIABLE( IMPOSED_ANGULAR_VELOCITY_Z )
    
    //For the DEM Application:
    KRATOS_REGISTER_VARIABLE(IMPOSED_VELOCITY_X_VALUE)
    KRATOS_REGISTER_VARIABLE(IMPOSED_VELOCITY_Y_VALUE)
    KRATOS_REGISTER_VARIABLE(IMPOSED_VELOCITY_Z_VALUE)
    KRATOS_REGISTER_VARIABLE(IMPOSED_ANGULAR_VELOCITY_X_VALUE)
    KRATOS_REGISTER_VARIABLE(IMPOSED_ANGULAR_VELOCITY_Y_VALUE)
    KRATOS_REGISTER_VARIABLE(IMPOSED_ANGULAR_VELOCITY_Z_VALUE)

    KRATOS_REGISTER_VARIABLE(   IS_INLET )
    KRATOS_REGISTER_VARIABLE(   IS_INTERFACE )
    KRATOS_REGISTER_VARIABLE(   IS_VISITED )
    KRATOS_REGISTER_VARIABLE( IS_EROSIONABLE )

    KRATOS_REGISTER_VARIABLE(   IS_STRUCTURE )
    KRATOS_REGISTER_VARIABLE(   IS_POROUS )
    KRATOS_REGISTER_VARIABLE(   IS_WATER )
    KRATOS_REGISTER_VARIABLE(   IS_FLUID )
    KRATOS_REGISTER_VARIABLE(   IS_BOUNDARY )
    KRATOS_REGISTER_VARIABLE(   IS_FREE_SURFACE )
    KRATOS_REGISTER_VARIABLE(   IS_AIR_EXIT )
    KRATOS_REGISTER_VARIABLE(   IS_LAGRANGIAN_INLET )
    KRATOS_REGISTER_VARIABLE(   IS_WATER_ELEMENT )


    KRATOS_REGISTER_VARIABLE(   IS_BURN )
    KRATOS_REGISTER_VARIABLE(   IS_DRIPPING )
    KRATOS_REGISTER_VARIABLE(   IS_PERMANENT )
    KRATOS_REGISTER_VARIABLE(   IS_WALL )

    KRATOS_REGISTER_VARIABLE(   Ypr ) //var name does not follow standard
    KRATOS_REGISTER_VARIABLE(   Yox )
    KRATOS_REGISTER_VARIABLE(   Yfuel )
    KRATOS_REGISTER_VARIABLE(   Hfuel )
    KRATOS_REGISTER_VARIABLE(   Hpr )
    KRATOS_REGISTER_VARIABLE(   Hpr1 )
    KRATOS_REGISTER_VARIABLE(   Hox )

    KRATOS_REGISTER_VARIABLE(   RADIATIVE_INTENSITY )
    KRATOS_REGISTER_VARIABLE(   RADIATIVE_INTENSITY_1 )
    KRATOS_REGISTER_VARIABLE(   RADIATIVE_INTENSITY_2 )
    KRATOS_REGISTER_VARIABLE(   RADIATIVE_INTENSITY_3 )
    KRATOS_REGISTER_VARIABLE(   RADIATIVE_INTENSITY_4 )
    KRATOS_REGISTER_VARIABLE(   RADIATIVE_INTENSITY_5 )
    KRATOS_REGISTER_VARIABLE(   RADIATIVE_INTENSITY_6 )
    KRATOS_REGISTER_VARIABLE(   RADIATIVE_INTENSITY_7 )
    KRATOS_REGISTER_VARIABLE(   RADIATIVE_INTENSITY_8 )
    KRATOS_REGISTER_VARIABLE(   RADIATIVE_INTENSITY_9 )
    KRATOS_REGISTER_VARIABLE(   RADIATIVE_INTENSITY_10 )
    KRATOS_REGISTER_VARIABLE(   RADIATIVE_INTENSITY_11 )
    KRATOS_REGISTER_VARIABLE(   RADIATIVE_INTENSITY_12 )
    KRATOS_REGISTER_VARIABLE(   RADIATIVE_INTENSITY_13 )
    KRATOS_REGISTER_VARIABLE(   RADIATIVE_INTENSITY_14 )
    KRATOS_REGISTER_VARIABLE(   RADIATIVE_INTENSITY_15 )
    KRATOS_REGISTER_VARIABLE(   RADIATIVE_INTENSITY_16 )
    KRATOS_REGISTER_VARIABLE(   RADIATIVE_INTENSITY_17 )
    KRATOS_REGISTER_VARIABLE(   RADIATIVE_INTENSITY_18 )
    KRATOS_REGISTER_VARIABLE(   RADIATIVE_INTENSITY_19 )
    KRATOS_REGISTER_VARIABLE(   RADIATIVE_INTENSITY_20 )
    KRATOS_REGISTER_VARIABLE(   RADIATIVE_INTENSITY_21 )
    KRATOS_REGISTER_VARIABLE(   RADIATIVE_INTENSITY_22 )
    KRATOS_REGISTER_VARIABLE(   RADIATIVE_INTENSITY_23 )
    KRATOS_REGISTER_VARIABLE(   RADIATIVE_INTENSITY_24 )

    KRATOS_REGISTER_VARIABLE(   rhoD )
    KRATOS_REGISTER_VARIABLE(   xi )
    KRATOS_REGISTER_VARIABLE(   a )
    KRATOS_REGISTER_VARIABLE(   b )


    KRATOS_REGISTER_VARIABLE(   IS_SLIP )

    //for Level Set application:
    KRATOS_REGISTER_VARIABLE(   IS_DIVIDED )

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( xi_c )


  }


}  // namespace Kratos.

// This define must be HERE
#undef DKRATOS_EXPORT_INTERFACE_2


