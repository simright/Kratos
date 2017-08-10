// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
typedef array_1d<double,3> Vector3;

    // Generalized eigenvalue problem
    KRATOS_CREATE_VARIABLE( int, BUILD_LEVEL )
    KRATOS_CREATE_VARIABLE( Vector, EIGENVALUE_VECTOR)
    KRATOS_CREATE_VARIABLE( Matrix , EIGENVECTOR_MATRIX )

    // Geometrical
    KRATOS_CREATE_VARIABLE( double, AREA )
    KRATOS_CREATE_VARIABLE( double, IT )
    KRATOS_CREATE_VARIABLE( double, IY )
    KRATOS_CREATE_VARIABLE( double, IZ )
    KRATOS_CREATE_VARIABLE( double, CROSS_AREA )
    KRATOS_CREATE_VARIABLE( double, MEAN_RADIUS )
    KRATOS_CREATE_VARIABLE( int,    SECTION_SIDES )
    KRATOS_CREATE_VARIABLE( Matrix , GEOMETRIC_STIFFNESS )

    //fusseder
    KRATOS_CREATE_VARIABLE( double, IZ_SENSITIVITY )

    // Truss generalized variables
    KRATOS_CREATE_VARIABLE(double, TRUSS_PRESTRESS_PK2)
    KRATOS_CREATE_VARIABLE(bool, TRUSS_IS_CABLE)

    // Beam generalized variables
    KRATOS_CREATE_VARIABLE(double, AREA_EFFECTIVE_Y)
    KRATOS_CREATE_VARIABLE(double, AREA_EFFECTIVE_Z)
    KRATOS_CREATE_VARIABLE(double, INERTIA_ROT_Y)
    KRATOS_CREATE_VARIABLE(double, INERTIA_ROT_Z)
    KRATOS_CREATE_VARIABLE(Vector, LOCAL_AXES_VECTOR)
    KRATOS_CREATE_VARIABLE(bool, LUMPED_MASS_MATRIX)
	KRATOS_CREATE_VARIABLE(Vector, LOCAL_INERTIA_VECTOR)

    // Shell generalized variables
    KRATOS_CREATE_VARIABLE( Matrix, SHELL_STRAIN )
    KRATOS_CREATE_VARIABLE( Matrix, SHELL_STRAIN_GLOBAL )
    KRATOS_CREATE_VARIABLE( Matrix, SHELL_CURVATURE )
    KRATOS_CREATE_VARIABLE( Matrix, SHELL_CURVATURE_GLOBAL )
    KRATOS_CREATE_VARIABLE( Matrix, SHELL_FORCE )
    KRATOS_CREATE_VARIABLE( Matrix, SHELL_FORCE_GLOBAL )
    KRATOS_CREATE_VARIABLE( Matrix, SHELL_MOMENT )
    KRATOS_CREATE_VARIABLE( Matrix, SHELL_MOMENT_GLOBAL )

    // Membrane variables
    KRATOS_CREATE_VARIABLE(Vector, MEMBRANE_PRESTRESS)

    // Cross section
    KRATOS_CREATE_VARIABLE( ShellCrossSection::Pointer, SHELL_CROSS_SECTION )
    KRATOS_CREATE_VARIABLE( int, SHELL_CROSS_SECTION_OUTPUT_PLY_ID )
    KRATOS_CREATE_VARIABLE( double, SHELL_CROSS_SECTION_OUTPUT_PLY_LOCATION )

    // Nodal stiffness for the nodal concentrated element
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( NODAL_STIFFNESS )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( NODAL_DAMPING_RATIO )

    // CONDITIONS
    /* Beam conditions */
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( POINT_MOMENT )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( LOCAL_POINT_MOMENT )
    /* Torque conditions */
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( POINT_TORQUE )

    // Adding the SPRISM EAS variables
    KRATOS_CREATE_VARIABLE(double, ALPHA_EAS);
    KRATOS_CREATE_VARIABLE(bool, EAS_IMP);
    KRATOS_CREATE_VARIABLE(bool, SPRISM_TL_UL);

    // Adding the SPRISM additional variables
    KRATOS_CREATE_VARIABLE(double, ANG_ROT);

    // Adding the Sprism number of transversal integration points
    KRATOS_CREATE_VARIABLE(int, NINT_TRANS);

    // Adding the SPRISM variable to deactivate the quadratic interpolation
    KRATOS_CREATE_VARIABLE(bool, QUAD_ON);

    // Additional strain measures
    KRATOS_CREATE_VARIABLE(Vector, HENCKY_STRAIN_VECTOR);
    KRATOS_CREATE_VARIABLE(Matrix, HENCKY_STRAIN_TENSOR);

    KRATOS_CREATE_VARIABLE(double, VON_MISES_STRESS ) 

    KRATOS_CREATE_VARIABLE(Matrix, REFERENCE_DEFORMATION_GRADIENT);
    KRATOS_CREATE_VARIABLE(double, REFERENCE_DEFORMATION_GRADIENT_DETERMINANT);
    
    // Rayleigh variables
    KRATOS_CREATE_VARIABLE(double,  RAYLEIGH_ALPHA )
    KRATOS_CREATE_VARIABLE(double,  RAYLEIGH_BETA )

    // Nodal load variables
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( POINT_LOAD )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( LINE_LOAD )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( SURFACE_LOAD )

    // Condition load variables
    KRATOS_CREATE_VARIABLE(Vector, POINT_LOADS_VECTOR )
    KRATOS_CREATE_VARIABLE(Vector, LINE_LOADS_VECTOR )
    KRATOS_CREATE_VARIABLE(Vector, SURFACE_LOADS_VECTOR )
    KRATOS_CREATE_VARIABLE(Vector, POSITIVE_FACE_PRESSURES_VECTOR )
    KRATOS_CREATE_VARIABLE(Vector, NEGATIVE_FACE_PRESSURES_VECTOR )


}
