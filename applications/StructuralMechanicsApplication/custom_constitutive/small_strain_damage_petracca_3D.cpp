// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo & Philip Kalkbrenner
//
// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_constitutive/small_strain_damage_petracca_3D.h"

#include "structural_mechanics_application_variables.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallStrainDamagePetracca3D::SmallStrainDamagePetracca3D()
    : ElasticIsotropic3D()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

SmallStrainDamagePetracca3D::SmallStrainDamagePetracca3D(const SmallStrainDamagePetracca3D& rOther)
    : ElasticIsotropic3D(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer SmallStrainDamagePetracca3D::Clone() const
{
    SmallStrainDamagePetracca3D::Pointer p_clone(new SmallStrainDamagePetracca3D(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

SmallStrainDamagePetracca3D::~SmallStrainDamagePetracca3D()
{
};

void SmallStrainDamagePetracca3D::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // To implement
}



} // Namespace Kratos
