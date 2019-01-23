//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Veronika Singer
//


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_conditions/particle_based_conditions/mpm_particle_base_load_condition.h"
#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{
    //******************************* CONSTRUCTOR ****************************************
    //************************************************************************************

    MPMParticleBaseLoadCondition::MPMParticleBaseLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry )
        : MPMBaseLoadCondition( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    MPMParticleBaseLoadCondition::MPMParticleBaseLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
        : MPMBaseLoadCondition( NewId, pGeometry, pProperties )
    {
    }

    //********************************* CREATE *******************************************
    //************************************************************************************

    Condition::Pointer MPMParticleBaseLoadCondition::Create(IndexType NewId,GeometryType::Pointer pGeom,PropertiesType::Pointer pProperties) const
    {
        return Kratos::make_shared<MPMParticleBaseLoadCondition>(NewId, pGeom, pProperties);
    }

    //************************************************************************************
    //************************************************************************************

    Condition::Pointer MPMParticleBaseLoadCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
    {
        return Kratos::make_shared<MPMParticleBaseLoadCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
    }

    //******************************* DESTRUCTOR *****************************************
    //************************************************************************************

    MPMParticleBaseLoadCondition::~MPMParticleBaseLoadCondition()
    {
    }

    //************************************************************************************
    //************************************************************************************


    double MPMParticleBaseLoadCondition::GetPointLoadIntegrationWeight()
    {
        return 1.0;
    }


} // Namespace Kratos


