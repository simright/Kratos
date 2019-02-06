//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra, Veronika Singer
//


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_conditions/particle_based_conditions/mpm_particle_point_load_condition.h"
#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{
    //******************************* CONSTRUCTOR ****************************************
    //************************************************************************************

    MPMParticlePointLoadCondition::MPMParticlePointLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry )
        : MPMParticleBaseLoadCondition( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    MPMParticlePointLoadCondition::MPMParticlePointLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
        : MPMParticleBaseLoadCondition( NewId, pGeometry, pProperties )
    {
    }

    //********************************* CREATE *******************************************
    //************************************************************************************

    Condition::Pointer MPMParticlePointLoadCondition::Create(IndexType NewId,GeometryType::Pointer pGeom,PropertiesType::Pointer pProperties) const
    {
        return Kratos::make_shared<MPMParticlePointLoadCondition>(NewId, pGeom, pProperties);
    }

    //************************************************************************************
    //************************************************************************************

    Condition::Pointer MPMParticlePointLoadCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
    {
        return Kratos::make_shared<MPMParticlePointLoadCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
    }

    //******************************* DESTRUCTOR *****************************************
    //************************************************************************************

    MPMParticlePointLoadCondition::~MPMParticlePointLoadCondition()
    {
    }

    //*************************COMPUTE CURRENT DISPLACEMENT*******************************
    //************************************************************************************
    /*
    This function convert the computed nodal displacement into matrix of (number_of_nodes, dimension)
    */
    Matrix& MPMParticlePointLoadCondition::CalculateCurrentDisp(Matrix & rCurrentDisp, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        GeometryType& rGeom = GetGeometry();
        const unsigned int number_of_nodes = rGeom.PointsNumber();
        const unsigned int dimension = rGeom.WorkingSpaceDimension();

        rCurrentDisp = ZeroMatrix(number_of_nodes, dimension);

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            const array_1d<double, 3 > & current_displacement  = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT);

            for ( unsigned int j = 0; j < dimension; j++ )
            {
                rCurrentDisp(i,j) = current_displacement[j];
            }
        }

        return rCurrentDisp;

        KRATOS_CATCH( "" )
    }


    void MPMParticlePointLoadCondition::UpdateMPC_Coord( GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);
        const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        const array_1d<double,3> & xg = this->GetValue(MPC_COORD);

        array_1d<double,3> delta_xg = ZeroVector(3);

        rVariables.N = this->MPMShapeFunctionPointValues(rVariables.N, xg);

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            if (rVariables.N[i] > 1e-16)
            {
                for ( unsigned int j = 0; j < dimension; j++ )
                {
                    delta_xg[j] += rVariables.N[i] * rVariables.CurrentDisp(i,j);
                }
            }
        }
        KRATOS_WATCH(rVariables.CurrentDisp);
        KRATOS_WATCH(delta_xg);
        // Update the Material Point Condition Position
        const array_1d<double,3>& new_xg = xg + delta_xg ;
        this -> SetValue(MPC_COORD,new_xg);
        KRATOS_WATCH(xg);
        KRATOS_WATCH(new_xg);

        KRATOS_CATCH( "" )
    }

    void MPMParticlePointLoadCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag
        )
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int block_size = this->GetBlockSize();

        // Resizing as needed the LHS
        const unsigned int matrix_size = number_of_nodes * block_size;

        if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
        {
            if ( rLeftHandSideMatrix.size1() != matrix_size )
            {
                rLeftHandSideMatrix.resize( matrix_size, matrix_size, false );
            }

            noalias( rLeftHandSideMatrix ) = ZeroMatrix(matrix_size,matrix_size); //resetting LHS
        }

        // Resizing as needed the RHS
        if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
        {
            if ( rRightHandSideVector.size( ) != matrix_size )
            {
                rRightHandSideVector.resize( matrix_size, false );
            }

            noalias( rRightHandSideVector ) = ZeroVector( matrix_size ); //resetting RHS
        }

        // Get imposed displacement and normal vector
        const array_1d<double, 3 > & xg_c = this->GetValue(MPC_COORD);
        const array_1d<double, 3 > & Point_Load = this->GetValue (POINT_LOAD);
        Matrix Nodal_Force = ZeroMatrix(3,number_of_nodes);


        // Prepare variables
        GeneralVariables Variables;

        // Calculating shape function
        Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg_c);

        // Here MP contribution in terms of force are added
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            for (unsigned int j = 0; j < dimension; j++)
            {
                Nodal_Force(j,i) = Variables.N[i] * Point_Load[j];
            }
        }

        for (unsigned int ii = 0; ii < number_of_nodes; ++ii)
        {
            const unsigned int base = ii*dimension;

            for(unsigned int k = 0; k < dimension; ++k)
            {
                rRightHandSideVector[base + k] += GetPointLoadIntegrationWeight() * Nodal_Force(k,ii);
            }
        }
        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************

    double MPMParticlePointLoadCondition::GetPointLoadIntegrationWeight()
    {
        return 1.0;
    }

    /**
     * Shape function values in given point. This method calculate the shape function
     * vector in given point.
     *
     * @param rPoint point which shape function values have to
     * be calculated in it.
     *
     * @return Vector of double which is shape function vector \f$ N \f$ in given point.
     *
    */
    Vector& MPMParticlePointLoadCondition::MPMShapeFunctionPointValues( Vector& rResult, const array_1d<double,3>& rPoint )
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        array_1d<double,3> rPointLocal = ZeroVector(3);
        rPointLocal = GetGeometry().PointLocalCoordinates(rPointLocal, rPoint);

        if (dimension == 2)
        {
            if (number_of_nodes == 3)
            {
                rResult.resize(3, false);

                rResult[0] = 1 - rPointLocal[0] - rPointLocal[1] ;
                rResult[1] = rPointLocal[0] ;
                rResult[2] = rPointLocal[1];
            }
            else if (number_of_nodes == 4)
            {
                rResult.resize(4, false);

                rResult[0] = 0.25 * (1 - rPointLocal[0]) * (1 - rPointLocal[1]) ;
                rResult[1] = 0.25 * (1 + rPointLocal[0]) * (1 - rPointLocal[1]) ;
                rResult[2] = 0.25 * (1 + rPointLocal[0]) * (1 + rPointLocal[1]) ;
                rResult[3] = 0.25 * (1 - rPointLocal[0]) * (1 + rPointLocal[1]) ;
            }
        }
        else if (dimension == 3)
        {
            if (number_of_nodes == 4)
            {
                rResult.resize(4, false);

                rResult[0] =  1.0-(rPointLocal[0]+rPointLocal[1]+rPointLocal[2]) ;
                rResult[1] = rPointLocal[0] ;
                rResult[2] = rPointLocal[1];
                rResult[3] = rPointLocal[2];
            }
            else if (number_of_nodes == 8)
            {
                rResult.resize(8, false);

                // Shape Functions (if the first node of the connettivity is the node at (-1,-1,-1))
                // NOTE: Implemented based on Carlos Felippa's Lecture on AFEM Chapter 11
                rResult[0] = 0.125 * (1 - rPointLocal[0]) * (1 - rPointLocal[1]) * (1 - rPointLocal[2]) ;
                rResult[1] = 0.125 * (1 + rPointLocal[0]) * (1 - rPointLocal[1]) * (1 - rPointLocal[2]) ;
                rResult[2] = 0.125 * (1 + rPointLocal[0]) * (1 + rPointLocal[1]) * (1 - rPointLocal[2]) ;
                rResult[3] = 0.125 * (1 - rPointLocal[0]) * (1 + rPointLocal[1]) * (1 - rPointLocal[2]) ;
                rResult[4] = 0.125 * (1 - rPointLocal[0]) * (1 - rPointLocal[1]) * (1 + rPointLocal[2]) ;
                rResult[5] = 0.125 * (1 + rPointLocal[0]) * (1 - rPointLocal[1]) * (1 + rPointLocal[2]) ;
                rResult[6] = 0.125 * (1 + rPointLocal[0]) * (1 + rPointLocal[1]) * (1 + rPointLocal[2]) ;
                rResult[7] = 0.125 * (1 - rPointLocal[0]) * (1 + rPointLocal[1]) * (1 + rPointLocal[2]) ;
            }
        }

        return rResult;

        KRATOS_CATCH( "" )
    }

    void MPMParticlePointLoadCondition::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
    {
        GeneralVariables Variables;
        this->UpdateMPC_Coord(Variables, rCurrentProcessInfo);

        //mFinalizedStep = true;

    }




} // Namespace Kratos


