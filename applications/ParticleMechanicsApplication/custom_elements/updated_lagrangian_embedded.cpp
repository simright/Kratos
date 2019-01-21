//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:   Bodhinanda Chandra
//


// System includes
#include <omp.h>
#include <sstream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "custom_elements/updated_lagrangian_embedded.hpp"
#include "utilities/divide_triangle_2d_3.h"
#include "utilities/divide_tetrahedra_3d_4.h"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "custom_utilities/particle_mechanics_math_utilities.h"
#include "particle_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianEmbedded::UpdatedLagrangianEmbedded( )
    : UpdatedLagrangian( )
{
    //DO NOT CALL IT: only needed for Register and Serialization!!!
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianEmbedded::UpdatedLagrangianEmbedded( IndexType NewId, GeometryType::Pointer pGeometry )
        : UpdatedLagrangian( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianEmbedded::UpdatedLagrangianEmbedded( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : UpdatedLagrangian( NewId, pGeometry, pProperties )
{
mFinalizedStep = true;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

UpdatedLagrangianEmbedded::UpdatedLagrangianEmbedded( UpdatedLagrangianEmbedded const& rOther)
    :UpdatedLagrangian(rOther)
{
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer UpdatedLagrangianEmbedded::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new UpdatedLagrangianEmbedded( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer UpdatedLagrangianEmbedded::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    UpdatedLagrangianEmbedded NewElement (NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    return Element::Pointer( new UpdatedLagrangianEmbedded(NewElement) );
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianEmbedded::~UpdatedLagrangianEmbedded()
{
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianEmbedded::CalculateElementalSystem( LocalSystemComponents& rLocalSystem,
        ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Initialize embedded variables
    this->InitializeEmbeddedVariables();

    // Create and initialize element variables (considering whether it is cut or not):
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags &ConstitutiveLawOptions=Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    // Auxiliary terms
    Vector volume_force;

    // Compute element kinematics B, F, DN_DX ...
    this->CalculateKinematics(Variables,rCurrentProcessInfo);

    // Set general variables to constitutivelaw parameters
    this->SetGeneralVariables(Variables,Values);

    // Calculate Material Response
    /* NOTE:
    The function below will call CalculateMaterialResponseCauchy() by default and then (may)
    call CalculateMaterialResponseKirchhoff() in the constitutive_law.*/
    mConstitutiveLawVector->CalculateMaterialResponse(Values, Variables.StressMeasure);

    /* NOTE:
    The material points will have constant mass as defined at the beginning.
    However, the density and volume (integration weight) are changing every time step.*/
    // Update MP_Density
    const double MP_Density = (GetProperties()[DENSITY]) / Variables.detFT;
    this->SetValue(MP_DENSITY, MP_Density);

    // The MP_Volume (integration weight) is evaluated
    const double MP_Volume = this->GetValue(MP_MASS)/this->GetValue(MP_DENSITY);
    this->SetValue(MP_VOLUME, MP_Volume);

    if ( rLocalSystem.CalculationFlags.Is(UpdatedLagrangian::COMPUTE_LHS_MATRIX) ) // if calculation of the matrix is required
    {
        // Contributions to stiffness matrix calculated on the reference configuration
        this->CalculateAndAddLHS ( rLocalSystem, Variables, MP_Volume );
    }

    if ( rLocalSystem.CalculationFlags.Is(UpdatedLagrangian::COMPUTE_RHS_VECTOR) ) // if calculation of the vector is required
    {
        // Contribution to forces (in residual term) are calculated
        volume_force  = this->CalculateVolumeForce( volume_force, Variables );
        this->CalculateAndAddRHS ( rLocalSystem, Variables, volume_force, MP_Volume );
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianEmbedded::InitializeEmbeddedVariables()
{
    const GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = GetGeometry().size();

    // Initialize mDistance as from obtained DISTANCE variables
    if (mDistance.size() != number_of_nodes)
        mDistance.resize(number_of_nodes, false);
    mDistance = ZeroVector(number_of_nodes);

    for (size_t i = 0; i < number_of_nodes; i++)
        mDistance[i] = rGeom[i].FastGetSolutionStepValue(DISTANCE);

    // Check whether the element is a SLIP element or not
    mIsSlip = this->Is(SLIP) ? true : false;

    // Clear previous positive and negative number of nodes and indices
    mNumPositiveNodes = 0;
    mNumNegativeNodes = 0;
    mPositiveIndices.clear();
    mNegativeIndices.clear();

    // Number of positive and negative distance function values
    for (size_t i = 0; i < number_of_nodes; ++i){
        // For positive distance
        if (mDistance[i] > 0.0) {
            mNumPositiveNodes++;
            mPositiveIndices.push_back(i);
        }
        // For negative distance
        else if (mDistance[i] < 0.0) {
            mNumNegativeNodes++;
            mNegativeIndices.push_back(i);
        }
    }

    if (this->IsCut()){
        // Initialize Geometry Splitter
        const unsigned int dimension       = rGeom.WorkingSpaceDimension();
        if (dimension == 2 && number_of_nodes == 3)
            mpGeometrySplitter = Kratos::make_shared<DivideTriangle2D3>(rGeom, mDistance);
        else if (dimension == 3 && number_of_nodes == 4)
            mpGeometrySplitter = Kratos::make_shared<DivideTetrahedra3D4>(rGeom, mDistance);
        else
            KRATOS_ERROR <<  "Dimension given is wrong: Something is wrong with the given dimension in function: InitializeEmbeddedVariables" << std::endl;

        // Perform the element splitting
        mpGeometrySplitter->GenerateDivision();
        mpGeometrySplitter->GenerateIntersectionsSkin();

        // If Splitted, then generate condensation matrix
        if(this->IsSplit()){
            Matrix p_matrix;
            this->SetCondensationMatrix(p_matrix, mpGeometrySplitter->GetEdgeIdsI(), mpGeometrySplitter->GetEdgeIdsJ(), mpGeometrySplitter->GetSplitEdges());
        }

    }
}

//************************************************************************************
//************************************************************************************

// Sets the condensation matrix to transform the subdivision values to entire element ones.
void UpdatedLagrangianEmbedded::SetCondensationMatrix(
    Matrix& rIntPointCondMatrix,
    const std::vector<int>& rEdgeNodeI,
    const std::vector<int>& rEdgeNodeJ,
    const std::vector<int>& rSplitEdges) {

    const GeometryType& rGeom = GetGeometry();
    const unsigned int n_edges = rGeom.EdgesNumber();
    const unsigned int n_nodes = rGeom.PointsNumber();

    // Resize condensation matrix
    rIntPointCondMatrix.resize(n_nodes + n_edges, n_nodes, false);

    // Initialize intersection points condensation matrix
    rIntPointCondMatrix = ZeroMatrix(n_nodes + n_edges, n_nodes);

    // Fill the original geometry points main diagonal
    for (unsigned int i = 0; i < n_nodes; ++i) {
        rIntPointCondMatrix(i,i) = 1.0;
    }

    // Compute the intersection points contributions
    unsigned int row = n_nodes;
    for (unsigned int id_edge = 0; id_edge < n_edges; ++id_edge) {
        // Check if the edge has an intersection point
        if (rSplitEdges[n_nodes+id_edge] != -1) {
            // Get the nodes that compose the edge
            const unsigned int edge_node_i = rEdgeNodeI[id_edge];
            const unsigned int edge_node_j = rEdgeNodeJ[id_edge];

            // Compute the relative coordinate of the intersection point over the edge
            const double aux_node_rel_location = std::abs (mDistance[edge_node_i] / (mDistance[edge_node_j] - mDistance[edge_node_i]));

            // Store the relative coordinate values as the original geometry nodes sh. function value in the intersections
            rIntPointCondMatrix(row, edge_node_i) = 1.0 - aux_node_rel_location;
            rIntPointCondMatrix(row, edge_node_j) = aux_node_rel_location;
        }
        row++;
    }
}

//************************************************************************************
//************************************************************************************

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianEmbedded::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    rSerializer.save("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.save("DeterminantF0",mDeterminantF0);
    rSerializer.save("InverseJ0",mInverseJ0);
    rSerializer.save("DeterminantJ0",mDeterminantJ0);

}

void UpdatedLagrangianEmbedded::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    rSerializer.load("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.load("DeterminantF0",mDeterminantF0);
    rSerializer.load("InverseJ0",mInverseJ0);
    rSerializer.load("DeterminantJ0",mDeterminantJ0);
}

} // Namespace Kratos

