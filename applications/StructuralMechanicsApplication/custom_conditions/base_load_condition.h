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

#if !defined(KRATOS_BASE_LOAD_CONDITION_3D_H_INCLUDED )
#define  KRATOS_BASE_LOAD_CONDITION_3D_H_INCLUDED

// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{
    
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION)  BaseLoadCondition
    : public Condition
{
public:

    ///@name Type Definitions
    ///@{
    
    // Counted pointer of BaseLoadCondition
    KRATOS_CLASS_POINTER_DEFINITION( BaseLoadCondition );

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    BaseLoadCondition()
    {};

    // Constructor using an array of nodes
    BaseLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry ):Condition(NewId,pGeometry)
    {};

    // Constructor using an array of nodes with properties
    BaseLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ):Condition(NewId,pGeometry,pProperties)
    {};

    // Destructor
    virtual ~BaseLoadCondition()
    {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * Called to initialize the element.
     * Must be called before any calculation is done
     */
    virtual void Initialize() override;

    /**
     * Called at the beginning of each solution step
     * @param rCurrentProcessInfo: the current process info instance
     */
    virtual void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo) override;

    /**
     * This is called for non-linear analysis at the beginning of the iteration process
     * @param rCurrentProcessInfo: the current process info instance
     */
    virtual void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * This is called for non-linear analysis at the beginning of the iteration process
     * @param rCurrentProcessInfo: the current process info instance
     */
    virtual void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;
    
    /**
     * Called at the end of eahc solution step
     * @param rCurrentProcessInfo: the current process info instance
     */
    virtual void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo) override;
    
    /**
     * Sets on rResult the ID's of the element degrees of freedom
     * @param rResult: The vector containing the equation id
     * @param rCurrentProcessInfo: The current process info instance
     */
    virtual void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo 
        ) override;
    
    /**
     * Sets on rElementalDofList the degrees of freedom of the considered element geometry
     * @param rElementalDofList: The vector containing the dof of the element
     * @param rCurrentProcessInfo: The current process info instance
     */
    virtual void GetDofList(
        DofsVectorType& ElementalDofList,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * Sets on rValues the nodal displacements
     * @param rValues: The values of displacements
     * @param Step: The step to be computed
     */
    virtual void GetValuesVector(
        Vector& rValues,
        int Step = 0 
        ) override;

    /**
     * Sets on rValues the nodal velocities
     * @param rValues: The values of velocities
     * @param Step: The step to be computed
     */
    virtual void GetFirstDerivativesVector(
        Vector& rValues,
        int Step = 0 
        ) override;
    
    /**
     * Sets on rValues the nodal accelerations
     * @param rValues: The values of accelerations
     * @param Step: The step to be computed
     */
    virtual void GetSecondDerivativesVector(
        Vector& rValues,
        int Step = 0 
        ) override;

    /**
     * This function provides a more general interface to the element. 
     * It is designed so that rLHSvariables and rRHSvariables are passed to the element thus telling what is the desired output
     * @param rLeftHandSideMatrices: container with the output left hand side matrices
     * @param rLHSVariables: paramter describing the expected LHSs
     * @param rRightHandSideVectors: container for the desired RHS output
     * @param rRHSVariables: parameter describing the expected RHSs
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix, 
        VectorType& rRightHandSideVector, 
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * This is called during the assembling process in order to calculate the elemental right hand side vector only
      * @param rRightHandSideVector: the elemental right hand side vector
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector, 
        ProcessInfo& rCurrentProcessInfo
        ) override;
        
    /**
      * This is called during the assembling process in order to calculate the elemental mass matrix
      * @param rMassMatrix: the elemental mass matrix
      * @param rCurrentProcessInfo: The current process info instance
      */
    virtual void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo 
        ) override;
    
    /**
      * This is called during the assembling process in order
      * to calculate the elemental damping matrix
      * @param rDampingMatrix: the elemental damping matrix
      * @param rCurrentProcessInfo: The current process info instance
      */
    virtual void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        ProcessInfo& rCurrentProcessInfo 
        ) override;

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    virtual int Check( const ProcessInfo& rCurrentProcessInfo ) override;
    
    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{
    
protected:
    
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{
    
    /**
     * This functions calculates both the RHS and the LHS
     * @param rLeftHandSideMatrix: The LHS
     * @param rRightHandSideVector: The RHS
     * @param rCurrentProcessInfo: The current process info instance
     * @param CalculateStiffnessMatrixFlag: The flag to set if compute the LHS
     * @param CalculateResidualVectorFlag: The flag to set if compute the RHS
     */
    virtual void CalculateAll(
        MatrixType& rLeftHandSideMatrix, 
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
        );
    
    /**
     * This functions computes the integration weight to consider
     * @param IntegrationPoints: The array containing the integration points
     * @param PointNumber: The id of the integration point considered
     * @param detJ: The determinant of the jacobian of the element
     */
    virtual double GetIntegrationWeight(
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const unsigned int PointNumber,
        const double detJ
        );
    
    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{



    ///@}
    ///@name Private Operators
    ///@{
    
    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
    }

    virtual void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
    }

}; // class BaseLoadCondition.

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.

#endif // KRATOS_BASE_LOAD_CONDITION_3D_H_INCLUDED  defined 
