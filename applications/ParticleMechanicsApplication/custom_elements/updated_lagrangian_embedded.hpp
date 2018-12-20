//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					    Kratos default license: kratos/license.txt
//
//  Main authors:   Bodhinanda Chandra
//


#if !defined(KRATOS_UPDATED_LAGRANGIAN_EMBEDDED_H_INCLUDED )
#define      KRATOS_UPDATED_LAGRANGIAN_EMBEDDED_H_INCLUDED

// System includes

// External includes
#include "custom_elements/updated_lagrangian.hpp"


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

/// Large Displacement Lagrangian Element for 3D and 2D geometries. (base class)

/**
 * Implements a Large Displacement Lagrangian definition for structural analysis.
 * This works for arbitrary geometries in 3D and 2D (base class)
 */

class UpdatedLagrangianEmbedded
    : public UpdatedLagrangian
{
public:

    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///StressMeasure from constitutive laws
    typedef ConstitutiveLawType::StressMeasure StressMeasureType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// Counted pointer of LargeDisplacementElement
    KRATOS_CLASS_POINTER_DEFINITION( UpdatedLagrangianEmbedded );
    ///@}

    /**
     * Parameters to be used in the Element as they are. Direct interface to Parameters Struct
     */



public:


    ///@name Life Cycle
    ///@{

    /// Empty constructor needed for serialization
    UpdatedLagrangianEmbedded();


    /// Default constructors
    UpdatedLagrangianEmbedded(IndexType NewId, GeometryType::Pointer pGeometry);

    UpdatedLagrangianEmbedded(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    UpdatedLagrangianEmbedded(UpdatedLagrangianEmbedded const& rOther);

    /// Destructor.
    virtual ~UpdatedLagrangianEmbedded();

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    /**
     * clones the selected element variables, creating a new one
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId,
                            NodesArrayType const& ThisNodes) const override;

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
    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{

    bool mIsSlip;
    std::size_t mNumPositiveNodes;
    std::size_t mNumNegativeNodes;

    std::vector< size_t > mPositiveIndices;
    std::vector< size_t > mNegativeIndices;

    Vector mDistance;

    ///@}
    ///@name Protected Operators
    ///@{

    /**
     * Calculates the elemental contributions
     * \f$ K^e = w\,B^T\,D\,B \f$ and
     * \f$ r^e \f$
     */
    void CalculateElementalSystem(LocalSystemComponents& rLocalSystem,
                                ProcessInfo& rCurrentProcessInfo) override;


    /**
     * Function to identify whether element is embedded or not and initialize necessary variables
     */
    virtual void InitializeEmbeddedVariables();


    /**
     * Initialize Element General Variables
     */
    void InitializeGeneralVariables(GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo) override;


    /**
     * Calculate Modified Jacobian in a given point considering cut element
     */
    virtual Matrix& MPMModifiedJacobian(Matrix& rResult, const array_1d<double,3>& rPoint);

    /**
     * Modified Jacobian in given point and given a delta position.
     * This method calculate jacobian matrix in given point and a given delta position.
     *
     * @param rPoint point which jacobians has to be calculated in it.
     * @return Matrix of double which is jacobian matrix \f$ J \f$ in given point and a given delta position.
     *
     * @see DeterminantOfJacobian
     * @see InverseOfJacobian
     */
    virtual Matrix& MPMModifiedJacobianDelta(Matrix& rResult, const array_1d<double,3>& rPoint, const Matrix& rDeltaPosition);

    /**
     * Modified Shape function values in given point.
     * This method calculate the shape function vector in given point.
     *
     * @param rPoint point which shape function values have to be calculated in it.
     * @return Vector of double which is shape function vector \f$ N \f$ in given point.
     */
    virtual Vector& MPMModifiedShapeFunctionPointValues(Vector& rResult, const array_1d<double,3>& rPoint);

    /**
     * Calculate Modified Shape Function gradient local Values in a given point in 3 dimension
     */
    virtual Matrix& MPMModifiedShapeFunctionsLocalGradients(Matrix& rResult);


    /**
     * Returns true if there is any cut in the element
     */
    bool IsCut() {
        return ((mNumPositiveNodes > 0) && (mNumNegativeNodes > 0));
    }

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{
    ///@}
    ///@name Protected Inquiry
    ///@{
    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

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

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class UpdatedLagrangian

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_UPDATED_LAGRANGIAN_EMBEDDED_H_INCLUDED  defined