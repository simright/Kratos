//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez and Riccardo Rossi
//

#include "potential_wall_condition.h"

namespace Kratos
{
///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer PotentialWallCondition<TDim, TNumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Condition::Pointer(Kratos::make_shared<PotentialWallCondition>(
        NewId, GetGeometry().Create(ThisNodes), pProperties));
    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer PotentialWallCondition<TDim, TNumNodes>::Create(
    IndexType NewId, Condition::GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(
        Kratos::make_shared<PotentialWallCondition>(NewId, pGeom, pProperties));
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer PotentialWallCondition<TDim, TNumNodes>::Clone(IndexType NewId,
                                                                  NodesArrayType const& rThisNodes) const
{
    Condition::Pointer pNewCondition =
        Create(NewId, GetGeometry().Create(rThisNodes), pGetProperties());

    pNewCondition->SetData(this->GetData());
    pNewCondition->SetFlags(this->GetFlags());

    return pNewCondition;
}

template <unsigned int TDim, unsigned int TNumNodes>
void PotentialWallCondition<TDim, TNumNodes>::Initialize()
{
    KRATOS_TRY;

    const array_1d<double, 3>& rNormal = this->GetValue(NORMAL);

    KRATOS_ERROR_IF(norm_2(rNormal) < std::numeric_limits<double>::epsilon())
        << "Error on condition -> " << this->Id() << "\n"
        << "NORMAL must be calculated before using this condition." << std::endl;

    if (mInitializeWasPerformed)
        return;

    mInitializeWasPerformed = true;

    const GeometryType& rGeom = this->GetGeometry();
    WeakPointerVector<Element> ElementCandidates;
    GetElementCandidates(ElementCandidates, rGeom);

    std::vector<IndexType> NodeIds, ElementNodeIds;
    GetSortedIds(NodeIds, rGeom);
    FindParentElement(NodeIds, ElementNodeIds, ElementCandidates);

    KRATOS_ERROR_IF(!mpElement.lock())
        << "error in condition # " << this->Id() << "\n"
        << "Condition cannot find parent element" << std::endl;
    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void PotentialWallCondition<TDim, TNumNodes>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                                                    ProcessInfo& rCurrentProcessInfo)
{
    VectorType RHS;
    CalculateLocalSystem(rLeftHandSideMatrix, RHS, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes>
void PotentialWallCondition<TDim, TNumNodes>::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                                     ProcessInfo& rCurrentProcessInfo)
{
    // TODO: improve speed
    Matrix tmp;
    CalculateLocalSystem(tmp, rRightHandSideVector, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes>
void PotentialWallCondition<TDim, TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TNumNodes)
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);
    if (rRightHandSideVector.size() != TNumNodes)
        rRightHandSideVector.resize(TNumNodes, false);
    rLeftHandSideMatrix.clear();

    array_1d<double, 3> An;
    if (TDim == 2)
        CalculateNormal2D(An);
    else
        CalculateNormal3D(An);

    const double density_infinity = rCurrentProcessInfo[DENSITY];

    const PotentialWallCondition& r_this = *this;
    const array_1d<double, 3>& v = r_this.GetValue(VELOCITY_INFINITY);
    const double value = density_infinity*inner_prod(v, An) / static_cast<double>(TNumNodes);

    for (unsigned int i = 0; i < TNumNodes; ++i)
        rRightHandSideVector[i] = value;
}

template <unsigned int TDim, unsigned int TNumNodes>
int PotentialWallCondition<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    int Check =
        Condition::Check(rCurrentProcessInfo); // Checks id > 0 and area > 0

    if (Check != 0)
    {
        return Check;
    }
    else
    {
        // Check that all required variables have been registered
        KRATOS_CHECK_VARIABLE_KEY(VELOCITY_POTENTIAL);
        KRATOS_CHECK_VARIABLE_KEY(AUXILIARY_VELOCITY_POTENTIAL);

        // Checks on nodes

        // Check that the element's nodes contain all required
        // SolutionStepData and Degrees of freedom
        for (unsigned int i = 0; i < this->GetGeometry().size(); ++i)
        {
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY_POTENTIAL,
                                                this->GetGeometry()[i]);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(AUXILIARY_VELOCITY_POTENTIAL,
                                                this->GetGeometry()[i]);

            return Check;
        }
    }
    return 0;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
void PotentialWallCondition<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType& rResult,
                                                               ProcessInfo& rCurrentProcessInfo)
{
    if (rResult.size() != TNumNodes)
        rResult.resize(TNumNodes, false);

    for (unsigned int i = 0; i < TNumNodes; i++)
        rResult[i] = GetGeometry()[i].GetDof(VELOCITY_POTENTIAL).EquationId();
}

template <unsigned int TDim, unsigned int TNumNodes>
void PotentialWallCondition<TDim, TNumNodes>::GetDofList(DofsVectorType& ConditionDofList,
                                                         ProcessInfo& CurrentProcessInfo)
{
    if (ConditionDofList.size() != TNumNodes)
        ConditionDofList.resize(TNumNodes);

    for (unsigned int i = 0; i < TNumNodes; i++)
        ConditionDofList[i] = GetGeometry()[i].pGetDof(VELOCITY_POTENTIAL);
}

template <unsigned int TDim, unsigned int TNumNodes>
void PotentialWallCondition<TDim, TNumNodes>::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    std::vector<double> rValues;
    ElementPointerType pElem = pGetElement();
    pElem->GetValueOnIntegrationPoints(PRESSURE, rValues, rCurrentProcessInfo);
    this->SetValue(PRESSURE, rValues[0]);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template <unsigned int TDim, unsigned int TNumNodes>
std::string PotentialWallCondition<TDim, TNumNodes>::Info() const
{
    std::stringstream buffer;
    this->PrintInfo(buffer);
    return buffer.str();
}

template <unsigned int TDim, unsigned int TNumNodes>
void PotentialWallCondition<TDim, TNumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "PotentialWallCondition" << TDim << "D #" << this->Id();
}

template <unsigned int TDim, unsigned int TNumNodes>
void PotentialWallCondition<TDim, TNumNodes>::PrintData(std::ostream& rOStream) const
{
    this->pGetGeometry()->PrintData(rOStream);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template <unsigned int TDim, unsigned int TNumNodes>
void PotentialWallCondition<TDim, TNumNodes>::CalculateNormal2D(array_1d<double, 3>& An) const
{
    const Geometry<Node<3>>& pGeometry = this->GetGeometry();

    An[0] = pGeometry[1].Y() - pGeometry[0].Y();
    An[1] = -(pGeometry[1].X() - pGeometry[0].X());
    An[2] = 0.00;
}

template <unsigned int TDim, unsigned int TNumNodes>
void PotentialWallCondition<TDim, TNumNodes>::CalculateNormal3D(array_1d<double, 3>& An) const
{
    const Geometry<Node<3>>& pGeometry = this->GetGeometry();

    array_1d<double, 3> v1, v2;
    v1[0] = pGeometry[1].X() - pGeometry[0].X();
    v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
    v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

    v2[0] = pGeometry[2].X() - pGeometry[0].X();
    v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
    v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

    MathUtils<double>::CrossProduct(An, v1, v2);
    An *= 0.5;
}

// serializer

template <unsigned int TDim, unsigned int TNumNodes>
void PotentialWallCondition<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
}

template <unsigned int TDim, unsigned int TNumNodes>
void PotentialWallCondition<TDim, TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
}

// private operators

template <unsigned int TDim, unsigned int TNumNodes>
inline Element::Pointer PotentialWallCondition<TDim, TNumNodes>::pGetElement() const
{
    KRATOS_ERROR_IF_NOT(mpElement.lock() != 0)
        << "No element found for condition #" << this->Id() << std::endl;
    return mpElement.lock();
}

template <unsigned int TDim, unsigned int TNumNodes>
void PotentialWallCondition<TDim, TNumNodes>::GetElementCandidates(
    WeakPointerVector<Element>& ElementCandidates, const GeometryType& rGeom) const
{
    for (SizeType i = 0; i < TDim; i++)
    {
        const WeakPointerVector<Element>& rNodeElementCandidates =
            rGeom[i].GetValue(NEIGHBOUR_ELEMENTS);
        for (SizeType j = 0; j < rNodeElementCandidates.size(); j++)
            ElementCandidates.push_back(rNodeElementCandidates(j));
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void PotentialWallCondition<TDim, TNumNodes>::GetSortedIds(std::vector<IndexType>& Ids,
                                                           const GeometryType& rGeom) const
{
    Ids.resize(rGeom.PointsNumber());
    for (SizeType i = 0; i < Ids.size(); i++)
        Ids[i] = rGeom[i].Id();
    std::sort(Ids.begin(), Ids.end());
}

template <unsigned int TDim, unsigned int TNumNodes>
void PotentialWallCondition<TDim, TNumNodes>::FindParentElement(
    std::vector<IndexType>& NodeIds,
    std::vector<IndexType>& ElementNodeIds,
    WeakPointerVector<Element> ElementCandidates)
{
    for (SizeType i = 0; i < ElementCandidates.size(); i++)
    {
        GeometryType& rElemGeom = ElementCandidates[i].GetGeometry();
        GetSortedIds(ElementNodeIds, rElemGeom);

        if (std::includes(ElementNodeIds.begin(), ElementNodeIds.end(),
                          NodeIds.begin(), NodeIds.end()))
        {
            mpElement = ElementCandidates(i);
            return;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template class instantiation

template class PotentialWallCondition<2, 2>;
template class PotentialWallCondition<3, 3>;

} // namespace Kratos