/*
==============================================================================
KratosShallowWaterApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

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
//   Last modified by:    Miguel Masó Sotomayor
//   Date:                May 2017
//   Revision:            1.2
//
//


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/move_shallow_water_particle_utility.h"
#include "custom_utilities/move_shallow_water_particle_hu_utility.h"
#include "custom_utilities/shallow_water_variables_utility.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"



namespace Kratos
{

namespace Python
{

  void  AddCustomUtilitiesToPython()
  {
    using namespace boost::python;

    //~ typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    //~ typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    //~ typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    class_< MoveShallowWaterParticleUtility<2> > ("MoveShallowWaterParticleUtility", init<ModelPart& , int >())
        .def("MountBin", &MoveShallowWaterParticleUtility<2>::MountBin)
        .def("MoveParticles", &MoveShallowWaterParticleUtility<2>::MoveParticles)
        .def("CorrectParticlesWithoutMovingUsingDeltaVariables", &MoveShallowWaterParticleUtility<2>::CorrectParticlesWithoutMovingUsingDeltaVariables)
        .def("PreReseed", &MoveShallowWaterParticleUtility<2>::PreReseed)
        .def("PostReseed", &MoveShallowWaterParticleUtility<2>::PostReseed)
        .def("ResetBoundaryConditions", &MoveShallowWaterParticleUtility<2>::ResetBoundaryConditions)
        .def("TransferLagrangianToEulerian",&MoveShallowWaterParticleUtility<2>::TransferLagrangianToEulerian)
        .def("CalculateVelOverElemSize", &MoveShallowWaterParticleUtility<2>::CalculateVelOverElemSize)
        .def("CalculateDeltaVariables", &MoveShallowWaterParticleUtility<2>::CalculateDeltaVariables)
        .def("CopyScalarVarToPreviousTimeStep", &MoveShallowWaterParticleUtility<2>::CopyScalarVarToPreviousTimeStep)
        .def("CopyVectorVarToPreviousTimeStep", &MoveShallowWaterParticleUtility<2>::CopyVectorVarToPreviousTimeStep)
        ;

    class_< MoveShallowWaterParticleHUUtility<2> > ("MoveShallowWaterParticleHUUtility", init<ModelPart& , int >())
        .def("MountBin", &MoveShallowWaterParticleHUUtility<2>::MountBin)
        .def("MoveParticles", &MoveShallowWaterParticleHUUtility<2>::MoveParticles)
        .def("CorrectParticlesWithoutMovingUsingDeltaVariables", &MoveShallowWaterParticleHUUtility<2>::CorrectParticlesWithoutMovingUsingDeltaVariables)
        .def("PreReseed", &MoveShallowWaterParticleHUUtility<2>::PreReseed)
        .def("PostReseed", &MoveShallowWaterParticleHUUtility<2>::PostReseed)
        .def("ResetBoundaryConditions", &MoveShallowWaterParticleHUUtility<2>::ResetBoundaryConditions)
        .def("TransferLagrangianToEulerian",&MoveShallowWaterParticleHUUtility<2>::TransferLagrangianToEulerian)
        .def("CalculateVelOverElemSize", &MoveShallowWaterParticleHUUtility<2>::CalculateVelOverElemSize)
        .def("CalculateDeltaVariables", &MoveShallowWaterParticleHUUtility<2>::CalculateDeltaVariables)
        .def("CopyScalarVarToPreviousTimeStep", &MoveShallowWaterParticleHUUtility<2>::CopyScalarVarToPreviousTimeStep)
        .def("CopyVectorVarToPreviousTimeStep", &MoveShallowWaterParticleHUUtility<2>::CopyVectorVarToPreviousTimeStep)
        .def("ComputeVelocity", &MoveShallowWaterParticleHUUtility<2>::ComputeVelocity)
        ;

    class_< ShallowWaterVariablesUtility > ("ShallowWaterVariablesUtility", init<ModelPart&>())
        .def("ComputeFreeSurfaceElevation", &ShallowWaterVariablesUtility::ComputeFreeSurfaceElevation)
        .def("ComputeVelocity", &ShallowWaterVariablesUtility::ComputeVelocity)
        .def("DryWetStateConservedVariables", &ShallowWaterVariablesUtility::DryWetStateConservedVariables)
        .def("DryWetStatePrimitiveVariables", &ShallowWaterVariablesUtility::DryWetStatePrimitiveVariables)
        ;
  }


}  // namespace Python.

} // Namespace Kratos
