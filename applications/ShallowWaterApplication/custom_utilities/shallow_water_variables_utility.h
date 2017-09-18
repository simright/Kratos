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
//   Last Modified by:    Miguel Masó Sotomayor
//   Date:                September 15th 2017
//   Revision:            1.0
//
//

#if !defined(KRATOS_SHALLOW_WATER_VARIABLES_UTILITY_H_INCLUDED)
#define  KRATOS_SHALLOW_WATER_VARIABLES_UTILITY_H_INCLUDED

// System includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"

// External includes 

// Project includes
#include "shallow_water_application.h"

namespace Kratos
{

    class ShallowWaterVariablesUtility
    {
    public:

        KRATOS_CLASS_POINTER_DEFINITION(ShallowWaterVariablesUtility);

        ShallowWaterVariablesUtility(ModelPart& model_part)
            : mr_model_part(model_part)  
        {
            KRATOS_TRY
            std::cout << "Initializing shallow water variables utility" << std::endl; 
            KRATOS_CATCH("")
        }

        ~ShallowWaterVariablesUtility()
        {}

        /**
         * This method computes the free surface elevation as HEIGHT + BATHYMETRY
         */
        void ComputeFreeSurfaceElevation()
        {
            KRATOS_TRY

            ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
            #pragma omp parallel for
            for(unsigned int i = 0; i < static_cast<unsigned int>(rNodes.size()); i++)
            {
                ModelPart::NodesContainerType::iterator inode = rNodes.begin() + i;
                inode->FastGetSolutionStepValue(FREE_SURFACE_ELEVATION) = inode->FastGetSolutionStepValue(HEIGHT) + inode->FastGetSolutionStepValue(BATHYMETRY);
            }

            KRATOS_CATCH("")
        }

        /** 
         * This method computes the velocity as the MOMENTUM / HEIGHT
         */
        void ComputeVelocity()
        {
            KRATOS_TRY

            ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
            #pragma omp parallel for
            for(unsigned int i = 0; i < static_cast<unsigned int>(rNodes.size()); i++)
            {
                ModelPart::NodesContainerType::iterator inode = rNodes.begin() + i;
                inode->GetSolutionStepValue(VELOCITY) = inode->FastGetSolutionStepValue(MOMENTUM) / inode->FastGetSolutionStepValue(HEIGHT);
            }

            KRATOS_CATCH("")
        }

        void DryWetStateConservedVariables()
        {
            KRATOS_TRY

            ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
            #pragma omp parallel for
            for(unsigned int i = 0; i < static_cast<unsigned int>(rNodes.size()); i++)
            {
                ModelPart::NodesContainerType::iterator inode = rNodes.begin() + i;
                if (inode->FastGetSolutionStepValue(HEIGHT) < 1e-3)
                {
                    inode->FastGetSolutionStepValue(HEIGHT)     = 1e-5;
                    inode->FastGetSolutionStepValue(MOMENTUM_X) = 0;
                    inode->FastGetSolutionStepValue(MOMENTUM_Y) = 0;
                }
            }
            KRATOS_CATCH("")
        }

        void DryWetStatePrimitiveVariables()
        {
            KRATOS_TRY

            ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
            #pragma omp parallel for
            for(unsigned int i = 0; i < static_cast<unsigned int>(rNodes.size()); i++)
            {
                ModelPart::NodesContainerType::iterator inode = rNodes.begin() + i;
                if (inode->FastGetSolutionStepValue(HEIGHT) < 1e-3)
                {
                    inode->FastGetSolutionStepValue(HEIGHT)     = 1e-5;
                    inode->FastGetSolutionStepValue(VELOCITY_X) = 0;
                    inode->FastGetSolutionStepValue(VELOCITY_Y) = 0;
                }
            }
            KRATOS_CATCH("")
        }


    protected:

    private:

        ModelPart& mr_model_part;

    }; // class ShallowWaterVariablesUtility

} // namespace Kratos

# endif // KRATOS_SHALLOW_WATER_VARIABLES_UTILITY_H_INCLUDED
