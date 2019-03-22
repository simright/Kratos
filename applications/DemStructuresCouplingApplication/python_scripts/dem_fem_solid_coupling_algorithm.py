from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import sys
import time as timer
import os
import weakref
import KratosMultiphysics as Kratos
from Kratos import Logger
import KratosMultiphysics.DEMApplication as Dem
sys.path.insert(0, '')
Logger.Print("Running under OpenMP........", label="DEM")
import DEM_procedures
import DEM_material_test_script
import KratosMultiphysics.SolversApplication
import KratosMultiphysics.ConstitutiveModelsApplication
import KratosMultiphysics.SolidMechanicsApplication
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem
import dem_structures_coupling_gid_output
from dem_fem_coupling_algorithm import Algorithm

class DemFemSolidCouplingAlgorithm(Algorithm):

    def __init__(self):
        self.model = Kratos.Model()

        import dem_main_script_ready_for_coupling_with_fem
        self.dem_solution = dem_main_script_ready_for_coupling_with_fem.Solution(self.model)
        self.dem_solution.coupling_algorithm = weakref.proxy(self)

        import dem_fem_solid_analysis
        structural_parameters_file_name = "ProjectParameters.json"

        self.structural_solution = dem_fem_solid_analysis.DemFemSolidAnalysis(self.model, structural_parameters_file_name)


if __name__ == "__main__":
    DemFemSolidCouplingAlgorithm().Run()
