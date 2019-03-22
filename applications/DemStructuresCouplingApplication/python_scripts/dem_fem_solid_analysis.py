# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Import system python modules
import sys

# Import kratos core and applications
import KratosMultiphysics
# import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

import KratosMultiphysics.DEMApplication as Dem
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem

sys.stdout.flush()

from solid_analysis import Solution

class DemFemSolidAnalysis(Solution):

    def __init__(self, model, project_parameters="ProjectParameters.json", file_name=None):

        super(DemFemSolidAnalysis,self).__init__(model,project_parameters,file_name)

        self.model = model
        self.project_parameters = project_parameters

    def _set_solution_step_variables(self):
        super(DemFemSolidAnalysis,self)._set_solution_step_variables()

        self._model.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        self._model.main_model_part.AddNodalSolutionStepVariable(DemFem.DEM_SURFACE_LOAD)
        self._model.main_model_part.AddNodalSolutionStepVariable(DemFem.BACKUP_LAST_STRUCTURAL_VELOCITY)
        self._model.main_model_part.AddNodalSolutionStepVariable(DemFem.BACKUP_LAST_STRUCTURAL_DISPLACEMENT)
        self._model.main_model_part.AddNodalSolutionStepVariable(DemFem.SMOOTHED_STRUCTURAL_VELOCITY)
        self._model.main_model_part.AddNodalSolutionStepVariable(Dem.DELTA_DISPLACEMENT)
        self._model.main_model_part.AddNodalSolutionStepVariable(Dem.DEM_PRESSURE)
        self._model.main_model_part.AddNodalSolutionStepVariable(Dem.DEM_NODAL_AREA)
        self._model.main_model_part.AddNodalSolutionStepVariable(Dem.ELASTIC_FORCES)
        self._model.main_model_part.AddNodalSolutionStepVariable(Dem.CONTACT_FORCES)
        self._model.main_model_part.AddNodalSolutionStepVariable(Dem.TANGENTIAL_ELASTIC_FORCES)
        self._model.main_model_part.AddNodalSolutionStepVariable(Dem.SHEAR_STRESS)
        self._model.main_model_part.AddNodalSolutionStepVariable(Dem.NON_DIMENSIONAL_VOLUME_WEAR)
        self._model.main_model_part.AddNodalSolutionStepVariable(Dem.IMPACT_WEAR)

    def FinalizeSolutionStep(self):
        clock_time = self._start_time_measuring()

        # Execution at the end of the solution step
        self._output.ExecuteFinalizeSolutionStep()

        # Processes to be executed at the end of the solution step
        self._processes.ExecuteFinalizeSolutionStep()

        # Execution at the end of the solution step
        self._model.ExecuteFinalizeSolutionStep()

        self._stop_time_measuring(clock_time, "Finalize Step", self.report)

if __name__ == "__main__":
    DemFemSolidAnalysis().Run()