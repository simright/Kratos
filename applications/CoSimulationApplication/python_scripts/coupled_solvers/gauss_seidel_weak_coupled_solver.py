from __future__ import print_function, absolute_import, division

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_base_coupled_solver import CoSimulationBaseCouplingSolver

def Create(model, cosim_solver_settings):
    return GaussSeidelWeakCouplingSolver(model, cosim_solver_settings)

class GaussSeidelWeakCouplingSolver(CoSimulationBaseCouplingSolver):
    def SolveSolutionStep(self):
        for solver_name in self.solver_names:
            solver = self.solvers[solver_name]
            self._SynchronizeInputData(solver, solver_name)
            solver.SolveSolutionStep()
            self._SynchronizeOutputData(solver, solver_name)

    def _Name(self):
        return self.__class__.__name__