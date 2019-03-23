from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
# import KratosMultiphysics.SolversApplication as KratosSolver
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

# Import the mechanical solver base class
from implicit_dynamic_solver import ImplicitMonolithicSolver

def CreateSolver(custom_settings, Model):
    return DemFemSolidSolver(Model, custom_settings)

class DemFemSolidSolver(ImplicitMonolithicSolver):

    def __init__(self, Model, custom_settings):

        # Construct the base solver.
        super(DemFemSolidSolver, self).__init__(Model, custom_settings)

    def GetComputingModelPart(self):
        return self.model_part

    def ComputeDeltaTime(self):
        return self.process_info[KratosMultiphysics.DELTA_TIME]

    def AdvanceInTime(self, current_time):
        # dt = self.ComputeDeltaTime()
        # new_time = current_time + dt
        # self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        # self.main_model_part.CloneTimeStep(new_time)
        return current_time

    def Predict(self):
        pass

    def SolveSolutionStep(self):
        converged = super(DemFemSolidSolver,self).Solve()
        self.process_info[KratosSolid.CONVERGENCE_ACHIEVED] = converged
        return converged
