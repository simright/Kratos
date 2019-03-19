from __future__ import print_function, absolute_import, division
import KratosMultiphysics as KM

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import CoSimulationAnalysis

import os, json

class CoSimulationTestCase(KratosUnittest.TestCase):
    '''This class is the basis for the testing the framework
    It can be used to test complete cases with the "CoSimulation-Analysis"
    '''

    def createTest(self, problem_dir_name, parameter_file_name):
        self.problem_dir_name = problem_dir_name

        with open(os.path.join(problem_dir_name, parameter_file_name + '_parameters.json'), 'r') as parameter_file:
            self.cosim_parameters = json.load(parameter_file)

        # To avoid many prints
        echo_level = self.cosim_parameters["problem_data"]["echo_level"]
        if (echo_level == 0):
            KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
        else:
            KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.INFO)

    def runTest(self):
        model = KM.Model()
        CoSimulationAnalysis(model, self.cosim_parameters).Run()
        kratos_utils.DeleteTimeFiles(self.problem_dir_name)

    # called only once for this class, opposed of tearDown()
    @classmethod
    def tearDownClass(cls):
        # Cleaning
        kratos_utils.DeleteDirectoryIfExisting("__pycache__")