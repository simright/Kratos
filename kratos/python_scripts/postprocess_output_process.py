from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return PostprocessOutputProcess(Model, settings["Parameters"])

# All the processes python processes should be derived from "Process"
class PostprocessOutputProcess(Process):
    """This process wraps all the alternative methods to postprocess the current problem

    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, Model, settings ):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """

        KratosMultiphysics.Process.__init__(self)

        # In case of GiD post process
        if settings["Parameters"]["postprocess_parameters"]["result_file_configuration"].Has("gidpost_flags"):
            model_part = Model[settings["Parameters"]["model_part_name"].GetString()]
            output_name = settings["Parameters"]["output_name"].GetString()
            postprocess_parameters = settings["Parameters"]["postprocess_parameters"]
            import gid_output_process
            self.post_process = gid_output_process.GiDOutputProcess(model_part, output_name, postprocess_parameters)
        elif settings["Parameters"]["postprocess_parameters"]["result_file_configuration"].Has("vtk_flags"): # In case of VTK post process
            vtk_params = KratosMultiphysics.Parameters('''{ }''')

            vtk_params.AddValue("model_part_name", settings["Parameters"]["model_part_name"])
            vtk_params.AddValue("file_format", settings["Parameters"]["postprocess_parameters"]["result_file_configuration"]["vtk_flags"]["file_format"])
            vtk_params.AddValue("output_control_type", settings["Parameters"]["postprocess_parameters"]["result_file_configuration"]["output_control_type"])
            vtk_params.AddValue("output_frequency", settings["Parameters"]["postprocess_parameters"]["result_file_configuration"]["output_frequency"])
            vtk_params.AddValue("nodal_solution_step_data_variables", settings["Parameters"]["postprocess_parameters"]["result_file_configuration"]["nodal_results"])
            vtk_params.AddValue("nodal_data_value_variables", settings["Parameters"]["postprocess_parameters"]["result_file_configuration"]["nodal_nonhistorical_results"])

            import vtk_output_process
            self.post_process = vtk_output_process.VtkOutputProcess(Model, vtk_params)

        elif settings["Parameters"]["postprocess_parameters"]["result_file_configuration"].Has("hdf5_flags"): # In case of HDF5 post process
            # TODO

    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.post_process.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed before starting the time loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.post_process.ExecuteInitializeSolutionStep()

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.post_process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.post_process.ExecuteFinalizeSolutionStep()

    def IsOutputStep(self):
        """ This method returns if the current step is an output

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        return self.post_process.IsOutputStep()

    def PrintOutput(self):
        """ This method is returns the post process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.post_process.PrintOutput()

    def ExecuteBeforeOutputStep(self):
        """ This method is executed right before the ouput process computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.post_process.ExecuteBeforeOutputStep()

    def ExecuteAfterOutputStep(self):
        """ This method is executed right after the ouput process computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.post_process.ExecuteAfterOutputStep()

    def ExecuteFinalize(self):
        """ This method is executed in order to finalize the current computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.post_process.ExecuteFinalize()
