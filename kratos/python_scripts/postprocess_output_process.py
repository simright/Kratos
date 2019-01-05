from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return PostprocessOutputProcess(Model, settings["Parameters"])

# All the processes python processes should be derived from "Process"
class PostprocessOutputProcess(Process):
    """This process wraps all the alternative methods to postprocess the current problem.
    The options are:
        - GiD postprocess
        - VTK postprocess
        - HDF5 postprocess

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

        self.auxiliar_process = [] # Here we add the processes for auxiliar executions before printint the result

        # In case of GiD post process
        if settings["Parameters"]["postprocess_parameters"]["result_file_configuration"].Has("gidpost_flags"):
            model_part = Model[settings["Parameters"]["model_part_name"].GetString()]
            output_name = settings["Parameters"]["output_name"].GetString()
            postprocess_parameters = settings["Parameters"]["postprocess_parameters"]
            import gid_output_process
            self.post_process = gid_output_process.GiDOutputProcess(model_part, output_name, postprocess_parameters)
        elif settings["Parameters"]["postprocess_parameters"]["result_file_configuration"].Has("vtk_flags"): # In case of VTK post process

            # Creating parameters
            vtk_params = KratosMultiphysics.Parameters('''{ }''')

            # model_part_name
            if settings["Parameters"].Has("model_part_name"):
                vtk_params.AddValue("model_part_name", settings["Parameters"]["model_part_name"])

            # file_format
            if settings["Parameters"]["postprocess_parameters"]["result_file_configuration"]["vtk_flags"].Has("file_format"):
                vtk_params.AddValue("file_format", settings["Parameters"]["postprocess_parameters"]["result_file_configuration"]["vtk_flags"]["file_format"])

            # output_control_type
            if settings["Parameters"]["postprocess_parameters"]["result_file_configuration"].Has("output_control_type"):
                vtk_params.AddValue("output_control_type", settings["Parameters"]["postprocess_parameters"]["result_file_configuration"]["output_control_type"])

            # output_frequency
            if settings["Parameters"]["postprocess_parameters"]["result_file_configuration"].Has("output_frequency"):
                vtk_params.AddValue("output_frequency", settings["Parameters"]["postprocess_parameters"]["result_file_configuration"]["output_frequency"])

            # nodal_solution_step_data_variables
            if settings["Parameters"]["postprocess_parameters"]["result_file_configuration"].Has("nodal_results"):
                vtk_params.AddValue("nodal_solution_step_data_variables", settings["Parameters"]["postprocess_parameters"]["result_file_configuration"]["nodal_results"])

            # nodal_data_value_variables
            if settings["Parameters"]["postprocess_parameters"]["result_file_configuration"].Has("nodal_nonhistorical_results"):
                vtk_params.AddValue("nodal_data_value_variables", settings["Parameters"]["postprocess_parameters"]["result_file_configuration"]["nodal_nonhistorical_results"])

            import vtk_output_process
            self.post_process = vtk_output_process.VtkOutputProcess(Model, vtk_params)

        elif settings["Parameters"]["postprocess_parameters"]["result_file_configuration"].Has("hdf5_flags"): # In case of HDF5 post process

            # Importing application
            try:
                import KratosMultiphysics.HDF5Application as HDF5Application
            except ImportError as e:
                raise Exception("Compile HDF5Application in order to use HDF5 postprocess")

            # Creating parameters
            hdf5_params = KratosMultiphysics.Parameters('''{ }''')

            # model_part_name
            if settings["Parameters"].Has("model_part_name"):
                hdf5_params.AddValue("model_part_name", settings["Parameters"]["model_part_name"])

            # file_name
            if settings["Parameters"].Has("output_name"):
                hdf5_params.AddValue("file_name", settings["Parameters"]["output_name"])

            # nodal_solution_step_data_settings
            if settings["Parameters"]["postprocess_parameters"]["result_file_configuration"].Has("nodal_results"):
                hdf5_params.AddEmptyValue("nodal_solution_step_data_settings")
                hdf5_params["nodal_solution_step_data_settings"].AddValue("list_of_variables", settings["Parameters"]["postprocess_parameters"]["result_file_configuration"]["nodal_results"])

            # nodal_data_value_settings
            if settings["Parameters"]["postprocess_parameters"]["result_file_configuration"].Has("nodal_nonhistorical_results"):
                hdf5_params.AddEmptyValue("nodal_data_value_settings")
                hdf5_params["nodal_data_value_settings"].AddValue("list_of_variables", settings["Parameters"]["postprocess_parameters"]["result_file_configuration"]["nodal_nonhistorical_results"])

            from single_mesh_temporal_output_process import Factory as HDF5Factory
            self.post_process = HDF5Factory(hdf5_params, Model)

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
        # First execute the auxiliar processes
        for process in self.auxiliar_process:
            process.Execute()

        # Calling PrintOutput
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
