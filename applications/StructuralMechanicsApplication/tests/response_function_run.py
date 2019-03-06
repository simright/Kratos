# Import Kratos core and apps
import os, shutil
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
import structural_response_function_factory
import KratosMultiphysics.kratos_utilities as kratos_utils
import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication


with open("response_function_tests/adjoint_strain_energy_response_parameters_3D_truss.json",'r') as parameter_file:
    parameters = KratosMultiphysics.Parameters( parameter_file.read())

model = KratosMultiphysics.Model()
response_function = structural_response_function_factory.CreateResponseFunction("dummy", parameters["kratos_response_settings"], model)

model_part_primal = response_function.primal_model_part
model_part_adjoint = response_function.adjoint_model_part

response_function.RunCalculation(calculate_gradient=True)

model_part_name = "primal_output_truss"
for name in os.listdir():
    if name.find(model_part_name) == 0:
        kratos_utils.DeleteFileIfExisting(name)