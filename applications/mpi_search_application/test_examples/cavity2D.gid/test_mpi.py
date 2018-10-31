from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import mpi
import fluid_only_var

#
#
# setting the domain size for the problem to be solved
domain_size = 2

#
#
# ATTENTION: here the order is important

from KratosMultiphysics import *

from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.TrilinosApplication import *
from KratosMultiphysics.MetisApplication import *
from KratosMultiphysics.MPISearchApplication import *


# defining a model part
model_part = ModelPart("FluidPart")
p_model_part = ModelPart("FluidPart")

# providing the variable list to the model part
import trilinos_fs_fluid_solver
trilinos_fs_fluid_solver.AddVariables(model_part)
trilinos_fs_fluid_solver.AddVariables(p_model_part)
model_part.AddNodalSolutionStepVariable(NORMAL)
p_model_part.AddNodalSolutionStepVariable(NORMAL)

# reading a model
gid_mode_flag = GiDPostMode.GiD_PostBinary
use_multifile = MultiFileFlag.MultipleFiles
deformed_print_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly

input_file_name = "cavity2D"
gid_io = GidIO(input_file_name, gid_mode_flag, use_multifile, deformed_print_flag, write_conditions)

model_part_io_fluid = ModelPartIO(input_file_name)

print("(" + str(mpi.rank) + "," + str(mpi.size) + ")" +"before performing the division")
number_of_partitions = mpi.size  # we set it equal to the number of processors
if mpi.rank == 0:
    print("(" + str(mpi.rank) + "," + str(mpi.size) + ")" +"start partition process")
    partitioner = MetisDivideInputToPartitionsProcess(model_part_io_fluid, number_of_partitions, domain_size);
    partitioner.Execute()

print("(" + str(mpi.rank) + "," + str(mpi.size) + ")" +"division performed")
mpi.world.barrier()

MPICommSetup = SetMPICommunicatorProcess(model_part)
MPICommSetup.Execute()

print("(" + str(mpi.rank) + "," + str(mpi.size) + ")" +"Comunicator Set")

my_input_filename = input_file_name + "_" + str(mpi.rank)
model_part_io_fluid = ModelPartIO(my_input_filename)
model_part_io_fluid.ReadModelPart(model_part)

# number_of_partitions = mpi.size #we set it equal to the number of processors
# print "number_of_partitions", number_of_partitions
# partitioner = MetisPartitioningProcess(model_part, gid_io, number_of_partitions, domain_size);
# partitioner.Execute()
# print "GetRank()",GetRank()

print(model_part)

# the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(3)
p_model_part.SetBufferSize(3)

# adding dofs
trilinos_fs_fluid_solver.AddDofs(model_part)
for node in model_part.Nodes:
    node.AddDof(DISPLACEMENT_X)
    node.AddDof(DISPLACEMENT_Y)
    node.AddDof(DISPLACEMENT_Z)

trilinos_fs_fluid_solver.AddDofs(p_model_part)

# creating a fluid solver object
fluid_solver = trilinos_fs_fluid_solver.IncompressibleFluidSolver(model_part, domain_size)
fluid_solver.laplacian_form = 1;
fluid_solver.predictor_corrector = True;
fluid_solver.vel_toll = 1e-3
fluid_solver.time_order = 2
fluid_solver.ReformDofAtEachIteration = False
fluid_solver.echo_level = 0

import PressureMultiLevelSolver
fluid_solver.pressure_linear_solver = PressureMultiLevelSolver.MultilevelLinearSolver(1e-4, 1000)


fluid_solver.Initialize()

# settings to be changed
Re = 100.0
nsteps = 100
output_step = 2

Dt = 1000.0 / Re
if(Dt > 0.1):
    Dt = 0.1
out = 0
for node in model_part.Nodes:
    node.SetSolutionStepValue(DENSITY, 0, 1.0)
    node.SetSolutionStepValue(VISCOSITY, 0, 1.0 / Re)


for elem in model_part.Elements:
    if (elem.Id < 200):
        elem.SetValue(SPLIT_ELEMENT, True)

ccc = ParallelFillCommunicator(model_part)
ccc.Execute()

print("*******************************************************************")

Comm = CreateCommunicator()
mesh_utility = TrilinosRefineMesh(model_part, Comm)
refine_on_reference = False
interpolate_internal_variables = False
mesh_utility.Local_Refine_Mesh(refine_on_reference, interpolate_internal_variables, domain_size)


print("-----------------------------------------------------------------")
ccc = ParallelFillCommunicator(model_part)
ccc.Execute()


mesh_utility.PrintDebugInfo()
aaa = BodyNormalCalculationUtils()
aaa.CalculateBodyNormals(model_part, 2)

mpi.world.barrier()

mesh_name = mpi.rank
gid_io.InitializeMesh(mesh_name);
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()
gid_io.Flush()

gid_io.InitializeResults(mesh_name, model_part.GetMesh())

BoxSize = 1;
bins_dynamic_mpi = BinsDynamicMpi(model_part, p_model_part, BoxSize)

print("----------------------------------------------------------")
bins_dynamic_mpi.MPIMultiSearchInRadiusTest()
# bins_dynamic_mpi.MPISingleSearchInRadiusTest()
# bins_dynamic_mpi.MultiSearchInRadiusTest()
print("----------------------------------------------------------")

# Count = 0
# for node in model_part.Nodes:
  # if (Count == 7):
        # print "----------------------------------------------------------"
        # bins_dynamic_mpi.MPISingleSearchInRadiusTest()
        # bins_dynamic_mpi.MultiSearchInRadiusTest()
        # print "----------------------------------------------------------"
  # Count = Count + 1;

# break_here

# for step in range(0,10):

    # time = Dt*step
    # model_part.CloneTimeStep(time)

    # if(mpi.rank == 0):
        # print time
    # print model_part.ProcessInfo()[TIME]

    # solving the fluid problem
    # if(step > 4):
        # fluid_solver.Solve()

# for node in model_part.Nodes:
# node.SetSolutionStepValue(PRESS_PROJ,0,zero);


    # print the results
    # if(out == output_step):
        # gid_io.WriteNodalResults(PARTITION_INDEX,model_part.Nodes,time,0)
        # gid_io.WriteNodalResults(NORMAL,model_part.Nodes,time,0)
        # gid_io.WriteNodalResults(PRESSURE,model_part.Nodes,time,0)
        # gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,time,0)
        # out = 0
    # out = out + 1

# for step in range(0,nsteps):

    # time = Dt*step
    # model_part.CloneTimeStep(time)

    # if(mpi.rank == 0):
        # print time
    # print model_part.ProcessInfo()[TIME]

    # solving the fluid problem
    # if(step > 4):
        # fluid_solver.Solve()

# for node in model_part.Nodes:
# node.SetSolutionStepValue(PRESS_PROJ,0,zero);


    # print the results
    # if(out == output_step):
        # gid_io.WriteNodalResults(PARTITION_INDEX,model_part.Nodes,time,0)
        # gid_io.WriteNodalResults(NORMAL,model_part.Nodes,time,0)
        # gid_io.WriteNodalResults(PRESSURE,model_part.Nodes,time,0)
        # gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,time,0)
        # out = 0
    # out = out + 1

# gid_io.FinalizeResults()