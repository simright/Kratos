from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
import time as timer
import random
import sys

#from exaqute.ExaquteTaskPyCOMPSs import *
#from exaqute.ExaquteTaskHyperLoom import *

#@ExaquteTask()
def f(serializer):

    before_load = timer.time()
    current_model = KratosMultiphysics.Model()

    #here doing the loading
    serializer.Load("Model", current_model)
    del current_model
    after_load = timer.time()
    print("Kratos loading from Serializer time = ",after_load-before_load)

#@ExaquteTask(returns=1)
def create_object():
    file_name = "restart_file"

    before_read = timer.time()

    model_part_name = "MainRestart"
    current_model = KratosMultiphysics.Model()
    model_part = current_model.CreateModelPart(model_part_name)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
    model_part_io = KratosMultiphysics.ModelPartIO("tube_3d_simple")
    model_part_io.ReadModelPart(model_part)

    model_part.Nodes[100].SetSolutionStepValue(KratosMultiphysics.VISCOSITY,0,100.0)

    after_read = timer.time()
    print("reading time = ",after_read-before_read)


    before_save = timer.time()

    serializer = KratosMultiphysics.StreamSerializer()
    serializer.Save("Model", current_model)
    
    after_save = timer.time()
    print("Kratos saving to Serializer time = ",after_save-before_save)
    return serializer

if __name__ == "__main__":
    n_nodes = sys.argv[1]
    ncores = 12 + (int(n_nodes) - 1) * 24
    nrepetitions = 4
    ncores = 1
    init_time = timer.time()
    object_to_share = create_object()
    #barrier()
    end_creation_time = timer.time() 

    for i in range(ncores * nrepetitions): 
        f(object_to_share) 
    
    #barrier()
    end_execution_time = timer.time()

    print("Amount of calls ", ncores, " with ", nrepetitions, " repetitions")
    print("Object creation time: ", end_creation_time - init_time)
    print(nrepetitions, " iterations time with ", n_nodes, " computational nodes (", ncores, " cores): ", end_execution_time - end_creation_time)
    print("Mean time per iteration: ", (end_execution_time - end_creation_time) / 4)