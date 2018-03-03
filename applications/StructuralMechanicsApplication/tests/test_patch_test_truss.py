from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

from math import sqrt, sin, cos, pi, exp, atan

class TestTruss3D2N(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def _add_dofs(self,mp):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)


    def _apply_material_properties(self,mp,dim):
        #define properties
        mp.GetProperties()[0].SetValue(KratosMultiphysics.YOUNG_MODULUS,210e9)
        mp.GetProperties()[0].SetValue(KratosMultiphysics.DENSITY,7850)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.CROSS_AREA,0.01)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.TRUSS_PRESTRESS_PK2,0)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.RAYLEIGH_ALPHA,0)
        mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.RAYLEIGH_BETA,0)

        g = [0,0,0]
        mp.GetProperties()[0].SetValue(KratosMultiphysics.VOLUME_ACCELERATION,g)

        cl = StructuralMechanicsApplication.LinearElastic3DLaw()
        mp.GetProperties()[0].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)


    def _apply_BCs(self,mp,which_dof):
        if (which_dof == 'xyz'):
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_X, True, mp.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Y, True, mp.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Z, True, mp.Nodes)
        if (which_dof == 'xz'):
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_X, True, mp.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Z, True, mp.Nodes)
        if (which_dof == 'yz'):
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Y, True, mp.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.DISPLACEMENT_Z, True, mp.Nodes)

    def _apply_Neumann_BCs(self,mp,which_dof,load_size_dir):
        if(which_dof == 'y'):
            KratosMultiphysics.VariableUtils().SetScalarVar(StructuralMechanicsApplication.
                POINT_LOAD_Y, load_size_dir, mp.Nodes)
            # for node in mp.Nodes:
            #     node.SetSolutionStepValue(StructuralMechanicsApplication.
            #     POINT_LOAD_Y,0,load_size_dir)
        if(which_dof == 'x'):
            KratosMultiphysics.VariableUtils().SetScalarVar(StructuralMechanicsApplication.
                POINT_LOAD_X, load_size_dir, mp.Nodes)
            # for node in mp.Nodes:
            #     node.SetSolutionStepValue(StructuralMechanicsApplication.
            #     POINT_LOAD_X,0,load_size_dir)

    def _solve_linear(self,mp):
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

        compute_reactions = True
        reform_step_dofs = True
        calculate_norm_dx = False
        move_mesh_flag = True
        strategy = KratosMultiphysics.ResidualBasedLinearStrategy(mp,
                                                                scheme,
                                                                linear_solver,
                                                                builder_and_solver,
                                                                compute_reactions,
                                                                reform_step_dofs,
                                                                calculate_norm_dx,
                                                                move_mesh_flag)
        strategy.SetEchoLevel(0)

        strategy.Check()
        strategy.Solve()

    def _solve_nonlinear(self,mp):
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        convergence_criterion = KratosMultiphysics.ResidualCriteria(1e-8,1e-8)
        convergence_criterion.SetEchoLevel(0)

        max_iters = 1000
        compute_reactions = True
        reform_step_dofs = True
        move_mesh_flag = True
        strategy = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(mp,
                                                                scheme,
                                                                linear_solver,
                                                                convergence_criterion,
                                                                builder_and_solver,
                                                                max_iters,
                                                                compute_reactions,
                                                                reform_step_dofs,
                                                                move_mesh_flag)
        strategy.SetEchoLevel(0)

        strategy.Check()
        strategy.Solve()

    def _solve_dynamic(self,mp):
        #define a minimal newton raphson solver
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedBossakDisplacementScheme(0.00)
        convergence_criterion = KratosMultiphysics.ResidualCriteria(1e-8,1e-8)
        convergence_criterion.SetEchoLevel(0)

        max_iters = 1000
        compute_reactions = True
        reform_step_dofs = True
        move_mesh_flag = True
        strategy = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(mp,
                                                                scheme,
                                                                linear_solver,
                                                                convergence_criterion,
                                                                builder_and_solver,
                                                                max_iters,
                                                                compute_reactions,
                                                                reform_step_dofs,
                                                                move_mesh_flag)
        strategy.SetEchoLevel(0)

        strategy.Check()
        strategy.Solve()


    def _check_results_linear(self,mp):
        #1.) check displacement result
        displacement_nodes = [mp.Nodes[1].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT),mp.Nodes[2].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT),mp.Nodes[3].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT)]
        for i in range(3):
            k = 0.00
            if (i == 1): k = -0.0026619856874997507
            self.assertAlmostEqual(displacement_nodes[1][i], k)
            self.assertAlmostEqual(displacement_nodes[0][i], 0)
            self.assertAlmostEqual(displacement_nodes[2][i], 0)

        #2.) check reactions
        reaction_nodes = [mp.Nodes[1].GetSolutionStepValue(
            KratosMultiphysics.REACTION),mp.Nodes[2].GetSolutionStepValue(
            KratosMultiphysics.REACTION),mp.Nodes[3].GetSolutionStepValue(
            KratosMultiphysics.REACTION)]

        Force_y = 1000000
        self.assertAlmostEqual(reaction_nodes[0][0], Force_y)
        self.assertAlmostEqual(reaction_nodes[2][0], -Force_y)
        self.assertAlmostEqual(reaction_nodes[0][1], Force_y/2)
        self.assertAlmostEqual(reaction_nodes[2][1], Force_y/2)

    def _check_results_nonlinear(self,mp,timestep,Force_i):
        ##node1
        node_temp = mp.Nodes[1]
        disp_temp = node_temp.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
        reac_temp = node_temp.GetSolutionStepValue(KratosMultiphysics.REACTION)
        #disp_y
        self.assertAlmostEqual(disp_temp[1], 0)
        #reaction_y
        reaction_y_node1 = Force_i*(-1)
        self.assertAlmostEqual(reac_temp[1],reaction_y_node1,6)
        #reaction_x
        reaction_x_node1 = [741464.9276510741,1485888.977636112,2233316.9009164227,
        2983794.615716884,3737369.2509119534,4494089.191516033,5254004.126397622,
        6017165.0983619075,6783624.556737246,7553436.412628534,8326656.0969987875,
        9103340.621764094,9883548.644087134,10667340.53408764,11454778.446184492,
        12245926.394319195,13040850.331311818,13839618.232645638,14642300.18497437,
        15448968.47968232,16259697.711865114,17074564.885111626,17893649.522510014,
        18717033.784337185,19544802.592936598,20377043.765315644,21213848.154068694,
        22055309.79726978,22901526.07803619,23752597.894554257,24608629.841397412,
        25469730.40309285,26336012.160943132,27207592.014241144,28084591.41712631,
        28967136.632448073,29855359.004160944,30749395.24993747,31649387.775854316,
        32555485.015236698,33467841.79396308,34386619.72480394,35311987.633672595,
        36244122.02100001,37183207.5618443,38129437.64879564,39083014.9822317,
        40044152.21308944,41013072.64397469,41990010.99523344,42975214.24348978,
        43968942.54123137,44971470.22722923,45983086.93901924,47004098.840344116,
        48034829.97843962,49075623.78836585,50126844.76436112,51188880.321473464,
        52262142.87466787,53347072.16731208,54444137.88663925,55553842.61067612,
        56676725.13951358,57813364.2740804,58964383.118212394,60130453.99549516,
        61312304.09186774,62510721.95949238,63726565.048343465,64960768.47141639,
        66214355.26005781,67488448.43149039,68784285.27627705,70103234.38657254,
        71446816.09699601,72816727.21376368,74214871.18642414,75643395.26295514,
        77104736.71277626,78601680.98043492,80137435.76669273,81715726.72014935,
        83340922.98755075,85018204.87282823,86753792.28015186,88555263.27634439,
        90432010.47500862,92395916.01811945,94462388.67855237,96652033.38549507,
        98993500.26523624,101528726.98334643,104323616.16359182,107493197.43582612,
        111276440.23647068,116390127.39236663,-62782528.388332605,-63351316.30823133,
        -63919034.598836,-64485690.945303164,-65051292.93836311,-65615848.075949684,
        -66179363.76479947,-66741847.3220098,-67303305.9765681,-67863746.8708426,
        -68423177.06204486,-68981603.5236578,-69539033.1468354,-70095472.74176757,
        -70650929.0390236,-71205408.69085957,-71758918.27250087,-72311464.28340018,
        -72863053.1484657,-73413691.21926463,-73963384.77520159,-74512140.02467461,
        -75059963.10620539,]
        self.assertAlmostEqual(reac_temp[0],reaction_x_node1[timestep])

        ##node2
        node_temp = mp.Nodes[2]
        disp_temp = node_temp.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
        reac_temp = node_temp.GetSolutionStepValue(KratosMultiphysics.REACTION)
        load_temp = node_temp.GetSolutionStepValue(StructuralMechanicsApplication.POINT_LOAD_Y)
        #pointLoad
        self.assertAlmostEqual(load_temp,Force_i)
        #reaction_x
        self.assertAlmostEqual(reac_temp[0],reaction_x_node1[timestep]*(-1))
        #displacement_y
        EA = 210e9*0.01
        L = sqrt(4+1)
        L3 = L*L*L
        P_i = ((EA/(2*L3))*(disp_temp*disp_temp +2*1*disp_temp)*(disp_temp+1))
        self.assertAlmostEqual(P_i,Force_i,1)

    def _check_pre_stress_output(self,mp,force,tolerance=9):
        for element in mp.Elements:
            out = element.CalculateOnIntegrationPoints(KratosMultiphysics.FORCE,mp.ProcessInfo)
            self.assertAlmostEqual(out[0][0],force,tolerance)
            self.assertAlmostEqual(out[0][1],0.00)
            self.assertAlmostEqual(out[0][2],0.00)


    def _check_results_dynamic(self,mp,time_i):

        #analaytical free-vibration node 3
        we1 = 7917.25
        we2 = 19113.94
        y1 = 1.4142*2.874e-5
        y2 = -1.4142*4.93107e-6
        test_disp_temp = y1*(1-cos(we1*time_i))-y2*(1-cos(we2*time_i))
        simulated_disp_temp = mp.Nodes[3].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT_X)

        self.assertAlmostEqual(simulated_disp_temp, test_disp_temp,6)

        #analaytical free-vibration node 2
        we1 = 7917.25
        we2 = 19113.94
        y1 = 1.000*2.874e-5
        y2 = 1.000*4.93107e-6
        test_disp_temp = y1*(1-cos(we1*time_i))-y2*(1-cos(we2*time_i))
        simulated_disp_temp = mp.Nodes[2].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT_X)

        self.assertAlmostEqual(simulated_disp_temp, test_disp_temp,6)

    def _check_results_cable(self,mp,Force_X):
        
        disp_u_2 = mp.Nodes[2].GetSolutionStepValue(
        KratosMultiphysics.DISPLACEMENT_X)
        r_u_1 = mp.Nodes[1].GetSolutionStepValue(
        KratosMultiphysics.REACTION_X)
        r_u_3 = mp.Nodes[3].GetSolutionStepValue(
        KratosMultiphysics.REACTION_X)


        self.assertAlmostEqual(disp_u_2, 0.022296019142446475,6)
        self.assertAlmostEqual(r_u_1, -Force_X,6)
        self.assertAlmostEqual(r_u_3, 0.00 ,4)

    def _check_results_dynamic_explicit_nonlinear(self,mp,time_i,time_step):
             
        simulated_disp_temp = mp.Nodes[2].GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT_Y)
        test_disp_temp = [-0.02187643575439285, -0.06198215720163387,
         -0.12627003172959766, -0.19804924254370843, -0.2631923777969079,
          -0.313031207094566, -0.34361350795228707, -0.35373182403184056,
           -0.34338092709445905, -0.31319057182012966, -0.26486546314308357,
            -0.2024755542645808, -0.13407164406937494, -0.07231351774628614,
             -0.03221019822561301, -0.02538203496105414, -0.05392477771156322,
              -0.10902565726302998, -0.17575172783220822, -0.23975893050044683,
               -0.2911859747582321, -0.32489277095205155, -0.33890901757885206,
                -0.33289145786308927, -0.30738557023931556, -0.26404183431169526,
                 -0.20665783067442736, -0.14257695438336992, -0.08331538499080232,
                  -0.04276029659491523, -0.03225900075675591]
 
        self.assertAlmostEqual(simulated_disp_temp, test_disp_temp[time_step])

    def _set_and_fill_buffer(self,mp,buffer_size,delta_time):
        # Set buffer size
        mp.SetBufferSize(buffer_size)

        # Fill buffer
        time = mp.ProcessInfo[KratosMultiphysics.TIME]
        time = time - delta_time * (buffer_size)
        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
        for size in range(0, buffer_size):
            step = size - (buffer_size -1)
            mp.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
            time = time + delta_time
            #delta_time is computed from previous time in process_info
            mp.CloneTimeStep(time)

        mp.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = False

    def test_truss3D2N_linear(self):
        dim = 3
        mp = KratosMultiphysics.ModelPart("solid_part")
        self._add_variables(mp)
        self._apply_material_properties(mp,dim)

        #create nodes
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,2.0,1.0,0.0)
        mp.CreateNewNode(3,4.0,0.0,0.0)
        #add dofs
        self._add_dofs(mp)
        #create condition
        mp.CreateNewCondition("PointLoadCondition3D1N",1,[2],mp.GetProperties()[0])

        #create submodelparts for dirichlet boundary conditions
        bcs_xyz = mp.CreateSubModelPart("Dirichlet_XYZ")
        bcs_xyz.AddNodes([1,3])
        bcs_xz = mp.CreateSubModelPart("Dirichlet_XZ")
        bcs_xz.AddNodes([2])

        #create a submodalpart for neumann boundary conditions
        bcs_neumann = mp.CreateSubModelPart("PointLoad3D_neumann")
        bcs_neumann.AddNodes([2])
        bcs_neumann.AddConditions([1])

        #create Element
        mp.CreateNewElement("TrussLinearElement3D2N", 1, [1,2], mp.GetProperties()[0])
        mp.CreateNewElement("TrussLinearElement3D2N", 2, [2,3], mp.GetProperties()[0])

        #apply boundary conditions
        Force_Y = -1000000.00
        self._apply_BCs(bcs_xyz,'xyz')
        self._apply_BCs(bcs_xz,'xz')
        self._apply_Neumann_BCs(bcs_neumann,'y',Force_Y)

        #solve + compare
        self._solve_linear(mp)
        self._check_results_linear(mp)

    def test_truss3D2N_nonlinear(self):
        dim = 3
        mp = KratosMultiphysics.ModelPart("solid_part")
        self._add_variables(mp)
        self._apply_material_properties(mp,dim)

        #create nodes
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,2.0,1.0,0.0)
        #add dofs
        self._add_dofs(mp)
        #create condition
        mp.CreateNewCondition("PointLoadCondition3D1N",1,[2],mp.GetProperties()[0])

        #create submodelparts for dirichlet boundary conditions
        bcs_xyz = mp.CreateSubModelPart("Dirichlet_XYZ")
        bcs_xyz.AddNodes([1])
        bcs_xz = mp.CreateSubModelPart("Dirichlet_XZ")
        bcs_xz.AddNodes([2])

        #create a submodalpart for neumann boundary conditions
        bcs_neumann = mp.CreateSubModelPart("PointLoad3D_neumann")
        bcs_neumann.AddNodes([2])
        bcs_neumann.AddConditions([1])

        #create Element
        mp.CreateNewElement("TrussElement3D2N", 1, [1,2], mp.GetProperties()[0])
        #apply constant boundary conditions
        self._apply_BCs(bcs_xyz,'xyz')
        self._apply_BCs(bcs_xz,'xz')

        #incrementally increase load -> nonlinear case
        Force_y = -37000000
        time_start = 0.00
        time_end = 0.05
        time_delta = 0.01
        time_i = time_start
        time_step = 0
        while (time_i < time_end):

            time_i += time_delta
            #apply non-constant boundary conditions
            Force_i = Force_y*time_i
            self._apply_Neumann_BCs(bcs_neumann,'y',Force_i)
            #solve + compare
            self._solve_nonlinear(mp)
            self._check_results_nonlinear(mp,time_step,Force_i)
            time_step += 1

    def test_truss3D2N_prestress_nonlinear_fix(self):
            dim = 3
            mp = KratosMultiphysics.ModelPart("solid_part")
            self._add_variables(mp)
            self._apply_material_properties(mp,dim)
            mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.TRUSS_PRESTRESS_PK2,1000000.00)
            #create nodes
            mp.CreateNewNode(1,0.0,0.0,0.0)
            mp.CreateNewNode(2,2.0,0.0,0.0)
            #add dofs
            self._add_dofs(mp)

            #create submodelparts for dirichlet boundary conditions
            bcs_xyz = mp.CreateSubModelPart("Dirichlet_XYZ")
            bcs_xyz.AddNodes([1])
            bcs_xz = mp.CreateSubModelPart("Dirichlet_XZ")
            bcs_xz.AddNodes([2])


            #create Element
            mp.CreateNewElement("TrussElement3D2N", 1, [1,2], mp.GetProperties()[0])

            #apply constant boundary conditions
            self._apply_BCs(bcs_xyz,'xyz')
            self._apply_BCs(bcs_xz,'xz')
            self._solve_nonlinear(mp)
            self._check_pre_stress_output(mp,10000.0)

    def test_truss3D2N_prestress_nonlinear_free(self):
            dim = 3
            mp = KratosMultiphysics.ModelPart("solid_part")
            self._add_variables(mp)
            self._apply_material_properties(mp,dim)
            mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.TRUSS_PRESTRESS_PK2,1000000.00)
            #create nodes
            mp.CreateNewNode(1,0.0,0.0,0.0)
            mp.CreateNewNode(2,2.0,0.0,0.0)
            #add dofs
            self._add_dofs(mp)

            #create submodelparts for dirichlet boundary conditions
            bcs_xyz = mp.CreateSubModelPart("Dirichlet_XYZ")
            bcs_xyz.AddNodes([1])
            bcs_xz = mp.CreateSubModelPart("Dirichlet_XZ")
            bcs_xz.AddNodes([2])


            #create Element
            mp.CreateNewElement("TrussElement3D2N", 1, [1,2], mp.GetProperties()[0])

            #apply constant boundary conditions
            self._apply_BCs(bcs_xyz,'xyz')
            self._apply_BCs(bcs_xz,'yz')
            self._solve_nonlinear(mp)
            self._check_pre_stress_output(mp,0.0,6)


    def test_truss3D2N_prestress_linear_fix(self):
            dim = 3
            mp = KratosMultiphysics.ModelPart("solid_part")
            self._add_variables(mp)
            self._apply_material_properties(mp,dim)
            mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.TRUSS_PRESTRESS_PK2,1000000.00)
            #create nodes
            mp.CreateNewNode(1,0.0,0.0,0.0)
            mp.CreateNewNode(2,2.0,0.0,0.0)
            #add dofs
            self._add_dofs(mp)

            #create submodelparts for dirichlet boundary conditions
            bcs_xyz = mp.CreateSubModelPart("Dirichlet_XYZ")
            bcs_xyz.AddNodes([1])
            bcs_xz = mp.CreateSubModelPart("Dirichlet_XZ")
            bcs_xz.AddNodes([2])


            #create Element
            mp.CreateNewElement("TrussLinearElement3D2N", 1, [1,2], mp.GetProperties()[0])

            #apply constant boundary conditions
            self._apply_BCs(bcs_xyz,'xyz')
            self._apply_BCs(bcs_xz,'xz')
            self._solve_linear(mp)
            self._check_pre_stress_output(mp,10000.0)

    def test_truss3D2N_prestress_linear_free(self):
            dim = 3
            mp = KratosMultiphysics.ModelPart("solid_part")
            self._add_variables(mp)
            self._apply_material_properties(mp,dim)
            mp.GetProperties()[0].SetValue(StructuralMechanicsApplication.TRUSS_PRESTRESS_PK2,1000000.00)
            #create nodes
            mp.CreateNewNode(1,0.0,0.0,0.0)
            mp.CreateNewNode(2,2.0,0.0,0.0)
            #add dofs
            self._add_dofs(mp)

            #create submodelparts for dirichlet boundary conditions
            bcs_xyz = mp.CreateSubModelPart("Dirichlet_XYZ")
            bcs_xyz.AddNodes([1])
            bcs_xz = mp.CreateSubModelPart("Dirichlet_XZ")
            bcs_xz.AddNodes([2])


            #create Element
            mp.CreateNewElement("TrussLinearElement3D2N", 1, [1,2], mp.GetProperties()[0])

            #apply constant boundary conditions
            self._apply_BCs(bcs_xyz,'xyz')
            self._apply_BCs(bcs_xz,'yz')
            self._solve_linear(mp)
            self._check_pre_stress_output(mp,0.0)

    def test_truss3D2N_dynamic(self):
        dim = 3
        mp = KratosMultiphysics.ModelPart("solid_part")
        self._add_variables(mp)
        self._apply_material_properties(mp,dim)

        #create nodes
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,0.5,0.0,0.0)
        mp.CreateNewNode(3,1.0,0.0,0.0)
        #add dofs
        self._add_dofs(mp)
        #create condition
        mp.CreateNewCondition("PointLoadCondition3D1N",1,[3],mp.GetProperties()[0])
        #create submodelparts for dirichlet boundary conditions
        bcs_xyz = mp.CreateSubModelPart("Dirichlet_XYZ")
        bcs_xyz.AddNodes([1])
        bcs_yz = mp.CreateSubModelPart("Dirichlet_YZ")
        bcs_yz.AddNodes([2,3])
        #create a submodalpart for neumann boundary conditions
        bcs_neumann = mp.CreateSubModelPart("PointLoad3D_neumann")
        bcs_neumann.AddNodes([3])
        bcs_neumann.AddConditions([1])
        #create Elements
        mp.CreateNewElement("TrussElement3D2N", 1, [1,2], mp.GetProperties()[0])
        mp.CreateNewElement("TrussElement3D2N", 2, [2,3], mp.GetProperties()[0])
        #apply constant boundary conditions
        Force_X = 100000
        self._apply_BCs(bcs_xyz,'xyz')
        self._apply_BCs(bcs_yz,'yz')
        self._apply_Neumann_BCs(bcs_neumann,'x',Force_X)

        #loop over time
        time_start = 0.00
        time_end = 0.000002
        time_delta = 0.000001
        time_i = time_start
        time_step = 0
        self._set_and_fill_buffer(mp,2,time_delta)


        while (time_i <= time_end):

            time_i += time_delta
            mp.CloneTimeStep(time_i)
            #solve + compare
            self._solve_dynamic(mp)
            self._check_results_dynamic(mp,time_i)
            time_step += 1

    def test_truss3D2N_cable(self):
            dim = 3
            mp = KratosMultiphysics.ModelPart("solid_part")
            self._add_variables(mp)
            self._apply_material_properties(mp,dim)

            #create nodes
            mp.CreateNewNode(1,0.0,0.0,0.0)
            mp.CreateNewNode(2,0.5,0.0,0.0)
            mp.CreateNewNode(3,1.0,0.0,0.0)
            #add dofs
            self._add_dofs(mp)
            #create condition
            mp.CreateNewCondition("PointLoadCondition3D1N",1,[2],mp.GetProperties()[0])
            #create submodelparts for dirichlet boundary conditions
            bcs_xyz = mp.CreateSubModelPart("Dirichlet_XYZ")
            bcs_xyz.AddNodes([1,3])
            bcs_yz = mp.CreateSubModelPart("Dirichlet_YZ")
            bcs_yz.AddNodes([2])
            #create a submodalpart for neumann boundary conditions
            bcs_neumann = mp.CreateSubModelPart("PointLoad3D_neumann")
            bcs_neumann.AddNodes([2])
            bcs_neumann.AddConditions([1]) 
            #create Elements
            mp.CreateNewElement("CableElement3D2N", 1, [1,2], mp.GetProperties()[0])
            mp.CreateNewElement("CableElement3D2N", 2, [2,3], mp.GetProperties()[0])
            #apply constant boundary conditions
            Force_X = 100000000
            self._apply_BCs(bcs_xyz,'xyz')
            self._apply_BCs(bcs_yz,'yz')
            self._apply_Neumann_BCs(bcs_neumann,'x',Force_X)

            self._solve_nonlinear(mp)  
            self._check_results_cable(mp,Force_X)

    def test_truss3D2N_dynamic_explicit_nonlinear(self):
                dim = 3
                mp = KratosMultiphysics.ModelPart("solid_part")
                self._add_variables(mp)
                _add_explicit_variables(mp)
                self._apply_material_properties(mp,dim)

                #create nodes
                mp.CreateNewNode(1,0.0,0.0,0.0)
                mp.CreateNewNode(2,2.0,1.0,0.0)
                #add dofs
                self._add_dofs(mp)
                #create condition
                mp.CreateNewCondition("PointLoadCondition3D1N",1,[2],mp.GetProperties()[0])
                #create submodelparts for dirichlet boundary conditions
                bcs_xyz = mp.CreateSubModelPart("Dirichlet_XYZ")
                bcs_xyz.AddNodes([1])
                bcs_yz = mp.CreateSubModelPart("Dirichlet_XZ")
                bcs_yz.AddNodes([2])
                #create a submodalpart for neumann boundary conditions
                bcs_neumann = mp.CreateSubModelPart("PointLoad3D_neumann")
                bcs_neumann.AddNodes([2])
                bcs_neumann.AddConditions([1])
                #create Elementsdb
                mp.CreateNewElement("CrBeamElement3D2N", 1, [1,2], mp.GetProperties()[0])
                #apply constant boundary conditions
                Force_Y = -24000000
                self._apply_BCs(bcs_xyz,'xyz')
                self._apply_BCs(bcs_yz,'xz')
                self._apply_Neumann_BCs(bcs_neumann,'y',Force_Y)

                #loop over time
                time_start = 0.00
                time_end = 0.012
                time_delta = 0.0004
                time_i = time_start
                time_step = 0
                self._set_and_fill_buffer(mp,2,time_delta)

                strategy_expl = _create_dynamic_explicit_strategy(mp)
                while (time_i <= time_end):
                    time_i += time_delta
                    mp.CloneTimeStep(time_i)
                    #solve + compare
                    strategy_expl.Solve()
                    self._check_results_dynamic_explicit_nonlinear(mp,time_i,time_step)
                    time_step += 1
            
                



def _add_explicit_variables(mp):
    mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.MIDDLE_VELOCITY)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE_RESIDUAL)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.RESIDUAL_VECTOR)
    mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.MIDDLE_ANGULAR_VELOCITY)
    mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.NODAL_INERTIA)
    mp.AddNodalSolutionStepVariable(KratosMultiphysics.MOMENT_RESIDUAL)

def _create_dynamic_explicit_strategy(mp):
    scheme = StructuralMechanicsApplication.ExplicitCentralDifferencesScheme(0.00,0.00,0.00)

    strategy = StructuralMechanicsApplication.ExplicitStrategy(mp,
                                        scheme,0,0,1)
    strategy.SetEchoLevel(0)
    return strategy 


        
if __name__ == '__main__':
    KratosUnittest.main()

