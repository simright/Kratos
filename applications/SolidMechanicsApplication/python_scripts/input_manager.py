from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
#import kratos core and applications
import KratosMultiphysics

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

#Base class to develop other solvers
class InputManager(object):
    """The base class for solid mechanic input parameters and materials.

    This class provides functions for seting parts to input parameters,
    and material parameters

    """
    def __init__(self, input_file):

        if os.path.isfile(input_file):
            parameters_file = open(input_file,'r')
            self.parameters = KratosMultiphysics.Parameters(parameters_file.read())
        else:
            raise Exception("Input file "+input_file+" does not exist")

        if(self.parameters.Has("input_settings")):
            custom_settings = self.parameters["input_settings"]
            self._set_custom_settings(custom_settings)
            
    #
    def GetProjectParameters(self):

        self.project_parameters = self.parameters
        
        if(self.parameters.Has("input_settings")):
            self._set_input_parts()

        #print(self.project_parameters.PrettyPrintJsonString())
        return self.project_parameters

    #
    def HasMaterialFile(self):

        if(self.parameters.Has("input_settings")):
            if(self.settings["materials_file_name"].GetString() != "None"):
                return True

        return False
    
    #
    def GetMaterialParameters(self):

        if os.path.isfile("Materials.json"):
            materials_file = open("Materials.json",'r')
            self.material_parameters = KratosMultiphysics.Parameters(materials_file.read())
            
        if(self.parameters.Has("input_settings")):
            self._set_material_parts()
            
        #print(self.material_parameters.PrettyPrintJsonString())
            
        return self.material_parameters
        

        
    #### Input manager internal methods ####

    #
    def _set_custom_settings(self, custom_settings):
                
        default_settings = KratosMultiphysics.Parameters("""
        {
            "parameters_file_name": "None",
            "model_parts":{
                "bodies_list"          : [],
                "domain_parts_list"    : [],
                "processes_parts_list" : []
            },
            "processes_parts":{
                "constraint_processes_parts" : [],
                "loads_processes_parts" : [],
                "output_processes_parts" : [],
                "check_processes_parts" : []
            },
            "materials_file_name" : "None",
            "material":{
                "material_ids"   : [],
                "material_parts" : []
            }
        }
        """)
        
        # Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

    #
    def _set_input_parts(self):
       
        if(self.settings["parameters_file_name"].GetString() != "None"):

            parameters_file = self.settings["parameters_file_name"].GetString()
            if os.path.isfile(parameters_file):
                parameters_file = open(parameters_file,'r')
                self.project_parameters = KratosMultiphysics.Parameters(parameters_file.read())
            else:
                raise Exception("Parameters file "+parameters_file+" does not exist")

            # model
            self._set_model_parts()
            
            # processes
            self._set_processes_parts()
   
            
    #
    def _set_model_parts(self):

        if( self.project_parameters["model_settings"].Has("bodies_list") == False ):
            if( self.settings["model_parts"].Has("bodies_list") ):
                if( self.settings["model_parts"]["bodies_list"].size() > 0 ):
                    self.project_parameters["model_settings"].AddValue("bodies_list", self.settings["model_parts"]["bodies_list"])
                
        if( self.project_parameters["model_settings"].Has("domains_parts_list") == False ):
            if( self.settings["model_parts"].Has("domain_parts_list") ):
                if( self.settings["model_parts"]["domain_parts_list"].size() > 0 ):
                    self.project_parameters["model_settings"].AddValue("domain_parts_list", self.settings["model_parts"]["domain_parts_list"])
                
        if( self.project_parameters["model_settings"].Has("processes_parts_list") == False ):
            if( self.settings["model_parts"].Has("processes_parts_list") ):
                if( self.settings["model_parts"]["processes_parts_list"].size() > 0 ):
                    self.project_parameters["model_settings"].AddValue("processes_parts_list", self.settings["model_parts"]["processes_parts_list"])

    #
    def _set_processes_parts(self):

        if( self.project_parameters.Has("constraints_process_list") ):
            processes_list = "constraints_process_list"
            parts_list     = "constraint_processes_parts"

            self._set_processes_type(processes_list,parts_list)
            
        if( self.project_parameters.Has("loads_process_list") ):
            processes_list = "loads_process_list"
            parts_list     = "loads_processes_parts"

            self._set_processes_type(processes_list,parts_list)

        if( self.project_parameters.Has("output_process_list") ):
            processes_list = "output_process_list"
            parts_list     = "output_processes_parts"

            self._set_processes_type(processes_list,parts_list)
            
        if( self.project_parameters.Has("check_process_list") ):
            processes_list = "check_process_list"
            parts_list     = "check_processes_parts"

            self._set_processes_type(processes_list,parts_list)
        
    #
    def _set_processes_type(self,processes_list,parts_list):
        
        parts = self.settings["processes_parts"][parts_list]
        if(parts.size() > 0):
            if(self.project_parameters[processes_list].size() == parts.size()):
                for i in range(0,parts.size()):
                    part_name = self.settings["model_parts"]["processes_parts_list"][parts[i].GetInt()].GetString()
                    process   = self.project_parameters[processes_list][i]
                    if( process.Has("Parameters") ):
                        process_parameters = process["Parameters"]
                        if( process_parameters.Has("model_part_name") ):
                            self.project_parameters[processes_list][i]["Parameters"]["model_part_name"].SetString(part_name)
                        else:
                            self.project_parameters[processes_list][i]["Parameters"].AddEmptyValue("model_part_name").SetString(part_name)
                    else:
                        if( process.Has("model_part_name") ):
                            self.project_parameters[processes_list][i]["model_part_name"].SetString(part_name)
                        else:
                            self.project_parameters[processes_list][i].AddEmptyValue("model_part_name").SetString(part_name)

            else:
                raise Exception(processes_list+" and "+parts_list+" do not have the same size")

    
    #
    def _set_material_parts(self):
        
        if(self.settings["materials_file_name"].GetString() != "None"):

            materials_file_name = self.settings["materials_file_name"].GetString()
            if os.path.isfile(materials_file_name):
                materials_file = open(materials_file_name,'r')
                self.material_parameters = KratosMultiphysics.Parameters(materials_file.read())
            else:
                raise Exception("Materials file "+materials_file_name+" does not exist")

            
            self._set_material_items()
            
    #
    def _set_material_items(self):

        materials_list = KratosMultiphysics.Parameters("""
        {
        }
        """)
        
        material = self.settings["material"]
        size = material["material_ids"].size()
        if( size > 0):
            materials_list.AddEmptyValue("material_models_list").SetVector(KratosMultiphysics.Vector(size))
            if( material["material_parts"].size() == size ):
                materials = self.material_parameters["material_models_list"]
                if(materials.size() >= size ):
                    for i in range(0,materials.size()):
                        material_parameters = materials[i]["Parameters"]
                        part_name = material["material_parts"][i].GetString()
                        if( material_parameters["properties_id"].GetInt() == material["material_ids"][i].GetInt() ):
                                                        
                            if( material_parameters.Has("model_part_name") ):
                                materials[i]["Parameters"]["model_part_name"].SetString(part_name)
                            else:
                                materials[i]["Parameters"].AddEmptyValue("model_part_name").SetString(part_name)

                            materials_list["material_models_list"].__setitem__(i, materials.__getitem__(i))

                                            
                else:
                    raise Exception("material_models_list size is too small")

                
                self.material_parameters = materials_list

            else:
                raise Exception("material_ids and material_parts do not have the same size")
