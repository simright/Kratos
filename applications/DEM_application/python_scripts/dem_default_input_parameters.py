from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import KratosMultiphysics
KratosMultiphysics.CheckForPreviousImport()

def GetDefaultInputParameters():
    
    default_settings = KratosMultiphysics.Parameters("""
        {
            "Dimension"                        : 3,
            "PeriodicDomainOption"             : false,
            "BoundingBoxOption"                : false,
            "AutomaticBoundingBoxOption"       : false,
            "BoundingBoxEnlargementFactor"     : 1.0,
            "BoundingBoxStartTime"             : 0.0,
            "BoundingBoxStopTime"              : 1000.0,
            "BoundingBoxMaxX"                  : 10.0,
            "BoundingBoxMaxY"                  : 10.0,
            "BoundingBoxMaxZ"                  : 10.0,
            "BoundingBoxMinX"                  : -10.0,
            "BoundingBoxMinY"                  : -10.0,
            "BoundingBoxMinZ"                  : -10.0,

            "dem_inlet_option"                 : true,
            "GravityX"                         : 0.0,
            "GravityY"                         : 0.0,
            "GravityZ"                         : -9.81,

            "EnergyCalculationOption"          : false,
            "PotentialEnergyReferencePointX"   : 0.0,
            "PotentialEnergyReferencePointY"   : 0.0,
            "PotentialEnergyReferencePointZ"   : 0.0,
           
            "VelocityTrapOption"               : false,
            "VelocityTrapMaxX"                 : 0.0,
            "VelocityTrapMaxY"                 : 0.0,
            "VelocityTrapMaxZ"                 : 0.0,
            "VelocityTrapMinX"                 : 0.0,
            "VelocityTrapMinY"                 : 0.0,
            "VelocityTrapMinZ"                 : 0.0,
            "RotationOption"                   : true,
            "CleanIndentationsOption"          : false,
            "RemoveBallsInEmbeddedOption"      : false,

            "DeltaOption"                      : "Absolute",
            "SearchTolerance"                  : 0.0,
            "CoordinationNumber"               : 10,
            "AmplifiedSearchRadiusExtension"   : 0.0,
            "MaxAmplificationRatioOfSearchRadius" : 10,
            "ModelDataInfo"                    : false,
            "VirtualMassCoefficient"           : 1.0,
            "RollingFrictionOption"            : false,
            "ComputeStressTensorOption"        : false,
            "GlobalDamping"                    : 0.0,
            "PoissonEffectOption"              : true,
            "ShearStrainParallelToBondOption"  : true,
            "DontSearchUntilFailure"           : false,
            "ContactMeshOption"                : false,
            "OutputFileType"                   : "Binary",
            "Multifile"                        : "multiple_files",
            "TestType"                         : "None",
            "ElementType"                      : "SphericPartDEMElement3D",

            "TranslationalIntegrationScheme"   : "Symplectic_Euler",
            "RotationalIntegrationScheme"      : "Direct_Integration",
            "AutomaticTimestep"                : false,
            "DeltaTimeSafetyFactor"            : 1.0,
            "MaxTimeStep"                      : 5e-5,
            "FinalTime"                        : 0.05,
            "ControlTime"                      : 4.0,
            "NeighbourSearchFrequency"         : 50,   
            
            "TestType"                         : false,
            "ConfinementPressure"              : 0.0,
            "LoadingVelocityTop"               : -0.10,
            "LoadingVelocityBot"               : 0.0,
            "MeshType"                         : "1",
            "MeshPath"                         : "0",
            "SpecimenLength"                   : 0.3,
            "SpecimenDiameter"                 : 0.15,
            "MeasuringSurface"                 : 0.01767145867644375,                                             

            "GraphExportFreq"                  : 1e-3,
            "VelTrapGraphExportFreq"           : 1e-3,
            "OutputTimeStep"                   : 1e-2,
            "PostBoundingBox"                  : false,
            "PostDisplacement"                 : false,
            "PostVelocity"                     : true,
            "PostTotalForces"                  : false,
            "PostRigidElementForces"           : false,
            "PostSkinSphere"                   : false,
            "PostPoissonRatio"                 : false,
            "PostRadius"                       : false,
            "PostAngularVelocity"              : false,
            "PostParticleMoment"               : false,
            "PostEulerAngles"                  : false,
            "PostRollingResistanceMoment"      : false,
            "PostCharacteristicLength"         : false,
            "PostElasticForces"                : false,
            "PostContactForces"                : false,
            "PostTangentialElasticForces"      : false,
            "PostShearStress"                  : false,
            "PostReactions"                    : false,
            "PostPressure"                     : false,
            "PostNonDimensionalVolumeWear"     : false,
            "PostNodalArea"                    : false,
            "PostStressStrainOption"           : false,
            "PostContactSigma"                 : false,
            "PostContactTau"                   : false,
            "PostLocalContactForce"            : false,
            "PostFailureCriterionState"        : false,
            "PostContactFailureId"             : false,
            "PostMeanContactArea"              : false,
            "PostRHS"                          : false,
            "PostDampForces"                   : false,
            "PostAppliedForces"                : false,
            "PostGroupId"                      : false,
            "PostExportId"                     : false,
                      
            "IntegrationScheme"                : "deprecated_key_since_6_december_2017",
            
            "problem_name" : "dummy_name.Provide_a_real_one"
            }""")
            
    return default_settings
