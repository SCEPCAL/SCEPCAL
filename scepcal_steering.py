from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import mm, GeV, MeV

SIM = DD4hepSimulation()

SIM.macroFile = ""
SIM.compactFile = ['install/share/compact/SCEPCAL.xml']

SIM.crossingAngleBoost = 0.0 # Lorentz boost for the crossing angle, in radian!
SIM.enableDetailedShowerMode = False
SIM.enableG4GPS = False
SIM.enableG4Gun = False
SIM.enableGun = False

# .stdhep, .slcio, .HEPEvt, .hepevt, .pairs, .hepmc, .hepmc.gz, .hepmc.xz, 
# .hepmc.bz2, .hepmc3, .hepmc3.gz, .hepmc3.xz, .hepmc3.bz2, .hepmc3.tree.root
SIM.inputFiles = ['wzp6_ee_ZZ_test_ecm240_10.stdhep']

SIM.numberOfEvents = 10

# .slcio, edm4hep.root, .root
SIM.outputFile = "wzp6_ee_ZZ_test_ecm240_test10_FG.slcio"
SIM.physicsList = None

# (1) VERBOSE, DEBUG, INFO, WARNING, ERROR, FATAL, ALWAYS (7)
SIM.printLevel = 3
SIM.runType = "batch"
SIM.skipNEvents = 0
SIM.steeringFile = None

SIM.vertexOffset = [0.0, 0.0, 0.0, 0.0]
SIM.vertexSigma = [0.0, 0.0, 0.0, 0.0]

SIM.action.calo = "ScepcalSDAction"
SIM.action.calorimeterSDTypes = ['ScepcalSD']
SIM.action.mapActions = { 'Scepcal' : "ScepcalSDAction" }
# SIM.action.mapActions = {}

SIM.action.tracker = None
# ('Geant4TrackerWeightedAction', {'HitPositionCombination': 2, 'CollectSingleDeposits': False})
SIM.action.trackerSDTypes = []

# SIM.field.delta_chord = 0.25
# SIM.field.delta_intersection = 0.001
# SIM.field.delta_one_step = 0.01
# SIM.field.eps_max = 0.001
# SIM.field.eps_min = 5e-05
# SIM.field.equation = "Mag_UsualEqRhs"
# SIM.field.largest_step = 10000.0
# SIM.field.min_chord_step = 0.01
# SIM.field.stepper = "ClassicalRK4"

SIM.filter.calo = ""
SIM.filter.filters = {
                        'geantino': {'name': 'GeantinoRejectFilter/GeantinoRejector', 'parameter': {}}, 
                        'edep1kev': {'name': 'EnergyDepositMinimumCut', 'parameter': {'Cut': 0.001}}, 
                        'edep0': {'name': 'EnergyDepositMinimumCut/Cut0', 'parameter': {'Cut': 0.0}}
                      }
SIM.filter.mapDetFilter = {}
SIM.filter.tracker = "edep1kev"

SIM.geometry.dumpGDML = ""
SIM.geometry.dumpHierarchy = 0
SIM.geometry.enableDebugElements = False
SIM.geometry.enableDebugMaterials = False
SIM.geometry.enableDebugPlacements = False
SIM.geometry.enableDebugReflections = False
SIM.geometry.enableDebugRegions = False
SIM.geometry.enableDebugShapes = False
SIM.geometry.enableDebugSurfaces = False
SIM.geometry.enableDebugVolumes = False
SIM.geometry.enablePrintPlacements = False
SIM.geometry.enablePrintSensitives = False

SIM.guineapig.particlesPerEvent = "-1"

SIM.hepmc3.Flow1 = "flow1"
SIM.hepmc3.Flow2 = "flow2"
SIM.hepmc3.useHepMC3 = True

SIM.inputConfig.userInputPlugin = []

SIM.lcio.mcParticleCollectionName = "MCParticle"
SIM.meta.eventNumberOffset = 0
## Event parameters to write in every event. Use C/F/I ids to specify parameter type. E.g parameterName/F=0.42 to set a float parameter
SIM.meta.eventParameters = []
SIM.meta.runNumberOffset = 0

SIM.output.geometry = 4
SIM.output.inputStage = 4
SIM.output.kernel = 4
SIM.output.part = 7
SIM.output.random = 6

SIM.outputConfig.forceDD4HEP = False
SIM.outputConfig.forceEDM4HEP = False
SIM.outputConfig.forceLCIO = False
SIM.outputConfig.userOutputPlugin = None

SIM.part.enableDetailedHitsAndParticleInfo = False
SIM.part.keepAllParticles = False
## Minimal distance between particle vertex and endpoint of parent after
##     which the vertexIsNotEndpointOfParent flag is set
SIM.part.minDistToParentVertex = 2.2e-14
## MinimalKineticEnergy to store particles created in the tracking region
SIM.part.minimalKineticEnergy = 1.0
SIM.part.printEndTracking = True
SIM.part.printStartTracking = True
SIM.part.saveProcesses = ['Decay']
# SIM.part.userParticleHandler = "Geant4TCUserParticleHandler"
SIM.part.userParticleHandler = None

SIM.physics.decays = False
SIM.physics.list = "FTFP_BERT"
SIM.physics.pdgfile = None
SIM.physics.rangecut = 0.7
SIM.physics.rejectPDGs = {1, 2, 3, 4, 5, 6, 3201, 3203, 4101, 4103, 21, 23, 24, 5401, 25, 
                          2203, 5403, 3101, 3103, 4403, 2101, 5301, 2103, 5303, 4301, 1103, 
                          4303, 5201, 5203, 3303, 4201, 4203, 5101, 5103, 5503}
SIM.physics.zeroTimePDGs = {17, 13, 15}

SIM.random.enableEventSeed = False
SIM.random.file = None
SIM.random.luxury = 1
SIM.random.replace_gRandom = True
SIM.random.seed = None
SIM.random.type = None

SIM.ui.commandsConfigure = []
SIM.ui.commandsInitialize = []
SIM.ui.commandsPostRun = []
SIM.ui.commandsPreRun = []
SIM.ui.commandsTerminate = []
