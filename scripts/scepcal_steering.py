from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import mm, GeV, MeV, keV
from math import pi

SIM = DD4hepSimulation()

SIM.compactFile = ['install/share/compact/SCEPCAL.xml']

SIM.crossingAngleBoost = 0.0 # Lorentz boost for the crossing angle, in radian!
SIM.enableDetailedShowerMode = False
SIM.enableG4GPS = False
SIM.enableG4Gun = False
SIM.enableGun = False

SIM.inputFiles = ['examples/wzp6_ee_ZZ_test_ecm240_1k.stdhep']
SIM.outputFile = 'examples/wzp6_ee_ZZ_test_ecm240_n1_cut0_BEonly.root'

SIM.macroFile = ""
SIM.steeringFile = None

SIM.numberOfEvents = 1
SIM.skipNEvents = 0

# (1) VERBOSE, DEBUG, INFO, WARNING, ERROR, FATAL, ALWAYS (7)
SIM.printLevel = 7
SIM.runType = "batch"

SIM.vertexOffset = [0.0, 0.0, 0.0, 0.0]
SIM.vertexSigma = [0.0, 0.0, 0.0, 0.0]

# SIM.action.calo = "Geant4CalorimeterAction"
SIM.action.calo = "ScepcalSDAction"

SIM.action.calorimeterSDTypes = ['ScepcalSD']

# SIM.action.mapActions = { 'Scepcal' : "ScepcalSDAction" }
SIM.action.mapActions = {}

SIM.action.tracker = ('Geant4TrackerWeightedAction', 
                         {'HitPositionCombination': 2, 'CollectSingleDeposits': False})
SIM.action.trackerSDTypes = ['tracker']

################################################################################
## Configuration for the magnetic field (stepper) 
################################################################################
# SIM.field.delta_chord = 0.25
# SIM.field.delta_intersection = 0.001
# SIM.field.delta_one_step = 0.01
# SIM.field.eps_max = 0.001
# SIM.field.eps_min = 5e-05
# SIM.field.equation = "Mag_UsualEqRhs"
# SIM.field.largest_step = 10000.0
# SIM.field.min_chord_step = 0.01
# SIM.field.stepper = "ClassicalRK4"

SIM.filter.filters = {
                        'geantino': {'name': 'GeantinoRejectFilter/GeantinoRejector', 'parameter': {}}, 
                        'edep1keV': {'name': 'EnergyDepositMinimumCut/1keV', 'parameter': {'Cut': 1.0*keV }}, 
                        'edep0': {'name': 'EnergyDepositMinimumCut/Cut0', 'parameter': {'Cut': 0.0}}
                      }

SIM.filter.calo = "edep1keV"

SIM.filter.mapDetFilter = {}
SIM.filter.tracker = "edep0"

# SIM.geometry.dumpGDML = ""
# SIM.geometry.dumpHierarchy = 0
# SIM.geometry.enableDebugElements = False
# SIM.geometry.enableDebugMaterials = False
# SIM.geometry.enableDebugPlacements = False
# SIM.geometry.enableDebugReflections = False
# SIM.geometry.enableDebugRegions = False
# SIM.geometry.enableDebugShapes = False
# SIM.geometry.enableDebugSurfaces = False
# SIM.geometry.enableDebugVolumes = False
# SIM.geometry.enablePrintPlacements = False
# SIM.geometry.enablePrintSensitives = False

SIM.guineapig.particlesPerEvent = "-1"

SIM.hepmc3.Flow1 = "flow1"
SIM.hepmc3.Flow2 = "flow2"
SIM.hepmc3.useHepMC3 = False

SIM.inputConfig.userInputPlugin = []

SIM.lcio.mcParticleCollectionName = "MCParticle"
SIM.meta.eventNumberOffset = 0

SIM.meta.eventParameters = []
SIM.meta.runNumberOffset = 0

SIM.output.geometry = 7
SIM.output.inputStage = 7
SIM.output.kernel = 7
SIM.output.part = 1
SIM.output.random = 7

# SIM.outputConfig.forceDD4HEP = False
# SIM.outputConfig.forceEDM4HEP = False
# SIM.outputConfig.forceLCIO = False

def exampleUserPlugin(dd4hepSimulation):
     from DDG4 import EventAction, Kernel
     dd = dd4hepSimulation  
     evt_root = EventAction(Kernel(), 'SCEGeant4Output2ROOT/' + dd.outputFile, True)
     evt_root.HandleMCTruth = True
     evt_root.Control = True
     output = dd.outputFile
     if not dd.outputFile.endswith(dd.outputConfig.myExtension):
          output = dd.outputFile + dd.outputConfig.myExtension
     evt_root.Output = output
     evt_root.enableUI()
     Kernel().eventAction().add(evt_root)
     return None

# SIM.outputConfig.userOutputPlugin = exampleUserPlugin
# SIM.outputConfig.myExtension = '.root'

SIM.part.enableDetailedHitsAndParticleInfo = False
SIM.part.keepAllParticles = False
SIM.part.minDistToParentVertex = 2.2e-14
SIM.part.minimalKineticEnergy = 1.0
SIM.part.printEndTracking = False
SIM.part.printStartTracking = False
SIM.part.saveProcesses = ['Decay']

# SIM.part.userParticleHandler = "Geant4TCUserParticleHandler"
SIM.part.userParticleHandler = ''


#### Particle Gun

SIM.gun.multiplicity = 1
SIM.gun.position = (0, 0, 0)
# SIM.gun.direction = (0, 1, 0)

# SIM.gun.isotrop = False
# SIM.gun.distribution = 'uniform'

# SIM.gun.energy = 10*GeV
# SIM.gun.particle = "gamma"

# SIM.gun.momentumMin = 0.0
# SIM.gun.momentumMax = 10000.0

# SIM.gun.phiMin = 0
# SIM.gun.phiMax = 2*pi

# SIM.gun.thetaMin = pi/2
# SIM.gun.thetaMax = pi/2


#### Physics

SIM.physicsList = None
SIM.physics.decays = False
SIM.physics.list = "FTFP_BERT"
SIM.physics.pdgfile = None
SIM.physics.rangecut = None

SIM.physics.rejectPDGs = {1, 2, 3, 4, 5, 6, 3201, 3203, 4101, 4103, 21, 23, 24, 5401, 25, 
                          2203, 5403, 3101, 3103, 4403, 2101, 5301, 2103, 5303, 4301, 1103, 
                          4303, 5201, 5203, 3303, 4201, 4203, 5101, 5103, 5503}

SIM.physics.zeroTimePDGs = {17, 11, 13, 15}

def setupCerenkovScint(kernel):
     from DDG4 import PhysicsList
     seq = kernel.physicsList()

     scint = PhysicsList(kernel, 'Geant4ScintillationPhysics/ScintillationPhys')
     scint.VerboseLevel = 0
     scint.TrackSecondariesFirst = True
     scint.enableUI()
     seq.adopt(scint)

     cerenkov = PhysicsList(kernel, 'Geant4CerenkovPhysics/CerenkovPhys')
     cerenkov.VerboseLevel = 0
     cerenkov.MaxNumPhotonsPerStep = 10
     cerenkov.MaxBetaChangePerStep = 10.0
     cerenkov.TrackSecondariesFirst = True
     cerenkov.enableUI()
     seq.adopt(cerenkov)

     ph = PhysicsList(kernel, 'Geant4OpticalPhotonPhysics/OpticalGammaPhys')
     ph.addParticleConstructor('G4OpticalPhoton')
     ph.VerboseLevel = 0
     ph.enableUI()
     seq.adopt(ph)

     return None

# SIM.physics.setupUserPhysics(setupCerenkovScint)


#### RNG

SIM.random.enableEventSeed = False
SIM.random.file = None
SIM.random.luxury = 1
SIM.random.replace_gRandom = True

# SIM.random.seed = None
SIM.random.seed = 2126136508

SIM.random.type = None