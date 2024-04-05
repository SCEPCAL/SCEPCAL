from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import mm, GeV, MeV
from math import pi

SIM = DD4hepSimulation()

SIM.compactFile = ['install/share/compact/SCEPCAL.xml']

# ddsim --steeringFile scepcal_steering.py -G --gun.position="0.,0.,0." --gun.direction "0 1 0" --gun.energy "1*GeV" --gun.particle="gamma" --gun.distribution=uniform -N 1 -O gamma_y0cm_1GeV_no_opt.root


SIM.crossingAngleBoost = 0.0 # Lorentz boost for the crossing angle, in radian!
SIM.enableDetailedShowerMode = False
SIM.enableG4GPS = False
SIM.enableG4Gun = False
SIM.enableGun = False

# .stdhep, .slcio, .HEPEvt, .hepevt, .pairs, .hepmc, .hepmc.gz, .hepmc.xz, 
# .hepmc.bz2, .hepmc3, .hepmc3.gz, .hepmc3.xz, .hepmc3.bz2, .hepmc3.tree.root
SIM.inputFiles = ['wzp6_ee_ZZ_test_ecm240_1k.stdhep']
# SIM.inputFiles = []

# Macro file to execute for runType 'run' or 'vis'
SIM.macroFile = ""

# number of events to simulate, used in batch mode
SIM.numberOfEvents = 1
SIM.skipNEvents = 0

# .slcio, edm4hep.root, .root
# SIM.outputFile = "photongun_1gev_n50_cli.root"
# SIM.outputFile = "magtest_noB_photongun240gev_n10.root"

SIM.physicsList = None

# (1) VERBOSE, DEBUG, INFO, WARNING, ERROR, FATAL, ALWAYS (7)
SIM.printLevel = 7

# The type of action to do in this invocation
# batch: just simulate some events, needs numberOfEvents, and input file or gun
# vis: enable visualisation, run the macroFile if it is set
# run: run the macroFile and exit
# shell: enable interactive session
SIM.runType = "batch"

SIM.steeringFile = None

SIM.vertexOffset = [0.0, 0.0, 0.0, 0.0]
SIM.vertexSigma = [0.0, 0.0, 0.0, 0.0]

# SIM.action.calo = "Geant4CalorimeterAction"
SIM.action.calo = "ScepcalSDAction"

SIM.action.calorimeterSDTypes = ['ScepcalSD']

# SIM.action.mapActions = { 'Scepcal' : "ScepcalSDAction" }
SIM.action.mapActions = {}

SIM.action.tracker = ('Geant4TrackerWeightedAction', {'HitPositionCombination': 2, 'CollectSingleDeposits': False})
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
                        'edep1kev': {'name': 'EnergyDepositMinimumCut', 'parameter': {'Cut': 0.001}}, 
                        'edep0': {'name': 'EnergyDepositMinimumCut/Cut0', 'parameter': {'Cut': 0.0}}
                      }

SIM.filter.calo = "edep0"

SIM.filter.mapDetFilter = {}
SIM.filter.tracker = "edep1kev"

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
## Event parameters to write in every event. Use C/F/I ids to specify parameter type. E.g parameterName/F=0.42 to set a float parameter
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
## Minimal distance between particle vertex and endpoint of parent after
##     which the vertexIsNotEndpointOfParent flag is set
SIM.part.minDistToParentVertex = 2.2e-14
## MinimalKineticEnergy to store particles created in the tracking region
SIM.part.minimalKineticEnergy = 1.0
SIM.part.printEndTracking = False
SIM.part.printStartTracking = False
SIM.part.saveProcesses = ['Decay']
# SIM.part.userParticleHandler = "Geant4TCUserParticleHandler"
SIM.part.userParticleHandler = ''


##  direction of the particle gun, 3 vector
##  shoot in the barrel region
SIM.gun.position = (0*mm, 0*mm, 0*mm)
SIM.gun.direction = (0, 1, 0)

# ## choose the distribution of the random direction for theta
# ##
# ##     Options for random distributions:
# ##
# ##     'uniform' is the default distribution, flat in theta
# ##     'cos(theta)' is flat in cos(theta)
# ##     'eta', or 'pseudorapidity' is flat in pseudorapity
# ##     'ffbar' is distributed according to 1+cos^2(theta)
# ##
# ##     Setting a distribution will set isotrop = True
# ##
# SIM.gun.distribution = 'uniform'

SIM.gun.energy = 10*GeV
# SIM.gun.particle = "gamma"

SIM.gun.multiplicity = 1

SIM.gun.momentumMin = 0.0
SIM.gun.momentumMax = 10000.0

# ##  isotropic distribution for the particle gun
# ##
# ##     use the options phiMin, phiMax, thetaMin, and thetaMax to limit the range of randomly distributed directions
# ##     if one of these options is not None the random distribution will be set to True and cannot be turned off!
# ##

SIM.gun.isotrop = False

# SIM.gun.phiMin = 0
# SIM.gun.phiMax = 2*pi

# SIM.gun.thetaMin = pi/2
# SIM.gun.thetaMax = pi/2

# ##  position of the particle gun, 3 vector. unit mm

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


SIM.random.enableEventSeed = False
SIM.random.file = None
SIM.random.luxury = 1
SIM.random.replace_gRandom = True
SIM.random.seed = None
SIM.random.type = None
