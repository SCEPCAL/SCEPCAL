from Gaudi.Configuration import *
from Configurables import ApplicationMgr
from GaudiKernel import SystemOfUnits as units

from Configurables import k4DataSvc
dataservice = k4DataSvc("EventDataSvc")

from Configurables import GenAlg, MomentumRangeParticleGun
pgun = MomentumRangeParticleGun("PGun",
  PdgCodes=[11], # electron
  MomentumMin = 3.*units.GeV, # GeV
  MomentumMax = 3.*units.GeV, # GeV
  ThetaMin = 0.1, # rad
  ThetaMax = 3.1, # rad
  PhiMin = 0.0, # rad
  PhiMax = 6.28, # rad
)

from Configurables import FlatSmearVertex
smearTool = FlatSmearVertex("VertexSmearingTool",
  yVertexMin = -36.42, # mm
  yVertexMax = -26.42, # mm
  zVertexMin = -52.135, # mm
  zVertexMax = -42.135, # mm
  beamDirection = 0 # 1, 0, -1
)

from Configurables import HepMCToEDMConverter
hepmc2edm = HepMCToEDMConverter("Converter")

gen = GenAlg("ParticleGun", SignalProvider=pgun, VertexSmearingTool=smearTool)

from Configurables import GeoSvc
geoservice = GeoSvc(
  "GeoSvc",
  detectors = [
    'file:/afs/cern.ch/work/f/fcetorel/private/work2/scepcal_key4hep/install/share/compact/SCEPCAL.xml'
  ]
)

from Configurables import SimG4Svc, SimG4OpticalPhysicsList

opticalPhysicsTool = SimG4OpticalPhysicsList("opticalPhysics", fullphysics="SimG4FtfpBert", SetScintillation=False)

from Configurables import SimG4SCEPCALActions
actionTool = SimG4SCEPCALActions("SimG4SCEPCALActions")

# Name of the tool in GAUDI is "XX/YY" where XX is the tool class name and YY is the given name
geantservice = SimG4Svc("SimG4Svc",
  physicslist = opticalPhysicsTool,
  actions = actionTool
)

from Configurables import SimG4Alg, SimG4PrimariesFromEdmTool
# next, create the G4 algorithm, giving the list of names of tools ("XX/YY")
edmConverter = SimG4PrimariesFromEdmTool("EdmConverter")

from Configurables import SimG4SaveSCEPCALHits, SimG4SaveSCEPCALMCTruth
saveSCEPCALTool = SimG4SaveSCEPCALHits("saveSCEPCALTool", readoutNames = ["SCEPCALreadout"])
saveMCTruthTool = SimG4SaveSCEPCALMCTruth("saveMCTruthTool") # need SimG4SCEPCALActions

geantsim = SimG4Alg("SimG4Alg",
  outputs = [
    "SimG4SaveSCEPCALHits/saveSCEPCALTool",
    "SimG4SaveSCEPCALMCTruth/saveMCTruthTool"
  ],
  eventProvider = edmConverter
)

from Configurables import PodioOutput
podiooutput = PodioOutput("PodioOutput", filename = "scepcal_lxplus.root")
podiooutput.outputCommands = ["keep *"]

from Configurables import RndmGenSvc, HepRndm__Engine_CLHEP__RanluxEngine_
rndmEngine = HepRndm__Engine_CLHEP__RanluxEngine_("RndmGenSvc.Engine",
  SetSingleton = True,
  UseTable = True,
  Column = 0, # 0 or 1
  Row = 123 # 0 to 214
)

rndmGenSvc = RndmGenSvc("RndmGenSvc",
  Engine = rndmEngine.name()
)

ApplicationMgr(
  TopAlg = [gen, hepmc2edm, geantsim, podiooutput],
  EvtSel = 'NONE',
  EvtMax = 100,
  # order is important, as GeoSvc is needed by SimG4Svc
  ExtSvc = [rndmEngine, rndmGenSvc, dataservice, geoservice, geantservice]
)
