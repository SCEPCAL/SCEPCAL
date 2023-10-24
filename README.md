# SCEPCAL
Repository for Gaudi simulation &amp; analysis of the DR Segmented Crystal ECAL. Forked from https://github.com/HEP-FCC/dual-readout.

## How-to
### Compile
After fetching the repository, do

```sh
source init_lcg.sh # init_hsf.sh for /cvmfs/sw.hsf.org/key4hep
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=<path_to_install_directory> ..
make -j4
make install
```

Note that to use the installed binary & library files, need to do following (assuming `$PWD=<path_to_install_directory>`)

```sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/lib64
export PYTHONPATH=$PYTHONPATH:$PWD/python
```

Everytime you log in:
```sh
source setupSCEPCAL.sh
```
To get the correct env.


### Generating events
Generating events relies on [k4Gen](https://github.com/HEP-FCC/k4Gen), generates primary particle(s) in `HepMC` format with either particle gun or `Pythia8` then converts it to `edm4hep`. Please refer to [k4Gen/options](https://github.com/HEP-FCC/k4Gen/tree/b7c735e401298a8c72915819dc0404a83f46a0fe/k4Gen/options) for example configurations.

### Running GEANT4 simulation
An example configuration is located at `SCEPCALsim/SCEPCALsimG4Components/test/runSCEPCALsim.py`. After modifying the configuration based on your needs, run

```sh
k4run runSCEPCALsim.py
```

Note `Gaudi` components specific to this simulation - `SimG4OpticalPhysicsList`, and `SimG4SCEPCALActions`.

```python3
from Configurables import SimG4Svc, SimG4OpticalPhysicsList
opticalPhysicsTool = SimG4OpticalPhysicsList("opticalPhysics", fullphysics="SimG4FtfpBert")

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

from Configurables import SimG4SaveSmearedParticles, SimG4SaveDRcaloHits, SimG4SaveDRcaloMCTruth
savePtcTool = SimG4SaveSmearedParticles("saveSmearedParticles")
saveSCEPCALTool = SimG4SaveSCEPCALHits("saveSCEPCALTool", readoutNames = ["SCEPCALreadout"])
saveMCTruthTool = SimG4SaveSCEPCALMCTruth("saveMCTruthTool") # need SimG4SCEPCALActions

geantsim = SimG4Alg("SimG4Alg",
  outputs = [
    "SimG4SaveSmearedParticles/saveSmearedParticles",
    "SimG4SaveSCEPCALHits/saveSCEPCALTool",
    "SimG4SaveSCEPCALMCTruth/saveMCTruthTool"
  ],
  eventProvider = edmConverter
)
```

Optical physics is NOT simulated by the `GEANT4` default physics list due to extensive compute, it is turned on separately in `SimG4OpticalPhysicsList` and the Cherenkov and scillation processes are configured, however in `SimG4SCEPCALSteppingAction` Cherenkov photons are killed on the first step after being counted.

`SimG4SCEPCALoActions` is responsible for initializing `SimG4SCEPCALSteppingAction`, which retrieves MC truth energy deposits inside each cell. The resulting MC-truth energy deposit and counted number of photoelectrons are stored separate `edm4hep` collections of type `SimCalorimeterHits`.

### Running on Condor

Go to `condor_submission` folder. 
Here you find a template of the gaudi config to run the SCEPCal Simulation `runSCEPCalsimTemplate.py` and the `submitJobs.py` scripts.
The `submitJobs.py` creates the gaudi configs to be run for each simulation job from the template and send the jobs on condor.
It also set up the right env to run SCEPCal simulation on condor.

Here an examples of usage:


```sh
python submitJobs.py \
--label MYLABEL \ 
--basedir MYSCEPCALDIR \
--cfg runSCEPCALsimTemplate.py \ # you can provide your own template
--jobs #JOBS \
--events #EVENTSperJOB \ 
--energy ENERGYOFGUN \
--pdg PDG \
--queue CONDORQUEUE \
--submit #to submit jobs, if not specified only dry run
```

In the `condor_submission` folder, it will create a folder named MYLABEL.
In this folder, three subdirectories are created:
- `config` contains the gaudi config created for each job from `runSCEPCalsimTemplate.py`;
- `condor` contains the .out, .err, .log of the condor jobs.
- `output` contains the .root produced.


