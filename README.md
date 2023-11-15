# SCEPCAL
Repository for Gaudi simulation &amp; analysis of the DR Segmented Crystal ECAL. Forked from https://github.com/HEP-FCC/dual-readout.

## How-to
### Compile
After fetching the repository, do

```sh
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
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

Adjust the files `SCEPCALsim.sh` and `SCEPCALsim.sub` to your configuration and then

```
condor_submit SCEPCALsim.sub
```

on lxplus. Note that the shell file sources a `SCEPCALenv.sh` file created by you to get the env variables to be able to use your local installation.

```
cd install/test
env | sed -e 's/^/export\ /' > SCEPCALenv.sh
chmod +x SCEPCALenv.sh
```

This will write your current environment to file, appending the string "export " at the beginning of each line to make it sourceable. This should be run in the same login session as when you ran `init_lcg.sh`. Since all environment variables, including ones irrelevant to the simulation, are printed by `env`, it is recommended to read through the generated file and remove env variables that are unneeded.
