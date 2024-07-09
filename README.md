# Segmented Crystal Electromagnetic Precision Calorimeter (SCEPCAL)

![SCEPCAL3d](https://github.com/SCEPCAL/SCEPCAL/blob/main/examples/SCEPCAL3d.png?raw=true)

Repository for full simulation and analysis.

Ref: [New perspectives on segmented crystal calorimeters for future colliders](https://arxiv.org/abs/2008.00338)

See also:

[Slides from US-FCC collaboration meeting (March 2024)](https://indico.mit.edu/event/876/contributions/2860/)

[Slides from CALOR (May 2024)](https://indico.cern.ch/event/1339557/contributions/5898506/)


## Instructions
### Compile and Install

**Note:** use lxplus9

```sh
git clone git@github.com:SCEPCAL/SCEPCAL.git
cd SCEPCAL

export TOP_DIR=$PWD

mkdir build
mkdir install

export MY_INSTALL_DIR=$PWD/install

source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh

cd build
cmake -DCMAKE_INSTALL_PREFIX=$MY_INSTALL_DIR -DPython_EXECUTABLE=$(which python) ..
make -j4
make install

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MY_INSTALL_DIR/lib64
export PYTHONPATH=$PYTHONPATH:$MY_INSTALL_DIR/python

```

### Event selection
Use the built-in `ddsim` particle gun (described further below) or generate events from your favorite MC program. `ddsim` currently accepts the following inputs:

```
.stdhep
.slcio
.HEPEvt
.hepevt
.pairs
.hepmc [.gz, .xz, .bz2]
.hepmc3 [.gz, .xz, .bz2, .tree.root]
```

Set your input and output files in `scripts/scepcal_steering.py`:

```python
SIM.inputFiles = ['examples/wzp6_ee_ZZ_test_ecm240_1k.stdhep']
SIM.outputFile = 'examples/wzp6_ee_ZZ_test_ecm240_n1_output.root'
```

See also: https://github.com/HEP-FCC/EventProducer/

### Running the simulation

Set all options in the steering file `scripts/scepcal_steering.py` and/or set them via the command line at runtime. Command line options override the steering file.

See `scripts/dd4hep_steering_template.py` for explanations of all options.

Once the steering file is configured as desired,

```sh
cd $TOP_DIR
ddsim scripts/scepcal_steering.py
```

### Simulation Options

#### Optical Physics

Disabled by default. Uncomment line 195 in `scripts/scepcal_steering.py` to enable optical physics if you are confident in your computing resources.

```python
# SIM.physics.setupUserPhysics(setupCerenkovScint)
```

#### Energy Cuts

You may change the energy cut if you wish. Energy deposits in the readout are accumulated per calorimeter cell; the energy cut here refers to the energy deposited per Geant4 particle step. Using an energy cut of 0 will work to show every single deposit but is not likely to be useful. See lines 56-62 of `scripts/scepcal_steering.py`. The default cut is set to 1 keV.

```python
... 'edep1keV': {'name': 'EnergyDepositMinimumCut/1keV', 'parameter': {'Cut': 1.0*keV }}
...
SIM.filter.calo = "edep1keV"
```

#### Random Seed

Change the random seed on line 204 of `scripts/scepcal_steering.py` for reproducibility. Default is `None`.

```python
# SIM.random.seed = None
SIM.random.seed = 2126136508
```

#### Example using command line options to run a particle gun input

```sh
ddsim --steeringFile scripts/scepcal_steering.py -G --gun.direction "1 1 0" --gun.energy "1*GeV" --gun.particle="gamma" -O gamma_1GeV_noopt_1MeVcut.root
```

`-G` and the following options tells the simulation to use the particle gun with the stated settings. Can additionally specify `--gun.position "0 0 0"`, etc. See the steering template file for more. **Note that a photon is called a gamma for the particle gun.**

`-O` sets the output filename.

`-N` sets the number of events (omitted here, set to 1 in steering file)

#### Changing the Geometry

Detector and crystal dimensions can be changed in `install/share/compact/SCEPCAL.xml`. This file is copied from `Detector/SCEPCAL/compact/SCEPCAL.xml` when you run `make install`, so any changes will be lost every time you run `make install`. The timing layer and barrel/endcap are written as separate subdetectors. **The timing layer is disabled by default.**

```xml
<detectors>
    <detector id="2"
                name="Scepcal_BE"
                type="SCEPCAL_BE" 
                readout="SCEPCAL_BE_readout"
                vis="scepcalAssemblyGlobalVis"
                sensitive="true">
        ...
    </detector>
    
<-- Disabled by default
    <detector id="3"
                name="Scepcal_TL" 
                type="SCEPCAL_TL" 
                readout="SCEPCAL_TL_readout"
                vis="scepcalAssemblyGlobalVis"
                sensitive="true">
        ...
    </detector>
-->
</detectors>
```

Uncomment the `<detector>` tag for the timing layer to enable it.

See the slides linked at the top for a brief description of the geometry parametrization and construction, or see it directly in `Detector/SCEPCAL/src/`.

#### Material Properties

Simulation results are heavily dependent on the material properties defined for the detector, here we use LYSO crystals for the timing layer and PbWO4 for the barrel/endcap. The material properties are defined in the compact XML file described above, however they are pulled from various sources and are by no means set in stone. Change the material properties as you like, or contact to collaborate about them.

### Running on Condor

Change the user directory paths in `scripts/SCEPCALsim.sh` and `scripts/SCEPCALsim.sub` to your own, and also change the path of the compact XML file in `scripts/scepcal_steering.py` to the full absolute file path rather than the relative path, e.g.:

```python
# SIM.compactFile = ['install/share/compact/SCEPCAL.xml']
SIM.compactFile = ['/eos/user/w/wochung/src/SCEPCAL/install/share/compact/SCEPCAL.xml']
```

and then run

```sh
condor_submit scripts/SCEPCALsim.sub
```

on lxplus9.

### Analysis

See `examples/SCEPCALanalysis.ipynb` for a simple example plotting calorimeter hits of type `edm4hep.SimCalorimeterHitData`:

![BEhits](https://github.com/SCEPCAL/SCEPCAL/blob/main/examples/BEhits.png?raw=true)

The plot above is from `examples/wzp6_ee_ZZ_test_ecm240_n1_cut0_BEonly.root`, which shows 1 event with an energy cut of 0, Barrel/Endcap only (i.e. no timing layer).

The files `examples/gamma_10MeV_d<xyz>_BE_noopt_1MeVcut.root` correspond to runs with the particle gun, shooting gammas at 10 MeV in the xyz vector direction, from the origin, Barrel/Endcap only (i.e. no timing layer), with no optical physics, with a 1 MeV energy cut applied. Reproduce these runs and plot them to verify the particle gun is shooting as expected, and that the detector is oriented as expected.

### Collaborate and Contact

This release is stable, but there is a lot of duplicated code in the detector constructors that need refactoring.

wonyongc@princeton.edu