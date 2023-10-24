#!/bin/python
import os
import glob
import math
from array import array
import sys
import time
import subprocess
import argparse
import time
import datetime

#print date
print("----------------------------------------------------------------------------------")
print(datetime.datetime.now())
print("----------------------------------------------------------------------------------")

#parse arguments
parser = argparse.ArgumentParser(description='Module characterization summary plots')
parser.add_argument("--label",       required=True,   type=str,   help="job label")
parser.add_argument("--cfg",         required=True,   type=str,   help="gaudi config template file")
parser.add_argument("--basedir",    required=True,   type=str,   help="base folder of your scepcal repo")
parser.add_argument("--jobs",      required=True,   type=int,   help="number of jobs")
parser.add_argument("--events",     required=True,   type=int,   help="number of events for each job")
parser.add_argument("--energy",     required=True,   type=str,   help="energy of the particle")
parser.add_argument("--pdg",        required=True,   type=int,   help="PDG number for particle gun")
parser.add_argument('--submit',              action='store_true',         default=False,      help='submit jobs')
parser.add_argument("--queue",                          type=str, default="workday",  help="condor queue: espresso, longlunch, workday...")  

args = parser.parse_args()


myJobsDir = args.basedir+"condor_submission/"+args.label
cfgDir    = myJobsDir+"/config/"
condorDir = myJobsDir+"/condor/"
outDir    = myJobsDir+"/output/"
cfglist = ""
# create folders to store the gaudi config file, .out,.err,.log from condor,  and  .root files
os.system("mkdir -p "+myJobsDir+" "+cfgDir+" "+condorDir+" "+outDir)



myConfString = "PDG%i_E%sGeV"%(args.pdg, args.energy)
for ijob in range(0,args.jobs):
   
  with open(args.cfg) as fi:
    contents = fi.read()
    replaced_contents = contents.replace("MYPARTICLE", str(args.pdg))
    replaced_contents = replaced_contents.replace("MYBASEDIR",str(args.basedir)) 
    replaced_contents = replaced_contents.replace("MYENERGY",str(args.energy)) 
    replaced_contents = replaced_contents.replace("MYNEVENTS",str(args.events)) 
    replaced_contents = replaced_contents.replace("MYJOBID",str(ijob)) 
    cfgfilename=cfgDir+"/runSCEPCALsim_%s_job%i.py"%(myConfString, ijob)
    cfglist += cfgfilename+"\n"
  with open(cfgfilename, "w") as fo:
    fo.write(replaced_contents)

with open (myJobsDir+"/condor_arg_many_jobs.txt", "w") as argfile:
  argfile.write(cfglist)
 
##### creates script #######
outScriptName=myJobsDir+"/SCEPCALsim.sh"
outScript = open(outScriptName,"w")
outScript.write("#!/bin/sh -e\n")
outScript.write('CONFIGFILE=$1\n')
#outScript.write('source %s/setupSCEPCAL.sh\n'%args.basedir)

outScript.write("source /cvmfs/sft.cern.ch/lcg/views/LCG_100/x86_64-centos7-gcc8-opt/setup.sh\n")
outScript.write("export PYTHIA8_ROOT_DIR=/cvmfs/sft.cern.ch/lcg/views/LCG_100/x86_64-centos7-gcc8-opt/\n")
outScript.write("source /cvmfs/sft.cern.ch/lcg/releases/LCG_100/DD4hep/01.16.01/x86_64-centos7-gcc8-opt/bin/thisdd4hep.sh\n")
outScript.write("source %s/init_k4.sh\n"%args.basedir)
outScript.write("source %s/init_sipm.sh\n"%args.basedir)

outScript.write("export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:%s/install/lib64\n"%args.basedir)
outScript.write("export PYTHONPATH=$PYTHONPATH:%s/install/python\n"%args.basedir)


outScript.write('cd %s\n'%outDir)
outScript.write('k4run $CONFIGFILE\n')
outScript.write('cd -\n')
outScript.close()

#generate condor multijob submitfile for each task
condorsubFilename=myJobsDir+"/submit_jobs.sub"
condorsub = open( condorsubFilename,"w")
condorsub.write("executable            = %s\n"%outScriptName)
condorsub.write("output                = %s$(ClusterId).$(ProcId).out\n"%condorDir)
condorsub.write("error                 = %s$(ClusterId).$(ProcId).err\n"%condorDir)
condorsub.write("log                   = %s$(ClusterId).log\n"%condorDir)
condorsub.write("+JobFlavour           = %s\n"%args.queue)
condorsub.write("queue arguments from %s/condor_arg_many_jobs.txt"%myJobsDir)
condorsub.close()

submit_command = "condor_submit "+condorsubFilename
print("SUBMIT COMMAND: "+submit_command)
#submit in case the option is given
if(args.submit):
    os.system(submit_command)
