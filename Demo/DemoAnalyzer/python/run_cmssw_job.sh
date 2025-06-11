#!/bin/bash

MODE=$1  # QCD / EWK / Interference
echo "DEBUG: MODE is '${MODE}'"

# Setup CMSSW environment without changing working directory
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh
cd /afs/cern.ch/user/p/pravesh/private/MYDEMOANALYZER/CMSSW_15_0_3/src
eval `scramv1 runtime -sh`

# Handle Condor scratch directory fallback for manual execution
if [ -z "$CONDOR_SCRATCH_DIR" ]; then
    echo "CONDOR_SCRATCH_DIR not set, using current directory."
    export CONDOR_SCRATCH_DIR=$(pwd)
fi

cd "$CONDOR_SCRATCH_DIR"  # Back to Condor environment
cd Demo/DemoAnalyzer/
# If the directory exists, remove it
#if [ -d "output" ]; then
#    echo "Directory 'output' exists. Removing it..."
#    rm -rf output
#fi

# Create a new empty output directory
#mkdir -p output
#echo "Created new 'output' directory."

# Define input/output for each mode
if [ "$MODE" == "QCD" ]; then
  INPUT="step1_0.root"
  OUTPUT="output/ZGJJtoNuNuGJJ_QCD.root"
elif [ "$MODE" == "EWK" ]; then
  INPUT="step1_0.root"
  OUTPUT="output/ZGJJtoNuNuGJJ_EWK.root"
elif [ "$MODE" == "Interference" ]; then
  INPUT="step1_0.root"
  OUTPUT="output/ZGJJtoNuNuGJJ_Interference.root"
else
  echo "Unknown mode: $MODE"
  exit 1
fi

# Diagnostic output
echo "Listing files in Condor scratch:"
ls -lh
echo "Working directory is: $(pwd)"

# Path to original config
CONFIG_SRC="/afs/cern.ch/user/p/pravesh/private/MYDEMOANALYZER/CMSSW_15_0_3/src/Demo/DemoAnalyzer/python/ConfFile_cfg.py"

# Ensure config file exists
if [ ! -f "$CONFIG_SRC" ]; then
  echo "Error: ConfFile_cfg.py not found at $CONFIG_SRC!"
  exit 1
fi

# Copy and modify config file
cp "$CONFIG_SRC" python/ConfFile_mode_cfg.py
FILELIST="/eos/user/p/pravesh/Condor/ZGJJtoNuNuGJJ_${MODE}/step1_0.root', \
'file:/eos/user/p/pravesh/Condor/ZGJJtoNuNuGJJ_${MODE}/step1_1.root', \
'file:/eos/user/p/pravesh/Condor/ZGJJtoNuNuGJJ_${MODE}/step1_2.root', \
'file:/eos/user/p/pravesh/Condor/ZGJJtoNuNuGJJ_${MODE}/step1_3.root', \
'file:/eos/user/p/pravesh/Condor/ZGJJtoNuNuGJJ_${MODE}/step1_4.root', \
'file:/eos/user/p/pravesh/Condor/ZGJJtoNuNuGJJ_${MODE}/step1_5.root', \
'file:/eos/user/p/pravesh/Condor/ZGJJtoNuNuGJJ_${MODE}/step1_6.root', \
'file:/eos/user/p/pravesh/Condor/ZGJJtoNuNuGJJ_${MODE}/step1_7.root', \
'file:/eos/user/p/pravesh/Condor/ZGJJtoNuNuGJJ_${MODE}/step1_8.root', \
'file:/eos/user/p/pravesh/Condor/ZGJJtoNuNuGJJ_${MODE}/step1_9.root"


sed -i "s|__INPUTFILE__|${FILELIST}|" python/ConfFile_mode_cfg.py
sed -i "s|__OUTPUTFILE__|${OUTPUT}|" python/ConfFile_mode_cfg.py

# Run CMSSW job
cmsRun python/ConfFile_mode_cfg.py
