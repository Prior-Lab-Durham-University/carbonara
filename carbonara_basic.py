#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 18:46:24 2026

@author: rob
"""

import CarbonaraDataTools as CDT
import numpy as np
import re
import fittingAnalysis as fa
import subprocess


# @title Run this to set the carbonara model up with a given run name
run_name = "humanSmarcal"
pdb_name = "pdbFiles/foldSmarcal.cif"
saxs_name = "saxsFiles/smrclcnc_a2.dat"
pae_name = "paeFiles/foldSmarcal.json"


subprocess.run([ "python", "setup_carbonara_allAtom.py", "-p", pdb_name, "-s", saxs_name, "-n", run_name ])


# for use with pae file (f for flexibility)

#!python setup_carbonara_allAtom.py -p $pdb_name -s $saxs_name -f $pae_name -n $run_name --alphaFoldFlex

#set up foxs script
from pathlib import Path
foxs_cmd = f'python3 {Path("external/pyFoXS/pyFoXS/foxs.py").resolve()}'

max_q=0.2

foxs_result = CDT.run_initial_foxs_check(
    pdb_name=pdb_name,
    saxs_name=saxs_name,
    foxs_cmd=foxs_cmd,
     max_q =max_q
)

CDT.toggle_startk("RunMe_"+run_name+".sh",max_q)

from pathlib import Path
from run_frontend import CarbonaraRunner

runner = CarbonaraRunner(run_name, foxs_cmd=foxs_cmd)
runner.start()
runner.start_monitor(threshold=2.5, every_s=10, defer_backmap_seconds=600)

