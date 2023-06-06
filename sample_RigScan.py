import pandas as pd
import os
import random
import math
import numpy as np
import time
import argparse
import sys
import configparser
import logging
import json
from rdkit import Chem
from rdkit.Chem import rdmolfiles, rdMolTransforms

import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
from decimal import Decimal
from decimal import ROUND_HALF_UP,ROUND_HALF_EVEN
import subprocess
from collections import Counter
import configparser

from util.RigScan import main as Scan

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO)

class main():
    def __init__(self, **args):
        self.param_config = args["config"]
        self._input = args["input_sdf"] ## initial sdf
        #self.prefix = _input.split(".")[0]
        self.mode = args["mode"]

        try:
            self.project_id = int(args["lbg_project_id"])
        except Exception as e:
            self.project_id = 10533


        config = configparser.ConfigParser()
        config.read(self.param_config)

        run_config = config["run_config"]

        try:
            self.platform = run_config["PLATFORM"]
        except Exception as e:
            self.platform = "local"
        else:
            if not self.platform in ["local", "lbg"]:
                self.platform = "local"

        try:
            self.nproc = int(run_config["NPROC"])
        except Exception as e:
            self.nproc = 32
        
        self._default_machine_type = {16: "c16_m128_cpu", 
                                      32: "c32_m256_cpu"}
        
        try:
            self.machine_type = self._default_machine_type[self.nproc]
        except Exception as e:
            self.machine_type = "c32_m256_cpu"
        
        
    
    def setup_local(self):
        cmd = f"python sample_RigScan.py --config {self.param_config} --input_sdf {self._input} --mode run_scan \n"
        with open("run.sh", "w+") as cc:
            cc.write("#!/bin/sh \n")
            cc.write("\n")
            cc.write(cmd)
            cc.write("\n")
        #logging.info("Prepared for local run, run by 'nohup bash run.sh &' to start running")

        return
    
    def setup_lbg(self):
        lbg_json = {
            "job_name": "rigid_dih_scan",
            "command": f"python sample_RigScan.py --config {self.param_config} --input_sdf {self._input} --mode run_scan",
            "log_file": "tmp_log",
            "backward_files": [],
            "program_id": f"{self.project_id}" ,
            "platform": "ali",
            "job_group_id": "",
            "disk_size": 128,
            "machine_type": f"{self.machine_type}",
            "job_type": "container",
            "image_name": "registry.dp.tech/dptech/prod-1364/dihscan:run0.0.3"
        }

        with open("input.json", "w+") as cc:
            json.dump(lbg_json, cc, indent=4)
        return 
    
    def setup(self):
        _dic_setup = {"local": self.setup_local, 
                      "lbg": self.setup_lbg}
        _dic_setup[self.platform]()
        return 
    
    def run_scan(self):
        Scan(config=self.param_config, 
             input_sdf=self._input).proceed()

        return 

    def option(self):
        _dic_run_mode = {"setup": self.setup, 
                         "run_scan": self.run_scan}
        
        if not self.mode in [mm for mm in _dic_run_mode.keys()]:
            logging.info(f"Wrong mode for running, should be from {[mm for mm in _dic_run_mode.keys()]}")
            return
        else:
            _dic_run_mode[self.mode]()

        return 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='rigid dih scan with g16')
    parser.add_argument('--config', type=str, required=True,
                        help='input config file')
    parser.add_argument('--input_sdf', type=str, required=True, 
                        help='input sdf file')
    parser.add_argument("--mode", type=str, required=True,
                        help="run mode, ['setup', 'run_scan']")
    parser.add_argument("--lbg_project_id", type=int, default=2273, 
                        help="if use lbg for calc define available project id, default is DPDH public 10533")
    
    args = parser.parse_args()

    main(input_sdf=args.input_sdf, 
         config=args.config, 
         mode=args.mode,
         lbg_project_id=args.lbg_project_id).option()   
