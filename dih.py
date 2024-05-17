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
#

logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO)

class main():
    def __init__(self, **args):
        self.param_config = args["config"]
        #self.run_type = args["run_type"]
        #self.prefix = _input.split(".")[0]
        #self.mode = args["mode"]

        try:
            self.do_submit = args["submit"]
        except Exception as e:
            self.do_submit = False
        
        try:
            self.do_col_result = args["col_result"]
        except Exception as e:
            self.do_col_result = False
        
        try:
            self.do_scan = args["scan"]
        except Exception as e:
            self.do_scan = False


        config = configparser.ConfigParser()
        config.read(self.param_config)

        general = config["GENERAL"]

        try:
            self._input = general["input_sdf"]
        except Exception as e:
            self._input = None
        else:
            if not os.path.isfile(self._input):
                self._input = None

        #try:
        #    self.platform = general["platform"]
        #except Exception as e:
        #    self.platform = "lbg"
        #else:
        #    if not self.platform in ["local", "lbg"]:
        #        self.platform = "lbg"
        
        try:
            self.nproc = int(general["nproc"])
        except Exception as e:
            self.nproc = 32
        
        self._default_machine_type = {16: "c16_m128_cpu", 
                                      32: "c32_m256_cpu"}
        
        try:
            self.machine_type = self._default_machine_type[self.nproc]
        except Exception as e:
            self.machine_type = "c32_m256_cpu"
        

        lbg = config["LBG"]
        try:
            self.project_id = lbg["lbg_project_id"]
        except Exception as e:
            self.project_id = None

        try:
            self.lbg_image = lbg["lbg_image"]
        except Exception as e:
            self.lbg_image = None
        
        self.main_dir = os.getcwd()
        
    """
    def setup_local(self):
        cmd = f"python dih.py --config {self.param_config} --input_sdf {self._input} --mode run_scan \n"
        with open("run.sh", "w+") as cc:
            cc.write("#!/bin/sh \n")
            cc.write("\n")
            cc.write(cmd)
            cc.write("\n")
        #logging.info("Prepared for local run, run by 'nohup bash run.sh &' to start running")

        return "run.sh"
    """
    
    def setup_lbg(self):
        if not (self.project_id and self.lbg_image):
            logging.info("Missing requiered information for LBG setting, setup failed")
            return None
        
        logging.info("Preparing ...")
        
        lbg_json = {
            "job_name": "DIH_scan",
            "command": f"python dih.py --config {self.param_config} --scan",
            "log_file": "tmp_log",
            "backward_files": [],
            "program_id": f"{self.project_id}" ,
            "platform": "ali",
            "job_group_id": "",
            "disk_size": 128,
            "machine_type": f"{self.machine_type}",
            "job_type": "container",
            "image_name": f"{self.lbg_image}"
        }

        with open("input.json", "w+") as cc:
            json.dump(lbg_json, cc, indent=4)
        return "input.json"
    
    def submit(self):
        if not self._input:
            logging.info("Input mol file doesn't exist, nothing to do")
            return
        
        setup_tag = self.setup_lbg()

        if not setup_tag:
            return 
        
        logging.info("Ready for job submission...")

        tracker = {}
        
        submission_cmd = f"lbg job submit -i {setup_tag} -p ./ > TRACK_JOB"

        (_, _) = subprocess.getstatusoutput(submission_cmd)   
            
        if os.path.isfile("TRACK_JOB") and os.path.getsize("TRACK_JOB"):
            with open("TRACK_JOB", "r+") as f:
                content = [ll.strip() for ll in f.readlines() if ll]
            
            if content:
                try:
                    JOB_ID = content[-1].split()[-1]
                except Exception as e:
                    JOB_ID = None
            else:
                JOB_ID = None
            
            if JOB_ID:
                tracker.setdefault(os.getcwd(), JOB_ID)
            else:
                logging.info(f"Submit failed, check in folder for more information")
        else:
            logging.info(f"Submit failed, check in folder for more information")
        
        if tracker:
            df_track = pd.DataFrame({"local_path": [kk for kk in tracker.keys()],
                                 "JOB_ID": [vv for vv in tracker.values()]})
        
            df_track.to_csv(os.path.join(self.main_dir,"submission_track_list"), index=None)
        else:
            logging.info("Failed")
            return 

        logging.info("Finish")

        return 

    def col_result(self):
        
        if not os.path.isfile("submission_track_list"):
            logging.info("No tracking information available, abort")
            return
        
        logging.info("Collecting results ....")

        df = pd.read_csv("submission_track_list", header=0)

        for idx, row in df.iterrows():
            os.chdir(row["local_path"])

            cmd = f"lbg job download {row['JOB_ID']}"
            (_, _) = subprocess.getstatusoutput(cmd)
            if os.path.isdir(str(row['JOB_ID'])):
                (_, _) = subprocess.getstatusoutput(f"mv {row['JOB_ID']}/* ./")
                (_, _) = subprocess.getstatusoutput(f"rm -rf {row['JOB_ID']}")
                #transfer += [os.path.join(row["local_path"], dd) for dd in os.listdir(row["local_path"]) if os.path.isdir(dd)]
                #should_delete.append(row["local_path"])
                os.system("rm -f submission_track_list TRACK_JOB fort.* lbg*.sh")
                logging.info("Done") 
            else:
                logging.info("lbg job collection failed, check TRACK_JOB for more detail")
                return 
            
            os.chdir(self.main_dir)
        
        return 
                    
    
    def scan(self):
    
        Scan(config=self.param_config, 
             input_sdf=self._input).run()

        return 

    def option(self):
        if self.do_scan:
            self.scan()
        
        elif self.do_submit:
            self.submit()
        
        elif self.do_col_result:
            self.col_result()
        else:
            logging.info("No actually command, terminate")

        return 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='rigid dih scan with g16')
    parser.add_argument('--config', type=str, required=True,
                        help='input config file')
    parser.add_argument('--submit', default=False, action='store_true', 
                        help='if add, do lbg job preparison and submit')
    parser.add_argument('--col_result', default=False, action='store_true', 
                        help='if add, do lbg job collection')
    parser.add_argument('--scan', default=False, action='store_true',
                        help="if add, will initiate dih scan")
    
    args = parser.parse_args()

    main(config=args.config,
         submit=args.submit,
         col_result=args.col_result,
         scan=args.scan).option()   
