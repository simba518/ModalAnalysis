#! /usr/bin/env python

import os

web_data_dir = "//10.76.0.107/Experiment/data"
web_app = "//10.76.0.107/modal_analysis"
local_data_dir = "./local_data"
local_app = "./modal_analysis"

if not os.path.exists(local_data_dir):
    os.system("mkdir "+local_data_dir)

os.system("mount -t cifs "+web_data_dir+" "+local_data_dir)

def modal_analysis(mass_filepath, num_modes=16):
    model_filepath = mass_filepath[0:-4]
    num_nodes_per_elem = 8
    if (model_filepath.find("_hex") < 0):
        num_nodes_per_elem = 4
    cmd = local_app+" "+model_filepath+"stiffness "+model_filepath+"mass "+str(num_modes)+" "+model_filepath+"off "+str(num_nodes_per_elem)
    print "running command: "+cmd
    os.system(cmd)

for dirname, dirnames, filenames in os.walk(local_data_dir):
    for filename in filenames:
        filepath = str(os.path.join(dirname, filename))
        if filepath.find(".svn") < 0 and filepath.find(".mass") > 0:
            modal_analysis(filepath)

os.system("umount "+local_data_dir)
