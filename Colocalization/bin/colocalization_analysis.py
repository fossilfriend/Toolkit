#!/usr/bin/env python

'''
uses python R wrapper to do colocalization analysis comparing two gwas datasets
takes a json config file

REQUIRES: python3.x

[
"regions": "full path to region file",
"conditions": [
{"name": "name of condition 1",
 "file": "full path to condition 1 gwas",
{<condition 2>}
]
]
'''


import argparse
import json

import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import rpy2.robjects.lib.ggplot2 as ggplot2
# import pandas.rpy.common as rcom

from os import path
from Toolkit.Utils.utils import warning, xstr, qw


def load_data(fileConfig, sep="\t"):
    ''' load data from the specified file path/update column names '''
    return pd.read_csv(fileConfig.file, sep=sep)
  

def parse_config():
    ''' parse JSON config file '''
    try:
        with open(args.config, 'r') as f:
            data = json.load(f)

    except (OSError, IOError) as e: # FileNotFoundError does not exist below python 3.3
        warning("Could not open config file " + args.config)
        die(e)

    except JSONDecodeError as e:
        warning("Error parsing config file")
        die(e)

    return data


def perform_colocalization():
    ''' compare data1 & data2 '''
    warning("Performing colocalization analysis")
    result = {}
    for row in regions:
        die(row.chr)


if __name__ == "__main__":
    """ perform colocationzation analysis """
    parser = argparse.ArgumentParser("Perform colocalization analysis comparing two gwas datasets")
    parser.add_argument('-c', '--config', help="full path to json config file", required=True)
    parser.add_argument('-o', '--outputFilePath', help="full output file path", required=True)
    parser.add_argument('--verbose', action='store_true')

    args = parser.parse_args()

    config = parse_config()

    condition1 = config.conditions[0]
    condition2 = config.conditions[1]

    warning("Loading Regions from", config.regions)
    regions = load_data(config.regions, sep=",")
    
    warning("Loading GWAS for", condition1.name)
    #data1 = load_data(condition1)

    warning("Loading GWAS for", condition2.name)
    #data2 = load_data(condition2)

    #base = importr('base')
    #coloc = importr('coloc')

    perform_colocalization()
   
