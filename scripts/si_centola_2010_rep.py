# Centola 2010 "in silico" 복제 시뮬레이션 (SI 섹션 L 데이터 생성, 2D 격자 네트워크)
"""
Script for "in silico" replication of Centola 2010
Reproduces data for SI Section L

Parameter settings
k = 6,8
i = 2
b = 1,2
sig_thresh_sd = 0, 0.1
2D lattices
One random seed
"""

import networkx as nx
import random
import numpy as np
import pandas as pd
import os
import math
import argparse
import helper_functions as hf

# command line arguments
parser = argparse.ArgumentParser(description="parameters to run network contagion model")

parser.add_argument('-r_start', '--r_start', type=int, required=True,
                    help='seed to start trials at')
                    
parser.add_argument('-tr', '--trials', type=int, required=True,
                    help='number of trials for each parameter run')
                    
parser.add_argument('-k', '--k', type=int, required=True,
                    help='network degree')
                    
parser.add_argument('-b', '--b', type=int, required=True, # mem
                    help='time of influence')
                    
parser.add_argument('-p1', '--p1', type=float, required=True,
                    help='below threshold adoption rate ')
                    
parser.add_argument('-sig_thresh_sd', '--sig_thresh_sd', type=float, required=True,
                    help='heterogeneity')

args = parser.parse_args()


# Define k and network types
# Hexagonal lattice if k = 6, Moore lattice if k = 8
k = args.k

if k == 6:
    G = nx.triangular_lattice_graph(16,16, periodic = True)
    mappingG = {n:i for i, n in enumerate(G.nodes())}
    G = nx.relabel_nodes(G, mappingG)
    G_name = "Hex"
    # 128 nodes
if k == 8:
    G = hf.moore_lattice(12,periodic= True, create_using=None)
    mappingG = {n:i for i, n in enumerate(G.nodes())}
    G = nx.relabel_nodes(G, mappingG)
    G_name = "Moore"
    # 144 nodes
    
universal_seed_strat = hf.seed_strat_random

# Run simulation
#for rand_net in range(10):
#G_name = str(rand_net)
for perc_rewire in [0,1]:
    hf.main(G, 
            G_name, 
            k = args.k, 
            model = hf.sirsic_wmem , 
            trials = args.trials, 
            r_start = args.r_start, 
            beta = args.b, 
            perc = perc_rewire, 
            rand_path = None, 
            seeds = 1, 
            seed_strat0 = universal_seed_strat,
            seed_strat = universal_seed_strat,
            sig_thresh_sd = args.sig_thresh_sd,
            thrshld = 2, 
            p1 = args.p1, 
            #p2_list = [0.6], 
            p2_list = np.arange(args.p1, 1.02, step = 0.02),
            cf_file = "/home/wan.a/complex_contagion_repo/output_data/centola_2010/cf_sim_" 
                        + str(args.k)+"k_"+str(args.b)+"T_"
                        +"2td_"+str(args.sig_thresh_sd)+"sd"+str(args.p1)+"p1"+".csv", 
            df_file = None,
            full_series = True,
            rand_network = "ES")


