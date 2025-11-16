# 메인 복잡한 전파 시뮬레이션 실행 스크립트 (모든 주요 그림 및 SI 섹션 B-I 데이터 생성)
"""
Script to run main model

Reproduces data for:
- All main figures
- SI sections B - I
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

parser.add_argument('-tr', '--trials', type=int, required=True,
                    help='number of trials for each parameter run')

parser.add_argument('-G_name', '--G_name', type=str, required=True,
                    help='"WS" (Watts strogatz ring lattice), "MR" (Moore),\
                    "HX" (Hex)')
                    
parser.add_argument('-n', '--n', type=int, required=True,
                    help='number of nodes')
                    
parser.add_argument('-k', '--k', type=int, required=True,
                    help='network degree')
                    
parser.add_argument('-b', '--b', type=int, required=True, # mem
                    help='time of influence')
                    
parser.add_argument('-td', '--thrshld', type=int, required=True, #i
                    help='number of neighbors to activate p2 value')
                    
parser.add_argument('-p1', '--p1', type=float, required=True,
                    help='below threshold adoption rate ')
                    
#parser.add_argument('-p2', '--p2', type=float, required=True,
                    #help='above threshold adoption rate ')
                    
parser.add_argument('-perc', '--perc', type=float, required=True,
                    help='percent rewiring')
                    
parser.add_argument('-seed_strat', '--seed_strat', type=str, required=True,
                    help='String indicating seeding strategy to be used, \
                    can be random, adj, or nbrs')
                    
parser.add_argument('-sig_thresh_sd', '--sig_thresh_sd', type=float, required=True,
                    help='heterogeneity')
                    
parser.add_argument('-rand', '--rand', type= str, required=True,
                    help='random network construction')

args = parser.parse_args()


# Define network degreen size, and network type
k = args.k
n_nodes_approx = args.n
    
if args.G_name == "WS":
    G = nx.watts_strogatz_graph(n_nodes_approx,k,0)
    
elif args.G_name == "MR" and k ==8:
    m = round((n_nodes_approx)**(0.5))
    G = hf.moore_lattice(m,periodic= True, create_using=None)
    mappingG = {n:i for i, n in enumerate(G.nodes())}
    G = nx.relabel_nodes(G, mappingG)
    
elif args.G_name == "HX" and k == 6:
    m = round((2*n_nodes_approx)**(0.5))
    G = nx.triangular_lattice_graph(m,m, periodic = True)
    mappingG = {n:i for i, n in enumerate(G.nodes())}
    G = nx.relabel_nodes(G, mappingG)
    
else:
    print("No graph for k specified")
    
if args.seed_strat == "random":
    universal_seed_strat = hf.seed_strat_random
if args.seed_strat == "adj":
    universal_seed_strat = hf.seed_strat_adj
if args.seed_strat == "one_nbrs":
    universal_seed_strat = hf.seed_strat_one_nbrs


# Run simulation
hf.main(G, 
        G_name = args.G_name, 
        k = args.k, 
        model = hf.sirsic_wmem , 
        trials = args.trials, 
        r_start = 0, 
        beta = args.b, 
        perc = args.perc, 
        rand_path = None, 
        seeds = args.thrshld, 
        seed_strat0 = universal_seed_strat,
        seed_strat = universal_seed_strat,
        sig_thresh_sd = args.sig_thresh_sd,
        thrshld = args.thrshld, 
        p1 = args.p1, 
        #p2_list = [args.p2], 
        p2_list = np.arange(args.p1, 1.02, step = 0.02),
        cf_file = "/home/wan.a/complex_contagion_repo/output_data/cf/cf_" +
                    str(args.G_name) + "_"+
                    str(args.k) +"k_"+
                    str(args.n) +"n_"+
                    str(args.b) +"T_"+
                    str(args.thrshld) +"td_"+
                    str(args.p1) +"p1_"+
                    str(args.perc) +"perc_"+
                    str(args.sig_thresh_sd) +"sd_"+
                    str(args.rand) +"rand"+
                    ".csv", 
        df_file = None,
        full_series = False,
        rand_network = args.rand)


