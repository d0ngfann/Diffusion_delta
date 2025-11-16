# Centola and Macy 2007 복제 시뮬레이션 (SI 섹션 K 데이터 생성, 누적 로지스틱 함수 사용)
"""
Script to replicate Centola and Macy 2007
Reproduces data for SI Section K

Uses cumulative logistic function parameterized by m instead of step function

Parameter settings
k = 8
i = 2
b = 0, 1
sig_thresh_sd = 0
1D and 2D lattices
"neighbor" seeding
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
                    
parser.add_argument('-b', '--b', type=int, required=True, 
                    help='time of influence')
                    
parser.add_argument('-perc', '--perc', type=float, required=True,
                    help='percent rewiring')
                    
parser.add_argument('-m', '--m', type=float, required=True,
                    help='slope')
                    
parser.add_argument('-seed_strat', '--seed_strat', type=str, required=True,
                    help='String indicating seeding strategy to be used, \
                    can be random, adj, or nbrs')
                    
parser.add_argument('-sig_thresh_sd', '--sig_thresh_sd', type=float, required=True,
                    help='heterogeneity')
                    
parser.add_argument('-rand', '--rand', type= str, required=True,
                    help='random network construction')

args = parser.parse_args()

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
if args.seed_strat == "nbrs":
    universal_seed_strat = hf.seed_strat_one_nbrs
  
cf_file = "/home/wan.a/complex_contagion_repo/output_data/macy_centola_07/cf_"+ str(args.G_name) + "_"+ str(args.k) +"k_"+str(args.n) +"n_"+str(args.b) +"T_"+str(args.perc) +"perc_"+str(args.sig_thresh_sd) +"sd_"+str(args.rand) +"rand"+".csv"
df_file = None
i = 100
prob_list =  [hf.get_mc_sigmoid(x, 2, args.m) for x in range(i)]

ddf, testcf =hf.run_simulation(G, 
                                G_name = args.G_name, 
                                model =hf.sirsic_wmem, 
                                trials = 100, 
                                r_start = 0, 
                                beta = args.b, 
                                prob_list = prob_list, 
                                perc =  args.perc, 
                                rand_path = None,
                                seeds = 2,
                                seed_strat0 = hf.seed_strat_one_nbrs, 
                                seed_strat = hf.seed_strat_one_nbrs, 
                                sig_thresh_sd = 0,
                                full_series = False, 
                                rand_network = args.rand)
            
testcf["m"] = args.m
testcf["p1"] = prob_list[1]
testcf["p2"] = 0.5
testcf["k"] = k
testcf["thrshld"] = 2
testcf["max_spread_norm"] = testcf["max_spread"]/testcf["n_nodes"]
testcf["G_name"] = args.G_name
            
testcf.columns = testcf.columns.astype(str)

if os.path.exists(cf_file):
    existing_cf = pd.read_csv(cf_file, sep=',')
    testcf = pd.concat([existing_cf, testcf], ignore_index = True)
testcf.to_csv(cf_file, sep=',', index=False)

if df_file is not None:
    if os.path.exists(df_file):
        existing_df = pd.read_csv(df_file, sep=',')
        ddf = pd.concat([existing_df, ddf], ignore_index = True)
    ddf.to_csv(df_file, sep=',', index=False)

