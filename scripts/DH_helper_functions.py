# 네트워크 전파 시뮬레이션 핵심 함수 모음 (네트워크 생성, SIR 확산 모델, 시딩 전략, 후처리)
import networkx as nx
import random
import numpy as np
import pandas as pd
import os
import math

###################
# Network Construction

def moore_lattice(m,periodic=False, create_using=None):
    """
    constructs moore lattice where every node has degree k
    m^2 = number of nodes
    """
    G = nx.empty_graph()

    rows = list(range(m))
    cols = list(range(m))

    row_pairs = []
    for i in rows[0:-1]:
        pair = (i,i+1)
        row_pairs.append(pair)

    col_pairs = []
    for i in cols[0:-1]:
        pair = (i,i+1)
        col_pairs.append(pair)

    G.add_nodes_from((i, j) for i in rows for j in cols)
    G.add_edges_from(((i, j), (pi, j)) for pi, i in row_pairs for j in cols)
    G.add_edges_from(((i, j), (i, pj)) for i in rows for pj, j in col_pairs)

    G.add_edges_from(((i, j), (i+1, pj)) for i in rows[0:-1] for pj, j in col_pairs)
    G.add_edges_from(((i+1, j), (i, pj)) for i in rows[0:-1] for pj, j in col_pairs)
    
    try:
        periodic_r, periodic_c = periodic
    except TypeError:
        periodic_r = periodic_c = periodic

    if periodic_r and len(rows) > 2:
        first = rows[0]
        last = rows[-1]
        G.add_edges_from(((first, j), (last, j)) for j in cols)
        G.add_edges_from(((first, j), (last, j+1)) for j in cols[0:-1])
        G.add_edges_from(((first, j+1), (last, j)) for j in cols[0:-1])

    if periodic_c and len(cols) > 2:
        first = cols[0]
        last = cols[-1]
        G.add_edges_from(((i, first), (i, last)) for i in rows)
        G.add_edges_from(((i, first), (i+1, last)) for i in rows[0:-1])
        G.add_edges_from(((i+1, first), (i, last)) for i in rows[0:-1])

    if periodic == True:
        G.add_edge((rows[0], cols[0]), (rows[-1], cols[-1]))
        G.add_edge((rows[0], cols[-1]), (rows[-1], cols[0]))

    # both directions for directed
    if G.is_directed():
        G.add_edges_from((v, u) for u, v in G.edges())
    return G

###################
# Seeding Strategies

def seed_strat_one_nbrs(G, n, seeds):
    """
    pick one random node and x other randomly chosen neighbors for that node
    """
    seed_node_index = random.randint(0, n-1)
    seed_node_nbr = list(G.neighbors(seed_node_index))
    random.shuffle(seed_node_nbr)
    seed_list = seed_node_nbr[0:seeds-1]
    seed_list.append(seed_node_index)
    return seed_list

def seed_strat_random(G, n, seeds):
    """
    pick x random seeds
    """
    seed_list = random.sample(list(range(n)),seeds )
    return seed_list

def seed_strat_adj(G, n, seeds):
    """
    pick x adjacent seeds
    """
    seed_list = list(range(0,seeds))
    return seed_list
    
###################
# Intermediary processing

def get_cum_freq(node_dict, timesteps):
    """
    Get Cummulative adoption over time.

    'node_dict' is a dictionary where key = node_id and 
    val = dictionary of node attributes.
    Produced from running the simulation.
    'timesteps' is the total number of timesteps for the simulation
    
    Returns dictionary where every key is timestep and 
    every value is the number of cummulative adopters at that time.
    """
    freq_dict = {i:0 for i in range(timesteps+1)}
    for key, val in node_dict.items():
        if val['adoption_time']== 'NA':
            continue
        freq_dict[val['adoption_time']]= freq_dict[val['adoption_time']]+1

    cum_freq = {0:freq_dict[0]}
    for i in range(1, timesteps+1):
            #print(freq_dict[val[1]])
        cum_freq[i]= freq_dict[i]+cum_freq[i-1]
    return cum_freq

def get_t_spread(cf_list, spread_amount,n):
    """
    Returns number of time steps it takes for simulation to saturate a specified
    amount (spread_amount) of the network
    """
    target_n = spread_amount*n
    if cf_list[-1] < target_n:
            i ="NA"
    else:
        i = 0
        while cf_list[i] < target_n:
            i = i+1
        
    return i
    
def get_mc_sigmoid(x, k, m):
    """
    Sigmoid function used in Centola and Macy 2007.
    """
    denom = 1+(math.exp((k -x)*m))
    p = 1/denom
    return p

###################
# Main simulation

def sirsic_wmem(G,beta, prob_list, 
                seeds,seed_strat0,seed_strat, perc, sig_thresh_sd, random_seed):
    """
    Models one diffusion process for a given network, adoption trajectory, 
    number of seeds, seeding strategy, rewiring level, and adoption heterogeneity
    Diffucion process runs until there are either no infective of no susceptible 
    individuals left.
    
    Parameters:
        G: networkx graph object
            Graph to diffuse behavior on
        beta: integer
            Time steps adopters remain influential for before going dormant. 
            Set beta = 0 for adopters to remain influential for the entire simulation
        prob_list: list of probabilisties corresponding to adoption trajectory
            Ex. prob_list = [0,0.1, 0.5, 0.5, 0.5,...] means p_0 = 0, p_1 = 0.1,
            p_2 = 0.5 and so on.
        seeds: integer
            NUmber of seeds to start simulation
        seed_strat0: function name of seeding strategy used when rewiring = 0 
            Can be seed_strat_one_nbrs, seed_strat_random, or seed_strat_adj
        seed_strat0: function name of seeding strategy used when rewiring > 0 
            Can be seed_strat_one_nbrs, seed_strat_random, or seed_strat_adj
        perc: Float between 0 and 1 inclusive
            Amount that original network is rewired. Used to determine seeding 
            strategy.
        sig_thresh_sd: Float between 0 and 1 inclusive
            Standard deviation of normal distribution centered at adoption prob.
            sig_thresh_sd = 0 means homogeneous within individual adoption
            sig_thresh_sd > 0 increases heterogeneity
        random_seed: integer
            sets the trial number and random seed
            
    Returns: dictionary of dictionries
        Keys are node ids
        Values are dictionaries with node level information including whether the
        node adopted, the time of adoption, how many signals from different 
        adopting neighbors, and a list of influential neighbors they are in contact with
    """
    random.seed(random_seed)
    
    # Initialize node dictionary
    node_list = list(G.nodes)

    node_dict = {n:{"adopted": 0,"adoption_time" : "NA", 
                    "signals": 0,  "influencers":[]} for n in node_list}
    n = len(node_list)
    
    # Pick seeding strategy based on levels of rewiring
    if perc == 0:
        seed_list = seed_strat0(G, n, seeds)
    else:
        seed_list = seed_strat(G, n, seeds)
    
    i_list = [seed_list]
    s_list = [list(G.nodes)]    
    for seed in seed_list:
        node_dict[seed]={"adopted": 1,"adoption_time" : 0, 
                         "signals": 0, "influencers":[]}
        s_list[0].remove(seed)
    
    # If beta = 0, adopted remain influential for the entire simulation
    # else, that are only influential for beta time steps
    if beta == 0:
        i_list1 = i_list
    else:
        i_list1 = i_list[-beta::] 
    #flatten list
    ii_list =  [item for sublist in i_list1 for item in sublist]

    t = 1
    while len(ii_list) > 0:
        new_i_list =[]
        s_remove = []
        for s in s_list[t-1]:
            for nb in list(G.neighbors(s)):
                if nb in ii_list:
                    if nb not in node_dict[s]['influencers']:
                        node_dict[s]['signals']=node_dict[s]['signals']+1
                    #transmission rule
                    node_dict[s]['influencers'].append(nb)
                    signal_count = node_dict[s]['signals']
                    # draw probability from normal dist. 
                    sig_thrshld = np.random.normal(prob_list[signal_count],  
                                                   scale=sig_thresh_sd, 
                                                   size=None)
                    x = np.random.uniform(0,1)
                    if x < sig_thrshld:
                        node_dict[s]['adopted']=1
                        node_dict[s]['adoption_time']=t
                        new_i_list.append(s)
                        s_remove.append(s)
                        break
        
        i_list.append(new_i_list)
        new_s_list = [i for i in s_list[t-1] if i not in s_remove]
        
        if len(new_s_list) < 1 :
            break
        random.shuffle(new_s_list)
        s_list.append(new_s_list)

        #get infected from last x timesteps
        if beta == 0:
            i_list1 = i_list
        else:
            i_list1 = i_list[-beta::] 
        
        ii_list =  [item for sublist in i_list1 for item in sublist]

        t = t+1
        if t > 10000:
            break

    return node_dict

###################
# Wrapper Function to run many simulations and post-process

def run_simulation(G, G_name, model, trials, 
                   r_start, beta, prob_list, perc, 
                   rand_path,seeds,seed_strat0,seed_strat, sig_thresh_sd,
                  full_series, rand_network):
    '''
    Runs diffusion model for specified number of trials and post processes. 
    
    Parameters:
        G, beta, prob_list, perc, rand_path, seeds, seed_strat0, seed_strat, 
        sig_thresh_sd are the same as in sirsic_wmem.
        G_name: string
            User defined name to identify what network type is used
        model: function
            Function name diffusion model is used. In this case, 
            model = sirsic_wmem
        trials: integer
            Number of simulations to run with the specified parameters
        r_start: integer
            Which trial number to start at. 
        full_series: Boolean
            If True, saves the adoption time series (list of the number of nodes
            adopting at the first time step, second, and so on) to a column in the
            output data frame.
            If False, does not save
        rand_network: string
            User defined name to identify what proecedure is used to create 
            random network if perc > 0.
    Returns:
        df: Node level data frame from each trial in one data frome
        cf: Trial level data frame aggregated across all nodes
    '''
    n_edges = G.number_of_edges()
    n_nodes = len(G)

    cum_freq_list = []
    df_list= []

    # for each trial
    for r in range(r_start, trials):
        # rewire the network if needed
        R = G.copy()
        if perc > 0:
            niter = perc*n_edges
            if rand_path is None:
                nx.connected_double_edge_swap(R, nswap=niter, seed=r)
            else:
                saved_name = rand_path+G_name+"_"+str(r)+".txt"
                if os.path.exists(saved_name):
                    R = nx.read_edgelist(saved_name, data= False, delimiter = ' ')
                else:
                    nx.connected_double_edge_swap(R, nswap=niter, seed=r)
                    nx.write_edgelist(R, saved_name, data = False, delimiter = ' ')
        # run diffusion
        node_dict = model(R, beta, prob_list, 
                seeds,seed_strat0,seed_strat, perc, sig_thresh_sd, r)
                
        # create node-level data frame from dictionary, add additional meta data
        node_df = pd.DataFrame.from_dict(node_dict, orient='index')
        node_df.index.name = 'node_id'
        node_df.reset_index(inplace=True)
        adoption_time = [i for i in list(node_df['adoption_time']) if i != "NA"]
        adoption_time = [int(i) for i in adoption_time]

        node_df["seed"] = r
        node_df["n_nodes"] = n_nodes
        node_df["n_edges"] = n_edges
        node_df["perc_rewire"] = perc
        node_df["time_of_influence"] = beta
        node_df["sig_thresh_sd"] = sig_thresh_sd
        if perc==0:
            node_df["seed_strat"]  = seed_strat0.__name__
        else:
            node_df["seed_strat"]  = seed_strat0.__name__
        df_list.append(node_df)
        
        # create trial-leve; data frame from dictionary, add additional meta data
        cum_freq = get_cum_freq(node_dict, max(4,max(adoption_time)))
        cf_list = [val for key,val in cum_freq.items()]
        
        cum_freq_dict = {
            "spread_time_0": cf_list[0],
            "spread_time_1":  cf_list[1],
            "spread_time_2":  cf_list[2],
            "spread_time_3":  cf_list[3],
           # "full_timeseries": cf_list,
            "max_spread": cf_list[-1],
            "time_to_spread": len(cf_list)-1,
            "time_to_60_spread": get_t_spread(cf_list, 0.6,n_nodes),
            "time_to_75_spread": get_t_spread(cf_list, 0.75,n_nodes),
            "time_to_90_spread": get_t_spread(cf_list, 0.9,n_nodes),
            "seed": r,
            "n_nodes": n_nodes,
            "n_edges": n_edges,
            "perc_rewire": perc,
            "n_seeds":seeds,
            "time_of_influence":beta,
            "sig_thresh_sd": sig_thresh_sd,
            #"seed_strat": seed_strat.__name__,
            "randomization_type": rand_network
        }

        if full_series == True:
            cum_freq_dict["full_timeseries"] = cf_list
        else:
            cum_freq_dict["full_timeseries"] = "Not Saved"

        if perc == 0:
            cum_freq_dict["seed_strat"] = seed_strat0.__name__
        else:
            cum_freq_dict["seed_strat"] = seed_strat.__name__
    
        cum_freq_list.append(cum_freq_dict)

    # add node-level or trial-level data to existing data frames
    cf = pd.DataFrame.from_dict(cum_freq_list, orient='columns')
    df = pd.concat(df_list)
    return df, cf

###################
# Main function

def main(G, G_name,k, model, trials, 
        r_start, beta,  perc, 
        rand_path,seeds,seed_strat0,seed_strat, sig_thresh_sd,
         thrshld, p1, p2_list, 
         cf_file, df_file, full_series, rand_network):
    """
    Main function to run simulation. Contains auxilliary functions sirsic_wmem 
    and run_simuation. Automatically saves node-level and trial level data 
    outputted from run_simulation to specified file names.
    
    Parameters:
        G, G_name,model, trials, r_start, beta,  perc, rand_path,seeds, 
        seed_strat0,seed_strat, sig_thresh_sd, full_series, rand_network 
        are the same as in sirsic_wmem and run_simulation.
        k: integer
            Degree of network (all nodes have the same degree)
        thrshld: integer
            Social reinforcement threshold. The number of exposures to different 
            neighbors required to adopt at p_2 instread of p_1
        p1:
            Base rate adoption
        p2_list: list
            List of socially reinforced adoption rates to loop through. 
            List length of one can be specified for just one p_2 value
        cf_file: string
            File name to save trial-level data
        df_file: string
            File name to save node-level data. If df_file = None, node-level data
            will not be saved.
    
    """
    n = len(G)
    for p2 in p2_list:
        prob_list = [0]+ [p1]*(thrshld - 1)+ [p2]*n

        ddf, testcf =run_simulation(G, G_name, model, trials, 
                   r_start, beta, prob_list, perc, 
                   rand_path,seeds,seed_strat0, seed_strat, sig_thresh_sd,
                                        full_series, rand_network)

        ddf["p1"] = round(p1, 4)
        ddf["p2"] = round(p2, 4)
        ddf["k"] = k
        ddf["thrshld"] = thrshld
        ddf["G_name"] = G_name
            
        testcf["p1"] = round(p1, 4)
        testcf["p2"] = round(p2, 4)
        testcf["k"] = k
        testcf["thrshld"] = thrshld
        testcf["max_spread_norm"] = testcf["max_spread"]/testcf["n_nodes"]
        testcf["G_name"] = G_name
            
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
