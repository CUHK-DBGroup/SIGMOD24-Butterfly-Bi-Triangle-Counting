import networkx as nx
import time
import random
import utils
import numpy as np
from networkx.algorithms import bipartite
# G = nx.davis_southern_women_graph()
G = nx.Graph()
L = 0 # set L size
R = 0 # set R size
edge_num =  0
G.add_nodes_from(list(range(1,L+1)), bipartite=0)
G.add_nodes_from(list(range(L+1,L+R+1)), bipartite=1)
E_list = []
edge_file = open("xx.txt") # input dataset from xx.txt
for i in edge_file.readlines():
    u,v = map(int,i.split())
    v += L
    edge_num += 1
    E_list.append((u,v))
    G.add_edge(u,v)

probs = []
for i in range(1, L+1):
    probs.append(len(G[i]))
    
sampler = utils.discreteSampler(probs)
bias = edge_num * edge_num
for i in range(1,L+1):
    bias -= len(G[i]) * len(G[i])
bias/=2

def Est_RMSE(est_res):
    k = len(est_res)
    res = 0
    for i in est_res:
        res += i
    res /= k
    res2 = 0
    for i in est_res:
        res2 += (i-res) * (i-res)
    res2 /= (k-1)
    return np.sqrt(res2)

def est_4_bf(sample_num):
    num_4_bf = 0
    for i in range(0,sample_num):
        u = sampler.generate()
        v = sampler.generate()
        while u == v:   
            u = sampler.generate()
            v = sampler.generate()
        inter_uv = len(set(G[u])&set(G[v]))
        if inter_uv>=2:
            num_4_bf += (inter_uv *(inter_uv-1) /2)/len(G[u])/len(G[v])
    return num_4_bf / sample_num * bias

def est_3_path(sample_num):
    num_3_path = 0
    for i in range(0,sample_num):
        u,v = random.choice(E_list)
        if G[u] and G[v]:
            num_3_path += (len(G[u])-1) * (len(G[v])-1)
    return num_3_path * edge_num / sample_num

time_start=time.time()
rnd_bf = 0 # set round
rnd_path = 0 # set round
est_coe = 4*est_4_bf(rnd_bf)/est_3_path(rnd_path)
true_coe = 1.0 # set true coe
time_end=time.time()
print('time cost',time_end-time_start,'s')
print((est_coe - true_coe)/true_coe)
# one can calculate the reasonable estimate of RMSE by the results of estimations