# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import sys
import networkx as nx
import networkx.algorithms.approximation as nx_approx
import pandas as pd
import numpy as np
#from numexpr3 import evaluate as ev
from scipy.spatial.distance import cdist, pdist, squareform
import scipy.integrate as integrate
import time
import math
import os
import csv

from numexpr import evaluate as ev
from tqdm import tqdm

def giulia_spatial_entropy(edgelist, nodelist):

    start = time.time()

    g = nx.read_edgelist(edgelist, delimiter=',')
    nodes = list(g.nodes())
    nodes.sort()

    PP = np.array(nx.to_numpy_matrix(g,nodelist=nodes))

    n = len(nodes)
    print("I'm here...Number:", str(n))
    distance = np.zeros((n,n))

    ens_id_expr_df = pd.read_csv(nodelist,names=['ens_id','expr_val'])
    ens_id_expr_df = ens_id_expr_df.set_index('ens_id')
    
    ens_id_expr_map = ens_id_expr_df.to_dict()['expr_val']

    for row in tqdm(range(n)): # first tqdm
        for col in range(n):
            #try:
            x = ens_id_expr_map[nodes[row]]
            y = ens_id_expr_map[nodes[col]]
            distance[row][col] = math.fabs(x-y)

            #except KeyError:
                #distance[row][col] = 0
                
                

    Nbin = int(math.sqrt(n)+1)

    print("before linear binning...")
    #linear binning
    mi,ma = np.min(distance[distance>0]),np.max(distance)
    limiti = np.linspace(mi,ma,Nbin+1)
    limiti[-1]+=0.1*ma
    limiti[0]-=0.1*mi

    b = np.searchsorted(limiti,distance)-1

    print("before massimo...")
    massimo = np.max(b)+1
   
    # BC gives how many links fall in each bin
    BC = [ np.sum(b[PP>0]==i)/2 for i in range(massimo) ]


    N=PP.shape[0];
    
    connectivity = np.sum(PP,1);
    avg_conn = np.mean(connectivity)

    print("before lagragian...")
    #Lagrangians
    z = connectivity/(math.sqrt(avg_conn*N))
    w = BC / (avg_conn*N)

    old_z = z
    old_w = w

    loops = 5  # CHANGE to 10000 again
        
    precision = 1E-5
    #precision = 1E-3 # CHANGED NOW!

    for idx in tqdm(range(loops)): # second tqdm
        bigW = w.take(b)

        for i in range(N):  bigW[i,i]=0.

        U = ev("bigW * z")
        UT = U.T    
        D = ev("(UT * z) + 1.")

        UD = ev("U/D")

        del D,U,UT
        


        for i in range(N):  UD[i,i]=0.

        z = connectivity / np.sum(UD,1)

        zt = z[:].reshape(N,1)
        D2 = ev("(z*zt) / ( bigW *(z*zt) + 1.)")

        B2 = np.array([np.sum(np.where(b==i,D2,0.)) for i in range(len(w))])/2.

        print("And calculating B2 AND D2 done!!!!! inside for loop out of 5")

        w = np.where( (BC!=0) & (B2!=0),BC/B2,0)
        rz= ev("abs(1.-z/old_z)")
        rw= ev("abs(1.-w/old_w)")
        rz[np.isinf(rz)]=0.
        rw[np.isinf(rw)]=0.
        rz[np.isnan(rz)]=0.
        rw[np.isnan(rw)]=0.

        if max(rz)<precision and max(rw)<precision:
            break

        old_z = z
        old_w = w
        
    print("JUST OUT OF THE BIG FORR LOOP...!")

    bigW = w.take(b)
    for i in range(N):  bigW[i,i]=0.

    z2 = bigW * np.outer(z,z)
    P = z2 / ( z2 + 1)
    Q=1.-P
    
    print("And calculations done!!!!!")
    S = -1*(np.sum(np.log(P**P) + np.log(Q**Q) ))/2.
#    S = -1*(np.sum((P*np.log(P)) + (Q*np.log(Q) )))/2.
    print("Done S  ....!!!!!")
    stop = time.time()

    output = dict()
    output['num_nodes'] = n
    output['num_edges'] = len(g.edges())
    output['giulia_spatial_entropy'] = S
    output['runtime'] = stop-start
    output['nodelist'] = nodelist
    output['edgelist'] = edgelist

    return output

#groups = ['FA_youngB','FA_oldB','FA_gerB']
group='FA_'

#for group in groups:

print(group)
#print(constant)

## have the filled input folders ready
#nodelist_dir = group + '_nodelists/'
#edgelist_dir = group + '_edgelists/'
nodelist_dir = 'nodelists/'
edgelist_dir = 'edgelists/'
   
##gets all the contents of "nodelist_dir" and gets the cellID
all_nodelists = os.listdir(nodelist_dir)
all_ids = [x.split('_')[0] for x in all_nodelists]

#number of cells in each nodelist folder
total = len(all_ids)

#creating a new csv file to write computed spatial entropy values and giving headers

for cell_id in all_ids:
    print('Need to complete: ' + str(total))
    count = 1
    with open(group +'_'+ cell_id+ '_spatial_entropy.csv', 'w') as csvfile:

        fieldnames = ['nodelist','edgelist','num_nodes','num_edges', \
                    'giulia_spatial_entropy','runtime']
    
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
    
        #total number of cells

    #for cell_id in all_ids:
        
         ## you should comment this out later!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        print('Working on: ' + str(count) + '/' + str(total))
    
    #creating one csv for one cell - nodelist AND edgelist
        edgelist = edgelist_dir + cell_id + '_edgelist_' + group.split('_')[0] + '_.csv'
        nodelist = nodelist_dir + cell_id + '_nodelist_' + group.split('_')[0] + '_.csv'

        myoutput=giulia_spatial_entropy(edgelist,nodelist)
        writer.writerow(myoutput)
        
        print(cell_id)
        print(myoutput)        
        print("ROW PRINTED!")
        
        count += 1

print('\n\n') 

#    
#        cell_id='FA-01'
#    
#     ## you should comment this out later!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#    print('Working on: ' + str(count) + '/' + str(total))
#
##creating one csv for one cell - nodelist AND edgelist
#    edgelist = edgelist_dir + cell_id + '_edgelist_' + group.split('_')[0] + '_.csv'
#    nodelist = nodelist_dir + cell_id + '_nodelist_' + group.split('_')[0] + '_.csv'
#
#    writer.writerow(giulia_spatial_entropy(edgelist,nodelist))
#    
#    print("ONE ROW PRINTED!")
#    count += 1
#
#print('\n\n')
