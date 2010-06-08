#!/usr/bin/env python

"""For each HCV protein, calcuate the likelyhood
   of the GO BP similarity between predictions
   and gold standard. Do this for H1H2 & H1.
"""
import sys, utils_stats, utils_graph, utils_humanVirus, random, os

hhe_file = sys.argv[1]
hhp_file = sys.argv[2]
background_file = sys.argv[3]
out_file = sys.argv[4]

# this takes a long time
# utils_stats.gene_set_go_sim(background_file, 'results/HPRD.ls.entrez.gosim')

hhe_vp2hp = utils_humanVirus.loadHHETargetPairs(hhe_file)
pred2vp2hp = utils_humanVirus.loadPredictions_predType2vp2hp(hhp_file)
all_hps = utils_graph.getNodes(background_file)

for pred_type in ('h1', 'h1h2'):
    for vp in pred2vp2hp[pred_type].keys():
        if hhe_vp2hp.has_key(vp):
            hhe = utils_graph.intersectLists([hhe_vp2hp[vp], all_hps]).keys()
            preds = pred2vp2hp[pred_type][vp].keys()
            go_pval = utils_stats.gene_set_go_sim_pval(preds, hhe,
                                                       'results/HPRD.ls.entrez.gosim')
            print('%s\t%s\t%.3f' %
                  (vp, pred_type, go_pval))
