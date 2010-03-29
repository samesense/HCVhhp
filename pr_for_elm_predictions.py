""" Given results file:
    vp elm hp predType,
    give precision and recall
    for each vp and all.
"""
import utils_scripting, utils_humanVirus, utils_graph, sys, utils_stats

req_args = ['niaid triplet file',
            'prediction file',
            'human proteins in study']
examples = ['../../Runs/Clustering.domain.s/all_niaid_triplets',
            '../../Runs/Conservation70_Cutoff.2_Window10',
            '../../Data/human.hprd.prosite']
utils_scripting.checkStart(sys.argv, req_args, examples, len(req_args), True)

hhe_vp2hp = utils_humanVirus.loadHHETargetPairs(sys.argv[1])
pred2vp2hp = utils_humanVirus.loadPredictions_predType2vp2hp(sys.argv[2])
all_hps = utils_graph.getNodes(sys.argv[3])
print 'START'
print 'Prediction Type\tVP\tHHE\tHHP\tMatch\tPrecsion\tRecall\tRandomPrecision\tPval'
for predtype in pred2vp2hp.keys():
    for vp in pred2vp2hp[predtype].keys():
        if hhe_vp2hp.has_key(vp):
            hhe = utils_graph.intersectLists([hhe_vp2hp[vp], all_hps])

            hhe_len = len(hhe.keys())
            preds = pred2vp2hp[predtype][vp]
            preds_len = len(preds.keys())
            match = utils_graph.intersectLists([hhe, preds])
            match_len = len(match.keys())
            precision = int(round(float(100)*float(match_len)/float(preds_len)))
            if hhe_len > 0:
                recall = int(round(float(100)*float(match_len)/float(hhe_len)))
            else:
                recall = 'NA'
            random_precision = int(round(float(100)*float(hhe_len)/float(len(all_hps.keys()))))
            if match_len != 0:
                pval = utils_stats.prob3(len(all_hps.keys()), preds_len, hhe_len, match_len)
            else:
                pval = 'No Matches'
            print predtype + '\t' + vp + '\t' \
                  + str(hhe_len) + '\t' \
                  + str(preds_len) + '\t' \
                  + str(match_len) + '\t' \
                  + str(precision) + '\t' \
                  + str(recall) + '\t' \
                  + str(random_precision) + '\t' \
                  + str(pval)
