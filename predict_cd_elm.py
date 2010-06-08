""" Given a list of ELM/ELM modules expected on H2,
     find H1 with the right CD.  
     Print H1, H2, H1 U H2, and H2_noRestrictions,
     which is every protein with the right ELM.
     Also print to error file that has
     vp h1 h2.
"""
import warnings
warnings.simplefilter('ignore')
import utils_scripting, sys, utils_motif, utils_humanVirus,utils_graph

req_args = ['virus elms/elm modules to use',
            'annotation type',
            'annotations of human elms/elm modules',
            'annotation type',
            'prosite human annotations',
            'annotation type',
            'network',
            'human proteins in this study',
            'version 2 gene id file',
            'elm 2 cd relations file',
            'out1',
            'out2']
examples = ['../../Data/ELM/HIV-1/Subtypes_B_C/HIV1.BC.70.conserved',
            'ELM',
            '../../Data/ELM/Human/human.website.elm',
            'ELM',
            '../../Data/ProfileScan/all.ProfileScan.scanHPRD.notNCBI',
            'ProfileScan',
            '../../Data/Network/Human/HPRD/hprd.intr',
            '../../Data/human.hprd.prosite',
            '../../Data/Network/Human/HPRD/version2entrezgeneid',
            '../../Data/Binding_Relations/ELM.ProfileScan.pairs',
            'some out 1',
            'some out 2']
utils_scripting.checkStart(sys.argv, req_args, examples, len(req_args), True)

virus_elm2protein = utils_motif.annotation2protein(sys.argv[1],
                                                   {sys.argv[2]:True})
study_hps = utils_graph.getNodes(sys.argv[8])
human_elm2protein = utils_motif.annotation2protein_forProteins(sys.argv[3],
                                                               {sys.argv[4]:True},
                                                               study_hps)
human_cd2protein = utils_motif.annotation2protein_forProteins(sys.argv[5],
                                                              {sys.argv[6]:True},
                                                              study_hps)
network = utils_graph.getEdges(sys.argv[7])
version2geneid = utils_humanVirus.get_version2entrez(sys.argv[9])
elm2cd = utils_humanVirus.get_elm2prosites(sys.argv[10])
outf1 = sys.argv[11]
outf2 = sys.argv[12]

vp_to_h1_to_h2 = {}
with open(outf1, 'w') as f:
    for elm in virus_elm2protein.keys():
        if human_elm2protein.has_key(elm):
            h2_noRestrictions = human_elm2protein[elm]
            h2 = {}
            h1 = {}
            h1_to_h2 = {}
            elms_in_module = elm.split(':')
            matching_cds = {}
            for e in elms_in_module:
                if elm2cd.has_key(e):
                    for cd in elm2cd[e].keys():
                        if human_cd2protein.has_key(cd):
                            matching_cds[cd] = True

            for hp in h2_noRestrictions.keys():
                for hp_neigh in network[hp].keys():
                    for cd in matching_cds.keys():
                        if human_cd2protein[cd].has_key(hp_neigh):
                            h1[hp_neigh] = True
                            h2[hp] = True
                            h1_to_h2[hp_neigh + ':' + hp] = True
            for vp in virus_elm2protein[elm].keys():
                for pred in h1.keys():
                    if version2geneid.has_key(pred):
                        f.write(vp + '\t' + elm + '\t'
                                + version2geneid[pred] + '\th1\n')
                for pred in h2.keys():
                    if version2geneid.has_key(pred):
                        f.write(vp + '\t' + elm + '\t'
                                + version2geneid[pred] + '\th2\n')
                for pred in utils_graph.unionLists([h1,h2]).keys():
                    if version2geneid.has_key(pred):
                        f.write(vp + '\t' + elm + '\t'
                                + version2geneid[pred] + '\th1h2\n')
                for pred in h2_noRestrictions.keys():
                    if version2geneid.has_key(pred):
                        f.write(vp + '\t' + elm + '\t' 
                                + version2geneid[pred] + '\th2All\n')
                for pair in h1_to_h2.keys():
                    if not vp_to_h1_to_h2.has_key(vp):
                        vp_to_h1_to_h2[vp] = {}
                    vp_to_h1_to_h2[vp][pair] = True

with open(outf2, 'w') as f:
    for vp in vp_to_h1_to_h2.keys():
        for pair in vp_to_h1_to_h2[vp]:
            [h1_gene, h2_gene] = pair.split(':')
            if version2geneid.has_key(h1_gene) and version2geneid.has_key(h2_gene):
                f.write(vp + '\t' + version2geneid[h1_gene]
                        + '\t' + version2geneid[h2_gene] + '\n')
            

