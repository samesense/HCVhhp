#---------------------------------------
#
# Author     Perry Evans
#            evansjp@mail.med.upenn.edu
# 2008
#---------------------------------------
"""
Convert HHP results in version #s to 
Entrez Gene IDs so they can be compared
to the NIAID database.
"""
import utils_scripting, sys, utils_graph

req_args = ['protein ls file to convert',
            'version2entrez translation file']
examples = ['../../Runs/Restrict_Domain_By_ELM_Restrict_ELM_By_Domain_Whole_Virus/H12',
            '../../Data/Network/Human/HPRD/version2entrezgeneid']
utils_scripting.checkStart(sys.argv, req_args, examples, 2, True)

translation_file = sys.argv[2]

proteins_version = utils_graph.getNodes(sys.argv[1])
version2entrez = {}
f_translation = open(translation_file)
for line in f_translation:
    [version, entrez_gene] = line.strip().split('\t')
    version2entrez[version]= entrez_gene
f_translation.close()
proteins_entrez = {}
for protein in proteins_version.keys():
    if version2entrez.has_key(protein):
        proteins_entrez[ version2entrez[protein] ] = True
for protein in proteins_entrez:
    print protein
#utils_graph.dumpNodes(output_file, proteins_entrez)
