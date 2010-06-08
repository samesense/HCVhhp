import utils_humanVirus, sys, utils_scripting, utils_graph

req_args = ['elm 2 cd file',
            'cd to name file',
            'directory for mylists',
            'outfile']
examples = ['../../Data/Binding_Relations/ELM.ProfileScan.pairs',
            '../../Data/prosite.id2name',
            '../../Data/ProfileScan/MyLists/'
            'outfile']
utils_scripting.checkStart(sys.argv, req_args, examples, len(req_args), True)

elm2prosite = utils_humanVirus.get_elm2prosites(sys.argv[1])
cd2nameFile = sys.argv[2]
mylistdir = sys.argv[3]
outfile = sys.argv[4]
cd2name = {}
f = open(cd2nameFile)
for line in f:
    [id, name] = line.strip().split('\t')
    cd2name[id] = name
f.close()

elms = elm2prosite.keys()
elms.sort()
with open(outfile, 'w') as f:
    f.write('ELM\tBinding PROSITE or Entrez Gene IDs\n')
    for elm in elms:
        for cd in elm2prosite[elm].keys():
            if cd2name.has_key(cd):
                f.write(elm + '\t' + cd2name[cd] + '\n')
            else:
                genes = utils_graph.getNodes(mylistdir + cd)
                genes_to_print = ''
                for gene in genes.keys():
                    genes_to_print = genes_to_print + gene + ';'
                f.write(elm + '\t' + genes_to_print.strip(';') + '\n')
