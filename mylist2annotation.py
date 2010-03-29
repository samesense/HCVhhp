""" Given a directory w/ hand 
    currated CD genes, convert
    the CD files there to 
    annotations in the prosite
    form.
"""
import utils_scripting, sys, os, utils_humanVirus

req_args = ['mylist directory',
            'translation file',
            'output file']
examples = ['../../Data/ProfileScan/MyLists/',
            'some translation file',
            'some out file']
utils_scripting.checkStart(sys.argv, req_args, examples, len(req_args), True)

mylist_dir = sys.argv[1]
trans = utils_humanVirus.get_entrez2version(sys.argv[2])
outfile = sys.argv[3]

with open(outfile, 'w') as outf:
    for file in os.listdir(mylist_dir):
        with open(mylist_dir + file) as f:
            for line in f:
                entrez = line.strip()
                if entrez in trans:
                    for version in trans[entrez]:
                        outf.write(version + '\t0\t0\t'
                                   + file + '\tseq\tProfileScan\n')

