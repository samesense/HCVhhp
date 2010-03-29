""" Given a directory w/ hand 
    currated CD genes, convert
    the CD files there to 
    annotations in the prosite
    form.
"""
import utils_scripting, sys, os, utils_humanVirus

req_args = ['mylist directory',
            'translation file']
examples = ['../../Data/ProfileScan/MyLists/',
            'some translation file']
utils_scripting.checkStart(sys.argv, req_args, examples, len(req_args), True)

mylist_dir = sys.argv[1]
trans = utils_humanVirus.entrez2version(sys.argv[2])

for file in os.listdir(mylist_dir):
    f = open(mylist_dir + file)
    for line in f:
        entrez = line.strip()
        for version in trans[entrez]:
            print(version + '\t0\t0\t'
                  + file + '\tseq\tProfileScan')
    f.close()
