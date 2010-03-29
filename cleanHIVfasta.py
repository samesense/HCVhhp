#---------------------------------------
#
# Author     Perry Evans
#            evansjp@mail.med.upenn.edu
# 2008
#---------------------------------------
"""
Remove the dashes and *s from the multiple alignment.
Print a fasta file.  Generate 2 maps, one from old positions to new,
and one from new positions to old. X's are left.  If an X causes a
gapped column, it should have been removed before this step.
Make maps from old2new and new2old sequences.  Maps are 1-indexed
b/c the ELM hits are.
"""
import sys, pickle
import utils_fasta, utils_scripting, utils_graph

req_args = ['multiple alignment to clean',
            'viral protein to add as suffix to seq names',
            'output file for protein ls',
            'new 2 old file',
            'old2new file']
examples = ['../../Data/FASTA/HIV-1/Subtypes_B_C/NEF.BC.dirty.fasta',
            'NEF',
            '../../Data/FASTA/HIV-1/Subtypes_B_C/NEF.proteins',
            '../../Data/FASTA/HIV-1/Mappings/NEF.BCtoclean.new2old.pickle',
            '../../Data/FASTA/HIV-1/Mappings/NEF.BCtoclean.old2new.pickle']
utils_scripting.checkStart(sys.argv, req_args, examples, len(req_args), True)

vp_suffix = sys.argv[2]
protein_ls_output_file = sys.argv[3]
new2old_file = sys.argv[4]
old2new_file = sys.argv[5]

def clean(seq):
    """ Remove dashes and stars. """

    return seq.replace('-', '').replace('*', '').lstrip().strip()

protein_ls = {}
new2old = {}
old2new = {}

old_fasta = utils_fasta.loadFASTA(sys.argv[1])
for protein in old_fasta.keys():
    new_seq = clean(old_fasta[protein])
    name = protein
    if name.split('.')[-1] != vp_suffix:
        name = name + '.' + vp_suffix
    utils_fasta.prettyPrint(name, new_seq)
    protein_ls[name] = True
    new2old[name] = {}
    old2new[name] = {}
    new_index = 0
    for old_index in xrange(len(old_fasta[protein])):
        old_residue = old_fasta[protein][old_index]
        if old_residue  == '-' or old_residue == '*':
            pass
        else:
            new2old[name][new_index+1] = old_index+1
            old2new[name][old_index+1] = new_index+1
            new_index += 1
f = open(new2old_file, 'w')
pickle.dump(new2old,f)
f.close()    
f = open(old2new_file, 'w')
pickle.dump(old2new,f)
f.close()    
utils_graph.dumpNodes(protein_ls_output_file, protein_ls)
