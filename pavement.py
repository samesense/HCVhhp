from paver.easy import *
"""Tasks for the project"""

@task
def HCV_seq_setup():
    """ Parse LANL HCV fasta files """
    
    hcv_proteins = ['ARFP', 'CORE', 'E1', 'E2',
                    'NS2', 'NS3', 'NS4A', 'NS4B',
                    'NS5A', 'NS5B', 'P7']:
    for hcv_protein in hcv_proteins:
        if hcv_protein == 'ARFP':
            suffix = 'F'
        else:
            suffix = hcv_protein
        sh('python cleanHIVfasta.py ' 
           + 'data/LANL/HCV_ALL_2008_'
           + hcv_protein
           + '_PRO.fasta '
           + suffix + ' '
           + 'data/'
           + suffix
           + '.proteins > '
           + 'data/'
           + suffix
           + '.clean.fasta')
    sh('cat data/*.clean.fasta '
       + '> data/HCV.fasta')
    sh("grep '>' "
       + 'data/HCV.fasta '
       + '> data/HCV.proteins')

@task
def hcv_elms():
    """Annotate HCV multiple alignments from LANL"""
    
    for hcv_protein in ('F', 'CORE', 'E1', 'E2',
                        'NS2', 'NS3', 'NS4A', 'NS4B', 
                        'NS5A', 'NS5B', 'P7'):
        sh('wget ' + protein + ' --output-document='
           + protein + 'fa')
