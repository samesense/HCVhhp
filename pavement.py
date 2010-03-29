from paver.easy import *
"""Tasks for the project"""

@task
def HCV_seq_setup():
    """Parse LANL HCV fasta files"""
    
    hcv_proteins = ['ARFP', 'CORE', 'E1', 'E2',
                    'NS2', 'NS3', 'NS4A', 'NS4B',
                    'NS5A', 'NS5B', 'P7']
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
           + '.proteins '
           + 'data/' + suffix + '.new2old '
           + 'data/' + suffix + '.old2new > '
           + 'data/'
           + suffix
           + '.clean.fasta')
    sh('cat data/*.clean.fasta '
       + '> data/HCV.fasta')
    sh("grep '>' "
       + 'data/HCV.fasta '
       + "| sed 's/>//g' "
       + '> data/HCV.proteins')

@task
def HCV_elms():
    """Annotate HCV multiple alignments from LANL"""
    
    sh('python matchELMpattern.py '
       + 'data/elm2pattern '
       + 'data/HCV.fasta '
       + '> data/HCV.elms')
    sh('python getConserved.py '
       + 'data/HCV.elms '
       + 'ELM '
       + '90 '
       + '1> data/HCV.conserved.90 '
       + '2> data/HCV.conservation')

def get_ls():
    """net_name, trans, ps_scan, net, elms, bg"""

    ls = [ ['OPHID',
            '../Thesis/Data/Network/Human/OPHID/swissProt2entrez',
            #'../Thesis/Data/ProfileScan/OPHID.fasta.ps_scan.parsed',
            'data/OPHID.ps_scan',
            '../Thesis/Data/Network/Human/OPHID/ophid.swiss',
            '../Thesis/Data/ELM/Human/OPHID/ophid.elms',
            'OPHID.bg'],
           ['STRING',
            '../Thesis/Data/Network/Human/STRING/ensp2entrez',
#            '../Thesis/Data/ProfileScan/STRING.human.ps_scan.parsed',
            'data/STRING.ps_scan'
            '../Thesis/Data/Network/Human/STRING/human.string.intr',
            '../Thesis/Data/ELM/Human/STRING/string.elms',
            'STRING.bg']]#,
           #['HPRD',
           # '../Thesis/Data/Network/Human/HPRD/version2entrezgeneid_new',
           # '../Thesis/Data/ProfileScan/STRING.human.ps_scan.parsed',
    #'data/HPRD.ps_scan',
    #'../Thesis/Data/Network/Human/HPRD/hprd_new.intr',
#'../Thesis/Data/ELM/Human/HPRD/hprd.elms',
#HPRD.bg] ]
    return ls

@task
def OPHID_domains():
    """Find CDs for HCV ELMs for OPHID"""
    sh('python scan_prosite.py '
       + 'data/elm2prosite '
       + '~/bioperry/Projects/Thesis/Data/FASTA/Human/ophid.fasta '
       + 'data/OPHID.fasta.hhp.ps_scan '
       + 'data/HCV.conserved.90')
    sh('python parse_prosite.py '
       + 'data/OPHID.fasta.hhp.ps_scan '
       + '> data/OPHID.ps_scan')

@task
def STRING_domains():
    """Find CDs for HCV ELMs for STRING"""
    sh('python scan_prosite.py '
       + 'data/elm2prosite '
       + '~/bioperry/Projects/Thesis/Data/FASTA/Human/STRING.fasta '
       + 'data/STRING.fasta.hhp.ps_scan '
       + 'data/HCV.conserved.90')
    sh('python parse_prosite.py '
       + 'data/STRING.fasta.hhp.ps_scan '
       + '> data/STRING.ps_scan')

@task
def HPRD_domains():
    """Find CDs for HCV ELMs for HPRD"""
    sh('python scan_prosite.py '
       + 'data/elm2prosite '
       + '~/bioperry/Projects/Thesis/Data/FASTA/Human/hprd_new.intr.fasta '
       + 'data/HPRD.fasta.hhp.ps_scan '
       + 'data/HCV.conserved.90')
    sh('python parse_prosite.py '
       + 'data/HPRD.fasta.hhp.ps_scan '
       + '> data/HPRD.ps_scan')

@task
def domains():
    """Merge PROSITE domains w/ myLists"""

    for net_name, trans, ps_scan, net, elms, bg in get_ls():
        sh('python mylist2annotation.py '
           + '../Thesis/Data/ProfileScan/MyLists/ '
           + trans + ' '
           + 'data/' + net_name + '.mylist')
        sh('cat data/'
           + net_name + '.mylist '
           + ps_scan 
           + '> data/'
           + net_name + '.ps_mylist')
        sh('cut -f 1 data/'
           + net_name + '.ps_mylist | sort -u > '
           + net_name + '.ps_mylist.ls')
        # restrict background to network genes
        # with an entrez translation
        sh('cut -f 1 '
           + net
           + ' | sort -u > 1')
        sh('cut -f 2 '
           + net
           + ' | sort -u > 2')
        sh('cat 1 2 | sort -u > 12')
        sh('cut -f 1 '
           + trans + ' | sort -u > '
           + net_name + '.ls_pre')
        sh('intersect.py '
           + net_name + '.ls_pre 12 > '
           + net_name + '.ls')
        sh('python convertVersion2EntrezGeneID2.py '
           + net_name + '.ls '
           + trans + ' > '
           + net_name + '.ls.entrez')
        # additional restriction on background
        # to genes with PROSITE matches
        sh('intersect.py '
           + net_name + '.ls '
           + net_name + '.ps_mylist.ls '
           + '> '
           + net_name + '.bg')
        sh('python convertVersion2EntrezGeneID2.py '
           + net_name + '.bg '
           + trans
           + '> ' + net_name + '.bg.entrez')
        sh('rm ' + net_name + '.ls_pre')
    sh('rm 1 2 12')

@task
def HCV_hhp():
    """Predict human interactors"""

    for net_name, trans, ps_scan, net, elms, bg in get_ls():
        sh('python predict_cd_elm.py ' 
           + 'data/HCV.conserved.90 ' 
           + 'ELM '
           + elms + ' '
           + 'ELM '
           + ps_scan + ' '
           + 'ProfileScan '
           + net + ' ' 
           + net_name + '.ls ' 
           + trans + ' ' 
           + 'data/elm2prosite ' 
           + 'results/' + net_name + '.hcv_hhp ' 
           + 'results/' + net_name + '.hcv_hhp.vp2h12h2.tab')
 
@task
def HCV_pr_rec():
    """Recall & precision for HHP"""

    for net_name, trans, ps_scan, net, elms, bg in get_ls():
        sh('python pr_for_elm_predictions.py '
           + '../Thesis/Data/Network/HCV/hcv.hhe '
           + 'results/' + net_name + '.hcv_hhp '
           + net_name + '.ls.entrez '
           + 'results/' + net_name + '.tab')
        
