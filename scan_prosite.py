### This does not work b/c I can't figure out the 
### directory stuff

import os, utils_humanVirus ,sys

elm2cd = utils_humanVirus.get_elm2prosites(sys.argv[1])
elms = {}
with open(sys.argv[4]) as f:
    for line in f:
        elm = line.strip().split('\t')[3]
        elms[elm] = True

prosites = {}
for elm in elms:
    if elm in elm2cd:
        for prosite in elm2cd[elm]:
            if prosite.find('.CD') == -1:
                prosites[prosite] = True

cat_line = ''
for p in prosites:
    os.system('perl /bin/ps_scan/ps_scan.pl -d /bin/ps_scan/prosite.dat -p '
              + p + ' ' + sys.argv[2] + ' -o pff > '
              + p + '.ps_scan')
    cat_line = cat_line + p + '.ps_scan '
os.system('cat ' + cat_line
          + '> ' + sys.argv[3])
os.system('rm ' + cat_line)
