#!/home/hkimlab/anaconda2/bin/python2.7

import os


def Make_ref():

    with open(r"./Input/Reference/Barcode.txt") as Barcode,\
        open(r"./Input/Reference/Target_region.txt") as Target,\
        open(r"./Input/Reference/Reference_sequence.txt") as Ref,\
        open('./Input/Reference/Consensus.fa', 'w') as Output:

        lName = [sBar.replace('\n', '') + ':' + sTar for sBar, sTar in zip(Barcode, Target)]

        for i, sRow in enumerate(Ref):
            Output.write('>'+lName[i]+sRow)

Make_ref()
