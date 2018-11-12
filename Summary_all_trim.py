#!/home/hkimlab/anaconda2/bin/python2.7

import os, sys
import subprocess as sp
from pdb import set_trace

sOutput_dir = sys.argv[1]
sProject    = sys.argv[2]
lRef_alt    = sys.argv[3].split(',')


def Concat_summary():

    sRef         = lRef_alt[0]
    sAlt         = lRef_alt[1]
    sSummary_dir = "{outdir}/Summary/Target".format(outdir=sOutput_dir)
    lHeader      = []
    lData        = []

    for sFile in os.listdir(sSummary_dir):
        if sRef + 'to' + sAlt in sFile:

            with open(sSummary_dir + '/' + sFile) as Input:
                for i, sRow in enumerate(Input):
                    if i == 0:
                        lCol = sRow.replace('\n', '').split('\t')
                        if lHeader:
                            for iCol_num in range(len(lHeader)):
                                if iCol_num > 6:

                                    if lHeader[iCol_num] == "" or lHeader[iCol_num] == " ":
                                        lHeader[iCol_num] = lCol[iCol_num]
                        else:
                            lHeader = lCol
                    else:
                        lData.append(sRow)
                #END: for
            #END: with
        #END: if
    #END: for

    if not os.path.isdir('{outdir}/Summary/Merge_target_result'.format(outdir=sOutput_dir)): os.mkdir('{outdir}/Summary/Merge_target_result'.format(outdir=sOutput_dir))

    with open('{outdir}/Summary/Merge_target_result/{project}_{ref}to{alt}_Summary.txt'.format(outdir=sOutput_dir,
                                                                                               project=sProject,
                                                                                               ref=sRef,
                                                                                               alt=sAlt), 'w') as Output:

        Output.write('\t'.join(lHeader) +'\n')
        for sData in lData:
            Output.write(sData)


if __name__ == '__main__':
    Concat_summary()
