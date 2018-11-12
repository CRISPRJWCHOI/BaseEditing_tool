#!/usr/bin/env python

import re
from pdb import set_trace
import sys
import os
import subprocess as sp
from datetime import datetime
from collections import OrderedDict

import numpy as np
from Bio import AlignIO

"""
variable prefix

s: string
l: list
i: int
f: float
d: dictionary
"""

sManual = """
Usage:

python2.7 ./indel_search_ver1.0.py splitted_input_1.fq splitted_input_2.fq reference.fa

splitted_input_1.fq : forward
splitted_input_2.fq : reverse

Total FASTQ(fq) lines / 4 = remainder 0.
"""
# iQual_cutoff = 20

if len(sys.argv) > 1:
    sForw_path       = sys.argv[1]
    sOG              = sys.argv[2]
    sOE              = sys.argv[3]
    sEG              = sys.argv[4]
    sEE              = sys.argv[5]
    sBarcode         = sys.argv[6]
    sRef             = sys.argv[7]
    sRef             = sRef[sRef.index(sBarcode):]  ## 'ACTG'<barcode>ACGACACACGCAT, leftside bases are redundant.
    lTarget_window   = sys.argv[8].split('-')
    lIndel_check_pos = sys.argv[9].split('-')
    lTarget_ref_alt  = sys.argv[10].split(',')
    sOutput_dir      = sys.argv[11]
    sSample_name     = sys.argv[12]
    sPAM_seq         = sys.argv[13]
    lPAM_pos         = sys.argv[14].split('-')
    lGuide_pos       = sys.argv[15].split('-')

else:
    print sManual
    sys.exit()

# index name, constant variable.
gNum_of_total = 0
gNum_of_ins = 1
gNum_of_del = 2
gNum_of_com = 3
gTotal_FASTQ = 4
gIns_FASTQ = 5
gDel_FASTQ = 6
gCom_FASTQ = 7
gINDEL_info = 8


def Open_Sequence_files():
    lSequence_forward = []
    with open(sForw_path) as fa_1:
        lSequence_forward = [sRow.replace('\n', '').upper() for sRow in fa_1]
    return lSequence_forward


# end: def


class BaseEdit_parser():

    def Calculate_BaseEdit_freq(self, lQuery_seq=[]):

        dRef = {}
        dResult = {}

        dRef[sBarcode] = (sRef)  # total matched reads, insertion, deletion, complex
        dResult[sBarcode] = [0, 0, 0, 0, [], [], [], [], [], [], []]

        # lRef   : [(ref_seq, ref_seq_after_barcode, barcode, barcode end pos, indel end pos, indel from barcode),(...)]
        # dResult = [# of total, # of ins, # of del, # of com, [total FASTQ], [ins FASTQ], [del FASTQ], [com FASTQ], info]
        iCount = 0


        for sQuery_seq_raw in lQuery_seq:

            iBarcode_matched = 0
            iNeedle_matched = 0
            iInsert_count = 0
            iDelete_count = 0
            iComplex_count = 0

            try:
                # Check the barcode pos and remove it.
                sQuery_seq_raw = sQuery_seq_raw.replace('\r', '')
                iBarcode_start_pos = sQuery_seq_raw.index(sBarcode)
                iBarcode_matched += 1

                sQuery_seq_with_barcode = sQuery_seq_raw[iBarcode_start_pos:]

                # _check = 0
                # if sQuery_seq_raw == 'TCTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTATCTCTATCAGCACACAAGCATGCAATCACCTTGGGTCCCAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA':
                #     _check = 1

                sRef_seq = r'<(echo -e ">{name}\n{seq}")'.format(name=sBarcode + '_ref', seq=sRef)
                sQuery_seq = r'<(echo -e ">{name}\n{seq}")'.format(name=sBarcode + '_query', seq=sQuery_seq_with_barcode)

                sNeedle_cmd = r"/bin/bash -c 'needle -filter {0} {1} -outfile stdout -gapopen {2} -gapextend {3} -endweight Y -endopen {4} -endextend {5}'".format(
                    sRef_seq, sQuery_seq, sOG, sOE, sEG, sEE)

                Needle_result = sp.Popen(sNeedle_cmd, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True,
                                         shell=True)

                lResult           = [Instance.seq._data for Instance in AlignIO.read(Needle_result.stdout, "emboss")]
                sRef_needle_ori   = lResult[0]
                sQuery_needle_ori = lResult[1]

                # if _check == 1:
                #     print(sRef_needle_ori)
                #     print(sQuery_needle_ori)
                #     set_trace()

                Needle_result.stdout.close()

                # detach forward ---, backward ---
                # e.g.    ref   ------AAAGGCTACGATCTGCG------
                #         query AAAAAAAAATCGCTCTCGCTCTCCGATCT
                # trimmed ref         AAAGGCTACGATCTGCG
                # trimmed qeury       AAATCGCTCTCGCTCTC
                iReal_ref_needle_start = 0
                iReal_ref_needle_end   = len(sRef_needle_ori)
                iRef_needle_len        = len(sRef_needle_ori)

                for i, sRef_nucle in enumerate(sRef_needle_ori):
                    if sRef_nucle in ['A', 'C', 'G', 'T']:
                        iReal_ref_needle_start = i
                        break

                for i, sRef_nucle in enumerate(sRef_needle_ori[::-1]):
                    if sRef_nucle in ['A', 'C', 'G', 'T']:
                        iReal_ref_needle_end = iRef_needle_len - (i + 1)
                        # forward 0 1 2  len : 3
                        # reverse 2 1 0,  len - (2 + 1) = 0
                        break

                sRef_needle = sRef_needle_ori[iReal_ref_needle_start:iReal_ref_needle_end + 1]
                if iReal_ref_needle_start:
                    sQuery_needle = sQuery_needle_ori[:iReal_ref_needle_end]
                sQuery_needle = sQuery_needle_ori[:len(sRef_needle)]
                # detaching completion

                # indel info making.
                iNeedle_match_pos_ref   = 0
                iNeedle_match_pos_query = 0
                iNeedle_insertion       = 0
                iNeedle_deletion        = 0

                lInsertion_in_read = []  # insertion result [[100, 1], [119, 13]]
                lDeletion_in_read  = []  # deletion result  [[97, 1], [102, 3]]

                # print 'sRef_needle', sRef_needle
                # print 'sQuery_needle', sQuery_needle
                for i, (sRef_nucle, sQuery_nucle) in enumerate(zip(sRef_needle, sQuery_needle)):

                    if sRef_nucle == '-':
                        iNeedle_insertion += 1

                    if sQuery_nucle == '-':
                        iNeedle_deletion += 1

                    if sRef_nucle in ['A', 'C', 'G', 'T']:
                        if iNeedle_insertion:
                            lInsertion_in_read.append([iNeedle_match_pos_ref, iNeedle_insertion])
                            iNeedle_insertion = 0
                        iNeedle_match_pos_ref += 1

                    if sQuery_nucle in ['A', 'C', 'G', 'T']:
                        if iNeedle_deletion:
                            lDeletion_in_read.append([iNeedle_match_pos_query, iNeedle_deletion])
                            iNeedle_match_pos_query += iNeedle_deletion
                            iNeedle_deletion = 0
                        iNeedle_match_pos_query += 1
                        # print 'sRef_needle', sRef_needle

                # print 'sQuery_needle', sQuery_needle
                # print 'lInsertion_in_read: onebase', lInsertion_in_read
                # print 'lDeletion_in_read: onebase', lDeletion_in_read
                # print 'i5bp_front_Indel_end', i5bp_front_Indel_end
                # print 'iIndel_end_from_barcode_pos', iIndel_end_from_barcode_pos

                lTarget_indel_result = []  # ['20M2I', '23M3D' ...]

                """
                ins case
                ...............................NNNNNNNNNNNNNN....NNNNNNNNNNNNNNNNNNN*NNNNNAGCTT
                """

                iCleavage_window_start = int(lIndel_check_pos[0])
                iCleavage_window_end = int(lIndel_check_pos[1]) - 1

                for iMatch_pos, iInsertion_pos in lInsertion_in_read:
                    if iCleavage_window_start <= iMatch_pos <= iCleavage_window_end:  # iMatch_pos is one base
                        iInsert_count = 1
                        lTarget_indel_result.append(str(iMatch_pos) + 'M' + str(iInsertion_pos) + 'I')
                """
                del case 1
                ...............................NNNNNNNNNNNNNN....NNNNNNNNNNNNNNNNNNNNN**NNNAGCTT
                del case 2
                ...............................NNNNNNNNNNNNNN....NNNNNNNNNNNNNNNNNNNNN**NNNNNCTT
                """
                for iMatch_pos, iDeletion_pos in lDeletion_in_read:

                    """
                    Insertion: 30M3I
                           ^
                    ACGT---ACGT
                    ACGTTTTACGT -> check this seq
                    Insertion just check two position

                    Deletion: 30M3D
                         ^
                    ACGTTTTACGT
                    ACGT---ACGT -> check this seq
                    But deletion has to includes overlap deletion.
                    """

                    if iMatch_pos <= iCleavage_window_end and iCleavage_window_start <= (iMatch_pos + iDeletion_pos):
                        iDelete_count = 1
                        lTarget_indel_result.append(str(iMatch_pos) + 'M' + str(iDeletion_pos) + 'D')

                if iInsert_count == 1 and iDelete_count == 1:
                    iComplex_count = 1
                    iInsert_count = 0
                    iDelete_count = 0

                    # """ test set
                    # print 'sBarcode', sBarcode
                    # print 'sTarget_region', sTarget_region
                    # print 'sRef_seq_after_barcode', sRef_seq_after_barcode
                    # print 'sSeq_after_barcode', sQuery_seq
                    # print 'iIndel_start_from_barcode_pos', iIndel_start_from_barcode_pos
                    # print 'iIndel_end_from_barcode_pos', iIndel_end_from_barcode_pos
                    # """

                    """
                    23M3I
                    23M is included junk_seq after barcode,

                    barcorde  junk   targetseq   others
                    *********ACCCT-------------ACACACACC
                    so should select target region.
                    If junk seq is removed by target region seq index pos.
                    """

                ## 8: indel info
                dResult[sBarcode][8].append(
                    [sRef, sQuery_seq_raw, lTarget_indel_result,
                     "", sRef_needle_ori, sQuery_needle_ori])  ## "" -> target seq, but this is not used this project.



            # end: try
            except ValueError as e:
                # print(e)
                continue

            # total matched reads, insertion, deletion, complex
            dResult[sBarcode][0] += iBarcode_matched
            dResult[sBarcode][1] += iInsert_count
            dResult[sBarcode][2] += iDelete_count
            dResult[sBarcode][3] += iComplex_count

            ## base editing frequency
            """
                   BaseEditPos : 0                                                    1                                  2
            [OrderedDict([('A',0),('C',0),('G',0),('T',0)]), OrderedDict([('A',0),('C',0),('G',0),('T',0)]), ...

            and sum the counts each position
            """

            if iInsert_count == 0 and iDelete_count == 0 and iComplex_count == 0:

                lBaseEdit = []
                iTarget_len = int(lTarget_window[1]) - int(lTarget_window[0]) + 1

                for i in range(iTarget_len):
                    lBaseEdit.append(OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]))

                iTarget_start = int(lTarget_window[0]) - 1
                iTarget_end = int(lTarget_window[1])

                """
                                       cleavage window start
                                        ^
                [barcode]ACGACGTACGACGT[cleavage]
                [barcode]ACGACGTACGACGT[cleavage]
                """

                iBase_edit_event = 0

                for i, tRef_Query_base in enumerate(zip(sRef_needle[iTarget_start: iTarget_end], sQuery_needle[iTarget_start: iTarget_end])):
                    sRef_base = tRef_Query_base[0]
                    sQuery_base = tRef_Query_base[1]

                    if sRef_base == '-' or sQuery_base == '-': continue

                    if sRef_base != sQuery_base and sQuery_base != 'N':
                        iBase_edit_event = 1
                        lBaseEdit[i][sQuery_base] += 1
                        # print(sQuery_needle)

                dResult[sBarcode][9].append(lBaseEdit)
                if iBase_edit_event == 1:
                    dResult[sBarcode][10].append([sRef, sQuery_seq_raw, lTarget_indel_result, [list(orderedDict.values()) for orderedDict in lBaseEdit], sRef_needle_ori, sQuery_needle_ori])
                # dResult[sBarcode] = [0, 0, 0, 0, [], [], [], [], [], [BaseEdit_freq_data]]

            iBarcode_matched = 0
            iInsert_count = 0
            iDelete_count = 0
            iComplex_count = 0
            # end: for sBarcode, lCol_ref
        # end: for lCol_FASTQ
        return dResult


def sClassify_file(sBarcode, sJudge): return ('{outdir}/result/{barcode}.{judge}.fastq'.format(outdir=sOutput_dir, barcode=sBarcode, judge=sJudge))


def Write_FASTQ_output(iFASTQ_kind, lValue, FASTQ_out):
    for lFASTQ in lValue[iFASTQ_kind]:
        FASTQ_out.write('\n'.join(lFASTQ) + '\n')


def Make_output(dResult):
    """
   {'TTTGGTGCACACACATATA': [6, 2, 2, 0, [], [], [], [], [['TATCTCTA..ref', 'GAGTCGGTG...query', [13M5D], '',
   'TTTGGTGCACACACATATAACTGGAACACAAAGCATAGACTGCGGGGCG------------------------------------------------------------',
   'TTTGGTGCACACACATATAACTGGAACACAAAGCATAGA-TGCGGGGCGTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA'],
   ['TTTGGTGCACACACATATAACTGGAACACAAAGCATAGACTGCGGGGCG', '', '', '',
   'TTTGGTGCACACACATATAACTGGAACACAAAGCATAGACTGCGGGGCG------------------------------------------------------------', ...
   [[OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 1)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)])],
   [OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 1)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)])]]]}
    """

    #set_trace()


    with open('{outdir}/result/{sample_name}_filtered_indel.txt'.format(outdir=sOutput_dir, sample_name=sSample_name), 'w') as Filtered,\
        open('{outdir}/result/{sample_name}_aligned_BaseEdit.txt'.format(outdir=sOutput_dir, sample_name=sSample_name), 'w') as Ref_Alt_edit:

        for sBarcode in dResult:
            for lAligned_indel_result in dResult[sBarcode][8]:  # 8 : indel list
                if lAligned_indel_result[2]:
                    Filtered.write('\t'.join(map(str, lAligned_indel_result)) + '\n')

            for lAligned_alt_result in dResult[sBarcode][10]:  # 10 : alt base list
                if lAligned_alt_result:
                    lAligned_alt_result[2] = str(lAligned_alt_result[2])
                    try:
                        Ref_Alt_edit.write('\t'.join(map(str, lAligned_alt_result)) + '\n')
                    except Exception:
                        set_trace()

        """
        lAligned_result
        ['TATCTCTATCAGCACACAAGCATGCAATCACCTTGGGTCCAAAGGTCC', 'TCTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTATCTCTATCAGCACACAAGCATGCAATCACCTTGGGTCAAAGGTCCAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAAT\r',
        ['38M1D'], '', 'TATCTCTATCAGCACACAAGCATGCAATCACCTTGGGTCCAAAGGTCC-----------------------------------------------------------------', 'TATCTCTATCAGCACACAAGCATGCAATCACCTTGGGT-CAAAGGTCCAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAAT']
        """

    dSelect_base = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    sTarget_ref = lTarget_ref_alt[0]
    sTarget_alt = lTarget_ref_alt[1]

    iTarget_base = dSelect_base[sTarget_alt]
    
    try:
        if not os.path.isdir('{outdir}/Summary/All'.format(outdir=sOutput_dir)):
            os.mkdir('{outdir}/Summary/All'.format(outdir=sOutput_dir))
        if not os.path.isdir('{outdir}/Summary/Target'.format(outdir=sOutput_dir)):
            os.mkdir('{outdir}/Summary/Target'.format(outdir=sOutput_dir))
    except OSError:
        pass

    for sBarcode, lValue in dResult.items():

        iBarcode_start_pos       = sRef.index(sBarcode)
        sRef_seq_without_barcode = sRef[iBarcode_start_pos+len(sBarcode):]

        llBaseEdit = lValue[9]
        lSum = []

        for i, lBaseEdit in enumerate(llBaseEdit):

            if not lSum:
                lSum = [[0, 0, 0, 0] for iQuery in range(len(lBaseEdit))]

            for j in range(len(lBaseEdit)):
                for k, iCount in enumerate(list(llBaseEdit[i][j].values())):
                    lSum[j][k] += iCount

        with open('{outdir}/Summary/All/{sample_name}_Summary.txt'.format(outdir=sOutput_dir, sample_name=sSample_name), 'w') as Summary, \
            open('{outdir}/Summary/Target/{sample_name}_{target}_Summary.txt'.format(outdir=sOutput_dir, sample_name=sSample_name, target=sTarget_ref + 'to' + sTarget_alt), 'w') as Target_summary:

            ## This Ref has barcode.
            sRef_target = sRef[int(lTarget_window[0]) - 1:int(lTarget_window[1])]

            iPAM_start    = int(lPAM_pos[0]) - 1
            iPAM_end      = int(lPAM_pos[1])
            iGuide_start  = int(lGuide_pos[0]) - 1
            iGuide_end    = int(lGuide_pos[1])
            iGuide_len    = iGuide_end - iGuide_start
            iBarcode_len  = len(sBarcode)

            """
            barcode Guide st,ed 
            <----><----------> NGG
            ACGTACGTACGTACGTACGTGGACG
            """
            #sRef_target[iPAM_start:iPAM_end] = sPAM_seq
            iWithout_target_len = len(sRef_target[iBarcode_len:iGuide_start])
            lWithout_target_pos = [-(i+1) for i in range(iWithout_target_len)][::-1]

            lWith_target_pos = [i + 1 for i in range(iGuide_len)]
            lAfter_PAM_pos   = [i + 1 for i in range(len(sRef) - iPAM_end + 1)]

            lPos_num = lWithout_target_pos + lWith_target_pos + list(sPAM_seq) + lAfter_PAM_pos
            lPos_annotated_ref = [str(i)+'.'+str(j) for i,j in zip(sRef_target, lPos_num)]
            ## ['A.-7', 'C.-6', 'A.-5', 'A.-4', 'G.-3', 'C.-2', 'A.-1', 'T.1', 'G.2', 'C.3', 'A.4', 'A.5', 'T.6', 'C.7', 'A.8', 'C.9', 'C.10', 'T.11', 'T.12', 'G.13', 'G.14',

            lMasked_pos_annotated_ref_target = []   ## '' '' '' A '' '' '' A A '' ''

            for sBase_pos in lPos_annotated_ref:
                sBase_only = sBase_pos.split('.')[0]
                if sBase_only != sTarget_ref:
                    lMasked_pos_annotated_ref_target.append(' ')
                else:
                    lMasked_pos_annotated_ref_target.append(sBase_pos)

            #set_trace()

            ## Making a header
            Summary.write("Sample\tBarcode\tRef\t# of Total\t# of Insertion\t# of Deletion\t# of Combination\t{refseq}\n".format(refseq='\t'.join(lPos_annotated_ref)))
            Target_summary.write("Sample\tBarcode\tRef\t# of Total\t# of Insertion\t# of Deletion\t# of Combination\t{refseq}\n".format(refseq='\t'.join(lMasked_pos_annotated_ref_target)))

            for i, lBase_count in enumerate(zip(*lSum)):  ## lBase_count [(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0)]

                if i == 0:
                    Summary.write("{sample}\t{bar}\t{ref}\t{NumTot}\t{NumIns}\t{NumDel}\t{NumCom}\t{BaseEditCount}\n".format(
                        sample=sSample_name, bar=sBarcode, ref=sRef_seq_without_barcode, NumTot=lValue[gNum_of_total], NumIns=lValue[gNum_of_ins], NumDel=lValue[gNum_of_del], NumCom=lValue[gNum_of_com],
                        BaseEditCount='\t'.join(map(str, lBase_count))))
                else:
                    Summary.write("\t\t\t\t\t\t\t{BaseEditCount}\n".format(BaseEditCount='\t'.join(map(str, lBase_count))))

            try:
                lTarget_base_count = zip(*lSum)[iTarget_base]
                lMasked_target_base_count = []  ## '' 20 '' 30 '' '' '' '' 20 ''

                for sMasked_ref, fCount in zip(lMasked_pos_annotated_ref_target, lTarget_base_count):

                    if sMasked_ref == ' ':
                        lMasked_target_base_count.append(' ')
                    else:
                        lMasked_target_base_count.append(fCount)

                Target_summary.write(("{sample}\t{bar}\t{ref}\t{NumTot}\t{NumIns}\t{NumDel}\t{NumCom}\t{BaseEditCount}\n".format(
                    sample=sSample_name, bar=sBarcode, ref=sRef_seq_without_barcode, NumTot=lValue[gNum_of_total], NumIns=lValue[gNum_of_ins], NumDel=lValue[gNum_of_del], NumCom=lValue[gNum_of_com],
                    BaseEditCount='\t'.join(map(str, lMasked_target_base_count)))))

            except IndexError:
                print('Null query: ', sForw_path)
                ## Null query base count is all zero.
                Target_summary.write(
                    ("{sample}\t{bar}\t{ref}\t{NumTot}\t{NumIns}\t{NumDel}\t{NumCom}\t{BaseEditCount}\n".format(
                        sample=sSample_name, bar=sBarcode, ref=sRef_seq_without_barcode, NumTot=lValue[gNum_of_total],
                        NumIns=lValue[gNum_of_ins], NumDel=lValue[gNum_of_del], NumCom=lValue[gNum_of_com],
                        BaseEditCount='\t'.join(['0'] * len(lPos_annotated_ref)))))


def Main():
    # Output: 1. Count information of matched barcode e.g. TACGATCTA\t# total\tins\t# del\t# com
    # Output: 2. classify FASTQ.    e.g. TAGAATATACACG.insertion.fastq
    try:
        if not os.path.isdir('{outdir}'.format(outdir=os.path.dirname(sOutput_dir))): os.mkdir('{outdir}'.format(outdir=os.path.dirname(sOutput_dir)))
        if not os.path.isdir('{outdir}'.format(outdir=sOutput_dir)): os.mkdir('{outdir}'.format(outdir=sOutput_dir))
        if not os.path.isdir('{outdir}/Summary/'.format(outdir=sOutput_dir)): os.mkdir('{outdir}/Summary/'.format(outdir=sOutput_dir))
        if not os.path.isdir('{outdir}/result/'.format(outdir=sOutput_dir)): os.mkdir('{outdir}/result/'.format(outdir=sOutput_dir))
        if not os.path.isdir('{outdir}/result/freq'.format(outdir=sOutput_dir)): os.mkdir('{outdir}/result/freq'.format(outdir=sOutput_dir))
        if not os.path.isdir('{outdir}/result/freq/freq_result'.format(outdir=sOutput_dir)): os.mkdir('{outdir}/result/freq/freq_result'.format(outdir=sOutput_dir))
    except OSError:
        pass

    print 'Program start : ' + sSample_name
    sStart_time = datetime.now()
    print sStart_time

    print 'File Open : ' + sSample_name
    lSequence_forward = Open_Sequence_files()

    Ins_BaesEdit = BaseEdit_parser()
    dResult_forward = Ins_BaesEdit.Calculate_BaseEdit_freq(lSequence_forward)

    print 'Calculate INDEL frequency : ' + sSample_name
    print 'Make output forward : ' + sSample_name
    Make_output(dResult_forward)

    sEnd_time = datetime.now()
    print sEnd_time
    print 'Program end : ' + sSample_name
# end: def Main


if __name__ == '__main__':
    Main()


