#!/usr/bin/env python

import os, sys
import multiprocessing as mp
import subprocess as sp
from datetime import datetime
import pdb
from optparse import OptionParser


class Path_info(object):

    def __init__(self, sProject, options):

        self.sProject         = sProject
        self.iCore            = options.multicore
        self.sGap_open        = options.gap_open

        self.sGap_extend      = options.gap_extend
        self.sEnd_open        = options.end_open
        self.sEnd_extend      = options.end_extend
        self.sTarget_window   = options.target_window
        self.sIndel_check_pos = options.indel_check_pos
        self.sTarget_ref_alt  = options.target_ref_alt

        self.sBarcode       = './Input/Reference/%s/Barcode.txt' % self.sProject
        self.sReference_seq = './Input/Reference/%s/Reference.txt' % self.sProject

        self.sRef_path   = './Input/Reference/%s/Reference.fa' % self.sProject  # reference e.g.Consensus.fa
        self.sOutput_dir = './Output/%s' % self.sProject

        self.sPAM_seq = options.PAM_seq
        self.sPAM_pos = options.PAM_pos
        self.sGuide_pos = options.Guide_pos
        #if os.path.isdir(self.sOutput_dir):
        #    print("Result file already exists. Please remove the output file and try again.\ncd Output\nmv Project_name tmp")
        #    sys.exit()

        if not os.path.isdir(self.sOutput_dir): os.mkdir(self.sOutput_dir)
        if not os.path.isdir(self.sOutput_dir + '/Summary'): os.mkdir(self.sOutput_dir + '/Summary')
        if not os.path.isdir(self.sOutput_dir + '/result'): os.mkdir(self.sOutput_dir + '/result')
        if not os.path.isdir(self.sOutput_dir + '/result/freq'): os.mkdir(self.sOutput_dir + '/result/freq')
        if not os.path.isdir(self.sOutput_dir + '/result/freq/freq_result'): os.mkdir(self.sOutput_dir + '/result/freq/freq_result')


class Single_node_controller(Path_info):

    def __init__(self, sProject, options):
        super(Single_node_controller, self).__init__(sProject, options)


    def Make_reference(self):

        with open(self.sBarcode) as Barcode, \
            open(self.sReference_seq) as Ref, \
            open(self.sRef_path, 'w') as Output:

            dBarcode = {}

            for sBarcode in Barcode:
                lBarcode = sBarcode.replace('\n','').replace('\r','').split(':')
                sBar_sample = lBarcode[0]
                sBarcode = lBarcode[1]
                dBarcode[sBar_sample] = sBarcode

            for sRef in Ref:
                lRef        = sRef.replace('\n','').replace('\r','').split(':')
                sRef_sample = lRef[0]
                sRef        = lRef[1]

                try:
                    sBarcode = dBarcode[sRef_sample]
                    Output.write(sRef_sample + '\t' +sBarcode + '\t' + sRef + '\n')
                except KeyError:
                    print('no matching')
                    print(sRef_sample,sRef)


    def Make_indel_searcher_CMD(self):

        lCmd = []

        with open(self.sRef_path) as Barcode_ref:

            for sBarcode_ref in Barcode_ref:
                lBarcode_ref = sBarcode_ref.replace('\n', '').replace('\r','').split('\t')
                sSample  = lBarcode_ref[0]
                sBarcode = lBarcode_ref[1]
                sRef     = lBarcode_ref[2]

                sForward_Query = './Input/Query/%s/%s.txt' % (self.sProject, sSample)
                sCmd = './BaseEdit_freq_ver1.0.py {forw} {GapO} {GapE} {EndO} {EndE} {barcode} {ref} {target_window} {indel_check_pos}'
                sCmd += ' {target_ref_alt} {outdir} {sample} {PAM_seq} {PAM_pos} {Guide_pos}'
                lCmd.append(sCmd.format(forw=sForward_Query, GapO=self.sGap_open, GapE=self.sGap_extend, EndO=self.sEnd_open,
                                        EndE=self.sEnd_extend, barcode=sBarcode, ref=sRef, target_window=self.sTarget_window,
                                        indel_check_pos=self.sIndel_check_pos, target_ref_alt=self.sTarget_ref_alt, outdir=self.sOutput_dir,
                                        sample=sSample, PAM_seq=self.sPAM_seq, PAM_pos=self.sPAM_pos, Guide_pos=self.sGuide_pos))

        return lCmd

    def Run_indel_freq_calculator(self):
        sp.call('./Indel_frequency_calculator.py {outdir}'.format(outdir=self.sOutput_dir), shell=True)
        sp.call('./Summary_all_trim.py {outdir}'.format(outdir=self.sOutput_dir), shell=True)


def Run_indel_searcher(sCmd):
    sp.call(sCmd, shell=True)


def Run_multicore(lCmd, iCore):
    for sCmd in lCmd:
        print(sCmd)

    p = mp.Pool(iCore)
    p.map_async(Run_indel_searcher, lCmd).get()
    p.close()


def Main():
    print 'BaseEdit program start: %s' % datetime.now()

    sCmd = "BaseEdit frequency analyzer\n./Run_BaseEdit_freq_ver1.0.py -t 15 -w 16-48 --indel_check_pos 39-40 --target_ref_alt A,T --PAM_seq NGG --PAM_pos 43-45 --Guide_pos 23-42"
    sCmd += " --gap_open 20 --gap_extend 1 --end_open 20 --end_extend 1\n\n"
    sCmd += "1: Barcode\n"
    sCmd += "2: Base target window\n"
    sCmd += "3: Indel check pos\n"
    sCmd += "4: PAM pos\n"
    sCmd += "5: Guide pos (without PAM)\n"
    sCmd += "TATCTCTATCAGCACACAAGCATGCAATCACCTTGGGTCCAAAGGTCC\n"
    sCmd += "<------1------><----------------2--------------->\n"
    sCmd += "                                     <3>  <4>   \n"
    sCmd += "                      <---------5-------->      \n\n"

    parser = OptionParser(sCmd)

    parser.add_option("-t", "--thread", default="1", type="int", dest="multicore", help="multiprocessing number")
    parser.add_option("--gap_open", default="20", type="int", dest="gap_open", help="gap open: 1~100")
    parser.add_option("--gap_extend", default="1", type="float", dest="gap_extend", help="gap extend: 0.0005~10.0")
    parser.add_option("--end_open", default="20", type="int", dest="end_open", help="end open: 1~100")
    parser.add_option("--end_extend", default="1", type="float", dest="end_extend", help="end extend: 0.0005~10.0")
    parser.add_option("-w", "--target_window", type="str", dest="target_window", help="a window size for target sequence : 20-48")
    parser.add_option("--indel_check_pos", type="str", dest="indel_check_pos", help="indel check position to filter : 39-40; insertion 39, deletion 39 & 40")
    parser.add_option("--target_ref_alt", type="str", dest="target_ref_alt", help="Ref 'A' is changed to Alt 'T': A,T")
    parser.add_option("--PAM_seq", type="str", dest="PAM_seq", help="PAM sequence: NGG, NGC ...")
    parser.add_option("--PAM_pos", type="str", dest="PAM_pos", help="PAM position range in the reference seqeunce : 43-45")
    parser.add_option("--Guide_pos", type="str", dest="Guide_pos", help="Guide position range in the reference seqeunce : 23-42")

    options, args = parser.parse_args()

    with open('Project_list.txt') as Project_list,\
        open('Rerun.log') as Rerun_log_check:

        lRerun_log_check = [sRow.replace('\n', '') for sRow in Rerun_log_check]

        for sRow in Project_list:
            if sRow[0] == '#': continue
            sProject = sRow.replace('\n', '').replace('\r', '').replace(' ', '')
            Ins_single = Single_node_controller(sProject, options)
            Ins_single.Make_reference()

            lCmd = Ins_single.Make_indel_searcher_CMD()
            ###print(lCmd[:5])
            Run_multicore(lCmd, options.multicore)

            sMerge_cmd = "./Summary_all_trim.py Output/{project} {project} {ref_alt}".format(project=sProject, ref_alt=options.target_ref_alt)
            print('BaseEdit merge command', sMerge_cmd)
            sp.call(sMerge_cmd, shell=True)
            
            if sProject in lRerun_log_check: continue ## avoid to run redundant.

            try:            
                sFile_check = len([sMerged_files_check for sMerged_files_check in os.listdir('./Output/%s' % sProject + '/Summary/Merge_target_result')])
            except OSError:
                print('Empty project : %s' % sProject)
                continue

            #print(sFile_check)
            if sFile_check == 0:
                print('No merge file: %s' % sProject)
                break
            
            iLeak_count = 0
            
            for sCmd in lCmd:
                lCol = sCmd.split()
                sOutput_dir  = lCol[11]
                sOutput_file = lCol[12]
                sFile_for_check = sOutput_dir + '/Summary/All/' + sOutput_file + '_Summary.txt'
                
                if not os.path.isfile(sFile_for_check):
                    iLeak_count = 1
                    print('Rerunning leaked file : %s, %s' % (sProject, sOutput_file))
                    sp.call(sCmd, shell=True)
            
            if iLeak_count == 1:
                ## Rewrite to include leaked files.
                sp.call(sMerge_cmd, shell=True)
                #Ins_single.Run_indel_freq_calculator()
            
            with open('Rerun.log', 'a') as Rerun_log:
                Rerun_log.write('%s\n' % sProject)

    print 'BaseEdit program end: %s' % datetime.now()


if __name__ == '__main__':
    Main()
