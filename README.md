# BaseEditing frequency measurement tool

### Prerequisite ###
    OS Requirements: Linux(Ubuntu 16.04)
    Install time: about 30 minutes
    
    - python2.7
    pip install numpy==1.12.1 pandas==0.20.1
    or
    - anaconda python2.7
    conda install numpy=1.12.1 pandas=0.20.1
    
    EMBOSS:6.6.0.0, needle
    http://emboss.sourceforge.net/download/
    

### Usage ###

    ./Project_list.txt
    Write the one project name per line
    
    /Input/Query/<project_name>/<sample_name>.txt
    /Input/Reference/<project_name>/Barcode.txt
    /Input/Reference/<project_name>/Reference.txt
    
    Test
    ./nohup.sh
    
    Detailed options
    python2.7 Run_BaseEdit_freq_ver1.0.py -h
    
    Output
    /Output/<project_name>/Summary_result/All
    Sample	Barcode	Ref	# of Total	# of Insertion	# of Deletion	# of Combination	C.-7	A.-6	T.-5	A.-4	G.-3	G.-2	A.-1	T.1	G.2	C.3	C.4	C.5	C.6	A.7	T.8	C.9	A.10	C.11	T.12	C.13	A.14	C.15	C.16	G.17	G.18	G.19	G.20	A.N	G.G	G.G	C.1	T.2	G.3
    Doench2014_900	TCACTCTCTGCTGCA CATAGGATGCCCCATCACTCACCGGGGAGGCTGAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA	9	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0

    /Output/<project_name>/Summary_result/Target
    Sample	Barcode	Ref	# of Total	# of Insertion	# of Deletion	# of Combination	 	A.-6	 	A.-4	 	 	A.-1	 	 	 	 	 	 	A.7	 	 	A.10	 	 	 	A.14	 	 	 	 	 	 	A.N	 	 	 	 	 
    Doench2014_900	TCACTCTCTGCTGCA	CATAGGATGCCCCATCACTCACCGGGGAGGCTGAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA	9	0	0	0	 	0	 	0	 	 	0	 	 	 	 	 	 	0	 	 	0	 	 	 	0	 	 	 	 	 	 	0	 	 	 	 	 

    /Output/<project_name>/Summary_result/Merge_target_result
    Sample	Barcode	Ref	# of Total	# of Insertion	# of Deletion	# of Combination	A.-7	A.-6	A.-5	A.-4	A.-3	A.-2	A.-1	A.1	A.2	A.3	A.4	A.5	A.6	A.7	A.8	A.9	A.10	A.11	A.12	A.13	A.14	A.15	A.16	A.17	A.18	A.19	A.20	A.N	 	 	A.1	A.2	A.3
    Doench2014_900	TCACTCTCTGCTGCA	CATAGGATGCCCCATCACTCACCGGGGAGGCTGAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA	9	0	0	0	 	0	 	0	 	 	0	 	 	 	 	 	 	0	 	 	0	 	 	 	0	 	 	 	 	 	 	0	 	 	 	 	 

    An additional base editing result
    vi Additional_BaseEdit_process_list.tsv
    column1: project name
    column2: other base editing
    column3: existing result file (/Output/<project_name>/Summary_result/Target)
    
    180903_split_hiseq_10_1_1   A,G 180903_split_hiseq_10_1_1_AtoT_Summary.txt
    180903_split_hiseq_10_1_2   A,C 180903_split_hiseq_10_1_2_AtoT_Summary.txt

    And run python program
    python2.7 Run_each_base_summary.py
