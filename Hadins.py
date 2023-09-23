# coding=utf-8
import sys
import os






def usage():
    sys.exit(
        '====================Hadins====================\n'
        'generate_feature:\n'
        'python Hadins.py generate_feature lbamfile sbamfile output process includecontig\n'
        'lbamfile is the path of the alignment file about the reference and the long read set\n'
        'sbamfile is the path of the alignment file about the reference and the short read set\n'
        'output is a folder which is used to store features\n'
        'process is the number of processes to use\n'
        'includecontig is the list of contig to preform detection.(default: [], all contig are used) \n'
        '----------------------------------------------------\n'
        'call_ins:\n'
        'python Hadins.py weight datapath lbamfile outvcfpath support includecontig\n'
        'weight is the path of the model weight\n'
        'datapath is a folder which is used to store features\n'
        'lbamfile is the path os the alignment file about the reference and the long read set\n'
        'outvcfpath is the path of vcf file\n'
        'support is min support reads\n'
        'includecontig is the list of contig to preform detection.(default: [], all contig are used)\n '
        '----------------------------------------------------\n'
        'eg:python Hadins.py generate_feature long_read.bam short_read.bam outpath 5 [12,13,14,15,16,17,18,19,20,21,'
        '22]\n'
        'eg:python Hadins.py call_ins ins.h5 datapath long_read.bam outvcfpath 10 [12,13,14,15,16,17,18,19,20,21,22]\n'
    )


model = sys.argv[1]
flag = 0  # 0代表正常，1代表异常
lens = len(sys.argv)
print(sys.argv)
if model == "generate_feature":
    print('hahaha')
    # 开始执行提取特征
    if lens not in [5, 6, 7]:
        flag = 1
        print("噢噢噢噢哦哦哦")
    else:
        print("start generate feature")
        bamfilepath_long = ''
        bamfilepath_short = ''
        outputpath = ''
        max_work = 5
        includecontig = []
        if lens == 7:
            bamfilepath_long, bamfilepath_short, outputpath, max_work, includecontig = sys.argv[2], sys.argv[3], \
                                                                                       sys.argv[4], sys.argv[5], [
                                                                                           str(contig) for contig in
                                                                                           eval(sys.argv[6])]
        elif lens == 6:
            bamfilepath_long, bamfilepath_short, outputpath, max_work, includecontig = sys.argv[2], sys.argv[3], \
                                                                                       sys.argv[4], sys.argv[5], []
        elif lens == 5:
            bamfilepath_long, bamfilepath_short, outputpath, max_work, includecontig = sys.argv[2], sys.argv[3], \
                                                                                       sys.argv[4], 10, []

        from exact_ls_feature import creat_data_longshort_by_pool
        creat_data_longshort_by_pool(bamfilepath_long, bamfilepath_short, outputpath, includecontig, int(max_work))

elif model == 'call_ins':
    if lens != 8:
        flag = 1
    else:
        ins_predict_weight, datapath, bamfilepath, outvcfpath, support, contig = \
            sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], [str(contig) for contig in
                                                                              eval(sys.argv[8])]
        import tensorflow as tf
        from model_predict import predict_funtion
        predict_funtion(ins_predict_weight, datapath, bamfilepath, outvcfpath, support, contig)
else:
    flag = 1

if flag == 1:
    usage()
