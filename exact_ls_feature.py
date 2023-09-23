# coding=utf-8
import pysam
import numpy as np
import pandas as pd
import multiprocessing as mp
import math
import os

WINDOW = 200
LENGTHS = 10000000


def decode_flag(Flag):
    signal = {1 << 2: 0, 1 >> 1: 1, 1 << 4: 2, 1 << 11: 3, 1 << 4 | 1 << 11: 4}
    return signal[Flag] if (Flag in signal) else 0


def c_pos(cigar, refstart):
    number = ''
    numlist = [str(i) for i in range(10)]
    readstart = False
    readend = False
    refend = False
    readloc = 0
    refloc = refstart
    for c in cigar:
        if (c in numlist):
            number += c
        else:
            number = int(number)
            if readstart == False and c in ['M', 'I', '=', 'X']:
                readstart = readloc
            if readstart != False and c in ['H', 'S']:
                readend = readloc
                refend = refloc
                break

            if c in ['M', 'I', 'S', '=', 'X']:
                readloc += number

            if c in ['M', 'D', 'N', '=', 'X']:
                refloc += number
            number = ''
    if readend == False:
        readend = readloc
        refend = refloc

    return refstart, refend, readstart, readend


def calculate_splitread(chr_name, bamfile):
    dada = []
    for read in bamfile.fetch(chr_name):

        # print(read.has_tag('SA'))
        if read.has_tag('SA'):
            code = decode_flag(read.flag)
            rawsalist = read.get_tag('SA').split(';')
            for sa in rawsalist[:-1]:
                sainfo = sa.split(',')
                # print(sainfo)
                tmpcontig, tmprefstart, strand, cigar = sainfo[0], int(sainfo[1]), sainfo[2], sainfo[3]
                if tmpcontig != chr_name:
                    continue
                if (strand == '-' and (code % 2) == 0) or (strand == '+' and (code % 2) == 1):
                    refstart_1, refend_1, readstart_1, readend_1 = read.reference_start, read.reference_end, read.query_alignment_start, read.query_alignment_end

                    refstart_2, refend_2, readstart_2, readend_2 = c_pos(cigar, tmprefstart)
                    a = readend_1 - readstart_2
                    b = refend_1 - refstart_2
                    if abs(b - a) < 30:
                        continue
                    # print(b-a)
                    if abs(b) < 2000:
                        if (b - a) > 50 and abs(b - a) < 200000:
                            dada.extend([refstart_2])
                            dada.extend([refend_1])
    data = pd.value_counts(dada)

    return data


def feature_extraction_long(bamfile_long, chr_name_long, start, end, mapq):
    ref_pos = []
    ins_count = []
    loci_clip_sm = []
    loci_clip_ms = []
    for read in bamfile_long.fetch(chr_name_long, start, end):
        aligned_length = read.reference_length
        if aligned_length is None:
            aligned_length = 0
        if (read.mapping_quality > mapq) and (aligned_length > 0):
            cigar = np.array(read.cigartuples)
            ref_pos += (read.get_reference_positions())
            ref_pos_start = read.reference_start + 1
            cigar_shape = cigar.shape
            for i in range(cigar_shape[0]):
                if cigar[i, 0] == 0:
                    ref_pos_start = cigar[i, 1] + ref_pos_start
                elif cigar[i, 0] == 7:
                    ref_pos_start = cigar[i, 1] + ref_pos_start
                elif cigar[i, 0] == 8:
                    ref_pos_start = cigar[i, 1] + ref_pos_start
                elif cigar[i, 0] == 2:
                    ref_pos_start = cigar[i, 1] + ref_pos_start
                elif cigar[i, 0] == 1:
                    if cigar[i, 1] >= 20:
                        ins_count += [ref_pos_start]
            if cigar[0, 0] == 4:
                loci_clip_sm.append(read.reference_start + 1)
            if cigar[-1, 0] == 4:
                loci_clip_ms.append(read.reference_end)
    ref_pos = np.array(ref_pos) + 1
    if len(ref_pos) == 0:
        ref_pos = np.array([0])
    ref_pos = pd.value_counts(ref_pos)
    ins_count = pd.value_counts(ins_count)
    loci_clip_sm = np.array(loci_clip_sm)
    loci_clip_ms = np.array(loci_clip_ms)

    loci_clip_sm = pd.value_counts(loci_clip_sm)
    loci_clip_ms = pd.value_counts(loci_clip_ms)

    return ref_pos, ins_count, loci_clip_sm, loci_clip_ms


def feature_extraction_short(bamfile_short, chr_name_short, start, end, mapq):
    ref_pos = []
    loci_clip_sm = []
    loci_clip_ms = []
    for read in bamfile_short.fetch(chr_name_short, start, end):
        aligned_length = read.reference_length
        if aligned_length is None:
            aligned_length = 0
        if (read.mapping_quality > mapq) and aligned_length > 0:
            cigar = np.array(read.cigartuples)
            ref_pos += (read.get_reference_positions())
            if cigar[0, 0] == 4:
                loci_clip_sm.append(read.reference_start + 1)
            if cigar[-1, 0] == 4:
                loci_clip_ms.append(read.reference_end)

    ref_pos = np.array(ref_pos) + 1
    if len(ref_pos) == 0:
        ref_pos = np.array([0])
    ref_pos = pd.value_counts(ref_pos)
    loci_clip_sm = np.array(loci_clip_sm)
    loci_clip_ms = np.array(loci_clip_ms)
    loci_clip_sm = pd.value_counts(loci_clip_sm)
    loci_clip_ms = pd.value_counts(loci_clip_ms)

    return ref_pos, loci_clip_sm, loci_clip_ms


def compute(cover_long, loci_ins_long, loci_clip_long_sm, loci_clip_long_ms, split_read_long, cover_short,
            loci_clip_short_sm, loci_clip_short_ms, start, end):
    s_e = np.arange(start, end)
    cover_long = cover_long.reindex(index=s_e).fillna(value=0).values.reshape(-1, 1)
    loci_ins_long = loci_ins_long.reindex(index=s_e).fillna(value=0).values.reshape(-1, 1)
    loci_clip_long_sm = loci_clip_long_sm.reindex(index=s_e).fillna(value=0).values.reshape(-1, 1)
    loci_clip_long_ms = loci_clip_long_ms.reindex(index=s_e).fillna(value=0).values.reshape(-1, 1)
    split_read_long = split_read_long.reindex(index=s_e).fillna(value=0).values.reshape(-1, 1)
    cover_short = cover_short.reindex(index=s_e).fillna(value=0).values.reshape(-1, 1)
    loci_clip_short_sm = loci_clip_short_sm.reindex(index=s_e).fillna(value=0).values.reshape(-1, 1)
    loci_clip_short_ms = loci_clip_short_ms.reindex(index=s_e).fillna(value=0).values.reshape(-1, 1)
    info = np.concatenate(
        [cover_long, loci_ins_long, loci_clip_long_sm, loci_clip_long_ms, split_read_long, cover_short,
         loci_clip_short_sm, loci_clip_short_ms], axis=1
    )
    return info


def normalization(x_data):
    oshape = x_data.shape
    x_data = x_data.reshape(-1, 8).astype('float32')
    x_data = x_data - x_data.mean(axis=0)
    x_data /= (np.sqrt(x_data.var(axis=0)) + 1e-10)
    return x_data.reshape(oshape)


def label_data(vcfpath, chr_name_long, start, end, index):
    if 'chr' in chr_name_long:
        chr_name_long = chr_name_long[3:]
    gold = []

    for rec in pysam.VariantFile(vcfpath).fetch():
        if rec.contig != chr_name_long:
            continue
        if rec.info['SVTYPE'] == 'INS':
            gold.append([rec.start, rec.stop, rec.stop - rec.start, 1])
    gold = (pd.DataFrame(gold).sort_values([0, 1]).values).astype('float64')
    y = []
    for rec in index:
        if ((gold[:, 1:2] > rec) & (gold[:, :1] < (rec + WINDOW))).sum() != 0:
            y.append((((gold[:, 1:2] >= rec) & (gold[:, :1] <= (rec + WINDOW))) * gold[:, 3:]).sum())
        else:
            y.append(0)

    return (np.array(y) > 0).astype('float32')


def creat_data_longshort(bamfile_long_path, bamfile_short_path, outputpath, contig, contiglength):
    """
    提取长短读数的ins特征   8个特征
    loci_coverage_long、loci_insertion_long、loci_lms_long, loci_lsm_long, loci_lts_long, loci_sc、loci_sms, loci_ssm
    :param bamfile_long_path:
    :param bamfile_short_path:
    :param outputpath:
    :param contig:
    :param contiglength:
    :return:
    """
    bamfile_short = pysam.AlignmentFile(bamfile_short_path, 'rb', threads=20)
    bamfile_long = pysam.AlignmentFile(bamfile_long_path, 'rb', threads=20)

    for ct in contig:
        print("{} is over ".format(ct))
        chr_name_long = ct
        chr_name_short = ct
        chr_length = contiglength[ct]
        ider = math.ceil(chr_length / LENGTHS)
        start = 0
        end = LENGTHS
        s = 0
        # print(f'chr_name_long:{chr_name_long}\tchr_length:{chr_length}\tider:{ider}')
        split_read_long = calculate_splitread(chr_name_long, bamfile_long)
        # print('split_read_long is over !!!')
        for n in range(ider):
            x_data = []
            index = []
            loci_cover_long, loci_ins_long, loci_clip_long_sm, loci_clip_long_ms = feature_extraction_long(bamfile_long,
                                                                                                           chr_name_long,
                                                                                                           start, end,
                                                                                                           0)
            if len(loci_cover_long) == 0:
                start = start + LENGTHS
                end = end + LENGTHS
                continue
            # print('feature_extraction_long is over !!!')
            loci_cover_short, loci_clip_short_sm, loci_clip_short_ms = feature_extraction_short(bamfile_short,
                                                                                                chr_name_short, start,
                                                                                                end, 0)
            # print('feature_extraction_short is over !!!')
            xx = compute(loci_cover_long, loci_ins_long, loci_clip_long_sm, loci_clip_long_ms, split_read_long,
                         loci_cover_short, loci_clip_short_sm, loci_clip_short_ms, start, end)
            # print('computer is over !!!')
            # print('xx.shape:{}'.format(xx.shape))
            xx_test = xx[:, :5].reshape(-1, 1000)
            # print('xx_test.shape:{}'.format(xx_test.shape))
            xx = xx.reshape(-1, 1600)
            # print('xx.shape:{}'.format(xx.shape))
            for k in range(xx_test.shape[0]):
                if xx_test[k].any() != 0:
                    x_data.append(xx[k])
                    index.append(s)
                s = s + WINDOW
            x_data = np.array(x_data)
            index = np.array(index)
            start = start + LENGTHS
            end = end + LENGTHS
            if len(x_data) == 0:
                continue
            # print('chr_name_long:{}\t x_data.shape:{}\tindex.shape:{}'.format(chr_name_long, x_data.shape, index.shape))
            x_data = normalization(x_data)
            # y_label = label_data('/home/yuyanan/data/vcf/hg002.ins_1-22.vcf.gz', chr_name_long, start - LENGTHS,
            #                      end - LENGTHS, index)
            # y_label = y_label.reshape(-1, 1)
            # x_data = np.concatenate([x_data, y_label], axis=1)
            feature_fold = os.path.join(outputpath, 'ont_48x')
            if not os.path.exists(feature_fold):
                os.mkdir(feature_fold)

            if 'chr' in chr_name_long:
                filename_data = os.path.join(feature_fold, f'{chr_name_long}_{start - LENGTHS}_{end - LENGTHS}.npy')
                filename_index = os.path.join(feature_fold,
                                              f'{chr_name_long}_{start - LENGTHS}_{end - LENGTHS}_index.npy')
            else:
                filename_data = os.path.join(feature_fold, f'chr{chr_name_long}_{start - LENGTHS}_{end - LENGTHS}.npy')
                filename_index = os.path.join(feature_fold,
                                              f'chr{chr_name_long}_{start - LENGTHS}_{end - LENGTHS}_index.npy')
            np.save(filename_data, x_data)
            np.save(filename_index, index)
            print(f'chrom:{chr_name_long} start:{start - LENGTHS} end:{end - LENGTHS} is exacted')



def creat_data_longshort_by_pool(bamfile_long_path, bamfile_short_path, outputpath, contig, max_work=5):
    bamfile = pysam.AlignmentFile(bamfile_long_path, 'rb', threads=20)
    contiglength = {}
    if len(contig) == 0:
        contigs = []
        for index in range(len(bamfile.get_index_statistics())):
            if bamfile.lengths[index] >= 40000000:
                contigs.append(bamfile.get_index_statistics()[index].contig)
                contiglength[bamfile.get_index_statistics()[index].contig] = bamfile.lengths[index]
        # print(contiglength)
    else:
        contigs = np.array(contig).astype(str)
        for index in range(len(bamfile.get_index_statistics())):
            if bamfile.lengths[index] >= 40000000:
                contiglength[bamfile.get_index_statistics()[index].contig] = bamfile.lengths[index]
    # print(contiglength)
    pool = mp.Pool(processes=max_work)
    for contig in contigs:
        pool.apply_async(creat_data_longshort,
                         args=(bamfile_long_path, bamfile_short_path, outputpath, [contig], contiglength))
    pool.close()
    pool.join()


"""
 ===========================================================================================================
 
 ===========================================================================================================
"""

# 82服务器
# bamfile_long_path = "/public_data/SV/long_read/HG002_PB_70x_RG_HP10XtrioRTG.bam"
# bamfile_short_path = "/public_data/SV/short_read/HG002.hs37d5.2x250.bam"
# outputpath = ''
# contig = ['1']
# max_work = 5

# 85服务器
# bamfile_long_path = "/public_data/SV/SV/long_read/HG002_GRCh37_ONT-UL_UCSC_20200508.phased.MD.bam"
# bamfile_short_path = "/public_data/SV/SV/short_read/HG002.hs37d5.2x250_downsample50.bam"
# outputpath = '/home/yuyanan/SV/ls_feature/ont_predict'
# contig = [str(i) for i in range(12, 23)]
# # contig = ['12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']
# max_work = 10
#
# creat_data_longshort_by_pool(bamfile_long_path, bamfile_short_path, outputpath, contig, max_work)
