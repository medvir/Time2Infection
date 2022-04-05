# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 10:13:44 2022

Input is a fastq file or a vcf file. Extract the variations in the PR and RT region of 
the HXB2 genome. Then calculate the APD score and ambiguity score based on paper
https://doi.org/10.1093/infdis/jiz094  for the third codon positions. 
If fastq file is given, it is filtered for low quality reads, mapped to HXB2 genome, VCF file 
is generated. 
If vcf file is given it should be based on HXB2 genome.     

@author: maryamzaheri
"""

import argparse
from argparse import RawTextHelpFormatter
import subprocess
import pandas as pd 
import os
import gzip
import shutil
import sys

SNV_MIN_COV = 100
MIN_ALLEL_FREQ = 0.01
substr_th = str(MIN_ALLEL_FREQ).replace('.', 'p')
SCORE_FILE = 'score_%s_PR_RT_3rdcodonPos.txt' %substr_th
start_pos = 2253  #start of PR
end_pos = 3870 #end of RT  # end of INT 5096  
PR_RT_len = end_pos - start_pos # HXB2 PR/RT  3870 - 2253 

def filter_low_quality(fastq_file):
    # check for read length, check for quality score, trim low quality bases
    # min_len 0.8 * 150 = 120 results to too many sequence loss, so I take 70% of the length. 
    hq_file_prefix = 'hq_R'
    
    prinseq_extra_param = '-ns_max_n 4 -min_qual_mean 30 -trim_qual_left 30 ' \
    '-trim_qual_right 30 -trim_qual_window 10'
    
    command = 'prinseq -fastq %s -out_format 3 -out_good %s -out_bad bad_q ' \
    '-min_len %d %s -log prinseq_log' %(fastq_file, hq_file_prefix, 105, prinseq_extra_param)
    subprocess.call(command, shell=True)
    
    hq_file = '%s.fastq' %hq_file_prefix
    
    return hq_file

def align(hq_fastq_file, ref):
    
    command = 'samtools faidx %s' %ref
    subprocess.call(command, shell=True)
    
    command = 'bwa index %s' %ref
    subprocess.call(command, shell=True)
    
    aligned_file = 'aln_ref.sam'
    command = 'bwa mem -t %d %s %s > %s' %(4, ref, hq_fastq_file, aligned_file)
    subprocess.call(command, shell=True)
    
    aligned_bam_file = 'aln_ref.bam'
    command = 'samtools view -h -F 4 -F 2048 %s > %s' %(aligned_file, aligned_bam_file)
    subprocess.call(command, shell=True)
    
    
    return aligned_bam_file

def generate_vcf(aln_bam_file, ref):
    
    indelqual_file = 'aligned_indelqual.bam'
    command = 'lofreq indelqual --dindel -f %s -o %s %s' % (ref, indelqual_file, aln_bam_file)
    subprocess.call(command, shell=True)
    
    indelqual_sorted_file = 'indelqual_sorted_file.bam'
    command = 'samtools sort -T tmp -o %s %s' %(indelqual_sorted_file, indelqual_file)
    subprocess.call(command, shell=True)
    
    command = 'samtools index %s' %indelqual_sorted_file
    subprocess.call(command, shell=True)
    
    snv_file = 'snv.vcf'
    command = 'lofreq call-parallel --pp-threads 4 --call-indels -f %s -o %s %s' %(ref, 'tmp.vcf', indelqual_sorted_file)
    subprocess.call(command, shell=True)
    
    command = 'lofreq filter -i tmp.vcf -o %s -v %d -a %f' %(snv_file, SNV_MIN_COV, MIN_ALLEL_FREQ)
    subprocess.call(command, shell=True)
    return snv_file

def parse_vcf(snv_file):
    
    col_names = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER', 'INFO']
    # read vcf file
    snv_df = pd.read_csv(snv_file, sep='\t', header=None, comment='#', names=col_names)
    
    pol_PRRT_snv_df = snv_df.loc[(snv_df['POS'] >= start_pos) & (snv_df['POS'] <= end_pos)] #3870
    
    pol_PRRT_snv_df['AF'] = pol_PRRT_snv_df['INFO'].str.extract(r'AF=([\d]*\.[\d]*);', expand=True).astype(float)
    
    # Drop positions with freq 1
    pol_snv_woOnes_df = pol_PRRT_snv_df[pol_PRRT_snv_df['AF'] != 1.00]
    
    # Drop 1st and 2nd codon positions
    pol_snv_woOnes_3rdcodon_df = pol_snv_woOnes_df[((pol_snv_woOnes_df['POS'] - (start_pos + 2))%3) == 0]
    
    return pol_snv_woOnes_df, pol_snv_woOnes_3rdcodon_df

def calc_amb_scores(snv_df):
    
    amb_score = 0
    for idx in range(snv_df.shape[0]):
        pointer2 = idx
        freq_ls = []
        while pointer2 < snv_df.shape[0] and snv_df.iloc[pointer2]['POS'] == snv_df.iloc[idx]['POS']:
            freq_ls.append(snv_df.iloc[pointer2]['AF'])
            pointer2 += 1
            
        freq_ls.append(1 - sum(freq_ls))
        
        max_freq = max(freq_ls)
        
        if 1 - max_freq > MIN_ALLEL_FREQ:
            amb_score += 1
            
    amb_score = (amb_score / PR_RT_len)
    
    return amb_score

def calc_APD(snv_df):
    APD = 0
    for idx in range(snv_df.shape[0]):
        pointer2 = idx
        freq_ls = []
        while pointer2 < snv_df.shape[0] and snv_df.iloc[pointer2]['POS'] == snv_df.iloc[idx]['POS']:
            freq_ls.append(snv_df.iloc[pointer2]['AF'])
            pointer2 += 1
            
        freq_ls.append(1 - sum(freq_ls))
        
        max_freq = max(freq_ls)
        diversity_contribution_sum = 0
        if 1 - max_freq > MIN_ALLEL_FREQ:
            for nc_freq in freq_ls:
                diversity_contribution_sum += (nc_freq * (1 - nc_freq))
            
        APD += diversity_contribution_sum
    
    APD = (APD / PR_RT_len)
    
    return APD

def main(input_file, ref_path):
    
    ref = os.path.abspath(ref_path)
    dir_path, file_name = os.path.split(os.path.abspath(input_file))
    os.chdir(dir_path)
    
    if file_name.endswith('.fastq') or file_name.endswith('.gz'):
        
        if file_name.endswith('.gz'):
            root, extension = os.path.splitext(file_name)
            #fastq_gz_file = fastq_file
            fastq_file = root
            with gzip.open(file_name, 'rb') as f_in:
                with open(fastq_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            fastq_file = file_name
        
        ## Filter for low-quality reads
        hq_fastq_file = filter_low_quality(fastq_file)
    
        # downsample if needed 
        ## Align to HXB2
        aligned_file = align(hq_fastq_file, ref)
        
        ## generate vcf file
        snv_file = generate_vcf(aligned_file, ref)
    elif file_name.endswith('.vcf'):
        snv_file = file_name
    
    else:
        sys.exit('The file extension is unknown (fastq, fastq.gz and vcf extensions are accepted).')
    
    pol_snv_woOnes_df, pol_snv_woOnes_3rdcodon_df = parse_vcf(snv_file)
    # Calculate APD and Ambigiuty score
    amb_score = calc_amb_scores(pol_snv_woOnes_df)
    APD = calc_APD(pol_snv_woOnes_3rdcodon_df)
    
    print(APD, amb_score)
    
    if os.path.exists(SCORE_FILE):
        line = '%s , %f , %f \n' %(file_name, APD, amb_score)
        with open(SCORE_FILE, "a+") as file_object:
            file_object.write(line)
    else:
        title = "input_file, APD, amb_score \n"
        line = '%s , %f , %f \n' %(file_name, APD, amb_score)
        with open(SCORE_FILE, "a+") as file_object:
            file_object.write(title)
            file_object.write(line)
    
    file_object.close()

    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='given the fastq file or snv file calculate ambiguity score and APD score in the score*.txt file.'
                                     , formatter_class=RawTextHelpFormatter)
    
    group = parser.add_mutually_exclusive_group()
    
    group.add_argument('-f','--file', dest='file', type=str, 
                       help='the fastq file name')
    group.add_argument('-s','--snv', dest='snv', type=str, 
                       help='the snv file with respect to HXB2')
    
    parser.add_argument('-r','--ref', dest='ref', type=str, 
                        default='HXB2.fasta',
                        help='the reference file name')

   
    args = parser.parse_args()
    
    if args.file:
        main(input_file=args.file, ref_path=args.ref )
    else:
        main(input_file=args.snv, ref_path=args.ref )


