#!/usr/bin/env python3
"""
Author: Kristy Joosten
Script to calculate the rpkm (reads per kilo base per million
transcripts) and list the top 10 most highly expressed genes.
usage:python3 Kallisto_parser.py [index_file] [fastq_file]
    index_file -- _io.TextIOWrapper, binary file that contains \
    index information.
    fastq_file -- _io.TextIOWrapper, file that contains the reads \
    in fastq format.
"""
# import statements
from sys import argv
import os.path
import subprocess

def run_kallisto(index_file, fastq_file):
    """run kallisto through command line

    Key arguments:
        index_file -- _io.TextIOWrapper, binary file that contains \
        index information.
        fastq_file -- _io.TextIOWrapper, file that contains the reads \
        in fastq format.
    Returns:
        out_kallisto -- _io.TextIOWrapper, lists for each transcript \
        its length, effective length, estimated counts and tpmw.
    """
    out_kallisto = "abundance.tsv"
    index_file = argv[1]
    fastq_file = argv[2]
    length = 200
    sd = 20
    if os.path.exists(out_kallisto):
        return out_kallisto
    cmd = "kallisto quant -i %s -o %s -l %s -s %s --single %s" \
    %(index_file, out_kallisto, length, sd, fastq_file)
    e = subprocess.check_call(cmd, shell=True)
    print('EXIT STATUS AND TYPE', e, type(e))
    return out_kallisto
    

def parse_kallisto(kallisto_file):
    """parses Kallisto file

    Key arguments:
        kallisto_file -- _io.TextIOWrapper, lists for each transcript \
        its length, effective length, estimated counts and tpm.
    Returns:
        eff_length -- list, effect length of each transcript
        sequence_id -- list, id of sequence of each transcript
        est_counts -- list, estimeted counts of each transcript
        tpm -- list, normalized expression fort transcripts per million
    """
    sequence_id = []
    eff_length = []
    est_counts = []
    tpm = []
    kalisto_list = kallisto_file.splitlines()
    kallisto_clean = kalisto_list[1:]
    kallisto_clean = [i.split('\t') for i in kallisto_clean] 
    for items in kallisto_clean:
        seq_id = items[0]
        eff = items[2]
        est = items[3]
        tpms = items[4]
        sequence_id.append(seq_id)
        eff_length.append(eff)
        est_counts.append(est)
        tpm.append(tpms)
    return eff_length, sequence_id, est_counts, tpm
 
def calculate_RPKM(est_counts, eff_length):
    """calculates the RPKM

    Key arguments:
        eff_length -- list, effect length of each transcript
        est_counts -- list, estimeted counts of each transcript
    returns:
        RPKM -- list,reads per kilo base per million transcripts
    """
    est_counts = [ float(element.strip()) for element in est_counts]
    eff_length = [ float(element.strip()) for element in eff_length]
    sum_est_counts = sum(est_counts)
    RPKM = [(count/(length*sum_est_counts))*10000000000 for count, \
    length in zip(est_counts, eff_length)]
    return RPKM

def calculate_TPM(RPKM):
    """ calculates the TPM
    
    Key arguments:
        RPKM -- list,reads per kilo base per million transcripts
    returns:
        TPM -- list, with the TPM values
    """
    sum_RPKM = sum(RPKM)
    TPM = [(rpkm/sum_RPKM)*1000000 for rpkm in RPKM]
    return TPM
    
def validate_TPM(TPM, tpm, sequence_id):
    """ compares the calculated TPM with TPM from tool
    
    key values:
        TPM -- list, with the TPM values
        tpm -- list, normalized expression fort transcripts per million
        sequence_ID -- list, id of sequence of each transcript
    returns:
        validation -- list, with values that differ more then 0.1 form \
        each other
     """
    tpm = [ float(element.strip()) for element in tpm]
    TPM = [round(num, 2) for num in TPM]
    tpm = [round(num, 2) for num in tpm]
    validation = []
    for i in range(len(TPM)): 
        if (TPM[i] - tpm[i]) > 0.1:
            validation.append([sequence_id[i], TPM[i]])
    return validation

def rpkm_gene(sequence_id, RPKM):
    gene_RPKM_dict = {}
    for i in range(len(RPKM)):
        try:
            prev_rpkm = gene_RPKM_dict[sequence_id[i].split('.')[0]]
            gene_RPKM_dict[sequence_id[i].split('.')[0]] = RPKM[i] + prev_rpkm
        except KeyError:
            gene_RPKM_dict[sequence_id[i].split('.')[0]] = RPKM[i]
    return gene_RPKM_dict

def create_output(gene_RPKM_dict, validation):
    """create output    
    key arguments:
    sequence_ID -- list, id of sequence of each transcript
    RPKM -- list,reads per kilo base per million transcripts
    Returns:
    none """
    sorted_gene_RPKM_dict = sorted(gene_RPKM_dict.items(), \
    key=lambda item: item[1])
    sorted_gene_RPKM_dict.reverse()
    top_5 = sorted_gene_RPKM_dict[:5]
    print("Given and calculated tpm values differ > 0.1 for:")
    for i in range(len(validation)):
        print(str(validation[i][0]))
    for i in range(len(top_5)):
        print(str(i+1)+"\t"+ str(top_5[i][0]) + "\t" + str(top_5[i][1]))
    return
    
def main():
    """Main function of this module"""
    index_file = argv[1]
    fastq_file = argv[2]
    out_kallisto = run_kallisto(open(index_file), open(fastq_file))
    kallisto_file = open("abundance.tsv/abundance.tsv", mode='r', \
    encoding='UTF-8')
    kallisto_file = kallisto_file.read()
    eff_length, sequence_id, est_counts, tpm = \
    parse_kallisto(kallisto_file)
    RPKM = calculate_RPKM(est_counts, eff_length)
    TPM = calculate_TPM(RPKM)
    validation = validate_TPM(TPM, tpm, sequence_id)
    gene_RPKM_dict = rpkm_gene(sequence_id, RPKM)
    create_output(gene_RPKM_dict,validation)

if __name__ == "__main__":
    main()
