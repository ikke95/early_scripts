#!/usr/bin/env python3
"""
Author: Kristy
Script to assess and trim the read quality of a fastq file.
usage: python3 Assess_qaulity_trimmer.py [fastq_file]
    fastq_file: string, sample information in fastq format.
"""
# import statements
from sys import argv
import os.path
import subprocess

def parse_file(fastq_file):
    """parses fastq file to get quality and fasta sequence

    Key arguments:
        fastq_file -- string, sample information in fastq format.
    Returns:
        quality_score -- list with strings, with quality score.
        fasta_sequence -- list with strings, with sequence in fasta \
        format.
    """
    quality_score = []
    fasta_sequence = []
    quality = False
    sequence = False
    for line in fastq_file:
        if line.startswith("+"):
            quality = True
        elif quality:
            quality_line = line.rstrip('\n')
            quality = False
            quality_score.append(quality_line)
        if line.startswith("@"):
            sequence = True
        elif sequence:
            sequence_line = line.rstrip('\n')
            sequence = False
            fasta_sequence.append(sequence_line)
    return quality_score, fasta_sequence
 
def calculate_length_and_quality(fasta_sequence, quality_score):
    """calculate read quality and read length for each sample

    Key arguments:
        quality_score -- list with strings, with quality score.
        fasta_sequence -- list with strings, with sequence in fasta \
        format.
    Returns:
    max_sequence, min_sequence, average_sequence, ascii_score, \
    average_quality_score
        max_sequence -- integer, size of longest sequence.
        min_sequence -- integer, size of smallest sequence.
        average_sequence -- integer, average size of sequence.
        ascii_score -- list of list with integers, with quality score.
        average_quality_score -- list with integers, with average \
        quality.
    """
    average_quality_score = []
    ascii_numbers = []
    ascii_score = []
    max_sequence = len(max(fasta_sequence, key=len))
    min_sequence = len(min(fasta_sequence, key=len))
    average_sequence = sum(len(sequence)for sequence in \
    fasta_sequence)/len(fasta_sequence)
    for lines in quality_score:
        for character in lines:
            ascii_numbers.append(ord(character)-64)
        ascii_score.append(ascii_numbers)
        ascii_numbers = []  
    for index in range(max_sequence):
        amount_of_reads=0
        total_score=0
        for score in ascii_score:
            if len(score) > index:
                total_score += score[index]
                amount_of_reads += 1
        average=total_score/amount_of_reads
        average_quality_score.append(average)
    return max_sequence, min_sequence, average_sequence, ascii_score, \
    average_quality_score
    
def quality_trimm(fastq_file):
    """uses fastq_quality_trimmer to trim the reads with a low quality.

    Key arguments:
        fastq_file -- string, sample information in fastq format.
    Returns:
        trimmed_fq -- string, trimmed sample information in fastq format
    """
    trimmed_fq = "trimmed_V2.fq"
    ASCII = 64
    if os.path.exists(trimmed_fq):
        return trimmed_fq
    cmd = 'fastq_quality_trimmer -t 30 -i %s -o %s  -Q %s'\
    %(fastq_file, trimmed_fq, ASCII)
    e = subprocess.check_call(cmd, shell=True)
    print('EXIT STATUS AND TYPE', e, type(e))
    return trimmed_fq

def compare_quality_trimmed_and_untrimmed(average_quality_score, \
trimmed_average_quality_score):
    """compares the quality of the trimmed and untrimmed reads and
    calculates the quality increase.

    Key arguments:
        average_quality_score -- list with integers, with average \
        quality.
        trimmed_average_quality_score -- list with integers, with \
        average quality of trimmed sequences.
    Returns:
        difference_in_quality -- list with integers, difference in \
        quality score between trimmed and untrimmed sequences.
    """
    difference_in_quality = []
    for element in range(len(average_quality_score)):
        difference = trimmed_average_quality_score[element] - \
        average_quality_score[element] 
        difference_in_quality.append(difference)
    return difference_in_quality

def main():
    """Main function of this module"""
    fastq_file = argv[1]
    quality_score, fasta_sequence = parse_file(open(fastq_file))
    max_sequence, min_sequence, average_sequence, ascii_score, \
    average_quality_score = calculate_length_and_quality\
    (fasta_sequence, quality_score)
    trimmed_fq = quality_trimm(fastq_file)
    trimmed_quality_score, trimmed_fasta_sequence = parse_file\
    (open(trimmed_fq))
    trimmed_max_sequence, trimmed_min_sequence, \
    trimmed_average_sequence, trimmed_ascii_score, \
    trimmed_average_quality_score = calculate_length_and_quality\
    (trimmed_fasta_sequence, trimmed_quality_score)
    difference_in_quality = compare_quality_trimmed_and_untrimmed\
    (average_quality_score, trimmed_average_quality_score)
    print("ORIGINAL: min=",min_sequence, "max=",max_sequence, \
    "avg=",average_sequence)
    print("TRIMMED: min=",trimmed_min_sequence, "max=",\
    trimmed_max_sequence, "avg=",trimmed_average_sequence)
    for index in range(len(average_quality_score)):
        output_list=[str(index + 1),str(average_quality_score[index]),\
         str(trimmed_average_quality_score[index]), \
         str(difference_in_quality[index])]
        print("\t".join(output_list))

if __name__ == "__main__":
    main()
