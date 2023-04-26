#! /usr/bin/python3.10
import sys
##unique 128-bit values popularly used to uniquely identify entities on the internet
import uuid
##gzip reduces the size of the named files using Lempel–Ziv coding (LZ77). Whenever possible, each file is replaced by one with the extension ‘.gz’
import gzip
import random
import string
## Matplotlib is a comprehensive library for creating static, animated, and interactive visualizations in Python
####matplotlib.pyplot is a collection of functions; Each pyplot function makes some change to a figure: e.g., creates a figure, creates a plotting area in a figure, plots some lines in a plotting area...
from matplotlib import pyplot as plt
import matplotlib
##subprocess lets you start new applications right from the Python program you are currently writing. So, if you want to run external programs from a git repository you can use subprocess in Python.
import subprocess as sb
import re
##The OS module in Python provides functions for interacting with the operating system
import os
##bio.Entrez rovides code to access NCBI over the WWW
###This module provides a number of functions like efetch (short for Entrez Fetch) which will return the data as a handle object (Handle objects enable more than one variable to refer to the same object.). 
###This is a standard interface used in Python for reading data from a file, or in this case a remote network connection, and provides methods like .read() or offers iteration over the contents line by line.
from Bio import Entrez
##Bio.SeqIO provides a simple uniform interface to input and output assorted sequence file formats (including multiple sequence alignments), but will only deal with sequences as SeqRecord objects. 
##There is a sister interface Bio.AlignIO for working directly with sequence alignment files as Alignment objects.
from Bio import SeqIO
#Create a commandline for the NCBI BLAST+ program blastp.
####NcbiblastnCommandline - Nucleotide-Nucleotide BLAST
from Bio.Blast.Applications import NcbiblastnCommandline 
###multiprocessing is a package that supports spawning processes using an API similar to the threading module. 
####The multiprocessing package offers both local and remote concurrency, effectively side-stepping the Global Interpreter Lock by using subprocesses instead of threads. 
###Due to this, the multiprocessing module allows the programmer to fully leverage multiple processors on a given machine. (???)
###A prime example of this is the Pool object which offers a convenient means of parallelizing 
###the execution of a function across multiple input values, distributing the input data across processes (data parallelism).
from multiprocessing import Pool
##This module provides classes and functions for comparing sequences. 
##It can be used for example, for comparing files, and can produce information about file differences in various formats, including HTML
##class difflib.SequenceMatcher: This is a flexible class for comparing pairs of sequences of any type, so long as the sequence elements are hashable.
from difflib import SequenceMatcher
#The tqdm module allows for the generation of progress bars in Python. (barra di caricamento)
from tqdm import tqdm
##Tempfile is a Python module used in a situation, where we need to read multiple files, change or access the data in the file, and gives output files based on the result of processed data.
##Each of the output files produced during the program execution was no longer needed after the program was done.
##In this situation, tempfiles are used to create temporary files so that next time we don’t have to find delete when our program is done with them
import tempfile
##The argparse module makes it easy to write user-friendly command-line interfaces. 
###The program defines what arguments it requires, and argparse will figure out how to parse those out of sys.argv. 
###The argparse module also automatically generates help and usage messages. 
import argparse
##Pandas is an open source library in Python. It provides ready to use high-performance data structures and data analysis tools. (?)
import pandas as pd
##Matplotlib has a number of built-in colormaps accessible via matplotlib.colormaps. However, we often want to create or manipulate colormaps in Matplotlib. This can be done using the class ListedColormap or LinearSegmentedColormap. 
# Seen from the outside, both colormap classes map values between 0 and 1 to a bunch of colors
##Colormap objects based on lookup tables using linear segments.
from matplotlib.colors import LinearSegmentedColormap
import csv
##The Environment for Tree Exploration (ETE) is a Python programming toolkit that assists in the recontruction, manipulation, analysis and visualization of phylogenetic trees 
###(although clustering trees or any other tree-like data structure are also supported).
##ETE’s ncbi_taxonomy module provides utilities to efficiently query a local copy of the NCBI Taxonomy database. 
# The class NCBITaxonomy offers methods to convert from taxid to names (and vice versa), to fetch pruned topologies connecting a given set of species, or to download rank, names and lineage track information.
from ete3 import NCBITaxa

##matplotlib.style.use Use Matplotlib style settings from a style specification.
matplotlib.style.use('ggplot')
#ktImportTaxonomy commonly used for the analysis and visualization of phylogenetic trees. 
#The "-tax" option is used to specify the type of taxonomic data that is being imported.For example, you might use "-tax ncbi" to specify that the taxonomic data is in NCBI Taxonomy format
#The "-o" option is used to specify the output file for the annotated tree. This option accepts the path and filename of the output file. 
KTIMPORTTAX = "ktImportTaxonomy -tax %s -o %s %s"
##Minimap2 is a software tool used for fast and accurate pairwise alignment of DNA or RNA sequences
##MINIMAP: ALIGNMENT NANOPORE READS AGAINST EACH OTHER
#"-a" specifies that the output should be in the SAM alignment format.
#"-x ava-ont" specifies the preset alignment mode to use, in this case, "ava-ont" which is optimized for aligning Oxford Nanopore Technologies (ONT) reads against each other.
#"-t" specifies the number of threads (or CPU cores) to use for the alignment process.
MINIMAP = "minimap2 -a -x ava-ont -t %s %s %s "
#MINIMAP_S: ALIGNMENT NANOPORE READS AGAINST A REFERENCE GENOME
#"-x map-ont" specifies the preset alignment mode to use, in this case,"map-ont" which is optimized for aligning Oxford Nanopore Technologies (ONT) reads against a reference genome or transcriptome.
#"--secondary=no" specifies that secondary alignments (i.e., alignments that are not the primary alignment for a given read) should not be included in the output.
MINIMAP_S = "minimap2 -a -x map-ont --secondary=no -t %s %s %s "
MINIMAP_SY = "minimap2 -a -x map-ont -t %s %s %s "
##Samtools sort is a software tool used for sorting sequencing reads in a SAM or BAM file by their alignment position on a reference genome
##Samtools sort offers several options for customizing the sorting process, such as specifying the number of threads to use for sorting, sorting by read name instead of alignment position
##Sorting the reads by their alignment position can be useful for several reasons. For example, it can make it easier to identify reads that overlap a specific genomic region of interest, or to identify variants that are present in a specific genomic region.
#"-o" specifies the output file path for the sorted BAM file.
##"-O BAM" specifies the output file format for the sorted file, in this case, "BAM"
SAMSORT = "samtools sort -o %s -O BAM"
SAMINDEX = "samtools index %s"
##Racon is a software tool for polishing genome assemblies using long read sequencing data. Its main function is to correct errors in an initial draft assembly using long-read sequencing data from PacBio or Oxford Nanopore Technologies (ONT) platforms.
##First, Racon maps the long reads back to the initial draft assembly using a fast and efficient alignment algorithm. Next, it identifies discrepancies between the draft assembly and the long reads, such as mismatches, insertions, and deletions. Finally, 
##Racon generates a new consensus sequence for each region of the assembly, based on the majority vote of the aligned long reads.
#"t": threads
#"-f" specifies the input file path for the FASTQ or FASTA format reads to be used for polishing the assembly. This option tells Racon which long read sequencing data to use to correct errors in the initial draft assembly.
RACON = "racon -t %s -f %s %s %s"
##porechop is a software tool used for adapter trimming and demultiplexing of Oxford Nanopore Technologies (ONT) sequencing data (the second one separating sequencing reads into different files based on their barcode tags). 
##The primary function of Porechop is to remove adapter sequences and barcode tags from raw ONT sequencing data
##"-i" specifies the input file path for the ONT sequencing data to be processed. "-t": threads -o" specifies the output file path for the processed sequencing data. 
PORECHOP = "porechop -i %s -t %s -o %s"
##Jellyfish is a software tool used for counting k-mers in DNA or RNA sequencing data. K-mers are substrings of length k that are extracted from sequencing reads
##it's a command-line program that reads FASTA and multi-FASTA files containing DNA sequences. It outputs its k-mer counts in a binary format (hash table), which can be translated into a human-readable text format using the "jellyfish dump" command
##The primary function of Jellyfish is to efficiently count the occurrences of k-mers in large sequencing datasets, to count the number of occurrences of each k-mer in the reads.
##The output of Jellyfish is a histogram of k-mer frequencies, showing the number of k-mers that occur at each frequency in the input reads. It's used for estimating genome size, assessing sequencing quality, or identifying repeat regions in a genome.
##"count" is the command to run the Jellyfish k-mer counting function.
##"-m 100" specifies the k-mer length to use for counting. In this case, k-mers of length 100 will be counted.
##"-s 100M" specifies the size of the hash table to use for counting (100 million) 
##"-t 10": threads
##"-o /dev/stdout" specifies the output file path for the counted k-mer frequencies. In this case, the output will be sent to the standard output stream (?) (stdout) instead of a file.
##"-C" specifies that Jellyfish should canonicalize the k-mers before counting. This means that each k-mer and its reverse complement will be counted as a single k-mer, reducing the memory usage and improving the accuracy of the k-mer frequency estimation.
JELLYFISH_COUNT= "jellyfish count -m 100 -s 100M -t 10 -o /dev/stdout -C %s"
##is used to write the k-mer count data to the standard output (stdout) stream, which is typically displayed in the terminal window or redirected to another program or file using pipes.
##The "/dev/stdin" part of the command specifies that the k-mer count data should be read from the standard input stream (?).
JELLYFISH_DUMP= "jellyfish dump -o %s /dev/stdin"
##spades is a popular software package for de novo genome assembly of bacterial, fungal, and small eukaryotic genomes from next-generation sequencing data.
##The main function of the SPAdes package is to take short reads generated by next-generation sequencing technologies and assemble them into longer contiguous sequences, called contigs.
##It first performs error correction on the reads to reduce sequencing errors, and then it constructs the De Bruijn graph, a data structure that represents the overlaps between reads. 
##The algorithm then simplify the graph and identify potential paths that correspond to genomic sequences. Finally, it generates contigs from these paths and tries to scaffold them into larger sequences using paired-end reads or other information.
#"-s": specifies the path to the input sequencing data in single-end mode. This means that the input reads are from a single end of the DNA fragments. You would need to replace the dash with the path to your input file, for example, "-s /path/to/reads.fastq"
#"-o": specifies the output directory where the assembly results will be saved
#"--only-assembler": specifies that SPAdes should only perform the assembly step and not perform any downstream analysis, such as gene prediction or functional annotation. 
SPADES = "spades -s %s -o %s --only-assembler"
##Nanopolish is specifically designed to analyze the raw signal data produced by nanopore sequencers to perform: Basecalling; Variant calling: Nanopolish can detect single nucleotide variants (SNVs), insertion-deletion variants (indels), and structural variants (SVs) in DNA or RNA sequencing data by comparing the raw signal data to a reference genome or draft assembly...
##"--consensus": This option tells Nanopolish to generate a consensus sequence from the input reads and the detected variants.
##"-o": This option specifies the output file where the variants will be written
#"-w" This option sets the size of the window used to call variants. A larger window size can improve sensitivity but decrease specificity.
#"-r": This option specifies the format of the file where there are the input reads
#"-b <alignments.bam>": This option specifies the input alignments in BAM format
#"-g": This option specifies the format of the file where the reference genome is stored
#"-t":number of threads
#"--min-candidate-frequency 0.1": This option sets the minimum frequency of a candidate variant in the input reads to be considered for calling. A lower frequency threshold can increase sensitivity but decrease specificity.
#"-p 1": This option enables polishing of the consensus sequence using the input reads and the detected variants.
#"--fix-homopolymers": This option enables a homopolymer correction algorithm that can improve accuracy in regions with long stretches of the same nucleotide.
NANOPOLISHV = "nanopolish variants --consensus -o %s -w %s -r %s -b %s -g %s -t %s --min-candidate-frequency 0.1 -p 1 " \
             "--fix-homopolymers"
#The "nanopolish index -d" command is used to create an index for a directory of raw nanopore sequencing data files in the fast5 format. 
#"nanopolish index": This is the Nanopolish subcommand used to create an index for the raw nanopore sequencing data files.
#"-d": This option specifies that the input is a directory of raw nanopore sequencing data files. After the "nanopolish index -d" command, you should specify the path to the directory containing the fast5 files you want to index
NANOPOLISHI = "nanopolish index -d %s %s "
#The "nanopolish vcf2fasta -g" command is used to generate a consensus sequence from a VCF file containing variants called by Nanopolish.
#"-g": This option specifies that the output should be in FASTA format and include only the consensus sequence.
NANOPOLISHVA = "nanopolish vcf2fasta -g %s %s"
BWAI = "bwa index %s"
#The "bwa mem -t" command is used to align short reads to a reference genome 
#-t: threads
BWA = "bwa mem -t %s %s %s"
#BCFtools is a suite of software tools for working with VCF and BCF files, which are commonly used in bioinformatics to store genomic variation data.
#The "bcftools mpileup -Ou -f" command is used to generate a pileup of mapped reads against a reference genome, which can be used for variant calling using BCFtools (?)
####---> Pileup format is a text-based format for summarizing the base calls of aligned reads to a reference sequence. This format facilitates visual display of SNP/indel calling and alignment. 
##"-Ou": This option specifies that the output should be in uncompressed BCF format and should be written to the standard output. The uncompressed format is used to save disk space and reduce processing time.
#"-f <ref.fa>": This option specifies the path to the reference genome in FASTA format
#After the "bcftools mpileup -Ou -f" command, you should specify the path to one or more input BAM files containing the mapped reads.
BCFTOOLS_MP = "bcftools mpileup -Ou -f %s %s"
#The "bcftools call -mv -o --ploidy 1" command is used to perform variant calling and output the results in VCF format
#"-mv": These options specify that both SNPs and indels should be called, and that only high-confidence variants should be output.
#"-o": This option specifies the path and name of the output VCF file.
#"--ploidy 1": This option specifies the ploidy of the sample being analyzed. A ploidy of 1 indicates that the sample is haploid, meaning that there is only one copy of each chromosome.
##After the "bcftools call -mv -o --ploidy 1" command, you should specify the path to the input BCF file containing the pileup of mapped reads.
BCFTOOLS_CALL = "bcftools call -mv -o %s --ploidy 1"
BCFTOOLS_IN = "bcftools index %s"
##The "bcftools consensus" command is used to generate a consensus sequence in FASTA format based on a given reference genome and a set of variants from a VCF file.
##"-f": This specifies the path and name of the reference genome in FASTA format.
##"-o": This specifies the path and name of the output consensus sequence in FASTA format.
BCFTOOLS_CO = "bcftools consensus -f %s -o %s %s"
#FreeBayes is a program used for variant discovery and genotyping in diploid genomes. 
# It uses a Bayesian inference model to identify variants from aligned sequencing data (e.g. BAM files) and outputs the discovered variants in a VCF format.
#The command "freebayes -f <reference.fa> -p 1" is used to call variants in a single sample using a specified reference genome.
##"-f" This specifies the path and name of the reference genome in FASTA format.
##"-p 1": This specifies the ploidy of the sample being analyzed. In this case, the sample is assumed to be diploid (ploidy=2), so the option "-p 1" is used to specify that only one individual is being analyzed.
FREEBAYES = "freebayes -f %s -p 1 %s"
#Krocus which can predict a sequence type directly from uncorrected long reads
##First of all you need MLST databases. usage: krocus_database_downloader [options] --species: Species to download (default: None)
KROCUS_D = "krocus_database_downloader --species %s"
#Bgzip is a block-based compression format that is widely used in bioinformatics to compress large genomic data files such as variant call format (VCF) files, BAM files, and FASTQ files.
BGZIP = "bgzip %s"
#The bgzip tool is typically used in conjunction with tabix, a tool for indexing Bgzip-compressed files, to efficiently query and retrieve subsets of data from large genomic data files. 
##The option -p vcf is used to specify that the input file is in VCF format. 
# When tabix is run with this option, it automatically detects the VCF file format and generates an index file (with extension .tbi) that allows for efficient querying of the VCF file.
TABIX = "tabix -p vcf %s"

mlst_database = "https://pubmlst.org/data/dbases.xml"

#initialize an instance of the NCBITaxa class using the NCBITaxa() constructor.
ncbi = NCBITaxa()


def get_desired_ranks(taxid, desired_ranks):
    ####To use the get_lineage() function from the NCBITaxa module, you first need to import the module and initialize an instance of the NCBITaxa class
    ###"ncbi.get_lineage(taxid)"" - retrieves the taxonomic lineage (i.e., the hierarchical relationships of taxonomic categories) for a given taxonomic identifier (ID)
    ##The output of this code will be a list of taxon IDs representing the lineage of organism identified with the taxon ID, starting with the root node and ending with the specified taxon ID
    lineage = ncbi.get_lineage(taxid)
    #- retrieves the taxonomic rank (i.e., the level of a taxonomic category in the hierarchy) for each taxonomic identifier in the lineage 
    lineage2ranks = ncbi.get_rank(lineage)
    ##with "dict" creates a dictionary (ranks2lineage) that maps each taxonomic rank to the corresponding taxonomic identifier in the lineage, 
    ###by iterating (with a for loop) over the items in lineage2ranks and swapping (scambiando) the key-value pairs (because we create the dictionary with rank like key and taxid like value, but in the for loop we invert the order of these two).
    ##items() --> it's used for dictionaries. It returns a view object that contains a list of key-value pairs as tuples(like a list but immutable). Each tuple in the list represents a key-value pair from the dictionary.
    ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
    ##for each taxonomic rank in desired_ranks (a list of taxonomic ranks to be included in the output, which we provide at the beginning as an argument of the function),
    ###creates a dictionary with keys that are formatted strings ({}_id) that represent each taxonomic rank in desired_ranks, 
    ###and values that are obtained from ranks2lineage for each rank, or <not present> if the rank is not found in ranks2lineage
    ###Start with a string that contains one or more placeholder markers that look like {} or {index}. These placeholders are where you want to insert values into the string.
    ###Call the .format() method on the string, passing one or more values as arguments. These values will be inserted into the placeholders in the order they are given. 
    
    return {'{}_id'.format(rank): ranks2lineage.get(rank, '<not present>') for rank in desired_ranks}

##This function takes a list of taxonomic identifiers (taxids), retrieves taxonomic information at the desired ranks for each identifier, and writes the results to a tab-separated CSV file.
##taxids: A list of taxonomic identifiers for which you want to retrieve taxonomic information.
##desired_ranks: A list of taxonomic ranks that you want to include in the output.
##path: A file path where the output will be written.
def main(taxids, desired_ranks, path):
    ##It opens a new file (cvsfile) at the specified path
    with open(path, 'w') as csvfile:
        ###It creates a list of column headers for the output CSV file by iterating over the desired_ranks list and formatting each rank as a string with the suffix _id. 
        # This list of headers is used to create a csv.DictWriter object, which will write rows of taxonomic information to the output file.
        ##It writes the column headers to the output file using the writeheader method of the csv.DictWriter object.
        fieldnames = ['{}_id'.format(rank) for rank in desired_ranks]
        writer = csv.DictWriter(csvfile, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        for taxid in taxids:
            ###For each taxonomic identifier in the taxids list, it retrieves taxonomic information at the desired ranks using the get_desired_ranks function.
            ##t writes the taxonomic information for the current taxonomic identifier to the output file as a new row, using the writerow method of the csv.DictWriter object.
            writer.writerow(get_desired_ranks(taxid, desired_ranks))

#@Gooey(optional_cols=2, program_name="Subparser Layout Demo")
def setting():
    ##argparse makes it easy to write user-friendly command-line interfaces, defines what arguments it requires, also automatically generates help and usage messages.
    ##argparse.ArgumentParser provides a way for developers to define and handle command-line arguments in their Python programs. It allows you to create a command-line interface with options and arguments that can be easily parsed and validated by the program.
    #once an ArgumentParser object is created, you can add arguments and options to it using methods such as add_argument()
    ##When the command-line arguments are passed to the program, argparse will automatically parse and validate them based on the definitions provided in the ArgumentParser object. 
    ##It will also generate a help message for the program based on the arguments and options defined.
    #####The prog parameter sets the name of the program that will be displayed in the help message.
    #The usage parameter sets the usage message that will be displayed in the help message. In this case, it includes the program name from prog and the string "[options]", indicating that the program expects some options to be passed in. 
    #After creating the ArgumentParser object, you would typically call the add_argument() method on it to define the options that the program expects
    parser = argparse.ArgumentParser(prog='quantify and detect pathogens', usage='%(prog)s [options]')
    ##-o and --output are two possible ways to specify the argument on the command line. 
    ##nargs="?" indicates that this argument can take zero or one value. If no value is specified, the default value "output" will be used.
    ##default="output" sets the default value for the argument to "output" if no value is specified.
    ##help="output name" provides a brief description of the argument to the user. This will be displayed if the user requests help documentation for the program. In this case, it tells the user that the argument specifies the output file name.
    parser.add_argument("-o","--output", nargs="?", default="output", help="output name")
    ##--barcode is a long form of the argument used to specify the barcode to be analyzed.
    ##required=True specifies that this argument is mandatory, meaning that the user must provide a value for it in order for the program to run. If the user does not provide a value for this argument, an error message will be displayed.
    parser.add_argument("--barcode", help="barcode to analyse", required=True)
    ##-f5 is a shorthand notation and --folder_fast5 is the long form.
    ###default="/data/sharedata/" sets the default value for the argument to "/data/sharedata/" if no value is specified.
    parser.add_argument("-f5","--folder_fast5", help="folder containing fast5 files",default="/data/sharedata/")
    ##nargs="?" indicates that this argument can take zero or one value.
    ##type=int specifies that the argument value should be converted to an integer.
    ##default="10" sets the default value for the argument to 10 if no value is specified.
    parser.add_argument("-t", "--threads", nargs="?", type = int, default="10", help="number of thres")
    #nargs="?" indicates that this argument can take zero or one value.
    #choices=['Plant_Bacteria_Funghi.nal',...] sets the list of allowed choices for the argument. The user must select one of these choices when running the program, otherwise an error will be raised. The allowed choices are a fasta file or a subset of the NCBI database.
    #default="Plant_Bacteria_Funghi.nal" sets the default value for the argument to "Plant_Bacteria_Funghi.nal" if no value is specified.
    parser.add_argument("-d", "--database", nargs="?", choices=['Plant_Bacteria_Funghi.nal', 'nt.nal', 'ITS_16S_18S_28S_LSU_SSU.nal',
                                                                'Metazoa_Organism_OR_metazoa_All_Fields_.nal','Xanthomonas_genomes.nal','curatedXylellaDatabase.nal'],
                        default="Plant_Bacteria_Funghi.nal",
                        help="database name; can be a fasta file or a subset of ncbi database; Default is Plant_Bacteria_Funghi")
    ##nargs="?" indicates that this argument can take zero or one value.
    ##default="/tmp/" sets the default value for the argument to "/tmp/" if no value is specified.
    parser.add_argument("-w", "--workdir", nargs="?", default="/tmp/")
    ##required=True specifies that the argument must be provided by the user when running the program. If no value is provided, an error will be raised.
    ##default="/data/sharedata/" sets the default value for the argument to "/data/sharedata/" if no value is specified.
    ##Note that there is a commented-out widget argument at the end, which suggests that this argument might be used in a graphical user interface (GUI) for selecting the folder interactively. 
    #However, it is commented.
    parser.add_argument("--folder_fastq", help='folder where fastq file are located', required=True, default="/data/sharedata/")#, widget='DirChooser')
    ##nargs="?" indicates that this argument can take zero or one value.
    ##default="" sets the default value for the argument to an empty string if no value is specified.
     ##required=True specifies that this argument is mandatory, meaning that the user must provide a value for it in order for the program to run. If the user does not provide a value for this argument, an error message will be displayed.
    parser.add_argument("-e", "--email", nargs="?", default="", help="email for blast search", required=True)
    ##this is to search the database directly on NCBI: ncbi search to retrieve GIs
    ##nargs="?" indicates that this argument can take zero or one value.
    ##default="" sets the default value for the argument to an empty string if no value is specified
    parser.add_argument("-s", "--search", nargs="?", default="", help="ncbi search to retrieve GIs")
    ##this specifies the minimum number of reads required to plot a species in the output figure as a fraction of the total mapped reads, which should be between 0 and 100. (barplot)
    parser.add_argument("-m", "--min",  default="10", help="minimum number of reads to plot a species as fraction of total mappd reads [0-100]")
    ##this specifies the "minimum identity to consider a valid alignment"
    #type=float specifies that the argument must be converted to a floating-point number.
    #default="90" sets the default value for the argument to "90" if no value is specified.
    parser.add_argument("-p", "--percentage", type = float, default="90", help="minimum identity to consider a valid alignment")
    ##this specifies the "minimum number of reads to sequence to consider a sample in the analysis"
    ##default="100" sets the default value for the argument to "100" if no value is specified.
    parser.add_argument("-mr", "--min_reads", type = float, default="100", help="minimum number of reads to sequence to consider a sample in the analysis")
    ##action='store_true' specifies that if the argument is present on the command line, its value should be set to True.
    ##default=False specifies that the default value of the argument is False if it is not present on the command line.
    parser.add_argument("-a", "--assemble", action='store_true', help="assembled-reads", default=False)
    ##action='store_true' specifies that if the argument is present on the command line, its value should be set to True.
    ##default=False specifies that the default value of the argument is False if it is not present on the command line.
    parser.add_argument("-ua", "--use_assembled", action='store_true', help="use assembled reads for discovery", default=False)
    parser.add_argument("-dh", "--database_history", action='store_true', help="", default=False)
    ##this allows you to correct the reads via racon
    parser.add_argument("-c", "--correct", action='store_true', help="correct or not reads via racon", default=False)
    parser.add_argument("-u", "--update", action='store_true', help="update database", default=False)#nargs="?", default="", required=True)
    #dest='verbose' specifies that the value of the argument will be stored in the variable named verbose.
    parser.add_argument('--verbose', help='be verbose', dest='verbose', action='store_true', default=False)
    ##Finally, you would call the parse_args() method on the ArgumentParser object to actually parse the command-line arguments and return an object containing the values of the options that were passed in.
    args = parser.parse_args()
    return args


def convert_color(test):
    values = dict()
    for index, letter in enumerate(string.ascii_letters):
        values[letter] = index + 1
    complete_colors = []
    for name_sp in test:
        if " " in name_sp:
            name_a = name_sp.split(" ")
            genera = name_a[0]
            del name_a[0]
            specie = "".join(name_a)
            name = name_sp
        else:
            genera = name_sp
            specie = "unknown"
            name  = name_sp + "unknown"
        genera = re.sub("[^a-z]+", "", genera.lower())
        specie = re.sub("[^a-z]+", "", specie.lower())
        name = re.sub("[^a-z]+", "", name.lower())
        genera = list(genera)
        specie = list(specie)
        name = list(name)
        rgb = [genera, name, specie]

        number_rgb = []
        for i in rgb:
            number  = sum([values[l] for l in i])
            random.seed(number)
            number_rgb.append(random.randint(0,255)/255)
        complete_colors.append([name_sp ,tuple(number_rgb)])
    return(complete_colors)

def blast(elm):
    blastn_cline = NcbiblastnCommandline(db=elm[1], query=elm[0], evalue=0.001, out=elm[2] ,outfmt = "6 qseqid sseqid bitscore sscinames pident evalue staxids qlen" )
    try:
        blastn_cline()
    except ValueError:
        print(blastn_cline)
    return(elm[2])

def common_start(*s):
   l = len(s)
   if l == 0:
       return None
   elif l == 1:
       return s[0]
   start = s[0]
   while start != "" and not all(ii.startswith(start) for ii in s[1:]):
       start = start[:-1]

   return start

def download_fasta(elm):
    search, retmax, retstart, email = elm
    Entrez.email = email
    handle = Entrez.esearch(db="nucleotide", term=search, retmax=retmax, retstart=retstart)
    record = Entrez.read(handle)
    records = str("\n".join(record["IdList"]))
    return (records)

def counpute_count(dict_match):
    genera_dict = {}
    species_dict = {}
    id, match_blast = dict_match
    if len(match_blast) > 1:
        substring_counts = {}
        length_array = int(len(match_blast))
        for i in range(0, length_array):
            for j in range(i + 1, length_array):
                string1 = match_blast[i]
                string2 = match_blast[j]
                match = SequenceMatcher(None, string1, string2).find_longest_match(0, len(string1), 0, len(string2))
                matching_substring = string1[match.a:match.a + match.size]
                if (matching_substring not in substring_counts):
                    substring_counts[matching_substring] = 1
                else:
                    substring_counts[matching_substring] += 1
        species = sorted(substring_counts.items(), key=lambda x: (x[1], x[0]), reverse=True)
        specie_id = True
        genera_id = True
        for organism in species:
            if " " in organism[0] and organism[0][0].isupper() and organism[0] in match_blast:
                fields = organism[0].split(" ")
                genera = fields[0]
                specie = organism[0]
                if genera != "" and specie != "" and specie_id:
                    if organism[0] in species_dict and specie_id:
                        specie_id = False
                        species_dict[organism[0]].append(id)
                        genera_dict[genera].append(id)
                        genera_id = False
                    else:
                        specie_id = False
                        species_dict[organism[0]] = [id]
                        genera_id = False
                        genera_dict[genera] = [id]
                elif genera != "" and specie == "" and genera_id:
                    if genera in genera_dict:
                        genera_dict[genera].append(id)
                        genera_id = False
                    else:
                        genera_id = False
                        genera_dict[genera] = [id]
    else:
        fields = match_blast[0].split(" ")
        genera = fields[0]
        species = match_blast[0]
        if genera in genera_dict:
            genera_dict[genera].append(id)
        else:
            genera_dict[genera] = [id]
        if species in species_dict:
            species_dict[species].append(id)
        else:
            species_dict[species] = [id]
    return([species_dict, genera_dict])

def racon(od, fastq, threads, fast5, cwd):
    for root, dirs, files in os.walk(os.path.join(cwd,fastq), topdown=False):
        reads_fastq = tempfile.NamedTemporaryFile(suffix=".fastq", delete=False)
        with open(reads_fastq.name, "w") as outfile:
            for filename in files:
                filename = os.path.join(cwd, fastq, filename)
                with gzip.open(filename, 'rb') as infile:
                    for line in infile:
                        outfile.write(line.decode())

    read_mapping=reads_fastq.name
    print("RUNNING PORECHOP")
    reads_porechop_1 = tempfile.NamedTemporaryFile(suffix=".fastq", delete=False)
    reads_porechop = "/tmp/" + uuid.uuid4().hex + ".fastq"
    porechop_cmd = PORECHOP % (read_mapping, threads, reads_porechop_1.name)
    ##sb.Popen serve per correre con python un programma di bash
    porechop = sb.Popen(porechop_cmd, shell=True, cwd=cwd, stderr=sb.PIPE, stdout=sb.PIPE)
    porechop.communicate()
    min = 1
    max = 200000
    with open(reads_porechop, "w") as fh:
        for record in SeqIO.parse(reads_porechop_1.name, "fastq"):
            if min < len(record.seq) < max:
                SeqIO.write(record, fh, "fastq")
    sam = tempfile.NamedTemporaryFile(suffix=".sam", delete=False)
    reads = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
    print("RUNNING MINIMAP")
    m = MINIMAP % (threads, reads_porechop, reads_porechop)
    print(m)
    minimap = sb.Popen(m, shell=True, cwd=cwd, stdout=sam, stderr=sb.PIPE)
    minimap.communicate()
    print("RUNNING RACON")
    r = RACON % (threads, reads_porechop, sam.name, reads_porechop)
    print(r)
    racon_cmd = sb.Popen(r, shell=True, cwd=cwd, stdout=reads, stderr=sb.PIPE)
    racon_cmd.communicate()
    sam = tempfile.NamedTemporaryFile(suffix=".sam")
    print("RUNNING MINIMAP")
    m = MINIMAP % (threads, reads.name, reads.name)
    print(m)
    minimap = sb.Popen(m, shell=True, cwd=cwd, stdout=sam, stderr=sb.PIPE)
    minimap.communicate()
    print("RUNNING RACON")
    r = RACON % (threads, reads.name, sam.name, reads.name)
    print(r)
    output = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
    jfc_out = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
    racon_cmd = sb.Popen(r, shell=True, cwd=cwd, stdout=output, stderr=sb.PIPE)
    racon_cmd.communicate()
    print("RUNNING JELLYFISH")
    jfc = JELLYFISH_COUNT % output.name
    jfd = JELLYFISH_DUMP % jfc_out.name
    print(jfc)
    print(jfd)
    jellyfishc = sb.Popen(jfc, shell=True, cwd=cwd, stdout=sb.PIPE)
    jellyfishd = sb.Popen(jfd, shell=True, cwd=cwd, stdin=jellyfishc.stdout)
    jellyfishd.communicate()
    count = 0
    kmer = "/tmp/" + uuid.uuid4().hex + ".fasta"
    with open(kmer, "w") as fh:
        for record in SeqIO.parse(jfc_out.name, "fasta"):
            if int(record.id) > 10 and len(record.seq) == 100:
                repetitive = 0
                while 10 >= repetitive:
                    repetitive += 1
                    count += 1
                    record.id = "kmer_" + str(count)
                    SeqIO.write(record, fh, "fasta")
    print(kmer)
    tmp_dir = tempfile.mkdtemp(dir="/tmp")
    print("RUNNING SPADES")
    spa = SPADES % (kmer, tmp_dir)
    print(spa)
    try:
        spade = sb.Popen(spa, shell=True, cwd=cwd, stderr=sb.PIPE, stdout=sb.PIPE)
        spade.communicate()
        assembled_contigs = os.path.join(tmp_dir, "contigs.fasta")
        bam = tempfile.NamedTemporaryFile(suffix=".bam", delete=False)
        print("RUNNING MINIMAP")
        m = MINIMAP_S % (threads, assembled_contigs, read_mapping)
        ss = SAMSORT % bam.name
        print(m)
        minimap = sb.Popen(m, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        samtools_sort = sb.Popen(ss, shell=True, cwd=cwd, stdin=minimap.stdout, stdout=sb.PIPE, stderr=sb.PIPE)
        samtools_sort.communicate()
        si = SAMINDEX % bam.name
        samtoos_index = sb.Popen(si, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        samtoos_index.communicate()
        ni = NANOPOLISHI % (os.path.join(fast5,fastq), reads_porechop)
        print(ni)
        nanopolish_index = sb.Popen(ni, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        nanopolish_index.communicate()
        regions = []
        print("RUNNING NANOPOLISH")
        dict_vcf_fasta = []
        for record in SeqIO.parse(assembled_contigs, "fasta"):
            region = record.id + ":0-" + str(len(record.seq))
            vcf = tempfile.NamedTemporaryFile(prefix="polished.", suffix=".vcf", delete=False)
            nv = NANOPOLISHV % (vcf.name, region, read_mapping, bam.name, assembled_contigs, threads)
            print(nv)
            nanopolish_var = sb.Popen(nv, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
            nanopolish_var.communicate()
            if os.path.isfile(vcf.name) and os.path.getsize(vcf.name) > 0:
                regions.append(vcf.name)
            else:
                dict_vcf_fasta.append(record)
        nva = NANOPOLISHVA % (assembled_contigs, " ".join(regions))
        print(nva)
        mlst_done = os.path.join(od, fastq + ".MLST.done.fasta")
        mlst = os.path.join(assembled_contigs + ".MLST.fasta")
        with open(mlst, "w") as fh:
            nanopolish_vcf = sb.Popen(nva, shell=True, cwd=cwd, stdout=fh, stderr=sb.PIPE)
        nanopolish_vcf.communicate()
        bam = tempfile.NamedTemporaryFile(suffix=".bam", delete=False)
        bi = BWAI % mlst
        bwa_index = sb.Popen(bi, shell=True, cwd="/tmp")
        bwa_index.communicate()
        bm = BWA % (threads, mlst, output.name)
        #shutil.copyfile(mlst, mlst_done) #os.path.join(cwd,output + ".fasta"))
        ss = SAMSORT % bam.name
        print(bm)
        bwa_mem = sb.Popen(bm, shell=True, cwd="/tmp", stdout=sb.PIPE, stderr=sb.PIPE)
        samtools_sort = sb.Popen(ss, shell=True, cwd=cwd, stdin=bwa_mem.stdout, stdout=sb.PIPE, stderr=sb.PIPE)
        samtools_sort.communicate()
        si = SAMINDEX % bam.name
        samtoos_index = sb.Popen(si, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        samtoos_index.communicate()
        vcf = tempfile.NamedTemporaryFile(suffix=".vcf", delete=False)
        bcfm = BCFTOOLS_MP % (mlst, bam.name) #FREEBAYES
        bcfc = BCFTOOLS_CALL % vcf.name
        print(bcfm)
        print(bcfc)
        bcftools_mpile = sb.Popen(bcfm, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        bcftools_call = sb.Popen(bcfc, shell=True, cwd=cwd, stdin=bcftools_mpile.stdout)
        bcftools_call.communicate()
        bcfi = BCFTOOLS_IN % vcf.name

        #BGZIP = "bgzip %s"
        #TABIX = "tabix -p vcf %s"
        bgzip = BGZIP % (vcf.name)
        vcf_gz = vcf.name + ".gz"
        tabix = TABIX % (vcf_gz)
        bgzip_run = sb.Popen(bgzip, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        bgzip_run.communicate()
        print(bgzip)
        tabix_run = sb.Popen(tabix, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        tabix_run.communicate()
        print(tabix)
        print(bcfi)
        bcftools_index = sb.Popen(bcfi, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        bcftools_index.communicate()
        bcfco = BCFTOOLS_CO % (mlst, mlst_done, vcf_gz)
        print(bcfco)
        #print(mlst_done)
        bcftools_call = sb.Popen(bcfco, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        bcftools_call.communicate()
        return(mlst_done, "A")
    except:
        return(output.name, "B")


def analysis():
    print("PROGRAM STARTED")
    args = setting()
    mem_threads = 0
    cwd = args.folder_fastq
    od = os.getcwd()
    Entrez.email = args.email
    #This is a Python code that uses the re module to substitute any character that is not a digit or a letter (uppercase or lowercase) with an underscore in the string contained in the variable "args.search".
    # re.sub: This is a method from the re module that substitutes a pattern in a string with a new string.
    #"[^0-9a-zA-Z]+": This is a regular expression pattern that matches any character that is not a digit (0-9) or a letter (a-zA-Z), one or more times (+). The caret (^) inside the square brackets means "not".
    #"_": This is the string that will replace any match of the regular expression pattern.
    #args.search: This is the string that will be searched and have its matches substituted.
    new_string = re.sub("[^0-9a-zA-Z]+", "_", args.search)
    name_plot = args.output + ".pdf"
    name_table = args.output + ".xlsx"
    figure_file = os.path.join(od, name_plot)
    name_table = os.path.join(od, name_table)
    os.environ["BLASTDB"] = "/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/tmp" ##have I this directory "/data2/blastdb"? I'm changing this with "/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI"
    blastdb = os.environ['BLASTDB']
    if args.search != "":
        #os.path.join(blastdb, new_string) is concatenating two path components - the variable blastdb and the variable new_string - 
        # using the appropriate path separator for the current operating system (/).
        #blastdb contains a path to a directory and new_string contains a string representing a filename so  the result will be a string representing the full path to the file
        database_name = os.path.join(blastdb, new_string)
        print("Using database %s" % database_name)
        database_name_check = os.path.join(blastdb, new_string + ".nal") ##it's adding .nal to database_name
        ##if a BLAST database with the name specified in the database_name_check (the same with .nal) does not exist or the args.update flag is set, 
        ##the code will proceed to download sequences from the NCBI database and create or update the BLAST database.
        if not os.path.exists(database_name_check) or args.update:
            #This line sends a query to the NCBI nucleotide database using the Entrez.esearch() function from the Biopython library. The query is specified by the args.search variable
            handle_n = Entrez.esearch(db="nucleotide", term=args.search)
            #This line reads the results of the query returned by the NCBI database and stores them in the record_n variable.
            record_n = Entrez.read(handle_n)
            ##retmax = 10000000: This line sets the maximum number of records to download in a single batch.
            retmax = 10000000
            #?
            retstart = 0
            #repeat = (str(int(record_n['Count'])/retmax)).split(".")[0] This line calculates the number of batches required to download all the records returned by the NCBI database.
            #The Count field in the record_n variable contains the total number of records that match the query. 
            #The repeat variable is calculated by dividing Count by retmax, converting the result to an integer, and rounding down to the nearest whole number.
            repeat = (str(int(record_n['Count'])/retmax)).split(".")[0]
            ##print FOUND and the number of the sequenced found with this specified database
            print("FOUND " + str(record_n['Count']) + " SEQUENCES\n")
            ##?????
            repeat = int(repeat) + 1
            #id_download = []: This line initializes an empty list that will be used to store the IDs of the records to be downloaded.
            id_download = []
            count = 0
            #while (count < repeat):: This line starts a loop that will run once for each batch of records to be downloaded.
            while (count < repeat):
                ##credo che con count +=1, ad ogni iterazione del loop, conta il numero di batch che abbiamo analizzato così che quando arriva all'ultimo batch, si blocca il loop (perchè count=repeat)
                count += 1
                #id_download.append([args.search,retmax, retstart, args.email]): This line appends a list of download parameters to the id_download list. 
                #Each list contains the search term (args.search), the maximum number of records to download (retmax), the starting record (retstart), and the email address of the user (args.email). 
                id_download.append([args.search,retmax, retstart, args.email])
                #The retstart variable is updated in each iteration of the loop to skip the records that have already been downloaded -> count: number of reiteration; retmax: numero di records per batch
                retstart = count * retmax
            ##result_list_tqdm = [] this line initializes an empty list
            result_list_tqdm = []
            ##if the number of batches required to download all the records returned by the NCBI database > 10, use 10 threads
            if repeat > 10:
                threads = 10
            ##if the batch are less then 10, the number of threads used is = to the number of the batch
            else:
                threads = repeat
            #with Pool(processes=threads) as pool:: This line initializes a multiprocessing pool with a number of worker processes equal to the threads variable.
            with Pool(processes=threads) as pool:
                ##for result in tqdm(pool.imap(func=download_fasta, iterable=id_download), total=len(id_download)):: This line starts a loop that downloads the sequences in each batch using the download_fasta() function and the parameters in the id_download list. 
                ##The pool.imap() function maps the download_fasta() function to the id_download list and distributes the batches to the worker processes in the multiprocessing pool. 
                ##The tqdm() function provides a progress bar for the downloads, and the total len of the progress bar is the len of the id list
                for result in tqdm(pool.imap(func=download_fasta, iterable=id_download), total=len(id_download)):
                    ##result_list_tqdm.append(result): This line appends the downloaded sequences to the result_list_tqdm list.
                    result_list_tqdm.append(result)
            ##the variable records join all the downloaded sequences in the list dividing them with a new line
            records = "\n".join(result_list_tqdm)
            #fp = tempfile.NamedTemporaryFile(dir="/tmp", delete= False): This line creates a named temporary file in the /tmp directory that will be used to store the downloaded sequences in FASTA format.
            fp = tempfile.NamedTemporaryFile(dir="/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/tmp", delete= False) ##here there is the directory dir="/tmp", TRY to change it with the same directory of the output?
            #fp.write(records.encode()): This line writes the downloaded sequences to the temporary file in binary format. ???
            fp.write(records.encode())
            ##creates a string cm that runs the blastdb_aliastool command with the following arguments (it creates the .n.gil database and an alias file with the .nal extention):
            #-db nt: sets the name of the database to "nt" (corrispond to all the sequences in the database nucleotide on NCBI) we change it with a smaller database that we have dowloaded previously 
            #-dbtype nucl: specifies that the database type is nucleotide
            #-gilist <path to tempfile>: specifies the path to a file containing a list of GI numbers to include in the database
            #-out <database name>: sets the output database name 
            ##--> l'output is the database ".n.gil" comprexed and an alias file with the same name but with the extension .nal
            cm = "/usr/local/bin/blastdb_aliastool -db /Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/tmp/xylella_fastidiosa_database.fasta -dbtype nucl -gilist %s -out %s" % (fp.name, database_name)
            print(cm)
            #Then the Popen() function from the subprocess module is used to run the cm command in a new process with a shell. 
            #Finally, the communicate() method is called on the blastalias object to wait for the process to finish and collect any output.
            blastalias = sb.Popen(cm, shell= True)
            blastalias.communicate()
    
    elif args.database != "":
        fasta_all = []
        if os.path.isabs(args.database):
            if args.database.endswith("fasta"):
                for rec in SeqIO.parse(args.database, "fasta"):
                    desc = rec.description
                    if 200 < len(rec.seq) < 10000:
                        single_elem = desc.split(" ")
                        if not "sp." in single_elem[2]:
                            if len(single_elem) >= 3:
                                species = " ".join([single_elem[1], single_elem[2]])
                                rec.description = species
                                fasta_all.append(rec)
                database_name = os.path.join(blastdb,args.database + "clean.fasta")
                SeqIO.write(fasta_all, database_name , "fasta")
                print("Using database %s" % database_name)

            elif args.database.endswith(".nal") or args.database.endswith(".Nal"):
                database_name = os.path.abspath(args.database)
                interm = database_name.split(".")[:-1]
                database_name = ".".join(interm)
                print("Using database %s" % database_name)
            else:
                print("DATABASE NOT FOUND")
        elif os.path.isfile(os.path.abspath(args.database)):
            database_name_orig = os.path.abspath(args.database)
            if database_name_orig.endswith("fasta"):
                for rec in SeqIO.parse(database_name_orig, "fasta"):
                    desc = rec.description
                    if 200 < len(rec.seq) < 10000:
                        single_elem = desc.split(" ")
                        if not "sp." in single_elem[2]:
                            if len(single_elem) >= 3:
                                species = " ".join([single_elem[1], single_elem[2]])
                                rec.description = species
                                fasta_all.append(rec)
                database_name = os.path.join(blastdb,args.database + "clean.fasta")
                print("Using database %s" % database_name)
                SeqIO.write(fasta_all, database_name , "fasta")
            elif database_name_orig.endswith(".nal") or database_name_orig.endswith(".Nal"):
                interm = database_name_orig.split(".")[:-1]
                database_name = ".".join(interm)
            else:
                print("DATABASE NOT FOUND")
        elif os.path.isfile(os.path.join(blastdb, args.database)):
            database_name_orig = os.path.join(blastdb, args.database)
            if database_name_orig.endswith("fasta"):
                for rec in SeqIO.parse(database_name_orig, "fasta"):
                    desc = rec.description
                    if 200 < len(rec.seq) < 10000:
                        single_elem = desc.split(" ")
                        if not "sp." in single_elem[2]:
                            if len(single_elem) >= 3:
                                species = " ".join([single_elem[1], single_elem[2]])
                                rec.description = species
                                fasta_all.append(rec)
                database_name = os.path.join(blastdb, args.database + "clean.fasta")
                print("Using database %s" % database_name)
                SeqIO.write(fasta_all, database_name, "fasta")
            elif database_name_orig.endswith(".nal") or database_name_orig.endswith(".Nal"):
                interm = database_name_orig.split(".")[:-1]
                database_name = ".".join(interm)
                print("Using database %s" % database_name)
            else:
                print("DATABASE NOT FOUND")
        else:
            print("DATABASE NOT FOUND")
    else:
        database_name = os.path.join(blastdb, "nt")
        print("Using database %s" % database_name)
    if "," in args.barcode:
        list_barcode = args.barcode.split(",")
    else:
        list_barcode = [args.barcode]
    dict_species_all = {}
    dict_table_all = {}
    number_species = []
    files_kt_list = []
    for barcode in list_barcode:
        ##For each barcode in list_barcode, the code constructs a path to the directory containing the barcode data (dir_fastq). 
        #If this directory exists, the code sets the result_blast flag to False, prints a message indicating which barcode is being processed, and initializes two empty lists (fastx_all and fastx_all_assembled) for storing data later in the code.
        dir_fastq = os.path.join(cwd, barcode) ###cwd= args.folder_fastq (la cartella dove c'è il documento fastq da analizzare)
        if os.path.exists(dir_fastq):
            result_blast = False
            print("EXECUTING BARCODE: %s" % barcode)
            fastx_all = []
            fastx_all_assembled = []
            ##IF I DECIDE TO MAKE THE ASSEMBLE
            if args.assemble or args.use_assembled:
                if args.folder_fast5 == "":
                    print("ASSEMBLY NOT DONE, PLEASE ADD FAST5 PARAMETER - RUNNING ALIGNMENT")
                    if args.use_assembled:
                        sys.exit("FAST5 NOT PASSED WHILE ASKED FOR ASSEMBLED READS FOR ALIGNEMENT")
                else:
                    barcode_corr, spade_t = racon(od, barcode, args.threads, args.folder_fast5, cwd)
                    for record in SeqIO.parse(barcode_corr, "fasta"):
                        fp = tempfile.NamedTemporaryFile(dir="/tmp", suffix=".fasta", mode="w", delete=False)
                        fo = tempfile.NamedTemporaryFile(dir="/tmp", mode="w", prefix=barcode, suffix=".blastn", delete=False)
                        SeqIO.write(record, fp, "fasta")
                        fastx_all_assembled.append([fp.name, database_name, fo.name])
                        # else:
                        #     print(barcode + " BARCODE NOT DONE WITH CONSENSUS SEQUENCES")
                        #     continue
            ##IF I DECIDE TO USE THE ASSEMBLED READS
            if args.use_assembled:
                if spade_t == "A":
                    print("USING ASSEMBLED READS FOR ALIGNEMENT TO DATABASE FOR BARCODE " + barcode)
                else:
                    print("USING NON-ASSEMBLED READS FOR ALIGNEMENT TO DATABASE FOR BARCODE " + barcode)
                fastx_all = fastx_all_assembled
            ##IF I DON'T DECIDE TO USE THE USED THE ASSEMBLED READS
            else:
                ##porechop is used for adapter trimming and demultiplexing of Oxford Nanopore Technologies (ONT) sequencing data
                print("RUNNING PORECHOP")
                reads_porechop = tempfile.NamedTemporaryFile(suffix=".fastq", delete=False)
                print(reads_porechop)

                porechop_cmd = PORECHOP % (barcode, args.threads, reads_porechop.name)
                print(porechop_cmd)
                print ("reads_porechop.name:" + reads_porechop.name)  #VALERIA PER CONTROLLO
                print("porechop command: " + porechop_cmd) #DANI PER CONTROLLO (?) ##this is the folder where it is stored /var/folders/c2/gghtq5tn77z8g_cgbf14n7nr0000gn/T/
                porechop = sb.Popen(porechop_cmd, shell=True, cwd=cwd, stderr=sb.PIPE, stdout=sb.PIPE) ###cwd= args.folder_fastq (la cartella dove c'è il documento fastq da analizzare)
                porechop.communicate()
                ##IF I DECIDE TO CORRECT THE READS
                if args.correct:
                    sam = tempfile.NamedTemporaryFile(suffix=".sam", delete=False)
                    reads = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
                    print("RUNNING MINIMAP")
                    m = MINIMAP % (args.threads, reads_porechop.name, reads_porechop.name)
                    print(m)
                    minimap = sb.Popen(m, shell=True, cwd=cwd, stdout=sam, stderr=sb.PIPE)
                    minimap.communicate()
                    print("RUNNING RACON")
                    r = RACON % (args.threads, reads_porechop.name, sam.name, reads_porechop.name)
                    print(r)
                    racon_cmd = sb.Popen(r, shell=True, cwd=cwd, stdout=reads, stderr=sb.PIPE)
                    racon_cmd.communicate()
                    if int(os.path.getsize(reads.name)) > 0:
                        sam = tempfile.NamedTemporaryFile(suffix=".sam")
                        print("RUNNING MINIMAP")
                        m = MINIMAP % (args.threads, reads.name, reads.name)
                        minimap = sb.Popen(m, shell=True, cwd=cwd, stdout=sam, stderr=sb.PIPE)
                        minimap.communicate()
                        print("RUNNING RACON")
                        r = RACON % (args.threads, reads.name, sam.name, reads.name)
                        output = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
                        racon_cmd = sb.Popen(r, shell=True, cwd=cwd, stdout=output, stderr=sb.PIPE)
                        racon_cmd.communicate()
                        count = 0
                        for record in SeqIO.parse(output.name, "fasta"):
                            count += 1
                            fp = tempfile.NamedTemporaryFile(dir="/tmp", suffix=".fasta", mode="w", delete=False)
                            fo = tempfile.NamedTemporaryFile(dir="/tmp", mode="w", prefix=barcode, suffix=".blastn", delete=False)
                            SeqIO.write(record, fp, "fasta")
                            fastx_all.append([fp.name, database_name,fo.name])
                    else:
                        count = 0
                        print(reads_porechop.name)
                        for record in SeqIO.parse(reads_porechop.name, "fastq"):
                            count += 1
                            fp = tempfile.NamedTemporaryFile(dir="/tmp", suffix=".fasta", mode="w", delete=False)
                            fo = tempfile.NamedTemporaryFile(dir="/tmp", mode="w", prefix=barcode, suffix=".blastn", delete=False)
                            SeqIO.write(record, fp, "fasta")
                            fastx_all.append([fp.name, database_name, fo.name])
                            if not barcode.endswith("raw"):
                                barcode = barcode + "raw"
                ##IF I DON'T DECIDE TO CORRECT THE READS
                else:
                    print(reads_porechop.name)
                    ##The following: SeqIO module to read a Fastq file, perform some processing on each record, and add the processed record to a list.
                    # The for loop iterates over the records in the Fastq file reads_porechop.name (the file with our sequences with the adapter removed and demultiplexed) using the SeqIO.parse() function with the file format specified as "fastq".
                    #reads_porechop= (the file with our sequences with the adapter removed and demultiplexed) in a temporary file format with the extension (.fastq)
                    for record in SeqIO.parse(reads_porechop.name, "fastq"):
                        #print("records: " + record) ##DANI PER CONTROLLO
                        #For each record, the code creates a temporary file (fp) to write the sequence demultiplexed and without adapter in fasta format.
                        fp = tempfile.NamedTemporaryFile(dir="/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/fp_fo", suffix=".fasta", mode="w", delete=False) ##here i changed the directory from "/tmp" to "/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/fp_fo"
                        #It also creates another temporary file (fo) with a name derived from the barcode (presumably for writing the BLAST output).
                        fo = tempfile.NamedTemporaryFile(dir="/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/fp_fo", mode="w", prefix=barcode, suffix=".blastn", delete=False) ##here i changed the directory from "/tmp" to "/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/fp_fo"
                        #The code then uses the SeqIO.write() function to write the record to the fp file in fasta format. 
                        SeqIO.write(record, fp, "fasta")
                        #It then appends a list [fp.name, database_name, fo.name] to the fastx_all list, 
                        #which contains the name of the fasta file, the name of the BLAST database to search against (database_name), and the name of the BLAST output file.
                        fastx_all.append([fp.name, database_name,fo.name]) #esempio del contenuto di questa lista: "[['/tmp/tmposk0hxgn.fasta', '/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/tmp/Xylella_NOT_uncultured_All_Fields_AND_100_5000_SLEN_AND_subsp_All_Fields_', '/tmp/barcode18p4h_1dy8.blastn']]"
                        if not barcode.endswith("raw"):
                            if barcode.endswith("/"):
                                barcode = barcode.split("/")[-2]+ "raw"
                            else:
                                barcode = barcode.split("/")[-1]+ "raw"
                if len(fastx_all) < int(args.min_reads):
                    print("Not analysing barcode " + barcode +". Less than 50 reads")
                    continue
            result_list = []
            print("RUNNING BLAST USING "  + str(args.threads) + " THREADS")
            ##This code block is responsible for running BLAST on the demultiplexed reads for each barcode.
            #First, an empty list is created to store the results (result_list). 
            # #Then, the code enters a with block where a multiprocessing.Pool object is created with a number of processes equal to args.threads.
            with Pool(processes=args.threads) as pool:
                #The Pool.imap() method is used to apply the blast() function (the luigi fuction that ask for tree element to be specified)  to each element in the fastx_all iterable (i.e., the list of input files for BLAST, composed by: the fil fp with the reads demulptiplexed and without adapter in a fasta format; the name of the database.n.gil without the extention, and the output name). 
                # ---> the function blast is the one created by Luigi: blast(elm): blastn_cline = NcbiblastnCommandline(db=elm[1], query=elm[0], evalue=0.001, out=elm[2] ,outfmt = "6 qseqid sseqid bitscore sscinames pident evalue staxids qlen" )
                #tqdm library is used to display a progress bar; the total lenght of the bar is the lenght of the fastx_all file.
                for result in tqdm(pool.imap(func=blast, iterable=fastx_all), total=len(fastx_all)): ##in questo caso stiamo specificando che elm (l'argomento della funzione blust creata da luigi) è fastx_all!!! che è una lista formata da tre elementi (ecco perchè elm [0], elm [1] etc)
                    #The results (files .blstn) are stored in the result_list variable
                    result_list.append(result)
            result_list_tqdm = []
            ##od = os.getcwd() --> The os module provides a way to interact with the operating system in Python. The getcwd() function from the os module returns a string representing the current working directory. The code od = os.getcwd() assigns the value returned by getcwd() to the variable od.
            #we are defining the name of the output file of blust for each of the analized barcode 
            name_blast = os.path.join(od, args.output + barcode + ".blast.txt")
            #write a new file that will be named with the name_blust defined before
            with open(name_blast, "w") as new_file:
                ##for each of the result name (the ones .blstn) open a list name data_single_file
                for name in result_list:
                    data_single_file = []
                    ##For each result name in the list, it opens the file using the with statement and assigns it to the variable f
                    with open(name) as f:
                        #It then loops through each line in the file using a for loop, writes the line to the new_file, and appends it to the data_single_file list.
                        for line in f:
                            new_file.write(line)
                            data_single_file.append(line)
                        #After reading all lines in the file, it writes a newline character to new_file.
                        new_file.write("\n")
                    #and then joins the data_single_file lists into a single string using the newline character as a separator and assigns it to the variable blast_out_single. 
                    blast_out_signle = "\n".join(data_single_file)
                    ##The blast_out_single string is then appended to the result_list_tqdm list.
                    result_list_tqdm.append(blast_out_signle)
                    #Finally, the code uses the os module's remove() function to delete the file that was just processed.
                    os.remove(name)
                    #print(result_list_tqdm) #VALERIA PER CONTROLLO
                    #### ^ The purpose of this code is to merge the contents of multiple files into a single file and to store the individual file contents as a list of strings in result_list_tqdm. The with statement is used to ensure that each file is closed after it has been read and written to new_file, and the os.remove() function is used to delete the original file after it has been processed to save space on the disk.
            #The loop over result_list_tqdm check if any lines of data are present in the result_list_tqdm list and if result_blast is not true (before there is written that, if the path to the fastq file exist the result_blast should be false)
            for line in result_list_tqdm:
                if line != "" and not result_blast:
                    result_blast = True
            #If no data is present (if result_blust is not true, like we said in the previous line) the loop continues
            if not result_blast:
                continue
            read_best_hit = {}
            ##the code below parses each line of output from a sequence alignment file, then sorting and selecting the highest-scoring alignment(s) for each read based on certain criteria. The selected alignment(s) are stored in a dictionary called read_best_hit.
            #The loop over result_list_tqdm processes each line of the file by splitting it into a list of elements using tab ('\t') as a separator.
            for output in result_list_tqdm:
                if output != "":
                    align_species = {}
                    #each element of the reslut_list_tqdm is split with a \n and this new format is stored in the variable align
                    align = output.split("\n") 
                    #each element in align is split with a \t
                    for single_align in align:
                        align_elm = single_align.split("\t")
                        #If the length of the single_align list is greater than 3, the program proceeds to the next step. Otherwise, it ignores the line.
                        if len(single_align) > 3:
                            align_species_score = {}
                            read = align_elm[0]
                            score = float(align_elm[2])
                            align_species_score[align_elm[4]] =  [align_elm[3]]
                            if score > 200 and float(align_elm[4]) >= args.percentage :
                                if score in align_species:
                                    for key in  align_species[score]:
                                        if align_elm[4] in align_species_score:
                                            align_species[score][key].append(align_elm[3])
                                        else:
                                            align_species[score][key] = [align_elm[3]]
                                else:
                                    align_species[score] = align_species_score
                    if align_species:
                        list_align_species = list(align_species.items())
                        list_align_species.sort(reverse=True)
                        align_species_ident = list_align_species[0][1]
                        list_align_species_b = list(align_species_ident.items())
                        list_align_species_b.sort(reverse=True)
                        align_species_ident_b = list_align_species_b[0][1]
                        read_best_hit[read] = align_species_ident_b
            dict_match = []
            for match in read_best_hit:
                dict_match.append([match, read_best_hit[match]])
            result_list_tqdm = []
            with Pool(processes=args.threads) as pool:
                for result in tqdm(pool.imap(func=counpute_count, iterable=dict_match), total=len(dict_match)):
                    result_list_tqdm.append(result)

            kt_barcode = tempfile.NamedTemporaryFile(suffix=".txt", prefix=barcode, delete=False, mode = "w")
            files_kt_list.append(kt_barcode.name)
            for value in result_list_tqdm:
                for key in value[0]:
                    read = value[0][key][0]
                    taxid = ncbi.get_name_translator([key])
                    #print("before")
                    #print(taxid)
                    #print("after")
                    if bool(taxid):
                        try:
                            #print(str(taxid[key][0]))
                            species_count = str(taxid[key][0]).split(";")[0]
                            kt_barcode.write(read + "\t" + species_count + "\n")
                        except:
                            print("NO TAXID FOR " + taxid[key][0])
                            continue
            species_dict = {}
            genera_dict = {}
            for result in result_list_tqdm:
                sp = list(result[0].keys())
                gen = list(result[1].keys())
                if len(sp) > 0:
                    if sp[0] in species_dict:
                        species_dict[sp[0]] = species_dict[sp[0]]+ 1
                    else:
                        species_dict[sp[0]] = 1
                if len(gen) > 0:
                    if gen[0] in genera_dict:
                        genera_dict[gen[0]] = [genera_dict[gen[0]][0] + 1]
                    else:
                        genera_dict[gen[0]] = [1]
            total_reads_mapped = sum(species_dict.values())
            min_value = total_reads_mapped * (float(args.min)/100)
            if min_value < 1 :
                if not args.assemble:
                    min_value = 1
            print("minimum number of reads used " + str(min_value))
            species_retain = {}
            specied_final = {}
            for key in species_dict:
                current_value = species_dict[key]
                if current_value >= min_value:#float(args.min):
                    species_retain[key] = species_dict[key]
            for key in species_retain:
                specied_final[key] = species_retain[key] / sum(species_retain.values()) * 100
                number_species.append(key)
            dict_species_all[barcode] = specied_final
            dict_table_all[barcode] = dict(sorted(species_dict.items(), key=lambda item : item[1]))
    # KTIMPORTTAX
        else:
            print("PATH FOR BARCODE " + barcode + " DO NOT EXISTS. CHECK BARCODE NAME OR PATH" )
    if len(files_kt_list) > 0:
        files_kt = " ".join(files_kt_list)
        k = KTIMPORTTAX % (blastdb, args.output + ".html", files_kt)
        print(k)
        ktimport = sb.Popen(k, shell=True, cwd=od)
        ktimport.communicate()

        species_pd = pd.DataFrame.from_dict(dict_species_all, orient='index')
        species_name_colors = [name for barcode in dict_species_all for name in dict_species_all[barcode]]
        species_name_colors = list(dict.fromkeys(species_name_colors))
        cmap = convert_color(species_name_colors)
        cmap2 = []
        for i in cmap:
            cmap2.append(i[1])
        species_abs = pd.DataFrame.from_dict(dict_table_all, orient='index')
        species_abs.transpose().to_excel(name_table)
        cmap1 = LinearSegmentedColormap.from_list("my_colormap", cmap2)
        species_pd.head()
        species_pd.fillna(0)
        f = plt.figure()
        species_pd.plot(kind="bar", stacked=True, colormap=cmap1, ax=f.gca())
        plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
        plt.title("Abundance species")
        plt.xlabel("barcode")
        plt.ylabel("Relative abundance(%)")
        #plt.show()
        plt.savefig(figure_file)
    else:
        print("NOTHING TO PLOT")




if __name__ == '__main__':
    analysis()