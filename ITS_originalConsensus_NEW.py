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
#The "-o" option is used to specify the output file for the annotated tree. This option accepts the path and filename of the output file, but here we specified only the name;
#the final %s is to specified the input.file 
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
##Racon is a software tool for polishing (migliorare) genome assemblies using long read sequencing data. Its main function is to correct errors in an initial draft assembly using long-read sequencing data from PacBio or Oxford Nanopore Technologies (ONT) platforms.
##First, Racon maps the long reads back to the initial draft assembly (I don't know if has to be an assembly) using a fast and efficient alignment algorithm. Next, it identifies discrepancies between the draft assembly and the long reads, such as mismatches, insertions, and deletions. Finally, 
##Racon generates a new consensus sequence for each region of the assembly, based on the majority vote of the aligned long reads.
#"t": threads
#-f, --fragment-correction perform fragment correction instead of contig polishing (overlaps file should contain dual/self overlaps!)
#?? "-f" specifies the input file path for the FASTQ or FASTA format reads to be used for polishing the assembly. This option tells Racon which long read sequencing data to use to correct errors in the initial draft assembly.
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
    ##there is a file blastn for each read
    blastn_cline = NcbiblastnCommandline(db=elm[1], query=elm[0], evalue=0.001, out=elm[2] ,outfmt = "6 qseqid sseqid bitscore sscinames pident evalue staxids qlen" ) ##qseqid: Query sequence ID. sseqid: Subject sequence ID. bitscore: Bit score of the alignment. sscinames: Scientific name of the subject organism. pident: Percent identity of the alignment. evalue: Expect value of the alignment. staxids: Taxonomy ID of the subject organism. qlen: Length of the query sequence.
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

#This function takes a dictionary dict_match as input, which contains the id and the Blast hits of a sequence.
def counpute_count(dict_match):
    genera_dict = {}
    species_dict = {}
    id, match_blast = dict_match
    #The function then checks if there are multiple Blast hits for the given sequence. 
    if len(match_blast) > 1:
        #If there are, it initializes an empty substring_counts dictionary to keep track of the counts of matching substrings between any two Blast hits. 
        substring_counts = {}
        length_array = int(len(match_blast))
        ##It then loops over all pairs of Blast hits and uses the SequenceMatcher function to find the longest common substring between them. 
        for i in range(0, length_array):
            for j in range(i + 1, length_array):
                string1 = match_blast[i]
                string2 = match_blast[j]
                ##SequenceMatcher is a Python standard library module that provides a class for comparing sequences. It compares pairs of input sequences (e.g. strings, lists) and returns a similarity ratio between 0 and 1, where 0 means completely different sequences and 1 means identical sequences.
                ##The find_longest_match method of SequenceMatcher finds the longest contiguous matching subsequence between two strings. --> It takes four arguments, the starting and ending indices of the first sequence and the starting and ending indices of the second sequence.
                #It returns a named tuple (the match object) containing the starting indices of the longest matching subsequence in each of the two sequences, as well as the length of the matching subsequence.
                match = SequenceMatcher(None, string1, string2).find_longest_match(0, len(string1), 0, len(string2))
                ##The match object returned by find_longest_match() contains information about the matching subsequence, including its starting position in string1 (match.a), its length (match.size), and its starting position in string2 (match.b).
                #We extract the matching substring from string1 using slicing: string1[match.a:match.a + match.size].
                #Finally, we store the matching substring in the matching_substring variable.
                matching_substring = string1[match.a:match.a + match.size]
                if (matching_substring not in substring_counts):
                    #The substring_counts dictionary is then updated with the counts of each matching substring.
                    substring_counts[matching_substring] = 1
                else:
                    substring_counts[matching_substring] += 1
        ##The function then sorts the substring_counts dictionary by the count of matching substrings in descending order, and then by the matching substring in ascending order.
        species = sorted(substring_counts.items(), key=lambda x: (x[1], x[0]), reverse=True)
        specie_id = True
        genera_id = True
        #It then loops over each matching substring in the sorted list and checks if it is a species name or a genus name. 
        for organism in species:
            #If the matching substring contains a space and starts with an uppercase letter (a convention for species names), it is considered a species name. 
            if " " in organism[0] and organism[0][0].isupper() and organism[0] in match_blast:
                #The function then extracts the genus and species names from the matching substring and checks if they are empty or not.
                fields = organism[0].split(" ")
                genera = fields[0]
                specie = organism[0]
                #If they are not empty, the function adds the sequence id to the corresponding lists in the species_dict and genera_dict dictionaries.
                if genera != "" and specie != "" and specie_id:
                    # If the genus or species name has already been added to the dictionary, the function does not add it again.
                    #If the matching substring contains a space but does not start with an uppercase letter, it is considered a genus name, and the function follows the same logic as for a species name but only updates the genera_dict dictionary.
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
                ##If the matching substring does not contain a space, it is considered a genus name, and the function adds the sequence id to the corresponding list in the genera_dict dictionary. 
                ##If the genus name has already been added to the dictionary, the function does not add it again.
                elif genera != "" and specie == "" and genera_id:
                    if genera in genera_dict:
                        genera_dict[genera].append(id)
                        genera_id = False
                    else:
                        genera_id = False
                        genera_dict[genera] = [id]
    ###If there is only one Blast hit for the given sequence, the function extracts the genus and species names from the Blast hit and adds the sequence id to the corresponding lists in the species_dict and genera_dict dictionaries.
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
    ###The function then returns a list containing the species_dict and genera_dict dictionaries.!
    return([species_dict, genera_dict])


##od = os.getcwd() --> The os module provides a way to interact with the operating system in Python. The getcwd() function from the os module returns a string representing the current working directory. The code od = os.getcwd() assigns the value returned by getcwd() to the variable od.
def racon(od, fastq, threads, fast5, cwd):
    #the function uses os.walk() to traverse a directory tree rooted at os.path.join(cwd,fastq) and process each FastQ file found.
    #os.walk() returns a generator that yields a 3-tuple (dirpath, dirnames, filenames) for each directory in the tree. 
    #The topdown argument is set to False so that the directory tree is traversed bottom-up.
    for root, dirs, files in os.walk(os.path.join(cwd,fastq), topdown=False):
        reads_fastq = tempfile.NamedTemporaryFile(suffix=".fastq", delete=False)
        #the function opens the temporary file in write mode
        with open(reads_fastq.name, "w") as outfile:
            for filename in files:
                #The filename variable is assigned the absolute path of the FastQ file by joining the cwd, fastq, and filename using the os.path.join() function.
                filename = os.path.join(cwd, fastq, filename)
                #For each file, the function opens the file using gzip.open() (assuming it's a gzipped FastQ file) and reads its contents line by line in binary mode ('rb').
                with gzip.open(filename, 'rb') as infile:
                    for line in infile:
                        # It then writes each line to a temporary file (reads_fastq) after decoding it from bytes to string using the decode() method.
                        # This is necessary because the contents of the FastQ file are in binary format, so they need to be converted to strings before they can be processed
                        #This process is repeated for each file in the directory tree, and the resulting temporary file contains all the FastQ reads from all the files in the directory tree.
                        outfile.write(line.decode())
    #reads_fastq is the file written above (fixed in the variable colled outfile) with the line of the fastq file decoded
    ##this file is stored in the variable read_mapping
    #This variable is later used as an argument to the PORECHOP command (a command-line tool for trimming adapters from Oxford Nanopore reads) to specify the input FastQ file.
    read_mapping=reads_fastq.name
    print("RUNNING PORECHOP")
    #reads_porechop_1 is a temporary file with a .fastq suffix that will be used to store the trimmed reads output by Porechop,
    reads_porechop_1 = tempfile.NamedTemporaryFile(suffix=".fastq", delete=False)
    #reads_porechop is a file path string that specifies the name and location of the final output FastQ file.
    #it generates a random file name for the final output FastQ file.
    #uuid.uuid4().hex generates a random UUID (Universally Unique IDentifier) version 4, which is a 128-bit value represented as a hexadecimal string.
    #This UUID is then concatenated with ".fastq" to create a file name with a .fastq suffix.
    reads_porechop = "/tmp/" + uuid.uuid4().hex + ".fastq" ##change /tmp/ with /Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/RACON
    #porechop_cmd string using the PORECHOP constant and the variables read_mapping, threads, and reads_porechop_1.name. The % operator is used to substitute these variables into the PORECHOP string (-input -threads - output)
    porechop_cmd = PORECHOP % (read_mapping, threads, reads_porechop_1.name)
    ##sb.Popen serve per correre con python un programma di bash (commento di luigi?)
    #The PORECHOP command is executed using sb.Popen(). The shell=True argument indicates that the command should be run in a shell environment. 
    #The cwd argument specifies the working directory where the command should be run, and the stderr and stdout arguments specify that the output and error messages should be captured in pipes
    porechop = sb.Popen(porechop_cmd, shell=True, cwd=cwd, stderr=sb.PIPE, stdout=sb.PIPE)
    #The porechop.communicate() method is called to start the porechop program and wait for it to finish. This method blocks the execution of the script until the command finishes running.
    porechop.communicate()
    min = 1
    max = 200000
    #then it open the files reads_porechop with a random UUID in a write mode and storing it in the variable fh
    with open(reads_porechop, "w") as fh:
        #the code reads from the intermediate output file (reads_porechop_1) in "fastq" format using SeqIO.parse() method
        for record in SeqIO.parse(reads_porechop_1.name, "fastq"):
            #For each record in the intermediate output file, if the length of the read sequence is between min and max (inclusive)
            if min < len(record.seq) < max:
                ##the record is written to the output file (reads_porechop) using SeqIO.write() method in "fastq" format
                SeqIO.write(record, fh, "fastq")
    sam = tempfile.NamedTemporaryFile(suffix=".sam", delete=False)
    reads = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
    print("RUNNING MINIMAP")
    #After trimming the reads and saving the output, the code runs Minimap2 to align the trimmed reads against each other (it is specified also because the second and the third element in the bracket are both "reads_porechop")
    m = MINIMAP % (threads, reads_porechop, reads_porechop)
    print(m)
    #The minimap variable is a subprocess created using subprocess.Popen() method to execute the command specified by MINIMAP.
    #the output is stored in the SAM file called "sam" created before
    minimap = sb.Popen(m, shell=True, cwd=cwd, stdout=sam, stderr=sb.PIPE)
    minimap.communicate()
    print("RUNNING RACON")
    #The code then runs Racon to polish the read alignments generated by Minimap2 (in this case we don't use a reference sequence (like a genome) to polish the reads but we align the reads against each other; in this case it's used like a read error-correction tool) 
    ##Racon do it, based on the majority vote of the aligned long reads
    r = RACON % (threads, reads_porechop, sam.name, reads_porechop)
    print(r)
    #The racon_cmd variable is a subprocess created using subprocess.Popen() method to execute the command specified by RACON
    ##The polished reads are saved in a Temporary file named reads, that with create before, which is in FASTA format.
    racon_cmd = sb.Popen(r, shell=True, cwd=cwd, stdout=reads, stderr=sb.PIPE)
    racon_cmd.communicate()
    ###The code then runs Minimap2 and Racon again to further polish the reads, and saves the polished reads in a file named output.
    sam = tempfile.NamedTemporaryFile(suffix=".sam")
    print("RUNNING MINIMAP")
    m = MINIMAP % (threads, reads.name, reads.name)
    print(m)
    minimap = sb.Popen(m, shell=True, cwd=cwd, stdout=sam, stderr=sb.PIPE)
    minimap.communicate()
    print("RUNNING RACON")
    r = RACON % (threads, reads.name, sam.name, reads.name)
    print(r)
    ##here it create the output file of rancon called "output" in a fasta format
    output = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
    #here it create the output file of jellifish called "jfc_out" in a fasta format
    jfc_out = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
    racon_cmd = sb.Popen(r, shell=True, cwd=cwd, stdout=output, stderr=sb.PIPE)
    racon_cmd.communicate()
    print("RUNNING JELLYFISH")
    ##The code then runs Jellyfish to count the k-mers in the polished reads
    ##K-mers are substrings of length k that are extracted from sequencing reads
    ##It's used for estimating genome size, assessing sequencing quality, or identifying repeat regions in a genome.
    ##It outputs its k-mer counts in a binary format (hash table), which can be translated into a human-readable text format using the "jellyfish dump" command
    jfc = JELLYFISH_COUNT % output.name
    jfd = JELLYFISH_DUMP % jfc_out.name
    print(jfc)
    print(jfd)
    #jellyfishc takes the output of JELLYFISH_COUNT and pipes it to jellyfishd, which converts the binary output to text. 
    jellyfishc = sb.Popen(jfc, shell=True, cwd=cwd, stdout=sb.PIPE)
    jellyfishd = sb.Popen(jfd, shell=True, cwd=cwd, stdin=jellyfishc.stdout)
    jellyfishd.communicate()
    count = 0
    ##Finally, the code creates a temporary file named kmer with a random file name in the /tmp/ directory to save the k-mers.
    kmer = "/tmp/" + uuid.uuid4().hex + ".fasta" #change /tmp/ with /Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/RACON
    #it opens this file in write mode, 
    with open(kmer, "w") as fh:
        #and iterates over the sequences in the jfc_out.name file (output of Jellyfish), selecting only those sequences that have an ID greater than 10 and a length of 100 nucleotides.
        for record in SeqIO.parse(jfc_out.name, "fasta"):
            if int(record.id) > 10 and len(record.seq) == 100:
                repetitive = 0
                #For each selected sequence, the code writes 10 copies of the sequence to the kmer file with a new ID starting with "kmer_1", "kmer_2", and so on.
                while 10 >= repetitive:
                    repetitive += 1
                    count += 1
                    record.id = "kmer_" + str(count)
                    SeqIO.write(record, fh, "fasta")
    print(kmer) #--> /tmp/c8319f6386ba418fb88e7e79668053b9.fasta
    #After generating the k-mer file, it creates a temporary directory in "/tmp" using tempfile.mkdtemp()
    tmp_dir = tempfile.mkdtemp(dir="/tmp/") #change /tmp/ with /Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/RACON
    print(tmp_dir) #valeria per controllo --> output: /tmp/tmpserqoo89 è una cartella
    print("RUNNING SPADES")
    #and runs the SPADES genome assembler using the k-mer file and the temporary directory as input. 
    #SPADES = "spades -s %s -o %s --only-assembler"
    #he main function of the SPAdes package is to take short reads generated by next-generation sequencing technologies and assemble them into longer contiguous sequences, called contigs.
    #the output it's stored in the tmp_dir
    spa = SPADES % (kmer, tmp_dir)
    print(spa) #--> output spades -s /tmp/c8319f6386ba418fb88e7e79668053b9.fasta -o /tmp/tmpserqoo89 --only-assembler
    try:
        # the sb.Popen() method is used to run SPADES. The stderr and stdout are redirected to sb.PIPE so that they can be captured and printed later.
        spade = sb.Popen(spa, shell=True, cwd=cwd, stderr=sb.PIPE, stdout=sb.PIPE)
        spade.communicate()
        #The assembled contigs file is created by joining the temporary directory and the file name "contigs.fasta" using the os.path.join() method.
        assembled_contigs = os.path.join(tmp_dir, "contigs.fasta")
        #A temporary BAM file is created using the tempfile.NamedTemporaryFile() method.
        bam = tempfile.NamedTemporaryFile(suffix=".bam", delete=False)
        print("RUNNING MINIMAP")
        #MINIMAP is then run using the MINIMAP_S command with the threads, assembled_contigs, and read_mapping files as input. --> ALIGNMENT NANOPORE READS AGAINST A REFERENCE GENOME (the contigs generated with spade)
        #MINIMAP_S = "minimap2 -a -x map-ont --secondary=no -t %s %s %s " --> 
        ##--> here it is mapping the reads to the assembled contigs (generated with spade) using MINIMAP
        m = MINIMAP_S % (threads, assembled_contigs, read_mapping)
        ##SAMSORT = "samtools sort -o %s -O BAM" --> the output is stored in the bam file
        ##Samtools sort is a software tool used for sorting sequencing reads in a SAM or BAM file by their alignment position on a reference genome
        ss = SAMSORT % bam.name
        print(m)
        #The output of MINIMAP is piped into SAMTOOLS_SORT and SAMTOOLS_INDEX using the sb.Popen() method
        minimap = sb.Popen(m, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        #samtools is a tool used for sorting (ordinare) sequencing reads in a SAM or BAM file by their alignment position on a reference genome
        ##the input file is the output of the minimap
        samtools_sort = sb.Popen(ss, shell=True, cwd=cwd, stdin=minimap.stdout, stdout=sb.PIPE, stderr=sb.PIPE)
        samtools_sort.communicate()
        # SAMINDEX: create the index used by samtools
        si = SAMINDEX % bam.name
        samtoos_index = sb.Popen(si, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        samtoos_index.communicate()
        #The "nanopolish index -d" command is used to create an index for a directory of raw nanopore sequencing data files in the fast5 format.  --> This option specifies that the input is a directory of raw nanopore sequencing data files.
        ni = NANOPOLISHI % (os.path.join(fast5,fastq), reads_porechop)
        print(ni)
        nanopolish_index = sb.Popen(ni, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        nanopolish_index.communicate()
        regions = []
        print("RUNNING NANOPOLISH")
        #Here, nanopolish is being run to call variants (SNPs and small indels) in the assembly.
        dict_vcf_fasta = []
        #The code iterates over the contigs assembled (with spade) starting from the reads and constructs a "region" string that specifies the contig ID and its length. 
        for record in SeqIO.parse(assembled_contigs, "fasta"):
            region = record.id + ":0-" + str(len(record.seq))
            vcf = tempfile.NamedTemporaryFile(prefix="polished.", suffix=".vcf", delete=False)
            ####Nanopolish is specifically designed to analyze the raw signal data produced by nanopore sequencers to perform: Basecalling; Variant calling: Nanopolish can detect single nucleotide variants (SNVs), insertion-deletion variants (indels), and structural variants (SVs) in DNA or RNA sequencing data by comparing the raw signal data to a reference genome or draft assembly...
            ##"--consensus": This option tells Nanopolish to generate a consensus sequence from the input reads and the detected variants. see above for all the variable specified!
            ### nanopolish is then used to call variants in the specified region using the aligned reads from the BAM file.
            #NANOPOLISHV = "nanopolish variants --consensus -o %s -w %s -r %s -b %s -g %s -t %s --min-candidate-frequency 0.1 -p 1 " \"--fix-homopolymers" --> -w" This option sets the size of the window used to call variants. "-r": This option specifies the format of the file where there are the input reads. "-b <alignments.bam>": This option specifies the input alignments in BAM format. "-g": This option specifies the format of the file where the reference genome is stored
            nv = NANOPOLISHV % (vcf.name, region, read_mapping, bam.name, assembled_contigs, threads)
            print(nv)
            nanopolish_var = sb.Popen(nv, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
            nanopolish_var.communicate()
            #The resulting VCF file is written to a temporary file "vcf.name", and the filename is added to a list if it exists and is non-empty.
            if os.path.isfile(vcf.name) and os.path.getsize(vcf.name) > 0:
                regions.append(vcf.name)
            else:
                #If the VCF file is empty, the contig is added to a list called dict_vcf_fasta
                dict_vcf_fasta.append(record)
        #Finally, nanopolish is called again to perform a consensus polish on the assembled contigs, using the VCF files for each contig as input.
        ##The "nanopolish vcf2fasta -g" command is used to generate a consensus sequence from a VCF file containing variants called by Nanopolish.
        #" ".join(regions): a string containing the paths to the VCF files for each contig joined by a space.
        nva = NANOPOLISHVA % (assembled_contigs, " ".join(regions))
        print(nva)
        ##the code here creates a path to the output file for the MLST analysis
        mlst_done = os.path.join(od, fastq + ".MLST.done.fasta")
        #Defines the path to the MLST file
        mlst = os.path.join(assembled_contigs + ".MLST.fasta")
        print (mlst) #VALERIA PER CONTROLLO ---> FILE VUOTO
        #Opens the MLST file for writing
        with open(mlst, "w") as fh:
            #executes a command to perform consensus polishing on the assembled contigs using the VCF files generated by Nanopolish.
            #The resulting consensus sequence is written to the MLST file (fh)
            nanopolish_vcf = sb.Popen(nva, shell=True, cwd=cwd, stdout=fh, stderr=sb.PIPE)
        nanopolish_vcf.communicate()
        #Creates a temporary file with a ".bam" extension to store the BAM file generated by BWA-MEM alignment
        bam = tempfile.NamedTemporaryFile(suffix=".bam", delete=False)
        ##BWAI = "bwa index %s" --> create the index for bwa
        ##Executes a command to build an index of the MLST file for use in the BWA-MEM alignment
        bi = BWAI % mlst
        bwa_index = sb.Popen(bi, shell=True, cwd="/tmp") #change /tmp with /Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/RACON
        bwa_index.communicate()
        ##Executes a command to perform the BWA-MEM alignment of the reads to the MLST file and pipes the output to a SAM file: BELOW "stdin=bwa_mem.stdout" in samtools_sort command
        bm = BWA % (threads, mlst, output.name)
        #shutil.copyfile(mlst, mlst_done) #os.path.join(cwd,output + ".fasta")) ((LUIGI))
        ####Executes a command to sort the SAM file and output the result to the BAM file
        ss = SAMSORT % bam.name
        print(bm)
        bwa_mem = sb.Popen(bm, shell=True, cwd="/tmp", stdout=sb.PIPE, stderr=sb.PIPE) #change /tmp with /Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/RACON
        samtools_sort = sb.Popen(ss, shell=True, cwd=cwd, stdin=bwa_mem.stdout, stdout=sb.PIPE, stderr=sb.PIPE)
        samtools_sort.communicate()
        #here the sorted BAM file is indexed using SAMtools
        si = SAMINDEX % bam.name
        samtoos_index = sb.Popen(si, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        samtoos_index.communicate()
        # Creates a temporary file to store the VCF output
        vcf = tempfile.NamedTemporaryFile(suffix=".vcf", delete=False)
        ##then code uses BCFtools mpileup to generate a binary VCF file using the indexed BAM file and the MLST reference sequence (consensus).
        ##The "BCFTOOLS_MP: bcftools mpileup -Ou -f s%" command is used to generate a pileup of mapped reads against a reference genome, which can be used for variant calling using BCFtools (?)
        ####---> Pileup format is a text-based format for summarizing the base calls of aligned reads to a reference sequence. This format facilitates visual display of SNP/indel calling and alignment. 
        #### After the "bcftools mpileup -Ou -f" command, you should specify the path to one or more input BAM files containing the mapped reads.
        bcfm = BCFTOOLS_MP % (mlst, bam.name) #FREEBAYES
        #the code onverts the binary VCF file to a text file using BCFtools_call and writes it to the temporary file created earlier.
        ##The "bcftools call -mv -o --ploidy 1" command is used to perform variant calling and output the results in VCF format
        #vcf.name specifies the path and name of the output VCF file.
        bcfc = BCFTOOLS_CALL % vcf.name
        print(bcfm)
        print(bcfc)
        bcftools_mpile = sb.Popen(bcfm, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        bcftools_call = sb.Popen(bcfc, shell=True, cwd=cwd, stdin=bcftools_mpile.stdout)
        bcftools_call.communicate()
        ##BCFTOOLS_IN = "bcftools index %s" --> it indexes the resulting VCF file using BCFtools index.
        bcfi = BCFTOOLS_IN % vcf.name
        ##BGZIP  is used to compress large genomic data files such as variant call format (VCF) files, BAM files, and FASTQ files
        ##The bgzip tool is typically used in conjunction with tabix, a tool for indexing Bgzip-compressed files, to efficiently query and retrieve subsets of data from large genomic data files.
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
        ##The "bcftools consensus" command is used to generate a consensus sequence in FASTA format based on a given reference genome and a set of variants from a VCF file.
        ##"-f": This specifies the path and name of the reference genome in FASTA format.
        ##"-o": This specifies the path and name of the output consensus sequence in FASTA format.
        #BCFTOOLS_CO = "bcftools consensus -f %s -o %s %s"
        ##--> the output it is written in the mlst_done file created before
        bcfco = BCFTOOLS_CO % (mlst, mlst_done, vcf_gz)
        print(bcfco)
        #print(mlst_done)
        bcftools_call = sb.Popen(bcfco, shell=True, cwd=cwd, stdout=sb.PIPE, stderr=sb.PIPE)
        bcftools_call.communicate()
        print("THIS IS MLST_DONE" + mlst_done) #valeria per controllo
        return(mlst_done, "A")
    except:
        print
        return(output.name, "B")
        print("THIS IS output" + output)


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
    os.environ["BLASTDB"] = "/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/tmp/taxdb" ##have I this directory "/data2/blastdb"? I'm changing this with "/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/tmp"
    ### ^ I have insert taxdb because we have seen that the path of bastdb has to be the same where this file is store (this file is a database with the taxonomy information of all the NCBI sequence)
    blastdb = os.environ['BLASTDB']
    if args.search != "":
        #os.path.join(blastdb, new_string) is concatenating two path components - the variable blastdb(a path) and the variable new_string(a file name) - 
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
            #? -->starting record(?)
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
                ##for result in tqdm(pool.imap(func=download_fasta, iterable=id_download), total=len(id_download)):: This line starts a loop that downloads the gi_ID of the downloaded sequences in each batch using the download_fasta() function and the parameters in the id_download list. 
                ##The pool.imap() function maps the download_fasta() function to the id_download list and distributes the batches to the worker processes in the multiprocessing pool. 
                ##The tqdm() function provides a progress bar for the downloads, and the total len of the progress bar is the len of the id list
                for result in tqdm(pool.imap(func=download_fasta, iterable=id_download), total=len(id_download)):
                    ##result_list_tqdm.append(result): This line appends the list of gi_IDs of the downloaded sequences to the result_list_tqdm list.
                    result_list_tqdm.append(result)
            ##the variable records join all the gi_ID of the downloaded sequences in the list dividing them with a new line
            records = "\n".join(result_list_tqdm)
            ### print("this are records" + records) #valeria PER CONTROLLO 
            #fp = tempfile.NamedTemporaryFile(dir="/tmp", delete= False): This line creates a named temporary file in the /tmp directory that will be used to store the gi_IDs of the downloaded sequences.
            fp = tempfile.NamedTemporaryFile(dir="/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/tmp", delete= False) ##here there is the directory dir="/tmp", TRY to change it with the same directory of the output?
            #fp.write(records.encode()): This line writes the gi_IDs downloaded sequences to the temporary file in binary format. ???
            fp.write(records.encode())
            ##creates a string cm that runs the blastdb_aliastool command with the following arguments (it creates the .n.gil database and an alias file with the .nal extention):
            #-db nt: sets the name of the database to "nt" (corrispond to all the sequences in the database nucleotide on NCBI) we change it with a smaller database that we have dowloaded previously 
            #-dbtype nucl: specifies that the database type is nucleotide
            #-gilist <path to tempfile>: specifies the path to a file containing a list of GI numbers to include in the database
            #-out <database name>: sets the output database name 
            ##--> l'output is the database ".n.gil" comprexed and an alias file with the same name but with the extension .nal
            cm = "/usr/local/bin/blastdb_aliastool -db /Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/tmp/final_xylella_fastidiosa_database.fasta -dbtype nucl -gilist %s -out %s" % (fp.name, database_name)
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
        ##FOR EACH BARCODE IN LIST_BARCODE, the code constructs a path to the directory containing the barcode data (dir_fastq). 
        dir_fastq = os.path.join(cwd, barcode) ###cwd= args.folder_fastq (la cartella dove c'è il documento fastq da analizzare)
        #If this directory exists, the code sets the result_blast flag to False, prints a message indicating which barcode is being processed, and initializes two empty lists (fastx_all and fastx_all_assembled) for storing data later in the code. 
        if os.path.exists(dir_fastq):
            result_blast = False
            print("EXECUTING BARCODE: %s" % barcode)
            fastx_all = [] 
            fastx_all_assembled = []
            ##IF I DECIDE TO MAKE THE ASSEMBLE OR TO USE THE ASSEMBLED READS
            if args.assemble or args.use_assembled:
                if args.folder_fast5 == "":
                    print("ASSEMBLY NOT DONE, PLEASE ADD FAST5 PARAMETER - RUNNING ALIGNMENT")
                    if args.use_assembled:
                        sys.exit("FAST5 NOT PASSED WHILE ASKED FOR ASSEMBLED READS FOR ALIGNEMENT")
                else:
                    ##If "folder_fast5" is not empty, it calls a function called "racon" with some parameters (od, barcode, args.threads, args.folder_fast5, cwd) and stores the result in "barcode_corr" and "spade_t".
                    barcode_corr, spade_t = racon(od, barcode, args.threads, args.folder_fast5, cwd)
                    #It then reads the "barcode_corr" (output of racon) file using Biopython's SeqIO module and for each record in the file, 
                    for record in SeqIO.parse(barcode_corr, "fasta"):
                        #it creates a temporary file using Python's "tempfile" module with a ".fasta" extension 
                        fp = tempfile.NamedTemporaryFile(dir="/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/fp_fo_assemble", suffix=".fasta", mode="w", delete=False) ##I've change the directory from /tmp/ to /Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/fp_fo_assemble
                        fo = tempfile.NamedTemporaryFile(dir="/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/fp_fo_assemble", mode="w", prefix=barcode, suffix=".blastn", delete=False) ##I've change the directory from /tmp/ to /Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/fp_fo_assemble
                        #and writes the record to it (the temporary fp file) in "fasta" format.
                        SeqIO.write(record, fp, "fasta")
                        #Finally, it appends a list of three items to "fastx_all_assembled" containing the paths of the two temporary files and the name of the "database_name".
                        ##this is useful only if you put also "use_assembled". In that case, this list become the fastx_all list used after. Otherwise, it is not used
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
            ##IF I DON'T DECIDE TO USE THE USED ASSEMBLED READS
            else:
                ##porechop is used for adapter trimming and demultiplexing of Oxford Nanopore Technologies (ONT) sequencing data
                print("RUNNING PORECHOP")
                reads_porechop = tempfile.NamedTemporaryFile(suffix=".fastq", delete=False) ###this should be the file output name of porechop with the path to this file: "/var/folders/c2.../T/" (because then, in PORECHOP % the third element is this)
                porechop_cmd = PORECHOP % (barcode, args.threads, reads_porechop.name)
                #print(porechop_cmd)
                print ("reads_porechop.name:" + reads_porechop.name)  #VALERIA PER CONTROLLO
                print("porechop command: " + porechop_cmd) #DANI PER CONTROLLO (?) ##this is the folder where it is stored /var/folders/c2/gghtq5tn77z8g_cgbf14n7nr0000gn/T/
                porechop = sb.Popen(porechop_cmd, shell=True, cwd=cwd, stderr=sb.PIPE, stdout=sb.PIPE) ###cwd= args.folder_fastq (la cartella dove c'è il documento fastq da analizzare i miei fastq pass)
                porechop.communicate()
            ##IF I DECIDE TO CORRECT THE READS
                ###The correction process involves using a program called "minimap" to align the reads against a reference genome, and then using "racon" to generate a consensus sequence by correcting the reads based on the reference.
                ###The resulting corrected reads are then used to generate a set of fasta files.
                ###If the corrected reads are empty, the original reads are used instead. 
                ###The script also creates temporary files to store intermediate results and writes the final fasta files to a list called "fastx_all". If the barcode does not end with "raw", it is appended with "raw".
                if args.correct:
                    print("START CORRECTING PROCESS")
                    #Two temporary files are created, sam and reads, with suffixes .sam and .fasta, respectively. 
                    #These files will be used to store the SAM alignment output from minimap and the corrected reads output from racon, respectively.
                    sam = tempfile.NamedTemporaryFile(suffix=".sam", delete=False)
                    reads = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
                    print("RUNNING MINIMAP")
                    #minimap2 performe an allignment of the reads demultiplexed and trimmered to each other
                    m = MINIMAP % (args.threads, reads_porechop.name, reads_porechop.name) #---> reads_porechop contiene qualcosa, ALLORA PERCHè MINIMAP è VUOTO?
                    print(m)
                    ##The minimap output is redirected to the sam temporary file.
                    minimap = sb.Popen(m, shell=True, cwd=cwd, stdout=sam, stderr=sb.PIPE)
                    minimap.communicate()
                    print("RUNNING RACON")
                    ##RACON = "racon -t %s -f %s %s %s"
                    ####"-f" specifies the input file path for the FASTQ or FASTA format reads to be used for polishing the assembly. This option tells Racon which long read sequencing data to use to correct errors in the initial draft assembly.
                    ## the sam.name and reads_porechop.name variables, which specify the name of the temporary file containing the SAM output from minimap and the input file containing the original reads, respectively.
                    r = RACON % (args.threads, reads_porechop.name, sam.name, reads_porechop.name)
                    print(sam.name) #valeria PER CONTROLLO output: /var/folders/c2/gghtq5tn77z8g_cgbf14n7nr0000gn/T/tmpjyda4hyx.sam 
                    print(r)
                    #The racon output is redirected to the reads temporary file.
                    racon_cmd = sb.Popen(r, shell=True, cwd=cwd, stdout=reads, stderr=sb.PIPE)
                    racon_cmd.communicate()
                    print(racon_cmd) #VALERIA PER CONTROLLO
                    print(reads.name) #VALERIA PER CONTROLLO
                    #This line checks if the reads file is not empty. The code performes two times the correction of the reads
                    if int(os.path.getsize(reads.name)) > 0:
                        #create a temporary file with the suffix ".sam" and assigns it to the variable sam
                        sam = tempfile.NamedTemporaryFile(suffix=".sam")
                        print("RUNNING MINIMAP")
                        ##minimap performe an allignment of the corrected reads output from racon to each other
                        m = MINIMAP % (args.threads, reads.name, reads.name)
                        ## The output is redirected to the sam file
                        minimap = sb.Popen(m, shell=True, cwd=cwd, stdout=sam, stderr=sb.PIPE)
                        minimap.communicate()
                        print("RUNNING RACON")
                        #the sam.name and reads_porechop.name variables, which specify the name of the temporary file containing the new SAM output from minimap and the input file containing the reads corrected in the first step like the original reads, respectively.
                        r = RACON % (args.threads, reads.name, sam.name, reads.name)
                        output = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
                        #the output file is called "output"
                        racon_cmd = sb.Popen(r, shell=True, cwd=cwd, stdout=output, stderr=sb.PIPE)
                        racon_cmd.communicate()
                        #This initializes a counter variable count to 0.
                        count = 0
                        ##This loops over each record in the output file, which contains the corrected reads in FASTA format.
                        for record in SeqIO.parse(output.name, "fasta"):
                            #count += 1: This increments the count variable by 1.
                            count += 1
                            #fp = This creates a temporary file with the suffix ".fasta" in the "/tmp" directory and assigns it to the variable fp.
                            fp = tempfile.NamedTemporaryFile(dir="/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/fp_fo", suffix=".fasta", mode="w", delete=False) ###here i changed the directory from "/tmp" to "/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/fp_fo"
                            #fo = This creates a temporary file with the prefix barcode and the suffix ".blastn" in the "/tmp" directory and assigns it to the variable fo.
                            fo = tempfile.NamedTemporaryFile(dir="/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/fp_fo", mode="w", prefix=barcode, suffix=".blastn", delete=False)  ###here i changed the directory from "/tmp" to "/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/fp_fo"
                            #This writes the current record to the fp file in FASTA format.
                            SeqIO.write(record, fp, "fasta")
                            #This appends a list containing the names of the fp file, the database_name variable, and the fo file to the fastx_all list. (the input file, the database the we search in NCBI and the output files)
                            fastx_all.append([fp.name, database_name,fo.name])
                    #if the reads file is empty
                    else:
                        print("THE CORRECT FILE IS EMPTY") ##valeria PER CONTROLLO
                        count = 0
                        print(reads_porechop.name)
                        #or each record, it creates two temporary files using tempfile.NamedTemporaryFile: fp with a suffix of .fasta and fo with a prefix of barcode and a suffix of .blastn. 
                        for record in SeqIO.parse(reads_porechop.name, "fastq"):
                            count += 1
                            fp = tempfile.NamedTemporaryFile(dir="/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/fp_fo", suffix=".fasta", mode="w", delete=False)  ###here i changed the directory from "/tmp" to "/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/fp_fo"
                            fo = tempfile.NamedTemporaryFile(dir="/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/fp_fo", mode="w", prefix=barcode, suffix=".blastn", delete=False)  ###here i changed the directory from "/tmp" to "/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/fp_fo"
                            SeqIO.write(record, fp, "fasta")
                            fastx_all.append([fp.name, database_name, fo.name])
                            ##If barcode does not end with "raw", it appends "raw" to it
                            if not barcode.endswith("raw"):
                                barcode = barcode + "raw"
                ##IF I DON'T DECIDE TO CORRECT THE READS
                else:
                    print(reads_porechop.name)
                    ##The following: SeqIO module to read a Fastq file, perform some processing on each record, and add the processed record to a list.
                    # The for loop iterates over the records in the Fastq file reads_porechop.name (the file with our sequences with the adapter removed and demultiplexed, that is stored in the folders .../T/) using the SeqIO.parse() function with the file format specified as "fastq".
                    #reads_porechop= (the file with our sequences with the adapter removed and demultiplexed) in a temporary file format with the extension (.fastq)
                    for record in SeqIO.parse(reads_porechop.name, "fastq"):
                        #print("records: " + record) ##DANI PER CONTROLLO
                        #For each record, the code creates a temporary file (fp) to write the sequence demultiplexed and without adapter in fasta format. --> There is a fasta file for each read of a barcode
                        fp = tempfile.NamedTemporaryFile(dir="/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/fp_fo", suffix=".fasta", mode="w", delete=False) ##here i changed the directory from "/tmp" to "/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/fp_fo"
                        #It also creates another temporary file (fo) with a name derived from the barcode (for writing the BLAST output).--> there is a blastn file for each read of a barcode
                        fo = tempfile.NamedTemporaryFile(dir="/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/fp_fo", mode="w", prefix=barcode, suffix=".blastn", delete=False) ##here i changed the directory from "/tmp" to "/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/fp_fo"
                        #The code then uses the SeqIO.write() function to write the record to the fp file in fasta format. 
                        SeqIO.write(record, fp, "fasta")
                        #It then appends a list [fp.name, database_name, fo.name] to the fastx_all list, 
                        #which contains the name of the fasta file, the name of the BLAST database to search against (database_name), and the name of the BLAST output file (.blastn) (one file for each read)
                        fastx_all.append([fp.name, database_name,fo.name]) #esempio del contenuto di questa lista: "[['/fp_fo/tmposk0hxgn.fasta', '/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/tmp/taxdb/Xylella_NOT_uncultured_All_Fields_AND_100_5000_SLEN_AND_subsp_All_Fields_', '/fp_fo/barcode18p4h_1dy8.blastn']]"
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
            ##Then, the code enters a with block where a multiprocessing.Pool object is created with a number of processes equal to args.threads.
            with Pool(processes=args.threads) as pool:
                #The Pool.imap() method is used to apply the blast() function (the luigi fuction that ask for tree element to be specified) to each element in the fastx_all listS iterable (i.e., the list of input files for BLAST, composed by: the file fp one for each read demulptiplexed and without adapter in a fasta format; the name of the database.n.gil without the extention, and the output name). 
                # ---> the function blast is the one created by Luigi: blast(elm): blastn_cline = NcbiblastnCommandline(db=elm[1], query=elm[0], evalue=0.001, out=elm[2] ,outfmt = "6 qseqid sseqid bitscore sscinames pident evalue staxids qlen" )
                #tqdm library is used to display a progress bar; the total lenght of the bar is the lenght of the fastx_all file.
                for result in tqdm(pool.imap(func=blast, iterable=fastx_all), total=len(fastx_all)): ##in questo caso stiamo specificando che elm (l'argomento della funzione blust creata da luigi) è fastx_all!!! che è una lista formata da tre elementi (ecco perchè elm [0], elm [1] etc)
                    #The results (files .blastn) are stored in the result_list variable --> there is a .blastn file for each read of the barcode
                    result_list.append(result)
            #print(result_list) #valeria PER CONTROLLO. => A list with the path and the name of the different .blastn files ['/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/fp_fo/barcode18nzqh9jla.blastn', '/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/fp_fo/barcode18rawvmacyc9o.blastn', '/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/fp_fo/barcode18raw5h0kzx5a.blastn', '/Users/it031376/Desktop/Sequentia_python/SCRIPT_LUIGI/fp_fo/barcode18rawwl56hp1w.blastn', 
            result_list_tqdm = []
            ##od = os.getcwd() --> The os module provides a way to interact with the operating system in Python. The getcwd() function from the os module returns a string representing the current working directory. The code od = os.getcwd() assigns the value returned by getcwd() to the variable od.
            #we are defining the name of the output file of blust for each of the analized barcode 
            name_blast = os.path.join(od, args.output + barcode + ".blast.txt")
            #write a new file that will be named with the name_blust defined before
            with open(name_blast, "w") as new_file:
                ##for each result name in the result_list (the .blastn file names) open a list named data_single_file
                for name in result_list:
                    data_single_file = []
                    ##For each result name (.blastn files) in the list, it opens the file using the with statement and assigns it to the variable f
                    with open(name) as f:
                        #It then loops through each line (the informations of each possible match that this specific read) in the file using a for loop, writes the line to the new_file, and appends it to the data_single_file list.
                        for line in f:
                            new_file.write(line) ##a questo punto ha scritto tutti i file di output "nome_barcode.blast.txt" 
                            data_single_file.append(line) ##li salva anche in una lista perchè serve per i passaggi successivi. una lista per ogni barcode
                        #After reading all lines in the file, it writes a newline character to new_file.
                        new_file.write("\n")
                    #print(data_single_file) #valeria per controllo --> several listS with the the content of "nome_barcode.blast.txt" without newline ex. output ['9bbce8e7-3a75-458f-831b-9c00f58088ed\tgi|2293497712|gb|MW940137.1|\t1072\tN/A\t94.828\t0.0\t0\t683\n', '9bbce8e7-3a75-458f-831b-9c00f58088ed\tgi|397790742|gb|JQ290498.1|\t1066\tN/A\t94.684\t0.0\t0\t683\n',
                    #and then joins the data_single_file listS (one list for each read) into a single string using the newline character as a separator of each line and assigns it to the variable blast_out_single. 
                    blast_out_signle = "\n".join(data_single_file)
                    ##print(blast_out_signle) valeria PER CONTROLLO --> they are the lines inside of the "nome_barcode.blast.txt" file separated by newline (I know because I printed it)
                    ##The blast_out_single string is then appended to the result_list_tqdm list.
                    result_list_tqdm.append(blast_out_signle)
                    #Finally, the code uses the os module's remove() function to delete the file that was just processed.
                    os.remove(name)
            #print(result_list_tqdm) #VALERIA PER CONTROLLO --> is a unique list with the lines inside of the "nome_barcode.blast.txt" with /n between the several line for one read and with a "," between the line of each new read ex. \t976\tN/A\t92.609\t0.0\t0\t683\n\n9bbce8e7-3a75-458f-831b-9c00f58088ed\tgi|284814416|gb|FJ610176.1|\t976\tN/A\t92.598\t0.0\t0\t683\n\n9bbce8e7-3a75-458f-831b-9c00f58088ed\tgi|313104579|gb|HM243599.1|\t970\tN/A\t92.464\t0.0\t0\t683\n', '']
                    #### ^ The purpose of this code is to merge the contents of multiple files into a single file and to store the individual file contents as a list of strings in result_list_tqdm. The with statement is used to ensure that each file is closed after it has been read and written to new_file, and the os.remove() function is used to delete the original file after it has been processed to save space on the disk.
           
            #If an element in the list result_list_tqdm is not empty and has not yet been flagged as a hit (result_blast is initially set to False), then result_blast is set to True. (before there is written that, if the path to the fastq file exist the result_blast should be false)
            for line in result_list_tqdm:
                if line != "" and not result_blast:
                    result_blast = True
            #If there are no hits in the result, the loop continues to the next result using the continue statement.
            #A flag in Python acts as a signal to the program to determine whether or not the program as a whole or a specific section of the program should run. In other words, you can set the flag to True and the program will run continuously until any type of event makes it False. Then the program, loop, or whatever you're using a flag for will stop.
            if not result_blast:
                continue
            #For each result with a hit, the code creates a dictionary called read_best_hit. This dictionary will eventually contain the best match for each query sequence.
            read_best_hit = {}
            ##the code below parses each line of output from a sequence alignment file, then sorting and selecting the highest-scoring alignment(s) for each read based on certain criteria. The selected alignment(s) are stored in a dictionary called read_best_hit.
            ##The code iterates through the lines of the result again, splitting each line by the tab character to separate the different fields of the output.
            #The loop over result_list_tqdm processes each line of the file by splitting it into a list of elements using tab ('\t') as a separator.
            for output in result_list_tqdm:
                if output != "":
                    align_species = {}
                    #each line of block of match informations of a sigle read is devided with /n in result_list_tqdm. With this function (split) we divided each of these lines in different element and this new format is stored in the variable align.
                    align = output.split("\n") 
                    #print(align) #esempio di un pezzo di align '9bbce8e7-3a75-458f-831b-9c00f58088ed\tgi|827353834|gb|KP282072.1|\t976\tN/A\t92.598\t0.0\t0\t683', '', '9bbce8e7-3a75-458f-831b-9c00f58088ed\tgi|827353836|gb|KP282073.1|\t976\tN/A\t92.598\t0.0\t0\t683', '', '9bbce8e7-3a75-458f-831b-9c00f58088ed\tgi|827353838|gb|KP282074.1|\t976\tN/A\t92.598\t0.0\t0\t683', '', '9bbce8e7-3a75-458f-831b-9c00f58088ed\tgi|284814410|gb|FJ610173.1|\t976\tN/A\t92.609\t0.0\t0\t683', '', '9bbce8e7-3a75-458f-831b-9c00f58088ed\tgi|284814416|gb|FJ610176.1|\t976\tN/A\t92.598\t0.0\t0\t683', '', '9bbce8e7-3a75-458f-831b-9c00f58088ed\tgi|313104579|gb|HM243599.1|\t970\tN/A\t92.464\t0.0\t0\t683', '']
                    #then, each element of a line in the block of match informations, that are divided with a /t, is divided in different element with this function and the result is store in the variable align_elm
                    for single_align in align:
                        align_elm = single_align.split("\t")
                        #If the length of the single_align (that are the elements of the variable "align") is greater than 3, the program proceeds to the next step. Otherwise, it ignores the line.
                        #--> The if statement checks that the line has at least four elements (query ID, subject ID, percent identity, and alignment length? (qseqid sseqid bitscore sscinames)), and that the percent identity meets a certain threshold (args.percentage) and that the score (bitscore?) is greater than 200.
                        if len(single_align) > 3:
                            #If these conditions are met, the code creates a dictionary called align_species_score with the perc. of identity (pident?) as the key and a list containing the scientific name (sscinames?) as the value. 
                            align_species_score = {}
                            read = align_elm[0]
                            score = float(align_elm[2])
                            align_species_score[align_elm[4]] =  [align_elm[3]]
                            #print(align_species_score) ##VALERIA PER CONTROLLO ex. result--> {'95.819': ['N/A']} each of the sscinames is N/A
                            ###if bitscore (the element 2 in the align_elm variable) is > 200 and pident (element 4 in the align:elm variable) meets the threshold of args.percentage
                            if score > 200 and float(align_elm[4]) >= args.percentage :
                                #If the score is already present as a key in the align_species dictionary, the code iterates over the keys of the value (which is a dictionary of species and their respective percent identity)
                                if score in align_species:
                                    #The align_species_score dictionary is then added to the align_species dictionary with the score of the alignment as the key and the align_species_score dictionary as a value.
                                    for key in  align_species[score]:
                                        #then checks whether the percent identity of the current alignment (align_elm[4]) is already present in the align_species_score dictionary. 
                                        if align_elm[4] in align_species_score:
                                            #If it is, the scientific name (align_elm[3]) is appended to the list of species names for that percent identity. 
                                            align_species[score][key].append(align_elm[3])
                                        else:
                                            #If not, a new key-value pair is added to the align_species dictionary, where the bitscore is the key and the list of species names is the value.
                                            align_species[score][key] = [align_elm[3]]
                                else:
                                    #If the score is not already present in the align_species dictionary, a new key-value pair is added, where the score is the key and the align_species_score dictionary (with the percent identity as key and the list of species names as value) is the value.
                                    align_species[score] = align_species_score
                    #print(align_species) #VALERIA PER CONTROLLO there are several dictionary ex. of a part of the results: {765.0: {'96.186': ['N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A']}, 760.0: {'95.975': ['N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A']}, 726.0: {'94.715': ['N/A']}, 721.0: {'94.503': ['N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A']}, 715.0: {'94.292': ['N/A', 'N/A', 'N/A']}, 710.0: {'94.080': ['N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A']}, 704.0: {'93.869': ['N/A']}, 699.0: {'93.658': ['N/A']}, 339.0: {'94.196': ['N/A']}, 274.0: {'95.882': ['N/A']}}
                    #The code then checks if any alignments were found for the query sequence 
                    if align_species:
                        #(align_species is not empty). If so, it sorts the align_species dictionary by score in descending order and selects the dictionary with the highest score 
                        list_align_species = list(align_species.items())
                        list_align_species.sort(reverse=True)
                        #(align_species_ident). It then sorts align_species_ident by percent identity in descending order and selects the list of percent identities with the highest value 
                        align_species_ident = list_align_species[0][1]
                        #(align_species_ident_b). This list of percent identities is then added to the read_best_hit dictionary with the query sequence ID (align_elm[0]) as the key.
                        list_align_species_b = list(align_species_ident.items())
                        list_align_species_b.sort(reverse=True)
                        align_species_ident_b = list_align_species_b[0][1]
                        read_best_hit[read] = align_species_ident_b
            #print(read_best_hit) #VALERIA PER CONTROLLO ex. of the end of read_best_hit dictionary --> '15d0e68e-b6b4-4143-9852-eadc946f211e': ['N/A', 'N/A', 'N/A'], '2ff01489-113e-41d8-a61b-20a45638e687': ['N/A'], '28b9fdfa-c2ba-4648-8507-762ceb1399e5': ['N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A'], '9bbce8e7-3a75-458f-831b-9c00f58088ed': ['N/A']}
            ####This code is taking a list of matches (read_best_hit) and creating a new list (dict_match) that contains the key-value pairs from read_best_hit. 
            dict_match = []
            ##match are the different key-value elements in read_best_hit dictionary.
            for match in read_best_hit:
                dict_match.append([match, read_best_hit[match]])
            result_list_tqdm = []
            ####Then it creates a pool of worker processes using the Python multiprocessing module 
            ###and uses the "imap" method to parallelize the "compute_count" function over the items in dict_match (=the dict_match is the object of luigi function "counpute_count" that returns a list containing the species_dict and genera_dict dictionaries of all of the Blast hits for the given sequence in dict_match). 
            with Pool(processes=args.threads) as pool:
                for result in tqdm(pool.imap(func=counpute_count, iterable=dict_match), total=len(dict_match)):
                    ###The results of each parallel computation are stored in a list called "result_list_tqdm".
                    result_list_tqdm.append(result)
            #print(result_list_tqdm) ##VALERIA PER CONTROLLO ex. first part of the result output: [[{}, {}], [{'N/A': ['fad34af2-cce1-41a3-8012-900420a63dc2']}, {'N/A': ['fad34af2-cce1-41a3-8012-900420a63dc2']}], [{'N/A': ['ebb84d75-16bd-4e79-8b39-052c9ff14dab']},
            #after this, the code creates a temporary file using the "tempfile" module with a unique name based on the barcode, and writes the output to this file. 
            kt_barcode = tempfile.NamedTemporaryFile(suffix=".txt", prefix=barcode, delete=False, mode = "w")
            files_kt_list.append(kt_barcode.name)

            #The output is generated by looping over the items in "result_list_tqdm", 
            for value in result_list_tqdm:
                #and for each item, looping over the keys in the first element of the item (which is a dictionary).
                for key in value[0]:
                    ##For each key, it looks up the corresponding taxonomic ID using the "get_name_translator" method from the BioPython library, 
                    ##and writes the read ID and taxonomic ID to the temporary file.
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
                            ## If the taxonomic ID is not found, the code prints a warning message and continues to the next iteration of the loop.
                            print("NO TAXID FOR " + taxid[key][0])
                            continue
            ##This code seems to be calculating the frequency of species and genera in the result_list_tqdm. 
            species_dict = {}
            genera_dict = {}
            ##The code loops through each dictionary in result_list_tqdm, extracts the species and genera keys, and updates the counts in species_dict and genera_dict.
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
            #print(species_dict) #valeria PER CONTROLLO ---> OUTPUT: {'N/A': 182}
            #print(genera_dict)#valeria PER CONTROLLO ----> OUTPUT: {'N/A': [182]}
            total_reads_mapped = sum(species_dict.values())
            #After counting the species and genera, it calculates the minimum number of reads needed for a species to be considered present, based on the args.min argument,
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
            ##and then calculates the final percentage of each species based on the retained species. 
            for key in species_retain:
                specied_final[key] = species_retain[key] / sum(species_retain.values()) * 100
                number_species.append(key)
            ## It then stores this information in dict_species_all and dict_table_all.
            dict_species_all[barcode] = specied_final
            dict_table_all[barcode] = dict(sorted(species_dict.items(), key=lambda item : item[1]))
            #print(dict_species_all) #valeria PER CONTROLLO OUTPUT: --> {'barcode18raw': {'N/A': 100.0}}
            #print(dict_table_all) #valeria PER CONTROLLO output: --> {'barcode18raw': {'N/A': 182}}
    
    # KTIMPORTTAX
        ###if the path of the fastq file (dir_fastq) doesn't exist:
        else:
            print("PATH FOR BARCODE " + barcode + " DO NOT EXISTS. CHECK BARCODE NAME OR PATH" )
    #print(files_kt_list) #valeria PER CONTROLLO --> OUTPUT: is the path to the folders with the file with our sequences with the adapter removed and demultiplexed (porechop output): ['/var/folders/c2/gghtq5tn77z8g_cgbf14n7nr0000gn/T/barcode18rawu8sx7o99.txt'] 
    
    ##This section of code generates a stacked bar plot of the relative abundance of species in each barcode, based on the output of the previous steps. 
    ##The plot is saved as a file with the name specified by figure_file.
    if len(files_kt_list) > 0:
        #This section of the code creates a single string, files_kt, by joining the paths in files_kt_list with a space character in between
        files_kt = " ".join(files_kt_list)
        #This string is then used to create a command k that is passed to subprocess.Popen().
        ##KTIMPORTTAX is the function "ktImportTaxonomy -tax %s -o %s %s", where the argument -tax is , and -o specify the output file for the annotated phylogenetic tree. This option accepts the path and filename of the output file (but er put only the name). e il secondo elemento %s è il file di input  
        ####in this case, blastdb is a path to the output folder (?)
        k = KTIMPORTTAX % ("/Users/it031376/local/bin/KronaTools-2.8.1/taxonomy", args.output + ".html", files_kt)
        print(k) # OUTPUT: --> ktImportTaxonomy -tax /Users/it031376/local/bin/KronaTools-2.8.1/taxonomy -o output.html /var/folders/c2/gghtq5tn77z8g_cgbf14n7nr0000gn/T/barcode18raw1fgyrkjc.txt
        ##ktimport is an instance of the Popen class, which starts a new process to execute the k command. 
        ####The shell=True argument allows us to use shell syntax in the k command. cwd=od specifies the working directory where the command should be run (THIS COULD BE THE PROBLEM!). 
        #####-->od = os.getcwd() --> The os module provides a way to interact with the operating system in Python. The getcwd() function from the os module returns a string representing the current working directory.
        ktimport = sb.Popen(k, shell=True, cwd=od)
        #he communicate() method is called on ktimport to wait for the command to finish executing before continuing with the rest of the code.
        ktimport.communicate()
        #This part creates a stacked bar chart using the Pandas DataFrame species_pd that is generated from the dictionary dict_species_all containing the relative abundance of each species for each barcode.
        species_pd = pd.DataFrame.from_dict(dict_species_all, orient='index')
        #First, the code creates a list of species names called species_name_colors, which is used to generate a color map for the plot.
        species_name_colors = [name for barcode in dict_species_all for name in dict_species_all[barcode]]
        species_name_colors = list(dict.fromkeys(species_name_colors))
        #Then, the code creates a linear segmented colormap called cmap1 using cmap2, a list of the colors corresponding to each species in species_name_colors
        cmap = convert_color(species_name_colors)
        cmap2 = []
        for i in cmap:
            cmap2.append(i[1])
        species_abs = pd.DataFrame.from_dict(dict_table_all, orient='index')
        species_abs.transpose().to_excel(name_table)
        cmap1 = LinearSegmentedColormap.from_list("my_colormap", cmap2)
        #The code then creates a new figure and plots the data in species_pd as a stacked bar chart using the color map cmap1. 
        #The legend is placed outside of the plot to the right, and the title, x-axis label, and y-axis label are set. Finally, the plot is saved as a file with the name specified by figure_file.
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
    ##If files_kt_list is empty, then the code prints "NOTHING TO PLOT" and does not generate a plot.
    else:
        print("NOTHING TO PLOT")




if __name__ == '__main__':
    analysis()