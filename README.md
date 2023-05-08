# Project Title

pipeline for Nanopore data analysis

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Installation

A step-by-step guide on how to install and run the project on a local machine.

## Usage

"Instructions and examples on how to use the project, including code snippets and/or screenshots."  
To run this script it's necessary to download a database from NCBI with the sequences with which to align our queries.  
To do this you have to: 
- click on this site: https://www.ncbi.nlm.nih.gov
- select the 'nucleotide' database and in the search bar at the top search for the database of interest
- then click on "send to" in the top right-hand corner
- select "Complete Record" -> "File", the format "FASTA" and tick the box with Show GI
![example: HOW TO DOWNLOAD NCBI DATABASES](https://github.com/dani-julian/Sapienza_environmetal_biology_Valeria/blob/869202326a3aefd852e83987837cbd2708a0bbad/example_NCBI.png)  
In order to have the database with the correct format to run the script, you must add the GI and scientific name informations of the organisms associated with the downloaded sequences.  
To do this:
- download the file with the taxonomic information of the sequences on NCBI taxdb.tar.gz from this link: https://ftp.ncbi.nlm.nih.gov/blast/db/
- save the file in the same directory of the downloaded NCBI database
- download from NCBI the information of the sequence ID and taxa ID and merge them into a txt file with two columns separated by a space:  
$ esearch -db nucleotide -query "name_of_your_database" | efetch -format acc > name_of_first_txt_file.txt  
$ esearch -db nucleotide -query "name_of_your_database" | efetch -format docsum | xtract -pattern DocumentSummary -element TaxId > name_of_second_txt_file.txt  
$ paste -d " " name_of_first_txt_file.txt name_of_second_txt_file.txt > definitive_txt_table.txt
- Build a BLAST database with your (local) sequences adding sequence ID and taxa ID information:
$ makeblastdb -in "name_of_your_database.fasta" -dbtype nucl -input_type fasta -parse_seqids -taxid_map definitive_txt_table.txt

## Contributing

Information on how to contribute to the project, including guidelines for code style and pull request formatting.

## License

This project is licensed under the [License Name] License - see the [LICENSE.md](LICENSE.md) file for details.

