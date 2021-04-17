# parse_gtf.py
This program has been designed to obtain information about transcripts from a
GTF file. The gtf files that have been used for this can be found on
the [Ensembl downloads website](https://uswest.ensembl.org/info/data/ftp/index.html).
This program took approximately 1 hour and 15 minutes to process the human gtf
file version 98 from Ensembl using a 2018 MacBook Pro and using the -t and -c options.

A comma separated value (.csv) file containing information about each transcript will be
generated. In addition, if you would like a separate file with the information
only about the longest transcript per gene or the transcript with the longest
exonic length, you can specify to have those files created with the -t or -c
options.

Example usage: python parse_gtf.py --gtf_file
Homo_sapiens.GRCh38.96.gtf --output_prefix HG38.96 -c

optional arguments:

  -h, --help

                        show this help message and exit

  -g GTF_FILE, --gtf_file GTF_FILE

                        This is the GTF file from which you would like to
                        parse information. If you do not provide a file, the
                        program will attempt to find a GTF file in your
                        present working directory.

  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX

                        This is the name of the output file that you would
                        like created. If you do not specify a prefix, the
                        prefix will be taken from the GTF file that you
                        provided.

  -t, --max_transcript_length

                        Specify this option if you would also like to make a
                        file with only the transcripts containing the longest
                        transcript length.

  -c, --max_cds

                        Specify this option if you would also like to make a
                        file with only the transcripts containing the longest
                        coding sequence (CDS).
