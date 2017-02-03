usage: 1.primer3InputMaker.py [-h] [-t TABLE] [-r REFERENCE]
                              [-s AMPLICON_SIZE] [-i INFILE] [-o OUTFILE]

optional arguments:
-h, --help            show this help message and exit
-t TABLE, --table TABLE
    chr start end sample refSample ref alt [whatever...]
-r REFERENCE, --reference REFERENCE
-s AMPLICON_SIZE, --amplicon-size AMPLICON_SIZE
-i INFILE, --infile INFILE
-o OUTFILE, --outfile OUTFILE

#Example commands
$ python 1.primer3InputMaker.py -t markerList.xls -r aSpecies_genome.fna -s 100-200
$ mkdir forrev
$ mv *.for forrev/
$ mv *.rev forrev/
$ python 2.primer3OutputParser.py primer3.output.100-200 primer3.input.info

