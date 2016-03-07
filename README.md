#PredZinc

##Description
PredZinc is a program for predicting zinc-binding sites in proteins from their
amino acid sequences. The program is witten in c/c++ and bash shell scripts.
Currently, PredZinc can be run on Linux and Windows (with Cygwin). 

PredZinc is copyrighted (c) to Nanjiang Shu, Structural Chemistry, Stockholm
University, Sweden and is free for academic use.

###Web-server
The web server of PredZinc is available at 
http://predzinc.bioshu.se


##Reference
Nanjiang Shu, Tuping Zhou and Sven Hovmoller (2008). "Prediction of
zinc-binding sites in proteins from sequence." Bioinformatics 24(6):
775-782. http://www.ncbi.nlm.nih.gov/pubmed/18245129

##Contact
Nanjiang Shu    Science for Life Laboratory, Stockholm

####    Email
nanjiang.shu@scilifelab.se

nanjiang.shu@mmk.su.se

##Installation:

Download the package by

    $ git clone https://github.com/nanjiangshu/predzinc

Enter the predzinc directory

    $ cd predzinc

Fetch the data set by 
    
    $ git lfs fetch
    $ git lfs checkout

If you don't have git-lfs installed, please installed follow the instructions
at https://git-lfs.github.com/

After that, run

    $ make 
    $ make install

Make sure that the NCBI nr database formatted for PSI-BLAST is installed. The
environmental variable BLASTDB points to the directory storing nr blast
database needed by PSI-BLAST

    $ export BLASTDB=path-storing-blast-nr-database


## Usage
```
usage: predzinc.sh  [-cpu INT] [-outpath DIR] [-blastdb STR] [-db STR]
                    [-pssm FILE] [-pssmfilelist LISTFILE]
                    [-showexample] [-verbose] [-version] [-not-clean] [-h]
                    [-outname STR]
                    FILE

Note: The supplied sequence FILE should be in Fasta format
      Do not use ';' in the rootname of the input file 
      since it will be used as field separator.        

Options:
 -cpu          INT   Set the number of cpu cores to be used to run
                     PSI-BLAST, (default: 1)
 -outpath      DIR   Output the result to the specified path, (default: ./)
 -outname      STR   Output name, (default: query)
 -blastdb      FILE  Database for psi-blast, (default: $BLASTDB/nr)
 -db           STR   Database for PredZinc, (default: passe)
 -pssm         FILE  Supply pssm file in PSI-BLAST -Q flag output format,
                     if supplied, PSI-BLAST will not run
 -pssmfilelist FILE  Supply a file with a list of pssm files for batch mode
                     prediction
 -not-clean          If supplied, the temporary files will be kept
 -verbose            Print verbose information 
 -version            Print version 
 -showexample        Print How to use with examples
 -h|--help           Print this help message and exit

Note that the rootname of the file should be <= 100 characters.

```

##Examples:
In the subfolder `test`

* Carry out the prediction by supplying sequence file (with one or more sequences) in FASTA format

        $ ../predzinc.sh test3seq.fa -outpath out1

* Carry out the prediction by supplying a single pssm file

        $ ../predzinc.sh --pssm 2AZ4A.pssm -outpath out2

* Carry out the prediction by supplying a list of pssm files

        $ ../predzinc.sh --pssmfilelist test.pssmfilelist -outpath out3

##Others
The blastpgp used in this version of PredZinc is 2.2.25.
If you want to use a different version of PSIBLAST, please copy the program
'blastpgp' to $PredZinc/bin and the corresponding matrix file to $PredZinc/data

    $ cp blastpgp $PredZinc/bin
    $ cp BLOSUM62 $PredZinc/data

The gist-svm used in this version of PredZinc is 2.1.1
If you want to use a different version of gist-svm, please copy the program 
'gist-svm-train', 'gist-classify' and 'gist-score-svm' to $PredZinc/bin

    $ cp gist-svm-train $PredZinc/bin
    $ cp gist-classify $PredZinc/bin
    $ cp gist-score-svm $PredZinc/bin

##Running time

It takes on average about 6 minutes to predict a protein sequence with 300
amino acids, when running on a single core with 2GHZ cpu and 4GB RAM. Most
time is taken by PSI-BLAST for building the sequence profile. When the
sequence profile is obtained, it takes about 30 seconds to get the
prediction result. With modern computes, the speed can be significantly faster.

##Result examples:

in the `test` folder, run

    $ ../predzinc.sh 2AZ4A.aa

or if the pssm file is already built, run

    $ ../predzinc.sh --pssm  2AZ4A.pssm

The result of the predict will be output to "query.predzinc.report"
as shown below. The explanation of the colums in the result file is as follows

    *Res       : three-letter amino acid code
    *SerialNo  : the residue number in sequence, starting from 1
    *Score     : predicted zinc-binding score for the residue, ranging from 0 to 1


```
Zinc-binding site prediction by PredZinc version 1.3 (c) Shu.

Reference: 
   Shu, N., Zhou, T. and Hovmoller, S. (2008) Prediction of zinc-binding
   sites in proteins from sequence, Bioinformatics, 24, 775-782.

Parameters:

Training set: 2727 PDB chains
blastpgp version: 2.2.25
blast nr database: 10,083,419 sequences         

//BEGIN query 1
Your input amino acid sequence (containing 429 aa) is
MESKAKTTVTFHSGILTIGGTVIEVAYKDAHIFFDFGTEFRPELDLPDDHIETLINNRLVPELKDLYDPR
LGYEYHGAEDKDYQHTAVFLSHAHLDHSRMINYLDPAVPLYTLKETKMILNSLNRKGDFLIPSPFEEKNF
TREMIGLNKNDVIKVGEISVEIVPVDHDAYGASALLIRTPDHFITYTGDLRLHGHNREETLAFCEKAKHT
ELLMMEGVSISFPEREPDPAQIAVVSEEDLVQHLVRLELENPNRQITFNGYPANVERFAKIIEKSPRTVV
LEANMAALLLEVFGIEVRYYYAESGKIPELNPALEIPYDTLLKDKTDYLWQVVNQFDNLQEGSLYIHSDA
QPLGDFDPQYRVFLDLLAKKDITFVRLACSGHAIPEDLDKIIALIEPQVLVPIHTLKPEKLENPYGERIL
PERGEQIVL

The following 5 residues were predicted as zinc-binding for protein "2az4_A" (with score >=  0.450), sorted by scores

Res SerialNo  Score

HIS       92  0.878
HIS       94  0.878
HIS      167  0.844
ASP       96  0.669
HIS      404  0.614

Prediction scores for the rest 84 selected residues, sorted by scores

Res SerialNo  Score

HIS      382  0.372
GLU      237  0.321
ASP      189  0.245
HIS       97  0.114
HIS      193  0.091
GLU      409  0.029
HIS      182  0.025
ASP       65  0.022
GLU      205  0.019
GLU      216  0.019
GLU      417  0.016
GLU      303  0.016
ASP      166  0.013
HIS       31  0.010
GLU      161  0.007
GLU      296  0.006
ASP      324  0.005
GLU      136  0.005
GLU      199  0.005
GLU      396  0.005
HIS       50  0.004
GLU       39  0.004
ASP       35  0.004
GLU      115  0.004
GLU      422  0.003
GLU      157  0.003
CYS      379  0.003
ASP      389  0.003
HIS       12  0.003
GLU      341  0.002
GLU      282  0.002
GLU      198  0.002
ASP       29  0.002
GLU       62  0.002
ASP       68  0.002
HIS      195  0.001
HIS      347  0.001
GLU      238  0.001
GLU      309  0.001
GLU       52  0.001
HIS       85  0.001
GLU      224  0.001
ASP       45  0.001
GLU      425  0.001
GLU      273  0.001
GLU       74  0.001
GLU       24  0.001
ASP       49  0.001
GLU      137  0.001
GLU      226  0.001
ASP      371  0.000
ASP      105  0.000
ASP      168  0.000
GLU      291  0.000
ASP      349  0.000
GLU      211  0.000
ASP      128  0.000
GLU      315  0.000
GLU      250  0.000
GLU      248  0.000
GLU      143  0.000
GLU      386  0.000
HIS       76  0.000
HIS      209  0.000
CYS      204  0.000
ASP      357  0.000
ASP      355  0.000
ASP       48  0.000
GLU      266  0.000
HIS      243  0.000
ASP      387  0.000
ASP      337  0.000
ASP      319  0.000
ASP      181  0.000
GLU       79  0.000
GLU       43  0.000
ASP       80  0.000
GLU      412  0.000
ASP       82  0.000
ASP      151  0.000
ASP      239  0.000
GLU        2  0.000
ASP      228  0.000
ASP      365  0.000
//End
```
