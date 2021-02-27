#!/bin/bash
# copyright(c) Nanjiang Shu, Structural Chemistry, Stockholm University, Sweden
# 2007-10-15, 
# Contact: nanjiang.shu@mmk.su.se
# 
# Description:
#   Predicting zinc-binding sites, i.e. which amino acids are binding to zinc,
#   from amino acid sequences
# 
# Version: 1.4
#
# Reference:
#   Prediction of zinc-binding sites in proteins from sequence
#   Nanjiang Shu, Tuping Zhou and Sven Hovmoller 
#   Division of Structural Chemistry, Stockholm University, SE-106 91 Stockholm, Sweden.
#
# usage: ./predzinc.sh [options] sequence-file-in-fasta-format
#
# Updated 2011-10-20

function AddAbsolutePath() #$path#{{{
{
    local var=$1
    if [ "${var:0:1}" != "/" ];then
        var=$PWD/$var # add the absolut path
    fi
    echo $var
    return 0
}
#}}}
function IsProgExist()#{{{
# usage: IsProgExist prog
# prog can be both with or without absolute path
{
    type -P $1 &>/dev/null || { echo "The program \"$1\" is required but it's not installed. Aborting $0" >&2; exit 1; }
}
#}}}
function IsPathExist()#{{{
# supply the effective path of the program 
{
    if ! test -d $1; then
        echo "Directory $1 does not exist. Aborting $0" >&2
        exit
    fi
}
#}}}

PREDZINC=`dirname $0`
PREDZINC=`AddAbsolutePath $PREDZINC`
export PREDZINC

LD_LIBRARY_PATH=$PREDZINC/lib:$LD_LIBRARY_PATH

export LD_LIBRARY_PATH
export BINPATH=$PREDZINC/bin
export DATADIR=$PREDZINC/data

export BLASTBIN=$BINPATH
export BLASTMAT=$DATADIR
blastbin=$BLASTBIN
blastdbname="uniref90.fasta"
gistbin=$PREDZINC/bin

usage="
`cat $PREDZINC/version`

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
 -blastdb      FILE  Database for PSI-BLAST, (default: uniref90.fasta)
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
"

examples="
Examples:

# In the subfolder \"test\" 
# Carry out the prediction by supplying sequence file (with one or more
# sequences) in Fasta format
    ../predzinc.sh test3seq.fa

# Carry out the prediction by supplying a single pssm file
    ../predzinc.sh --pssm 2AZ4A.pssm

# Carry out the prediction by supplying a list of pssm files
    ../predzinc.sh --pssmfilelist test.pssmfilelist
"
debugoptioninfo="
 --tmpdir  DIR     Set the directory to store temporary files, 
                   by default it will e created under /tmp by mktemp
"
MIN_LENGTH_SEQ=10
MAX_LENGTH_SEQ=10000
SEQ_MODE_AASEQ=0
SEQ_MODE_PSSM=1

VERSION="PREDZINC version 1.4, 2011-10-11"
if [ -s $PREDZINC/version ]; then 
    VERSION=`cat $PREDZINC/version`
fi
function CheckSequence() #$file $mode #{{{
{
    local file=$1
    local mode=$2
    local lengthseq=0
    local cntCHDE=0
    if [ $mode -eq $SEQ_MODE_PSSM ]; then
        lengthseq=`awk 'function isnum(x){return(x==x+0)}BEGIN{cnt=0}{if (isnum(substr($1,1,1))) {cnt+=1}}END{print cnt}'  $file ` 
        cntCHDE=` awk 'function isnum(x){return(x==x+0)}BEGIN{cnt=0}{if (isnum(substr($1,1,1))) {if( index("CHDE",$2) != 0 ){cnt+=1;} } } END{print cnt}' $file ` 
    elif [ $mode -eq $SEQ_MODE_AASEQ ] ; then 
        lengthseq=`awk 'BEGIN{cnt=0}/^[^>]/{str=$0;sub(/ \r\t\n/,"",str); cnt+=length(str);}END{print cnt}'  $file ` 
        cntCHDE=` awk 'BEGIN{cnt=0}/^[^>]/{str=$0;sub(/ \r\t\n/,"",str); num=gsub(/[CHDEchde]/,".",str); cnt+=num;}END{print cnt}'  $file ` 
    else
        echo "Unrecognized sequence mode for $file. Ignore" >&2
        echo "1"
        return
    fi

    if [ $lengthseq -lt  $MIN_LENGTH_SEQ -o $lengthseq -gt $MAX_LENGTH_SEQ ];  then
        echo "Length of the sequence for $file ($lengthseq) is out of range ($MIN_LENGTH_SEQ - $MAX_LENGTH_SEQ). Ignore."  >&2
        echo "1"
        return
    elif [ $cntCHDE -lt  0 ];  then
        echo "The sequence in $file has no cysteins, histidines, aspartates and glutamates (CHDE). This protein is unlikely to be zinc-binding. Ignore."  >&2
        echo "1"
        return
    else
        echo "0"
        return 
    fi
}
#}}}
function PrintHelp()#{{{
{
    echo "$usage"
}
#}}}
function PrintVersion()#{{{
{
    echo "$VERSION"
}
#}}}
function CheckErrMsg() #$errFile#{{{
{
    local errFile=$1
    if [ -s "$errFile" ]; then 
        cat $errFile >&2 
        echo "1"
    else 
        echo "0"
    fi  
}
#}}}

function BuildProfileFromAASeq() #seqfile outpath seqfilelist pssmfilelist modmfilelist#{{{
{
    local seqfile=$1
    local outpath=$2

    local seqfilelist=$3
    local pssmfilelist=$4
    local modmfilelist=$5

    local basename=`basename $seqfile`
    local rootname=${basename%.*}
    local pssmfile=$outpath/${rootname}.pssm
    local blastfile=$outpath/${rootname}.blast
    local chkfile=$outpath/${rootname}.chk
    echo "Running PSI-BLAST for $seqfile..."
    if [ $isPrintVerboseInfo -eq 1 ]; then
        echo "$blastbin/blastpgp -a $numCPU -d "$blastdbname"  -j 3 -h 0.001 -i $seqfile -o $blastfile -m 9 -Q $pssmfile -C $chkfile "
    fi

    $blastbin/blastpgp -a $numCPU -d "$blastdbname"  -j 3 -h 0.001 -i $seqfile -o $blastfile -m 9 -Q $pssmfile -C $chkfile
    if [ ! -s $chkfile ]; then
        return
    fi
    if [ $isPrintVerboseInfo -eq 1 ]; then
        echo "$BINPATH/chk2mtx.sh -d $outpath $chkfile $seqfile -nc 2> $errFile"
        $BINPATH/chk2mtx.sh -d $outpath $chkfile $seqfile -nc  2> $errFile
    else 
        $BINPATH/chk2mtx.sh -d $outpath $chkfile $seqfile -nc -q 2> $errFile
    fi
    st=`CheckErrMsg $errFile`
    if [ $st != 0 ]; then
        return
    fi
    # output the modm file in ASCII format with floating type, for
    # the program create_svm_vector
    mtxfile=$outpath/${rootname}.mtx
    if [ $isPrintVerboseInfo -eq 1 ]; then
        echo "$BINPATH/mtx2modm -d $outpath --pssm $mtxfile 1> /dev/null 2> $errFile "
    fi
    $BINPATH/mtx2modm -d $outpath --pssm $mtxfile 1> /dev/null 2> $errFile 
    st=`CheckErrMsg $errFile`
    if [ $st != 0 ]; then
        return
    fi
    modmfile=$outpath/${rootname}.modm

    if [ -s $modmfile -a -s $pssmfile ]; then 
        echo $modmfile >> $modmfilelist
        echo $pssmfile >> $pssmfilelist
        echo $seqfile >> $seqfilelist
    fi
}
#}}}
function BuildProfileFromPSSM() #pssmfile outpath seqfilelist pssmfilelist modmfilelist#{{{
{
    local pssmfile=$1
    local outpath=$2
    local seqfilelist=$3
    local pssmfilelist=$4
    local modmfilelist=$5
    local basename=`basename $pssmfile`
    local rootname=${basename%.*}
    if [ $isPrintVerboseInfo -eq 1 ]; then
        echo "$BINPATH/pssm2modm $pssmfile -d $outpath -t log --typeprofile 2 1> /dev/null 2> $errFile"
    fi 
    $BINPATH/pssm2modm $pssmfile -d $outpath -t log --typeprofile 2 1> /dev/null 2> $errFile
    st=`CheckErrMsg $errFile`
    if [ $st != 0 ]; then
        return
    fi
    modmfile=$outpath/${rootname}.modm
    seqfile=$outpath/${rootname}.aa
    $BINPATH/pssm2fasta.sh $pssmfile -o $seqfile
    if [ -s $seqfile -a $modmfile ]; then
        echo $modmfile >> $modmfilelist
        echo $pssmfile >> $pssmfilelist
        echo $seqfile >> $seqfilelist
    fi
}
#}}}
function Predzinc() # $outpath $seqfilelist $pssmfilelist modmfilelist predzincParaFile#{{{
{
    local outpath=$1
    local seqfilelist=$2
    local pssmfilelist=$3
    local modmfilelist=$4
    local predzincParaFile=$5 
    local basename=`basename $seqfilelist`
    local rootname=${basename%.*}
    #--------------------------------------------------------------------
    #SVM prediction using single-site vectors#{{{
    echo "SVM prediction using single-site vectors..."
    local resListArray=(C  H  D  E)
    local resNameListArray=(Cysteins  Histidines  Aspartates  Glutamates)
    for (( i=0; i <4; i++)); do
        local resList=${resListArray[$i]}
        local resName=${resNameListArray[$i]}
        local train_single_zn=$DATADIR/$dbname/train_single_${resList}
        local test_single_zn=$TMPDIR/${rootname}_single_${resList}
        #create vectors for svm prediction
        if [ $isPrintVerboseInfo -eq 1 ]; then
            echo "$BINPATH/create_svm_vector --mode pred -ns 1 --not-use-shcr -a $resList -s -1 -c -1 -l $modmfilelist --vector $test_single_zn.vector"
        fi 
        $BINPATH/create_svm_vector --mode pred -ns 1 --not-use-shcr -a $resList -s -1 -c -1 -l $modmfilelist --vector $test_single_zn.vector
        if [ -s $test_single_zn.vector ]; then 
            $gistbin/gist-classify -train $train_single_zn.vector -learned $train_single_zn.weight -test $test_single_zn.vector  1> $test_single_zn.predict  2> $errFile
        else
            echo "vector file empty for $resName" >&2
        fi
    done

    local test_single_zn=$TMPDIR/${rootname}_single_zn
    cat $TMPDIR/${rootname}_single_[CHDE].predict | sort -t\; -k1,1 -k5,5n > $test_single_zn.predict

    #--------------------------------------------------------------------
    #SVM prediction using pair-based vectors
    echo "SVM prediction using pair-based vectors..."
    local resList=CHDE
    local train_pair_zn=$DATADIR/$dbname/train_pair_${resList}
    local test_pair_zn=$TMPDIR/${rootname}_pair
    local metalTranMatrixDir=$DATADIR/$dbname/metalTranMatrix

    if [ $isPrintVerboseInfo -eq 1 ]; then
        echo "$BINPATH/create_svm_vector -ns 2  -s -1 --mode pred -a $resList -l $modmfilelist --vector $test_pair_zn.vector --metaltrans $metalTranMatrixDir"
    fi
    $BINPATH/create_svm_vector -ns 2  -s -1 -c 7.0 --win-pair 100 --win-min-pair 20 --win3 150 --win4 300 --mode pred -a $resList -l $modmfilelist --vector $test_pair_zn.vector --metaltrans $metalTranMatrixDir

    if [ -s $test_pair_zn.vector ]; then
        #carry out prediction on pair-based vectors by gist svm program
        if [ $isPrintVerboseInfo -eq 1 ]; then
            echo "$gistbin/gist-classify -train $train_pair_zn.vector   -learned $train_pair_zn.weight -test $test_pair_zn.vector   1> $test_pair_zn.predict 2> $errFile"
        fi
        $gistbin/gist-classify -train $train_pair_zn.vector   -learned $train_pair_zn.weight -test $test_pair_zn.vector   1> $test_pair_zn.predict    2> $errFile
    else 
        echo "vector file empty for pair vector" >&2
    fi

    #--------------------------------------------------------------------
    local test_gate_zn=$TMPDIR/${rootname}_gate_zn
    if [ -f $test_single_zn.predict -a -f $test_pair_zn.predict ]; then
        #using gating network 
        if [ $isPrintVerboseInfo -eq 1 ] ; then 
            echo "$BINPATH/gating_gistPred -i $test_single_zn.predict $test_pair_zn.predict -o $test_gate_zn.predict"
        fi
        $BINPATH/gating_gistPred -i $test_single_zn.predict $test_pair_zn.predict -o $test_gate_zn.predict
    fi 
    #--------------------------------------------------------------------#}}}

    #Homology-based prediction
    echo "homology-based prediction..."
    local fragMatchFile=$TMPDIR/${rootname}.fragbin
    local bkFreqFile=$DATADIR/$dbname/background_freq.txt
    mkdir -p $TMPDIR/Qij
    $BINPATH/pssm2Qij -list $pssmfilelist -d $TMPDIR/Qij --wb --bkfile $bkFreqFile --sad -a AVLIPFMKRHGSTCYNEWDQ --newseq --idtype 1 --typeprofile 1 1> /dev/null 2> $errFile
    st=`CheckErrMsg $errFile`
    # for the test chain, the Qij files are used as modm and fragacc bin, so that
    # no matter what the merging ratio is, the Qij matrix is used. 2009-04-07 
    mkdir -p $TMPDIR/modm $TMPDIR/fragacc
    for file in $(find $TMPDIR/Qij -name "*.Qijbin");do
        local bname1=`basename $file`
        local rtname1=${bname1%.*}
        ln -s $file $TMPDIR/modm/${rtname1}.modmbin
        ln -s $file $TMPDIR/fragacc/${rtname1}.fragaccbin
    done

    local trainIDListFile=$DATADIR/$dbname/train.idlist
    local testIDListFile=$TMPDIR/${rootname}.idlist
    find $TMPDIR/Qij/ -name "*.Qijbin" | $BINPATH/rootname > $testIDListFile

    #local trainQijPath=$DATADIR/$dbname/Qij
    local trainQijPath=$DATADIR/$dbname/train_Qij
    local testQijPath=$TMPDIR/Qij
    local qijformat=1
    #local frag_acc_path=$DATADIR/$dbname/fragacc/
    local frag_acc_path=$DATADIR/$dbname/train_fragacc
    local test_frag_acc_path=$TMPDIR/fragacc
    #local modmpath=$DATADIR/$dbname/modm
    local modmpath=$DATADIR/$dbname/train_modm
    local test_modmpath=$TMPDIR/modm
    local modmformat=1
    local bkfreqfile=$DATADIR/$dbname/background_freq.txt
    local parafile=$DATADIR/$dbname/control_parameters.txt
    local resultpath=$TMPDIR/fragmatch
    local dbtype=1
    mkdir -p $resultpath
    if [ $isPrintVerboseInfo -eq 1 ] ; then 
        echo "$BINPATH/search_new -dbtype $dbtype --rb --wb --train $trainIDListFile --test $testIDListFile --train-qij $trainQijPath --test-qij $testQijPath --qijformat $qijformat --fragacc $frag_acc_path --modm $modmpath --test-modm $test_modmpath --test-fragacc $test_frag_acc_path --modmformat $modmformat --bkfile $bkfreqfile --para $parafile --result $resultpath 1> /dev/null 2> $errFile"
    fi
    $BINPATH/search_new -dbtype $dbtype --rb --wb --train $trainIDListFile --test $testIDListFile --train-qij $trainQijPath --test-qij $testQijPath --qijformat $qijformat --fragacc $frag_acc_path --modm $modmpath --test-modm $test_modmpath --test-fragacc $test_frag_acc_path --modmformat $modmformat --bkfile $bkfreqfile --para $parafile --result $resultpath 1> /dev/null 2> $errFile
    st=`CheckErrMsg $errFile`

    local trainSeqLenFile=$DATADIR/$dbname/train.seqlen
    local testSeqLenFile=$TMPDIR/${rootname}.seqlen
    local fragmatchfilelist=$TMPDIR/${rootname}.fragmatchfilelist

    find $resultpath -name "*.fragbin" > $fragmatchfilelist
    local dumpedfastaseq=$TMPDIR/dumpedseq.fa
    cat `cat $seqfilelist`  > $dumpedfastaseq
    $BINPATH/getseqlen.py $dumpedfastaseq -i $dumpedfastaseq -o $testSeqLenFile

    if [ $isPrintVerboseInfo -eq 1 ]; then
        echo "$BINPATH/postscan-frag-search --outpath $resultpath -l $fragmatchfilelist --rb --len-file $trainSeqLenFile  --test-len-file $testSeqLenFile 1> /dev/null 2> $errFile"
    fi
    $BINPATH/postscan-frag-search --outpath $resultpath -l $fragmatchfilelist --rb --len-file $trainSeqLenFile  --test-len-file $testSeqLenFile 1> /dev/null 2> $errFile

    st=`CheckErrMsg $errFile`

    local proportion=0.45
    local max_homopro=3
    local resList=CHDE
    local cutoff_homoScore=6
    local weightOnHomo=1.0
    local homo_search_new=$resultpath
    local metalProFile=$DATADIR/$dbname/revised-closeMetalPro.dat
    local ssbondProFile=$DATADIR/$dbname/train.ssbond

    local test_final_zn=$TMPDIR/${rootname}_final_zn
    #normalize the score to 0~1

    if [ $isPrintVerboseInfo -eq 1 ]; then
        echo "$BINPATH/znpred-postscan --metal $metalProFile --ssbond $ssbondProFile -a $resList --homo $homo_search_new  --gate  $test_gate_zn.predict --homo-format 1 -o $test_final_zn.predict > /dev/null 2> $errFile"
    fi
    $BINPATH/znpred-postscan --metal $metalProFile --ssbond $ssbondProFile -a $resList --homo $homo_search_new  --gate  $test_gate_zn.predict --homo-format 1 -o $test_final_zn.predict > /dev/null 2> $errFile 
    st=`CheckErrMsg $errFile`

    # report the results of prediction
    if [ $isPrintVerboseInfo -eq 1 ]; then
        echo "$BINPATH/report_zincpred $test_final_zn.predict  -seqfile $dumpedfastaseq --ishtml yes --parafile $predzincParaFile -o $test_final_zn.report 1> /dev/null 2> $errFile"
    fi
    $BINPATH/report_zincpred $test_final_zn.predict -seqfile $dumpedfastaseq --ishtml yes --parafile $predzincParaFile -o $test_final_zn.report 1> /dev/null 2> $errFile
    st=`CheckErrMsg $errFile`

    if [ -f $test_final_zn.report  -a -f $test_final_zn.report.html ] ; then 
        /bin/cp -f $test_final_zn.predict $outpath/$outname.predzinc.predict
        /bin/cp -f $test_final_zn.report $outpath/$outname.predzinc.report
        /bin/cp -f $test_final_zn.report.html $outpath/$outname.predzinc.report.html
        if [ -d $TMPDIR/psiblast ]; then 
            #/bin/cp -f $TMPDIR/psiblast/*.mtx $outpath
            /bin/cp -f $TMPDIR/psiblast/*.pssm $outpath
            /bin/cp -f $TMPDIR/psiblast/*.blast $outpath
        fi
        echo "PredZinc succeeded. The results have been output to
        $outpath/$outname.predzinc.report
        $outpath/$outname.predzinc.report.html"
    else
        echo "PredZinc failed. No results have been output. Check the error message above."
    fi
}
#}}}

if [ $# -lt 1 ]; then
#    echo "too few arguments" >&2
    PrintHelp
    exit 1
fi

TMPDIR=
isClean=1
outpath=./
outname=query
seqFile=
pssmFile=
pssmFileListFile=
numCPU=1
dbname=passe

isPrintVerboseInfo=0

cntInputFileArgument=0 # count the number of input file argument supplied, note that only one of the four optional argument (seqFile, pssmFile, seqFileListFile, pssmFileListFile) can be supplied

isNonOptionArg=false
while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        idList="$idList $1"
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h|--help) PrintHelp; exit 0;;
            -showexample|--showexample) echo "$examples"; exit 0;;
            -version|--version) PrintVersion; exit 1;;
            -pssm|--pssm) pssmFile=$2;shift;;
            -pssmfilelist|--pssmfilelist) pssmFileListFile=$2;shift;;
            -blastdb|--blastdb) blastdbname=$2;shift;;
            -db|--db|-dbname|--dbname) dbname=$2;shift;;
            -outpath|--outpath) outpath=$2;shift;;
            -outname|--outname) outname=$2;shift;;
            -cpu|--cpu) numCPU=$2;shift;;
            -nc|-not-clean|--not-clean) isClean=0;;
            -verbose|--verbose) isPrintVerboseInfo=1;;
            -tmpdir|--tmpdir) TMPDIR=$2;shift;;
            -q|-quiet|--quiet) isQuiet=true;;
            -*) echo "Error! Wrong argument: $1"; exit 20;;
        esac
    else
        seqFile=$1
    fi
    shift
done


#check for necessary programs and files
#-----------------------------------------------------------
IsProgExist $blastbin/blastpgp
IsPathExist $DATADIR
#-----------------------------------------------------------

if [ $isPrintVerboseInfo -eq 1 ]; then
    echo "BLASTMAT=$BLASTMAT"
    echo "BLASTBIN=$BLASTBIN"
    echo "BLASTDB=$BLASTDB"
fi

mkdir -p $outpath

# create the temporary dirs
if [ "$TMPDIR" == "" ]; then
    TMPDIR=$(mktemp -d /tmp/tmpdir.predzinc.XXXXXXXXX) || { echo "Failed to create temp dir" >&2; exit 1; }   
else
    mkdir -p $TMPDIR
fi
errFile=$TMPDIR/predzinc.err

tmpSeqFileListFile0=$TMPDIR/query.0.seqfilelist
tmpPSSMFileListFile0=$TMPDIR/query.0.pssmfilelist
tmpMODMFileListFile0=$TMPDIR/query.0.modmfilelist

tmpSeqFileListFile1=$TMPDIR/query.1.seqfilelist
tmpPSSMFileListFile1=$TMPDIR/query.1.pssmfilelist
tmpMODMFileListFile1=$TMPDIR/query.1.modmfilelist

# create predzincParaFile, parameters
blast_version=`$blastbin/blastpgp - | awk '/^blastpgp/{print $2}'`
trainingSet="PDB20 (April 2004), 2727 unique chains"
blastdb_para=`$blastbin/blastpgp -i $DATADIR/1letterseq.fa -d $blastdbname 2>/dev/null| grep "^Database\|sequences.*letters" | sed 's/Database//g' | tr '\n' ' '`
predzincParaFile=$TMPDIR/predzincparafile.txt
echo "Training set: $trainingSet" >$predzincParaFile
echo "blastpgp version: $blast_version" >>$predzincParaFile
echo "blast database: $blastdb_para" >>$predzincParaFile

OS=`uname -o`
case $OS in 
	Cygwin*)
	blastdbname=`cygpath -w "$blastdbname" `
	export BLASTMAT=`cygpath -w "$BLASTMAT"`
	;;
esac


if [  "$seqFile" != "" ] ; then #the input is aa sequence
    if [ -s "$seqFile" ]; then
        numSeq=`grep "^>" $seqFile | wc -l`
        if [ $numSeq -gt 0 ]; then 
            echo "$numSeq sequences recognized in the input file $seqFile"
            $BINPATH/splitfasta.py $seqFile -outpath $TMPDIR/seq -q
            cntValidSeq=0
            for file in $(find $TMPDIR/seq -name "*.aa"); do
                v=`CheckSequence $file $SEQ_MODE_AASEQ` 
                if [ $v -eq 0 ]; then 
                    ((cntValidSeq++))
                    echo "$file" >> $tmpSeqFileListFile0 
                fi
            done

            if [ $cntValidSeq -gt 0 ]; then

                if [ -f "$blastdbname.phr"  -o -f "$blastdbname.00.phr" ]; then 
                    echo "Using blastdb under current directory."
                    export BLASTDB=$PWD
                elif [ -f "$DATADIR/blastdb/$blastdbname.phr"  -o  -f "$DATADIR/blastdb/$blastdbname.00.phr" ]; then 
                    echo "Using blastdb at $DATADIR/blastdb."
                    export BLASTDB=$DATADIR/blastdb
                elif [ "$BLASTDB" != "" -a  \( -f "$BLASTDB/$blastdbname.phr" -o  -f "$BLASTDB/$blastdbname.00.phr" \) ]; then 
                    echo "Env BLASTDB is set and using blastdb at $BLASTDB/$blastdbname."
                else
                    echo "Fatal! Could not find blastdb $blastdbname" >&2 
                    exit
                fi

                echo "$cntValidSeq sequences satisfy the requirement and will be processed further..."
                echo "Running PSI-BLAST for $cntValidSeq sequences."
                echo "Be patient with PSI-BLAST for building sequence profiles..."
                echo
                mkdir -p $TMPDIR/psiblast
                for file in $(cat $tmpSeqFileListFile0); do
                    BuildProfileFromAASeq $file $TMPDIR/psiblast $tmpSeqFileListFile1 $tmpPSSMFileListFile1 $tmpMODMFileListFile1 
                done
            else
                echo "No sequences satisfy the requirement. Ignore..."
            fi
        else
            echo "Sequence file $seqFile format error (not Fasta). Ignore!" >&2
        fi
    else
        echo "Sequence file $seqFile does not exists or is empty. Ignore!" >&2
    fi
fi

if [ "$pssmFile" != ""  ]; then
    if [ ! -s "$pssmFile"  ]; then
        echo "PSSM file $pssmFile does not exist or empty. Ignore."  >&2 
    else 
        echo "1 PSSM file \"$pssmFile\" recognized."
        echo $pssmFile >> $tmpPSSMFileListFile0
    fi
fi

if [ "$pssmFileListFile" != "" ] ; then
    if [ ! -s "$pssmFileListFile"  ]; then
        echo "PSSM list file $pssmFileListFile does not exist or empty. Ignore."  >&2 
    else 
        cntfile=0
        for file in $(cat $pssmFileListFile); do 
            echo $file >> $tmpPSSMFileListFile0
            ((cntfile++))
        done 
        echo "$cntfile PSSM files recognized from the list file \"$pssmFileListFile\"."
    fi
fi

if [ -s $tmpPSSMFileListFile0 ]; then
    for file in $(cat $tmpPSSMFileListFile0); do
        v=`CheckSequence $file $SEQ_MODE_PSSM` 
        cntValidSeq=0
        mkdir -p $TMPDIR/pssm
        if [ $v == 0 ]; then 
            ((cntValidSeq++))
            BuildProfileFromPSSM $file $TMPDIR/pssm $tmpSeqFileListFile1 $tmpPSSMFileListFile1 $tmpMODMFileListFile1 
        fi
    done
fi

NFile=`cat $tmpPSSMFileListFile1 | wc -l`

if [ $NFile -le 0 ] ; then 
    echo "Valid number of sequences is zero. Exit." >&2
    exit 1
fi
echo "Run PredZinc for $NFile sequence(s)"
echo
Predzinc $outpath $tmpSeqFileListFile1 $tmpPSSMFileListFile1 $tmpMODMFileListFile1 $predzincParaFile

#clean up intermediate output files
if [ $isClean -eq 1  ]; then
    rm -rf $TMPDIR
else
    echo "Temporary files are kept at $TMPDIR"
fi
