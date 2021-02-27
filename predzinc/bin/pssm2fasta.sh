#!/bin/bash
#extract fasta sequences from pssm files

usage="
usage:   pssm2fasta.sh FILE [FILE ...]

Extract sequences in Fasta format from PSSM files output by blastpgp

Options:
  -o     FILE  output to file, default to stdout
  -l     FILE  set the idListFile
  -q           quiet mode
  -h|--help    print this help message and exit

Created 2009-05-26, updated 2011-10-11, Nanjiang

Examples:
  pssm2fasta.sh test.pssm
"
function PrintHelp()
{
    echo "$usage"
}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

isQuiet=false
outfile=
idListFile=
idList=

function PSSM2Fasta() #$file#{{{
{
    local file=$1
    local basename=`basename "$1"`
    local rootname=${basename%.*}
    echo ">$rootname"
    awk '{if (NF>40) printf "%c", $2} END{printf "\n"}' $file 
}
#}}}
isNonOptionArg=false
while [ "$1" != "" ]; do
    if [ "$isNonOptionArg" == "true" ]; then 
        idList="$idList $1"
        isNonOptionArg=false
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) PrintHelp; exit;;
            -o|--o|-outfile|--outfile) outfile=$2;shift;;
            -l|--l|--list) idListFile=$2;shift;;
            -q|--q|--quiet) isQuiet=true;;
            -f|--outformat) outFormat=$2;shift;;
            -*) echo "Error! Wrong argument: $1" >&2; exit;;
        esac
    else
        idList="$idList $1"
    fi
    shift
done

tmpfile=$(mktemp /tmp/tmp.pssm2fasta.XXXXXXXXX) || { echo "Failed to create temp file" >&2; exit 1; }  


(
for file in ${idList}; do 
    PSSM2Fasta $file
done

if [ -f "$idListFile" ]; then
    for file in $(cat $idListFile ) ; do
        PSSM2Fasta $file
    done
fi
) >$tmpfile

if [ "$outfile" != "" ]; then 
    /bin/mv -f $tmpfile $outfile
else 
    /bin/cat $tmpfile
fi

/bin/rm -f $tmpfile
