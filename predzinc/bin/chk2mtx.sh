#!/bin/bash 
# generate the mtx file from chk of the blastpgp output by makemat program

usage="
usage: chk2mtx.sh [options] chkfile seqfile
  shell script for converting chkfile to mtx file make sure that the
  chkfile and seqfile have the same rootname. not batch mode

  one time one file

options:
  -d DIR    Output directory (default: ./)
  -S FLOAT  Scaling factor to avoid round-off problems (default: 100.0)
  -q        quiet mode, do not print any error message
  -nc       not clean

Created 2007-01-17, updated 2011-09-27, Nanjiang Shu
"

function PrintHelp()
{
    echo "$usage"
}
function AddAbsolutePath()#{{{
{
    local oldPath=$1
    if [ "${oldPath:0:1}" == "/" ]; then
        newPath=$oldPath
    else
        newPath=$PWD/$oldPath
    fi
    echo "$newPath"
}
#}}}
if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

BINPATH=`dirname $0`
blastbin=$BINPATH
export BLASTMAT=$BINPATH/../data
chkfile=
seqfile=
OUTDIR=./
isQuietMode=false
scale=100.0
isClean=1


while [ "$1" != "" ]; do
    case $1 in 
        -h|--help) PrintHelp; exit;;
        -d) OUTDIR=$2; shift;;
        -S) scale=$2;shift;;
        -q) isQuietMode=true;;
        -nc) isClean=0;;
        *) 
        if [ "$chkfile" == "" ]; then
            chkfile=$1 
        else
            seqfile=$1
        fi
    esac
    shift
done

if [ ! -d $OUTDIR ]; then
    mkdir $OUTDIR
fi

chkfile=`AddAbsolutePath $chkfile` 
seqfile=`AddAbsolutePath $seqfile`
blastbin=`AddAbsolutePath $blastbin` 

basename=`basename "$chkfile"`
rootname=${basename%.*}
filepath=`dirname $chkfile`

snfile=$OUTDIR/$rootname.sn
pnfile=$OUTDIR/$rootname.pn

echo "$seqfile" > $snfile
echo "$chkfile" > $pnfile

currDir=$PWD

cd $OUTDIR
if [ "$isQuietMode" == "true" ]; then
    $blastbin/makemat -S $scale -P $rootname > /dev/null 2>&1
else
    echo "$blastbin/makemat -S $scale -P $rootname"
    $blastbin/makemat -S $scale -P $rootname
    echo -e "$OUTDIR/$rootname.mtx output"
fi

if [ "$OUTDIR" != "$filepath" ]; then
    /bin/cp -f $filepath/$rootname.mtx $OUTDIR
fi

# clean
if [ $isClean -eq 1 ]; then 
    rm -f $snfile
    rm -f $pnfile
    rm -f $rootname.aux
    rm -f $rootname.mn
fi

cd $currDir
