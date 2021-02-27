#!/bin/bash
# convert the ascii format modmfile to binary format, the output file name is
# $filename+bin

function PrintHelp()
{
    echo "Usage: txtmodm2binary.sh [options] txtmodmfile | -l modmfilelist-file"
    echo "  convert the ascii format modmfile to the binary format, plus bin"
    echo "  to the end of the original filename"
    echo "options:"
    echo "  --sad          : indicating the modm file containing SAD information"
    echo "  -a str         : supply the forced alphabet"
    echo "  --outpath path : setting the output path, default outpath is the"
    echo "                 : same as the original modmfile"
    echo "  -h|--help      : print this help message and exit"
    echo 
    echo " Created 2007-11-03, last modified 2007-11-03, Nanjiang Shu"
}

if [ $# -lt 1 ]; then
    PrintHelp
    exit
fi

outpath=
isHaveSAD=false
modmfilelist_file=
rtnameprg=`echo $0 | sed 's/.*\///g' | sed 's/\.[^\.]*$//g'`

TMPFILE=`mktemp  /tmp/$rtnameprg.XXXXXX`
if [ $? -ne 0 ]; then
    echo "$0: Can't create temp file, exiting..."
    exit 1
fi

((i=0))
while [ "$1" != "" ]; do 
    case $1 in 
        -h|--help) PrintHelp; exit;;
        --sad) isHaveSAD=true;;
        -a) alphabet=$2;shift;;
        --outpath) outpath=$2;shift;;
        -l)modmfilelist_file=$2;shift;;
        *) echo $1>> $TMPFILE
    esac
    shift
done

if [ "$modmfilelist_file" == "" ];   then
    modmfilelist_file=$TMPFILE
fi

BASECMD="$BINPATH/test-readmodm --mode 0 -b --notime -l $modmfilelist_file"
APPCMD=

if [ "$alphabet" != "" ];then
    APPCMD="$APPCMD -a $alphabet"
fi
if [ "$outpath" != "" ]; then
    APPCMD="$APPCMD --outpath $outpath"
fi
if [ "$isHaveSAD" == "true" ]; then
    APPCMD="$APPCMD --sad"
fi

CMD="$BASECMD $APPCMD"
$CMD

