#!/usr/bin/python
# split multiple sequences in the fasta file in individual files
# one sequence per file
import sys,re,os;

BLOCK_SIZE=100000;

usage="""
Usage:  splitfasta.py [Options] [-i] fastafile
Options:
  -i       <file> : input file
  -outpath <dir>  : output the result to dir, default=./
  -ext     <str>  : set the file extension, default = aa
  -q              : quiet mode, do not write any messages
  -nameseq        : name the splitted files by $rootname_i, i = 0,1,2...
                  : this is useful when it is the filename extracted from the 
                  : annotation line having a duplicated name
  -h|--help       : print this help message and exit

Created 2010-10-22, updated 2010-10-22, Nanjiang

Examples:
    splitfasta.py example.fasta 
    splitfasta.py example.fasta -outpath outdir
"""

def PrintHelp():
    print usage;

def GetSeqIDFromAnnotation(line):#{{{
# get the ID from the annotation line of the fasta  file
# updated 2010-08-24
    seqID = "";
    line = line.lstrip('>').split()[0]; #get the first word after '>'
    if line and line.find('|') >= 0:# if the annotation line has |, e.g. >sp|P0AE31|ARTM_ECOL6 Arginine ABC transporter permease
        strs = line.split("|");
        if (strs[0] == "sp" or strs[0] == "lcl" or strs[0] == "tr") : seqID = strs[1];
        else                 : seqID = strs[0];
    else:
        seqID=line;
    return seqID;#}}}

def OutputSplittedSeq(seqWithAnno,cntSeq, rootname, outpath):#{{{
    begseq=seqWithAnno.find("\n");
    seq=seqWithAnno[begseq:];
    seqID=GetSeqIDFromAnnotation(seqWithAnno[0:begseq]);
#   seq=seq.replace('\n','').replace(' ','');
#   length=len(seqWithAnno[seqWithAnno.find("\n"):].replace('\n','').replace(' ',''));
    if seqID == "" or isNameFileSequentially:
        seqID= rootname + "_%d"%cntSeq;
    outfile=outpath+os.sep+ seqID + "." + extension;
    try:
        fpout = open(outfile,"w");
        fpout.write("%s"%seqWithAnno[0:begseq]);
        fpout.write("%s\n" % seq);
        fpout.close();
        if not isQuiet:
            print >> sys.stdout, "%d\t%s output"%(cntSeq, outfile);
    except:
        print >> sys.stderr ,"Failed to open %s for write"%outfile;
        exit(1);


#}}}
def SplitFasta(inFile, outpath):#{{{
# The faster version
    rootname=os.path.basename(os.path.splitext(inFile)[0]);
    isFirstSeq=True;
    cntSeq = 0;
    fpin = open(inFile, "r");
    buff = fpin.read(BLOCK_SIZE);
    brokenSeqWithAnnoLine=""; ##for the annotation line broken by BLOCK read
    while buff:
        beg=0;
        end=0;
        while 1:
            if brokenSeqWithAnnoLine:
                if brokenSeqWithAnnoLine[len(brokenSeqWithAnnoLine)-1] == "\n":
                    end=buff.find(">");
                else:
                    end=buff.find("\n>");
                if end >= 0:
                    seqWithAnno = brokenSeqWithAnnoLine + buff[0:end];
                    OutputSplittedSeq(seqWithAnno,cntSeq, rootname, outpath);
                    brokenSeqWithAnnoLine = "";
                    cntSeq += 1;
                    beg=end;
                else:
                    brokenSeqWithAnnoLine += buff;
                    break;

            beg=buff.find(">",beg);
            end=buff.find("\n>",beg+1);
            if beg >= 0:
                if end >=0:
                    seqWithAnno=buff[beg:end];
                    OutputSplittedSeq(seqWithAnno, cntSeq, rootname, outpath);
                    cntSeq += 1;
                    beg=end;
                else:
                    brokenSeqWithAnnoLine=buff[beg:];
                    break;
            else:
                break;

        buff = fpin.read(BLOCK_SIZE);
    fpin.close();
    if brokenSeqWithAnnoLine:
        seqWithAnno=brokenSeqWithAnnoLine;
        OutputSplittedSeq(seqWithAnno, cntSeq, rootname, outpath);
        cntSeq += 1;
        
    return cntSeq;

#}}}

if __name__ == '__main__' :
    # Check argv
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp();
        sys.exit(1);

    outpath="./";
    inFile="";
    extension="aa";
    isQuiet=False;
    isNameFileSequentially=False;

    i = 1;
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            isNonOptionArg=False;
            i = i + 1;
        elif sys.argv[i] == "--":
            isNonOptionArg=True;
            i = i + 1;
        elif sys.argv[i][0] == "-":
            if sys.argv[i] ==  "-h" or  sys.argv[i] == "--help":
                PrintHelp();
                sys.exit(0);
            elif sys.argv[i] == "-i" or sys.argv[i] == "--infile":
                inFile=sys.argv[i+1];
                i = i + 2;
            elif sys.argv[i] == "-ext" or sys.argv[i] == "--ext":
                extension=(sys.argv[i+1]);
                i = i + 2;
            elif sys.argv[i] == "-outpath" or sys.argv[i] == "--outpath":
                outpath=sys.argv[i+1];
                i = i + 2;
            elif sys.argv[i] == "-q" or sys.argv[i] == "--q" or sys.argv[i] == "--quiet":
                isQuiet=True;
                i = i + 1;
            elif sys.argv[i] == "-nameseq" or sys.argv[i] == "--nameseq":
                isNameFileSequentially=True;
                i = i + 1;
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i]);
                sys.exit(1);
        else:
            inFile=sys.argv[i];
            i+=1;
           

    if inFile == "":
        print >> sys.stderr,"Error! Input file not set.";

    os.system("mkdir -p %s"%outpath);

    try :
        numSeq = SplitFasta(inFile,outpath);
        if not isQuiet:
            print >> sys.stdout,"%d sequences are splitted and output to %s"%(numSeq, outpath);

    except :
        print >>sys.stderr, "except for the input file: %s" % inFile;
        raise ;
