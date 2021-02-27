#!/usr/bin/env python
import sys,os;

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
def ReadSingleFasta(inFile):#{{{
# return seqID, annotation, aaSeq
# the leading '>' of the annotation is removed
    seqID="";
    aaSeq="";
    annotation="";
    try:
        fpin = open(inFile, "r");
        lines = fpin.readlines();
        fpin.close();
        for line in lines:
            if line[0] == ">":
                seqID = GetSeqIDFromAnnotation(line);
                annotation = line.lstrip(">").strip();
            else:
                aaSeq+=line.strip();
    except IOError: 
        print >> sys.stderr, "Except for the input file ", inFile, "in the function ReadSingleFasta";
    return (seqID, annotation, aaSeq);
#}}}
def ReadFasta(inFile):#{{{
#read in all fasta sequences, just the original sequence, do not varify
#return idList, annotationList, seqList
    idList=[];
    annotationList=[];
    seqList=[];

    fpin = open(inFile, "r");
    lines = fpin.readlines();
    fpin.close();
    i=0;
    while i < len(lines):
        line = lines[i];
        if line[0] == ">":
            seqID = GetSeqIDFromAnnotation(line);
            annoLine = line.lstrip('>').strip();
            seq="";
            i += 1;
            while i < len(lines) and lines[i][0] != '>':
                seq+=lines[i].strip();
                i=i+1;
            idList.append(seqID);
            annotationList.append(annoLine);
            seqList.append(seq);

    return (idList, annotationList, seqList);
#}}}
def coverage(a1,b1,a2,b2):#{{{
# return the coverage of two intervals
# a1, b1, a2, b2 are integers
# when the return value <=0, it means there is no coverage
    return (min(b1,b2)-max(a1,a2));
#}}}

def isnumeric(s):#{{{
# determine whether a string is numeric value
    try:
        i = float(s);
        return True;
    except ValueError, TypeError:
        return False;
#}}}
def isnumeric_extended(lit):#{{{
#including complex, hex, binary and octal numeric literals
    # Handle '0'
    if lit == '0': 
        return True;
    # Hex/Binary
    litneg = lit[1:] if lit[0] == '-' else lit
    if litneg[0] == '0':
        if litneg[1] in 'xX':
            try:
                v=int(lit,16)
                return True;
            except ValueError, TypeError:
                return False;
        elif litneg[1] in 'bB':
            try:
                v=int(lit,2);
                return True;
            except ValueError, TypeError:
                return False;
        else:
            try:
                v=int(lit,8)
                return True;
            except ValueError:
                pass
 
    # Int/Float/Complex
    try:
        v=int(lit)
        return True;
    except ValueError:
        pass
    try:
        v=float(lit)
        return True;
    except ValueError:
        pass
    try:
        v=complex(lit)
        return True;
    except ValueError:
        return False;
#}}}

def confirm(prompt=None, resp=False):#{{{
#     """prompts for yes or no response from the user. Returns True for yes and
#     False for no.
# 
#     'resp' should be set to the default value assumed by the caller when
#     user simply types ENTER.
# Example :
# if confirm(prompt='Create Dir', resp=True) == True:
#     print "run"
# else:
#     print "ignore"
#     """
    if prompt is None:
        prompt = 'Confirm'

    if resp:
        prompt = '%s [%s]|%s: ' % (prompt, 'y', 'n')
    else:
        prompt = '%s [%s]|%s: ' % (prompt, 'n', 'y')

    while True:
        ans = raw_input(prompt)
        if not ans:
            return resp
        if ans not in ['y', 'Y', 'n', 'N']:
            print 'please enter y or n.'
            continue
        if ans == 'y' or ans == 'Y':
            return True
        if ans == 'n' or ans == 'N':
            return False
#}}}
