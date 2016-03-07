/*
 * =====================================================================================
 * 
 *        Filename:  postscan-frag-search.cpp
 *     Description:  calculate the homology score of proteins by scanning the
 *     			     dot plots
 *         Version:  1.0
 *         Created:  11/06/2007 05:01:43 PM CET
 *        Compiler:  g++
 *          Author:  Nanjiang Shu (Shu), nanjiang@struc.su.se
 *         Company:  Structural Chemistry, Stockholm Univesity
 * 
 * =====================================================================================
 */
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "myfunc.h"
#include "mytemplate.h"
#include "mypro.h"
#include "array.h"
using namespace std;

// Change Log/*{{{*/
/*****************************************************************************
 * ChangeLog 2007-07-16:
 * adding new algorithm for scanning the pair dot plot
 * 
 * ChangeLog 2007-10-11
 *      change the default sequence length list file to $DATADIR/pdbaa.seqlen,
 *      update the pdbaa.seqlen when update pdbaa
 *
 * ChangeLog 2007-10-26
 *      add options 
 *          --topn topN
 *          --cs cutoff_rawper
 *      set the size of id to a larger value
 * ChangeLog 2007-11-18
 *      add option
 *          --nm, the multidomain detection functionality can be turned off
 *      also, if the multidomain functionality is off, the matches with great
 *      different in size will not be considered as homology, that in, the
 *      score will be normalized by avgLength
 * ChangeLog 2008-01-15
 *      add weight for calculating the homology score
 *          scheme 1: using the fragment-fragment score output in the frag file
 *          scheme 2: defining different regions in the protien depend on the core
 *                    region of the family
 *      function for identifying homology domains, the matched dot must be at
 *      near the diagnol if the domain length are comparable
 * ChangeLog 2008-01-16
 *      add option: --ignore-non-diag
 *          for domains with similar length, ignore the non main diagnol
 *          matches. This option should not be used when multiple domains exist
 *          in the sequences.
 *      add option: --cfragscore
 * ChangeLog 2008-01-25
 *      On 2008-01-24, the frag file has been change to 5 columns, the original
 *      fifth column containing the index of sequnces has been removed.
 *      Therefore, the ReadInPair function should be modified to cope with the
 *      new format of frag file
 * ChangeLog 2009-06-09 
 *      the option --test-len-file added
 * ChangeLog 2009-06-25
 * 		add the option --rb
 * 		add the function to auto detect the types of frag file
 * ChangeLog 2010-01-06 
 * 		add the function --begin --end
 ****************************************************************************/
/*}}}*/
char changeLog[] = /*{{{*/\
"ChangeLog 2007-07-16:                                                           \n\
    adding new algorithm for scanning the pair dot plot                        \n\
                                                                                 \n\
ChangeLog 2007-10-11                                                            \n\
   change the default sequence length list file to $DATADIR/pdbaa.seqlen,     \n\
   update the pdbaa.seqlen when update pdbaa                                  \n\
                                                                                \n\
ChangeLog 2007-10-26                                                            \n\
   add options                                                                \n\
      --topn topN                                                            \n\
      --cs cutoff_rawper                                                     \n\
   set the size of id to a larger value                                       \n\
ChangeLog 2007-11-18                                                            \n\
   add option                                                                 \n\
      --nm, the multidomain detection functionality can be turned off        \n\
   also, if the multidomain functionality is off, the matches with great      \n\
   different in size will not be considered as homology, that in, the         \n\
   score will be normalized by avgLength                                      \n\
ChangeLog 2008-01-15                                                            \n\
   add weight for calculating the homology score                              \n\
      scheme 1: using the fragment-fragment score output in the frag file    \n\
      scheme 2: defining different regions in the protien depend on the core \n\
                region of the family                                         \n\
   function for identifying homology domains, the matched dot must be at      \n\
   near the diagnol if the domain length are comparable                       \n\
ChangeLog 2008-01-16                                                            \n\
   add option: --ignore-non-diag                                              \n\
      for domains with similar length, ignore the non main diagnol           \n\
      matches. This option should not be used when multiple domains exist    \n\
      in the sequences.                                                      \n\
   add option: --cfragscore                                                   \n\
ChangeLog 2008-01-25                                                            \n\
   On 2008-01-24, the frag file has been change to 5 columns, the original    \n\
   fifth column containing the index of sequnces has been removed.            \n\
   Therefore, the ReadInPair function should be modified to cope with the     \n\
   new format of frag file                                                    \n\
ChangeLog 2008-04-12                                                            \n\
   A memory leak fixed. The memory leak was caused when return the function    \n\
   before freed the memory. Valgrind checked on 2008-04-12.  \n\
      ";/*}}}*/

// Note/*{{{*/
/*****************************************************************************
 * Note 2008-01-16
 *  using weight scheme 1 get worse result than weight scheme 0 (no weight).
 *  Moreover, the more sophisticated weighting scheme considering the
 *  distribution of each frag file is worse than using a very simple empirical
 *  weighting deriving scheme, w = (fragscore - 45000) / 7500
 *
 ****************************************************************************/
///*}}}*/
#define EMPTY -1
#define FILLED 1000

#define LOWER_CORNER 0
#define UPPER_CORNER 1

#define DATATYPE_DOT_SCORE double

#ifndef FRAGFORMAT
#define FRAGFORMAT
#define FRAGFORMAT_TUPING 0  /*six columns*/
#define FRAGFORMAT_NANJIANG 1 /*five columns*/
#endif 

int fragformat = FRAGFORMAT_NANJIANG; /*default frag format 2009-06-25*/
bool isReadBinaryFragfile = false; /*set whether to read the binary frag file*/
double  cutoff_rawper = -1.0;
double  cutoff_postper = 5.0;
float cutoff_fragscore = 0.0; /*cutoff_fragscore should be derived from the disribution of fragscore*/
int topN = 20;
int scale_full = 16;
int scale_half = 8;
bool isFindingMultiDomain = true;/*whether to detect multi domains, this option is turned on by default. However, if the dataset is already the SCOP domain, this option shold be turned off, 2007-11-18, Nanjiang*/  
bool isIgnoreNonDiag = false;/*whether to ignore the non main diagnol matches, 2008-01-16, Nanjiang*/  

int beginID = 0;       /*in order to run the program simultaneously, setting the begin and end position of the idlist to run*/
int endID = 0x7FFFFFFF; /*by default, running all ids in the idListFile, 2007-11-16 */

#define SIZE_ID 30
/*****************************************************************************
 * NOTE: printing out the beg end alignment for the target sequence and the
 * candidate sequence, the position is from 0, that is
 * [begPos, endPos] = [0, 12], means aaSeqIndex >= 0 and aaSeqIndex < 12 is
 * with the segment, 2007-07-19
 ****************************************************************************/
int max_num_diag = 2; /*maxium number of diagnol to be used for each plot, which means a two domain polypeptide*/

double DotSegLengthScore[LONGEST_SHAPE] ; /*the raw score of dots with respect to the length of the continuous segment*/
double DotSegDistWeight[LONGEST_SHAPE] ; /* the raw weights for two segments seprated by N distance*/
double DiagDistWeight[LONGEST_SHAPE] ; /* the raw weights for diagnol seprated by N distance*/

struct DotSegment
{
    int beg;
    int end;
    int size;
};
struct Rectangle
{
    int x1;
    int y1;
    int x2;
    int y2;
};

void PrintHelp()
{
    fprintf(stdout,"Usage: postscan-frag-search [options] frag-search-file\n");
    fprintf(stdout,"options:\n");
    fprintf(stdout,"  -l filelist         : frag-search-file list\n");
    fprintf(stdout,"  --method 0|1|2      : method for frag search, 0 - tuping's old algorithm,  default = 2\n");
    fprintf(stdout,"                      : method 1 -- my new algorithm, but still use BuildPairMatrix\n");
    fprintf(stdout,"                      : method 2 -- new algorithm, not using BuildPairMatrix\n");
    fprintf(stdout,"  --len-file file     : list for sequence length for the candidate sequences\n");
    fprintf(stdout,"  --test-len-file file: list for sequence length for the test sequences\n");
    fprintf(stdout,"                      : if not set, used the same lengthFile as for the candidate sequences\n");
    fprintf(stdout,"  --outpath path      : output path for postscaned frag search file, default = ./\n");
    fprintf(stdout,"  --ext extension     : file extension for postscaned frag search file, default = fragpost\n");
    fprintf(stdout,"  --topn int          : output only topN hits, default = 20\n");
    fprintf(stdout,"  --cs double         : cutoff for the raw hit percentage , default = max(10.0, pivatValue)\n");
    fprintf(stdout,"  --cps double        : cutoff for the postcaned score of hits, default = 5.0\n");
    fprintf(stdout,"  --cfragscore int    : cutoff for the fragscore, only dots with fragscore>cutoff are used, default = 0\n");
    fprintf(stdout,"  --nm                : turn off the multi-domain functionality\n");
    fprintf(stdout,"  --ignore-non-diag   : for domains with similar length, ignore the non main diagnol diags\n");
    fprintf(stdout,"  --weight 0|1|2|3    : set weight scheme, 0 -- no weight\n");
    fprintf(stdout,"                      :                    1 -- weighted by the frag-frag score\n");
    fprintf(stdout,"                      :                    2 -- weighted by the core region, if so, extra information should be supplied\n");
    fprintf(stdout,"                      :                    3 -- both 1 and 2, \n");
    fprintf(stdout,"  --rb                : read the binary frag file\n");
    fprintf(stdout,"  --begin   <int>     : start number of ids in the idListFile to run database build\n");
    fprintf(stdout,"  --end     <int>     : end number of ids in the idListFile to run database build\n");
    fprintf(stdout,"  --changelog         : print the change log and exit\n");
    fprintf(stdout,"  -h|--help           : print this help message and exit\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"Created on 2007-07-01, updated 2010-01-07, Nanjiang Shu\n");
}
void PrintVerboseHelp() { }

double GetDotSegScoreLength(int n)/*{{{*/
    // get the raw score for each dot according to the length of the segment
{
    //if (n == 0)
    //{
        //return 0.0;
    //}
    //else if (n <= 2)
    //{
        //return double(n/2);
    //}
    //else
    //{
        //return pow(double(n), 0.6);
    //}
    return double (n);
}/*}}}*/
double GetDotSegDistWeight(int n)/*{{{*/
    // get the weight for two segments separeted by N residues
{
    return pow(double(n+1), 0.5);
}/*}}}*/
double GetDiagDistWeight(int n)/*{{{*/
    // get the raw weight for two diagnols separated by N distance
{
    return pow(double(n+1), 0.5);
}/*}}}*/
int AnalyzeFragScore(float *fragScore, int num, float &meanFragScore, float &medianFragScore, float &minFragScore, float &maxFragScore, float &cutoff_fragscore)/*{{{*/
/*****************************************************************************
 * Analyze the frag score and get the parameters for cutoff_fragscore and
 * normalizing the fragscore to weight
 * minFragScore and maxFragScore are not the real minFragScore or maxFragScore,
 * but by ruling out top 0.1% and bottom 0.1%
 ****************************************************************************/
{
    int i;
    float pmin = 0.001; /*since the distribution of fragscore is not normal distribution, but something like extrem distribution with a long tail, use a larger value for the tail, 2008-01-16, Nanjiang*/
    float pmax = 0.01;
    float pcutoff = 0.05;
    Array1D <float> tmpFragScore_1darray(num);
    float* tmpFragScore = tmpFragScore_1darray.array1D;
    for(i = 0; i< num; i ++) { tmpFragScore[i] = fragScore[i]; }
    QuickSort(tmpFragScore, 0, num-1);
    meanFragScore = Average(tmpFragScore, 0, num-1);
    medianFragScore = tmpFragScore[num/2];
    minFragScore = tmpFragScore[Integer(num*pmin)];
    maxFragScore = tmpFragScore[Integer(num*(1.0-pmax)-1)];
    cutoff_fragscore = min(tmpFragScore[Integer(num*pcutoff)], float(45000.0)); /*45000 is the average number*/
    return 0;
}
/*}}}*/
void NormalizeFragScore(float* fragScore, float mean, float min, float max, int num)/*{{{*/
/*scale the allFragScore to 0~2, so that it can be used as weight, 2008-01-16,Nanjiang*/
{
    int i;
    float halfwidth = (max-min)/2.0;
    for(i = 0; i < num; i ++)
    {
        //fragScore[i] = (fragScore[i] - mean)/halfwidth;
        //if(fragScore[i] < 0.0) { fragScore[i] = 0.0; }
        //if(fragScore[i] > 2.0) { fragScore[i] = 2.0; }
        if(fragScore[i] >= max) 
        {
            fragScore[i] = 1.0 + (fragScore[i] - max) / (halfwidth);
        }
        else
        {
            fragScore[i] = 1.0;
        }
    }
}
/*}}}*/
int RemoveFragByCutoffScore(int *allPosTar, int *allPosCan, float *allFragScore, char **allProIDCan, int &numAllPair, float cutoff_fragscore)/*{{{*/
{   
    int i;
    Array1D <int> tmpAllPosCan_1darray(numAllPair);
    Array1D <int> tmpAllPosTar_1darray(numAllPair);
    Array1D <float> tmpAllFragScore_1darray(numAllPair);
    Array2D <char> tmpAllProIDCan_2darray(numAllPair, SIZE_ID+1);
    int *tmpAllPosCan = tmpAllPosCan_1darray.array1D;
    int *tmpAllPosTar = tmpAllPosTar_1darray.array1D;
    float *tmpAllFragScore = tmpAllFragScore_1darray.array1D;
    char **tmpAllProIDCan = tmpAllProIDCan_2darray.array2D;
    for(i= 0; i < numAllPair; i ++)
    {
        tmpAllPosTar[i] = allPosTar[i];
        tmpAllPosCan[i] = allPosCan[i];
        tmpAllFragScore[i] = allFragScore[i];
        my_strcpy(tmpAllProIDCan[i],allProIDCan[i],SIZE_ID);
    }
    int cntPair = 0;
    for(i = 0; i< numAllPair; i ++)
    {
        if(tmpAllFragScore[i] >= cutoff_fragscore)
        {
            allPosTar[cntPair] = tmpAllPosTar[i];
            allPosCan[cntPair] = tmpAllPosCan[i];
            allFragScore[cntPair] = tmpAllFragScore[i];
            my_strcpy(allProIDCan[cntPair],tmpAllProIDCan[i],SIZE_ID);
            cntPair ++;
        }
    }
    numAllPair = cntPair;
    return cntPair;
}
/*}}}*/
//int GetMatchedPart(int **pair_matrix, int posDiag, int corner, int length_target, int length_candidate, int &begPosTar, int &endPosTar, int &begPosCan, int &endPosCan)[>{{{<]
/*****************************************************************************
 * given dot line with highest score in the pair_matrix, give out the region
 * where candidate sequence and the target sequence aligned
 * begPos start from 1
 ****************************************************************************/
//{
//    double score;
//    int i,j,k ;
//    int numDot = 0;
//    Array1D <int> dot_1darray (length_target+length_candidate);
//    int *dot = dot_1darray.array1D;
//    [>get the 1D line for the main diagnol<]
//    numDot = GetDiag(dotMain, pair_matrix, diagPos, length_target, length_candidate, corner);

//    Array1D <DotSegment> dotSegMain_1darray(numDotMain);
//    DotSegment *dotSegMain = dotSegMain_1darray.array1D;
//    for(i = 0 ; i < numDotMain; i ++) { InitDotSegment(&dotSegMain[i]); }
//    int numDotSegMain = GetDotSeg(dotMain, numDotMain, dotSegMain);



//}[>}}}<]
int GetDiagPos(int x, int y, int sizeX, int sizeY, int &posDiag, int &cornerDiag)/*{{{*/
/*****************************************************************************
 * give the dot on the 2D plot (x,y) and sizeX and sizeY, posDiag and
 * cornerDiag
 * return the size of the diagnol sizeDiag
 ****************************************************************************/
{
    posDiag = abs(x-y);
    cornerDiag = ((x > y) ? UPPER_CORNER : LOWER_CORNER);
    if(cornerDiag == UPPER_CORNER) 
    { 
        return (y + min(sizeX - x, sizeY-y)); 
    }
    else 
    { 
        return (x + min(sizeY - y, sizeX-x)); 
    }
}
/*}}}*/
int FindDiagIndex(int x,  int y, int *indexDiag, int *indexDiagCorner, int numDiag)/*{{{*/
/*****************************************************************************
 * give a dot (x,y) on the 2D plot and a set of diagnols (defined by the
 * indexDiag and corner, e.g, indexDiag[0] = 1, indexDiagCorner[0]=UPPER_CORNER
 * means the diagnol just above the main diagnol), return the index of the
 * diagnol in (0, numDiag-1) if found or -1 if not found
 ****************************************************************************/
{
    int i;
    int idxDiag = -1;
    int posDiag = abs(x-y);
    int cornerDiag = (x > y) ? UPPER_CORNER : LOWER_CORNER;
    for(i = 0 ; i < numDiag; i ++)
    {
        if(indexDiag[i] == posDiag && indexDiagCorner[i] == cornerDiag)
        {
            idxDiag = i;
            break;
        }
    }
    return idxDiag;
}
/*}}}*/
void InitDiag(DATATYPE_DOT_SCORE *diag, int sizeDiag)/*{{{*/
{
    for(int i = 0 ; i < sizeDiag; i ++)
    {
        diag[i] = EMPTY;
    }
}
/*}}}*/
void InitDotSegment(DotSegment *p)/*{{{*/
{
    p -> beg = 0;
    p -> end = 0;
    p -> size = 0;
}/*}}}*/
bool IsWithinRect(int x,int y, Rectangle *rect, int numRect)/*{{{*/
{
    for(int i = 0 ; i < numRect; i ++)
    {
        if(x >= rect[i].x1 && y >= rect[i].y1 && x < rect[i].x2 && y < rect[i].y2)
        {
            return true;
        }
    }
    return false;
}/*}}}*/
int GetDiag(DATATYPE_DOT_SCORE *dot, int **pair_matrix, int diagPos, int length_target, int length_candidate, int corner)/*{{{*/
{
    int i,j;
    int cntDot = 0;
    if(corner == LOWER_CORNER)
    {
        i = diagPos;
        j = 0;
    }
    else  /*UPPER_CORNER*/
    {
        i = 0;
        j = diagPos;
    }
    cntDot = 0;
    while(i < length_target && j < length_candidate)
    {
        dot[cntDot] = DATATYPE_DOT_SCORE(pair_matrix[i][j]);
        i++; j ++;
        cntDot ++;
    }
        
    return cntDot;
}/*}}}*/
void SetDiag(DATATYPE_DOT_SCORE *dot, int **pair_matrix, int diagPos, int length_target, int length_candidate, int corner)/*{{{*/
{
    int i,j;
    int cntDot = 0;
    if(corner == LOWER_CORNER)
    {
        i = diagPos;
        j = 0;
    }
    else  /*UPPER_CORNER*/
    {
        i = 0;
        j = diagPos;
    }
    cntDot = 0;
    while(i < length_target && j < length_candidate)
    {
        pair_matrix[i][j] = int (dot[cntDot]) ;
        i++; j ++;
        cntDot ++;
    }
}/*}}}*/
int GetDotSeg(DATATYPE_DOT_SCORE *dot, int numDot, DotSegment *dotSeg)/*{{{*/
/*****************************************************************************
 * Get the continuous segments of dots in the line
 * for a diagnol like this
 *   xx  xxxx xx  xxxxxxx  xxx
 * there will be 5 segments and the lenght of segments are 2, 4, 2, 7 and 3.
 ****************************************************************************/
{
    int cntDotSeg = 0;
    int i = 0 ;
    while(i < numDot)
    {
        if(dot[i] == EMPTY)
        {
            i ++;
        }
        else
        {
            int j = 1;
            if (i + j < numDot)      /*2010-01-07 */
            {
                while(dot[i+j] > 0)
                {
                    j ++;
                    if (i+j >= numDot)      /*2010-01-07 */
                    {
                        break;
                    }
                }
            }
            for(int k = 0 ; k < j ; k ++)
            {
                dot[i+k] = scale_full * DotSegLengthScore[j]; 
                /*set the raw score of each dots as the func(length) of the
                 * segment , DotSegLengthScore is a global array storing the 
                 * score with respect to the length of the segments 2008-04-24*/
            }
            dotSeg[cntDotSeg].beg = i;
            dotSeg[cntDotSeg].end = i+j;
            dotSeg[cntDotSeg].size = j ;
            cntDotSeg ++;
            i += j;
        }
    }
    return cntDotSeg;
}/*}}}*/
int SetDiagDotScore(DATATYPE_DOT_SCORE *diag, float *fragScore, int sizeDiag, int numFilledDotDiag, int weightScheme = 0)/*{{{*/
/*****************************************************************************
 * set the score of each dot on the diagnol
 ****************************************************************************/
{
    int i,j;
    Array1D <DotSegment> dotSeg_1darray(numFilledDotDiag);
    DotSegment *dotSeg = dotSeg_1darray.array1D;
    //for(i = 0 ; i < numFilledDotDiag; i ++) { InitDotSegment(&dotSeg[i]); }
    int numDotSeg = GetDotSeg(diag, sizeDiag, dotSeg);
    for(i = 0 ; i < numDotSeg; i ++)
    {
        double scoreSum = 0.0;
        for(j = 0 ; j < numDotSeg ; j ++)
        {
            double score = 0.0;
            int dist = 0;
            if(j == i) { continue; } /*ignore the i segment itself*/
            
            if(j < i)
            {
                dist = dotSeg[i].beg - dotSeg[j].end;
            }
            else if(j > i)
            {
                dist = dotSeg[j].beg - dotSeg[i].end;
            }
            score = DotSegLengthScore[dotSeg[j].size] / DotSegDistWeight[dist] * scale_half;  /*changed 2008-04-24, DiagDistScore is a global variable*/
            scoreSum += score;
        }
        double totalScore = scoreSum + DotSegLengthScore[dotSeg[i].size] * scale_full;
        for(j = dotSeg[i].beg; j < dotSeg[i].end; j ++)
        {
            if(weightScheme == 0)
            {
                diag[j] = Integer (totalScore);
            }
            else if (weightScheme ==1)
            {
                diag[j] = Integer (totalScore*fragScore[j]); /*multiplied by fragScore[j] which rangs from 0 to 2, 2008-01-16, Nanjiang*/
            }
        }
    }
    return numDotSeg;
}
/*}}}*/

void Map1D_2D(int diagPos, int pos1D, int &x, int &y, int corner)/*{{{*/
/*****************************************************************************
 * given pos1D, return the position of the dot on 2D plot (x,y)
 ****************************************************************************/
{
    if(corner == LOWER_CORNER) /*but fixed 2007-07-19, x, y mixed before*/
    {
        x = pos1D ;
        y = pos1D + diagPos;
    }
    else
    {
        x = pos1D + diagPos;
        y = pos1D ;
    }
}/*}}}*/
int Map2D_1D(int x, int y, int corner)/*{{{*/
/*****************************************************************************
 * return the position of the dot on 2D plot (x,y), return the position of the
 * dot on the diagnol,  pos1D
 ****************************************************************************/
{
    if(corner == LOWER_CORNER)
    {
        return x ; /*pos1D*/
    }
    else
    {
        return y ; /*pos1D*/ 
    }
}/*}}}*/

int GetNumSatisfiedDot(int **pair_matrix, int diagPos, int length_target, int length_candidate, int corner)/*{{{*/
{
    int i,j;
    int numDot = 0;
    Array1D <DATATYPE_DOT_SCORE> dot_1darray (length_target+length_candidate);
    DATATYPE_DOT_SCORE *dot = dot_1darray.array1D;

    /*get the diagnol line, store in 1D array dot*/
    numDot = GetDiag(dot, pair_matrix, diagPos, length_target, length_candidate, corner);
    /* now checking the 1D array dot[0..numDot-1]*/
    Array1D <DotSegment> dotSeg_1darray(numDot);
    DotSegment *dotSeg = dotSeg_1darray.array1D;
    for(i = 0 ; i < numDot; i ++) { InitDotSegment(&dotSeg[i]); }

    int numDotSeg = GetDotSeg(dot, numDot, dotSeg);

    for(i = 0 ; i < numDotSeg; i ++)
    {
        double scoreSum = 0.0;
        for(j = 0 ; j < numDotSeg ; j ++)
        {
            double score = 0.0;
            int dist = 0;
            if(j == i) { continue; }
            
            if(j < i)
            {
                dist = dotSeg[i].beg - dotSeg[j].end;
            }
            else if(j > i)
            {
                dist = dotSeg[j].beg - dotSeg[i].end;
            }
            score = dotSeg[j].size / pow(double(dist+1), 0.6) * scale_half;
            scoreSum += score;
        }
        double totalScore = scoreSum + dotSeg[i].size * scale_full;
        for(j = dotSeg[i].beg; j < dotSeg[i].end; j ++)
        {
            dot[j] = Integer (totalScore);
        }
    }

    SetDiag(dot, pair_matrix, diagPos, length_target, length_candidate, corner);

    int cntFilledDot = 0;
    int maxSegSize = MIN_INT;
    for(i = 0 ; i < numDotSeg; i ++)
    {
        cntFilledDot += dotSeg[i].size;
        if(dotSeg[i].size > maxSegSize)
        {
            maxSegSize = dotSeg[i].size;
        }
    }

    if(maxSegSize < 3 && (double (cntFilledDot) / numDot < 0.1))
    {
        return 0;
    }
    else
    {
        return cntFilledDot;
    }
}/*}}}*/
int GetRect(Rectangle *rect, int diagPos, DotSegment *dotSeg, int numDotSeg, int numDot, int rowSize, int colSize, int corner)/*{{{*/
{
    int i;
    int cntRect  = 0 ;
    if(dotSeg[0].beg > 0)
    {
        rect[cntRect].x1 = 0;
        rect[cntRect].y1 = 0;
        Map1D_2D(diagPos, dotSeg[0].beg, rect[cntRect].x2, rect[cntRect].y2, corner);
        cntRect ++;
    }

    for(i = 1 ; i < numDotSeg; i ++)
    {
        Map1D_2D(diagPos, dotSeg[i-1].end, rect[cntRect].x1, rect[cntRect].y1, corner);
        Map1D_2D(diagPos, dotSeg[i].beg, rect[cntRect].x2, rect[cntRect].y2, corner);
        cntRect ++;
    }

    if(dotSeg[numDotSeg-1].end < numDot)
    {
        Map1D_2D(diagPos, dotSeg[numDotSeg-1].end, rect[cntRect].x1, rect[cntRect].y1, corner);
        rect[cntRect].x2 = colSize;
        rect[cntRect].y2 = rowSize;
        cntRect ++;
    }

    return cntRect;
}/*}}}*/
    
double GetDiagScore2(int pivat, DATATYPE_DOT_SCORE **diag, int *indexDiag, int *indexDiagCorner, int *sizeDiag, int ** indexFilledDotDiag, int *numFilledDotDiag,int numDiag,  int sizeX, int sizeY, int &begX, int &endX, int &begY, int &endY)/*{{{*/
/*****************************************************************************
 * calculate the score with the diagonal line at the diagPos as the main line
 * return the score and the start and end position for both candidate sequence
 * and target sequence
 ****************************************************************************/
{
    int i,k ;

    Array1D <DotSegment> dotSegMain_1darray(numFilledDotDiag[pivat]);
    DotSegment *dotSegMain = dotSegMain_1darray.array1D;
    //for(i = 0 ; i < numDotMain; i ++) { InitDotSegment(&dotSegMain[i]); }
    int numDotSegMain = GetDotSeg(diag[pivat], sizeDiag[pivat], dotSegMain);

    /*first get the begin and end position for the aligned seuqneces according
     * to the main diagnol line*/
    int diagPosMain = indexDiag[pivat];
    int cornerMain = indexDiagCorner[pivat];
    int numFilledDotDiagMain = numFilledDotDiag[pivat];
    int numDotMain = sizeDiag[pivat];
    Map1D_2D(diagPosMain, dotSegMain[0].beg, begX, begY, cornerMain);
    Map1D_2D(diagPosMain, dotSegMain[numDotSegMain-1].end, endX, endY, cornerMain);

    
    Array1D <Rectangle> rectMain_1darray(numFilledDotDiagMain+4);
    Rectangle *rectMain = rectMain_1darray.array1D;

    int numRect = GetRect(rectMain, diagPosMain, dotSegMain, numDotSegMain, numDotMain, sizeY, sizeX, cornerMain);

    double scoreSum  = 0.0;
    if (numFilledDotDiagMain < 3)
    {
        scoreSum = 0.0;
    }
    else
    {
        /*first calculate the score of the main diagnol*/
        for (k = 0 ; k < numFilledDotDiagMain; k ++)
        {
            scoreSum += double (diag[pivat][indexFilledDotDiag[pivat][k]]);
        }
        int iDiag = 0;
        int x,y;
        /*add up the score at other diagnols*/
        for(iDiag = 0 ; iDiag < numDiag; iDiag ++)
        {
            int diagPosTmp = indexDiag[iDiag] ;
            int cornerTmp = indexDiagCorner[iDiag];
            if(cornerMain ==  cornerTmp && diagPosMain == diagPosTmp) { continue; }

            //if (numFilledDotDiag[iDiag] < 3){continue;}

            int dist = 0;
            double avgLen = double (sizeX+sizeY) /2.0;
            if(cornerMain == cornerTmp)
            { dist = abs(diagPosTmp - diagPosMain); }
            else 
            { dist = diagPosTmp + diagPosMain -1; }

            //i = diagPosTmp;    [>i is the iterator for row, sizeRow = sizeY<]
            //j = 0;             [>j is the iterator for column, sizeColumn = sizeX<] 
            double tempScale = DiagDistWeight[dist];
            for(k = 0; k < numFilledDotDiag[iDiag]; k ++)
            {
                DATATYPE_DOT_SCORE diagScore =  diag[iDiag][indexFilledDotDiag[iDiag][k]] ;
                if (diagScore != EMPTY)
                {
                    Map1D_2D(diagPosTmp, indexFilledDotDiag[iDiag][k], x, y, cornerTmp);
                    if(IsWithinRect(x,y, rectMain, numRect))
                    {
                        scoreSum += diagScore / tempScale;
                        if(diagScore >= 3 * scale_full && dist/avgLen < 0.5 && (diagScore/scale_full/((dist+1)/avgLen)) >= 30.0 )
                        {
                            if(x < begX) { begX = x; }
                            if(x > endX) { endX = x; }
                            if(y < begY) { begY = y; }
                            if(y > endY) { endY = y; }
                        }
                    }
                }
            }
                
            //while(i < sizeY && j < sizeX)
            //{
            //    if(pair_matrix[i][j] != EMPTY)
            //    {
            //        if(IsWithinRect(j,i, rectMain, numRect))
            //        {
            //            scoreSum += pair_matrix[i][j] / pow(double(dist+1), 0.5);
            //            if(pair_matrix[i][j] >= 3 * scale_full && dist/avgLen < 0.5 && (pair_matrix[i][j]/scale_full/((dist+1)/avgLen)) >= 30.0 )
            //            {
            //                if(j < begPosCan) { begPosCan = j; }
            //                if(j > endPosCan) { endPosCan = j; }
            //                if(i < begPosTar) { begPosTar = i; }
            //                if(i > endPosTar) { endPosTar = i; }
            //            }
            //        }
            //    }
            //    i++; j ++;
            //}
        }

        /*extend the aligned region along the diagnol by 5% of the sequence
         * length
         * function:= 1.6*x^0.5 
         * which result slope = 0.05 when sequence length = 250
         *                    = 0.14 when sequence length = 30*/

        int max_enxtend_length = Integer(1.6*pow(double(min(sizeX, sizeY)), 0.5));
        for(i = 0 ; i < max_enxtend_length; i ++)
        {
            if(begX <= 0 || begY <= 0)
            { break; }
            begX --;
            begY --;
        }
        for(i = 0 ; i < max_enxtend_length; i ++)
        {
            if(endY >= sizeY || begX >= sizeX)
            { break; }
            endY ++;
            endX ++;
        }

    }
    return scoreSum;
}
/*}}}*/

double GetDiagScore(int diagPos,int corner, int *indexDiagLowerCorner, int numDiagLowerCorner, int *indexDiagUpperCorner,int numDiagUpperCorner, int **pair_matrix, int length_target, int length_candidate, int &begPosTar, int &endPosTar, int &begPosCan, int &endPosCan) /*{{{*/
/*****************************************************************************
 * calculate the score with the diagonal line at the diagPos as the main line
 * return the score and the start and end position for both candidate sequence
 * and target sequence
 ****************************************************************************/
{
    int i,j,k ;
    int numDotMain = 0;
    Array1D <DATATYPE_DOT_SCORE> dotMain_1darray (length_target+length_candidate);
    DATATYPE_DOT_SCORE *dotMain = dotMain_1darray.array1D;
    /*get the 1D line for the main diagnol*/
    numDotMain = GetDiag(dotMain, pair_matrix, diagPos, length_target, length_candidate, corner);

    Array1D <DotSegment> dotSegMain_1darray(numDotMain);
    DotSegment *dotSegMain = dotSegMain_1darray.array1D;
    for(i = 0 ; i < numDotMain; i ++) { InitDotSegment(&dotSegMain[i]); }
    int numDotSegMain = GetDotSeg(dotMain, numDotMain, dotSegMain);


    /*first get the begin and end position for the aligned seuqneces according
     * to the main diagnol line*/
    Map1D_2D(diagPos, dotSegMain[0].beg, begPosCan, begPosTar, corner);
    Map1D_2D(diagPos, dotSegMain[numDotSegMain-1].end, endPosCan, endPosTar, corner);

    
    Array1D <Rectangle> rectMain_1darray(numDotSegMain+4);
    Rectangle *rectMain = rectMain_1darray.array1D;

    int numRect = GetRect(rectMain, diagPos, dotSegMain, numDotSegMain, numDotMain, length_target, length_candidate, corner);
    

    double scoreSum  = 0.0;
    /*first calculate the score of the main diagnol*/
    for (k = 0 ; k < numDotMain; k ++)
    {
        if(dotMain[k] != EMPTY)
        {
            scoreSum += double (dotMain[k]);
        }
    }

    int iDiag;
    /*first calculate the score of the diagnols in the lower corner*/
    for(iDiag = 0 ; iDiag < numDiagLowerCorner; iDiag ++)
    {
        if(corner == LOWER_CORNER && indexDiagLowerCorner[iDiag] == diagPos) { continue; }
        int diagPosTmp = indexDiagLowerCorner[iDiag] ;
        int dist = 0;
        double avgLen = double (length_candidate+length_target) /2.0;
        if(corner == LOWER_CORNER)
        { dist = abs(diagPosTmp - diagPos); }
        else 
        { dist = diagPosTmp + diagPos -1; }

        i = diagPosTmp;    /*i is the iterator for row, sizeRow = length_target*/
        j = 0;             /*j is the iterator for column, sizeColumn = length_candidate*/ 
        while(i < length_target && j < length_candidate)
        {
            if(pair_matrix[i][j] != EMPTY)
            {
                if(IsWithinRect(j,i, rectMain, numRect))
                {
                    scoreSum += pair_matrix[i][j] / pow(double(dist+1), 0.5);
                    if(pair_matrix[i][j] >= 3 * scale_full && dist/avgLen < 0.5 && (pair_matrix[i][j]/scale_full/((dist+1)/avgLen)) >= 30.0 )
                    {
                        if(j < begPosCan) { begPosCan = j; }
                        if(j > endPosCan) { endPosCan = j; }
                        if(i < begPosTar) { begPosTar = i; }
                        if(i > endPosTar) { endPosTar = i; }
                    }
                }
            }
            i++; j ++;
        }
    }

    for(iDiag = 0 ; iDiag < numDiagUpperCorner; iDiag ++)
    {
        if(corner == UPPER_CORNER && indexDiagUpperCorner[iDiag] == diagPos) { continue; }
        int diagPosTmp = indexDiagUpperCorner[iDiag] ;
        int dist = 0;
        if(corner == UPPER_CORNER)
        { dist = abs(diagPosTmp - diagPos); }
        else 
        { dist = diagPosTmp + diagPos -1; }

        double avgLen = double (length_candidate+length_target) /2.0;
        i = 0;
        j = diagPosTmp;
        while(i < length_target && j < length_candidate)
        {
            if(pair_matrix[i][j] != EMPTY)
            {
                if(IsWithinRect(j,i, rectMain, numRect) )
                {
                    scoreSum += pair_matrix[i][j]/ pow(double(dist+1), 0.5);
                    if(pair_matrix[i][j] >= 3 * scale_full && dist/avgLen < 0.5 && (pair_matrix[i][j]/scale_full/((dist+1)/avgLen)) >= 30.0 )
                    {
                        if(j < begPosCan) { begPosCan = j; }
                        if(j > endPosCan) { endPosCan = j; }
                        if(i < begPosTar) { begPosTar = i; }
                        if(i > endPosTar) { endPosTar = i; }
                    }
                }
            }
            i++; j ++;
        }
    }

    /*extend the aligned region along the diagnol by 5% of the sequence
     * length
     * function:= 1.6*x^0.5 
     * which result slope = 0.05 when sequence length = 250
     *                    = 0.14 when sequence length = 30*/

    int max_enxtend_length = Integer(1.6*pow(double(min(length_target, length_candidate)), 0.5));
    for(i = 0 ; i < max_enxtend_length; i ++)
    {
        if(begPosTar <= 0 || begPosCan <= 0)
        { break; }
        begPosTar --;
        begPosCan --;
    }
    for(i = 0 ; i < max_enxtend_length; i ++)
    {
        if(endPosTar >= length_target || begPosCan >= length_candidate)
        { break; }
        endPosTar ++;
        endPosCan ++;
    }

    return scoreSum;
}/*}}}*/
int ReadInSeqLenth(const char * seqlenFile, char **proIDList, int *proLength)/*{{{*/
/*****************************************************************************
 * read in the length list file for each chains, in the format of
 * ID1 length1
 * ID2 length2
 ****************************************************************************/
{
    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;
    FILE *fpin = fopen(seqlenFile, "r");
    checkfilestream(fpin, seqlenFile, "r", true);
    bool isSorted = true;//assume the seqlenFile is sorted ascendingly according to the proID
    int cntID = 0;

    char id[SIZE_ID+1]  ="";
    int seqlen = 0;
    char idFormer[SIZE_ID+1] = "";

    while((linesize = fgetline(fpin, line ,maxline )) != EOF)
    {
        if(linesize <= 0 ) continue;
        if (sscanf(line, "%s %d", id, &seqlen) == 2)
        {
            my_strcpy(proIDList[cntID], id, SIZE_ID);
            proLength[cntID] = seqlen;
            if(strcmp(idFormer, id) > 0) 
            {
                isSorted = false;
            }
            my_strcpy(idFormer, id, SIZE_ID);
            cntID ++;
        }
    }
    fclose(fpin);

    int numID = cntID;
    if(!isSorted) // if the seqlenFile is not sorted ascendingly according to the proID, sort it
    {

        Array2D <char> tmpProIDList_2darray(numID, SIZE_ID+1);
        char **tmpProIDList = tmpProIDList_2darray.array2D;
        Array1D <int> tmpProLength_1darray(numID);
        int *tmpProLength = tmpProLength_1darray.array1D;
        int i = 0 ;
        for(i = 0 ; i < numID; i ++)
        {
            my_strcpy(tmpProIDList[i], proIDList[i], SIZE_ID);
            tmpProLength[i] = proLength[i];
        }
        Array1D <int> idx_1darray(numID);
        int *idx = idx_1darray.array1D;
        for(i = 0 ; i < numID; i ++)
        {
            idx[i] = i;
        }

        QuickSort_String(idx, proIDList, 0 , numID-1);
        for(i = 0 ; i < numID; i ++)
        {
            my_strcpy(proIDList[i], tmpProIDList[idx[i]], SIZE_ID);
            proLength[i] = tmpProLength[idx[i]];
        }
    }

    return numID;
}
/*}}}*/
void PrintPairMatrix(int *pair_matrix, int length_target, int length_candidate)/*{{{*/
{
    int i,j;
    int apos = 0 ;
    char printchar = ' ';
    for(i = 0 ; i < length_target; i ++)
    {
        for(j = 0 ; j < length_candidate; j ++)
        {
            apos = i * length_candidate + j;
            if(pair_matrix[apos] == EMPTY)
            {
                printchar = ' ';
            }
            else
            {
                printchar = 'X';
            }
            fprintf(stdout,"%c", printchar);
        }
        fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
}
/*}}}*/
void PrintPairMatrix(int **pair_matrix, int length_target, int length_candidate)/*{{{*/
{
    int i,j;
    char printchar = ' ';
    for(i = 0 ; i < length_target; i ++)
    {
        for(j = 0 ; j < length_candidate; j ++)
        {
            if(pair_matrix[i][j] == EMPTY)
            {
                printchar = ' ';
            }
            else
            {
                printchar = 'X';
            }
            fprintf(stdout,"%c", printchar);
        }
        fprintf(stdout,"\n");
    }
    fprintf(stdout,"\n");
}
/*}}}*/
int CheckFragFileFormat(const char *fragfile)/*{{{*/
{
    int linesize;
    int maxline = 500;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;
    FILE *fpin = fopen(fragfile, "r");
    int fragformat = 0;
    if (checkfilestream(fpin,fragfile, "r", true) != -1)
    {
        linesize = fgetline(fpin, line, maxline);
        char *tmpline = strtrim(line);
        char *pch; 
        pch = strtok (tmpline, WHITE_SPACE);
        int cntField=0;
        while (pch != NULL)
        {
            cntField ++;
            pch = strtok(NULL, WHITE_SPACE);
        }
#ifdef DEBUG
        fprintf(stdout,"fragfile=%s, line=%s, cntField=%d\n",fragfile, line, cntField );
#endif
        if (cntField == 5)
        {
            fragformat = FRAGFORMAT_NANJIANG;
        }
        else if (cntField == 6)
        {
            fragformat = FRAGFORMAT_TUPING;
        }
        else
        {
            fprintf(stderr,"Warning! Unrecognized format for the fragfile: %s, line = %s\n", fragfile, line);
        }
    }
    fclose(fpin);
    return fragformat;
}/*}}}*/
int ReadInFragParameter(char *file, int &fragformat, int &numID, int &maxSizeID, int &length, int &totalFragCan)/*{{{*/
/* read in the fragformat
 * totalFragCan: the total number of candidate fragments for all target positions
 * */
{
    FILE *fpin = fopen(file,"rb");
    if (  fpin == NULL  )
	{
		return -1;
	}
	size_t nread = 0;/*return value for fread*/
    /*read the fragFileType*/
    nread=fread(&fragformat, sizeof(int), 1, fpin); /*read the fragFileType*/ 
	if (nread != 1)
	{ fprintf(stderr,"fread error! File:%s, function:%s, var:%s\n",file,"ReadInFragParameter","fragformat"); }
    nread=fread(&numID, sizeof(int), 1, fpin); /*read the number of unique ids*/
	if (nread != 1)
	{ fprintf(stderr,"fread error! File:%s, function:%s, var:%s\n",file,"ReadInFragParameter","numID"); }
    nread=fread(&maxSizeID, sizeof(int),1,fpin);/*read the maximal size of ids*/
	if (nread != 1)
	{ fprintf(stderr,"fread error! File:%s, function:%s, var:%s\n",file,"ReadInFragParameter","maxSizeID"); }
    nread=fread(&length, sizeof(int), 1, fpin);  /*read length of the target sequence*/ 
	if (nread != 1)
	{ fprintf(stderr,"fread error! File:%s, function:%s, var:%s\n",file,"ReadInFragParameter","length"); }
    nread=fread(&totalFragCan, sizeof(int), 1, fpin);   /*read the total number of candidate fragments*/ 
	if (nread != 1)
	{ fprintf(stderr,"fread error! File:%s, function:%s, var:%s\n",file,"ReadInFragParameter","totalFragCan"); }
    fclose(fpin);
    return fragformat;
}/*}}}*/
int ReadInPair(const char *fragPairFile, int *posTar = NULL, int *posCan = NULL, float *fragScore = NULL, char **proIDCan = NULL)/*{{{*/
/*****************************************************************************
 * ChangeLog 2008-01-25
 * the format of frag file has been change to five columns
 *      seqTarget rank idCan seqCan score
 ****************************************************************************/
{
	fragformat = CheckFragFileFormat(fragPairFile);

    FILE *fpin = fopen(fragPairFile, "r");
    checkfilestream(fpin, fragPairFile, "r");
    int status_sscanf = 0;

    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;
    int pos1, pos2;
    float score;
    int rank;
    char proid[100] = "";
	int idxChain = 0; /*this is only for the frag file of Tuping's format, the index of the chain in the chain id list 2009-06-25 */
    int cntPair = 0;
    while((linesize = fgetline(fpin, line, maxline)) !=EOF)
    {
        if(linesize  <= 0) continue;
		if (fragformat == FRAGFORMAT_NANJIANG)
		{
			status_sscanf = sscanf(line,"%d %d %s %d %f", &pos1, &rank, proid, &pos2, &score);
			if (status_sscanf != 5)
			{
				fprintf(stderr, "Error fragformat! fragfile=%s, line=%s, fragformat=%d\n", fragPairFile, line, fragformat);
			}
			else
			{
				/*read in the frag-frag score as well, Nanjiang, 2008-01-15 */   /*cutoff_fragscore added 2008-01-16, Nanjiang*/ 
				if( posTar != NULL) posTar[cntPair] = pos1;
				if (posCan != NULL) posCan[cntPair] = pos2;
				if (fragScore != NULL) fragScore[cntPair] = score;
				if (proIDCan != NULL) { my_strcpy(proIDCan[cntPair], proid, SIZE_ID);}
				cntPair ++;
			}
		}
		else 
		{
			status_sscanf = sscanf(line,"%d %d %s %d %d %f", &pos1, &rank, proid, &idxChain, &pos2, &score);
			if (status_sscanf != 6)
			{
				fprintf(stderr, "Error fragformat! fragfile=%s, line=%s, fragformat=%d\n", fragPairFile, line, fragformat);
			}
			else
			{
				if( posTar != NULL) posTar[cntPair] = pos1;
				if (posCan != NULL) posCan[cntPair] = pos2;
				if (fragScore != NULL) fragScore[cntPair] = score;
				if (proIDCan != NULL) { my_strcpy(proIDCan[cntPair], proid, SIZE_ID);}
				cntPair ++;
			}
		}
    }
    fclose(fpin);
    return cntPair;
}/*}}}*/
int ReadInPairBinary(const char *binfragPairFile, int *posTar = NULL, int *posCan = NULL, float *fragScore = NULL, char **proIDCan = NULL)/*{{{*/
/*****************************************************************************
 * Read in pairs from binary frag files
 * 2009-06-25
 ****************************************************************************/
{
    /*check whether the input file is a vaid fragbinary file*/
    int lenFileName = strlen(binfragPairFile);
    if(strcmp (&(binfragPairFile[lenFileName-3]), "bin" ) != 0)
    {
        fprintf(stderr,"Error, the binfragfile \"%s\" is not end with \"bin\"\n", binfragPairFile);
        return -1;
    }

    int numID = 0;
    int maxSizeID = 0;
    int length = 0;
    int totalFragCan = 0;

    ReadInFragParameter((char*)binfragPairFile, fragformat, numID, maxSizeID, length, totalFragCan);
    Array2D <char> idList_2darray(numID, maxSizeID+1);
    idList_2darray.Init('\0');
    char **idList = idList_2darray.array2D;

    Array1D <short> tmpposTar_1darray(length);
    Array1D <short> tmpnumCan_1darray(length);
    tmpposTar_1darray.Init(0);
    tmpnumCan_1darray.Init(0);
    short *tmpposTar = tmpposTar_1darray.array1D;
    short *tmpnumCan = tmpnumCan_1darray.array1D;
    FragCanShort5 *tmpfragCan5 = NULL;
    FragCanShort6 *tmpfragCan6 = NULL;
    if (fragformat == FRAGFORMAT_NANJIANG)
    { tmpfragCan5 = new FragCanShort5[totalFragCan+1]; }
    else 
    { tmpfragCan6 = new FragCanShort6[totalFragCan+1]; }
    
    if (fragformat == FRAGFORMAT_NANJIANG)
    {
        GetBinaryFragShort5(binfragPairFile, fragformat, idList, numID, maxSizeID, length, tmpposTar, tmpnumCan, totalFragCan, tmpfragCan5);
    }
    else
    {
        GetBinaryFragShort6(binfragPairFile, fragformat, idList, numID, maxSizeID, length, tmpposTar, tmpnumCan, totalFragCan, tmpfragCan6);
    }

    int cntPair = 0;
	int i,j;
    for (i = 0; i < length; i ++)
	{
		for (j = 0; j < tmpnumCan[i]; j ++)
		{
			if (fragformat == FRAGFORMAT_NANJIANG)
			{
				if( posTar != NULL) posTar[cntPair] = tmpposTar[i];
				if (posCan != NULL) posCan[cntPair] = tmpfragCan5[cntPair].posCan;
				if (fragScore != NULL) fragScore[cntPair] = tmpfragCan5[cntPair].score;
				if (proIDCan != NULL) { my_strcpy(proIDCan[cntPair], idList[tmpfragCan5[cntPair].idxInner], SIZE_ID);}
			}
			else
			{
				if( posTar != NULL) posTar[cntPair] = tmpposTar[i];
				if (posCan != NULL) posCan[cntPair] = tmpfragCan6[cntPair].posCan;
				if (fragScore != NULL) fragScore[cntPair] = tmpfragCan6[cntPair].score;
				if (proIDCan != NULL) { my_strcpy(proIDCan[cntPair], idList[tmpfragCan6[cntPair].idxInner], SIZE_ID);}
			}
			cntPair ++;
		}
	}

    /*free memory*/
    if (fragformat == FRAGFORMAT_NANJIANG)
    { delete [] tmpfragCan5; }
    else 
    { delete [] tmpfragCan6; }

    return cntPair;
}/*}}}*/

void  BuildPairMatrix(int *pair_matrix, int *posTar, int *posCan, int numPair, int length_target, int length_candidate)/*{{{*/
{
    int i,j;
    for(i = 0 ; i < length_target * length_candidate; i ++)
    {
        pair_matrix[i] = EMPTY;
    }

    for(j = 0 ; j < numPair; j ++) // matrix is 
        //    ----------------- posCan, for candidate
        //    |
        //    |
        //    |
        //    posTar, for target
    {
        i = posTar[j] * length_candidate + posCan[j];
        pair_matrix[i] = FILLED;
    }
}
/*}}}*/
void  BuildPairMatrix(int **pair_matrix, int *posTar, int *posCan, int numPair, int length_target, int length_candidate)/*{{{*/
{
    int i,j;
    for(i = 0 ; i < length_target ; i ++)
    {
        for(j = 0 ; j < length_candidate; j ++)
        {
            pair_matrix[i][j] = EMPTY;
        }
    }

    for(j = 0 ; j < numPair; j ++) // matrix is 
        //    ----------------- posCan, for candidate
        //    |
        //    |
        //    |
        //    posTar, for target
    {
        pair_matrix[posTar[j]][posCan[j]]= FILLED;
    }
}
/*}}}*/
int Treat_pair_Percent_homology(int *Datapair, int Lengthtar, int Lengthcan)/*{{{*/
{
    int Numbre_Pair, ik,ij, ip, Len2_beg, Len1_end, Len_cal_tar,Len_cal_can,Npos,nbeg;
    int Len_Segment, Npair_last_can,Npair_last_tar;
    BYTE byte0, byte1 , *Point_can,*Point_tar, Number_seg_less500;
    int Xtar_seg[500], Ycan_seg[500], XYLen_seg[500], Number_seg ;//suppose the number of consecutive segments < 500
//    int NSquare, Low_end_tar, Low_end_can, High_end_tar, High_end_can, NBoundary;
    int Max_Numbre_Pair, Keep_Segment, Icode, Num_Data, Max_Point_Line,Point_Line , Max_line_threshhold;
    int Num_blank, BBLANK, Num_Check;
    //CString xs;



    //Datapair[]--; -1 for no point, 1000 for point, 0-500 for the number of consecutive segments
    byte0 = BYTE (0);
    byte1 = BYTE (1);
    Numbre_Pair = 2000;
    Keep_Segment = 3;
    BBLANK = 3;
    Max_Point_Line = 0;
    Point_Line = 0;
    Max_line_threshhold = 5;//percent
    Num_Data = Lengthtar*Lengthcan;
    //keep only those pairs that are at least 3 in consecutive
    Number_seg_less500 = TRUE;
    while (   Number_seg_less500  )
    {
        Number_seg_less500 = FALSE;//suppose one time
        Len2_beg = Lengthcan - Keep_Segment;//first position--Lengthcan - 1, at leats 3 segment
        Number_seg = 0;
        for (ik=Len2_beg; ik>=0; ik--) //from top-left to origin
        {
            Len_cal_tar = 0;
            Len_Segment = 0;
            Point_Line = 0;
            Num_blank = 0;
            Num_Check = 0;
            for (ij=ik; ij<Lengthcan; ij++)
            {
                Npos = Len_cal_tar*Lengthcan + ij;
                if (  Datapair[Npos] >= 1000  )
                {
                    Len_Segment++;
                    Point_Line++;
                    Num_blank = 0;
                }
                else
                {
                    Num_blank++;
                    if  (   (Num_blank>=BBLANK) || (Len_cal_tar==(Lengthtar-1))   )
                    {
                        if (  Len_Segment < Keep_Segment  )//thrown away the points which are spread-- Len_Segment points
                        {
                            for (ip=0; ip<Num_Check; ip++)
                            {
                                Npos = (Len_cal_tar-ip)*Lengthcan + ij-ip;
                                if (  (Npos>=0) && (Npos<Num_Data)  )
                                {
                                    Datapair[Npos] = -1;
                                }
                                else
                                {
                                    //xs.Format("%4d %4d 1-1",Lengthtar,Lengthcan);
                                    ////AfxMessageBox(xs);
                                }
                            }
                        }
                        else if (  Len_Segment >= Keep_Segment  )
                        {
                            for (ip=0; ip<Num_Check; ip++)
                            {
                                Npos = (Len_cal_tar-ip-1)*Lengthcan + ij-ip-1;
                                if (  (Npos>=0) && (Npos<Num_Data)  )
                                {
                                    if (  Datapair[Npos] >= 0  )
                                    {
                                        Datapair[Npos] = Number_seg;
                                    }
                                }
                                else
                                {
                                    //xs.Format("%4d %4d 1-2",Lengthtar,Lengthcan);
                                    ////AfxMessageBox(xs);
                                }
                            }
                            Xtar_seg[Number_seg] = Len_cal_tar - Len_Segment;
                            Ycan_seg[Number_seg] = ij - Len_Segment;
                            XYLen_seg[Number_seg] = Len_Segment;
                            Number_seg++;
                            if (  Number_seg >= 500  )
                            {
                                ////AfxMessageBox("Number_seg >= 500 in CCommonClass::Treat_pair_Percent_homology(BYTE *Datapair, int Lengthtar, int Lengthcan)");
                                //return Numbre_Pair ;
                                ik = -1;
                                Keep_Segment = Keep_Segment + 3;
                                Number_seg_less500 = TRUE;
                                Len_Segment = 0;
                                for (ip=0; ip<Num_Data; ip++)
                                {
                                    if (  Datapair[ip] >= 0   )
                                    {
                                        Datapair[ip] = 1000;
                                    }
                                }
                                break;//ij
                            }
                        }
                        Num_Check = 0;
                        Num_blank = 0;
                        Len_Segment = 0;
                    }
                }
                Len_cal_tar++; 
                Num_Check++;
                if (  Len_cal_tar >= Lengthtar  )
                {
                    break;
                }
            }
            if (  Point_Line > Max_Point_Line )
            {
                Max_Point_Line = Point_Line ;
            }
        }



        Len1_end = Lengthcan - 2;
        for (ik=Len1_end; ik>=0; ik--)  //from origin to right domn
        {
            Len_cal_can = 0;
            Len_Segment = 0;
            Point_Line = 0;
            Num_blank = 0;
            Num_Check = 0;
            nbeg = Lengthcan - 1 - ik;
            for (ij=nbeg; ij<Lengthtar; ij++)
            {
                Npos = ij*Lengthcan + Len_cal_can;
                if (  Datapair[Npos] >=1000  )
                {
                    Len_Segment++;
                    Point_Line++;
                    Num_blank = 0;
                }
                else
                {
                    Num_blank++;
                    if (   (Num_blank>=BBLANK) || (Len_cal_can==(Lengthcan-1))   )
                    {
                        if (  Len_Segment < Keep_Segment  )//thrown away the points which are spread-- Len_Segment points
                        {
                            for (ip=0; ip<Num_Check; ip++)
                            {
                                Npos = (ij-ip)*Lengthcan + Len_cal_can-ip;
                                if (  (Npos>=0) && (Npos<Num_Data)  )
                                {
                                    Datapair[Npos] = -1;
                                }
                                else
                                {
                                    //xs.Format("%4d %4d  2-1",Lengthtar,Lengthcan);
                                    //AfxMessageBox(xs);
                                }
                            }
                        }
                        else if (  Len_Segment >= Keep_Segment  )
                        {
                            for (ip=0; ip<Num_Check; ip++)
                            {
                                Npos = (ij-ip-1)*Lengthcan + Len_cal_can-ip-1;
                                if (  (Npos>=0) && (Npos<Num_Data)  )
                                {
                                    if (  Datapair[Npos] >= 0  )
                                    {
                                        Datapair[Npos] = Number_seg;
                                    }
                                }
                                else
                                {
                                    //xs.Format("%4d %4d 2-2",Lengthtar,Lengthcan);
                                    //AfxMessageBox(xs);
                                }
                            }
                            Xtar_seg[Number_seg] = ij - Len_Segment;
                            Ycan_seg[Number_seg] = Lengthcan - Len_Segment;
                            XYLen_seg[Number_seg] = Len_Segment;
                            Number_seg++;
                            if (  Number_seg >= 500  )
                            {
                                ////AfxMessageBox("Number_seg >= 500 in CCommonClass::Treat_pair_Percent_homology(BYTE *Datapair, int Lengthtar, int Lengthcan)");
                                //return Numbre_Pair;
                                ik = -1;
                                Keep_Segment = Keep_Segment + 3;
                                Number_seg_less500 = TRUE;
                                Len_Segment = 0;
                                for (ip=0; ip<Num_Data; ip++)
                                {
                                    if (  Datapair[ip] >= 0   )
                                    {
                                        Datapair[ip] = 1000;
                                    }
                                }
                                break;//ij
                            }
                        }
                    }
                    Len_Segment = 0;
                    Num_blank = 0;
                    Num_Check = 0;
                }
                Len_cal_can++;
                Num_Check++;
                if (  Len_cal_can >= Lengthcan  )
                {
                    break;
                }
            }
            if (  Point_Line > Max_Point_Line )
            {
                Max_Point_Line = Point_Line ;
            }
        }
    }




    Point_can = new BYTE [Lengthcan+1];
    Point_tar = new BYTE [Lengthtar+1];
    Max_Numbre_Pair = 0;
    Number_seg = 0;


    ////treatment of consecutive segments--to seach the maxium of the matches, begin with a segment, then merge all others[>{{{<]
    //for (ik=0; ik<Number_seg; ik++)//begin--seed
    //{
    //    for (ij=0; ij<Lengthcan; ij++)
    //    {
    //        Point_can[ij] = byte0;
    //    }
    //    for (ij=0; ij<Lengthtar; ij++)
    //    {
    //        Point_tar[ij] = byte0;
    //    }
    //    //keep the seed
    //    for (ij=Ycan_seg[ik]; ij<Ycan_seg[ik]+XYLen_seg[ik]; ij++)
    //    {
    //        Point_can[ij] = byte1;
    //    }
    //    for (ij=Xtar_seg[ik]; ij<Xtar_seg[ik]+XYLen_seg[ik]; ij++)
    //    {
    //        Point_tar[ij] = byte1;
    //    }
    //    Low_end_can = Ycan_seg[ik];
    //    Low_end_tar = Xtar_seg[ik];
    //    High_end_can = Ycan_seg[ik] + XYLen_seg[ik] -1;
    //    High_end_tar = Xtar_seg[ik] + XYLen_seg[ik] -1;

    //    //merge all other segments
    //    //check the low-end and search for the closest segment to j-th segment on the left-down------low-end
    //    NSquare = 1;
    //    NBoundary = __min(Low_end_can, Low_end_tar);
    //    while ( NSquare <= NBoundary )
    //    {
    //        Number_seg_tmp = -1;
    //        Npos = (Low_end_tar-NSquare)*Lengthcan + Low_end_can-NSquare;//diagonal
    //        if  (  (Datapair[Npos]>=0) && (Datapair[Npos]<500)  )  
    //        {
    //            Number_seg_tmp = Datapair[Npos];
    //            tar_NSquare = Low_end_tar - NSquare;
    //            can_NSquare = Low_end_can - NSquare;
    //            Icode = 0;
    //        }
    //        if (   (Number_seg_tmp==-1) && (NSquare>=2)   )//parallel-up
    //        {
    //            Npos = (Low_end_tar-NSquare)*Lengthcan + Low_end_can-1;//parallel-up
    //            if  (  (Datapair[Npos]>=0) && (Datapair[Npos]<500)  )  
    //            {
    //                Number_seg_tmp = Datapair[Npos];
    //                tar_NSquare = Low_end_tar - NSquare;
    //                can_NSquare = Low_end_can - 1;
    //                Icode = 1;
    //            }
    //        }
    //        if   (   (Number_seg_tmp==-1) && (NSquare>=2)   )//vertical-right
    //        {
    //            Npos = (Low_end_tar-1)*Lengthcan + Low_end_can-NSquare;//vertical-right
    //            if  (  (Datapair[Npos]>=0) && (Datapair[Npos]<500)  )  
    //            {
    //                Number_seg_tmp = Datapair[Npos];
    //                tar_NSquare = Low_end_tar-1;
    //                can_NSquare = Low_end_can-NSquare;
    //                Icode = 2;
    //            }
    //        }
    //        if   (   (Number_seg_tmp==-1) && (NSquare>=2)   )//vertical-left
    //        {
    //            for (ip=0; ip<NSquare-2; ip++)
    //            {
    //                Npos = (Low_end_tar-NSquare)*Lengthcan + Low_end_can-2-ip;//vertical-left
    //                if  (  (Datapair[Npos]>=0) && (Datapair[Npos]<500)  ) 
    //                {
    //                    Number_seg_tmp = Datapair[Npos];
    //                    tar_NSquare = Low_end_tar-NSquare;
    //                    can_NSquare = Low_end_can-2-ip;
    //                    Icode = 3;
    //                    break;
    //                }
    //            }
    //        }

    //        if   (   (Number_seg_tmp==-1) && (NSquare>=2)   )//parallel-down
    //        {
    //            for (ip=0; ip<NSquare-2; ip++)
    //            {
    //                Npos = (Low_end_tar-2-ip)*Lengthcan + Low_end_can-NSquare;//parallel-down
    //                if  (  (Datapair[Npos]>=0) && (Datapair[Npos]<500)  ) 
    //                {
    //                    Number_seg_tmp = Datapair[Npos];
    //                    tar_NSquare = Low_end_tar-2-ip;
    //                    can_NSquare = Low_end_can-NSquare;
    //                    Icode = 4;
    //                    break;
    //                }
    //            }
    //        }
    //        if (  Number_seg_tmp >= 0  )
    //        {
    //            Low_end_can = Ycan_seg[Number_seg_tmp];
    //            Low_end_tar = Xtar_seg[Number_seg_tmp];
    //            for (ip=Ycan_seg[Number_seg_tmp]; ip<Ycan_seg[Number_seg_tmp]+XYLen_seg[Number_seg_tmp]; ip++)
    //            {
    //                if (  ip > can_NSquare  )
    //                {
    //                    break;
    //                }
    //                Point_can[ip] = byte1;
    //            }
    //            for (ip=Xtar_seg[Number_seg_tmp]; ip<Xtar_seg[Number_seg_tmp]+XYLen_seg[Number_seg_tmp]; ip++)
    //            {
    //                if (  ip > tar_NSquare  )
    //                {
    //                    break;
    //                }
    //                Point_tar[ip] = byte1;
    //            }
    //            Number_seg_tmp = -1;
    //            NSquare = 0;
    //            NBoundary = __min(Low_end_can, Low_end_tar);
    //        }
    //        NSquare++;
    //    }

    //    //check the high-end and search for the closest segment to j-th segment on the left-down------high-end
    //    NSquare = 1;
    //    NBoundary = __min(Lengthcan-High_end_can-1, Lengthtar-High_end_tar-1);
    //    while ( NSquare <= NBoundary )
    //    {
    //        Number_seg_tmp = -1;
    //        Npos = (High_end_tar+NSquare)*Lengthcan + High_end_can+NSquare;//diagonal
    //        if  (  (Datapair[Npos]>=0) && (Datapair[Npos]<500)  )  
    //        {
    //            Number_seg_tmp = Datapair[Npos];
    //            tar_NSquare = High_end_tar + NSquare;
    //            can_NSquare = High_end_can + NSquare;
    //        }
    //        if (   (Number_seg_tmp==-1) && (NSquare>=2)   )//parallel-down
    //        {
    //            Npos = (High_end_tar+NSquare)*Lengthcan + High_end_can+1;//parallel-down
    //            if  (  (Datapair[Npos]>=0) && (Datapair[Npos]<500)  )  
    //            {
    //                Number_seg_tmp = Datapair[Npos];
    //                tar_NSquare = High_end_tar + NSquare;
    //                can_NSquare = High_end_can + 1;
    //            }
    //        }
    //        if   (   (Number_seg_tmp==-1) && (NSquare>=2)   )//vertical-left
    //        {
    //            Npos = (High_end_tar+1)*Lengthcan + High_end_can+NSquare;//vertical-left
    //            if  (  (Datapair[Npos]>=0) && (Datapair[Npos]<500)  )  
    //            {
    //                Number_seg_tmp = Datapair[Npos];
    //                tar_NSquare = High_end_tar + 1;
    //                can_NSquare = High_end_can + NSquare;
    //            }
    //        }
    //        if   (   (Number_seg_tmp==-1) && (NSquare>=2)   )//vertical-right
    //        {
    //            for (ip=0; ip<NSquare-2; ip++)
    //            {
    //                Npos = (High_end_tar+NSquare)*Lengthcan + High_end_can+2+ip;//vertical-right
    //                if  (  (Datapair[Npos]>=0) && (Datapair[Npos]<500)  ) 
    //                {
    //                    Number_seg_tmp = Datapair[Npos];
    //                    tar_NSquare = High_end_tar + NSquare;
    //                    can_NSquare = High_end_can + 2 + ip;
    //                    break;
    //                }
    //            }
    //        }

    //        if   (   (Number_seg_tmp==-1) && (NSquare>=2)   )//parallel-up
    //        {
    //            for (ip=0; ip<NSquare-2; ip++)
    //            {
    //                Npos = (High_end_tar+2+ip)*Lengthcan + High_end_can+NSquare;//parallel-up
    //                if  (  (Datapair[Npos]>=0) && (Datapair[Npos]<500)  ) 
    //                {
    //                    Number_seg_tmp = Datapair[Npos];
    //                    tar_NSquare = High_end_tar + 2 + ip;
    //                    can_NSquare = High_end_can + NSquare;
    //                    break;
    //                }
    //            }
    //        }
    //        if (  Number_seg_tmp >= 0  )
    //        {
    //            High_end_can = Ycan_seg[Number_seg_tmp] + XYLen_seg[Number_seg_tmp] -1;
    //            High_end_tar = Xtar_seg[Number_seg_tmp] + XYLen_seg[Number_seg_tmp] -1;
    //            for (ip=Ycan_seg[Number_seg_tmp]; ip<Ycan_seg[Number_seg_tmp]+XYLen_seg[Number_seg_tmp]; ip++)
    //            {
    //                if (  ip < can_NSquare  )
    //                {
    //                    continue;
    //                }
    //                Point_can[ip] = byte1;
    //            }
    //            for (ip=Xtar_seg[Number_seg_tmp]; ip<Xtar_seg[Number_seg_tmp]+XYLen_seg[Number_seg_tmp]; ip++)
    //            {
    //                if (  ip < tar_NSquare  )
    //                {
    //                    continue;
    //                }
    //                Point_tar[ip] = byte1;
    //            }
    //            Number_seg_tmp = -1;
    //            NSquare = 0;
    //            NBoundary = __min(Lengthcan-High_end_can-1, Lengthtar-High_end_tar-1);
    //        }
    //        NSquare++;
    //    }

    //    Npair_last_can = 0;
    //    for (ij=0; ij<Lengthcan; ij++)
    //    {
    //        if (  Point_can[ij] == byte1  )
    //        {
    //            Npair_last_can++;
    //        }
    //    }
    //    Npair_last_tar = 0;
    //    for (ij=0; ij<Lengthtar; ij++)
    //    {
    //        if (  Point_tar[ij] == byte1  )
    //        {
    //            Npair_last_tar++;
    //        }
    //    }
    //    Numbre_Pair = __min(Npair_last_can,Npair_last_tar);
    //    if (  Numbre_Pair > Max_Numbre_Pair  )
    //    {
    //        Max_Numbre_Pair = Numbre_Pair;
    //    }
    //}[>}}}<]


    for (ij=0; ij<Lengthcan; ij++)
    {
        Point_can[ij] = byte0;
    }
    for (ij=0; ij<Lengthtar; ij++)
    {
        Point_tar[ij] = byte0;
    }


    for (ij=0; ij<Lengthtar; ij++)
    {
        for (ip=0; ip<Lengthcan; ip++)
        {
            Npos = ij*Lengthcan + ip;
            if  (  Datapair[Npos] >= 0  ) 
            {
                Point_can[ip] = byte1;
                Point_tar[ij] = byte1;
            }
        }
    }

    Npair_last_can = 0;
    for (ij=0; ij<Lengthcan; ij++)
    {
        if (  Point_can[ij] == byte1  )
        {
            Npair_last_can++;
        }
    }

    Npair_last_tar = 0;
    for (ij=0; ij<Lengthtar; ij++)
    {
        if (  Point_tar[ij] == byte1  )
        {
            Npair_last_tar++;
        }
    }

    Numbre_Pair = min(Npair_last_can,Npair_last_tar);
    Max_Numbre_Pair = Numbre_Pair;
    Len_Segment = min(Lengthcan,Lengthtar);//shorter sequence
    Icode = Max_Point_Line*100/Len_Segment;
    if (  Icode > Max_line_threshhold  )
    {
        //Max_Numbre_Pair = Numbre_Pair + (Icode-Max_line_threshhold)*Len_Segment;

    }

    delete [] Point_can;
    delete [] Point_tar;

    return Max_Numbre_Pair;
}
/*}}}*/
int  ScanPairDotPlot(int **pair_matrix, int length_target, int length_candidate, int *posMaxDiag, int *corner, double *scoreDiag, int *begPosTar, int *endPosTar, int *begPosCan, int *endPosCan, int &numMaxDiag)/*{{{*/
/*****************************************************************************
 * given the 2D pair_matrix matrix, find the best diagnol line represeting the
 * homologous relationship between the target protein sequence and the candidate
 * protein sequence
 * 2007-07-16, Nanjiang Shu
 ****************************************************************************/
{
    int i,j ;
    /*first scan every diagnol line and record the line number */

    /*-----------------------------------------------------------------------------
     *  pair matrix
     *  +-----------------------+ 0-length_candidate
     *  | x    x    upper corner|
     *  |  x    x               |
     *  |   x    x              |
     *  |    x                  |
     *  |     x                 |
     *  |      x                |
     *  |       x               |
     *  |        x              |
     *  |    x    x             |
     *  |     x                 |
     *  | lower corner          |
     *  +-----------------------+
     *  0-length_target
     *-----------------------------------------------------------------------------*/
    Array1D <int> indexDiagUpperCorner_1darray(length_candidate);
    Array1D <int> indexDiagLowerCorner_1darray(length_target);
    int *indexDiagUpperCorner = indexDiagUpperCorner_1darray.array1D;
    int *indexDiagLowerCorner = indexDiagLowerCorner_1darray.array1D;

    int cntDiag = 0 ;
    cntDiag = 0;
    for(i = 0 ;  i < length_target; i ++)/*first scan diagnol lines in the lower corner*/
    {
        if(GetNumSatisfiedDot(pair_matrix, i, length_target, length_candidate, LOWER_CORNER) > 0)
        {
            indexDiagLowerCorner[cntDiag] = i;
            cntDiag ++;
        }
    }
    int numDiagLowerCorner = cntDiag;

    cntDiag = 0;
    for(j = 1 ;  j < length_candidate; j ++)/*then scan diagnol lines in the upper corner, j starting from 1, since the main diagnol has been used*/
    {
        if(GetNumSatisfiedDot(pair_matrix, j, length_target, length_candidate, UPPER_CORNER) > 0)
        {
            indexDiagUpperCorner[cntDiag] = j;
            cntDiag ++;
        }
    }
    int numDiagUpperCorner = cntDiag;

    if(numDiagUpperCorner+numDiagLowerCorner < 1)
    {
        numMaxDiag = 0;
        return numMaxDiag;
    }

#ifdef DEBUG_PRINT_MATRIX_SCORE
    fprintf(stdout,"print the matrix score for Method, numDiag = %d\n", numDiagLowerCorner+numDiagUpperCorner);

    Array1D <int> idxDiagPos_1darray(numDiagUpperCorner+numDiagLowerCorner);
    int *idxDiagPos = idxDiagPos_1darray.array1D;
    for(i = 0 ; i < numDiagLowerCorner; i ++) { idxDiagPos[i] = i; }
    QuickSort_index(idxDiagPos, indexDiagLowerCorner, 0 , numDiagLowerCorner-1);

    for(i = 0 ; i < numDiagLowerCorner; i ++)
    {
        int x,y;
        x = 0; 
        y = indexDiagLowerCorner[idxDiagPos[i]];

        while(x < length_candidate && y < length_target)
        {
            if(pair_matrix[y][x] != EMPTY)
            {
                fprintf(stdout,"diag = %3d, LowerCorner, dotScore[%3d] = %3d\n", indexDiagLowerCorner[idxDiagPos[i]], x, pair_matrix[y][x]);
            }
            x ++;
            y ++;
        }
    }
    for(i = 0 ; i < numDiagUpperCorner; i ++) { idxDiagPos[i] = i; }
    QuickSort_index(idxDiagPos, indexDiagUpperCorner, 0 , numDiagUpperCorner-1);
    for(i = 0 ; i < numDiagUpperCorner; i ++)
    {
        int x,y;
        x = indexDiagUpperCorner[idxDiagPos[i]]; 
        y = 0;

        while(x < length_candidate && y < length_target)
        {
            if(pair_matrix[y][x] != EMPTY)
            {
                fprintf(stdout,"diag = %3d, UpperCorner, dotScore[%3d] = %3d\n", indexDiagUpperCorner[idxDiagPos[i]], y, pair_matrix[y][x]);
            }
            x ++;
            y ++;
        }
    }
#endif

    /* for each selected diagnol line, calculating the score*/
    Array1D <double> score_1darray(numDiagLowerCorner+numDiagUpperCorner);
    double *score = score_1darray.array1D;

    Array1D <int> tmp_begPosTar_1darray(numDiagLowerCorner+numDiagUpperCorner);
    Array1D <int> tmp_endPosTar_1darray(numDiagLowerCorner+numDiagUpperCorner);
    Array1D <int> tmp_begPosCan_1darray(numDiagLowerCorner+numDiagUpperCorner);
    Array1D <int> tmp_endPosCan_1darray(numDiagLowerCorner+numDiagUpperCorner);
    int *tmp_begPosTar = tmp_begPosTar_1darray.array1D;
    int *tmp_endPosTar = tmp_endPosTar_1darray.array1D;
    int *tmp_begPosCan = tmp_begPosCan_1darray.array1D;
    int *tmp_endPosCan = tmp_endPosCan_1darray.array1D;


    for(i = 0 ; i < numDiagLowerCorner; i ++)
    {
        score[i] = GetDiagScore(indexDiagLowerCorner[i],LOWER_CORNER, indexDiagLowerCorner,numDiagLowerCorner,indexDiagUpperCorner,numDiagUpperCorner, pair_matrix, length_target, length_candidate, tmp_begPosTar[i], tmp_endPosTar[i], tmp_begPosCan[i], tmp_endPosCan[i]) ;
    }

    for(j = 0 ; j < numDiagUpperCorner; j ++)
    {
        score[j+numDiagLowerCorner] = GetDiagScore(indexDiagUpperCorner[j],UPPER_CORNER, indexDiagLowerCorner,numDiagLowerCorner, indexDiagUpperCorner,numDiagUpperCorner, pair_matrix, length_target, length_candidate, tmp_begPosTar[j+numDiagLowerCorner], tmp_endPosTar[j+numDiagLowerCorner], tmp_begPosCan[j+numDiagLowerCorner], tmp_endPosCan[j+numDiagLowerCorner]) ;
    }

    Array1D <int> idxScore_1darray(numDiagLowerCorner+numDiagUpperCorner);
    int *idxScore = idxScore_1darray.array1D;
    for(i = 0 ; i < numDiagUpperCorner+numDiagLowerCorner; i ++)
    {
        idxScore[i] = i;
    }
    QuickSort_index(idxScore, score, 0, numDiagLowerCorner+numDiagUpperCorner-1, DESCENDING);

    //double  maxScore = 0.0;
    
    int cntMaxDiag = 0;
    scoreDiag[cntMaxDiag] = score[idxScore[0]];
    if(idxScore[0] < numDiagLowerCorner)
    { 
        posMaxDiag[cntMaxDiag] = indexDiagLowerCorner[idxScore[0]]; 
        corner[cntMaxDiag] = LOWER_CORNER;
    }
    else
    { 
        posMaxDiag[cntMaxDiag] = indexDiagUpperCorner[idxScore[0]-numDiagLowerCorner];
        corner[cntMaxDiag] = UPPER_CORNER;
    }
    begPosTar[cntMaxDiag] = tmp_begPosTar[idxScore[0]];
    endPosTar[cntMaxDiag] = tmp_endPosTar[idxScore[0]];
    begPosCan[cntMaxDiag] = tmp_begPosCan[idxScore[0]];
    endPosCan[cntMaxDiag] = tmp_endPosCan[idxScore[0]];
    cntMaxDiag ++;

    //double lengthProp =   min(length_candidate, length_target) / double(max(length_target, length_candidate));
    double avgLen = double (length_candidate+length_target) /2.0;
    
    if (isFindingMultiDomain)
    {
        if(max(length_target, length_candidate) > 150) /*do not consider multidomain if both sequence are less than 150 aa long*/
        {
            for(i = 1 ; cntMaxDiag < max_num_diag && i < (numDiagUpperCorner+numDiagLowerCorner); i ++)
            {
                if(score[idxScore[i]] / score[idxScore[0]] < 0.2)
                {
                    break; /*if the line with maxScore is much higher than the next one, consider as single domain*/
                }
                bool isUniqDiag = true;
                for(j = 0; j < cntMaxDiag; j ++)
                {
                    scoreDiag[cntMaxDiag] = score[idxScore[i]];
                    if(idxScore[i] < numDiagLowerCorner)
                    { 
                        posMaxDiag[cntMaxDiag] = indexDiagLowerCorner[idxScore[i]]; 
                        corner[cntMaxDiag] = LOWER_CORNER;
                    }
                    else
                    { 
                        posMaxDiag[cntMaxDiag] = indexDiagUpperCorner[idxScore[i]-numDiagLowerCorner];
                        corner[cntMaxDiag] = UPPER_CORNER;
                    }
                    begPosTar[cntMaxDiag] = tmp_begPosTar[idxScore[i]];
                    endPosTar[cntMaxDiag] = tmp_endPosTar[idxScore[i]];
                    begPosCan[cntMaxDiag] = tmp_begPosCan[idxScore[i]];
                    endPosCan[cntMaxDiag] = tmp_endPosCan[idxScore[i]];

                    /*now check if the line cntMaxDiag is uniq to the previous added
                     * lines*/
                    int dist;
                    if(corner[cntMaxDiag] == corner[j])
                    { dist = abs(posMaxDiag[cntMaxDiag] - posMaxDiag[j]); }
                    else 
                    { dist = posMaxDiag[cntMaxDiag] + posMaxDiag[j] -1; }
                    if(dist/avgLen < 0.3)/* the distance between this line and the previous lines should be big enought to be considered as unique*/
                    {
                        isUniqDiag = false;
                        break;
                    }
                }
                if(isUniqDiag)
                {
                    cntMaxDiag ++;
                }
            }
        }
    }
    
    numMaxDiag = cntMaxDiag;

    
    /*unify scoreDiag*/
    for(i = 0 ; i < numMaxDiag; i ++)
    {
        //if(length_target < 150 && length_candidate < 150 && lengthProp >= 0.7) // that should be single domain
        //{
        //    if(posMaxDiag[i] /avgLen > 0.3)         [>this solve 1EWWA-1TBN<]
        //    {
        //        //scoreDiag[i] /= (avgLen / posMaxDiag[i]);
        //        avgLen -= posMaxDiag[i]; [>if the alignment shift a lot from the main diagnol, decrease the average length,2007-11-18, Nanjiang, the code before is probably wrong<]
        //    }
        //}
        if(numMaxDiag <=1)
        {
            scoreDiag[i] /= (avgLen);
        }
        else
        {
            scoreDiag[i] /= double(min(length_candidate,length_target));
        }
    }
    return numMaxDiag ;
}
/*}}}*/
int  ScanPairDotPlot_Method2(int *posTar, int *posCan, float *fragScore, int numPair, int length_target, int length_candidate, int *posMaxDiag, int *corner, double *scoreDiag, int *begPosTar, int *endPosTar, int *begPosCan, int *endPosCan, int &numMaxDiag, int weightScheme=0)/*{{{*/
/*****************************************************************************
 * given the posTar and posCan and numPair, calculate the diagnol with the
 * highest similarity score. return the position of the diagnol and the score
 * 2007-11-12, Nanjiang
 * Note that in this algorithm, the BuildPairMatrix subroutine is not run
 ****************************************************************************/
{
    int i,j ;
    /*first scan every diagnol line and record the line number */
    /*-----------------------------------------------------------------------------
     *  pair matrix
     *  +-----------------------+ 0-length_candidate
     *  | x    x    upper corner|
     *  |  x    x               |
     *  |   x    x              |
     *  |    x                  |
     *  |     x                 |
     *  |      x                |
     *  |       x               |
     *  |        x              |
     *  |    x    x             |
     *  |     x                 |
     *  | lower corner          |
     *  +-----------------------+
     *  0-length_target
     *-----------------------------------------------------------------------------*/
    int MAX_NUM_DIAGNOL = min(length_target+length_candidate, numPair);
    Array1D <int> indexDiag_1darray(MAX_NUM_DIAGNOL);
    Array1D <int> indexDiagCorner_1darray(MAX_NUM_DIAGNOL);
    int *indexDiag = indexDiag_1darray.array1D;
    int *indexDiagCorner = indexDiagCorner_1darray.array1D;

    Array1D <DATATYPE_DOT_SCORE*> diag_1darray(MAX_NUM_DIAGNOL);
    Array1D <float*> fragScoreDiag_1darray(MAX_NUM_DIAGNOL);
    Array1D <int*> indexFilledDotDiag_1darray(MAX_NUM_DIAGNOL);
    Array1D <int> sizeDiag_1darray(MAX_NUM_DIAGNOL);
    Array1D <int> numFilledDotDiag_1darray(MAX_NUM_DIAGNOL);
    DATATYPE_DOT_SCORE **diag = diag_1darray.array1D; /*diag[i] is a 1d array, recording the filling in the diagnol */
    float **fragScoreDiag = fragScoreDiag_1darray.array1D; /*fragScore of each diagnol*/
    int **indexFilledDotDiag = indexFilledDotDiag_1darray.array1D; 
    int *sizeDiag = sizeDiag_1darray.array1D;
    int *numFilledDotDiag = numFilledDotDiag_1darray.array1D;
    
    /*each diagnol is defined by the index and corner, say, 12, UPPER_CORNER*/

    int cntDiag = 0 ;
    cntDiag = 0;
    for(i = 0 ;  i < numPair; i ++)/*add all dots to different diagnols*/
    {
        int idxDiag = 0;
        int pos = 0;
        if((idxDiag = FindDiagIndex(posCan[i], posTar[i], indexDiag, indexDiagCorner, cntDiag)) == -1) /*create a new diagnol*/
        {
            sizeDiag[cntDiag] = GetDiagPos(posCan[i], posTar[i], length_candidate, length_target, indexDiag[cntDiag], indexDiagCorner[cntDiag]);
            assert(sizeDiag[cntDiag] > 0);
            diag[cntDiag] = new DATATYPE_DOT_SCORE [sizeDiag[cntDiag]];
            fragScoreDiag[cntDiag] = new float [sizeDiag[cntDiag]];
            InitDiag(diag[cntDiag], sizeDiag[cntDiag]);
            idxDiag = cntDiag;
            numFilledDotDiag[idxDiag] = 0;
            cntDiag ++;
        }
        /*add this dot to the diagnol*/
        pos = Map2D_1D(posCan[i],posTar[i], indexDiagCorner[idxDiag]);
        diag[idxDiag][pos] = FILLED;
        fragScoreDiag[idxDiag][pos] = fragScore[i]; /*input the fragScore for dot i to fragScoreDiag, 2008-01-16, Nanjiang*/
        numFilledDotDiag[idxDiag] ++;
    }
    int numDiag = cntDiag;

    for (i = 0 ; i < numDiag; i ++)
    {
        indexFilledDotDiag[i] = new int[numFilledDotDiag[i]];
        int cnt = 0;
        for(j = 0; j < sizeDiag[i]; j ++)
        {
            if(diag[i][j] != EMPTY)
            {
                indexFilledDotDiag[i][cnt] = j;
                cnt ++;
            }
        }
        if (cnt != numFilledDotDiag[i])
        {
            fprintf(stderr,"cnt = %d, numFilledDotDiag[%d] = %d\n", cnt, i, numFilledDotDiag[i]);
            assert(cnt == numFilledDotDiag[i]);
        }
    }

    int numValidDiag = 0; /*number of valid diagnols, some diagnols which have too few filled dots are removed*/
    Array1D <DATATYPE_DOT_SCORE*> validDiag_1darray(numDiag);
    Array1D <float*> fragScoreValidDiag_1darray(numDiag);
    Array1D <int> indexValidDiag_1darray(numDiag);
    Array1D <int> indexValidDiagCorner_1darray(numDiag);
    DATATYPE_DOT_SCORE **validDiag = validDiag_1darray.array1D; /*validDiag and fragScoreValidDiag do not allocate memory to the second level of the array*/
    float **fragScoreValidDiag = fragScoreValidDiag_1darray.array1D;
    int *indexValidDiag = indexValidDiag_1darray.array1D;
    int *indexValidDiagCorner = indexValidDiagCorner_1darray.array1D;
    Array1D <int> sizeValidDiag_1darray(numDiag);
    Array1D <int*> indexFilledDotValidDiag_1darray(numDiag);
    Array1D <int> numFilledDotValidDiag_1darray(numDiag);
    int * sizeValidDiag = sizeValidDiag_1darray.array1D;
    int ** indexFilledDotValidDiag = indexFilledDotValidDiag_1darray.array1D;
    int * numFilledDotValidDiag = numFilledDotValidDiag_1darray.array1D;
    int cntValidDiag = 0;
    for(i = 0; i < numDiag; i ++)
    {

        if (numFilledDotDiag[i] < 2 )
        {
            continue;
            //double(numFilledDotDiag[i]) / sizeDiag[i] < 0.1
        }
        /*checking if the maxSegSize > 3*/
        j = 0;
        int maxSegSize = 1;
        while(j < numFilledDotDiag[i] -1 )
        {
            if(indexFilledDotDiag[i][j+1] - indexFilledDotDiag[i][j] == 1)
            {
                maxSegSize ++;
            }
            else 
            {
                maxSegSize = 1;
            }
            if (maxSegSize >= 3) break;
            j ++;
        }
        if (maxSegSize >= 3 || double(numFilledDotDiag[i]) / sizeDiag[i] >= 0.1 )
        {
            validDiag[cntValidDiag] = diag[i];
            fragScoreValidDiag[cntValidDiag] = fragScoreDiag[i]; /*added 2008-01-16, Nanjiang*/
            indexValidDiag[cntValidDiag] = indexDiag[i];
            indexValidDiagCorner[cntValidDiag] = indexDiagCorner[i];
            indexFilledDotValidDiag[cntValidDiag] = indexFilledDotDiag[i];
            sizeValidDiag[cntValidDiag] = sizeDiag[i];
            numFilledDotValidDiag[cntValidDiag] = numFilledDotDiag[i];
            cntValidDiag ++;
        }
    }
    numValidDiag = cntValidDiag;
    

    if( numValidDiag < 1)
    {
        numMaxDiag = 0;
//        return numMaxDiag; /*this is a bug, if return here, the memory will not be freed, 2008-04-12*/
    }
    else/*{{{*/
    {
        /*initialize the score of dot on each diagnol*/
        for (i = 0; i < numValidDiag; i ++)
        {
            SetDiagDotScore(validDiag[i], fragScoreValidDiag[i], sizeValidDiag[i], numFilledDotValidDiag[i], weightScheme);
        }

#ifdef DEBUG_PRINT_MATRIX_SCORE
        Array1D <int> idxDiagPos_1darray(numValidDiag);
        int *idxDiagPos = idxDiagPos_1darray.array1D;
        for(i = 0 ; i < numValidDiag; i ++) { idxDiagPos[i] = i; }
        QuickSort_index(idxDiagPos, indexValidDiag, 0 , numValidDiag-1);
        fprintf(stdout,"print the matrix score for Method, numDiag = %d\n", numValidDiag);
        for(i = 0 ; i < numValidDiag; i ++)
        {
            if(indexValidDiagCorner[idxDiagPos[i]] == LOWER_CORNER)
            {
                for(j = 0; j < sizeValidDiag[idxDiagPos[i]]; j ++)
                {
                    if(validDiag[idxDiagPos[i]][j] != EMPTY)
                    {
                        fprintf(stdout,"diag = %3d, %s, dotScore[%3d] = %3d\n", indexValidDiag[idxDiagPos[i]], ((indexValidDiagCorner[idxDiagPos[i]] == UPPER_CORNER) ? "UpperCorner" : "LowerCorner"), j, validDiag[idxDiagPos[i]][j]);
                    }
                }
            }
        }
        for(i = 0 ; i < numValidDiag; i ++)
        {
            if(indexValidDiagCorner[idxDiagPos[i]] == UPPER_CORNER)
            {
                for(j = 0; j < sizeValidDiag[idxDiagPos[i]]; j ++)
                {
                    if(validDiag[idxDiagPos[i]][j] != EMPTY)
                    {
                        fprintf(stdout,"diag = %3d, %s, dotScore[%3d] = %3d\n", indexValidDiag[idxDiagPos[i]], ((indexValidDiagCorner[idxDiagPos[i]] == UPPER_CORNER) ? "UpperCorner" : "LowerCorner"), j, validDiag[idxDiagPos[i]][j]);
                    }
                }
            }
        }
#endif
        /* for each selected diagnol line, calculating the score*/
        Array1D <double> score_1darray(numValidDiag);
        double *score = score_1darray.array1D;
        for(i = 0; i < numValidDiag; i ++)
        {
            score[i] = 0.0;
        }

        Array1D <int> tmp_begPosTar_1darray(numValidDiag);
        Array1D <int> tmp_endPosTar_1darray(numValidDiag);
        Array1D <int> tmp_begPosCan_1darray(numValidDiag);
        Array1D <int> tmp_endPosCan_1darray(numValidDiag);
        int *tmp_begPosTar = tmp_begPosTar_1darray.array1D;
        int *tmp_endPosTar = tmp_endPosTar_1darray.array1D;
        int *tmp_begPosCan = tmp_begPosCan_1darray.array1D;
        int *tmp_endPosCan = tmp_endPosCan_1darray.array1D;

        double lengthProp =   min(length_candidate, length_target) / double(max(length_target, length_candidate));
        double avgLen = double (length_candidate+length_target) /2.0;


        for(i = 0 ; i < numValidDiag; i ++)
        {
            double diagShift = 0.0;
            if(indexValidDiagCorner[i] == UPPER_CORNER) { diagShift = indexValidDiag[i] / double(length_candidate); }
            else { diagShift = indexValidDiag[i] / double(length_target); }

            if(isIgnoreNonDiag && diagShift >= 0.30 && (abs(length_target-length_candidate) <= 25 || lengthProp >= 0.7))
            {
                score[i] = -5000000.0; 
            }
            else
            {
                score[i] = GetDiagScore2(i, validDiag, indexValidDiag, indexValidDiagCorner, sizeValidDiag, indexFilledDotValidDiag, numFilledDotValidDiag, numValidDiag, length_candidate, length_target, tmp_begPosCan[i], tmp_endPosCan[i], tmp_begPosTar[i], tmp_endPosTar[i]) ;
            }
        }


        Array1D <int> idxScore_1darray(numValidDiag);
        int *idxScore = idxScore_1darray.array1D;
        for(i = 0 ; i < numValidDiag; i ++) { idxScore[i] = i; }
        QuickSort_index(idxScore, score, 0, numValidDiag-1, DESCENDING);

        //double  maxScore = 0.0;
        int cntMaxDiag = 0;
        scoreDiag[cntMaxDiag] = score[idxScore[0]];
        posMaxDiag[cntMaxDiag] = indexValidDiag[idxScore[0]]; /*bug fixed, using validDiag*/
        corner[cntMaxDiag] = indexValidDiagCorner[idxScore[0]];

        begPosTar[cntMaxDiag] = tmp_begPosTar[idxScore[0]];
        endPosTar[cntMaxDiag] = tmp_endPosTar[idxScore[0]];
        begPosCan[cntMaxDiag] = tmp_begPosCan[idxScore[0]];
        endPosCan[cntMaxDiag] = tmp_endPosCan[idxScore[0]];
        cntMaxDiag ++;


        if (isFindingMultiDomain) /*finding multidomain is run only when isFindingMultiDomain = true, 2007-11-18, Nanjiang*//*{{{*/
        {
            if(max(length_target, length_candidate) > 150) /*do not consider multidomain if both sequence are less than 150 aa long*/
            {
                for(i = 1 ; cntMaxDiag < max_num_diag && i < (numValidDiag); i ++)
                {
                    if(score[idxScore[i]] / score[idxScore[0]] < 0.2)
                    {
                        break; /*if the line with maxScore is much higher than the next one, consider as single domain*/
                    }
                    bool isUniqDiag = true;
                    for(j = 0; j < cntMaxDiag; j ++)
                    {
                        scoreDiag[cntMaxDiag] = score[idxScore[i]];
                        posMaxDiag[cntMaxDiag] = indexValidDiag[idxScore[i]]; 
                        corner[cntMaxDiag] = indexValidDiagCorner[idxScore[i]];
                        begPosTar[cntMaxDiag] = tmp_begPosTar[idxScore[i]];
                        endPosTar[cntMaxDiag] = tmp_endPosTar[idxScore[i]];
                        begPosCan[cntMaxDiag] = tmp_begPosCan[idxScore[i]];
                        endPosCan[cntMaxDiag] = tmp_endPosCan[idxScore[i]];

                        /*now check if the line cntMaxDiag is uniq to the previous added
                         * lines*/
                        int dist;
                        if(corner[cntMaxDiag] == corner[j])
                        { dist = abs(posMaxDiag[cntMaxDiag] - posMaxDiag[j]); }
                        else 
                        { dist = posMaxDiag[cntMaxDiag] + posMaxDiag[j] -1; }
                        if(dist/avgLen < 0.3)/* the distance between this line and the previous lines should be big enought to be considered as unique*/
                        {
                            isUniqDiag = false;
                            break;
                        }
                    }
                    if(isUniqDiag)
                    {
                        cntMaxDiag ++;
                    }
                }
            }
        }/*}}}*/

        numMaxDiag = cntMaxDiag;
        /*unify scoreDiag*/
        for(i = 0 ; i < numMaxDiag; i ++)
        {
            //        if(length_target < 150 && length_candidate < 150 && lengthProp >= 0.7) // that should be single domain
            //        {
            //            if(posMaxDiag[i] /avgLen > 0.3)         [>this solve 1EWWA-1TBN<]
            //            {
            //                //scoreDiag[i] /= (avgLen / posMaxDiag[i]);
            ////                avgLen -= posMaxDiag[i]; [>if the alignment shift a lot from the main diagnol, decrease the average length,2007-11-18, Nanjiang, the code before is probably wrong<]
            //            }
            //        }
            if(numMaxDiag <=1)
            {
                scoreDiag[i] /= (avgLen);
            }
            else
            {
                scoreDiag[i] /= double(min(length_candidate,length_target));
            }
            /*give extra points to the domains with similar length, 2007-11-18*/
            if( abs(length_target-length_candidate) <= 25 || lengthProp >= 0.7)
            {
                if(max(length_candidate, length_target) < 60)
                {
                    scoreDiag[i] *= (1.5 - (abs(length_target-length_candidate)) / avgLen);
                }
                else if (lengthProp >= 0.7)
                {
                    scoreDiag[i] *= (lengthProp+0.5 );
                }
            }
        }
    }/*}}}*/

    //free memory
    for(i = 0 ; i < numDiag; i ++)
    {
        delete [] diag[i];
        delete [] fragScoreDiag[i];
        delete [] indexFilledDotDiag[i];
    }

    return numMaxDiag ;
}
/*}}}*/
//int  ScanPairDotPlot_Method3(int *posTar, int *posCan, int numPair, int length_target, int length_candidate, int *posMaxDiag, int *corner, double *scoreDiag, int *begPosTar, int *endPosTar, int *begPosCan, int *endPosCan, int &numMaxDiag)[>{{{<]
/*****************************************************************************
 * given the posTar and posCan and numPair, calculate the diagnol with the
 * highest similarity score. return the position of the diagnol and the score
 * 2007-11-12, Nanjiang
 * Note that in this algorithm, the BuildPairMatrix subroutine is not run
 ****************************************************************************/
//{
//    int i,j ;
//    [>first scan every diagnol line and record the line number <]
    /*-----------------------------------------------------------------------------
     *  pair matrix
     *  +-----------------------+ 0-length_candidate
     *  | x    x    upper corner|
     *  |  x    x               |
     *  |   x    x              |
     *  |    x                  |
     *  |     x                 |
     *  |      x                |
     *  |       x               |
     *  |        x              |
     *  |    x    x             |
     *  |     x                 |
     *  | lower corner          |
     *  +-----------------------+
     *  0-length_target
     *-----------------------------------------------------------------------------*/
//    int MAX_NUM_DIAGNOL = min(length_target+length_candidate, numPair);
//    Array1D <int> indexDiag_1darray(MAX_NUM_DIAGNOL);
//    Array1D <int> indexDiagCorner_1darray(MAX_NUM_DIAGNOL);
//    int *indexDiag = indexDiag_1darray.array1D;
//    int *indexDiagCorner = indexDiagCorner_1darray.array1D;

//    Array1D <int*> diagDotPos_1darray(MAX_NUM_DIAGNOL);
//    Array1D <int> diagNumFilledDot_1darray(MAX_NUM_DIAGNOL);
//    Array1D <int> sizeDiag_1darray(MAX_NUM_DIAGNOL);
//    int **diagDotPos = diagDotPos_1darray.array1D; [>diag[i] is a 1d array, recording the filling in the diagnol <]
//    int *sizeDiag = sizeDiag_1darray.array1D;
//    int *diagNumFilledDot = diagNumFilledDot_1darray.array1D;
    
//    [>each diagnol is defined by the index and corner, say, 12, UPPER_CORNER<]

//    int cntDiag = 0 ;
//    cntDiag = 0;
//    for(i = 0 ;  i < numPair; i ++)[>add all dots to different diagnols<]
//    {
//        int idxDiag = 0;
//        int pos = 0;
//        if((idxDiag = FindDiagIndex(posCan[i]. posTar[i], indexDiag, indexDiagCorner, cntDiag)) == -1) [>create a new diagnol<]
//        {
//            sizeDiag [cntDiag] = GetDiagIndex(posCan[i], postTar[i], length_candidate, length_target, indexDiag[cntDiag], indexDiagCorner[cntDiag]);
//            assert(sizeDiag[cntDiag] > 0);
//            diagDotPos[cntDiag] = new int [sizeDiag[cntDiag]];
//            //InitDiag(diag[cntDiag], sizeDiag[cntDiag]);
//            idxDiag = cntDiag;
//            cntDiag ++;
//        }
//        [>add this dot to the diagnol<]
//        pos = Map2D_1D(posCan[i],posTar[i], indexDiagCorner[idxDiag]);
//        diag[idxDiag][pos] = FILLED;
//    }
//    int numDiag = cntDiag;

//    if( numDiag < 1)
//    {
//        numMaxDiag = 0;
//        return numMaxDiag;
//    }

//    [> for each selected diagnol line, calculating the score<]
//    Array1D <double> score_1darray(numDiagLowerCorner+numDiagUpperCorner);
//    double *score = score_1darray.array1D;

//    Array1D <int> tmp_begPosTar_1darray(numDiagLowerCorner+numDiagUpperCorner);
//    Array1D <int> tmp_endPosTar_1darray(numDiagLowerCorner+numDiagUpperCorner);
//    Array1D <int> tmp_begPosCan_1darray(numDiagLowerCorner+numDiagUpperCorner);
//    Array1D <int> tmp_endPosCan_1darray(numDiagLowerCorner+numDiagUpperCorner);
//    int *tmp_begPosTar = tmp_begPosTar_1darray.array1D;
//    int *tmp_endPosTar = tmp_endPosTar_1darray.array1D;
//    int *tmp_begPosCan = tmp_begPosCan_1darray.array1D;
//    int *tmp_endPosCan = tmp_endPosCan_1darray.array1D;


//    for(i = 0 ; i < numDiagLowerCorner; i ++)
//    {
//        score[i] = GetDiagScore(indexDiagLowerCorner[i],LOWER_CORNER, indexDiagLowerCorner,numDiagLowerCorner,indexDiagUpperCorner,numDiagUpperCorner, pair_matrix, length_target, length_candidate, tmp_begPosTar[i], tmp_endPosTar[i], tmp_begPosCan[i], tmp_endPosCan[i]) ;
//    }

//    for(j = 0 ; j < numDiagUpperCorner; j ++)
//    {
//        score[j+numDiagLowerCorner] = GetDiagScore(indexDiagUpperCorner[j],UPPER_CORNER, indexDiagLowerCorner,numDiagLowerCorner, indexDiagUpperCorner,numDiagUpperCorner, pair_matrix, length_target, length_candidate, tmp_begPosTar[j+numDiagLowerCorner], tmp_endPosTar[j+numDiagLowerCorner], tmp_begPosCan[j+numDiagLowerCorner], tmp_endPosCan[j+numDiagLowerCorner]) ;
//    }

//    Array1D <int> idxScore_1darray(numDiagLowerCorner+numDiagUpperCorner);
//    int *idxScore = idxScore_1darray.array1D;
//    for(i = 0 ; i < numDiagUpperCorner+numDiagLowerCorner; i ++)
//    {
//        idxScore[i] = i;
//    }
//    QuickSort_index(idxScore, score, 0, numDiagLowerCorner+numDiagUpperCorner-1, DESCENDING);

//    double  maxScore = 0.0;
    
//    int cntMaxDiag = 0;
//    scoreDiag[cntMaxDiag] = score[idxScore[0]];
//    if(idxScore[0] < numDiagLowerCorner)
//    { 
//        posMaxDiag[cntMaxDiag] = indexDiagLowerCorner[idxScore[0]]; 
//        corner[cntMaxDiag] = LOWER_CORNER;
//    }
//    else
//    { 
//        posMaxDiag[cntMaxDiag] = indexDiagUpperCorner[idxScore[0]-numDiagLowerCorner];
//        corner[cntMaxDiag] = UPPER_CORNER;
//    }
//    begPosTar[cntMaxDiag] = tmp_begPosTar[idxScore[0]];
//    endPosTar[cntMaxDiag] = tmp_endPosTar[idxScore[0]];
//    begPosCan[cntMaxDiag] = tmp_begPosCan[idxScore[0]];
//    endPosCan[cntMaxDiag] = tmp_endPosCan[idxScore[0]];
//    cntMaxDiag ++;

//    double lengthProp =   min(length_candidate, length_target) / double(max(length_target, length_candidate));
//    double avgLen = double (length_candidate+length_target) /2.0;
    
//    if(max(length_target, length_candidate) > 150) [>do not consider multidomain if both sequence are less than 150 aa long<]
//    {
//        for(i = 1 ; cntMaxDiag < max_num_diag && i < (numDiagUpperCorner+numDiagLowerCorner); i ++)
//        {
//            if(score[idxScore[i]] / score[idxScore[0]] < 0.2)
//            {
//                break; [>if the line with maxScore is much higher than the next one, consider as single domain<]
//            }
//            bool isUniqDiag = true;
//            for(j = 0; j < cntMaxDiag; j ++)
//            {
//                scoreDiag[cntMaxDiag] = score[idxScore[i]];
//                if(idxScore[i] < numDiagLowerCorner)
//                { 
//                    posMaxDiag[cntMaxDiag] = indexDiagLowerCorner[idxScore[i]]; 
//                    corner[cntMaxDiag] = LOWER_CORNER;
//                }
//                else
//                { 
//                    posMaxDiag[cntMaxDiag] = indexDiagUpperCorner[idxScore[i]-numDiagLowerCorner];
//                    corner[cntMaxDiag] = UPPER_CORNER;
//                }
//                begPosTar[cntMaxDiag] = tmp_begPosTar[idxScore[i]];
//                endPosTar[cntMaxDiag] = tmp_endPosTar[idxScore[i]];
//                begPosCan[cntMaxDiag] = tmp_begPosCan[idxScore[i]];
//                endPosCan[cntMaxDiag] = tmp_endPosCan[idxScore[i]];

                /*now check if the line cntMaxDiag is uniq to the previous added
                 * lines*/
//                int dist;
//                if(corner[cntMaxDiag] == corner[j])
//                { dist = abs(posMaxDiag[cntMaxDiag] - posMaxDiag[j]); }
//                else 
//                { dist = posMaxDiag[cntMaxDiag] + posMaxDiag[j] -1; }
//                if(dist/avgLen < 0.3)[> the distance between this line and the previous lines should be big enought to be considered as unique<]
//                {
//                    isUniqDiag = false;
//                    break;
//                }
//            }
//            if(isUniqDiag)
//            {
//                cntMaxDiag ++;
//            }
//        }
//    }
    
//    numMaxDiag = cntMaxDiag;

    
//    [>unify scoreDiag<]
//    for(i = 0 ; i < numMaxDiag; i ++)
//    {
//        if(length_target < 150 && length_candidate < 150 && lengthProp >= 0.7) // that should be single domain
//        {
//            if(posMaxDiag[i] /avgLen > 0.3)         [>this solve 1EWWA-1TBN<]
//            {
//                scoreDiag[i] /= (avgLen / posMaxDiag[i]);
//            }
//        }
//        if(numMaxDiag <=1)
//        {
//            scoreDiag[i] /= (avgLen);
//        }
//        else
//        {
//            scoreDiag[i] /= double(min(length_candidate,length_target));
//        }
//    }
//    return numMaxDiag ;
//}
//[>}}}<]
int PostScanFragSearch(const char *fragPairFile, const char *postFragPairFile, char **proIDListCan, int *proLengthCan, int numProCan, char **proIDListTar, int *proLengthTar, int numProTar, int method, int weightScheme)/*{{{*/
{
    int i,j ;
#ifdef DEBUG_TIME
    clock_t start, finish;
    double duration;
    int tmp;
    start = clock();
#endif
    char rtname [MAX_PATH+1] = "";
    char idTar[SIZE_ID+1] = "";
    char idCan[SIZE_ID+1] = "";
    rootname(fragPairFile, rtname);
    my_strcpy(idTar, rtname, SIZE_ID);
    int idx = 0;

    int length_target = 0; 
    int length_candidate = 0; 

    idx = BinarySearch_String(idTar, proIDListTar, numProTar);
    if(idx != -1)
    {
        length_target = proLengthTar[idx];
    }
    else
    {
        fprintf(stderr,"idTar=%s, not found in proIDListTar\n", idTar);
        assert(idx >= 0);
    }

    //int linesize;
    //int maxline = 300;
    //Array1D <char> line_1darray(maxline+1);
    //char *line = line_1darray.array1D;


	int max_num_pair = 0;
	if (isReadBinaryFragfile)
	{
		int tmpi1 = 0;
		ReadInFragParameter((char*)fragPairFile, tmpi1, tmpi1, tmpi1, tmpi1, max_num_pair);
	}
	else
	{
		max_num_pair = fgetlinecnt(fragPairFile);
	}
    Array1D <int> allPosCan_1darray(max_num_pair);
    Array1D <int> allPosTar_1darray(max_num_pair);
    Array1D <float> allFragScore_1darray(max_num_pair);
    Array2D <char> allProIDCan_2darray(max_num_pair, SIZE_ID+1);
    int *allPosCan = allPosCan_1darray.array1D;
    int *allPosTar = allPosTar_1darray.array1D;
    float *allFragScore = allFragScore_1darray.array1D;
    char **allProIDCan = allProIDCan_2darray.array2D;

#ifdef DEBUG_TIME
    start = clock();
#endif

    int numAllPair = 0;
	if (isReadBinaryFragfile)
	{
		numAllPair = ReadInPairBinary(fragPairFile, allPosTar, allPosCan, allFragScore, allProIDCan); /*read in the binary frag file 2009-06-25, Nanjiang*/
	}
	else
	{
		numAllPair = ReadInPair(fragPairFile, allPosTar, allPosCan, allFragScore, allProIDCan); /*read allFragScore as well, 2008-01-15, Nanjiang*/
	}
    if(weightScheme == 1) /*weightScheme == 1, using the weight from frag-frag scores*/
    {
        float meanFragScore = 0.0;
        float medianFragScore = 0.0;
        float minFragScore = 0.0;
        float maxFragScore = 0.0;
        AnalyzeFragScore(allFragScore, numAllPair, meanFragScore, medianFragScore, minFragScore, maxFragScore, cutoff_fragscore);
#ifdef DEBUG_WEIGHT
        fprintf(stdout,"%s: meanFragScore= %.0f , minFragScore= %.0f , maxFragScore= %.0f , cutoff_fragscore= %.0f, numAllPair=%d\n", idTar, meanFragScore, minFragScore, maxFragScore, cutoff_fragscore, numAllPair);
        return 0;
#endif
//        RemoveFragByCutoffScore(allPosTar, allPosCan, allFragScore, allProIDCan, numAllPair, cutoff_fragscore);
        NormalizeFragScore(allFragScore, meanFragScore, minFragScore, maxFragScore, numAllPair); /*scale the allFragScore to 0~2, so that it can be used as weight, 2008-01-16,Nanjiang*/
    }

#ifdef DEBUG_TIME
    finish = clock();
    duration = double(finish-start)  /double(CLOCKS_PER_SEC);
    fprintf(stdout,"ReadInPair %s cost %lf seconds\n", fragPairFile, duration);
#endif


#ifdef DEBUG
    FILE *fpt1 = fopen("/tmp/t1.txt", "w");
    for(int it1 = 0 ; it1 < numAllPair; it1 ++)
    {
        fprintf(fpt1, "%d %s %d\n", allPosTar[it1], allProIDCan[it1], allPosCan[it1]);
    }
    fclose(fpt1);
#endif

#ifdef DEBUG_TIME
    start = clock();
#endif
    Array2D <char> uniqProIDCan_2darray(numAllPair, SIZE_ID+1);
    Array1D <int>  subTotalProIDCan_1darray(numAllPair);
    char ** uniqProIDCan = uniqProIDCan_2darray.array2D;
    int *subTotalProIDCan = subTotalProIDCan_1darray.array1D;
    int numUniqProIDCan = Grouping(allProIDCan, numAllPair, uniqProIDCan, subTotalProIDCan, SIZE_ID); // grouping the candidate segment according to the proID of candicate
    
#ifdef DEBUG_PRINT_ID
    for(i = 0 ; i < numUniqProIDCan; i ++)
    {
        fprintf(stdout, "%-4d %s count = %3d\n", i, uniqProIDCan[i], subTotalProIDCan[i] );
    }
    fprintf(stdout,"\n");
#endif

#ifdef DEBUG_TIME
    finish = clock();
    duration = double(finish-start)  /double(CLOCKS_PER_SEC);
    fprintf(stdout,"Grouping %s cost %lf seconds\n", fragPairFile, duration);
#endif

#ifdef DEBUG_TIME
    start = clock();
#endif
    Array1D <int*> indexUniqProInAllProID_1darray(numUniqProIDCan);
    int **indexUniqProInAllProID = indexUniqProInAllProID_1darray.array1D;
    Array1D <int> cntUniqProInAllProID_1darray(numUniqProIDCan);
    int *cntUniqProInAllProID = cntUniqProInAllProID_1darray.array1D;
    for(i = 0; i < numUniqProIDCan; i ++)
    {
        indexUniqProInAllProID[i] = new int [subTotalProIDCan[i]];
        cntUniqProInAllProID[i]  = 0;
    }
    for(i = 0 ; i < numAllPair; i ++)
    {
        idx = BinarySearch_String(allProIDCan[i], uniqProIDCan, numUniqProIDCan);
        if (idx >= 0)
        {
            indexUniqProInAllProID[idx][cntUniqProInAllProID[idx]] = i;
            cntUniqProInAllProID[idx] ++;
        }
    }
#ifdef DEBUG_TIME
    finish = clock();
    duration = double(finish-start)  /double(CLOCKS_PER_SEC);
    fprintf(stdout,"Getting indexUniqProInAllProID %s cost %lf seconds\n", fragPairFile, duration);
#endif

#ifdef DEBUG_TIME
    start = clock();
#endif
    /*determine the raw score*/
    Array1D <double> rawPerArray_1darray(numUniqProIDCan+1);
    double *rawPerArray = rawPerArray_1darray.array1D;
    for(i = 0 ; i < numUniqProIDCan; i ++)// determine the cutoff_rawper
    {
        idx = BinarySearch_String(uniqProIDCan[i], proIDListCan, numProCan);
        if(idx != -1)
        {
            length_candidate = proLengthCan[idx];
        }
        else
        {
            fprintf(stderr,"idCan=%s, not found in proIDListCan\n", uniqProIDCan[i]);
            assert(idx >= 0);
        }
        rawPerArray[i] = (subTotalProIDCan[i] / double (length_candidate + length_target)) * 2 * 100;
    }
#ifdef DEBUG_PRINT
    fprintf(stdout,"rawPerArray, numUniqProIDCan = %d\n", numUniqProIDCan);
    for(i = 0 ;i < numUniqProIDCan; i ++)
    {
        fprintf(stdout,"rawPer[%d] = %lf subTotalProIDCan = %d  %s -- %s\n", i, rawPerArray[i], subTotalProIDCan[i], idTar, uniqProIDCan[i]);
    }
#endif
    QuickSort(rawPerArray, 0, numUniqProIDCan-1, DESCENDING);

#ifdef DEBUG
    int a[] = { 0,2 ,3, 32, 2,3};
    QuickSort(a, 0, 5, DESCENDING);
    QuickSort(a, 0, 5, ASCENDING);
    int tmpidx[] = { 0,1,2,3,4,5,6,7,8 };
    QuickSort_index(tmpidx, a, 0, 5, DESCENDING);
    for(i = 0 ; i < 6;  i ++) { printf("%d\n", a[tmpidx[i]]); }
#endif
    
    int pivatIndex = min(numUniqProIDCan-1, 50); // the raw score of the 50th and the last one
    if (cutoff_rawper == -1.0)
    {
        cutoff_rawper = max(rawPerArray[pivatIndex], 10.0);
    }

#ifdef DEBUG_TIME
    finish = clock();
    duration = double(finish-start)  /double(CLOCKS_PER_SEC);
    fprintf(stdout,"Determining the pivatRawScore %s cost %lf seconds\n", fragPairFile, duration);
#endif
    Array2D <char> domainID_2darray(numUniqProIDCan*max_num_diag, SIZE_ID+1);
    char ** domainID = domainID_2darray.array2D;

    Array1D <int> begPosTar_1darray(numUniqProIDCan*max_num_diag);
    Array1D <int> endPosTar_1darray(numUniqProIDCan*max_num_diag);
    Array1D <int> begPosCan_1darray(numUniqProIDCan*max_num_diag);
    Array1D <int> endPosCan_1darray(numUniqProIDCan*max_num_diag);
    int *begPosTar = begPosTar_1darray.array1D;
    int *endPosTar = endPosTar_1darray.array1D;
    int *begPosCan = begPosCan_1darray.array1D;
    int *endPosCan = endPosCan_1darray.array1D;

    Array1D <double> postPer_1darray(numUniqProIDCan*max_num_diag);
    double *postPer = postPer_1darray.array1D;  // percentage after post scan
    

    int cntDomain = 0;

#ifdef DEBUG_TIME
    start = clock();
#endif
    int maxSeqLength = max_element( proLengthCan, 0, numProCan-1);
    for( i = 0; i < maxSeqLength; i ++) 
    { 
        DotSegLengthScore[i] = GetDotSegScoreLength(i); 
        DotSegDistWeight[i] = GetDotSegDistWeight(i); 
        DiagDistWeight[i] = GetDiagDistWeight(i); 
    } /*added 2008-04-24*/

    for(i = 0 ; i < numUniqProIDCan; i ++)// re-canculate the number of fragment for each candiate protein
    {
        idx = BinarySearch_String(uniqProIDCan[i], proIDListCan, numProCan);
        if(idx != -1)
        {
            length_candidate = proLengthCan[idx];
        }
        else
        {
            fprintf(stderr,"idCan=%s, not found in proIDListCan\n", uniqProIDCan[i]);
            assert(idx >= 0);
        }

        double rawPer;
        rawPer = (subTotalProIDCan[i] / double (length_candidate + length_target)) * 2 * 100;

#ifdef DEBUG_PRINT_ID
        fprintf(stdout, "%-4d %s count = %3d, rawPer = %.2lf, cutoff_rawper = %.2lf\n", i, uniqProIDCan[i], subTotalProIDCan[i], rawPer, cutoff_rawper );
#endif
        if(rawPer > cutoff_rawper) 
        {
            int numPair = subTotalProIDCan[i];
            my_strcpy(idCan, uniqProIDCan[i], SIZE_ID);

            int *pair_matrix_1d = NULL;
            int **pair_matrix_2d = NULL;
            if(method == 0)
            {
                pair_matrix_1d = new int[length_target*length_candidate];
            }
            else if(method == 1)
            {
                pair_matrix_2d = Create2DArray(pair_matrix_2d, length_target, length_candidate);
            }


            Array1D <int> posTar_1darray(numPair);
            Array1D <int> posCan_1darray(numPair);
            Array1D <float> fragScore_1darray(numPair);
            int *posTar = posTar_1darray.array1D; // position for target, the main sequence
            int *posCan = posCan_1darray.array1D; // position for the candidate, the found homo candidate
            float *fragScore = fragScore_1darray.array1D; // fragment-fragment profile comparison score

            //int cntPair = 0;
            //for(j = 0 ; j < numAllPair; j ++)
            //{
            //    if(strcmp(idCan, allProIDCan[j]) == 0)
            //    {
            //        if(cntPair > numPair-1)
            //        {
            //            fprintf(stderr,"idCan = %s, subTotalProIDCan[i] = %d, numAllPair = %d, j = %d\n", idCan, numPair, numAllPair,j );
            //            assert(cntPair < numPair);
            //        }
            //        posTar[cntPair] = allPosTar[j];
            //        posCan[cntPair] = allPosCan[j];
            //        cntPair ++;
            //    }
            //}

            //if(cntPair != numPair)
            //{
            //    fprintf(stderr,"idCan = %s, cntPair = %d, numPair = %d, numAllPair = %d\n", idCan, cntPair, numPair, numAllPair);
            //    assert(cntPair == numPair);
            //}
            for(j = 0 ; j < numPair; j ++)
            {
                posTar[j] = allPosTar[indexUniqProInAllProID[i][j]]; 
                posCan[j] = allPosCan[indexUniqProInAllProID[i][j]]; 
                if(weightScheme == 1)
                {
                    fragScore[j] = allFragScore[indexUniqProInAllProID[i][j]];
                }
            }

            if(method == 0)/*{{{*/
            {
                BuildPairMatrix(pair_matrix_1d, posTar, posCan, numPair, length_target, length_candidate);
                int revised_numPair = 0 ;
                revised_numPair = Treat_pair_Percent_homology(pair_matrix_1d, length_target, length_candidate);
                //            fprintf(stdout,"id = %s, numPair = %d, revised_numPair = %d\n",idCan,  numPair, revised_numPair);
                //            fprintf(stdout,"\n");
                //            PrintPairMatrix(pair_matrix, length_target, length_candidate);
                postPer[cntDomain] = double (Integer(revised_numPair / double (length_candidate+length_target) * 2 * 100));
                my_strcpy(domainID[cntDomain], uniqProIDCan[i], SIZE_ID);
                cntDomain++;
            }/*}}}*/
            else if (method == 1)/*{{{*/
            {
#ifdef DEBUG_TIME
    start = clock();
#endif
                BuildPairMatrix(pair_matrix_2d, posTar, posCan, numPair, length_target, length_candidate);
#ifdef DEBUG_TIME
    finish = clock();
    duration = double(finish-start)  /double(CLOCKS_PER_SEC);
    fprintf(stdout,"BuildPairMatrix length_candidate = %d, length_target = %d, cost %lf seconds\n", length_candidate, length_target, duration);
#endif

#ifdef DEBUG
                fprintf(stdout,"before scanning, %s (%d) vs %s (%d) \n", idTar, length_target, idCan, length_candidate);
                fprintf(stdout,"\n");
                PrintPairMatrix(pair_matrix_2d, length_target, length_candidate);
                fprintf(stdout,"\n");
#endif 
                Array1D <double> score_1darray(max_num_diag);
                double *score = score_1darray.array1D;

                Array1D <int> posMaxDiag_1darray(max_num_diag);
                Array1D <int> corner_1darray(max_num_diag);
                int *posMaxDiag = posMaxDiag_1darray.array1D;
                int *corner = corner_1darray.array1D;

                Array1D <int> tmp_begPosTar_1darray(max_num_diag);
                Array1D <int> tmp_endPosTar_1darray(max_num_diag);
                Array1D <int> tmp_begPosCan_1darray(max_num_diag);
                Array1D <int> tmp_endPosCan_1darray(max_num_diag);
                int *tmp_begPosTar = tmp_begPosTar_1darray.array1D;
                int *tmp_endPosTar = tmp_endPosTar_1darray.array1D;
                int *tmp_begPosCan = tmp_begPosCan_1darray.array1D;
                int *tmp_endPosCan = tmp_endPosCan_1darray.array1D;

                int numDiag = 0 ;
#ifdef DEBUG_TIME
    start = clock();
#endif

#ifdef DEBUG_PRINT_MATRIX_SCORE
    fprintf(stdout,"Target %s (%d) -- Candidate %s (%d)\n", idTar, length_target, idCan, length_candidate);
#endif
                numDiag = ScanPairDotPlot(pair_matrix_2d, length_target, length_candidate, posMaxDiag, corner, score, tmp_begPosTar, tmp_endPosTar, tmp_begPosCan, tmp_endPosCan, numDiag);

#ifdef DEBUG_TIME
    finish = clock();
    duration = double(finish-start)  /double(CLOCKS_PER_SEC);
    fprintf(stdout,"ScanPairDotPlot, cost %lf seconds\n", duration);
#endif
                //            fprintf(stdout,"id = %s, numPair = %d, revised_numPair = %d\n",idCan,  numPair, revised_numPair);
                //            fprintf(stdout,"\n");
#ifdef DEBUG
                fprintf(stdout,"after scanning, %s (%d) vs %s (%d) \n", idTar, length_target, idCan, length_candidate);
                fprintf(stdout,"\n");
                PrintPairMatrix(pair_matrix_2d, length_target, length_candidate);
                fprintf(stdout,"\n");
#endif 
                for(int k = 0 ;  k < numDiag;  k++)
                {
                    my_strcpy(domainID[cntDomain], uniqProIDCan[i], SIZE_ID);
                    postPer[cntDomain] = score[k];
                    begPosTar[cntDomain] = tmp_begPosTar[k];
                    endPosTar[cntDomain] = tmp_endPosTar[k];
                    begPosCan[cntDomain] = tmp_begPosCan[k];
                    endPosCan[cntDomain] = tmp_endPosCan[k];
                    cntDomain++;
                }
            }/*}}}*/
            else if (method == 2)/*{{{*/
                /*method 2, do not build the pair_matrix, create the diagnol
                 * directly. build matrix is a waste of time when there are
                 * only a few dots*/
            {
                Array1D <double> score_1darray(max_num_diag);
                double *score = score_1darray.array1D;

                Array1D <int> posMaxDiag_1darray(max_num_diag);
                Array1D <int> corner_1darray(max_num_diag);
                int *posMaxDiag = posMaxDiag_1darray.array1D;
                int *corner = corner_1darray.array1D;

                Array1D <int> tmp_begPosTar_1darray(max_num_diag);
                Array1D <int> tmp_endPosTar_1darray(max_num_diag);
                Array1D <int> tmp_begPosCan_1darray(max_num_diag);
                Array1D <int> tmp_endPosCan_1darray(max_num_diag);
                int *tmp_begPosTar = tmp_begPosTar_1darray.array1D;
                int *tmp_endPosTar = tmp_endPosTar_1darray.array1D;
                int *tmp_begPosCan = tmp_begPosCan_1darray.array1D;
                int *tmp_endPosCan = tmp_endPosCan_1darray.array1D;

#ifdef DEBUG_PRINT_MATRIX_SCORE
    fprintf(stdout,"Target %s (%d) -- Candidate %s (%d)\n", idTar, length_target, idCan, length_candidate);
#endif
                int numDiag = 0 ;
                numDiag = ScanPairDotPlot_Method2(posTar, posCan, fragScore, numPair, length_target, length_candidate, posMaxDiag, corner, score, tmp_begPosTar, tmp_endPosTar, tmp_begPosCan, tmp_endPosCan, numDiag, weightScheme);

                for(int k = 0 ;  k < numDiag;  k++)
                {
                    my_strcpy(domainID[cntDomain], uniqProIDCan[i], SIZE_ID);
                    postPer[cntDomain] = score[k];
                    begPosTar[cntDomain] = tmp_begPosTar[k];
                    endPosTar[cntDomain] = tmp_endPosTar[k];
                    begPosCan[cntDomain] = tmp_begPosCan[k];
                    endPosCan[cntDomain] = tmp_endPosCan[k];
                    cntDomain++;
                }
            }/*}}}*/

            if(method == 0)
            {
                delete [] pair_matrix_1d;
            }
            if(method == 1)
            {
                Delete2DArray(pair_matrix_2d, length_target);
            }
        }
        else  // neglect candiates with rawPer less than 10
        {
            postPer[cntDomain] = 0.0;
            my_strcpy(domainID[cntDomain], uniqProIDCan[i], SIZE_ID);
            begPosTar[cntDomain] = 0;
            endPosTar[cntDomain] = 0;
            begPosCan[cntDomain] = 0;
            endPosCan[cntDomain] = 0;
            cntDomain++;
        }
    }
    int numDomain = cntDomain;

#ifdef DEBUG_TIME
    finish = clock();
    duration = double(finish-start)  /double(CLOCKS_PER_SEC);
    fprintf(stdout,"Dot plot scanning %s cost %lf seconds\n", fragPairFile, duration);
#endif
#ifdef DEBUG_TIME
    start = clock();
#endif

    //sort postPer descendingly
    Array1D <int> indexPostPer_1darray(numDomain);
    int *indexPostPer = indexPostPer_1darray.array1D;
    for(j = 0; j < numDomain; j ++) { indexPostPer[j] = j; }
    QuickSort_index(indexPostPer, postPer, 0, numDomain-1, DESCENDING);

    int numPostPer10  =0;

    for(i = 0 ; i < numDomain; i ++)
    {
        if(postPer[i]  >= cutoff_postper)
            numPostPer10 ++;
    }

    FILE *fpout = fopen(postFragPairFile, "w");
    fprintf(fpout,"%s= %d\n","NumDistribution10", numPostPer10);
    int NN = min(topN, numPostPer10);

    fprintf(fpout,"#%-3s %-8s %6s %-8s %6s %6s %6s %-8s %6s %6s %6s\n", "Num", "idCan", "score", "idTar", "lenTar", "begTar", "endTar", "idCan", "lenCan", "begCan","endCan");
    for(i = 0 ; i < NN; i ++) // 2007-07-12, to keep  the format compatible with Res_$Id.txt(tuping's output), output the number 
    {
        idx = BinarySearch_String(domainID[indexPostPer[i]], proIDListCan, numProCan);
        if(idx != -1)
        {
            length_candidate = proLengthCan[idx];
        }
        else
        {
            fprintf(stderr,"idCan=%s, not found in proIDListCan\n", domainID[indexPostPer[i]]);
            assert(idx >= 0);
        }
        fprintf(fpout,"%4d %-8s %6.3lf %-8s %6d %6d %6d %-8s %6d %6d %6d\n", i, domainID[indexPostPer[i]], postPer[indexPostPer[i]], 
                idTar, length_target, begPosTar[indexPostPer[i]], endPosTar[indexPostPer[i]],
                domainID[indexPostPer[i]], length_candidate, begPosCan[indexPostPer[i]], endPosCan[indexPostPer[i]]
                );
    }
    fclose(fpout);

    for(i = 0; i < numUniqProIDCan; i ++)
    {
        delete [] indexUniqProInAllProID[i];
    }

#ifdef DEBUG_TIME
    finish = clock();
    duration = double(finish-start) /double(CLOCKS_PER_SEC);
    fprintf(stdout,"Post processing %s cost %lf seconds\n", fragPairFile, duration);
#endif

    return NN;
}

void test()/*{{{*/
{
    int *pair_matrix = NULL;
    int length_target = 154; //1AFA1
    int length_candidate = 63; //1CI6A

    Array1D <int> pair_matrix_1darray(length_target*length_candidate);
    pair_matrix = pair_matrix_1darray.array1D;

    char fragPairFile[MAX_PATH+1] = "t1.txt";

    Array1D <int> posTar_1darray(LONGEST_SHAPE);
    Array1D <int> posCan_1darray(LONGEST_SHAPE);
    int *posTar = posTar_1darray.array1D; // position for target, the main sequence
    int *posCan = posCan_1darray.array1D; // position for the candidate, the found homo candidate

    int numPair = ReadInPair(fragPairFile, posTar, posCan);
    BuildPairMatrix(pair_matrix, posTar, posCan, numPair, length_target, length_candidate);

    PrintPairMatrix(pair_matrix, length_target, length_candidate);

    int revised_numPair = 0 ;
    revised_numPair = Treat_pair_Percent_homology(pair_matrix, length_target, length_candidate);
    fprintf(stdout,"numPair = %d, revised_numPair = %d\n", numPair, revised_numPair);
    fprintf(stdout,"\n");
    PrintPairMatrix(pair_matrix, length_target, length_candidate);


}/*}}}*/
/*}}}*/

int main(int argc, char** argv)/*{{{*/
{
#ifdef DEBUG_TIME_MAIN
    clock_t start, finish;
    double duration;
#endif
    bool isNonOptionArg = false;
    if(argc < 2) 
    {
        PrintHelp();
        return 0;
    }
    int i,j;
    int method = 2;
    int weightScheme = 0; /*2008-01-15, Nanjiang*/
    char outpath[MAX_PATH+1] = "./";
    char extension[MAX_PATH+1] = "fragpost";

    char seqlenFileCan[MAX_PATH+1] = "/misc/casiodata2/pdbaa.seqlen";
    char seqlenFileTar[MAX_PATH+1] = "";
	bool isTestSeqLenFileSet = false;

    char listfile[MAX_PATH+1] = "";
    double value = 0.0;
    const char control_option[] = "aqs"; //options which control the program, and does not take parameters
    bool isAll = false;
    bool isQuiet = false;
    bool isSingle = false;

    set <string> filenamelist_set;
    set <string> ::iterator iss ;

    i = 1;
    while(i < argc)/*{{{*/
    {
        if(argv[i][0] == '-' && !isNonOptionArg) //options
        {
            isNonOptionArg = false;
            if(IsInCharSet(argv[i][1], control_option))//if argv[i][1] is in control_option, it might be used as -aqs
            {
                for(j = 1 ; j < int (strlen(argv[i])); j++)
                {
                    switch (argv[i][j])
                    {
                        case 'a': isAll = true; break;
                        case 'q': isQuiet = true; break;
                        case 's': isSingle = true; break;
                        default : fprintf(stderr,"Invalid option, non-control option '%c' can be used together with contorl-option\n", argv[i][j]); return -1;
                    }
                }
                i ++;
            }
            else if(strcmp(argv[i],"-h") == 0 ||strcmp(argv[i],"--help")==0 )
            {
                PrintHelp(); 
                return 0;
            }
            else if(strcmp(argv[i],"--changelog") == 0)
            {
                fprintf(stdout,"%s\n", changeLog);
                return 0;
            }
            else if(strcmp(argv[i],"-H") == 0 )
            {
                PrintVerboseHelp();
                return 0;
            }
            else if( (strcmp(argv[i],"--nm") == 0) )  
            {
                isFindingMultiDomain = false;
                i++;
            }
            else if( (strcmp(argv[i],"--rb") == 0) )  
            {
                isReadBinaryFragfile = true;
                i++;
            }
            else if( (strcmp(argv[i],"--ignore-non-diag") == 0) )  
            {
                isIgnoreNonDiag = true;
                i++;
            }
            else if( (strcmp(argv[i],"--outpath") == 0) )  
            {
                if( ( i = option_parser_filename(argc, argv, i, outpath)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--ext") == 0) )  
            {
                if( ( i = option_parser_filename(argc, argv, i, extension)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"-l") == 0) )  
            {
                if( ( i = option_parser_filename(argc, argv, i, listfile)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--len-file") == 0) )  
            {
                if( ( i = option_parser_filename(argc, argv, i, seqlenFileCan)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--test-len-file") == 0) )  
            {
                if( ( i = option_parser_filename(argc, argv, i, seqlenFileTar)) == -1)
                    return -1;
				isTestSeqLenFileSet = true;
            }
            else if (strcmp(argv[i], "--begin") == 0)
            {
                if( ( i = option_parser_numeric(argc, argv, i, beginID, true, 0, 1000000)) == -1)
                    return -1;
            }
            else if (strcmp(argv[i], "--end") == 0)
            {
                if( ( i = option_parser_numeric(argc, argv, i, endID, true, 0, 1000000)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--method") == 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, method, true, 0, 2)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--weight") == 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, weightScheme, true, 0, 4)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--topn") == 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, topN, true, 0, 10000)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--cs") == 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, cutoff_rawper, true, 0.0, 1000.0)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--cfragscore") == 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, cutoff_fragscore, true, float(0.0), float(10000000.0))) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--cps") == 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, cutoff_postper, true, 0.0, 1000.0)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--value") == 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, value, true, 0.0, 15.0)) == -1)
                    return -1;
            }
            else if (strcmp(argv[i], "--") == 0)//next item is non option argument
            {
                isNonOptionArg = true;
                i ++;
                continue;
            }
            else
            {
                fprintf(stderr,"Error! Invalid argument '%s'\n", argv[i]);
                return -1;
            }
        }
        else //non-option argument
        {
            filenamelist_set.insert(argv[i]);
            i ++;
        }
    }/*}}}*/

    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;

    if(strcmp(listfile,"") == 0 && filenamelist_set.size() == 0)
    {
        fprintf(stderr, "Error! Neither listfile nor fragPairFile set in the argument list \n");
        return -1;
    }
    else if(strcmp(listfile, "") != 0)
    {
        FILE *fplist;
        fplist = fopen(listfile, "r");
        checkfilestream(fplist, listfile, "r");
        while((linesize = fgetline(fplist, line ,maxline )) != EOF)
        {
            if(linesize <= 0 ) continue;
            filenamelist_set.insert(line);
        }
        fclose(fplist);
    }

    VerifyFolder(outpath);

    //get sequence length file
    int max_num_id_can = fgetlinecnt(seqlenFileCan);
    Array2D <char> proIDListCan_2darray(max_num_id_can, SIZE_ID+1);
    char **proIDListCan = proIDListCan_2darray.array2D;
    Array1D <int> proLengthCan_1darray(max_num_id_can);
    int *proLengthCan = proLengthCan_1darray.array1D;
    int numProCan = ReadInSeqLenth(seqlenFileCan, proIDListCan, proLengthCan);


	if (!isTestSeqLenFileSet)
	{
		my_strcpy(seqlenFileTar, seqlenFileCan, MAX_PATH); /*added 2009-06-25*/
	}
    int max_num_id_tar = fgetlinecnt(seqlenFileTar);
    Array2D <char> proIDListTar_2darray(max_num_id_tar, SIZE_ID+1);
    char **proIDListTar = proIDListTar_2darray.array2D;
    Array1D <int> proLengthTar_1darray(max_num_id_tar);
    int *proLengthTar = proLengthTar_1darray.array1D;
    int numProTar = ReadInSeqLenth(seqlenFileTar, proIDListTar, proLengthTar);

    char fragPairFile[MAX_PATH+1] = "";
    char postFragPairFile[MAX_PATH+1] = "";
    char rtname[MAX_PATH+1] = "";

    int cntfile = 0;
    for(iss = filenamelist_set.begin(); iss != filenamelist_set.end(); iss ++)
    {
        cntfile ++;
        if (cntfile <(beginID+1)|| cntfile > endID)
        {
            continue; /*run only the test idTar from begin to end, 2007-11-16, Nanjiang Shu*/
        }
        my_strcpy(fragPairFile, (*iss).c_str(), MAX_PATH);
        rootname(fragPairFile, rtname);
        sprintf(postFragPairFile, "%s/%s.%s", outpath, rtname, extension);

#ifdef DEBUG_TIME_MAIN
        start = clock();
#endif

        PostScanFragSearch(fragPairFile, postFragPairFile, proIDListCan, proLengthCan , numProCan, proIDListTar, proLengthTar, numProTar, method, weightScheme);
#ifdef DEBUG_TIME_MAIN
        finish = clock();
        duration = double(finish-start)  /double(CLOCKS_PER_SEC);
        fprintf(stdout,"PostFragPairFile cost %lf seconds\n", duration);
#endif

        fprintf(stdout,"%d \t %s output\n", cntfile, postFragPairFile);
    }

    return 0;
}
/*}}}*/
