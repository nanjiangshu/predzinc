#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include "array.h"
#include "mytemplate.h"
#include "mypro.h"
#include "myfunc.h"

using namespace std;
/*Valgrind checked
 * 2011-10-14 01:13:11 Friday Week 41  Valgrind, no error no memory leaking.
 *
 * */
#if defined(_Windows) || defined(__WINDOWS__) || \
    defined(__WIN32__) || defined(WIN32) || \
defined(__WINNT__) || defined(__NT__)
#   ifndef WINDOWS
#       define WINDOWS
#   endif
#endif

#ifndef INLINE
#  if __GNUC__
#    define INLINE extern inline
#  else
#    define INLINE inline
#  endif
#endif

#ifndef QIJ_NAME_FORMAT
#define QIJ_NAME_FORMAT
#define QIJ_FORMAT_TUPING 0
#define QIJ_FORMAT_NANJIANG 1
#endif

#ifndef MODM_NAME_FORMAT
#define MODM_NAME_FORMAT
#define MODM_FORMAT_TUPING 0
#define MODM_FORMAT_NANJIANG 1
#endif

#ifndef FRAGACC_NAME_FORMAT
#define FRAGACC_NAME_FORMAT
#define FRAGACC_FORMAT_TUPING 0
#define FRAGACC_FORMAT_NANJIANG 1
#endif

#ifndef FRAGFORMAT
#define FRAGFORMAT
#define FRAGFORMAT_TUPING 0
#define FRAGFORMAT_NANJIANG 1
#endif

#ifndef PROFILESCORETYPE
#define PROFILESCORETYPE
#define PROFILESCORETYPE_TUPING 0
#define PROFILESCORETYPE_NANJIANG 1
#endif

#define TRAINING_SET 0
#define TEST_SET 1

#undef SIZE_ID
#define SIZE_ID 100

#define INIT_FRAGSCORE -1234567
/*for debugging*/


#define DATATYPE_LOGPER float  /*using int instead of float save 10% computational time, however, with O3 option in g++, float is even faster than int*/
#define LOG_PER_SCALE 1 /*change this one to a larger value if use DATATYPE_LOGPER as int*/
#define PER_TEN_SCALE 1
#define FRAGSCORE_SCALE     (1.0/(PER_TEN_SCALE * LOG_PER_SCALE))
int profileScoreType = PROFILESCORETYPE_TUPING; /*default profileScoreType using the floating non scaled value 2009-06-15*/
int fragformat = FRAGFORMAT_NANJIANG; /*default fragformat is using 5 columns, 2009-06-23*/


int SCore_Sample = 100; /*topN to calculate the segment-segment score*/
int Save_Sample = 100; 
int Consens_Sample = 5; 
int NPer_Frag_Database = 45; /*ratio in percentage on FracAcc Matrix
see the formula profileMerged[j] = ratio*matFrag[i][j]+(1.0-ratio)*matMODM[i][j];*/
int Nsearch_beg = 0;
int beginID = 0;       /*in order to run the program simultaneously, setting the begin and end position of the idlist to run*/
int endID = 0x7FFFFFFF; /*by default, running all ids in the idListFile, 2007-11-16 */

bool isReadBinaryFile = false;
bool isWriteBinaryFile = false;
bool isPredictSecShape = false; /*by default, do not predict the secondary structures*/

int mergeSide = 0; /*mergeSide 0 on training set, and 1 on test set*/
int ratioScheme = 0;/*ratioScheme, 0 for fixed ratio and 1 for varied ratio*/
float speedRate = 0.0; /*set the speed up rate, if speedRate == 0, do not speed up. 2008-01-31, Nanjiang*/
int8 typeProfile=1; /*ProfileSADByte*/
int dsspMapMethod=0; /*this variable is actually not used for search_new*/

int dbtype=1; /*0 - each id has a file, 1 - all files dumped*/

#undef NUM_20_AA
#define NUM_20_AA 20
//int idtype = 0; [>default idtype = 0, meaning standardized 5 character chain identifier<]

const int blosum62[][NUM_BLOSUM] =/*{{{*/
{
    {  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4},
    { -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4},
    { -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4},
    { -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4},
    {  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4},
    { -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4},
    { -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
    {  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4},
    { -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4},
    { -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4},
    { -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4},
    { -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4},
    { -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4},
    { -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4},
    { -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4},
    {  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4},
    {  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4},
    { -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4},
    { -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4},
    {  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4},
    { -2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4},
    { -1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
    {  0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4},
    { -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1}
};
/*}}}*/
const int AAS_Code[] = /*{{{*/
{
    0         , //A--ALA
    20        , //B--Non
    13        , //C--CYS
    18        , //D--ASP
    16        , //E--GLU
    5         , //F--PHE
    10        , //G--GLY
    9         , //H--HIS
    3         , //I--ILE
    20        , //J--Non
    7         , //K--LYS
    2         , //L--LEU
    6         , //M--MET
    15        , //N--ASN
    20        , //O--Non
    4         , //P--PRO
    19        , //Q--GLN
    8         , //R--ARG
    11        , //S--SER
    12        , //T--THR
    20        , //U--Non
    1         , //V--VAL
    17        , //W--TRP
    20        , //X--Non
    14        , //Y--TYR
    20          //Z--Non
};
/*}}}*/
const int Shape_Code[] = /*{{{*/
{
    5        , //A--A
    8        , //B-- -
    8        , //C-- -
    8        , //D-- -
    8        , //E-- -
    8        , //F-- -
    7        , //G-- G
    8        , //H-- -
    8        , //I-- -
    8        , //J-- -
    4        , //K-- K
    8        , //L-- -
    8        , //M-- -
    8        , //N-- -
    8        , //O-- -
    8        , //P-- -
    8        , //Q-- -
    1        , //R-- R
    0        , //S-- S
    6        , //T-- T
    2        , //U-- U
    3        , //V-- V
    8        , //W-- -
    8        , //X-- -
    8        , //Y-- -
    8          //Z-- -
};/*}}}*/
const double  bkFreq[] =/*{{{*/  // ordered by AAAlphabet_Tuping AVLIPFMKRHGSTCYNEWDQ
{
    0.0788 ,  // A
    0.0673 ,  // V
    0.0965 ,  // L
    0.0590 ,  // I
    0.0483 ,  // P
    0.0396 ,  // F
    0.0238 ,  // M
    0.0593 ,  // K
    0.0540 ,  // R
    0.0229 ,  // H
    0.0696 ,  // G
    0.0683 ,  // S
    0.0540 ,  // T
    0.0150 ,  // C
    0.0303 ,  // Y
    0.0413 ,  // N
    0.0667 ,  // E
    0.0114 ,  // W
    0.0534 ,  // D
    0.0395    // Q
};
/*}}}*/
char usage[]="\n\
usage: search_new [options] --train trainIDListFile --test testIDListFile\n\
Note: when dbtype is set to 1, train-qij, train-modm and train-fragacc \n\
      are set as the filepath for dbname, see my_extractdb.py\n\
options:\n\
  -dbtype  0|1        Set the database type, (default  1)\n\
                      0 - each id has a file, 1 - all files are dumped in a single one.\n\
  --result   STR      Path for result files, (default: ./). FINAL OUTPUT\n\
  --rb                Read the binary format file matrix, *.Qijbin, *.modmbin, *.fragaccbin\n\
  --wb                Write the binary format file, *.fragbin and *.secpredbin file \n\
  --train-qij STR     Path for Qij files of the training set\n\
  --test-qij  STR     Path for Qij files of the test set\n\
  --qijformat 0|1     Format of qijfile name,  (default: 1)\n\
                      0 - Qijmatrix_$id.txt, 1 - $id.Qij\n\
  --train-modm STR    Path for modmatrix files\n\
  --test-modm  STR    Path for test modmatrix files\n\
  --modmformat 0|1    Format of modm files, (default: 1)\n\
                      0 - modmatrix_$id.txt, 1 - $id.modm\n\
  --train-fragacc STR Path for result acc files\n\
  --test-fragacc  STR Path for result acc files of test chains\n\
  --fragaccformat 0|1 Format of modm files, (default:1)\n\
                      0: frag_$id.txt, 1 - $id.fragacc\n\
  --matrix   STR      Substitution matrix file, if not set, use blosum62\n\
  --bkfile   STR      Background amino acid composition in PERCENTAGE, if not set, use bkFreq*100\n\
  --para     STR      File for control parameters\n\
  --NPer     INT      NPer_Frag_Database, default = 45\n\
  --mergeside    0|1  Merging side for fragMat and Qij, 0: on training set, 1: on test set, default = 0\n\
  --ratioscheme  0|1  How to set the ratio to merge the fragMat and Qij, (default: 0)\n\
                      0: ratio set by --Nper, 1: auto ratio depends on the profile score2\n\
  --begin    INT      Start number of ids in the idListFile to run database build\n\
  --end      INT      End number of ids in the idListFile to run database build\n\
                      if begin and end is not set, running all IDs in the\n\
                      idListFile, id[endID] will not be run that is, --begin\n\
                      = 10, --end = 12, will run the 11th and 12th items\n\
  -w|--fragsize  INT  Set the size of fragment, (default: 9)\n\
  -N|--topN      INT  Set the maximal number of high-scoring fragments to keep, (default: 100)\n\
  -s|--speed   FLOAT  Set the speed ratio, [0-1], default = 0, means do not speed up\n\
  --scoretype   0|1   Set the score type for profile-profile score. (default: 1).\n\
                      0: score will be scaled by int (score*15+5000)\n\
  --fragformat  0|1   Set the type for output frag file, (default: 1).    \n\
                      0: six columns, 1: five colummns\n\
  -h|--help           Print this help message and exit\n\
\n\
Created on 2007-06-06, updated 2011-10-13, Nanjiang Shu\n\
";

void PrintHelp()
{
    fprintf(stdout,"%s",usage);
}
void PrintVerboseHelp() { }

int ReadInIDList(const char* file, char ** idList, int longestlinelength=0)/*{{{*/
{
    if (longestlinelength <=0){
        int maxRawIDCnt;
        maxRawIDCnt = fgetlinecnt(file, longestlinelength, false);
        if (maxRawIDCnt < 0){
            return -1;
        }
    }
    FILE *fpin = fopen(file,"r");
    if (fpin != NULL){
        int linesize;
        int maxline = longestlinelength;
        Array1D <char> line_1darray(maxline+1);
        char *line = line_1darray.array1D;
        int cnt=0;
        while((linesize = fgetline(fpin, line ,maxline)) != EOF) //idList, main interation
        {
            if(linesize <= 0) { continue; }
            Array2D <char> tmpstrs_2darray(2, linesize+1);
            char ** tmpstrs = tmpstrs_2darray.array2D;
            int status_sscanf = 0;
            status_sscanf = sscanf(line, "%s %s", tmpstrs[0], tmpstrs[1]);
            if (status_sscanf == 1){
                strcpy(idList[cnt], tmpstrs[0]);
            }else if (status_sscanf == 2){
                strcpy(idList[cnt], tmpstrs[1]);
            }else {
                fprintf(stderr,"Wrong idlist line \"%s\". Ignore.\n",line);
                continue;
            }
            cnt ++;
        }
        fclose(fpin);
        return cnt;
    }
    else{
        fprintf(stderr,"Failed to open file \"%s\". Return -1.\n", file);
        return -1;
    }
}/*}}}*/

INLINE int GetChainCodeAAS(char aa)/*{{{*/
{
    int tint = 0;
    tint = aa - 'A';
    if (  (tint>=0) && (tint<=25)  )
    {
        return AAS_Code[tint];
    }
    else
    {
        return 20;
    }
}/*}}}*/
INLINE DATATYPE_LOGPER Prob_score(DATATYPE_LOGPER *Per_ten_one, DATATYPE_LOGPER *Log_per, DATATYPE_LOGPER* Log_Back_Comp, int* pMatTar, int *pMatCan)/*{{{*/
/*****************************************************************************
 *  calculate the prob_score which is sum(Ai*log(Bi/Pi)+Bi*log(Ai/Pi)), for i =  0 to 20
 *  this function is used to speed up the calculation
 *  the program must be compiled with -O3 options to achieve the fastest
 *  executable
 *  but actually with -O3 options, this procedure is not faster, do not used it
 *  it is probably because in each row, 5 arrays are used, they are not
 *  contiguous in memory
 *  2008-02-26, Nanjiang
 ****************************************************************************/
{
    return (
        (Per_ten_one[pMatTar[0]]*(Log_per[pMatCan[0]]-Log_Back_Comp[0]) + Per_ten_one[pMatCan[0]]*(Log_per[pMatTar[0]]-Log_Back_Comp[0]))+
        (Per_ten_one[pMatTar[1]]*(Log_per[pMatCan[1]]-Log_Back_Comp[1]) + Per_ten_one[pMatCan[1]]*(Log_per[pMatTar[1]]-Log_Back_Comp[1]))+
        (Per_ten_one[pMatTar[2]]*(Log_per[pMatCan[2]]-Log_Back_Comp[2]) + Per_ten_one[pMatCan[2]]*(Log_per[pMatTar[2]]-Log_Back_Comp[2]))+
        (Per_ten_one[pMatTar[3]]*(Log_per[pMatCan[3]]-Log_Back_Comp[3]) + Per_ten_one[pMatCan[3]]*(Log_per[pMatTar[3]]-Log_Back_Comp[3]))+
        (Per_ten_one[pMatTar[4]]*(Log_per[pMatCan[4]]-Log_Back_Comp[4]) + Per_ten_one[pMatCan[4]]*(Log_per[pMatTar[4]]-Log_Back_Comp[4]))+
        (Per_ten_one[pMatTar[5]]*(Log_per[pMatCan[5]]-Log_Back_Comp[5]) + Per_ten_one[pMatCan[5]]*(Log_per[pMatTar[5]]-Log_Back_Comp[5]))+
        (Per_ten_one[pMatTar[6]]*(Log_per[pMatCan[6]]-Log_Back_Comp[6]) + Per_ten_one[pMatCan[6]]*(Log_per[pMatTar[6]]-Log_Back_Comp[6]))+
        (Per_ten_one[pMatTar[7]]*(Log_per[pMatCan[7]]-Log_Back_Comp[7]) + Per_ten_one[pMatCan[7]]*(Log_per[pMatTar[7]]-Log_Back_Comp[7]))+
        (Per_ten_one[pMatTar[8]]*(Log_per[pMatCan[8]]-Log_Back_Comp[8]) + Per_ten_one[pMatCan[8]]*(Log_per[pMatTar[8]]-Log_Back_Comp[8]))+
        (Per_ten_one[pMatTar[9]]*(Log_per[pMatCan[9]]-Log_Back_Comp[9]) + Per_ten_one[pMatCan[9]]*(Log_per[pMatTar[9]]-Log_Back_Comp[9]))+
        (Per_ten_one[pMatTar[10]]*(Log_per[pMatCan[10]]-Log_Back_Comp[10]) + Per_ten_one[pMatCan[10]]*(Log_per[pMatTar[10]]-Log_Back_Comp[10]))+
        (Per_ten_one[pMatTar[11]]*(Log_per[pMatCan[11]]-Log_Back_Comp[11]) + Per_ten_one[pMatCan[11]]*(Log_per[pMatTar[11]]-Log_Back_Comp[11]))+
        (Per_ten_one[pMatTar[12]]*(Log_per[pMatCan[12]]-Log_Back_Comp[12]) + Per_ten_one[pMatCan[12]]*(Log_per[pMatTar[12]]-Log_Back_Comp[12]))+
        (Per_ten_one[pMatTar[13]]*(Log_per[pMatCan[13]]-Log_Back_Comp[13]) + Per_ten_one[pMatCan[13]]*(Log_per[pMatTar[13]]-Log_Back_Comp[13]))+
        (Per_ten_one[pMatTar[14]]*(Log_per[pMatCan[14]]-Log_Back_Comp[14]) + Per_ten_one[pMatCan[14]]*(Log_per[pMatTar[14]]-Log_Back_Comp[14]))+
        (Per_ten_one[pMatTar[15]]*(Log_per[pMatCan[15]]-Log_Back_Comp[15]) + Per_ten_one[pMatCan[15]]*(Log_per[pMatTar[15]]-Log_Back_Comp[15]))+
        (Per_ten_one[pMatTar[16]]*(Log_per[pMatCan[16]]-Log_Back_Comp[16]) + Per_ten_one[pMatCan[16]]*(Log_per[pMatTar[16]]-Log_Back_Comp[16]))+
        (Per_ten_one[pMatTar[17]]*(Log_per[pMatCan[17]]-Log_Back_Comp[17]) + Per_ten_one[pMatCan[17]]*(Log_per[pMatTar[17]]-Log_Back_Comp[17]))+
        (Per_ten_one[pMatTar[18]]*(Log_per[pMatCan[18]]-Log_Back_Comp[18]) + Per_ten_one[pMatCan[18]]*(Log_per[pMatTar[18]]-Log_Back_Comp[18]))+
        (Per_ten_one[pMatTar[19]]*(Log_per[pMatCan[19]]-Log_Back_Comp[19]) + Per_ten_one[pMatCan[19]]*(Log_per[pMatTar[19]]-Log_Back_Comp[19]))
        );
//    for (im=0;im<20; im++) [>this part cost 75% of the computational time<]/*{{{*/ /*OBSOLETE CODE 2008-02-26*/
//    {
//        // Removed code[>{{{<]
//        //if (   (matMerged[iTar][im]<back_comp[im]) && (matAllTrain[iksub][a0+im]<back_comp[im])  )
//        //{
//        //    param1 = Under_Avereage[im][matMerged[iTar][im]][matAllTrain[iksub][a0+im]];
//        //}
//        //else if (   (matMerged[iTar][im]>=back_comp[im]) && (matAllTrain[iksub][a0+im]>=back_comp[im])  )
//        //{
//        //    param1 = Per_ten_one[matMerged[iTar][im]]*(Log_per[matAllTrain[iksub][a0+im]]-Log_Back_Comp[im]) + Per_ten_one[matAllTrain[iksub][a0+im]]*(Log_per[matMerged[iTar][im]]-Log_Back_Comp[im]);
//        //}
//        //else
//        //{
//        //    param1 = Per_ten_one[matMerged[iTar][im]]*(Log_per[matAllTrain[iksub][a0+im]]-Log_Back_Comp[im]) + Per_ten_one[matAllTrain[iksub][a0+im]]*(Log_per[matMerged[iTar][im]]-Log_Back_Comp[im]);
//        //}[>}}}<]
//        [>calculate PICASO3 score, Nanjiang,2008-01-13<]

//        //param1 = Per_ten_one[matMerged[iTar][im]]*(Log_per[matAllTrain[iksub][a0+im]]-Log_Back_Comp[im]) + Per_ten_one[matAllTrain[iksub][a0+im]]*(Log_per[matMerged[iTar][im]]-Log_Back_Comp[im]);
//        value +=(Per_ten_one[pMatTar[im]]*(Log_per[pMatCan[im]]-Log_Back_Comp[im]) + Per_ten_one[pMatCan[im]]*(Log_per[pMatTar[im]]-Log_Back_Comp[im]));
//        //value +=(MA_MB[pMatTar[im]][pMatCan[im]] - MA_P[pMatTar[im]][im] + MA_MB[pMatCan[im]][pMatTar[im]] - MA_P[pMatCan[im]][im]);[>this is even slower<]

//#ifdef DEBUG_CHECK_FRAGSCORE
//        fprintf(stdout,"M1[%2d][%2d] (%3d)-- M2[%2d][%2d] (%3d) = %.2f\n", iTar, im,matMerged[iTar][im], jCan, im, matAllTrain[iksub][a0+im], param1);
//#endif [>DEBUG_CHECK_FRAGSCORE<]


//#ifdef DEBUG_MATRIX
//        //fprintf(stdout,"%3d, %3d : log2(pTar[%2d]) = %.2f, log2(pCan[%2d]) = %.2f, pTar[%2d] = %d, pCan [%2d] = %d, log2(P[%2d]) = %.2f, value=%.3f\n", iTar, jCan, im, Per_ten_one[matMerged[iTar][im]], im, Log_per[matAllTrain[iksub][a0+im]], im, matMerged[iTar][im], im, matAllTrain[iksub][a0+im], im, Log_Back_Comp[im], value);
//#endif
//    }/*}}}*/
}/*}}}*/
int WriteFragFile(const char* fragFile, int fragformat, int lengthTar, int Save_Sample, int fragSize, DATATYPE_LOGPER **Candidate_Score, int ** Candidate_Score_iksub, int **Candidate_Score_iklen, char **idListTrain )/*{{{*/
{
    int i,j;
    int  halfFragSize = fragSize/2;
    FILE *fpFrag = fopen(fragFile,"w");
    checkfilestream(fpFrag, fragFile, "w", true);
    for (i=0; i<lengthTar-fragSize; i++)/*i is the iterator of the target position*/
    {
        for (j=0; j<Save_Sample; j++)
        {
            if (  Candidate_Score[i][j] <= INIT_FRAGSCORE  )
            {
                continue;
            }
            if (fragformat == FRAGFORMAT_NANJIANG)
            {

                fprintf(fpFrag, "%d %d %s %d %d\n",i + halfFragSize,j, idListTrain[Candidate_Score_iksub[i][j]], 
                        Candidate_Score_iklen[i][j] + halfFragSize,Integer(Candidate_Score[i][j]));/*do not output the useless information: the index of chain in the chainIDList; 2008-01-25, Nanjiang*/
            }
            else /*(fragformat == FRAGFORMAT_TUPING)*/
            {
                fprintf(fpFrag, "%d %d %s %d %d %d\n",i ,j, idListTrain[Candidate_Score_iksub[i][j]], Candidate_Score_iksub[i][j],
                        Candidate_Score_iklen[i][j], Integer(Candidate_Score[i][j]));/*2009-06-23*/

            }
        }
    }
    fclose(fpFrag);
    return 0;
}/*}}}*/
int WriteBinaryFragFile(const char* fragFile, int fragformat, int lengthTar, int Save_Sample, int fragSize, DATATYPE_LOGPER **Candidate_Score, int ** Candidate_Score_iksub, int **Candidate_Score_iklen, char **idListTrain )/*{{{*/
{
    int halfFragSize = fragSize/2;
    set <string> idList_set; /*storing the actually used ids in the output frag file*/
    set <string> :: iterator iss;
    idList_set.clear();

    FragCanShort5 *fragCan5 = NULL;
    FragCanShort6 *fragCan6 = NULL;
    if (fragformat == FRAGFORMAT_NANJIANG)
    { fragCan5 = new FragCanShort5[lengthTar*Save_Sample+1]; }
    else 
    { fragCan6 = new FragCanShort6[lengthTar*Save_Sample+1]; }

    Array2D <char> idFragCanAll_2darray(lengthTar*Save_Sample, SIZE_ID+1); /*the maximum size of id is at least smaller than maxLineSize - 10*/
    Array1D <short> posTar_1darray(lengthTar*Save_Sample);
    Array1D <short> numCan_1darray(lengthTar*Save_Sample);
    posTar_1darray.Init(0);
    numCan_1darray.Init(0);
    char **idFragCanAll = idFragCanAll_2darray.array2D;
    short *posTar = posTar_1darray.array1D;
    short *numCan = numCan_1darray.array1D;

    int i,j;

    int cntCanTotal = 0;
    int cntTar = 0;
    unsigned int maxSizeID = 0;
    for (i=0; i<lengthTar-fragSize; i++)/*i is the iterator of the target position*/
    {
        if (fragformat == FRAGFORMAT_NANJIANG) 
        { posTar[i] = i+halfFragSize;}
        else 
        { posTar[i] = i;}
        numCan[i] = Save_Sample;

        for (j=0; j<Save_Sample; j++)
        {
            if (  Candidate_Score[i][j] <= INIT_FRAGSCORE  )
            {
                continue;
            }
            my_strcpy(idFragCanAll[cntCanTotal], idListTrain[Candidate_Score_iksub[i][j]], SIZE_ID);
            idList_set.insert(idListTrain[Candidate_Score_iksub[i][j]]);
            if (strlen(idFragCanAll[cntCanTotal]) > maxSizeID) 
            {
                maxSizeID = strlen(idFragCanAll[cntCanTotal]);
            }
            if (fragformat == FRAGFORMAT_NANJIANG)
            {
                fragCan5[cntCanTotal].posCan =  short (Candidate_Score_iklen[i][j] + halfFragSize);
                fragCan5[cntCanTotal].score =  float(Candidate_Score[i][j]);
            }
            else /*(fragformat == FRAGFORMAT_TUPING)*/
            {
                fragCan6[cntCanTotal].idxOuter = short (Candidate_Score_iksub[i][j] );
                fragCan6[cntCanTotal].posCan =  short (Candidate_Score_iklen[i][j] );
                fragCan6[cntCanTotal].score =  float(Candidate_Score[i][j]);
            }
            cntCanTotal ++;
        }
        cntTar ++;
    }

    int numID = idList_set.size();
    Array2D <char> idList_2darray(numID, maxSizeID+1);
    char **idList = idList_2darray.array2D;
    
    i = 0;
    for(iss = idList_set.begin(); iss != idList_set.end(); iss ++)
    {
        my_strcpy(idList[i], (*iss).c_str(), maxSizeID);
        i ++;
    }
    /*find the inner index for ids*/
    int idxInner = 0;
    for ( i = 0; i < cntCanTotal; i ++)
    {
        if((idxInner = BinarySearch_String(idFragCanAll[i], idList, numID)) == -1)
        {
            fprintf(stderr,"Error! id = %s not found in idList\n", idFragCanAll[i]);
            exit(1);
        }
        if (fragformat == FRAGFORMAT_NANJIANG)
        {
            fragCan5[i].idxInner  = short (idxInner);
        }
        if (fragformat == FRAGFORMAT_TUPING)
        {
            fragCan6[i].idxInner  = short (idxInner);
        }
    }
    if (fragformat == FRAGFORMAT_NANJIANG)
    {
        WriteBinaryFragShort5(fragFile, fragformat, idList, numID, maxSizeID, cntTar, posTar, numCan, cntCanTotal,fragCan5 );
    }
    else
    {
        WriteBinaryFragShort6(fragFile, fragformat, idList, numID, maxSizeID, cntTar, posTar, numCan, cntCanTotal, fragCan6);
    }


    if (fragformat == FRAGFORMAT_NANJIANG)
    { delete [] fragCan5; }
    else 
    { delete [] fragCan6; }
    return 0;
}/*}}}*/
int PredictSecondaryStructure(const char* idTar, char *aaSeqTar,  DATATYPE_LOGPER **Candidate_Score, int **Candidate_Score_iksub, int **Candidate_Score_iklen, int lengthTar, int fragSize, int *lengthListAllTrain, char ** idListTrain, char **secSeqAllTrain, char **shapeSeqAllTrain, int numChainTrain, const char * resultpath)/*{{{*/
/*****************************************************************************
 * predict the secondary structure from the pre-built segment candidate table
 * 2008-01-20, Nanjiang
 * idTar = ""; name of the target chain
 ****************************************************************************/
{
    int i,j ;
    int im = 0;
    int HSR = 0;
    int CCvalue = 0; 
    int halfFragSize = 0 ;
    int  HSRsum1;
    int  NSUM;
    int  SRUV[8];
    char CNameHSR[]                     = "HSRR";
    char CShape[]                       = "SRUVKATG-";
    int  NumDistribution10;
    int Pred_HSR[4];
    int Pred_SHAPE[10];
    int NPer[21];
    int  target_class;
    int NTG_COM;
    int GGroup[5];
    int Numberofgroup;
    int Weight_per;
    float Sumvalue;

    halfFragSize = fragSize/2;
    NTG_COM = 2;
    NTG_COM = Consens_Sample/15; /*??why*/

    Array1D <int> Candidate_Distribution_1darray(numChainTrain+1);
    Array1D <int> Candidate_Distribution10_1darray(numChainTrain+1);
    Array1D <int> Distribution10_chain_1darray(numChainTrain+1);
    int *Candidate_Distribution = Candidate_Distribution_1darray.array1D;  /*Candidate_Distribution is the number of valid high-scoring fragment in the table for each candidate chain, Nanjiang, 2008-01-13*/ 
    int *Candidate_Distribution10 = Candidate_Distribution10_1darray.array1D;/*Candidate_Distribution10 stores the percentages of candidate chains with numDot/length >=10% */
    int *Distribution10_chain = Distribution10_chain_1darray.array1D; /*Distribution10_chain is the index to chains with numDot /length > 10%, 2008-01-20, Nanjiang*/

    for (i=0; i<numChainTrain; i++) /*i is the iterator for chains, Nanjiang, 2008-01-13 */
    {
        Candidate_Distribution[i] = 0;
        Candidate_Distribution10[i] = 0;
        Distribution10_chain[i] = 0;
    }

    for (i=0; i<lengthTar-fragSize; i++)
    {
        for (j=0; j<SCore_Sample; j++)
        {
            if (  Candidate_Score[i][j] <= INIT_FRAGSCORE  ) /*ignore the items with Candidate_Score = INIT_FLOAT*/
            {
                continue;
            }
            Candidate_Distribution[Candidate_Score_iksub[i][j]]++; /*Candidate_Distribution is the number of valid high-scoring fragment in the table for each candidate chain, Nanjiang, 2008-01-13*/
        }
    }

    NumDistribution10 = 0;
    for (i=0; i<numChainTrain; i++) 
    {
        int percentFrag = Candidate_Distribution[i]*2*100/(lengthTar+lengthListAllTrain[i]); /*percentFrag is the normalized frequency of high-scoring fragments for each candidate chain in the fragment table, Nanjiang, 2008-01-13, HSR changed to percentFrag*/
#ifdef DEBUG
        //check which value is uninitialized
        fprintf(stdout,"percentFrag=%d, Candidate_Distribution[%d]=%d, Sumunit = %d, lengthTar=%d, lengthListAllTrain=%d\n", percentFrag, i, Candidate_Distribution[i], Sumunit, lengthTar, lengthListAllTrain[i]);
#endif
        if (  percentFrag >= 10  )
        {
            Candidate_Distribution10[NumDistribution10] = percentFrag;
            Distribution10_chain[NumDistribution10] = i;
            NumDistribution10++;
        }
        /*map the percentFrag to Candidate_Distribution10, Nanjiang,
         * 2008-01-13, ask tuping, why doing so*/
        if (  percentFrag >= 60  )
        {
            Candidate_Distribution[i] = 22;
        }
        else if (   (percentFrag<60) && (percentFrag>=50)  )
        {
            Candidate_Distribution[i] = 20;
        }
        else if (  (percentFrag<50) && (percentFrag>=40)  )
        {
            Candidate_Distribution[i] = 18;
        }
        else if (  (percentFrag<40) && (percentFrag>=30)  )
        {
            Candidate_Distribution[i] = 16;
        }
        else if (  (percentFrag<30) && (percentFrag>=20)  )
        {
            Candidate_Distribution[i] = 14;
        }
        else if (  (percentFrag<20) && (percentFrag>=15)  )
        {
            Candidate_Distribution[i] = 12;
        }
        else if (  (percentFrag<15) && (percentFrag>=10)  )
        {
            Candidate_Distribution[i] = 11;
        }
        else
        {
            Candidate_Distribution[i] = 10;
        }
    }

    Array2D <int> Pred_Result_2darray(lengthTar+1, 8);
    Array2D <int> Pred_Resulthsr_2darray(lengthTar+1,3);
    int **Pred_Result = Pred_Result_2darray.array2D;
    int **Pred_Resulthsr = Pred_Resulthsr_2darray.array2D;

    /*init*/
    for (i=0; i<lengthTar; i++)
    {
        for (j=0; j<8; j++) { Pred_Result[i][j] = 0; }
        for (j=0; j<3; j++) { Pred_Resulthsr[i][j] = 0; }
    }

    for (i=0; i<lengthTar-fragSize; i++)
    {
        //extract result
        for (im=0; im<fragSize; im++)//for fragSize positions of fragSize-long fragment
        {
            for (j=0; j<Consens_Sample; j++)
            {
                if (  Candidate_Score[i][j] <= MIN_FLOAT )
                {
                    continue;
                }
                //shape of candidate
                CCvalue = shapeSeqAllTrain[Candidate_Score_iksub[i][j]][Candidate_Score_iklen[i][j]+im];
                HSR = CCvalue - 'A';
                if  (  (HSR>=0) && (HSR<=25)  )
                {
                    HSR = Shape_Code[HSR];
                }
                else
                {
                    HSR = 8;
                }
                if (  HSR <= 7 )
                {
                    Pred_Result[i+im][HSR] = Pred_Result[i+im][HSR] + Candidate_Distribution[Candidate_Score_iksub[i][j]];
                }

                //HSR
                CCvalue = secSeqAllTrain[Candidate_Score_iksub[i][j]][Candidate_Score_iklen[i][j]+im];
                if  (  CCvalue == 'H'   )
                {
                    HSR = 0;
                }
                else if (  CCvalue == 'S'  )
                {
                    HSR = 1;
                }
                else
                {
                    HSR = 2;
                }
                if (  HSR <= 2 )
                {
                    Pred_Resulthsr[i+im][HSR] = Pred_Resulthsr[i+im][HSR] + Candidate_Distribution[Candidate_Score_iksub[i][j]];
                }
            }
        }
    }

    for (i=0; i<lengthTar; i++)
    {
        for (j=0; j<8; j++) { Pred_Result[i][j] = Pred_Result[i][j]/10; }
        for (j=0; j<3; j++) { Pred_Resulthsr[i][j] = Pred_Resulthsr[i][j]/10; }
    }

    for (i=0; i<8; i++) { Pred_SHAPE[i] = 0; }
    for (i=0; i<3; i++) { Pred_HSR[i] = 0; }

    /*write the secondary structure prediction result*//*{{{*/
    char secPredFile[MAX_PATH+1] = "";/*output file storing the secondary structure prediction*/
    sprintf(secPredFile,"%s/%s.secpred",resultpath,idTar);
    FILE *fpSecPred = fopen(secPredFile,"w");
    //shape and HSR
    //for (j=0; j<1800;j++)
    for(j = 0 ; j < lengthTar; j ++)
    {
        NSUM = 0;
        for (im=0; im<8;im++)
        { 
            NSUM = NSUM + Pred_Result[j][im];
            NPer[im] = Pred_Result[j][im];
        }
        target_class = NSUM/Consens_Sample;
        if (  target_class <= halfFragSize ) /*halfFragSize = 4, if the Pred_Result are all zero, NSUM = 0, target_class = 0, this will neglect non predicted positions, 2008-01-20, Nanjiang, ???*/
        { 
            continue;
        }
        NTG_COM = NSUM/15;
        GGroup[0] = NPer[0] + NPer[1] + NPer[2] + NPer[3];
        GGroup[1] = NPer[4] + NPer[5];
        GGroup[2] = NPer[6] + NPer[7] + NTG_COM;
        if (  (GGroup[2]>=GGroup[0]) && (GGroup[2]>=GGroup[1])   )//GT
        {
            if (  NPer[7] >= (NPer[6]-NTG_COM)  )
            {
                Numberofgroup = 7;
            }
            else
            {
                Numberofgroup = 6;
            }
            Weight_per = (GGroup[2]-NTG_COM)*100/NSUM;
        }
        else if (  (GGroup[0]>=GGroup[1]) && (GGroup[0]>=GGroup[2])   )//SRUV
        {
            Numberofgroup = 0;
            HSR = NPer[0];
            for (im=0; im<4; im++)
            {
                SRUV[im] = NPer[im];
            }
            SRUV[2] = SRUV[2] + NTG_COM;//U
            SRUV[3] = SRUV[3] + NTG_COM;//V
            for (im=1; im<4; im++)
            { 
                if (  SRUV[im] >= HSR  )
                {
                    HSR = SRUV[im];
                    Numberofgroup = im;
                }
            }
            Weight_per = GGroup[0]*100/NSUM;
        }
        else// (  (GGroup[1]>GGroup[0]) && (GGroup[1]>GGroup[2])   ), K A
        {
            if (  NPer[4] >= (NPer[5]-NTG_COM)  )
            {
                Numberofgroup = 4;
            }
            else
            {
                Numberofgroup = 5;
            }
            Weight_per = GGroup[1]*100/NSUM;
        }
        Pred_SHAPE[Numberofgroup]++;



        //HSR
        NPer[10] = Pred_Resulthsr[j][0];
        NPer[11] = Pred_Resulthsr[j][1];
        NPer[12] = Pred_Resulthsr[j][2];
        NTG_COM = (NPer[10]+NPer[11]+NPer[12])/15;
        NPer[11] = NPer[11] + NTG_COM;
        HSRsum1 = NPer[10] + NPer[11] + NPer[12];
        if  (   (NPer[11]>=NPer[10]) && (NPer[11]>=NPer[12])   )
        {
            HSR = 1;
            Sumvalue = float ( NPer[11]*100.0/HSRsum1 );
        }
        else if (   (NPer[10]>=NPer[11]) && (NPer[10]>=NPer[12])   )
        {
            HSR = 0;
            Sumvalue = float ( NPer[10]*100.0/HSRsum1 );
        }
        else
        {
            HSR = 2;
            Sumvalue = float ( NPer[12]*100.0/HSRsum1 );
        }
        Pred_HSR[HSR]++;

        fprintf(fpSecPred, "%4d %c  ",j,aaSeqTar[j]);
        fprintf(fpSecPred, " %3d  %3d  %3d  ",NPer[10],NPer[11],NPer[12]);
        for (im=0; im<8; im++)
        {
            fprintf(fpSecPred, "%3d ",NPer[im]);
        }
        fprintf(fpSecPred, "   %c %3d   %c %5.0f\n",CShape[Numberofgroup],Weight_per, CNameHSR[HSR],Sumvalue);
    }

    for (j=0; j<NumDistribution10;j++) /*order the Candidate_Distribution10 by descending order*//*{{{*/
    {
        for (im=j+1; im<NumDistribution10;im++)
        {
            if (  Candidate_Distribution10[im] > Candidate_Distribution10[j] )
            {
                HSR = Candidate_Distribution10[j];
                Candidate_Distribution10[j] = Candidate_Distribution10[im];
                Candidate_Distribution10[im] = HSR;

                HSR = Distribution10_chain[j];
                Distribution10_chain[j] = Distribution10_chain[im];
                Distribution10_chain[im] = HSR;
            }
        }
    }/*}}}*/

    int topN = min(10,NumDistribution10); /*topN candidate chains to output*/
    fprintf(fpSecPred,"NumDistribution10=%5d\n",NumDistribution10);
    for (j=0; j<topN;j++)
    {
        fprintf(fpSecPred,"%2d  %5s  %3d\n",j,idListTrain[Distribution10_chain[j]],Candidate_Distribution10[j]);
    }
    fclose(fpSecPred);/*}}}*/

    /*write the prediction list file*//*{{{*/
    char predictPercentListFile[MAX_PATH+1] = "";
    sprintf(predictPercentListFile,"%s/Test_resultper.txt",resultpath);
    FILE *fpPredList = fopen(predictPercentListFile,"a");
    fprintf(fpPredList,"%5s  ",idTar );
    for (im=0; im<3; im++)
    {
        fprintf(fpPredList,"%4d ",	Pred_HSR[im]);
    }
    fprintf(fpPredList,"     ");
    for (im=0; im<8; im++)
    {
        fprintf(fpPredList,"%4d ",	Pred_SHAPE[im]);
    }
    fprintf(fpPredList,"\n");
    fclose(fpPredList);/*}}}*/

    return 0;
}/*}}}*/

template <class T> int GetBinaryMODM_nofseek(BYTE *data, char *alphabet, int &length, T *profile, double *parameter, int8 &typeProfile)/*{{{*/
//Read in binary modm by given the FILE handler
//created 2011-10-13
{
#ifdef DEBUG_TIME_GETBINARYMODM_NOFSEEK
    clock_t start, finish;
    double duration;
    start = clock();
#endif

    int sizeAlphabet = 0;
    //long numByteReadIn=0;
    BYTE* pData = data;
    memcpy(&typeProfile, pData, sizeof(int8));
    pData += 1;
    memcpy(&sizeAlphabet, pData, sizeof(int));
    pData += sizeof(int);
    memcpy(&length, pData, sizeof(int));
    pData += sizeof(int);
    memcpy(alphabet, pData, sizeof(char)*sizeAlphabet);
    pData += sizeof(char)*sizeAlphabet;
    memcpy(parameter, pData, sizeof(double)*8);
    pData += sizeof(double)*8;
    memcpy(profile, pData, sizeof(T)*length);
#ifdef DEBUG_TIME_GETBINARYMODM_NOFSEEK
    finish = clock();
    duration = double(finish-start)  /double(CLOCKS_PER_SEC);
    fprintf(stdout,"GetBinaryMODM_nofseek %lf\n", duration);
#endif
    return 0;
} /*}}}*/

int ReadInDatabase_dumpedfile_nofseek(const char *id, map<string,dbindex>&dbindexfragacc, map<string,dbindex>&dbindexqij, map<string,dbindex>&dbindexmodm, vector <BYTE*> &dataList_fragacc, vector<BYTE*>&dataList_qij, vector<BYTE*>&dataList_modm, int **matFrag, int **matQij, int **matMODM, int **matMerged, float *matFragScore1, float *matFragScore2,  float *matQijScore1, float *matQijScore2, float *matMODMScore1, float *matMODMScore2, char *aaSeq, int *digitAASeq, char *shapeSeq, char *dsspSecSeq, int type_dataset, bool isReadBinaryFile, int NPer_Frag_Database, int ratioScheme, int mergeSide, int dsspMapMethod, int8 typeProfile, set <string> &testIDList_DEBUG_MERGESIDE)/*{{{*/
/*****************************************************************************
 * read in the dataset for a chain
 * including FragAcc matrix, Qij matrix and modm matrix
 * the type_dataset and isReadBinaryFile should be given
 * return the length of the chain if successful
 * return -1 for failure
 ****************************************************************************/
{
    int i,j;
    int lengthSeqMODM = 0;
    int lengthSeqQij = 0;
    int lengthFragacc = 0;
    Array1D <int> sumProfileFrag_1darray(LONGEST_SEQ);
    int *sumProfileFrag = sumProfileFrag_1darray.array1D; /*sum of the profile at each position for fragacc matrix*/
#ifdef NEGLECT_READ_DATABASE
    if (0){
#endif
    if (dbindexfragacc.find(id) == dbindexfragacc.end()){
        fprintf(stderr,"Can not find id %s in db %s. Ignore this id.\n",id, "fragacc" );
        return -1;
    }
    if (dbindexqij.find(id) == dbindexqij.end()){
        fprintf(stderr,"Can not find id %s in db %s. Ignore this id.\n",id, "Qij");
        return -1;
    }
    if (dbindexmodm.find(id) == dbindexmodm.end()){
        fprintf(stderr,"Can not find id %s in db %s. Ignore this id.\n",id, "modm");
        return -1;
    }

    BYTE *pData= NULL;

#ifdef DEBUG_READ_BIN_FILE
    clock_t start0,finish0;
    double duration0;
    start0=clock();
#endif
    dbindex *pDBindex= NULL;
    if (isReadBinaryFile == true) { /*isReadBinaryFile == true*/ /*{{{*/
        Array1D <ProfileSADByte> profile_1darray(LONGEST_SEQ);
        ProfileSADByte *profile = profile_1darray.array1D;
        char alphabet[50+1] = "";
        int length = 0; /*number of profiles read in*/
        Array1D <double> parameter_1darray(8);
        parameter_1darray.Init(0.0);
        double *parameter = parameter_1darray.array1D;

        int aaSeqIndex = 0;

        /*2007-11-11, Nanjiang, changed to read in Qij first and then if fragacc exist, replace it with fragacc*/
        /*read in Qij matrix file, for those residue positions does not have fragacc profile, using the Qij profile instead*/
#ifdef NEGLECT_READ_QIJ
        if(0){
#endif
        pDBindex = &(dbindexqij[id]) ;
#ifdef DEBUG_DBINDEX
        fprintf(stdout,"id=%s Qij %2d %8ld %8ld, dataList[0]=%s\n", id, pDBindex->dbfileindex, pDBindex->offset, pDBindex->size, dataList_qij[0]);
#endif
        pData = dataList_qij[ pDBindex->dbfileindex]+pDBindex->offset ;
        if (GetBinaryMODM_nofseek(pData, alphabet, length, profile, parameter, typeProfile) < 0 ) {
            fprintf(stderr, "Reading qij matrix failed for %s. offset=%ld, size=%ld. __LINE__=%d\n", id, dbindexqij[id].offset, dbindexqij[id].size, __LINE__);
            return -1;
        }
        lengthSeqQij = length ;
#ifdef NEGLECT_READ_QIJ
        }
#endif
#ifdef NEGLECT_COPY_QIJ
        if(0){
#endif
        for(i = 0 ; i < length; i ++) {
            aaSeqIndex = profile[i].aaSeqIndex-1;
            if (aaSeqIndex < 0) {
                fprintf(stderr,"Error aaSeqIndex (=%d) < 0, qij %s : %d\n", aaSeqIndex, id, i);
                return -1;
            }
            for (j=0; j<20; j++) {matQij[aaSeqIndex][j] = int(profile[i].p[j]); }
            matQijScore1[aaSeqIndex] = profile[i].score1; /*change i to aaSeqIndex, 2008-02-18, Nanjiang*/
            matQijScore2[aaSeqIndex] = profile[i].score2;
        }
#ifdef NEGLECT_COPY_QIJ
        }
#endif
        /*read in frag acc matrix file, replace the Psi_blosum_Frag matrix with fragacc profiles */

#ifdef NEGLECT_READ_FRAGACC
        if(0){
#endif
        if ((ratioScheme == 0 && NPer_Frag_Database == 0)|| mergeSide != type_dataset )
            /*changed 2009-10-22, if mergeSide != type_dataset, also do not
             * need to read in FragMatFile, then replace it by qijfile*/
        {
            pDBindex = &(dbindexqij[id]) ;
#ifdef DEBUG_DBINDEX
            fprintf(stdout,"id=%s Qij %2d %8ld %8ld\n", id, pDBindex->dbfileindex, pDBindex->offset, pDBindex->size);
#endif
            pData = dataList_qij[ pDBindex->dbfileindex]+pDBindex->offset ;
            if (GetBinaryMODM_nofseek(pData, alphabet, length, profile, parameter, typeProfile ) < 0) {
                fprintf(stderr, "Reading qij matrix failed for %s. offset=%ld, size=%ld. __LINE__=%d\n", id, dbindexqij[id].offset, dbindexqij[id].size, __LINE__);
                return -1;
            }
        } else {
            pDBindex = &(dbindexfragacc[id]) ;
#ifdef DEBUG_DBINDEX
            fprintf(stdout,"id=%s fragacc %2d %8ld %8ld\n", id, pDBindex->dbfileindex, pDBindex->offset, pDBindex->size);
#endif
            pData=dataList_fragacc[pDBindex->dbfileindex] +  pDBindex->offset;
            if (GetBinaryMODM_nofseek(pData, alphabet, length, profile, parameter, typeProfile ) < 0 ) {
                fprintf(stderr, "Reading fragacc matrix failed for %s. offset=%ld, size=%ld\n", id, dbindexfragacc[id].offset, dbindexfragacc[id].size);
                return -1;
            }
        }
        lengthFragacc = length;
#ifdef NEGLECT_READ_FRAGACC
        }
#endif

#ifdef NEGLECT_COPY_FRAGACC
       if(0){
#endif
        for(i=0; i<lengthSeqQij; i++) { sumProfileFrag[i] = 0; } /*init sumProfileFrag*/
        for(i = 0 ; i < lengthFragacc; i ++) {
            aaSeqIndex = profile[i].aaSeqIndex-1;
            if (aaSeqIndex < 0) {
                fprintf(stderr,"Error aaSeqIndex[%d] =%d < 0, fragacc for id %s\n", i, aaSeqIndex, id);
                assert(aaSeqIndex >= 0);
            }
            for (j=0; j<20; j++) {    
                matFrag[aaSeqIndex][j] =  int (profile[i].p[j]);
                sumProfileFrag[aaSeqIndex] +=  matFrag[aaSeqIndex][j];
            }
            matFragScore1[aaSeqIndex] = profile[i].score1; /*bug fixed, i changed to aaSeqIndex, 2008-02-18, Nanjiang*/
            matFragScore2[aaSeqIndex] = profile[i].score2;
        }
#ifdef NEGLECT_COPY_FRAGACC
       }
#endif

#ifdef NEGLECT_READ_MODM
        if(0){
#endif
        /*read in normal MODM matrix file*/
        pDBindex = &(dbindexmodm[id]) ;
#ifdef DEBUG_DBINDEX
        fprintf(stdout,"id=%s modm %2d %8ld %8ld\n", id, pDBindex->dbfileindex, pDBindex->offset, pDBindex->size);
#endif
        pData=dataList_modm[ pDBindex->dbfileindex] +  pDBindex->offset;
        if (GetBinaryMODM_nofseek(pData, alphabet, length, profile, parameter, typeProfile) < 0 ) {
            fprintf(stderr, "Reading modm matrix failed for %s. offset=%ld, size=%ld\n", id, dbindexmodm[id].offset, dbindexmodm[id].size);
            return -1;
        }

        lengthSeqMODM = length;
#ifdef NEGLECT_READ_MODM
        }
#endif

#ifdef NEGLECT_COPY_MODM
        if(0){
#endif
        for(i = 0 ; i < lengthSeqMODM; i ++) {
            aaSeqIndex = profile[i].aaSeqIndex-1;
            if (aaSeqIndex < 0) {
                fprintf(stderr,"Error aaSeqIndex (=%d) < 0, modm %s : %d\n", aaSeqIndex, id, i);
                assert(aaSeqIndex >= 0);
            }
            for (j=0; j<20; j++) {matMODM[aaSeqIndex][j] = int(profile[i].p[j]); }
            matMODMScore1[aaSeqIndex] = profile[i].score1;
            matMODMScore2[aaSeqIndex] = profile[i].score2;

            aaSeq[aaSeqIndex]   = profile[i].aa;   /*amino acids*/
            shapeSeq[aaSeqIndex] = profile[i].shape; /*shape strings*/
            dsspSecSeq[aaSeqIndex] = profile[i].dsspSec; /*dssp secondary structure*/
        }
#ifdef NEGLECT_COPY_MODM
        }
#endif
    }/*}}}*/
#if DEBUG_READ_BIN_FILE
    finish0=clock();
    duration0 = double(finish0-start0)  /double(CLOCKS_PER_SEC);
    fprintf(stdout,"Read bin data for for id %s costs %lf seconds\n", id, duration0);
#endif

#ifdef NEGLECT_READ_DATABASE
    }
#endif

#ifdef DEBUG_DATA_HANDEL
    clock_t start2, finish2;
    double duration2;
    start2=clock();
#endif

#ifdef NEGLECT_DATAHANDEL
    if(0){
#endif

    if (lengthSeqMODM != lengthSeqQij) {
        fprintf(stderr, "Warning! lengthSeqMODM(%d) not equal to lengthSeqQij(%d)\n", lengthSeqMODM , lengthSeqQij);
        fprintf(stderr, "lengthSeqMODM (%s) = %d\n", id , lengthSeqMODM);
        fprintf(stderr, "lengthSeqQij (%s) = %d\n", id, lengthSeqQij);
        /*2009-06-10, this will not affect the result, so no need for
         * assertion, just give the warning signal*/
    }
    aaSeq[lengthSeqMODM] = '\0';
    shapeSeq[lengthSeqMODM] = '\0';
    dsspSecSeq[lengthSeqMODM] = '\0';
    for(i = 0; i < lengthSeqMODM; i ++) { /*map 8-state dsspsec to 3-state, H,G-->Helix, E-->strand, others-->Random coil*/ 
        dsspSecSeq[i] = DSSP8to3(dsspSecSeq[i]); 
    }
    TreatAllZeroFij(matMODM, lengthSeqMODM, matMODMScore1, aaSeq, AAAlphabet_Tuping);
#ifdef DEBUG_READ_MODM
    fprintf(stdout,"modm for id %s\n", id);
    for(i = 0; i < lengthSeqMODM; i ++) {
        WriteMODMProfile(i+1, aaSeq[i], matMODM[i], 0,0, AAAlphabet_Tuping, stdout);
    }
    fprintf(stdout,"\n");
#endif

    int daa = 0;
    for(i = 0 ; i < lengthSeqMODM; i++) {
        /*get digit amino acid 0-19*/
        daa = aaSeq[i] - 'A';
        if (  (daa>=0) && (daa<=25)  ) { daa = AAS_Code[daa]; }
        else { daa = 20; }
        digitAASeq[i] = daa;
    }
    if(mergeSide == type_dataset) {/*{{{*/
        /*profile combination of Frag- and psi-blast*/
        for(i = 0; i < lengthSeqMODM; i ++) {
            double ratio = 0.0;
            double profileMerged[NUM_20_AA+1];
            if ( sumProfileFrag[i] >= 10  ) /*if the frag profile is valid, isAllZeroProfile = false*/
            {
                if(ratioScheme == 0){/*ratioScheme = 0, merge fragacc and modm by specified ratio, NPer_Frag_Database*/
                    ratio = NPer_Frag_Database/100.0;
                } else if(ratioScheme == 1) {/*ratioScheme = 1, merge fragacc and modm according to the score2, 2008-01-18, Nanjiang*/
                    ratio = GetMergeRatio1(i, matFragScore2, matMODMScore2, lengthSeqMODM);  /*bug fixed, matFragScore*/
                } else if (ratioScheme == 2) {
                    ratio = GetMergeRatio2(i, matFragScore2, matMODMScore2, lengthSeqMODM); 
                } else if (ratioScheme == 3) {
                    ratio = GetMergeRatio3(i, matFragScore2, matMODMScore2, lengthSeqMODM); 
                }

#ifdef DEBUG_MERGESIDE
                if (testIDList_DEBUG_MERGESIDE.find(id) !=  testIDList_DEBUG_MERGESIDE.end() )  {

                    fprintf(stdout,"id=%s, aa = %c, No = %d, mergeSide=%d, ratioScheme=%d, type_dataset=%d, ratio = %.3lf\n", id, aaSeq[i],i+1, mergeSide, ratioScheme, type_dataset, ratio);
                }
#endif
                double sum = 0.0;
                for (j=0; j<20; j++) {
                    profileMerged[j] = ratio*matFrag[i][j]+(1.0-ratio)*matMODM[i][j];
                    sum += profileMerged[j];
                }
                for (j=0; j<20; j++) {
                    matMerged[i][j] = Integer(profileMerged[j]*100.0/sum);/*normalize the sum of the profile to 100, Nanjiang, 2008-01-15*/
                }
#ifdef DEBUG_MERGESIDE /*test the mergeside and read in matrix*/
                if (testIDList_DEBUG_MERGESIDE.find(id) !=  testIDList_DEBUG_MERGESIDE.end() )  {
                    fprintf(stdout,"Profile for chain %s, three lines, matFrag, matMODM, matMerged\n", id);
                    WriteMODMProfile(i+1, aaSeq[i], matFrag[i], matFragScore1[i],matFragScore2[i], AAAlphabet_Tuping, stdout);
                    WriteMODMProfile(i+1, aaSeq[i], matMODM[i], matMODMScore1[i],matMODMScore2[i], AAAlphabet_Tuping, stdout);
                    //                WriteMODMProfile(i+1, aaSeq[i], profileMerged, 0.11,0.11, AAAlphabet_Tuping, stdout);
                    WriteMODMProfile(i+1, aaSeq[i], matMerged[i], 0,0, AAAlphabet_Tuping, stdout);
                }
#endif
            } else {/*if all components are zero in frag profiles, use the Qij profiles */
                for (j=0; j<20; j++) { matMerged[i][j] = matQij[i][j]; }
            }
        }/*}}}*/
    } else {/*do not merge the fragacc and qij on the training set, take the qij profile directly, 2008-01-18, Nanjiang*//*{{{*/
#ifdef DEBUG_MERGESIDE
        if (testIDList_DEBUG_MERGESIDE.find(id) !=  testIDList_DEBUG_MERGESIDE.end() )  {
            fprintf(stdout,"id=%s, mergeSide=%d, type_dataset=%d, do not merge\n", id, mergeSide, type_dataset );
            fprintf(stdout,"Profile for chain %s, matQij\n", id);
            WriteMODMProfile(i+1, aaSeq[i], matQij[i], matQijScore1[i],matQijScore2[i], AAAlphabet_Tuping, stdout);
        }
#endif
        for(i = 0; i < lengthSeqMODM; i ++) {
            for(j = 0; j < NUM_20_AA; j ++) { matMerged[i][j] = matQij[i][j];}
        }
    }/*}}}*/

#ifdef NEGLECT_DATAHANDEL
    }
#endif
    
#ifdef DEBUG_DATA_HANDEL
    finish2=clock();
    duration2 = double(finish2-start2)  /double(CLOCKS_PER_SEC);
    fprintf(stdout,"Data handel for for id %s costs %lf seconds. length= %d\n", id, duration2, lengthSeqMODM);
#endif
    return lengthSeqMODM;
}/*}}}*/
int GetDataList(vector <BYTE*> &dataList, string dbname, int maxIndex)/*{{{*/
{
    int i;
    string filename;
    struct stat st;
    int filesize=0;
    for (i=0; i<= maxIndex; i++){
        filename=dbname + int2string(i) + string(".db");
        stat(filename.c_str(), &st);
        filesize = st.st_size;
        BYTE *data = NULL;
        Array1D <BYTE> data_1darray(filesize+1);
        FILE *fp = fopen(filename.c_str(),"r");
        if (filesize > 0){
            if (fp != NULL){
                int nread = fread(data_1darray.array1D, sizeof(BYTE), filesize, fp);
                if (nread == filesize) {
                    data = new BYTE[nread];
                    memcpy(data, data_1darray.array1D, nread);
                } else {
                    fprintf(stderr,"fread error for dbfile %s in GetDataList\n",filename.c_str() );
                    exit(1);
                }
                fclose(fp);
            }else{
                fprintf(stderr,"can not open dbfile %s\n",filename.c_str() );
                exit(1);
            }
        }else{
            fprintf(stderr,"filesize < 0, for dbfile %s\n",filename.c_str() );
            exit(1);
        }
        if (data == NULL){
            fprintf(stderr,"data == NULL for dbfile %s\n",filename.c_str() );
            exit(1);
        }   else{
            dataList.push_back(data);
        }
    }
    return 0;
}/*}}}*/
int Read_database_SHU_nofseek(int dbtype, const char *trainIDListFile, int NPer_Frag_Database, int ratioScheme, const char *train_Qijpath,const char *train_modmpath, const char *train_fragaccpath,  const char *resultpath,  char **idListTrain, char **aaSeqAllTrain,char **shapeSeqAllTrain,char **secSeqAllTrain,int **matAllTrain,int **digitAASeqAllTrain, int *lengthListAllTrain ,int dsspMapMethod /*= 0*/, bool isReadBinaryFile /*= true*/, int8 typeProfile  /*=1*/, set <string> &testIDList_DEBUG_MERGESIDE)/*{{{*/
/*Modified by Nanjiang, using the function ReadInDatabase*/
{
    fprintf(stdout,"Read in matrix for the training set...\n");
#ifdef DEBUG_TIME
    clock_t start, finish;
    double duration;
    start = clock();
#endif
    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;

    char id[SIZE_ID+1] = "";
    Array2D <int> matQij_2darray(LONGEST_SEQ, NUM_20_AA);
    Array2D <int> matFrag_2darray(LONGEST_SEQ, NUM_20_AA);
    Array2D <int> matMODM_2darray(LONGEST_SEQ, NUM_20_AA);
    Array2D <int> matMerged_2darray(LONGEST_SEQ, NUM_20_AA);
    int ** matQij = matQij_2darray.array2D; /*Qij profile*/
    int ** matFrag= matFrag_2darray.array2D; /*fragacc profile*/
    int ** matMODM = matMODM_2darray.array2D; /*modm profile, modm and Qij might be the same, 2008-01-18, Nanjiang*/
    int ** matMerged = matMerged_2darray.array2D;

    matQij_2darray.Init(0);
    matFrag_2darray.Init(0);
    matMODM_2darray.Init(0);
    matMerged_2darray.Init(0);

    Array1D <float> matFragScore1_1darray(LONGEST_SEQ);
    Array1D <float> matFragScore2_1darray(LONGEST_SEQ);
    Array1D <float> matQijScore1_1darray(LONGEST_SEQ);
    Array1D <float> matQijScore2_1darray(LONGEST_SEQ);
    Array1D <float> matMODMScore1_1darray(LONGEST_SEQ);
    Array1D <float> matMODMScore2_1darray(LONGEST_SEQ);
    float * matFragScore1 = matFragScore1_1darray.array1D;
    float * matFragScore2 = matFragScore2_1darray.array1D;
    float * matQijScore1 = matQijScore1_1darray.array1D;
    float * matQijScore2 = matQijScore2_1darray.array1D;
    float * matMODMScore1 = matMODMScore1_1darray.array1D;
    float * matMODMScore2 = matMODMScore2_1darray.array1D;

    Array1D <int> digitAASeq_1darray(LONGEST_SEQ);
    Array1D <char> aaSeq_1darray(LONGEST_SEQ+1);
    Array1D <char> shapeSeq_1darray(LONGEST_SEQ+1);
    Array1D <char> dsspSecSeq_1darray(LONGEST_SEQ+1);
    int *digitAASeq = digitAASeq_1darray.array1D;
    char *aaSeq = aaSeq_1darray.array1D;
    char *shapeSeq = shapeSeq_1darray.array1D;
    char *dsspSecSeq = dsspSecSeq_1darray.array1D;

    int qijformat =  QIJ_FORMAT_NANJIANG;
    int modmformat = MODM_FORMAT_NANJIANG;
    int fragaccformat = FRAGACC_FORMAT_NANJIANG;
    int type_dataset = TRAINING_SET;
    int mergeSide = TRAINING_SET;

    map <string, dbindex> dbindexfragacc;
    map <string, dbindex> dbindexqij;
    map <string, dbindex> dbindexmodm;
    int maxDBIndexNumber_fragacc = 0;
    int maxDBIndexNumber_qij = 0;
    int maxDBIndexNumber_modm = 0;
    vector <BYTE*> dataList_fragacc;
    vector <BYTE*> dataList_qij;
    vector <BYTE*> dataList_modm;
    string filename;

#ifdef DEBUG_TIME
    clock_t start1, finish1;
    double duration1;
    start1 = clock();
#endif
    if (dbtype==1){ /*read index file of the dumped database, added 2011-10-13*/
        if (ReadDatabaseIndex(string(train_fragaccpath), dbindexfragacc, maxDBIndexNumber_fragacc) == -1 ){
            fprintf(stderr,"Read db fragacc failed for %s\n", train_fragaccpath);
            exit(1);
        }
        if (ReadDatabaseIndex(string(train_Qijpath), dbindexqij, maxDBIndexNumber_qij) == -1 ){
            fprintf(stderr,"Read db qij failed for %s\n", train_Qijpath);
            exit(1);
        }
        if (ReadDatabaseIndex(string(train_modmpath), dbindexmodm, maxDBIndexNumber_modm) == -1 ){
            fprintf(stderr,"Read db modm failed for %s\n", train_modmpath);
            exit(1);
        }
        GetDataList( dataList_fragacc, string(train_fragaccpath), maxDBIndexNumber_fragacc);
        
        GetDataList( dataList_qij, string(train_Qijpath), maxDBIndexNumber_qij);
        GetDataList( dataList_modm, string(train_modmpath), maxDBIndexNumber_modm);
    }
    
#ifdef DEBUG_TIME
    finish1 = clock();
    duration1 = double(finish1-start1)  /double(CLOCKS_PER_SEC);
    fprintf(stdout,"Read index file cost %lf seconds\n", duration1);
#endif

#ifdef DEBUG_TIME_NOFSEEK
    if(1){
        clock_t start, finish;
        double duration;
        start = clock();

        int longestlinelength=0;
        int maxRawIDCnt = fgetlinecnt(trainIDListFile, longestlinelength, false);
        Array2D<char> idList_2darray(maxRawIDCnt, longestlinelength+1);
        char ** idList = idList_2darray.array2D;
        int numChainTrain=ReadInIDList(trainIDListFile, idList);
        Array1D <char> alphabet(100);
        int length;
        Array1D <ProfileSADByte> profile(1000);
        Array1D <double> parameter(10);
        int8 typeProfile;
        BYTE *pData;
        dbindex *pDBindex;
        for (int i=0; i<numChainTrain;i++) {
            pDBindex = &(dbindexqij[idList[i]]) ;
            pData = dataList_qij[pDBindex->dbfileindex]+pDBindex->offset ;
            GetBinaryMODM_nofseek(pData, alphabet.array1D, length, profile.array1D, parameter.array1D, typeProfile);
        }
        for (int i=0; i<numChainTrain;i++) {
            pDBindex = &(dbindexfragacc[idList[i]]) ;
            pData = dataList_fragacc[pDBindex->dbfileindex]+pDBindex->offset ;
            GetBinaryMODM_nofseek(pData, alphabet.array1D, length, profile.array1D, parameter.array1D, typeProfile);
        }
        for (int i=0; i<numChainTrain;i++) {
            pDBindex = &(dbindexmodm[idList[i]]) ;
            pData = dataList_modm[pDBindex->dbfileindex]+pDBindex->offset ;
            GetBinaryMODM_nofseek(pData, alphabet.array1D, length, profile.array1D, parameter.array1D, typeProfile);
        }

        finish=clock();
        duration = double(finish-start)  /double(CLOCKS_PER_SEC);
        fprintf(stdout,"Run GetBinaryMODM_nofseek %d times cost %lf seconds\n", numChainTrain*3, duration);
    }
#endif 

    FILE *fpTrainList = fopen(trainIDListFile,"r");
    char OpenName[MAX_PATH+1] = "";
    sprintf(OpenName,"%s/subname_list.txt",resultpath);
    FILE *wfp = fopen(OpenName,"w");
    int lengthSeq = 0;

    int cntID = 0;
    while((linesize = fgetline(fpTrainList, line ,maxline)) != EOF) //idList, main interation
    {
        if(linesize <= 0) { continue; }
        Array2D <char> tmpstrs_2darray(2, linesize+1);
        char ** tmpstrs = tmpstrs_2darray.array2D;
        int status_sscanf = 0;
        
        status_sscanf = sscanf(line, "%s %s", tmpstrs[0], tmpstrs[1]);
        if (status_sscanf == 1){
            my_strcpy(id, tmpstrs[0],SIZE_ID);
        }else if (status_sscanf == 2){
            my_strcpy(id, tmpstrs[1],SIZE_ID);
        }else {
            fprintf(stderr,"Wrong idlist line \"%s\". Ignore.\n",line);
            continue;
        }
        matFragScore1_1darray.Init(0.0);
        matFragScore2_1darray.Init(0.0);
        matQijScore1_1darray.Init(0.0);
        matQijScore2_1darray.Init(0.0);
        matMODMScore1_1darray.Init(0.0);
        matMODMScore2_1darray.Init(0.0);



        if (dbtype==0){
            lengthSeq = ReadInDatabase(id, train_fragaccpath, train_Qijpath,train_modmpath, qijformat, modmformat, fragaccformat, matFrag, matQij, matMODM, matMerged, matFragScore1,matFragScore2, matQijScore1,matQijScore2,matMODMScore1,matMODMScore2, aaSeq, digitAASeq, shapeSeq, dsspSecSeq, type_dataset, isReadBinaryFile, NPer_Frag_Database, ratioScheme, mergeSide, dsspMapMethod, typeProfile);
        } else if (dbtype == 1){ /*added 2011-10-17 */
#ifdef DEBUG_TIME_EACHID
            clock_t start2, finish2;
            double duration2;
            start2=clock();
#endif
            lengthSeq = ReadInDatabase_dumpedfile_nofseek(id, dbindexfragacc, dbindexqij, dbindexmodm, dataList_fragacc, dataList_qij, dataList_modm, matFrag, matQij, matMODM, matMerged, matFragScore1,matFragScore2, matQijScore1,matQijScore2,matMODMScore1,matMODMScore2, aaSeq, digitAASeq, shapeSeq, dsspSecSeq, type_dataset, isReadBinaryFile, NPer_Frag_Database, ratioScheme, mergeSide, dsspMapMethod, typeProfile, testIDList_DEBUG_MERGESIDE);
#ifdef DEBUG_TIME_EACHID
            finish2=clock();
            duration2 = double(finish2-start2)  /double(CLOCKS_PER_SEC);
            fprintf(stdout,"Run ReadInDatabase_dumpedfile_nofseek for id %s cost %lf seconds\n", id, duration2);
#endif

        } else {
            fprintf(stderr,"dbtype = %d has not been implemented yet. Exit.\n", dbtype);
            exit(1);
        }
        if (lengthSeq < 0)
        {
            fprintf(stderr,"%s: ReadInDatabase error, lengthSeq = %d\n", id, lengthSeq);
            continue;
        }
#ifdef DEBUG_MATRIX
        fprintf(stdout,"MODM Matrix for %s (%d)\n", id, lengthSeq);
        for(i = 0 ; i < lengthSeq ; i++) {
            WriteMODMProfile(i+1, aaSeq[i], matMODM[i], 2,2, AAAlphabet_Tuping, stdout);
        }
        fprintf(stdout, "\n");
#endif

#ifdef DEBUG_MERGESIDE /*test the mergeside and read in matrix*/
        if (testIDList_DEBUG_MERGESIDE.find(id) !=  testIDList_DEBUG_MERGESIDE.end() )  {
            fprintf(stdout,"trainning set, mergeside=%d, id=%s\n", mergeSide, id);
            for(int i = 0; i < lengthSeq; i ++) {
                WriteMODMProfile(i+1, aaSeq[i], matMerged[i], matMODMScore1[i],matMODMScore2[i], AAAlphabet_Tuping, stdout);
            }
            fprintf(stdout,"\n");
        }
#endif

        idListTrain[cntID] = new char[SIZE_ID];
        aaSeqAllTrain[cntID] = new char[lengthSeq+3];
        shapeSeqAllTrain[cntID] = new char[lengthSeq+3];
        secSeqAllTrain[cntID] = new char[lengthSeq+1];
        matAllTrain[cntID] = new int [lengthSeq*20+3];
        digitAASeqAllTrain[cntID] = new int [lengthSeq+3];

        /*get aaSeq, shString, dsspSec as well as numeric aaSeq */
        int i;
        for (i=0; i<lengthSeq; i++)
        {
            aaSeqAllTrain[cntID][i] = aaSeq[i];
            digitAASeqAllTrain[cntID][i] = digitAASeq[i];
            shapeSeqAllTrain[cntID][i] = shapeSeq[i];
            secSeqAllTrain[cntID][i] = dsspSecSeq[i];
            int startPos = i*20;
            int j;
            for (j=0; j<20; j++) {
                matAllTrain[cntID][startPos+j] = matMerged[i][j];
            }
        }
        aaSeqAllTrain[cntID][lengthSeq] = '\0';
        shapeSeqAllTrain[cntID][lengthSeq] = '\0';
        secSeqAllTrain[cntID][lengthSeq] = '\0';
        lengthListAllTrain[cntID] = lengthSeq;
        my_strcpy(idListTrain[cntID], id, SIZE_ID);

        fprintf(wfp,"%4d %5s\n",cntID,id);

        cntID ++;
    }
    if (dbtype==1){
        int i;
        for (i=0;i<=maxDBIndexNumber_fragacc;i++){ 
            if (dataList_fragacc[i] != NULL) {
                delete [] dataList_fragacc[i];
            }
        }
        for (i=0;i<=maxDBIndexNumber_qij;i++){ 
            if (dataList_qij[i] != NULL) {
                delete [] dataList_qij[i];
            }
        }
        for (i=0;i<=maxDBIndexNumber_modm;i++){ 
            if (dataList_modm[i] != NULL) {
                delete [] dataList_modm[i];
            }
        }
    }
    fprintf(stdout,"Read training set finished.\n");
    fclose(fpTrainList);
    fclose(wfp);

#ifdef DEBUG_TIME
    finish = clock();
    duration = double(finish-start)  /double(CLOCKS_PER_SEC);
    fprintf(stdout,"Read in matrix for the training set cost %lf seconds\n", duration);
#endif
    return cntID; /*return the number of chains actually read in*/
}/*}}}*/

int Read_database_SHU(int dbtype, const char *trainIDListFile, int NPer_Frag_Database, int ratioScheme, const char *train_Qijpath,const char *train_modmpath, const char *train_fragaccpath,  const char *resultpath,  char **idListTrain, char **aaSeqAllTrain,char **shapeSeqAllTrain,char **secSeqAllTrain,int **matAllTrain,int **digitAASeqAllTrain, int *lengthListAllTrain ,int dsspMapMethod /*= 0*/, bool isReadBinaryFile /*= true*/, int8 typeProfile  /*=1*/, set <string> &testIDList_DEBUG_MERGESIDE)/*{{{*/
/*Modified by Nanjiang, using the function ReadInDatabase*/
{
    fprintf(stdout,"Read in matrix for the training set...\n");
#ifdef DEBUG_TIME
    clock_t start, finish;
    double duration;
    start = clock();
#endif
    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;

    char id[SIZE_ID+1] = "";
    Array2D <int> matQij_2darray(LONGEST_SEQ, NUM_20_AA);
    Array2D <int> matFrag_2darray(LONGEST_SEQ, NUM_20_AA);
    Array2D <int> matMODM_2darray(LONGEST_SEQ, NUM_20_AA);
    Array2D <int> matMerged_2darray(LONGEST_SEQ, NUM_20_AA);
    int ** matQij = matQij_2darray.array2D; /*Qij profile*/
    int ** matFrag= matFrag_2darray.array2D; /*fragacc profile*/
    int ** matMODM = matMODM_2darray.array2D; /*modm profile, modm and Qij might be the same, 2008-01-18, Nanjiang*/
    int ** matMerged = matMerged_2darray.array2D;

    matQij_2darray.Init(0);
    matMODM_2darray.Init(0);
    matFrag_2darray.Init(0);
    matMerged_2darray.Init(0);

    Array1D <float> matFragScore1_1darray(LONGEST_SEQ);
    Array1D <float> matFragScore2_1darray(LONGEST_SEQ);
    Array1D <float> matQijScore1_1darray(LONGEST_SEQ);
    Array1D <float> matQijScore2_1darray(LONGEST_SEQ);
    Array1D <float> matMODMScore1_1darray(LONGEST_SEQ);
    Array1D <float> matMODMScore2_1darray(LONGEST_SEQ);
    matFragScore1_1darray.Init(0.0);
    matFragScore2_1darray.Init(0.0);
    matQijScore1_1darray.Init(0.0);
    matQijScore2_1darray.Init(0.0);
    matMODMScore1_1darray.Init(0.0);
    matMODMScore2_1darray.Init(0.0);
    float * matFragScore1 = matFragScore1_1darray.array1D;
    float * matFragScore2 = matFragScore2_1darray.array1D;
    float * matQijScore1 = matQijScore1_1darray.array1D;
    float * matQijScore2 = matQijScore2_1darray.array1D;
    float * matMODMScore1 = matMODMScore1_1darray.array1D;
    float * matMODMScore2 = matMODMScore2_1darray.array1D;



    Array1D <int> digitAASeq_1darray(LONGEST_SEQ);
    Array1D <char> aaSeq_1darray(LONGEST_SEQ+1);
    Array1D <char> shapeSeq_1darray(LONGEST_SEQ+1);
    Array1D <char> dsspSecSeq_1darray(LONGEST_SEQ+1);
    int *digitAASeq = digitAASeq_1darray.array1D;
    digitAASeq_1darray.Init(0);
    aaSeq_1darray.Init('\0');
    shapeSeq_1darray.Init('\0');
    dsspSecSeq_1darray.Init('\0');
    char *aaSeq = aaSeq_1darray.array1D;
    char *shapeSeq = shapeSeq_1darray.array1D;
    char *dsspSecSeq = dsspSecSeq_1darray.array1D;

    int qijformat =  QIJ_FORMAT_NANJIANG;
    int modmformat = MODM_FORMAT_NANJIANG;
    int fragaccformat = FRAGACC_FORMAT_NANJIANG;
    int type_dataset = TRAINING_SET;
    int mergeSide = TRAINING_SET;

    map <string, dbindex> dbindexfragacc;
    map <string, dbindex> dbindexqij;
    map <string, dbindex> dbindexmodm;
    int maxDBIndexNumber_fragacc = 0;
    int maxDBIndexNumber_qij = 0;
    int maxDBIndexNumber_modm = 0;
    vector <FILE*> fpList_fragacc;
    vector <FILE*> fpList_qij;
    vector <FILE*> fpList_modm;
    string filename;

#ifdef DEBUG_TIME
    clock_t start1, finish1;
    double duration1;
    start1 = clock();
#endif
    if (dbtype==1){ /*read index file of the dumped database, added 2011-10-13*/
        if (ReadDatabaseIndex(string(train_fragaccpath), dbindexfragacc, maxDBIndexNumber_fragacc) == -1 ){
            fprintf(stderr,"Read db fragacc failed for %s\n", train_fragaccpath);
            exit(1);
        }
        if (ReadDatabaseIndex(string(train_Qijpath), dbindexqij, maxDBIndexNumber_qij) == -1 ){
            fprintf(stderr,"Read db qij failed for %s\n", train_Qijpath);
            exit(1);
        }
        if (ReadDatabaseIndex(string(train_modmpath), dbindexmodm, maxDBIndexNumber_modm) == -1 ){
            fprintf(stderr,"Read db modm failed for %s\n", train_modmpath);
            exit(1);
        }
        GetDBFPList( fpList_fragacc, string(train_fragaccpath), maxDBIndexNumber_fragacc);
        GetDBFPList( fpList_qij, string(train_Qijpath), maxDBIndexNumber_qij);
        GetDBFPList( fpList_modm, string(train_modmpath), maxDBIndexNumber_modm);
    }
    
#ifdef DEBUG_TIME
    finish1 = clock();
    duration1 = double(finish1-start1)  /double(CLOCKS_PER_SEC);
    fprintf(stdout,"Read index file cost %lf seconds\n", duration1);
#endif

    FILE *fpTrainList = fopen(trainIDListFile,"r");
    char OpenName[MAX_PATH+1] = "";
    sprintf(OpenName,"%s/subname_list.txt",resultpath);
    FILE *wfp = fopen(OpenName,"w");
    int lengthSeq = 0;

    int cntID = 0;
    while((linesize = fgetline(fpTrainList, line ,maxline)) != EOF) //idList, main interation
    {
        if(linesize <= 0) { continue; }
        Array2D <char> tmpstrs_2darray(2, linesize+1);
        char ** tmpstrs = tmpstrs_2darray.array2D;
        int status_sscanf = 0;
        
        status_sscanf = sscanf(line, "%s %s", tmpstrs[0], tmpstrs[1]);
        if (status_sscanf == 1){
            my_strcpy(id, tmpstrs[0],SIZE_ID);
        }else if (status_sscanf == 2){
            my_strcpy(id, tmpstrs[1],SIZE_ID);
        }else {
            fprintf(stderr,"Wrong idlist line \"%s\". Ignore.\n",line);
            continue;
        }

        if (dbtype==0){
            lengthSeq = ReadInDatabase(id, train_fragaccpath, train_Qijpath,train_modmpath, qijformat, modmformat, fragaccformat, matFrag, matQij, matMODM, matMerged, matFragScore1,matFragScore2, matQijScore1,matQijScore2,matMODMScore1,matMODMScore2, aaSeq, digitAASeq, shapeSeq, dsspSecSeq, type_dataset, isReadBinaryFile, NPer_Frag_Database, ratioScheme, mergeSide, dsspMapMethod, typeProfile);
        } else if (dbtype == 1){ /*added 2011-10-17 */
            lengthSeq = ReadInDatabase_dumpedfile(id, dbindexfragacc, dbindexqij, dbindexmodm, fpList_fragacc, fpList_qij, fpList_modm, matFrag, matQij, matMODM, matMerged, matFragScore1,matFragScore2, matQijScore1,matQijScore2,matMODMScore1,matMODMScore2, aaSeq, digitAASeq, shapeSeq, dsspSecSeq, type_dataset, isReadBinaryFile, NPer_Frag_Database, ratioScheme, mergeSide, dsspMapMethod, typeProfile);

        } else {
            fprintf(stderr,"dbtype = %d has not been implemented yet. Exit.\n", dbtype);
            exit(1);
        }
        if (lengthSeq < 0)
        {
            fprintf(stderr,"%s: ReadInDatabase error, lengthSeq = %d\n", id, lengthSeq);
            continue;
        }
#ifdef DEBUG_MATRIX
        fprintf(stdout,"MODM Matrix for %s (%d)\n", id, lengthSeq);
        for(i = 0 ; i < lengthSeq ; i++) {
            WriteMODMProfile(i+1, aaSeq[i], matMODM[i], 2,2, AAAlphabet_Tuping, stdout);
        }
        fprintf(stdout, "\n");
#endif

#ifdef DEBUG_MERGESIDE /*test the mergeside and read in matrix*/
        if (testIDList_DEBUG_MERGESIDE.find(id) !=  testIDList_DEBUG_MERGESIDE.end() )  {
            fprintf(stdout,"trainning set, mergeside=%d, id=%s\n", mergeSide, id);
            for(int i = 0; i < lengthSeq; i ++) {
                WriteMODMProfile(i+1, aaSeq[i], matMerged[i], matMODMScore1[i],matMODMScore2[i], AAAlphabet_Tuping, stdout);
            }
            fprintf(stdout,"\n");
        }
#endif

        idListTrain[cntID] = new char[SIZE_ID];
        aaSeqAllTrain[cntID] = new char[lengthSeq+3];
        shapeSeqAllTrain[cntID] = new char[lengthSeq+3];
        secSeqAllTrain[cntID] = new char[lengthSeq+1];
        matAllTrain[cntID] = new int [lengthSeq*20+3];
        digitAASeqAllTrain[cntID] = new int [lengthSeq+3];

        /*get aaSeq, shString, dsspSec as well as numeric aaSeq */
        int i;
        for (i=0; i<lengthSeq; i++)
        {
            aaSeqAllTrain[cntID][i] = aaSeq[i];
            digitAASeqAllTrain[cntID][i] = digitAASeq[i];
            shapeSeqAllTrain[cntID][i] = shapeSeq[i];
            secSeqAllTrain[cntID][i] = dsspSecSeq[i];
            int startPos = i*20;
            int j;
            for (j=0; j<20; j++)
            {
                matAllTrain[cntID][startPos+j] = matMerged[i][j];
            }
        }
        aaSeqAllTrain[cntID][i] = '\0';
        shapeSeqAllTrain[cntID][i] = '\0';
        secSeqAllTrain[cntID][i] = '\0';
        lengthListAllTrain[cntID] = lengthSeq;
        my_strcpy(idListTrain[cntID], id, SIZE_ID);

        fprintf(wfp,"%4d %5s\n",cntID,id);

        cntID ++;
    }
    if (dbtype==1){
        int i;
        for (i=0;i<=maxDBIndexNumber_fragacc;i++){ fclose(fpList_fragacc[i]); }
        for (i=0;i<=maxDBIndexNumber_qij;i++){ fclose(fpList_qij[i]); }
        for (i=0;i<=maxDBIndexNumber_modm;i++){ fclose(fpList_modm[i]); }
    }
    fprintf(stdout,"Read training set finished.\n");
    fclose(fpTrainList);
    fclose(wfp);

#ifdef DEBUG_TIME
    finish = clock();
    duration = double(finish-start)  /double(CLOCKS_PER_SEC);
    fprintf(stdout,"Read in matrix for the training set cost %lf seconds\n", duration);
#endif
    return cntID; /*return the number of chains actually read in*/
}/*}}}*/

int search_frag(const char *testIDListFile, const char *trainIDListFile, const char *test_Qijpath, const char *train_Qijpath,int qijformat, const char * test_modmpath, const char *train_modmpath, int modmformat, const char *test_fragaccpath, const char *train_fragaccpath, int fragaccformat, const char *resultpath, int **subMatrix, double *back_comp , int fragSize , set <string> &testIDList_DEBUG_MERGESIDE)/*{{{*/
/*****************************************************************************
 * Derived from Tuping's program
 ****************************************************************************/
{
#ifdef DEBUG_TIME
    clock_t start, finish;
    double duration;
#endif
    //char id[SIZE_ID+1] = ""; [>name of the chain<]
    char idTar[SIZE_ID+1] = ""; /*name of the target chain in the test set*/

    //int tmpi = 0; [>temporary variable for integer<]
    //float tmpf = 0.0; [>temporary variable for floating value<]
    int   i,j;
    int   im;
    int   TContinue;
    int numChainTrain = 0; /**/


    Array1D <int> digitAASeq_1darray(LONGEST_SEQ);
    Array1D <char> aaSeq_1darray(LONGEST_SEQ+1);
    Array1D <char> shapeSeq_1darray(LONGEST_SEQ+1);
    Array1D <char> dsspSecSeq_1darray(LONGEST_SEQ+1);
    int *digitAASeq = digitAASeq_1darray.array1D;
    char *aaSeq = aaSeq_1darray.array1D;
    char *shapeSeq = shapeSeq_1darray.array1D;
    char *dsspSecSeq = dsspSecSeq_1darray.array1D;

    int   iksub;
    DATATYPE_LOGPER   ScoreSum = 0; /*score between two fragments*/
    DATATYPE_LOGPER   ScoreTtemp;
    int   Nsearhigh;


    int linesize;
    int maxline = 1000;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;

    //DATATYPE_LOGPER param1;
    DATATYPE_LOGPER value;
    DATATYPE_LOGPER Log_per[105]; /*logarithm value for calculating the profile profile score, log2(x)*/
    DATATYPE_LOGPER Per_ten_one[105];
    DATATYPE_LOGPER Log_Back_Comp[21];

    Array2D <int> matQij_2darray(LONGEST_SEQ, NUM_AA);
    Array2D <int> matFrag_2darray(LONGEST_SEQ, NUM_AA);
    Array2D <int> matMODM_2darray(LONGEST_SEQ, NUM_AA);
    Array2D <int> matMerged_2darray(LONGEST_SEQ, NUM_AA);
    int ** matQij = matQij_2darray.array2D; /*Qij profile*/
    int ** matFrag= matFrag_2darray.array2D; /*fragacc profile*/
    int ** matMODM = matMODM_2darray.array2D; /*modm profile, modm and Qij might be the same, 2008-01-18, Nanjiang*/
    int ** matMerged = matMerged_2darray.array2D;

    matQij_2darray.Init(0);
    matMODM_2darray.Init(0);
    matMerged_2darray.Init(0);

    Array1D <float> matFragScore1_1darray(LONGEST_SEQ);
    Array1D <float> matFragScore2_1darray(LONGEST_SEQ);
    Array1D <float> matQijScore1_1darray(LONGEST_SEQ);
    Array1D <float> matQijScore2_1darray(LONGEST_SEQ);
    Array1D <float> matMODMScore1_1darray(LONGEST_SEQ);
    Array1D <float> matMODMScore2_1darray(LONGEST_SEQ);
    float * matFragScore1 = matFragScore1_1darray.array1D;
    float * matFragScore2 = matFragScore2_1darray.array1D;
    float * matQijScore1 = matQijScore1_1darray.array1D;
    float * matQijScore2 = matQijScore2_1darray.array1D;
    float * matMODMScore1 = matMODMScore1_1darray.array1D;
    float * matMODMScore2 = matMODMScore2_1darray.array1D;

    Array2D <DATATYPE_LOGPER> MA_MB_2darray(105,105);
    Array2D <DATATYPE_LOGPER> MA_P_2darray(105,20);
    DATATYPE_LOGPER ** MA_MB = MA_MB_2darray.array2D;
    DATATYPE_LOGPER ** MA_P = MA_P_2darray.array2D;
    DATATYPE_LOGPER ma, mb, pp;

    for (i=0; i<103; i++) /*Log_per*/
    {
        for(j = 0; j < 103; j ++)
        {
            ma = DATATYPE_LOGPER ( (i + 0.25) * PER_TEN_SCALE ); 
            mb = DATATYPE_LOGPER ( log(j+0.25)/M_LN2 * LOG_PER_SCALE);/*use log(2,x), 2 based logarithm, thus the result is scaled in bit, Nanjiang, 2008-01-15*/
            MA_MB[i][j] = ma*mb;
        }
    }

    for(i = 0 ; i < 103; i ++)
    {
        for (j=0; j<20; j++)
        { 
            ma = DATATYPE_LOGPER ( (i + 0.25) * PER_TEN_SCALE ); 
            pp = DATATYPE_LOGPER ( log(back_comp[j])/M_LN2 *LOG_PER_SCALE);//see Per_ten_one[i]
            MA_P[i][j] = ma*pp;
        }
    }

    //candidate
    int  halfFragSize = fragSize/2;
    

    for (i=0; i<103; i++) /*Log_per*/
    {
        Per_ten_one[i] = DATATYPE_LOGPER ( (i + 0.25) * PER_TEN_SCALE ); 
        Log_per[i] = DATATYPE_LOGPER ( log(i+0.25)/M_LN2 * LOG_PER_SCALE);/*use log(2,x), 2 based logarithm, thus the result is scaled in bit, Nanjiang, 2008-01-15*/
    }

    for (i=0; i<20; i++)
    { 
        Log_Back_Comp[i] = DATATYPE_LOGPER ( log(back_comp[i])/M_LN2 *LOG_PER_SCALE);//see Per_ten_one[i]
#ifdef DEBUG_CHECK_FRAGSCORE
        fprintf(stdout,"back_comp[%d]=%.3f\n", i, back_comp[i]);
#endif
    }
    
    int maxNumIDTrain = fgetlinecnt(trainIDListFile); /*first get the maximal number of ids in the training set, the ids in trainIDListFile should be one per line, 2008-01-20, Nanjiang*/
    
    /*allocate the memory for variables*/
    Array1D <char*> idListTrain_1darray(maxNumIDTrain+1);
    Array1D <char*> aaSeqAllTrain_1darray(maxNumIDTrain+1);
    Array1D <char*> shapeSeqAllTrain_1darray(maxNumIDTrain+1);
    Array1D <char*> secSeqAllTrain_1darray(maxNumIDTrain+1);
    Array1D <int*> digitAASeqAllTrain_1darray(maxNumIDTrain+1);
    Array1D <int*> matAllTrain_1darray(maxNumIDTrain+1);
    Array1D <int> lengthListAllTrain_1darray(maxNumIDTrain+1);

    char **idListTrain        = idListTrain_1darray.array1D;        /*list of chain names for all chains in the training set */
    char **aaSeqAllTrain      = aaSeqAllTrain_1darray.array1D;      /*amino acid sequences for all chains in the training set*/
    char **shapeSeqAllTrain   = shapeSeqAllTrain_1darray.array1D;   /*shape strings for all chains in the training set*/
    char **secSeqAllTrain     = secSeqAllTrain_1darray.array1D;     /*secondary structure sequences for all chains in the training set*/
    int  **digitAASeqAllTrain = digitAASeqAllTrain_1darray.array1D; /*digit encoded (A-->0, ...)amino acid sequences for all chains in the training set */
    int  **matAllTrain        = matAllTrain_1darray.array1D;        /*merged PSSMs for all chains in the training set */
    int   *lengthListAllTrain = lengthListAllTrain_1darray.array1D; /*length of sequences for all chains in the training set*/
    for (i=0; i<maxNumIDTrain+1; i++) { 
        idListTrain[i]        = NULL;
        aaSeqAllTrain[i]      = NULL;
        shapeSeqAllTrain[i]   = NULL;
        secSeqAllTrain[i]     = NULL;
        digitAASeqAllTrain[i] = NULL;
        matAllTrain[i]        = NULL;
        lengthListAllTrain[i] = 0;
    }

/*==Frame: read in training set=====================
 * Read in the profiles for the training set, the profiles are stored in
 * matMerge, which is a combination of Frag matrices and MODM matrices, 
 * For positions where Frag matices are not valid, use Qij profiles instead*/
   numChainTrain = Read_database_SHU(dbtype, trainIDListFile, NPer_Frag_Database, ratioScheme,  train_Qijpath, train_modmpath, train_fragaccpath,  resultpath,  idListTrain, aaSeqAllTrain,shapeSeqAllTrain, secSeqAllTrain, matAllTrain, digitAASeqAllTrain, lengthListAllTrain, dsspMapMethod, isReadBinaryFile , typeProfile, testIDList_DEBUG_MERGESIDE);
   //numChainTrain = Read_database_SHU_nofseek(dbtype, trainIDListFile, NPer_Frag_Database, ratioScheme,  train_Qijpath, train_modmpath, train_fragaccpath,  resultpath,  idListTrain, aaSeqAllTrain,shapeSeqAllTrain, secSeqAllTrain, matAllTrain, digitAASeqAllTrain, lengthListAllTrain, dsspMapMethod, isReadBinaryFile , typeProfile, testIDList_DEBUG_MERGESIDE);
/*===End of reading the training set ==================================================*/


    /*==Frame: read in test set, search the chain in the test set in the training set============ */
    int lengthTar = 0; /*length of the target sequence (in the test set)*/

    FILE *fpTestList = NULL;
    fpTestList = fopen(testIDListFile, "r");
    checkfilestream(fpTestList, testIDListFile, "r");

    int cnttestfile = 0 ;
    fprintf(stdout, "Search frag for the test set...\n");
    while((linesize = fgetline(fpTestList, line ,maxline)) != EOF)// idList for test set/*{{{*/
    {
        if (linesize <=0) {continue;}
        cnttestfile ++;
        if (linesize <= 0 || cnttestfile <(beginID+1)|| cnttestfile > endID)
        {
            continue; /*run only the test idTar from begin to end, 2007-11-16, Nanjiang Shu*/
        }
        my_strcpy(idTar, line, SIZE_ID);
        
        fprintf(stdout, "%-5d search frag for the target %s...\n", cnttestfile, idTar);

        /*===Read in matrix for the test chain ===============*/
        int type_dataset = TEST_SET;
        int lengthSeq = 0;
        lengthSeq = ReadInDatabase(idTar, test_fragaccpath, test_Qijpath, test_modmpath, qijformat, modmformat, fragaccformat, matFrag, matQij, matMODM, matMerged, matFragScore1,matFragScore2, matQijScore1,matQijScore2,matMODMScore1,matMODMScore2, aaSeq, digitAASeq, shapeSeq, dsspSecSeq, type_dataset, isReadBinaryFile,NPer_Frag_Database, ratioScheme, mergeSide, dsspMapMethod , typeProfile);
        if(lengthSeq < 0) {
            continue;
        }
        /*===End read in matrix for the test chain ==========*/
#ifdef DEBUG_MERGESIDE /*test the mergeside and read in matrix*/
        if (testIDList_DEBUG_MERGESIDE.find(id) !=  testIDList_DEBUG_MERGESIDE.end() )  {
            fprintf(stdout,"test set,  mergeside=%d, id=%s\n", mergeSide, idTar);
            for(i = 0; i < lengthSeq; i ++) {
                WriteMODMProfile(i+1, aaSeq[i], matMerged[i], matMODMScore1[i],matMODMScore2[i], AAAlphabet_Tuping, stdout);
            }
            fprintf(stdout,"\n");
        }
#endif
        //int iArray_Two_seq_score; [>index to the profile-profile-score<] 
        /*===Frame: Database searching ===============*/
        lengthTar = lengthSeq;
        Array1D <int> cntCandidate_1darray(lengthTar);
        Array2D <DATATYPE_LOGPER> Candidate_Score_2darray(lengthTar, SCore_Sample+1);
        Array2D <int> Candidate_Score_iksub_2darray(lengthTar, SCore_Sample+1);
        Array2D <int> Candidate_Score_iklen_2darray(lengthTar, SCore_Sample+1);
        int *cntCandidate = cntCandidate_1darray.array1D;
        DATATYPE_LOGPER **Candidate_Score = Candidate_Score_2darray.array2D;/*storing the fragment-fragment score for each candidate fragment*/
        int **Candidate_Score_iksub = Candidate_Score_iksub_2darray.array2D;/*the position in the candidate chain*/
        int **Candidate_Score_iklen = Candidate_Score_iklen_2darray.array2D;/*the position in the target chain*/
        /*init*/
        cntCandidate_1darray.Init(0);
        for (i=0; i<lengthTar; i++)
        {
            for (j=0; j<SCore_Sample+1; j++) { Candidate_Score[i][j] = INIT_FRAGSCORE; }
        }

#ifdef DEBUG_TIME
        start = clock();
#endif
        int iTar = 0; /*iTar is the iterator for the target sequence*/
        int jCan = 0; /*jCan is the iterator for the candidate sequence*/
        for (iksub=0; iksub<numChainTrain; iksub++)/*{{{*/ /*iksub is the iterator for the chains in the training set. Sumunit is the numProTrain, the number of proteins in the training set*/
        {    
#ifdef DEBUG_MATRIX
            fprintf(stdout,"matrices target %s (%d) -- candidate %s (%d)\n", idTar, lengthTar, idListTrain[iksub], lengthListAllTrain[iksub]);
#endif
            if (strcmp(idTar, idListTrain[iksub]) == 0 )
            {
                continue; /*ignore the target chain itself. this is used for leave-one-out crossvalidation, 2008-02-18, Nanjiang*/
            }
            //***first calculate the m X n for profile-profile comparison for two sequences
            int lengthCan = lengthListAllTrain[iksub]; /*length of the current candicate sequence*/
            //int Num_two_seq = lengthCan*lengthTar + 1; [> = m * n, m,n are sequence lengths<]
            Array2D <DATATYPE_LOGPER> Two_seq_score_2darray(lengthTar, lengthCan);
            Two_seq_score_2darray.Init(INIT_FRAGSCORE);
            DATATYPE_LOGPER **Two_seq_score = Two_seq_score_2darray.array2D;
            //for (i=0; i<Num_two_seq; i++) { Two_seq_score[i] = INIT_FRAGSCORE; }

            for (iTar=0; iTar<lengthTar; iTar++) /*iTar is the iterator for the target sequence*/
            {
                if(isPredictSecShape)/*{{{*/   /*only when predicting shape string and secondary structures, these non standard amino acid and no shape positions should be removed*/ 
                { 
                    if (   (digitAASeq[iTar]>=20) || (shapeSeq[iTar]=='-')   )
                    { 
                        continue;
                    }
                }/*}}}*/
                int posCanStart = 0; /*the start position for calculating the profile-profile score*/
                int posCanEnd = lengthCan;/*the end position (< END) for calculating the profile-profile score*/

                if(speedRate > 0.0) /*if the speed up is set, apply these code, 2008-01-31, Nanjiang*/
                {
                    //double lengthProp =   min(lengthTar, lengthCan) / double(max(lengthCan, lengthTar));
                    //double avgLen = double (lengthTar+lengthCan) /2.0;
//                    if(lengthProp >= 0.5 || (lengthCan < 150 && lengthTar < 150)) 
                    {
                        int shift = Integer (min(lengthCan, lengthTar)*(speedRate));/*2009-06-18, changed, see the excel file debug_ratio3.xls for more details*/
                        posCanStart = max(iTar - lengthTar + shift, 0);
                        posCanEnd = min(iTar + lengthCan - shift , lengthCan);
                    }
                }
#ifdef DEBUG_SPEEDRATE
                fprintf(stdout,"0 %s iTar=%d lenTar %d speedrate=%.2f pS %d  pE %d  lenCan %d\n", idListTrain[iksub], iTar, lengthTar, speedRate, posCanStart, posCanEnd, lengthCan);
#endif

                //int startPosForTwoSeqScore = iTar * lengthCan; [>for speeding up<]
                //int a0 = 0;
                int *pMatTar = NULL;
                int *pMatCan = NULL;
                pMatTar = matMerged[iTar];
                for (jCan=posCanStart; jCan<posCanEnd; jCan++) /*jCan is the iterator for the candidate sequence*/
                {
                    if(isPredictSecShape)/*{{{*/
                    {
                        if  (   (digitAASeqAllTrain[iksub][jCan]>=20) || (shapeSeqAllTrain[iksub][jCan]=='-')   )
                        {
                            continue;
                        }
                    }/*}}}*/
                    //score
                    //a0 = jCan*20;[>the 2D array stored as 1D array,a0 is the start position, since each profile has 20 numbers, <]
                    pMatCan = matAllTrain[iksub]+jCan*20;
//                    value = Prob_score(Per_ten_one, Log_per, Log_Back_Comp, pMatTar, pMatCan); /*this is not faster*/
                    value = 0; 
                    for (im=0;im<20; im++) /*this part cost 60% of the computational time*/
                    //for(im = 20; im --;) /*this does not speed up the program either*/
                    {
                        /*calculate PICASO3 score, Nanjiang,2008-01-13*/
                        //param1 = Per_ten_one[matMerged[iTar][im]]*(Log_per[matAllTrain[iksub][a0+im]]-Log_Back_Comp[im]) + Per_ten_one[matAllTrain[iksub][a0+im]]*(Log_per[matMerged[iTar][im]]-Log_Back_Comp[im]);
                        value +=(Per_ten_one[pMatTar[im]]*(Log_per[pMatCan[im]]-Log_Back_Comp[im]) + Per_ten_one[pMatCan[im]]*(Log_per[pMatTar[im]]-Log_Back_Comp[im]));
                        //value +=(MA_MB[pMatTar[im]][pMatCan[im]] - MA_P[pMatTar[im]][im] + MA_MB[pMatCan[im]][pMatTar[im]] - MA_P[pMatCan[im]][im]);[>this is even slower<]

#ifdef DEBUG_CHECK_FRAGSCORE
                        fprintf(stdout,"M1[%2d][%2d] (%3d)-- M2[%2d][%2d] (%3d) = %.2f\n", iTar, im,matMerged[iTar][im], jCan, im, matAllTrain[iksub][a0+im], param1);
#endif /*DEBUG_CHECK_FRAGSCORE*/


#ifdef DEBUG_MATRIX
                        //fprintf(stdout,"%3d, %3d : log2(pTar[%2d]) = %.2f, log2(pCan[%2d]) = %.2f, pTar[%2d] = %d, pCan [%2d] = %d, log2(P[%2d]) = %.2f, value=%.3f\n", iTar, jCan, im, Per_ten_one[matMerged[iTar][im]], im, Log_per[matAllTrain[iksub][a0+im]], im, matMerged[iTar][im], im, matAllTrain[iksub][a0+im], im, Log_Back_Comp[im], value);
#endif
                    }
                    //iArray_Two_seq_score = iTar*lengthCan + jCan;
                    //Two_seq_score[iArray_Two_seq_score] = value;
                    if (profileScoreType == PROFILESCORETYPE_TUPING)
                    {
                        Two_seq_score[iTar][jCan] = DATATYPE_LOGPER ( 5000 + value*15 );
                    }
                    else /*default*/
                    {
                        Two_seq_score[iTar][jCan] = DATATYPE_LOGPER (value *  FRAGSCORE_SCALE );
                    }

                    //Two_seq_score[iTar][jCan] = value;

#ifdef DEBUG_CHECK_FRAGSCORE
                    fprintf(stdout,"subtotal M1[%2d] -- M2[%2d] = %.1f\n", iTar, jCan, Two_seq_score[iTar][jCan]);
#endif /*DEBUG_CHECK_FRAGSCORE*/
#ifdef DEBUG_MATRIX
                    fprintf(stdout,"Two_seq_score[%3d][%3d] = %.1f\n", iTar, jCan, Two_seq_score[iTar][jCan]);
#endif
                }
            }
            //****first
#ifdef DEBUG_MATRIX
            fprintf(stdout,"matrix %s\n", idTar);
            for(i = 0; i < lengthTar; i ++)
            {
                WriteMODMProfile(i+1, aaSeq[i], matMerged[i], 0,0, AAAlphabet_Tuping, stdout);
            }
            fprintf(stdout,"\n");
            fprintf(stdout,"matrix %s\n", idListTrain[iksub]);
            for(i = 0 ; i < lengthListAllTrain[iksub]; i++)
            {
                WriteMODMProfile(i+1, aaSeqAllTrain[iksub][i], &matAllTrain[iksub][i*20], 1,1, AAAlphabet_Tuping, stdout);
            }
#endif
            /*calculate score for fragSize-long fragment
            * using diag, so that the computational complex is unreated to the fragSize, 2008-02-26
            * 
            *   |------------
            *   |\  \
            *   | \  \
            *   |  \  \
            */
            for (iTar=0; iTar<lengthTar-fragSize+1; iTar++) /*bug fixed , +1, so that the last residue is calculated as well, 2008-01-30, Nanjiang iTar is the iterator for the position of target sequence*//*{{{*/
            {
                if(isPredictSecShape)/*{{{*/
                {
                    TContinue = 1;
                    for (j=0; j<fragSize; j++)
                    {
                        if (   digitAASeq[iTar+j] >= 20    )
                        {  
                            TContinue = 0;
                            break;
                        }
                    }
                    if (  shapeSeq[iTar+halfFragSize] == '-'  ) /*for the
                                                                 homology detection, the requirement of valid shape simbol is not
                                                                 necessary, this is not shape string prediction, 2008-01-30, Nanjiang*/
                    {
                        //TContinue = 0;
                    }
                    if (  TContinue == 0  )
                    {
                        continue;
                    }
                }/*}}}*/
                int posCanStart = 0; /*the start position for calculating the profile-profile score*/ 
                int posCanEnd = lengthCan - fragSize +1 ;/*the end position (< END) for calculating the profile-profile score*/
                if(speedRate > 0.0) /*if the speed up is set, apply these code, 2008-01-31, Nanjiang*/
                {
                    //double lengthProp =   min(lengthTar, lengthCan) / double(max(lengthCan, lengthTar));
                    //double avgLen = double (lengthTar+lengthCan) /2.0;
//                    if(lengthProp >= 0.5 || (lengthCan < 150 && lengthTar < 150)) 
                    {
                        int shift = Integer (min(lengthCan, lengthTar)*(speedRate));/*2009-06-18, changed, see the excel file debug_ratio3.xls for more details*/
                        posCanStart = max(iTar - lengthTar + shift, 0);
                        posCanEnd = min(iTar + lengthCan - shift , lengthCan - fragSize + 1);
                    }
                }
#ifdef DEBUG_SPEEDRATE
                fprintf(stdout,"1 %s iTar=%d lenTar %d speedrate=%.2f pS %d  pE %d  lenCan %d\n", idListTrain[iksub], iTar, lengthTar, speedRate, posCanStart, posCanEnd, lengthCan);
#endif
                for (jCan=posCanStart; jCan < posCanEnd; jCan++)/*jCan is the iterator for the candidate sequence*/
                {
                    if(isPredictSecShape)/*{{{*/
                    {
                        //check AAS and shape
                        TContinue = 1;
                        int end = jCan + fragSize; /*end position*/
                        for (j=jCan; j< end; j++)
                        {
                            if  (   (digitAASeqAllTrain[iksub][j]>=20) || (shapeSeqAllTrain[iksub][j]=='-')   )
                            {
                                TContinue = 0;
                                break;
                            }
                        }
                        if (  TContinue == 0  )
                        {
                            continue;
                        }
                    }/*}}}*/
                    //score
                    ScoreSum = DATATYPE_LOGPER(0.0); /* store the score between two fragments*/
                    for (j=0; j<fragSize; j++)
                    {
                        //int iArray_Two_seq_score = (iTar+j)*lengthCan + jCan + j;
#ifdef WITH_ASSERT
                        assert (Two_seq_score[iTar+j][jCan+j] != INIT_FRAGSCORE);
#endif
                        ScoreSum += Two_seq_score[(iTar+j)][jCan+j];
                    }

                    //check
                    if (  ScoreSum <=  Candidate_Score[iTar][0]  ) /*ignore those fragment-fragment scores lower than the minimal Candidate_Score, Nanjiang, 2008-01-13 */
                    {
                        continue;
                    }
                    if(cntCandidate[iTar] < SCore_Sample)
                    {
                        cntCandidate[iTar] ++;
                    }
                    if(cntCandidate[iTar] < SCore_Sample)
                    {
                        Candidate_Score[iTar][cntCandidate[iTar]] = ScoreSum; /*kick the smallest Candidate_Score out, Nanjiang, 2008-01-13*/
                        Candidate_Score_iksub[iTar][cntCandidate[iTar]] = iksub; /*chain index for the high-scoring fragment, Nanjiang, 2008-01-13 */
                        Candidate_Score_iklen[iTar][cntCandidate[iTar]] = jCan; /*sequence index for the high-scoring fragment, Nanjiang, 2008-01-13 */ 
                    }
                    else
                    {
                        Candidate_Score[iTar][0] = ScoreSum; /*kick the smallest Candidate_Score out, Nanjiang, 2008-01-13*/
                        Candidate_Score_iksub[iTar][0] = iksub; /*chain index for the high-scoring fragment, Nanjiang, 2008-01-13 */
                        Candidate_Score_iklen[iTar][0] = jCan; /*sequence index for the high-scoring fragment, Nanjiang, 2008-01-13 */ 

                        //put the smallest score in Candidate_Score[iTar][0]
                        ScoreTtemp = Candidate_Score[iTar][0];
                        Nsearhigh = 0;
                        for (j=1; j<SCore_Sample; j++)
                        {
                            if ( Candidate_Score[iTar][j] < ScoreTtemp )
                            {
                                ScoreTtemp = Candidate_Score[iTar][j];
                                Nsearhigh = j; /*Nsearhigh store the index for the smallest Candidate_Score, Nanjiang, 2008-01-13*/
                            }
                        }
                        if (  Nsearhigh >= 1  ) /*if the index for the smallest Candidate_Score is not zero, Nanjiang, 2008-01-13 */
                        {
                            Candidate_Score[iTar][Nsearhigh] = Candidate_Score[iTar][0];
                            Candidate_Score[iTar][0] = ScoreTtemp;

                            Swap(Candidate_Score_iksub[iTar][0], Candidate_Score_iksub[iTar][Nsearhigh]);  /*swap chain index */
                            Swap(Candidate_Score_iklen[iTar][0], Candidate_Score_iklen[iTar][Nsearhigh]); /*swap chain index */

                        }
                    }
                }//length of chain
            }//length of target/*}}}*/
        }//subunits of database/*}}}*/
#ifdef DEBUG_TIME
        finish = clock();
        duration = double(finish-start)  /double(CLOCKS_PER_SEC);
        fprintf(stdout,"searching database cost %lf seconds\n", duration);  /*99% of time is used by searching database procedures*/
#endif
        /*===End: database searching finished===============*/

        /*Sort by descending order of Candidate_Score, computational complexity, about 300*100*100 = 3e6, Nanjiang, 2008-01-13*//*{{{*/
        for (iTar=0; iTar<lengthTar-fragSize+1; iTar++) /*this cost <1% time*/
        {
            //Array1D <int> idx_score_1darray(SCore_Sample);
            //int* idx_score = idx_score_1darray.array1D;
            //for(j= 0; j<SCore_Sample; j++){ idx_score[i] = j; }
            //QuickSort_index(idx_score, Candidate_Score[iTar], 0, SCore_Sample-1, DESCENDING);
            //Array1D <DATATYPE_LOGPER> tScore_1darray(SCore_Sample);
            //Array1D <int> tiksub_1darray(SCore_Sample);
            //Array1D <int> tiklen_1darray(SCore_Sample);
            //for(j = 0; j < SCore_Sample; j ++)
            //{
            //    tScore_1darray.array1D[j] = Candidate_Score[iTar][idx_score[j]];
            //    tiksub_1darray.array1D[j] = Candidate_Score_iksub[iTar][idx_score[j]];
            //    tiklen_1darray.array1D[j] = Candidate_Score_iklen[iTar][idx_score[j]];
            //}
            //for(j = 0; j < SCore_Sample; j ++)
            //{
            //     Candidate_Score[iTar][j] = tScore_1darray.array1D[j] ;
            //     Candidate_Score_iksub[iTar][j] = tiksub_1darray.array1D[j] ;
            //     Candidate_Score_iklen[iTar][j] = tiklen_1darray.array1D[j] ;
            //}

            for (j=0; j<SCore_Sample; j++)
            {
               for (im=j+1; im<SCore_Sample; im++)
               {
                   if (  Candidate_Score[iTar][im] > Candidate_Score[iTar][j]  )
                   {
                       Swap(Candidate_Score[iTar][im], Candidate_Score[iTar][j]);  /*swap frag-frag score*/
                       Swap(Candidate_Score_iksub[iTar][im], Candidate_Score_iksub[iTar][j]); /*swap chain index*/
                       Swap(Candidate_Score_iklen[iTar][im], Candidate_Score_iklen[iTar][j]); /*swap candidate position index*/
                   }
               }
            }
        }/*}}}*/

        /*Write out frag file, the table including all high-scoring fragments*//*{{{*/
        char fragFile[MAX_PATH+1] = ""; /*frag file, storing matched fragment pair from which the dot plot can be drawn*/
        if (fragformat == FRAGFORMAT_NANJIANG)
        {
            if(isWriteBinaryFile)
            { sprintf(fragFile,"%s/%s.fragbin",resultpath,idTar); /*2007-11-11, Nanjiang, change the name of fragFile from idTar.txt to idTar.frag*/ }
            else
            { sprintf(fragFile,"%s/%s.frag",resultpath,idTar); /*2007-11-11, Nanjiang, change the name of fragFile from idTar.txt to idTar.frag*/ }
        }
        else /*(fragformat == FRAGFORMAT_TUPING)*/
        {
            if(isWriteBinaryFile)
            { sprintf(fragFile,"%s/%s.txtbin",resultpath,idTar); /*2009-06-23 */ }
            else
            { sprintf(fragFile,"%s/%s.txt",resultpath,idTar); /*2009-06-23 */ }
        }

        if (isWriteBinaryFile)
        {
            WriteBinaryFragFile(fragFile,fragformat, lengthTar, Save_Sample, fragSize, Candidate_Score, Candidate_Score_iksub, Candidate_Score_iklen, idListTrain);
        }
        else
        {
            WriteFragFile(fragFile, fragformat, lengthTar, Save_Sample, fragSize, Candidate_Score, Candidate_Score_iksub, Candidate_Score_iklen, idListTrain);
        }
        /*}}}*/

        /*predict secondary structures and shape strings*/ 
        if(isPredictSecShape) /*execute the code only when predict secondary structure option is enabled, Nanjiang, 2008-01-13*/ 
        {
            PredictSecondaryStructure(idTar, aaSeq,  Candidate_Score, Candidate_Score_iksub, Candidate_Score_iklen, lengthTar, fragSize, lengthListAllTrain,  idListTrain, secSeqAllTrain, shapeSeqAllTrain, numChainTrain, resultpath);
        }

        fprintf(stdout, "  result file output to %s\n", resultpath);

    }/*}}}*/
    fclose(fpTestList);//unclosed file stream detected, added 2007-06-20
/*===End of reading and predicting on the test set ==================================================*/


    /*free memory*/
    for (i=0; i<maxNumIDTrain+1; i++)
    {
        if (aaSeqAllTrain[i] != NULL) {delete [] aaSeqAllTrain[i];}
        if (shapeSeqAllTrain[i] != NULL){delete [] shapeSeqAllTrain[i];}
        if (secSeqAllTrain[i] != NULL ){delete [] secSeqAllTrain[i];} //memory leak detected, should be delete, added 2007-06-20 
        if (matAllTrain[i] != NULL) {delete [] matAllTrain[i];}
        if (digitAASeqAllTrain[i] != NULL){delete [] digitAASeqAllTrain[i];}
        if (idListTrain[i] != NULL){delete [] idListTrain[i];}
    }

    return 0;   
}
/*}}}*/

int main(int argc, char** argv)/*{{{*/
{
    // char *testIDList_DEBUG_MERGESIDE[] = { "2CH5A", "1F7UA", "1UX5A", "2O0MA", "2HC1A", "1YKWA" };
    //char *testIDList_DEBUG_MERGESIDE[] = { "1YKWA" , "1WY3A", "2G50A"};
    // char *testIDList_DEBUG_MERGESIDE[] = { "d1fvka2", "d1bed_1"};
    //char *testIDList_DEBUG_MERGESIDE[] = { "1D1", "1OS6A", "1M1QA"}
    //char *testIDList_DEBUG_MERGESIDE[] = { "19HCA", "1OS6A", "1M1QA"};
    set <string> testIDList_DEBUG_MERGESIDE;
    testIDList_DEBUG_MERGESIDE.insert("1D2TA") ;
    testIDList_DEBUG_MERGESIDE.insert("1E9GA");

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

    const char control_option[] = ""; //options which control the program, and does not take parameters
    int fragSize = 9; /*the size of fragment, default = 9*/

    // char test_Qijpath[MAX_PATH+1] = "D:\\database\\Par_Qijnanjiang\\";
    // char trainIDListFile[MAX_PATH+1] = "D:\\database\\train_nanjinag.txt";
    // char CParameter_conser_control[MAX_PATH+1] = "D:\\Par_QijOther\\control_parameters.txt";
    // char train_modmpath[MAX_PATH+1] = "D:\\database\\modpssm\\modmatrix_";
    // char train_fragaccpath[MAX_PATH+1] = "D:\\database\\Fragacc_compnanjiang\\frag_";
    // char train_Qijpath[MAX_PATH+1] = "D:\\database\\Par_Qijnanjiang\\Qijmatrix_";
    // char resultpath[MAX_PATH+1] = "D:\\database\\Result_nanjiang\\";
    // char Back_composition[MAX_PATH+1] = "D:\\Par_QijOther\\composition_nr.txt";

    char testIDListFile[MAX_PATH+1] = "";
    char trainIDListFile[MAX_PATH+1] = "";

    int qijformat = QIJ_FORMAT_NANJIANG;
    int modmformat = MODM_FORMAT_NANJIANG;
    int fragaccformat = MODM_FORMAT_NANJIANG; /*added 2009-06-09*/

    char test_Qijpath[MAX_PATH+1] = "/misc/casiodata3/wk/passe/Par_Qijnanjiang"; //path of Qij files for the test set
    char train_Qijpath[MAX_PATH+1] = "/misc/casiodata3/wk/passe/Par_Qijnanjiang";  //path of Qij files for the training set
    char CParameter_conser_control[MAX_PATH+1] = "";
    char test_modmpath[MAX_PATH+1] = "";
    char test_fragaccpath[MAX_PATH+1] = "";


    char train_modmpath[MAX_PATH+1] = "/misc/casiodata3/wk/passe/modpssm"; //modm files, be careful about the alphabet order
    char train_fragaccpath[MAX_PATH+1] = "/misc/casiodata3/wk/passe/Fragacc_compnanjiang"; //output from database_build, 
    char resultpath[MAX_PATH+1] = "/misc/casiodata3/wk/passe/try1/result";

    char bkfreqfile[MAX_PATH+1] = "";
    char submatfile[MAX_PATH+1] = "";

    i = 1;
    while(i < argc)/*{{{*/
    {
        if(argv[i][0] == '-' && !isNonOptionArg) //options
        {
            isNonOptionArg = false;
            if(IsInCharSet(argv[i][1], control_option))//if argv[i][1] is in control_option, it might be used as -aqs
            {
                for(j = 1 ; j < int(strlen(argv[i])); j++)
                {
                    switch (argv[i][j])
                    {
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
            else if(strcmp(argv[i],"-H") == 0 )
            {
                PrintVerboseHelp();
                return 0;
            } else if (string(argv[i])== "--dbtype" || string(argv[i])== "-dbtype") {
                if( ( i = option_parser_numeric(argc, argv, i, dbtype, true, 0, 1)) == -1){
                    return -1;
                }
            }
            else if (strcmp(argv[i], "--rb") == 0)
            {
                isReadBinaryFile = true;
                i ++;
            }
            else if (strcmp(argv[i], "--wb") == 0)
            {
                isWriteBinaryFile = true;
                i ++;
            }
            else if( (strcmp(argv[i],"--test") == 0) )  
            {
                if( ( i = option_parser_filename(argc, argv, i, testIDListFile)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--train") == 0) )  
            {
                if( ( i = option_parser_filename(argc, argv, i, trainIDListFile)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--matrix") == 0) )  
            {
                if( ( i = option_parser_filename(argc, argv, i, submatfile)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--qijformat") == 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, qijformat, true, 0, 1)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--fragaccformat") == 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, fragaccformat, true, 0, 1)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--bkfile") == 0) )  
            {
                if( ( i = option_parser_filename(argc, argv, i, bkfreqfile)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--para") == 0) )  
            {
                if( ( i = option_parser_filename(argc, argv, i, CParameter_conser_control)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--test-qij") == 0) )  
            {
                if( ( i = option_parser_filename(argc, argv, i, test_Qijpath)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--train-qij") == 0) )  
            {
                if( ( i = option_parser_filename(argc, argv, i, train_Qijpath)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--result") == 0) || string(argv[i]) == "--outpath")  
            {
                if( ( i = option_parser_filename(argc, argv, i, resultpath)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--modm") == 0) || string(argv[i])== "--train-modm")  
            {
                if( ( i = option_parser_filename(argc, argv, i, train_modmpath)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--test-modm") == 0) )  
            {
                if( ( i = option_parser_filename(argc, argv, i, test_modmpath)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--modmformat") == 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, modmformat, true, 0, 1)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--fragacc") == 0) || string(argv[i])=="--train-fragacc")
            {
                if( ( i = option_parser_filename(argc, argv, i, train_fragaccpath)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--test-fragacc") == 0) )  
            {
                if( ( i = option_parser_filename(argc, argv, i, test_fragaccpath)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--NPer") == 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, NPer_Frag_Database, true, 0, 100)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--scoretype") == 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, profileScoreType, true, 0, 1)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--fragformat") == 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, fragformat, true, 0, 1)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "-N") == 0) || strcmp(argv[i], "--topN")==0)  
            {
                if( ( i = option_parser_numeric(argc, argv, i, SCore_Sample, true, 0, 1000)) == -1)
                    return -1;
                Save_Sample = SCore_Sample;
            }
            else if( (strcmp(argv[i], "-w") == 0) || strcmp(argv[i], "--fragsize")==0)  
            {
                if( ( i = option_parser_numeric(argc, argv, i, fragSize, true, 0, 100)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "-s") == 0) || strcmp(argv[i], "--speed")==0)  
            {
                if( ( i = option_parser_numeric(argc, argv, i, speedRate, true, float(0.0), float(1.0))) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--mergeside") == 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, mergeSide, true, 0, 1)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--ratioscheme") == 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, ratioScheme, true, 0, 10)) == -1)
                    return -1;
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
            i ++;
        }
    }/*}}}*/

    if (strcmp(test_modmpath, "") == 0) /*Added 2009-06-09*/ {
        my_strcpy(test_modmpath, train_modmpath, MAX_PATH);
    }
    if (strcmp(test_fragaccpath, "") == 0) {
        my_strcpy(test_fragaccpath, train_fragaccpath, MAX_PATH);
    }

    VerifyFolder(resultpath);

    Array2D <int> subMatrix_2darray(NUM_BLOSUM, NUM_BLOSUM);
    Array1D <double > back_comp_1darray(NUM_BLOSUM);
    int **subMatrix  = subMatrix_2darray.array2D;
    double *back_comp = back_comp_1darray.array1D;
    char alphabet[NUM_BLOSUM+1] = "";

    if(strcmp(submatfile,"") == 0)
    {
        for(i = 0 ; i < NUM_BLOSUM; i++)
            for(j = 0 ; j  < NUM_BLOSUM ; j++)
                subMatrix[i][j] = blosum62[i][j];
    }
    else
    {
        GetPatternTransMatrix(submatfile, subMatrix, alphabet);
    }
    //transfer the order of matrix
    char blosum_aa_alphabet[] = "ARNDCQEGHILKMFPSTWYV";
    ReorderMatrix(subMatrix,  blosum_aa_alphabet, AAAlphabet_Tuping);

    //background composition--back_comp
    if(strcmp(bkfreqfile,"")==0)
    {
        for(i = 0 ; i < NUM_20_AA; i++) back_comp[i] = bkFreq[i]*100.0;
    }
    else
    {
        GetBkFreq(bkfreqfile,back_comp, AAAlphabet_Tuping, NUM_20_AA);
        for(i = 0 ; i < NUM_20_AA; i ++) { back_comp[i] *=  100.0; }// convert to percentage
    }

    //conser_parameters.txt---search number: sample number:
    if(strcmp(CParameter_conser_control,"") != 0)
    {
        FILE *fp = fopen(CParameter_conser_control,"r");
        if (fp != NULL){
            while (  fscanf(fp,"%d %d %d %d %d %d\n", &SCore_Sample, &Save_Sample, &Consens_Sample, &NPer_Frag_Database, &Nsearch_beg , &Nsearch_beg) != EOF  ) {
            }
            fclose(fp);
        } else{
            fprintf(stderr,"Can not open file para_control_file \"%s\"\n", CParameter_conser_control);
        }
    }
    //==test 2009-06-09 
    //fprintf(stdout,"SCore_Sample = %d\n", SCore_Sample);
    //fprintf(stdout,"Save_Sample = %d\n", Save_Sample);
    //fprintf(stdout,"Consens_Sample = %d\n", Consens_Sample);
    //fprintf(stdout,"NPer_Frag_Database = %d\n", NPer_Frag_Database);
    //==test 2009-06-09 

#ifdef DEBUG_TIME_MAIN
    start = clock();
#endif
    search_frag(testIDListFile, trainIDListFile, test_Qijpath, train_Qijpath, qijformat, test_modmpath, train_modmpath, modmformat, test_fragaccpath, train_fragaccpath, fragaccformat, resultpath, subMatrix, back_comp, fragSize, testIDList_DEBUG_MERGESIDE);

#ifdef DEBUG_TIME_MAIN
        finish = clock();
        duration = double(finish-start)  /double(CLOCKS_PER_SEC);
        fprintf(stdout,"search_frag cost %lf seconds\n", duration);
#endif


    return 0;
}
