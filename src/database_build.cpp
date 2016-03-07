#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include "array.h"
#include "mytemplate.h"
#include "myfunc.h"
#include "mypro.h"

#ifndef QIJ_NAME_FORMAT
#define QIJ_NAME_FORMAT
#define QIJ_FORMAT_TUPING 0
#define QIJ_FORMAT_NANJIANG 1
#endif
#undef MAX_SIZE_ID
#define MAX_SIZE_ID 100

#undef MAX_NUM_CHAIN
#define MAX_NUM_CHAIN 7000

#ifndef INLINE
#  if __GNUC__
#    define INLINE extern inline
#  else
#    define INLINE inline
#  endif
#endif

#define DATATYPE_DIGITAASEQ int8          /*change the datatype from int to int8 does not save computational time, 2008-02-15, Nanjiang*/
#define DATATYPE_DIGITSHAPESEQ int8
#define DATATYPE_WATERACC int8
#define DATATYPE_ENCODEWATERACC int8
#define DATATYPE_MODE int8

/* ==========format of the fragFile (*.dbfrag) ====================
 * fragFile store the fragAcc profile for each sliding 9-long fragment,
 * the format of the fragfile (*.dbfrag) is
 * _______________________________________________________________
 * column 1    : index of the frag, [0, lengthSeq - 9]
 * column 2    : index within the frag, repeating from 0 to 8
 * column 3-22 : profile for each position in the 9-long fragment
 * column 23   : number of candidates found for this fragment
 * column 24   : average score for all candidates to this fragment
 * --------------------------------------------------------------
 * 2008-02-14 
 * ****************************************************************/

//int idtype = 0; [>default idtype = 0<] /*idtype option removed 2008-02-13 */

int FRAG_WINDOW_SIZE = 9;
int MAX_NUM_SAMEAA_FRAG = 1; /*maximal number of fragments with the same amino acid sequence allowed to build the database*/
bool isReadBinaryFile = false;
bool isWriteBinaryFile = false;
bool isCreateTriple = false; /*whether to create the triple frequency of shape string, do not read from pre-created file*/

int beginID = 0;       /*in order to run the program simultaneously, setting the begin and end position of the idlist to run*/
int endID = 0x7FFFFFFF; /*by default, running all ids in the idListFile, 2007-11-16 */


//ChangeLog /*{{{*/
/*****************************************************************************
 * ChangeLog 2007-10-22
 * option --idtype added
 *
 * ChangLog 2007-11-02
 * the file namning scheme for fragFileName and fragAccFileName changed, using
 * the extension "dbfrag" and "fragacc" to distinguish these two files
 *
 * ChangeLog 2007-11-03 
 * output binary file also, reading binary file is about 100 times faster then
 * reading the ascii file
 *
 * ChangeLog 2007-11-16
 * add option --create-triple,  the triple frequency of the shape string can be
 * created according to the database.
 *
 * ChangeLog 2007-11-16
 * add the functionality to check the redundency of the fragments, if there are
 * many fragment sharing the same amino acid sequence, use only
 * MAX_NUM_SAMEAA_FRAG = 1 of them to build the profile
 *
 * add options --begin --end 
 * Note: in building the triple frequency, U shape was grouped into S, V to R
 * and G to T
 *  S (+U)
 *  R (+V)
 *  K
 *  A
 *  T (+G)
 *
 * ChangeLog 2008-02-13
 *   many variables has been re-written. 
 *   option --idtype removed. All ids should be of exact file rootname. 
 * ChangeLog 2008-02-15
 *   bug fixed on database_build verion before 2008-02-08, in which NSample is
 *   not updated when (NAAS <=100), see database_build-2008-02-08.cpp
 *   However, this bug will not cause any error when the database is large so
 *   that (NAAS) is always > 100. 
 *
 *   speed up the database_build by at least 10 times. 
 *   the result has compared with version 2008-02-08, the same
 *
 *   later, all integer division was changed to first using the double
 *   calculation, then round to integer. This will cause some round off
 *   difference
 *
 * ChangeLog 2009-04-27
 *   add the debugging code #ifdef DEBUG_WRITE_FRAG
 *   which print out all matched fragments for each target fragment.
 * ChangeLog 2009-11-08
 *   the format of the binary MODM file changed. 
 *   and the --wb function added
 ****************************************************************************//*}}}*/
/*****************************************************************************
 * valgrind checked, no leaks are possible, 2008-02-14 , Nanjiang Shu
 ****************************************************************************/
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

const char  shapeAlphabet_8s[] = "SRUVKATG "; /*8-state shape string alphabet*/
const char  shapeAlphabet_5s[] = "SRKAT ";  /*5-state shape string alphabet*/

//const DATATYPE_DIGITAASEQ AAS_Code[] = [>{{{<]
//{
    //0         , //A--ALA
    //20        , //B--Non
    //13        , //C--CYS
    //18        , //D--ASP
    //16        , //E--GLU
    //5         , //F--PHE
    //10        , //G--GLY
    //9         , //H--HIS
    //3         , //I--ILE
    //20        , //J--Non
    //7         , //K--LYS
    //2         , //L--LEU
    //6         , //M--MET
    //15        , //N--ASN
    //20        , //O--Non
    //4         , //P--PRO
    //19        , //Q--GLN
    //8         , //R--ARG
    //11        , //S--SER
    //12        , //T--THR
    //20        , //U--Non
    //1         , //V--VAL
    //17        , //W--TRP
    //20        , //X--Non
    //14        , //Y--TYR
    //20          //Z--Non
//};
/*}}}*/
const DATATYPE_DIGITSHAPESEQ Shape_Code[] = /*{{{*/
{
    /*the shapeAlphabet_8s = "SRUVKATG "*/
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
const DATATYPE_DIGITSHAPESEQ Shape_convert[] = /*map 8-shape to 5-shape*//*{{{*/
{
    0, /* S --> S */     
    1, /* R --> R */     
    0, /* U --> S */     
    1, /* V --> R */     
    2, /* K --> K */     
    3, /* A --> A */     
    4, /* T --> T */     
    4, /* G --> T */     
    8  /* ' ' --> ' ' */     
};/*}}}*/
const DATATYPE_DIGITSHAPESEQ Shape_compare[] = /*map 5-shape to 3-shape, for comparing shape symbols to build shape blocks*//*{{{*/
{   /*(S,R) in one group, (A,K) in one group and (T) in the third group */
    0, /* S --> 0 */
    0, /* R --> 0 */ 
    1, /* K --> 1 */
    1, /* A --> 1 */
    2  /* T --> 3 */ 
};/*}}}*/

void PrintHelp()
{
    fprintf(stdout,"Usage: build-database [options] idListFile\n");
    fprintf(stdout,"    Note: one id per line in the idListFile\n");
    fprintf(stdout,"options:\n");
    fprintf(stdout,"  -o|--out file    : output the result to outfile, default = stdout\n");
    fprintf(stdout,"  --rb             : read the binary format file matrix, *.Qijbin, *.modmbin, *.fragaccbin\n");
    fprintf(stdout,"  --wb             : write the binary format file, *.fragbin and *.secpredbin file \n");
    fprintf(stdout,"  --qijpath path   : path for Qij files\n");
    fprintf(stdout,"  --qijformat 0|1  : format of qijfile name, 0 - Qijmatrix_$id.txt, 1- $id.Qij, default = 0\n");
    fprintf(stdout,"  --result  path   : path for result files, default  = ./\n");
    fprintf(stdout,"  --resultacc path : path for result acc files, FINAL OUTPUT\n");
    fprintf(stdout,"  --triple file    : shapeTriple file,  default = $WORKDIR/passe/Shape_triple3829.txt\n");
    fprintf(stdout,"  --create-triple  : create the triple frequency, do not read from the Shape_triple3829 file\n");
    fprintf(stdout,"  --matrix file    : substitution matrix file, default = blosum62\n");
    //fprintf(stdout,"  --idtype   0|1   : 0 -- standardized id, 1 -- scop seq id\n");
    fprintf(stdout,"  --begin   <int>  : start number of ids in the idListFile to run database build\n");
    fprintf(stdout,"  --end     <int>  : end number of ids in the idListFile to run database build\n");
    fprintf(stdout,"                   : if begin and end is not set, running all IDs in the idListFile, id[endID] will not be run\n");
    fprintf(stdout,"                   : that is, --begin = 10, --end = 12, will run the the 11th to 12th item items\n");
    fprintf(stdout,"  -h|--help        : print this help message and exit\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"Created on 2007-06-06, Updated on 2009-11-08, Nanjiang Shu\n");
}
void PrintVerboseHelp() { }

void ReorderMatrix(int **M, const char *alphabet_from, const char *alphabet_to)/*{{{*/
{
    int sizeAlphabet = strlen(alphabet_from);
    assert(strlen(alphabet_from) == strlen(alphabet_to));
    Array2D <int> Mtmp_2darray(sizeAlphabet, sizeAlphabet);
    int **Mtmp = Mtmp_2darray.array2D;
    int i,j;
    for(i = 0 ; i < sizeAlphabet; i++) 
        for(j = 0 ; j < sizeAlphabet; j++)
            Mtmp[i][j] = M[i][j] ;

    Array1D <int> mapindex_1darray(sizeAlphabet);
    int *mapindex = mapindex_1darray.array1D;
    for(i = 0 ; i < sizeAlphabet; i++)
    {
        mapindex[i] = Char2Digit(alphabet_to[i],  alphabet_from);
    }

    for(i = 0 ; i < sizeAlphabet; i++) 
        for(j = 0 ; j < sizeAlphabet; j++)
        {
            M[i][j] = Mtmp[mapindex[i]][mapindex[j]];
        }
}
/*}}}*/
int ReadInTripleShapeFreq(const char *shapeTripleFile, int *** freqTriple, int numShapeState = 5)/*{{{*/
/*****************************************************************************
 * Read in the triple-shape frequency from the file
 * format of the shapeTripleFile
 * code1 code2 code3 triple-shape frequency
 * 0 0 0 SSS 851
 * 0 0 1 SSR 328
 * 0 0 2 SSK 35
 *  
 * 2008-02-13 
 ****************************************************************************/
{
    int sumFreq = 0; /*total frequency of occupancies of triple-shapes*/
    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;

    int i,j,k;
    int freq;
    char triple[100+1] = "";
    FILE *fpTriple = fopen(shapeTripleFile,"r");
    if(checkfilestream(fpTriple, shapeTripleFile, "r") == -1)
    {
        return -1;
    }
    int cntLine = 0;
    while ((linesize = fgetline(fpTriple,line, maxline))!= EOF)
    {
        if(linesize > 0 && line[0] != '#' ) /*this should be && instead of ||, bug fixed 2009-07-27*/
        {
            if (sscanf(line,"%d %d %d %s %d",&i,&j,&k,triple,&freq ) == 5)
            {
                if(i>=0 && i <numShapeState && j>= 0 && j < numShapeState && k >= 0&& k <numShapeState)
                {
                    freqTriple[i][j][k] = freq;
                    sumFreq += freq;
                }
                else
                {
                    fprintf(stderr,"Error! shapeCode overflow for file %s at line %d:\n'%s'\n", shapeTripleFile, cntLine+1, line);
                    return -1;
                }
            }
            else
            {
                fprintf(stderr,"Error! shapeTripleFile %s does not contain 5 items at line %d\n'%s'\n", shapeTripleFile, cntLine+1, line);
                return -1;
            }
        }
        cntLine ++;
    }
    fclose(fpTriple);

    /*normalize to ten-thousand percentage*/
    for (i=0; i<numShapeState; i++)
    {
        for (j=0; j<numShapeState; j++)
        {
            for (k=0; k<numShapeState; k++)
            {
                freqTriple[i][j][k] = Integer (freqTriple[i][j][k]*10000.0/sumFreq);
                //freqTriple[i][j][k] = int (freqTriple[i][j][k]*10000.0/sumFreq + 0.5);
            }
        }
    }

    return EXIT_SUCCESS;
}/*}}}*/
int CreateShapeTripleFreq(DATATYPE_DIGITSHAPESEQ **digitShapeAll, int *lengthListAll, int numChain, char **idListAll, int ***freqTriple, DATATYPE_DIGITSHAPESEQ numShapeState, const char* tripleFreqOutFile)/*{{{*/
/*****************************************************************************
 * given the digitShapeAll, lengthListAll, numChain and numChain
 * calculated the frequency of triple-shapes and store them in freqTriple
 ****************************************************************************/
{
    int i,j,k;
    for(i = 0; i < numChain; i ++)
    {
        int lenShapeString = lengthListAll[i];
        if(lenShapeString <= 0)
        {
            fprintf(stderr,"Error! %s, lenShapeString = %d <= 0\n", idListAll[i], lenShapeString);
        }
        for(j = 0; j < lenShapeString-2; j ++)
        {
            if(digitShapeAll[i][j] < numShapeState && digitShapeAll[i][j+1] < numShapeState && digitShapeAll[i][j+2] < numShapeState)
            {
                freqTriple[digitShapeAll[i][j]][digitShapeAll[i][j+1]][digitShapeAll[i][j+2]] ++;
            }
        }
    }

    /*normalize the freqTriple*/
    int sumTripleFreq = 0;
    for(i = 0; i < numShapeState; i ++) 
    { for(j = 0; j < numShapeState ; j ++) 
        { for(k = 0; k < numShapeState; k ++) 
            { sumTripleFreq += freqTriple[i][j][k]; } 
        } 
    }

    for (i=0; i<numShapeState; i++)
    {
        for (j=0; j<numShapeState; j++)
        {
            for (k=0; k<numShapeState; k++)
            {
                freqTriple[i][j][k] = Integer (freqTriple[i][j][k]*10000.0/sumTripleFreq);
                //freqTriple[i][j][k] = int(freqTriple[i][j][k]*10000.0/sumTripleFreq+0.5);
            }
        }
    }

    /*Print normalized freqTriple to file*/
    FILE *fpwTriple = fopen(tripleFreqOutFile, "w");
    checkfilestream(fpwTriple, tripleFreqOutFile, "w");
    fprintf(fpwTriple, "#triple shape string frequency in ten thousandth\n");

    char tripleShString[5] = "";
    for(i = 0; i < numShapeState; i ++)
    {
        for(j = 0; j < numShapeState ; j ++)
        {
            for(k = 0; k < numShapeState; k ++)
            {
                sprintf(tripleShString, "%c%c%c", shapeAlphabet_5s[i], shapeAlphabet_5s[j], shapeAlphabet_5s[k]);
                fprintf(fpwTriple,"%d %d %d %s %d\n", i,j,k,tripleShString, freqTriple[i][j][k]);
            }
        }
    }
    fclose(fpwTriple);

    return EXIT_SUCCESS;
}/*}}}*/
void WriteBinaryFragAccFile(const char *outfile, char *aaSeq, char *shapeSeq,  DATATYPE_WATERACC  *waterAcc, char *secSeq, int **Number_blosum, int *Number_chain, int *Score_all, int seqLength )/*{{{*/
/*****************************************************************************
 * Write out the binary fragacc file
 ****************************************************************************/
{
    int i, j ;
    /*profile format: int8*/ 
    Array1D <ProfileSADByte> profileSADByte_1darray(seqLength);
    ProfileSADByte *profileSADByte = profileSADByte_1darray.array1D;
    int8 typeProfile = 1;

    string alphabet = "";
    alphabet = AAAlphabet_Tuping;
    Array1D <double> parameter(8);
    parameter.Init(0.0);

    int cntValidRes = 0;
    for(i = 0; i < seqLength; i ++)
    {
        if (  Number_chain[i] <= 5  ) { continue; }
        profileSADByte[cntValidRes].aaSeqIndex = short(i+1);
        profileSADByte[cntValidRes].aa = aaSeq[i];
        profileSADByte[cntValidRes].shape = shapeSeq[i];
        profileSADByte[cntValidRes].waterAcc = int8(waterAcc[i]);
        profileSADByte[cntValidRes].dsspSec = secSeq[i];

        for (j=0; j<NUM_20_AA; j++) /*j is usually used as iterator on profile columns*/
        {
            profileSADByte[cntValidRes].p[j] = Number_blosum[i][j];
        }
        profileSADByte[cntValidRes].score1 = float(Number_chain[i]/100.0);
        profileSADByte[cntValidRes].score2 = float(Score_all[i]/10.0);

        cntValidRes ++;
    }
    WriteBinaryMODM(outfile, alphabet.c_str(), cntValidRes, profileSADByte, parameter.array1D, typeProfile);
}
/*}}}*/

INLINE int IsSADValid(DATATYPE_DIGITAASEQ *digitAASeq, DATATYPE_DIGITSHAPESEQ *digitShapeSeq, DATATYPE_ENCODEWATERACC *encodeWaterAcc, int startPosFrag, int begFrag, int endFrag)/*{{{*/
/*check if the residue position contains valid amino acid, shape symbol and
 * waterAcc
 * startPosFrag: the stat position of the fragment
 * begFrag and endFrag are the index within the fragment, [0,9]*/
{
    int i;
    int beg = startPosFrag + begFrag;
    int end = startPosFrag + endFrag;
    for (i=beg; i<end; i++)
    {
        if (  (digitAASeq[i]>=20) || (digitShapeSeq[i]>=5) || (encodeWaterAcc[i] == INIT_WATERACC)  ) 
        {
            return (i-startPosFrag);
        }
    }
    return -1;
}/*}}}*/
INLINE bool IsSameShape(DATATYPE_DIGITSHAPESEQ *digitShapeTar, DATATYPE_DIGITSHAPESEQ *digitShapeCan, int posTar, int posCan, int begFrag, int endFrag, DATATYPE_MODE shapeCompareMode)/*{{{*/
/*check if the shape strings of two fragment are the same
 * digitShapeTar: shape string of the target chain
 * begFrag and endFrag are the start and end position within the fragment, that
 * is from 0 to 9
 * return -1 if the frag-shape are the same
 * return the index of the first position where the shapes are not the same
 * */
{
    int i;   /**/
    if(shapeCompareMode == 0)
    {
        for (i=begFrag; i<endFrag; i++) /*compare only the shape strings from nbeg to nend (exclude nend)*/
        {
            if (  digitShapeTar[posTar+i] != digitShapeCan[posCan+i]  ) { return false; }
        }
    }
    else /*use non-strict comparing*/
    {
        for (i=begFrag; i<endFrag; i++) /*compare only the shape strings from nbeg to nend (exclude nend)*/
        {
            if (  Shape_compare[digitShapeTar[posTar+i]] != Shape_compare[digitShapeCan[posCan+i]]  ) { return false; }
        }
    }
    return true;
}/*}}}*/
INLINE DATATYPE_DIGITAASEQ DigitAA_tuping(char aa)/*encode amino acid into digit value*//*{{{*/
{
    int tmp = 0;
    tmp = aa - 'A';
    if  (  (tmp>=0) && (tmp<=25)  )
    {
        return AAS_Code[tmp];
    }
    else
    {
        return  20;
    }
}/*}}}*/
INLINE DATATYPE_DIGITSHAPESEQ DigitShape_tuping(char shape) /*encode shape symbol into digit value*//*{{{*/
{
    int tmp = 0;
    tmp = shape - 'A';
    if  (  (tmp>=0) && (tmp<=25)  )
    {
        return Shape_Code[tmp];
    }
    else
    {
        return 8;
    }
}/*}}}*/
void EncodeWaterAcc(int *origWaterAcc, int *encodeWaterAcc, int length)/*{{{*/
/*****************************************************************************
 * Encode [0-9] waterAcc again in to [0-4]
 * 2008-02-14, Nanjiang
 ****************************************************************************/
{
    int i;
    int acc;   /*original water acc*/
    int enAcc; /*encode water acc*/
    for(i = 0; i < length; i ++)
    {
        acc = origWaterAcc[i];
        if((acc >= 0 ) && ( acc <= 9))
        {
            if (  acc == 0  ) { enAcc = 0; }      /*0     */
            else if (  acc <= 3  ) { enAcc = 1; } /*0-6   */
            else if (  acc <= 6  ) { enAcc = 2; } /*6--21 */
            else if (  acc <= 8  ) { enAcc = 3; } /*21--41*/
            else { enAcc = 4; }
        }
        else
        {
            enAcc = INIT_WATERACC;  //invalid water accessibility, (Nanjiang)
        }
        encodeWaterAcc[i] = enAcc ;
    }
}/*}}}*/
int ReadInDatabase(const char *id, const char *qijpath, int qijformat, char *aaSeq, char *shapeSeq, int *waterAcc, char *dsspSec, bool isReadBinaryFile)/*{{{*/
{
    int i;
    char qijFileName[MAX_PATH+1] = "";
    int lengthSeq = 0;
    if (! isReadBinaryFile)/*{{{*/
    {
        GetQijFilePath(id, qijFileName, qijpath, qijformat, isReadBinaryFile);
        lengthSeq = ReadInProfile(qijFileName, NULL, aaSeq, shapeSeq, waterAcc, dsspSec, NULL, NULL, NULL, LONGEST_SEQ);
        if(lengthSeq < 0) 
        { 
            fprintf(stderr,"Read file '%s' failed\n", qijFileName);
            return -1; 
        }
    }/*}}}*/
    else /*isReadBinaryFile == true*//*{{{*/
    {
        /*added 2009-11-08*/
        int8 typeProfile = 0;
        int sizeAlphabet = 0;
        int seqLength = 0;

        /*read in Qij matrix file*/
        GetQijFilePath(id, qijFileName, qijpath, qijformat, isReadBinaryFile);

        typeProfile = GetBinaryMODMPara(qijFileName, sizeAlphabet, seqLength, typeProfile);
        Array1D <char> alphabet_1darray(sizeAlphabet+1);
        char *alphabet = alphabet_1darray.array1D;
        double parameter[8];

        if (typeProfile == 0)
        {
            Array1D <ProfileSAD> profile_1darray(LONGEST_SEQ);
            ProfileSAD *profile = profile_1darray.array1D;
            int aaSeqIndex = 0;
            /*read in Qij matrix file*/
            if (GetBinaryMODM(qijFileName, alphabet, lengthSeq, profile, parameter, typeProfile) == -1 )
            {
                fprintf(stderr, "Can not open file %s\n", qijFileName);
                return -1;
            }
            for(i = 0 ; i <lengthSeq; i ++)
            {
                aaSeqIndex           = profile[i].aaSeqIndex-1;
                aaSeq[aaSeqIndex]    = profile[i].aa;           /* amino acid sequence */
                shapeSeq[aaSeqIndex] = profile[i].shape;        /* shape strings       */
                waterAcc[aaSeqIndex] = profile[i].waterAcc;     /* water accessibility */
                dsspSec[aaSeqIndex]  = profile[i].dsspSec;      /* dssp 8-state secondary structure*/
            }
        }
        else if (typeProfile == 1)
        {
            Array1D <ProfileSADByte> profile_1darray(LONGEST_SEQ);
            ProfileSADByte *profile = profile_1darray.array1D;
            int aaSeqIndex = 0;
            /*read in Qij matrix file*/
            GetQijFilePath(id, qijFileName, qijpath, qijformat, isReadBinaryFile);
            if (GetBinaryMODM(qijFileName, alphabet, lengthSeq, profile, parameter, typeProfile) == -1 )
            {
                fprintf(stderr, "Can not open file %s\n", qijFileName);
                return -1;
            }
            for(i = 0 ; i <lengthSeq; i ++)
            {
                aaSeqIndex           = profile[i].aaSeqIndex-1;
                aaSeq[aaSeqIndex]    = profile[i].aa;           /* amino acid sequence */
                shapeSeq[aaSeqIndex] = profile[i].shape;        /* shape strings       */
                waterAcc[aaSeqIndex] = profile[i].waterAcc;     /* water accessibility */
                dsspSec[aaSeqIndex]  = profile[i].dsspSec;      /* dssp 8-state secondary structure*/
            }
        }
        else
        {
            fprintf(stderr,"Error! The format of the binary Qij file %s is with typeProfile = %d, should be 0 or 1\n", qijFileName, typeProfile);
        }
    }/*}}}*/
    return lengthSeq;
}/*}}}*/

int Build_Database(const char *idListFile, const char *qijpath, int qijformat, const char *shape_triple_file, const char * resultpath, const char *result_acc_path , int **subMatrix)/*{{{*/
/*****************************************************************************
 * build substitution matrix based on the shape strings fragment and filtered
 * by the amino acid sequence and water accessibility
 * derived from Tuping's program
 *
 * For each chain in idListFile, search the shape string of this chain to all
 * chains (including itself), the fragacc matrix is built from the 9-residue
 * width shape string blocks.
 *
 * last updated on 2008-02-13
 ****************************************************************************/
{
    int ik /*iterator when level above i is needed*/;
    int i;/*i is used as the iterator for sequence positions or other general first order iterator*/
    int j;/*j is used as the iterator for profile columns or other second order iterator*/
    int status_sscanf = 0;/*store the return value of function sscanf()*/ 
    //int   SUM;    /*don't know what to do with this variable, Nanjiang,2008-02-13 */ 
    //int  *NShapedis; /*don't know what to do with this variable, Nanjiang,2008-02-13 */
    char id[MAX_SIZE_ID+1] = "";

    int linesize;
    int maxline = 1000;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;

    int iTar;       /*iterator for the target sequence                             */
    int iFragTar;   /*iterator for the fragment position of the target sequence    */
    int iCan;       /*iterator for the candidate sequence                          */
    int iFragCan;   /*iterator for the fragment position of the candidate sequence */
    //int Shape5_table[5][5];[>not used, this table is a substitution table when comparing shape symbols<]

    int HIGH = 100; /*maximal number of candidates to be used for high frequency triple-shapes */
    int LOW  = 200; /*maximal number of candidates to be used for low frequency triple-shapes  */

    int NSamplehigh = HIGH; /* 100,use less samples when the frequency of triple-shape is high, 2008-02-14, Nanjiang */
    int NSamplelow  = LOW;  /* 200,                                                                                  */
    int maxSample   = 0;    /*maximal number of candidates to be used                                                */

    int numShapeState = 5;
    Array3D <int> triple_dist_3darray(numShapeState,numShapeState,numShapeState);
    triple_dist_3darray.Init(0);
    int ***freqTriple = triple_dist_3darray.array3D; /*frequency of triple-shapestrings, shape strings are mapped to 5-state*/

    Array1D <char> aaSeq_1darray(LONGEST_SEQ+1);
    Array1D <char> shapeSeq_1darray(LONGEST_SEQ+1);
    Array1D <char> dsspSecSeq_1darray(LONGEST_SEQ+1);
    Array1D <int> waterAcc_1darray(LONGEST_SEQ+1);
    Array1D <int> encodeWaterAcc_1darray(LONGEST_SEQ+1); 
    char *aaSeq          = aaSeq_1darray.array1D;          /*amino acid sequence                          */
    char *shapeSeq       = shapeSeq_1darray.array1D;       /*shape strings                                */
    char *dsspSecSeq     = dsspSecSeq_1darray.array1D;     /*dssp secondary structure sequence            */
    int  *waterAcc       = waterAcc_1darray.array1D;       /*here watecAcc is mapped to 0-9 in SAD format */
    int  *encodeWaterAcc = encodeWaterAcc_1darray.array1D; /*encoded waterAcc, [0,4]                      */


    int maxNumIDTrain = fgetlinecnt(idListFile); /*first get the maximal number of ids in the training set, the ids in idListFile should be one per line, 2008-01-20, Nanjiang*/

    /*allocate the memory for variables*/
    Array1D <char*> idListAll_1darray(maxNumIDTrain+1);
    Array1D <char*> aaSeqAll_1darray(maxNumIDTrain+1);
    Array1D <char*> shapeSeqAll_1darray(maxNumIDTrain+1);
    Array1D <char*> dsspSecAll_1darray(maxNumIDTrain+1);
    Array1D <DATATYPE_DIGITAASEQ*> digitAASeqAll_1darray(maxNumIDTrain+1);
    Array1D <DATATYPE_DIGITSHAPESEQ*> digitShapeAll_1darray(maxNumIDTrain+1);
    Array1D <DATATYPE_WATERACC*> waterAccAll_1darray(maxNumIDTrain+1);
    Array1D <DATATYPE_ENCODEWATERACC*> encodeWaterAccAll_1darray(maxNumIDTrain+1);
    Array1D <int> lengthListAll_1darray(maxNumIDTrain+1);


    char **idListAll         = idListAll_1darray.array1D;         /*list of chain names for all chains in the training set  */
    char **aaSeqAll          = aaSeqAll_1darray.array1D;          /*amino acid sequences for all chains in the training set */
    char **shapeSeqAll       = shapeSeqAll_1darray.array1D;       /*shape strings for all chains in the training set        */
    char **dsspSecAll        = dsspSecAll_1darray.array1D;        /*dssp secondary structure sequence                       */
    DATATYPE_DIGITAASEQ  **digitAASeqAll     = digitAASeqAll_1darray.array1D;     /*digti encoded amino acid sequences for all chains       */
    DATATYPE_DIGITSHAPESEQ  **digitShapeAll     = digitShapeAll_1darray.array1D;     /*digit encoded shape strings for all chains              */
    DATATYPE_WATERACC  **waterAccAll       = waterAccAll_1darray.array1D;       /* water accessibility for all chains, [0,9]              */
    DATATYPE_ENCODEWATERACC  **encodeWaterAccAll = encodeWaterAccAll_1darray.array1D; /*encoded water accessibility for all chains, [0,4]       */
    int   *lengthListAll     = lengthListAll_1darray.array1D;     /*length of sequences for all chains in the training set  */

    //initialization/*{{{*/
    for (i=0; i<maxNumIDTrain+1; i++) 
    { 
        idListAll[i]          = NULL;
        aaSeqAll[i]           = NULL;
        shapeSeqAll[i]        = NULL;
        dsspSecAll[i]         = NULL;
        digitAASeqAll[i]      = NULL;
        digitShapeAll[i]      = NULL;
        waterAccAll[i]        = NULL;
        encodeWaterAccAll[i] = NULL;
        lengthListAll[i]      = 0;
    }
    //for (i=0; i<5; i++)/*{{{*/ /*OBSOLETE CODE*/
    //{
    //    for (j=0; j<5; j++)
    //    {
    //        Shape5_table[i][j] = 0;
    //    }
    //    if (  i == j  )
    //    {
    //        Shape5_table[i][j] = 10;
    //    }
    //    else if (  (i<=1) && (j<=1)  )
    //    {
    //        Shape5_table[i][j] = 5;
    //    }
    //    else if (  (i>=2) && (i<=3) && (j>=2) && (j<=3)  )
    //    {
    //        Shape5_table[i][j] = 5;
    //    }
    //}/*}}}*/
    //SUM = 125;
    //NShapedis = new int[SUM+1];
    //for (i=0; i<SUM; i++) { NShapedis[i] = 0; }
/*}}}*/

    if(!isCreateTriple)
    {
        if (ReadInTripleShapeFreq(shape_triple_file, freqTriple, numShapeState) == -1)
        {
            return -1;
        }
    }

    int cntChain = 0; /*count the number of chains in the database*/
    FILE *fpIDList = fopen(idListFile,"r");
    if(checkfilestream(fpIDList, idListFile ,"r") == -1)
    {
        return -1;
    }
    fprintf(stdout,"Read in Qij files...\n"); //debug
    while((linesize = fgetline(fpIDList,line, maxline))!= EOF  )/*{{{*/ /*readin Qij file for the database*/
    {
        if(linesize <= 0 || line[0] == '#') { continue; }
        my_strcpy(id, strtrim(line), MAX_SIZE_ID);

        int lengthSeq  = ReadInDatabase(id, qijpath, qijformat, aaSeq, shapeSeq, waterAcc, dsspSecSeq, isReadBinaryFile);

        EncodeWaterAcc(waterAcc, encodeWaterAcc, lengthSeq);

        lengthListAll[cntChain] = lengthSeq;
        idListAll[cntChain]         = new char[MAX_SIZE_ID+1];
        aaSeqAll[cntChain]          = new char[lengthSeq+1];
        shapeSeqAll[cntChain]       = new char[lengthSeq+1];
        dsspSecAll[cntChain]        = new char[lengthSeq+1];
        digitAASeqAll[cntChain]     = new DATATYPE_DIGITAASEQ [lengthSeq+1];
        digitShapeAll[cntChain]     = new DATATYPE_DIGITSHAPESEQ [lengthSeq+1];
        waterAccAll[cntChain]       = new DATATYPE_WATERACC [lengthSeq+1];
        encodeWaterAccAll[cntChain] = new DATATYPE_ENCODEWATERACC [lengthSeq+1];

#ifdef DEBUG
        fprintf(stdout,"lengthListAll[%d] = %d\n", cntChain, lengthSeq);
#endif

        //chain IDs
        my_strcpy(idListAll[cntChain], id, MAX_SIZE_ID);
        for (i=0; i<lengthSeq; i++)/*{{{*/
        {
            aaSeqAll[cntChain][i]      = aaSeq[i];
            shapeSeqAll[cntChain][i]   = shapeSeq[i];
            dsspSecAll[cntChain][i]    = dsspSecSeq[i];
            digitAASeqAll[cntChain][i] = DigitAA_tuping(aaSeq[i]);
            digitShapeAll[cntChain][i] = Shape_convert[DigitShape_tuping(shapeSeq[i])];
            waterAccAll[cntChain][i]   = waterAcc[i];
            encodeWaterAccAll[cntChain][i] = encodeWaterAcc[i];

        }/*}}}*/
        aaSeqAll[cntChain][lengthSeq] = '\0';
        shapeSeqAll[cntChain][lengthSeq] = '\0';
        dsspSecAll[cntChain][lengthSeq] = '\0';
        cntChain++;
//      fprintf(stdout,"%d %s result output to %s\n",cntChain, id, resultpath, fragFileName);
    }/*}}}*/
    fclose(fpIDList);
    /* ============== Read in Qij files Finished =============*/

    int numChain = cntChain; /*number of chains in the database*/
    if(isCreateTriple)/*added 2007-11-16, create triple frequency with different dtabases*//*{{{*/
    {
        char tripleFreqOutFile[MAX_PATH+1] = "";
        char rtname[MAX_PATH+1] = "";
        rootname(idListFile, rtname);
        sprintf(tripleFreqOutFile, "%s.triple.freq", rtname) ;
        CreateShapeTripleFreq(digitShapeAll, lengthListAll, numChain, idListAll, freqTriple, numShapeState, tripleFreqOutFile);
    }/*}}}*/

    /* ============= Building frag-acc matrix ============== */
    fprintf(stdout,"building frag-acc matrix for %d protein chains...\n", cntChain);
    //make fragment profiles from PDB
    char aaFragTar[100] = ""; /*amino acid sequence of the fragment for the target sequence, 2007-11-16*/
    //char aaFragCan[100] = ""; [>amino acid sequence of the fragment for the candidate sequence<] 
    char idTar[MAX_SIZE_ID+1] = ""; /*id of the target sequence*/
    char idCan[MAX_SIZE_ID+1] = ""; /*id of the candidate sequence*/
    char shapeFragTar[100] = "";/*shape strings of the fragment*/
    int NPer = 0; /*triple shape frequency*/
    int NEXT = 0; /*determining the range of fragment to be compared based on the background freq of the triple-shape*/
    int nbeg = 0; /* start position for a fragment actually used in comparison*/
    int nend = 0; /*end position for a fragment actually used in comparison*/
    int NScore = 0;/*score for a candidate fragment*/   

    int High_Score[300];
    int Pos_Chain[300];
    int Pos_length[300];

    int lengthSeqTar = 0;
    int lengthSeqCan = 0;

    endID = min(endID, cntChain);

    for (iTar = beginID; iTar < endID; iTar ++)/*{{{*/ /*adding the beginID and endID for running only part of the sequence in the database, 2007-11-16, nanjiang, */
        //for (iTar=2001; iTar<cntChain; iTar++)
    {
#ifdef DEBUG
        fprintf(stdout,"cntChain =%d ... lengthlist [%d]= %d\n", cntChain, iTar, lengthListAll[iTar]);
#endif                            
        my_strcpy(idTar, idListAll[iTar], MAX_SIZE_ID);
        lengthSeqTar = lengthListAll[iTar];
        DATATYPE_DIGITAASEQ *pDigitAASeqTar         = digitAASeqAll[iTar];
        DATATYPE_DIGITSHAPESEQ *pDigitShapeTar      = digitShapeAll[iTar];
        DATATYPE_ENCODEWATERACC *pEncodeWaterAccTar = encodeWaterAccAll[iTar];

        char fragFileName[MAX_PATH+1] = "";
        sprintf (fragFileName,"%s/%s.dbfrag",resultpath,idTar); /*dbfrag, file type for the frag file generated in database-build*/
        FILE *fpFrag = fopen(fragFileName,"w");


#ifdef DEBUG_WRITE_FRAG /*2009-04-26*/
        char fragBlockFileName[MAX_PATH+1] = ""; /*write out the frag blocks*/
        sprintf (fragBlockFileName,"%s/%s.fragblock",resultpath,idTar); /*fragBlockFile, */
        FILE *fpFragBlock = fopen(fragBlockFileName,"w");
#endif

        for(iFragTar=0; iFragTar<lengthSeqTar-9; iFragTar++)/*{{{*//*iFragTar is the iterator for the fragment position of the target sequence */
        {
            if(IsSADValid(pDigitAASeqTar, pDigitShapeTar, pEncodeWaterAccTar, iFragTar, 0, FRAG_WINDOW_SIZE) >= 0)
            {
                continue;
            }

            /*get the amino acid of the fragment*/
            int cntSameAAFrag = 0; /*couter for the matches of fragment with the same amino acid sequence, 2007-11-16.
            If there are many selected fragments having the same amino acid sequence, it's probably that the redundency of sequnece existing.
            To avoid of the potential bias of the profile caused by these repeated fragment, used only one of them in FRAGACC profile building*/
            my_strcpy(aaFragTar, aaSeqAll[iTar]+iFragTar, FRAG_WINDOW_SIZE); /*get the aaFrag for the target sequence*/
            my_strcpy(shapeFragTar, shapeSeqAll[iTar]+iFragTar, FRAG_WINDOW_SIZE); /*get the aaFrag for the target sequence*/

            //extension of shape similarity--length of match

            //NEXT is the parameter to control for length of shape comparison
            NPer = freqTriple[pDigitShapeTar[iFragTar+3]][pDigitShapeTar[iFragTar+4]][pDigitShapeTar[iFragTar+5]];
            NSamplehigh = HIGH;
            if (  NPer >= 500  )/*HHH, SSS, >= 5%*/
            {
                NEXT = 0;
                NSamplehigh = LOW; /*use more samples for HHH and SSS, 2008-02-14, Nanjiang */
            }
            else if (  NPer >= 200  )/*SSS, >= 2%*/ { NEXT = 1; }
            else if (  NPer >= 100  )/*SSR,RSS,SRS,AAK,RAA,SRR,SRA,RRR,RRS,RRA,ASS, >= 1%*/
            { NEXT = 2; }
            else
            { NEXT = 3; }

            /* if NPer >= 50, use exact shape comparison, 
             * else use a less restrict shape string comparison scheme, Shape_compare 2008-02-14, Nanjiang */
            /*difference between Nper>=50 and NPer <50
             *
             * 1. NSamplehigh (in Nper>=50), NSamplelow (in NPer<50)
             * 2. use NEXT to set the range of fragment to be compared (in
             *    NPer>=50), set exactly 3-6 as the range to compared (in Nper<50)
             * 3. use [1,8] when comparing waterAcc (in NPer>=50), use [2,7] (in Nper<50)
             *
             * */
            int fragBegAcc = 0;/*start position in fragment when calculating the score for waterAcc*/
            int fragEndAcc = 0;/*end position in fragment when calculating the score for waterAcc*/ 
            DATATYPE_MODE shapeCompareMode = 0; /*0 for exact comparing, 1 for not strict comparing*/
            if(NPer>=50) { maxSample = NSamplehigh; fragBegAcc = 1; fragEndAcc = 8; shapeCompareMode = 0;}
            else { maxSample = NSamplelow; fragBegAcc = 2; fragEndAcc = 7; shapeCompareMode = 1;}

            for (ik=0; ik<maxSample; ik++) { High_Score[ik] = -100; }

            int cntHighScoreFrag = 0; /*count the high scoring fragments, maximum is maxSample*/
            
            for (iCan=0; iCan<numChain; iCan++)/*Note that the sequence itself is also used in building the profile, 2007-11-16,Nanjiang*//*{{{*/
            {
                lengthSeqCan = lengthListAll[iCan];
                strcpy(idCan, idListAll[iCan] );

                DATATYPE_DIGITAASEQ *pDigitAASeqCan = digitAASeqAll[iCan];
                DATATYPE_DIGITAASEQ *pDigitShapeCan = digitShapeAll[iCan];
                DATATYPE_ENCODEWATERACC *pEncodeWaterAccCan = encodeWaterAccAll[iCan];
                char *pAASeqCan = aaSeqAll[iCan];

                int begFragCan = 0; 
                int endFragCan = lengthSeqCan-9;
                int isFormerSADValid = 0; /*keep the status of the former SAD checking, if isFormerSADValid <0, meaning the SAD of the former fragment is valid, non negative value stores the index of the first non-valid SAD position in the former fragment [0,8] */
                int isSADValid = 0;/*the SAD status of the current candidate fragment*/
                for (iFragCan=begFragCan; iFragCan<endFragCan; iFragCan++) /*this loop cost >95% of the time*/
                {
                    //continue;
#ifdef DEBUG_CHECK_FRAG
                    if(strcmp(idCan,"d12asa_") == 0 && iFragTar == 27 && iFragCan == 151)
                    {
                        fprintf(stdout,"%s:%d\n",idCan, iFragCan );
                    }
#endif
                    //my_strcpy(aaFragCan, aaSeqAll[iCan]+iFragCan, FRAG_WINDOW_SIZE); /*get the aaFrag for the candidate sequence, calling my_strcpy in the deepest loop is slow, especially when the copying string is short*/
//                    for(i=0;i<FRAG_WINDOW_SIZE;i++) { aaFragCan[i] = aaSeqAll[iCan][iFragCan+i]; } aaFragCan[i] = '\0';
                    //if(strcmp(aaFragTar, aaFragCan) == 0) { cntSameAAFrag ++; }
                    //if(strncmp(aaFragTar, aaSeqAll[iCan]+iFragCan, FRAG_WINDOW_SIZE) == 0) { cntSameAAFrag ++; }

                    /*do not call any function in the deepest loop*/
                    for(i = 0; i < FRAG_WINDOW_SIZE; i ++)
                    {
                        if(aaFragTar[i] != pAASeqCan[i+iFragCan])
                        { break; }
                    }
                    if(i == FRAG_WINDOW_SIZE) { cntSameAAFrag ++; }

                    if(cntSameAAFrag > MAX_NUM_SAMEAA_FRAG)
                    { continue; /*2007-11-16, Nanjiang, If there are many selected fragments having the same amino acid sequence, it's probably that the redundency of sequnece existing.
                         * To avoid of the potential bias of the profile caused by these repeated fragment, used only one of them in FRAGACC profile building*/
                    }

                    if(isFormerSADValid == 0) /*if it is the 0th position in the former fragment is not SAD valid, check the whole fragment*/
                    {
                       isSADValid = IsSADValid(pDigitAASeqCan, pDigitShapeCan, pEncodeWaterAccCan, iFragCan, 0, FRAG_WINDOW_SIZE);
                    }
                    else
                    {
                       if(isFormerSADValid < 0)/*this cost 25% of time originally, and this procedure save 25% of time: if isFormerSADValid < 0, meaning the SAD status of the former fragment is valid, so check only the 8th position in the fragment*/
                       {
                           i = iFragCan + FRAG_WINDOW_SIZE-1;
                           isSADValid = (  (pDigitAASeqCan[i]>=20) || (pDigitShapeCan[i]>=5) || (pEncodeWaterAccCan[i] == INIT_WATERACC)  ) ?(i-iFragCan): -1; 
                           //isSADValid = IsSADValid(pDigitAASeqCan, pDigitShapeCan, pEncodeWaterAccCan, iFragCan,FRAG_WINDOW_SIZE-1, FRAG_WINDOW_SIZE);
                       }
                       else /*if it is >= 1 position in the former fragment is not SAD valid, then the current fragment is of course not SAD valid, but have to decrement the index*/
                       {
                           isSADValid = isFormerSADValid -1;
                       }
                    }
                    isFormerSADValid = isSADValid;
                    if(isSADValid >= 0) /*if the current fragment is not SAD valid*/
                    { continue; } 

                    /*================== this simpler code seems to be faster*//*{{{*/
                    //if(isFormerSADValid < 0)
                    //{
                        //isSADValid = IsSADValid(pDigitAASeqCan, pDigitShapeCan, pEncodeWaterAccCan, iFragCan,FRAG_WINDOW_SIZE-1, FRAG_WINDOW_SIZE); 
                    //}
                    //else
                    //{
                        //isSADValid = IsSADValid(pDigitAASeqCan, pDigitShapeCan, pEncodeWaterAccCan, iFragCan,0, FRAG_WINDOW_SIZE);
                    //}
                    //isFormerSADValid = isSADValid;
                    //if(isSADValid >= 0) { continue; }
                    /*================*//*}}}*/
                    
                    //compare shapestring of fragment to check the similarity
                    /* NEXT:   start position within the fragment 
                     * 9-NEXT: end position within the fragment
                     * if HHH or SSS compare nine-long fragment, the comparison length is varied depends on the frequency of the central triple-shape in the 9-long fragment, 2008-02-14, Nanjiang*/


                    /*IsSameShape can not used the same method as used in
                     * IsSADValid to speed up, since in IsSADValid, checking is
                     * only on the candidate sequnece. but in IsSameShape, the
                     * shpape strings shifts, so the new comparison is totally
                     * different from the originally comparison, 2008-02-15,
                     * Nanjiang*/

                    /*do not call function, input the code here, to speed up */
                    if(shapeCompareMode == 0)
                    {
                        for (i=NEXT; i<9-NEXT; i++) /*compare only the shape strings from nbeg to nend (exclude nend)*/
                        {
                            if (  pDigitShapeTar[iFragTar+i] != pDigitShapeCan[iFragCan+i]  ) { break; }
                        }
                    }
                    else /*use non-strict comparing*/
                    {
                        for (i=NEXT; i<9-NEXT; i++) /*compare only the shape strings from nbeg to nend (exclude nend)*/
                        {
                            if (  Shape_compare[pDigitShapeTar[iFragTar+i]] != Shape_compare[pDigitShapeCan[iFragCan+i]]  ) { break;}
                        }
                    }
                    if(i < 9-NEXT) { continue; }


                    /*after here, it is not computational expensive; 2008-02-15, Nanjiang*/
                    //if (!IsSameShape(pDigitShapeTar, pDigitShapeCan, iFragTar, iFragCan, NEXT, 9-NEXT, shapeCompareMode))
                    //{ continue; }

                    //compare aas of fragment to check the similarity
                    NScore = 0;
                    for (i=0; i<FRAG_WINDOW_SIZE; i++)/*this is the deepest for loop, think about the speed , 2008-02-15, Nanjiang*/
                    {
                        NScore += subMatrix[pDigitAASeqTar[iFragTar+i]][pDigitAASeqCan[iFragCan+i]];
                    }
                    if (   NScore < 2  )/*the sum of BLOSUM score of the 9-long fragment should be >= 2*/
                    { continue; }

                    //access surface
                    int NScore_ACC = 0; /*score form the waterAcc*/
                    int diffAcc = 0; /*difference in waterAcc, ranged from 0 to 4*/
                    for (i=fragBegAcc; i<fragEndAcc; i++) /*i don't know why 1=<i<8, what about the 9th residue, 2008-02-14, Nanjiang, this is the deepest for loop, think about the speed*/
                    {
                        diffAcc = abs(pEncodeWaterAccTar[iFragTar+i]-pEncodeWaterAccCan[iFragCan+i]);
                        if (  diffAcc == 0 ) { NScore_ACC += 10; }
                        else if (  diffAcc == 1 ) { NScore_ACC += 5; }
                    }
                    //NScore += NScore_ACC/10; [>aaScore + waterAcc_Score/10<]
                    NScore += int(NScore_ACC/10.0+0.5); /*aaScore + waterAcc_Score/10*/
                    // NScore += Integer(NScore_ACC/10.0);
                    if (NScore < High_Score[0]) 
                    { 
//#ifdef DEBUG_FRAG
//                        fprintf(stderr,"NScore(%d) < High_Score[0](%d)\n", NScore, High_Score[0]);
//#endif
                        continue; 
                    }  /*High_Score[0] store the lowest score*/

                    if(cntHighScoreFrag < maxSample)
                    {
                        cntHighScoreFrag ++;   
                    }
                    if(cntHighScoreFrag < maxSample)/*when cntHighScoreFrag is < maxSample, just add this fragment sequentially from the second place in the array, 2008-02-15, Nanjiang, the first one, High_Score[0] is -100, always smallest*/
                    {
                        High_Score[cntHighScoreFrag] = NScore;
                        Pos_Chain[cntHighScoreFrag]  = iCan;
                        Pos_length[cntHighScoreFrag] = iFragCan;
                    }
                    else
                    {
//#ifdef DEBUG_FRAG
//                        fprintf(stderr,"cntHighScoreFrag(%d) >= maxSample(%d)\n", cntHighScoreFrag, maxSample);
//#endif
                        High_Score[0] = NScore;
                        Pos_Chain[0]  = iCan;
                        Pos_length[0] = iFragCan;
                        int low_score = NScore; /*low_score: used to find the smallest score in array High_Score*/
                        int idxLowScore = 0;
                        for (ik=1; ik<maxSample; ik++)/*this is the deepest for loop, think about the speed, 2008-02-15, Nanjiang*/
                        {
                            if ( High_Score[ik] < low_score )
                            {
                                low_score = High_Score[ik];
                                idxLowScore = ik; /*record the position of the lowest score in the array High_Score*/
                            }
                        }
                        if (  idxLowScore > 0 )
                        {
                            Swap(High_Score[0], High_Score[idxLowScore]); /*score*/
                            Swap(Pos_Chain[0], Pos_Chain[idxLowScore]);   /*position of the train_chain*/
                            Swap(Pos_length[0],Pos_length[idxLowScore]);  /*position of the target seq*/
                        }
                    }
                }
            }/*}}}*/

            /* ======== calculate the profile ============ *//*{{{*/
            int numValidCan = 0; /*number of valid candidates for each position*/
            Array2D <int> fragProfile_2darray(FRAG_WINDOW_SIZE, NUM_20_AA);
            fragProfile_2darray.Init(0);
            int **fragProfile = fragProfile_2darray.array2D; /*9x20 profile for each fragment*/
            /*when calculating the profile, the difference between >=50 and <50  are
             * 1. NSamplehigh in (>=50) and NSamplelow in (<50)
             * 2. when (<50), an extra code is used to further eliminate the candidate to 100
             * question, since "HHH" and "SSS" used NSamplehigh = 200, for
             * those fragments, the number of candidate can be over 100. ???
             * 2008-02-14, Nanjiang */

            int cntValidCan = 0;/*the number of fragment*/
            for (ik=0; ik<maxSample; ik++)
            {
                if (  High_Score[ik] < 0 ) { continue; }
                High_Score[cntValidCan] = High_Score[ik];
                Pos_Chain[cntValidCan]  = Pos_Chain[ik];
                Pos_length[cntValidCan] = Pos_length[ik];
                cntValidCan++;
            }
            numValidCan = cntValidCan;

            if ( NPer<50 &&  numValidCan > 100 )  /*if there are more than 100 candidates, higher similarity is required*/
            {
                cntValidCan = 0;
                for (ik=0; ik<numValidCan; ik++)
                {
                    nbeg = Pos_length[ik] + 2;
                    nend = Pos_length[ik] + 7 ;
                    if(!IsSameShape(digitShapeAll[iTar], digitShapeAll[Pos_Chain[ik]], iFragTar, Pos_length[ik], 2, 7, shapeCompareMode) )
                    {
                        continue;
                    }
                    High_Score[cntValidCan] = High_Score[ik];
                    Pos_Chain[cntValidCan]  = Pos_Chain[ik];
                    Pos_length[cntValidCan] = Pos_length[ik];
                    cntValidCan++;
                }
                numValidCan = cntValidCan;
            }

            /*sum up frequency*/
            for (ik=0; ik<numValidCan; ik++)
            {
                nbeg = Pos_length[ik];
                nend = Pos_length[ik] + FRAG_WINDOW_SIZE ;
                for (i= nbeg; i < nend; i++)
                {
                    fragProfile[i-nbeg][digitAASeqAll[Pos_Chain[ik]][i]]++;
                }
            }
#ifdef DEBUG_WRITE_FRAG /*2009-04-26*/ 
            /*print target fragment*/
            char aaFragCan[20] ="";
            char shapeFragCan[20] ="";
            fprintf(fpFragBlock,"%-9s %4s %10s %4d %s %s %s\n", "Target","cnt" , idTar, iFragTar,  aaFragTar, shapeFragTar, "Score");
            /*print all candidate fragments to this fragment*/
            for (ik=0; ik<numValidCan; ik++)
            {
                //nbeg = Pos_length[ik];
                //nend = Pos_length[ik] + FRAG_WINDOW_SIZE ;
                my_strcpy(aaFragCan, aaSeqAll[Pos_Chain[ik]]+Pos_length[ik], FRAG_WINDOW_SIZE); 
                my_strcpy(shapeFragCan, shapeSeqAll[Pos_Chain[ik]]+Pos_length[ik], FRAG_WINDOW_SIZE); 

                fprintf(fpFragBlock,"%-9s %4d %-10s %4d %s %s %6d\n", "Candidate",ik, idListAll[Pos_Chain[ik]], Pos_length[ik], aaFragCan, shapeFragCan, High_Score[ik] );
            }
            fprintf(fpFragBlock,"\n");
#endif
            //convert to percentage
            for (i=0; i<FRAG_WINDOW_SIZE; i++)
            {
                for (j=0; j<NUM_20_AA; j++)
                {
                    //fragProfile[i][j] = Integer( fragProfile[i][j]*100.0/numValidCan);
                    fragProfile[i][j] = int( fragProfile[i][j]*100.0/numValidCan+0.5);
                }
            }
            /*========== END calculating profiles ================*//*}}}*/
#ifdef DEBUG_FRAG
            /*print candidates*/
            fprintf(stdout,"frag of %s at %d has %d candidates\n", idTar,iFragTar , numValidCan);
            for(ik = 0; ik < numValidCan; ik ++)
            {
                fprintf(stdout,"candidate %d: %s %d score = %d\n", ik, idListAll[Pos_Chain[ik]], Pos_length [ik], High_Score[ik]);
            }
#endif

            //write to file
            if (  numValidCan <= 4  )  
            { 
#ifdef DEBUG_DBFRAG
                fprintf(stderr,"iFragTar=%d, numValidCan = %d <=4: frag=%s\n",iFragTar, numValidCan, shapeFragTar );
#endif
                continue; 
            }  /*neglect the position with <=4 candidates*/
            QuickSort(High_Score, 0, numValidCan-1);

            int scorePerPos = 0; /*score per position, which to some extent represents the quality of the fragacc profile at this position*/
            for (ik=0; ik<numValidCan; ik++)
            {
                scorePerPos += High_Score[ik];
            }
            //scorePerPos = Integer(scorePerPos/double(numValidCan));
            //scorePerPos = scorePerPos/numValidCan;
            scorePerPos = int(scorePerPos*1.0/numValidCan+0.5);
            for (ik=0; ik<FRAG_WINDOW_SIZE; ik++)
            {
                fprintf(fpFrag,"%4d %2d  ",iFragTar,ik);
                for (j=0; j<NUM_20_AA; j++)
                {
                    fprintf(fpFrag,"%3d ",fragProfile[ik][j]);
                }
                fprintf(fpFrag,"%3d %2d\n",numValidCan,scorePerPos);
            }
        }/*}}}*/
        fclose(fpFrag);   
#ifdef DEBUG_WRITE_FRAG /*2009-04-26*/
        fclose(fpFragBlock);
#endif

        Array1D <int> cntProfilePerPos_1darray(lengthSeqTar);
        Array1D <int> Number_chain_1darray(lengthSeqTar);
        Array1D <int> Score_all_1darray(lengthSeqTar);
        Array2D <int> Number_blosum_2darray(lengthSeqTar, NUM_20_AA);
        int  *cntProfilePerPos = cntProfilePerPos_1darray.array1D; /*the number of profiles for each residue position          */
        int  *Number_chain     = Number_chain_1darray.array1D;     /*sum of the number of candidates for each residue position */
        int  *Score_all        = Score_all_1darray.array1D;        /*sum of the similarity score for each residue position     */
        int **Number_blosum    = Number_blosum_2darray.array2D;    /*average profile for each residue position             */

        /*initialization*/
        for (i=0; i<lengthSeqTar; i++) /*i is usually used as the iterator for sequnece position or the general outest iterator*/
        {
            cntProfilePerPos[i] = 0;
            Number_chain[i]     = 0;
            Score_all[i]        = 0;
            for (j=0; j<20; j++) { Number_blosum[i][j] = 0; }
        }

        /*-----------read in fragFile --------------*/
        sprintf(fragFileName, "%s/%s.dbfrag",resultpath,idTar);
        fpFrag = fopen(fragFileName,"r");
        if( checkfilestream(fpFrag, fragFileName, "r") == -1) { continue; }
        int idxFrag      = 0; /*index of the 9-long fragment, [0,lengthSeqTar-9] */
        int posResInFrag = 0; /*index of the residue within each fragment, [0,8] */
        int aaSeqIndex   = 0; /*index of the residue, [0, lengthSeqTar-1]        */
        int NNPer[21]; /*storing the profile from fragFile*/
        int NSample = 0; /*number of candidates found for each residue position, read from dbfrag file*/
        aaSeqIndex = 0;
        while((linesize = fgetline(fpFrag,line, maxline))!= EOF  ) //read in Frag file/*{{{*/
        {
            status_sscanf= sscanf(line,"%d%d %d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d",
                    &idxFrag, &posResInFrag,&NNPer[0],&NNPer[1],&NNPer[2],&NNPer[3],&NNPer[4],&NNPer[5],&NNPer[6],&NNPer[7],&NNPer[8],&NNPer[9],
                    &NNPer[10],&NNPer[11],&NNPer[12],&NNPer[13],&NNPer[14],&NNPer[15],&NNPer[16],&NNPer[17],&NNPer[18],&NNPer[19],&NSample,&NScore);
            if(status_sscanf == 24)
            {
                aaSeqIndex = idxFrag+posResInFrag;
                for (j=0; j<20; j++)
                {
                    Number_blosum[aaSeqIndex][j] += NNPer[j]; /*add up the profile at the same residue position*/
                }
                cntProfilePerPos[aaSeqIndex]++; /*what does cntProfilePerPos mean here ??, 2007-11-03 */
                Number_chain[aaSeqIndex] += NSample; /*add up the number of candidates found for the fragment this residue belongs to*/
                Score_all[aaSeqIndex] = Score_all[aaSeqIndex] + NScore;
            }
            else
            {
                fprintf(stderr,"Error! dbfrag file %s is of wrong format at line \n%s", fragFileName, line);
            }
        }/*}}}*/
        fclose(fpFrag);

        for (i=0; i < lengthSeqTar; i++)/*get the average value*/ 
        {
            if (  cntProfilePerPos[i] <= 0 ) { continue; }
            for (j=0; j<NUM_20_AA; j++)
            {
                //Number_blosum[i][j] = Integer ( Number_blosum[i][j]/double(cntProfilePerPos[i]));
                Number_blosum[i][j] =  int(Number_blosum[i][j]*1.0/cntProfilePerPos[i]+0.5);
            }
            //Number_chain[i] = Integer(Number_chain[i]/double(cntProfilePerPos[i]));
            //Score_all[i] = Integer(Score_all[i]/double(cntProfilePerPos[i]));
            //Number_chain[i] = Number_chain[i]/cntProfilePerPos[i];
            //Score_all[i] = Score_all[i]/cntProfilePerPos[i];
            Number_chain[i] = int (Number_chain[i]/double(cntProfilePerPos[i])+0.5);
            Score_all[i] = int(Score_all[i]/double(cntProfilePerPos[i])+0.5);
        }

        /*write fragacc file*/
        char fragAccFileName[MAX_PATH+1] = "";
        if (isWriteBinaryFile)
        {   /*added 2009-11-08*/
            sprintf(fragAccFileName, "%s/%s.fragaccbin",result_acc_path,idTar); //added 2009-11-08 
            WriteBinaryFragAccFile(fragAccFileName, aaSeqAll[iTar], shapeSeqAll[iTar], waterAccAll[iTar], dsspSecAll[iTar], Number_blosum, Number_chain , Score_all, lengthSeqTar);

        }
        else
        {
            sprintf(fragAccFileName, "%s/%s.fragacc",result_acc_path,idTar); //frag_id.txt is the matrix file generated based on shape string and filtered by water accessibility, this is the final output, 2007-11-02, have the filename changed, using the extension to distinguish different file types

            FILE *fpFragAcc = fopen(fragAccFileName,"w");

            for(i = 0; i < lengthSeqTar; i ++)
            {
                if (  Number_chain[i] <= 5  ) { continue; }
                fprintf(fpFragAcc,"%5d %-2c %c%1d%c",i+1,aaSeqAll[iTar][i],shapeSeqAll[iTar][i], waterAccAll[iTar][i], dsspSecAll[iTar][i]);
                for (j=0; j<NUM_20_AA; j++) /*j is usually used as iterator on profile columns*/
                {
                    fprintf(fpFragAcc,"%4d",Number_blosum[i][j]);
                }
                fprintf(fpFragAcc,"  %4.2f %4.2f\n",Number_chain[i]/100.0,Score_all[i]/10.0);
            }
            fclose(fpFragAcc);
        }

        fprintf(stdout,"%d \t %s result output to %s\n",iTar, idTar, fragFileName);
    }/*}}}*/

    //free memory
    //delete [] NShapedis;
    for (i=0; i<maxNumIDTrain+1; i++)
    {
        if (aaSeqAll[i]           != NULL){delete [] aaSeqAll[i];          }
        if (shapeSeqAll[i]        != NULL){delete [] shapeSeqAll[i];       }
        if (dsspSecAll[i]         != NULL){delete [] dsspSecAll[i];        }
        if (digitAASeqAll[i]      != NULL){delete [] digitAASeqAll[i];     }
        if (digitShapeAll[i]      != NULL){delete [] digitShapeAll[i];     }
        if (waterAccAll[i]        != NULL){delete [] waterAccAll[i];       }
        if (encodeWaterAccAll[i]  != NULL){delete [] encodeWaterAccAll[i]; }
        if (idListAll[i]          != NULL){delete [] idListAll[i];         }
    }

    return 0;
}
/*}}}*/
int main(int argc, char** argv)/*{{{*/
{
    bool isNonOptionArg = false;
    if(argc < 2) 
    {
        PrintHelp();
        return 0;
    }
    int i,j;
    char submatfile[MAX_PATH+1] = "";
    char outfile[MAX_PATH+1] = "";
    double value = 0.0;
    const char control_option[] = ""; //options which control the program, and does not take parameters
    //bool isAll = false;
    //bool isQuiet = false;
    //bool isSingle = false;

    //cdatabase = "d:\\pred\\chains3829_list.txt";
    //idListFile = "I:\\pred\\chains_5246.txt";
    //qijpath = "D:\\Par_Qij3829";
    //qijpath = "D:\\Par_Qij5246";
    //shape_triple_file = "D:\\pred\\Shape_triple3829.txt";

    char idListFile[MAX_PATH+1] = "/misc/casiodata3/wk/passe/passe.idlist"; //file storing chain id list, e.g. 1UEOA 
    char qijpath[MAX_PATH+1] = "/misc/casiodata3/wk/passe/Par_Qijnanjiang"; //path for Qijfile, note that it should be in the format of Tuping,
    int qijformat = QIJ_FORMAT_TUPING;
    //note both the third column and the alphabet order is AAAlphabet_Tuping = AVLIPFMKRHGSTCYNEWDQ
    char shape_triple_file[MAX_PATH+1] = "/misc/casiodata3/wk/passe/Shape_triple3829.txt"; //frequency of triple shape string
    char resultpath[MAX_PATH+1] = "/misc/casiodata3/wk/passe/result_passe";        //intermediate result
    char result_acc_path[MAX_PATH+1] = "/misc/casiodata3/wk/passe/frag_acc_passe"; //result matrix file really needed

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
                        //case 'a': isAll = true; break;
                        //case 'q': isQuiet = true; break;
                        //case 's': isSingle = true; break;
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
            }
            else if( (strcmp(argv[i],"-o") == 0) || (strcmp(argv[i], "--out") == 0))  
            {
                if( ( i = option_parser_filename(argc, argv, i, outfile)) == -1)
                    return -1;
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
            else if (strcmp(argv[i], "--create-triple") == 0)
            {
                isCreateTriple = true;
                i ++;
            }
            else if( (strcmp(argv[i],"--matrix") == 0) )  
            {
                if( ( i = option_parser_filename(argc, argv, i, submatfile)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--triple") == 0) )  
            {
                if( ( i = option_parser_filename(argc, argv, i, shape_triple_file)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--qijpath") == 0) )  
            {
                if( ( i = option_parser_filename(argc, argv, i, qijpath)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--result") == 0) )  
            {
                if( ( i = option_parser_filename(argc, argv, i, resultpath)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--resultacc") == 0) )  
            {
                if( ( i = option_parser_filename(argc, argv, i, result_acc_path)) == -1)
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
            //else if (strcmp(argv[i], "--idtype") == 0)
            //{
            //    if( ( i = option_parser_numeric(argc, argv, i, idtype, true, 0, 1)) == -1)
            //        return -1;
            //}
            else if( (strcmp(argv[i], "--value") == 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, value, true, 0.0, 15.0)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--qijformat") == 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, qijformat, true, 0, 2)) == -1)
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
            my_strcpy(idListFile, argv[i],MAX_PATH);
            i ++;
        }
    }/*}}}*/

    VerifyFolder(resultpath);
    VerifyFolder(result_acc_path);


    Array2D <int> subMatrix_2darray(NUM_BLOSUM, NUM_BLOSUM);
    char alphabet[NUM_BLOSUM+1] = "";
    int **subMatrix  = subMatrix_2darray.array2D;
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

//     //debug
//     //print out matrix
//     for(i = 0 ; i  < 20; i++)
//     {
//         for(j = 0 ; j < 20; j++)
//         {
//             fprintf(stdout,"%4d ", subMatrix[i][j]);
//         }
//         fprintf(stdout,"\n");
//     }
//     fprintf(stdout,"\n");
// 
//     return 0;


    FILE *fpout = NULL;
    if(strcmp(outfile,"") == 0 || strcasecmp(outfile, "stdout") == 0)
    {
        fpout = stdout;
    }
    else
    {
        fpout = fopen(outfile, "w");
        checkfilestream(fpout, outfile,"w");
    }

    Build_Database(idListFile, qijpath, qijformat, shape_triple_file, resultpath, result_acc_path, subMatrix);


    if(fpout != NULL && fpout != stdout) fclose(fpout);

    return 0;
}
/*}}}*/
