/*
 * =====================================================================================
 *        Filename:  pssm2Qij.cpp
 *     Description:  convert weighted percentage matrix output by blastpgp -Q
 *                   flag to Qij, estimated percentages
 *         Version:  1.0
 *         Created:  09/01/2007 05:01:43 PM CET
 *        Compiler:  g++
 *          Author:  Nanjiang Shu (Shu), nanjiang@struc.su.se
 *         Company:  Structural Chemistry, Stockholm Univesity
 * =====================================================================================
 */

#include <cstdio>
#include <cmath>
#include <set>
#include <cstring>
#include "mytemplate.h"
#include "array.h"
#include "myfunc.h"
#include "mypro.h"

//ChangeLog/*{{{*/
/*****************************************************************************
 * ChangeLog: 2007-06-06  
 * The alpha (Nc-1) and beta (pseducounts) by default is not set to a constant
 * value any more. It will be derived from the weighted percentages and
 * pssm.score2.
 *
 * alpha = (number of non zero weighted percentages) / length of the sequence
 * beta  = max( (10 - score2/0.2) ,  0 )
 *
 * thus beta varies in each row
 *
 * formula for calculating beta can be adjusted
 *
 * ChangeLog: 2007-06-13
 * GetBkFreq checking if it is percentage or sum = 1
 * Add --ext so that the extension of qijfile can be set
 * Add -a    function, so that the alphabet order of the pssm can be set
 * Add -r    function, so that the default format for pssm is integer
 *           this is because Qij values here are derived weighted percentages,
 *           which is already rounded down to integer
 *
 * ChangeLog 2007-06-17
 * adding output shape string, water accessibility and dssp secondary
 * structure in tuping's format
 *
 * ChangeLog 2007-06-18
 * using different scale to map raw water accessibility from DSSP to 0~9 grade
 *
 *
 * ChangeLog 2007-10-22
 * add --idtype option, for idtype = 1, all parameter files are $id.seqmap
 *
 * ChangeLog 2007-11-01
 * output binary file also, reading binary file is about 100 times faster then
 * reading the ascii file
 *
 * ChangeLog 2008-01-15
 * add: normalize the Qij to have a sum of 100 on each position
 * add option: --not-hot: if --not-hot not set, treat all zero profiles. That
 * is, for all zero profiles, set the fij of the amino acid on this position
 * to 100
 *
 * ChangeLog 2008-07-06
 *   1. add option --pssmfilelist to solve the "argument list too long" error that
 *   might occur
 *   2. add the requirement of the seqmap file in the help message
 * ChangeLog 2008-07-07 
 *    add the option --seqmap to set the user defined seqmap path
 *    the argument parser changed to the standard one
 * ChangeLog 2009-07-27
 *    the function isTreatAllZeroProfile is updated
 * ChangeLog 2009-11-03
 *    add ProfileSADShort, so that the binary format matrix file can be even
 *    smaller
 ****************************************************************************/ /*}}}*/  

/*  BLOSUM 62 */
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
const double  bkFreq[] =/*{{{*/ // ordered by BLOSUM1D_alphabet ARNDCQEGHILKMFPSTWYV
{
    0.0780,  //  A
    0.0510,  //  R
    0.0450,  //  N
    0.0540,  //  D
    0.0190,  //  C
    0.0430,  //  Q
    0.0630,  //  E
    0.0740,  //  G
    0.0220,  //  H
    0.0520,  //  I
    0.0900,  //  L
    0.0570,  //  K
    0.0220,  //  M
    0.0390,  //  F
    0.0520,  //  P
    0.0710,  //  S
    0.0590,  //  T
    0.0130,  //  W
    0.0320,  //  Y
    0.0640   //  V
};
/*}}}*/
const int maxWaterAcc[] =/*{{{*/
{
    181,   //  A    
    326,   //  R    
    244,   //  N    
    243,   //  D    
    183,   //  C    
    267,   //  Q    
    286,   //  E    
    156,   //  G    
    254,   //  H    
    230,   //  I    
    245,   //  L    
    313,   //  K    
    260,   //  M    
    271,   //  F    
    204,   //  P    
    234,   //  S    
    217,   //  T    
    309,   //  W    
    304,   //  Y    
    216    //  V    
};/*}}}*/
const double waterAccScale[] = /*{{{*/
{
    0.0,  // 0
    0.02, // 1
    0.04, // 2
    0.06, // 3
    0.11, // 4
    0.16, // 5
    0.21, // 6  
    0.31, // 7
    0.41, // 8    
    1.00  // 9
};/*}}}*/

// const char AAAlphabet_Tuping[] = "AVLIPFMKRHGSTCYNEWDQ" ;
//const char *rescodes = "ARNDCQEGHILKMFPSTWYVBZX";
//const char *ncbicodes= "XAXCDEFGHIKLMNPQRSTVWXYXXX";

#undef NUM_PARAMETER
#define NUM_PARAMETER 8

#undef NUM_AA
#define NUM_AA 20

bool isQuietMode      = false;
bool isOutputInteger  = true;
bool isPseducountsSet = false;
bool isNcSet          = false;
bool isOutputSAD      = false; //whether to output shape string, water accessibility (0-9) and dssp secondary structure at the third column
bool isTreatAllZeroProfile = true; /*whether treat the all zero profile, if so, for all zero profiles, set the value of the current fij to 100*/
bool isNewSequence    = false; //whether the supplied pssm is derived from a new sequence whose structural information is not available

int8 typeProfile = 1; /*profile type, default = 1
                        0 - float ProfileSAD
                        1 - int8, ProfileSADByte
                        2 - float, Profile
                        3 - int8, ProfileByte*/

// scaling of water accessbility/*{{{*/
/*****************************************************************************
 * the acess surface is divided into 10 groups according to the percentage of the exposure surface
 * to water
 * 0 -- no       surface acccesible to water, takes up 14.1%;
 * 1 -- 0--2%    surface acccesible to water, takes up 10.9%;
 * 2 -- 2--4%    surface acccesible to water, takes up 6.4%;
 * 3 -- 4--6%    surface acccesible to water, takes up 4.9%;
 * 4 -- 6--11%   surface acccesible to water, takes up 9.6%;
 * 5 -- 11--16%  surface acccesible to water, takes up 8.4%;
 * 6 -- 16--21%  surface acccesible to water, takes up 7.7%;
 * 7 -- 21--31%  surface acccesible to water, takes up 13.9;
 * 8 -- 31--41%  surface acccesible to water, takes up 11.3%;
 * 9 -- 41--100% surface acccesible to water, takes up12.7%;
 *
 * while in search new, these 0~9 water accessbility are again mapped to five
 * grades, 
 * 0        --> 0
 * 1,2,3    --> 1
 * 4,5,6    --> 2
 * 7,8      --> 3
 * 9        --> 5
 ****************************************************************************/
/*}}}*/
void PrintHelp()
{
    fprintf(stdout,"Usage: pssm2Qij [OPTIONS] pssmfile\n");
    fprintf(stdout," Note that if --sad option is enabled, seqmap file should be available\n");
    fprintf(stdout,"OPTIONS: \n");
    fprintf(stdout,"  -list pssmfilelist    : set the list of pssm files\n");
    fprintf(stdout,"  -d outpath            : set the output path, default = ./\n");
    fprintf(stdout,"  --wb                  : output modm matrix file in binary format\n");
    fprintf(stdout,"  -c <real>             : pseducounts, if not set, calculated by max( (10 - score2/0.2),  0)\n");
    fprintf(stdout,"  -N <real>             : average number of independent residues in columns, if not set, calculated by \n");
    fprintf(stdout,"                        : (number of non zero weighted percentages) / length of the sequence + 1\n");
    fprintf(stdout,"  --maxpse <real>       : max_pseducount, default = 10\n");
    fprintf(stdout,"  --step   <real>       : step_score2, default = 0.2\n");
    fprintf(stdout,"                        : max_pseducount and step_score2 are used when pseducounts are not set, then\n");
    fprintf(stdout,"                        : pseducounts[i] = max (max_pseducount - score2[i]/step_score2, 0)\n");
    fprintf(stdout,"  -l  <real>            : lambda-std, natural scalue, default=0.3216, or read from the file\n");
    fprintf(stdout,"  -a  <string>          : supply the alphabet other than BLUSOM alphabet for outputing the pssm\n");
    fprintf(stdout,"  -r                    : whether output Qij in real format %%.2f, default is integer\n");
    fprintf(stdout,"  --bkfile file         : file storing background frequency\n");
    fprintf(stdout,"  --matrix file         : substitution matrix file, default is BLOSUM62\n");
    fprintf(stdout,"  --ext  string         : extension for qijfile, default = Qij\n");
    fprintf(stdout,"  --sad                 : output shape shtring, wateracc, dssp secondary structure, then the alphabet = AAAlphabet_Tuping\n");
    fprintf(stdout,"  --idtype 0 | 1        : 0 -- standardized id (5 chars), 1 -- exact id, default = 1\n");
    fprintf(stdout,"  --maxacc file         : file for maximal wateracc for each amino acid\n");
    fprintf(stdout,"  --newseq              : if new sequence, structural information is not available\n");
    fprintf(stdout,"  --accscale file       : file for scalue of wateracc\n");
    fprintf(stdout,"  --not-hot             : do not treat the all zero profiles\n");
    fprintf(stdout,"  --seqmap path         : set the seqmap path, default = $DATADIR/seqmap\n");
    fprintf(stdout,"  --typeprofile 0|1|2|3 : set the data type of profile, default = 1\n");
    fprintf(stdout,"                        : 0 and 1, with SAD,    0 -- for float type, 1 -- for short type\n");
    fprintf(stdout,"                        : 2 and 3, without SAD, 2 -- for float type, 3 -- for short type\n");
    fprintf(stdout,"  -q                    : quiet mode, do not output report\n");
    fprintf(stdout,"  -h|--help             : print this help message and exit\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"For Binary format file, the default extension: Qijbin\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"Created on 2007-01-09, updated 2009-11-08, Nanjiang Shu\n");
}
/*****************************************************************************
 * it is a smart way to use enum, in that case, ALA, ARG can be used in stead
 * of 0, 1, 2, which is more user readable
 ****************************************************************************/
enum aacodes // in blosum alphabet order/*{{{*/
{
    ALA, ARG, ASN, ASP, CYS,
    GLN, GLU, GLY, HIS, ILE,
    LEU, LYS, MET, PHE, PRO,
    SER, THR, TRP, TYR, VAL,
    UNK
};
/*}}}*/
//int GetBkFreq(const char *infile, double *P, const char *alphabet)[>{{{<]
/*****************************************************************************
 * read in background amino acid composition content from file, store them to
 * P according to the alphabet
 * 
 * the returned composition is unified, that is sum(P[i])  = 1.00;
 *
 *   format of the infile      
 *   ---------------------------------------
 *   # annotation line         
 *   A 0.078                   
 *   C 0.089                   
 * 
 * 2007-06-13
 ****************************************************************************/
//{
    //int i;
    //FILE *fpin;
    //fpin = fopen(infile, "r");
    //checkfilestream(fpin, infile, "r");

    //int linesize;
    //int maxline = 300;
    //Array1D <char> line_1darray(maxline+1);
    //char *line = line_1darray.array1D; 

    //i = 0 ; 
    //char ch;
    //double f;
    //int daa;
    //f_neglect_comment(fpin);
    //while((linesize = fgetline(fpin, line, maxline)) != EOF)
    //{
        //if(sscanf(line,"%c %lf", &ch, &f) == 2)
        //{
            //daa = Char2Digit(ch, alphabet);
            //if(daa >= 0)
            //{
                //P[daa] = f;
                //i ++;
            //}
        //}
    //}
    //fclose(fpin);
    //int numAA = i;
    //if( numAA != NUM_AA)
    //{
        //fprintf(stderr, "Warning! number of frequencies given in background frequency file %s is not equals to %d\n", infile, NUM_AA);
    //}
    
    //check if it is the percentage
    //double sum = 0.0; for(i = 0; i < numAA ; i ++) sum+= P[i];
    //if(fabs(sum - 100.0) < 1.0) // if is percentage, convert to frequency
    //{ for(i = 0 ; i < numAA; i ++) P[i]  /= 100.0; }
    //else if (fabs(sum-1.0) > 0.1)
    //{
        //fprintf(stderr, "Warning! the given background frequency is not unified in file %s\n", infile, NUM_AA);
    //}
    //return numAA;
//}
/*}}}*/
int GetMaxWaterAcc(const char *infile, int *acc, const char *alphabet)/*{{{*/
/*****************************************************************************
 * read in maximal water accessibility for each amino acid 
 * 
 *   format of the infile      
 *   ---------------------------------------
 *   # annotation line         
 *   A 182
 *   C 111
 * 
 * 2007-06-17
 ****************************************************************************/
{
    int i;
    FILE *fpin;
    fpin = fopen(infile, "r");
    checkfilestream(fpin, infile, "r");

    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D; 

    i = 0 ; 
    char ch;
    int a;
    int daa;
    f_neglect_comment(fpin);
    while((linesize = fgetline(fpin, line, maxline)) != EOF)
    {
        if(sscanf(line,"%c %d", &ch, &a) == 2)
        {
            daa = Char2Digit(ch, alphabet);
            if(daa >= 0)
            {
                acc[daa] = a;
                i ++;
            }
        }
    }
    fclose(fpin);
    int numAA = i;
    if( numAA != NUM_AA)
    {
        fprintf(stderr, "Warning! number of frequencies given in maxAcc file %s is not equals to %d\n", infile, NUM_AA);
    }
    
    return numAA;
}
/*}}}*/
int GetWaterAccScale(const char *infile, double *accScale)/*{{{*/
/*****************************************************************************
 * read in scale for mapping raw water accessibility to 0~9 grade
 * the scale should be ordered ascendingly
 * 2007-06-18
 ****************************************************************************/
{
    FILE *fpin;
    fpin = fopen(infile, "r");
    checkfilestream(fpin, infile, "r");

    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D; 

    int  cntBin = 0 ; 
    double d;
    f_neglect_comment(fpin);
    while((linesize = fgetline(fpin, line, maxline)) != EOF)
    {
        if(sscanf(line,"%lf", &d) == 1)
        {
            accScale[cntBin] = d;
            cntBin ++;
        }
    }
    fclose(fpin);
    return cntBin;
}
/*}}}*/
int ScaleWaterAcc(int rawWaterAcc, int maxAcc, double *accScale, int numBin)/*{{{*/
/*****************************************************************************
 * scalue raw water accessibility from dssp (which is the surface area of the
 * molecule exposed to water molecule in the unit of A**2) to 10 grade such
 * that the number of residues are generally equaly distributed in each grade
 *
 * rawWaterAcc : rawWaterAcc from dssp in the unit of A**2
 * maxAcc      : maximal rawWaterAcc for this amino acid
 * accScale    : scale for percentage of water accessibility, ascendingly ordered, e.g. 0.0 0.02 0.03 0.4
 ****************************************************************************/
{
    int i ; 
    double percentage = 0.0; //percentage = rawWaterAcc / maxAcc
    percentage = rawWaterAcc / double(maxAcc+1);
    for(i = 0 ; i < numBin; i ++)
    {
        if(percentage <= waterAccScale[i])
            return i;
    }
    return 0;
}/*}}}*/
void normalizeVector(double *v, int size, double normScale = 100.0)/*{{{*/
{
    int i;
    double sum = 0;
    for(i=0;i<size;i++) { sum += v[i]; }
    if(fabs(sum)>1e-6) { for(i=0;i<size;i++) { v[i] = v[i]/sum * normScale; } }
}
/*}}}*/
               
int GetNumIndependentRes(int **fij,int length, double *alpha)/*{{{*/
/*****************************************************************************
 * Get the independent number of residues in each column from weighted
 * percentages
 ****************************************************************************/
{
    int i,j ;
    int cnt = 0; // count of the non zero weighted percentages
    for(i = 0 ; i < length; i ++)
    {
        for (j = 0 ; j < NUM_AA; j ++)
        {
            if(fij[i][j] != 0) cnt ++;
        }
    }

    double a = cnt / double(length);
    for(i = 0 ; i < length; i ++) alpha[i] = a;
    return length;
}
/*}}}*/
int GetPseducount(double *score2, int length, double *beta, double max_pseducount = 10.0, double step_score2 = 0.2)/*{{{*/
/*****************************************************************************
 * Get pseducounts in each column from the score2 in PSSM matrix 
 * with the formula
 * pseducounts[i] = max(max_pseducount - score2[i] / step_score2, 0)
 ****************************************************************************/
{
    int i ;
    for(i = 0 ; i < length; i ++) 
    {
        beta[i] = max(max_pseducount - score2[i]/step_score2, 0.0);
    }
    return length;
}
/*}}}*/
int fij2Qij(int **fij, int length, double* alpha, double *beta, double lambda_std, double *P, double **Sij, double **Qij)/*{{{*/
/*****************************************************************************
 * fij        : weighted percentages
 * length     : length of the sequnece
 * alpha      : number of independent residues for each column, they can be different for each column 
 * beta       : pseducounts for each column
 * lambda_std : standard lambda, 
 * Pi         : background amino acid composition
 * Sij        : substitution matrix, e.g. BLOSUM62
 * formula    : Qkj = (alpha * fkj + beta * Pj * sum (fki * exp(lambda *Sij))) / (alpha + beta)
 *
 * Qij        : output value
 ****************************************************************************/
{
    int i,j,k;
    Array2D <double> eij_2darray(NUM_BLOSUM, NUM_BLOSUM);
    double ** eij = eij_2darray.array2D; // eij = exp(lambda_std * Sij)
    //calculating eij
    for(i = 0 ; i < NUM_AA ; i ++)
    {
        for(j = 0 ; j < NUM_AA; j++)
            eij[i][j] = exp(lambda_std * double (Sij[i][j]));
    }

    for(k = 0 ;  k < length ; k ++) //k: iterator for row
    {
        for(j = 0 ; j < NUM_AA; j ++)
        {
            double sum = 0 ; 
            for(i = 0 ; i < NUM_AA; i++)
                sum += fij[k][i]*eij[i][j];

            Qij[k][j] = (alpha[k]* fij[k][j] + beta[k] * P[j] * sum) /(alpha[k] + beta[k]);
        }
        /*normalize Qij profile to 100.0*/
        normalizeVector(Qij[k], NUM_AA, 100.0); /*added 2008-01-15, Nanjiang*/
    }
    return 0;
}
/*}}}*/
void WritePSSM(const char *outfile, char *aaSeq, int length, double **pssm, double* parameter, double *score1, double *score2, int* pssmIndex, char *outAlphabet, bool isOutputInteger, double *alpha, double *beta, double max_pseducount, double step_score2, Chain *pSeqmapChain)/*{{{*/
/*****************************************************************************
 * Write out the pssm
 ****************************************************************************/
{
    FILE *fpout;
    fpout = fopen(outfile, "w");
    checkfilestream(fpout, outfile, "w");
    int i,j;
    fprintf(fpout,"#Num: amino acid sequence number\n");
    fprintf(fpout,"#AA: 1 letter amino acid\n");
    fprintf(fpout,"# Qij, i.e. estimated frequencies are converted from fij (observed frequency)\n");

    fprintf(fpout,"# alpha (number of independent residues in column) = %.1lf\n", alpha[0]);

    if(isPseducountsSet)
        fprintf(fpout,"# beta (pseducounts) = %.1lf\n", beta[0]);
    else
        fprintf(fpout,"# beta (pseducounts) is estimated by max(%.1lf - score2[i]/%.1lf, 0)\n", max_pseducount, step_score2) ;

    if(isOutputSAD)
    {
        fprintf(fpout,"# SAD : shape string, water accessibility, secondary structure in dssp\n");
    }

    fprintf(fpout,"# \n");

    fprintf(fpout,"%-5s %2s","Num","AA");
    if(isOutputSAD)
    {
        fprintf(fpout," %3s", "SAD");
    }
    if(isOutputInteger)
        for(j = 0 ; j < NUM_AA ; j++) fprintf(fpout,"%4c",outAlphabet[j]);
    else
        for(j = 0 ; j < NUM_AA ; j++) fprintf(fpout," %5c",outAlphabet[j]);

    fprintf(fpout, "\n");
    for(i = 0 ; i < length ; i ++)
    {
        fprintf(fpout,"%5d %c",i+1,aaSeq[i]);
        if(isOutputSAD)
        {
            fprintf(fpout, " %c%1d%c", pSeqmapChain->shString[i], pSeqmapChain->waterAcc[i], pSeqmapChain->secStruc[i]);
        }
        if(isOutputInteger)
            for(j = 0 ; j < NUM_AA ; j++) fprintf(fpout,"%4d",Integer(pssm[i][pssmIndex[j]]));
        else
            for(j = 0 ; j < NUM_AA ; j++) fprintf(fpout," %5.2lf",pssm[i][pssmIndex[j]]);
        if(score1 != NULL)
            fprintf(fpout,"  %4.2f", score1[i]);
        if(score2!= NULL)
            fprintf(fpout," %4.2f", score2[i]);
        fprintf(fpout,"\n");
    }
    fprintf(fpout,"PSI Ungapped         %.4lf      %.4lf \n", parameter[4], parameter[5]   );
    fprintf(fpout,"PSI Gapped           %.4lf      %.4lf \n", parameter[6], parameter[7]   );

    fclose(fpout);
}
/*}}}*/
void WriteBinaryPSSMFile(const char *outfile, char *aaSeq, int length, double **pssm, double* parameter, double *score1, double *score2, int* pssmIndex, char *outAlphabet, Chain *pSeqmapChain)/*{{{*/
/*****************************************************************************
 * Write out the pssm
 ****************************************************************************/
{
    int i, j ;
    if (isOutputSAD)
    {
        if (typeProfile == 0)
        {   /*profile format: float*/
            Array1D <ProfileSAD> profile_1darray(length);
            ProfileSAD *profile = profile_1darray.array1D;
            int sizeAlphabet = strlen(outAlphabet);
            for(i = 0; i < length; i ++)
            {
                profile[i].aaSeqIndex = i+1;
                profile[i].aa = aaSeq[i];
                profile[i].shape = pSeqmapChain->shString[i];
                profile[i].waterAcc = pSeqmapChain->waterAcc[i];
                profile[i].dsspSec = pSeqmapChain->secStruc[i];
                for(j = 0 ; j < sizeAlphabet; j ++)
                {
                    profile[i].p[j] = float(pssm[i][pssmIndex[j]]);
                }
                if(score1 != NULL) { profile[i].score1 = float(score1[i]); }
                else { profile[i].score1 = 0.0; }
                if(score2 != NULL) { profile[i].score2 = float(score2[i]); }
                else { profile[i].score2 = 0.0; }
            }
                WriteBinaryMODM(outfile, outAlphabet, length, profile, parameter, typeProfile);
        }
        else
        {   /*profile format: int8*/ 
            Array1D <ProfileSADByte> profileSADByte_1darray(length);
            ProfileSADByte *profileSADByte = profileSADByte_1darray.array1D;

            int sizeAlphabet = strlen(outAlphabet);
            for(i = 0; i < length; i ++)
            {
                profileSADByte[i].aaSeqIndex = short(i+1);
                profileSADByte[i].aa = aaSeq[i];
                profileSADByte[i].shape = pSeqmapChain->shString[i];
                profileSADByte[i].waterAcc = int8(pSeqmapChain->waterAcc[i]);
                profileSADByte[i].dsspSec = pSeqmapChain->secStruc[i];
                for(j = 0 ; j < sizeAlphabet; j ++)
                {
                    profileSADByte[i].p[j] = Integer(pssm[i][pssmIndex[j]]);
                }
                if(score1 != NULL) { profileSADByte[i].score1 = float(score1[i]); }
                else { profileSADByte[i].score1 = 0.0; }
                if(score2 != NULL) { profileSADByte[i].score2 = float(score2[i]); }
                else { profileSADByte[i].score2 = 0.0; }
            }
#ifdef DEBUG_BINPROFILE
            for (i=0;i<length;i++)
            {
                fprintf(stdout,"%4d %-2c",profileSADByte[i].aaSeqIndex,profileSADByte[i].aa);
                fprintf(stdout, " %c%1d%c", profileSADByte[i].shape, profileSADByte[i].waterAcc, profileSADByte[i].dsspSec);
                for(j = 0 ; j < 20 ; j++) fprintf(stdout,"%4d",profileSADByte[i].p[j]);
                fprintf(stdout,"  %4.2f", profileSADByte[i].score1);
                fprintf(stdout," %4.2f", profileSADByte[i].score2);
                fprintf(stdout,"\n");
            }
#endif
            WriteBinaryMODM(outfile, outAlphabet, length, profileSADByte, parameter, typeProfile);
        }
    }
    else
    {
        if (typeProfile == 2)
        {
            Array1D <Profile> profile_1darray(length);
            Profile *profile = profile_1darray.array1D;
            int sizeAlphabet = strlen(outAlphabet);
            for(i = 0; i < length; i ++)
            {
                profile[i].aaSeqIndex = i+1;
                profile[i].aa = aaSeq[i];
                for(j = 0 ; j < sizeAlphabet; j ++)
                {
                    profile[i].p[j] = float(pssm[i][pssmIndex[j]]);
                }
                if(score1 != NULL) { profile[i].score1 = float(score1[i]); }
                else { profile[i].score1 = 0.0; }
                if(score2 != NULL) { profile[i].score2 = float(score2[i]); }
                else { profile[i].score2 = 0.0; }
            }
            WriteBinaryMODM(outfile, outAlphabet, length, profile, parameter, typeProfile);
        }
        else
        {
            Array1D <ProfileByte> profile_1darray(length);
            ProfileByte *profile = profile_1darray.array1D;
            int sizeAlphabet = strlen(outAlphabet);
            for(i = 0; i < length; i ++)
            {
                profile[i].aaSeqIndex = short (i+1);
                profile[i].aa = aaSeq[i];
                for(j = 0 ; j < sizeAlphabet; j ++)
                {
                    profile[i].p[j] = Integer(pssm[i][pssmIndex[j]]);
                }
                if(score1 != NULL) { profile[i].score1 = float(score1[i]); }
                else { profile[i].score1 = 0.0; }
                if(score2 != NULL) { profile[i].score2 = float(score2[i]); }
                else { profile[i].score2 = 0.0; }
            }
            WriteBinaryMODM(outfile, outAlphabet, length, profile, parameter, typeProfile);
        }
    }
}
/*}}}*/
bool IsValidAlphabet(char *alphabet)/*{{{*/
/*****************************************************************************
 * check if the alphabet is valid amino acid alphabet
 ****************************************************************************/
{
    int i;
    char stdAlphabet[] = "ARNDCQEGHILKMFPSTWYV";
    if(strlen(alphabet) < NUM_AA)
        return false;
    else 
    {
        for(i =0 ; i < NUM_AA; i ++)
        {
            if(!IsInCharSet(alphabet[i], stdAlphabet))
                return false;
        }
    }
    return true;
}/*}}}*/

int main(int argc, char** argv)/*{{{*/
{
    if (argc < 2)
    {
        PrintHelp();
        return -1;
    }

    //double scale_wateracc = 10; // scale water accessibility to 0-9
    int maxAcc[NUM_BLOSUM] ;
    
    int i = 0;
    char outpath[MAX_PATH+1] = "./";
    char qijfileext[MAX_PATH+1] = "Qij"; //default qijfile extension is Qij
    char outAlphabet[NUM_BLOSUM+1] = ""; //if outAlphabet not set, using the default BLOSUM alphabet
    double pseducounts = 9.0;
    double Nc  = 8.0;

    int idtype = 1; /*default idtype = 1, changed 2009-07-15*/

    double max_pseducount = 10.0; // max_pseducount and step_score2 are set for the formula to estimate pseducounts in each row with the formula
    double step_score2 = 0.2;     // pseducounts = max( max_pseducount - score2[i]/step_score2, 0)
    
    double lambda_std = 0.3216;
    char bkfreqfile[MAX_PATH+1] = "";
    char maxaccfile[MAX_PATH+1] = ""; // file list the maximal water accessibility for each amino acid
    char accscalefile[MAX_PATH+1] = "";// file list the scale for mapping water accessibility to 0~9
    char submatfile[MAX_PATH+1] = "";
    bool isOutputBinaryFile = false;

    char pssmFileList[MAX_PATH+1] =""; /*added 2008-07-06*/

    set <string> filenamelist_set;
    set <string> ::iterator iss;
    int errmsg = 0;

    char seqmappath[MAX_PATH+1] = "";
    char datadir[MAX_PATH+1] = "";
    GetDataDir(datadir);
    sprintf(seqmappath, "%s/%s", datadir, "seqmap");

    bool isNonOptionArg = false;

    i = 1;
    while(i < argc)/*{{{*/
    {
        if(argv[i][0] == '-' && !isNonOptionArg) //options
        {
            if(strcmp(argv[i], "-h") == 0|| strcmp(argv[i], "--help") == 0)
            {
                PrintHelp();
                errmsg = 1;
                break;
            }
            else if (strcmp(argv[i], "-d") == 0)
            {
                if( ( i = option_parser_filename(argc, argv, i, outpath)) == -1)
                    return -1;
            }
            else if (strcmp(argv[i], "--seqmap") == 0)
            {
                if( ( i = option_parser_filename(argc, argv, i, seqmappath)) == -1)
                    return -1;
            }
            else if (strcmp(argv[i], "-list") == 0)
            {
                if( ( i = option_parser_filename(argc, argv, i, pssmFileList)) == -1)
                    return -1;
            }
            else if (strcmp(argv[i], "-c") == 0)
            {
                if( ( i = option_parser_numeric(argc, argv, i, pseducounts, true, 0.0, 1e3)) == -1)
                    return -1;
                isPseducountsSet = true;
            }
            else if (strcmp(argv[i], "--newseq") == 0)
            {
                isNewSequence = true;
                i ++;
            }
            else if (strcmp(argv[i], "--maxpse") == 0)
            {
                if( ( i = option_parser_numeric(argc, argv, i, max_pseducount, true, 0.0, 1e3)) == -1)
                    return -1;
            }
            else if (strcmp(argv[i], "--typeprofile") == 0)
            {
                if( ( i = option_parser_numeric(argc, argv, i, typeProfile, true, int8(0), int8(3))) == -1)
                    return -1;
            }
            else if (strcmp(argv[i], "--wb") == 0)
            {
                isOutputBinaryFile = true;
                i ++;
            }
            else if (strcmp(argv[i], "--not-hot") == 0)
            {
                isTreatAllZeroProfile = false;
                i ++;
            }
            else if (strcmp(argv[i], "--step") == 0)
            {
                if( ( i = option_parser_numeric(argc, argv, i, step_score2, true, 0.0, 1e2)) == -1)
                    return -1;
            }
            else if (strcmp(argv[i], "-l") == 0)
            {
                if( ( i = option_parser_numeric(argc, argv, i, lambda_std, true, 0.0, 1e2)) == -1)
                    return -1;
            }
            else if (strcmp(argv[i], "--idtype") == 0)
            {
                if( ( i = option_parser_numeric(argc, argv, i, idtype, true, 0, 1)) == -1)
                    return -1;
            }
            else if (strcmp(argv[i], "-N") == 0)
            {
                if( ( i = option_parser_numeric(argc, argv, i, Nc, true, 0.0, 1e3)) == -1)
                    return -1;
                isNcSet = true;
            }
            else if (strcmp(argv[i], "--sad") == 0)
            {
                isOutputSAD = true;
                i ++;
            }
            else if (strcmp(argv[i], "--alphabet") == 0 || strcmp(argv[i], "-a") == 0)
            {
                my_strcpy(outAlphabet,argv[i+1], NUM_BLOSUM);
                if(!IsValidAlphabet(outAlphabet))
                {
                    errmsg = -1;
                    break;
                }
                i += 2;
            }
            else if (strcmp(argv[i], "--bkfile") == 0)
            {
                if( ( i = option_parser_filename(argc, argv, i, bkfreqfile)) == -1)
                    return -1;
            }
            else if (strcmp(argv[i], "--maxacc") == 0)
            {
                if( ( i = option_parser_filename(argc, argv, i, maxaccfile)) == -1)
                    return -1;
            }
            else if (strcmp(argv[i], "--accscale") == 0)
            {
                if( ( i = option_parser_filename(argc, argv, i, accscalefile)) == -1)
                    return -1;
            }
            else if (strcmp(argv[i], "--ext") == 0)
            {
                if( ( i = option_parser_filename(argc, argv, i, qijfileext)) == -1)
                    return -1;
            }
            else if (strcmp(argv[i], "--matrix") == 0)
            {
                if( ( i = option_parser_filename(argc, argv, i, submatfile)) == -1)
                    return -1;
            }
            else if (strcmp(argv[i], "-r") == 0)
            {
                isOutputInteger = false;
                i ++;
            }
            else if (strcmp(argv[i], "-q") == 0)
            {
                isQuietMode = true;
                i ++;
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
        else
        {
            filenamelist_set.insert(argv[i]);
            i ++;
        }
    }/*}}}*/


    VerifyFolder(outpath);
    if(errmsg != 0)
    {
        return errmsg;
    }

    /*check consistency, 2009-11-05 */
    if (isOutputSAD && typeProfile >= 2)
    {
        fprintf(stderr,"Erro! the typeProfile should be set to 0 or 1 when SAD is enabled.\n");
        return -1;
    }
    if ((! isOutputSAD)&& typeProfile < 2)
    {
        fprintf(stderr,"Erro! the typeProfile should be set to 2 or 3 when SAD is disabled.\n");
        return -1;
    }

    if (strcmp(pssmFileList, "") != 0)
    {
        int linesize;
        int maxline = MAX_PATH+1;
        Array1D <char> line_1darray(maxline+1);
        char *line = line_1darray.array1D; 
        FILE *fplist = fopen(pssmFileList,"r");
        if(checkfilestream(fplist, pssmFileList,"r") != -1 )
        {
            while((linesize = fgetline(fplist, line, maxline)) != EOF)
            {
                filenamelist_set.insert(line);
            }
        }
        else
        {
            return -1;
        }
        fclose(fplist);
    }

    //int numfile = filenamelist_set.size();

    Array2D <double> Sij_2darray(NUM_BLOSUM, NUM_BLOSUM);
    Array1D <double > P_1darray(NUM_BLOSUM);
    double **Sij  = Sij_2darray.array2D;
    double *P = P_1darray.array1D;
    char alphabet[NUM_BLOSUM+1] = "";
    int j;

    if(strcmp(submatfile,"") == 0)
    {
        for(i = 0 ; i < NUM_BLOSUM; i++)
            for(j = 0 ; j  < NUM_BLOSUM ; j++)
                Sij[i][j] = blosum62[i][j];
    }
    else { GetPatternTransMatrix(submatfile, Sij, alphabet); }

    if(strcmp(bkfreqfile,"")==0) { for(i = 0 ; i < NUM_AA; i++) P[i] = bkFreq[i]; }
    else 
    { 
        GetBkFreq(bkfreqfile,P, BLOSUM1D_alphabet, NUM_AA); 
    }

#ifdef DEBUG_CONST_CHAR
    printf("BLOSUM1D_alphabet=%s\n", BLOSUM1D_alphabet);
#endif

    if(strcmp(maxaccfile,"")==0) { for(i = 0 ; i < NUM_AA; i++) maxAcc[i] = maxWaterAcc[i]; }
    else { GetMaxWaterAcc(maxaccfile,maxAcc, BLOSUM1D_alphabet); }

    double *accScale = NULL;
    int numScale = sizeof(waterAccScale)/sizeof(double);
    if(strcmp(accscalefile,"")==0)
    {
        numScale = sizeof(waterAccScale)/sizeof(double);
        accScale = new double [numScale];
        for(i = 0 ; i < numScale; i++) accScale[i] = waterAccScale[i];
    }
    else
    {
        int linecnt = fgetlinecnt(accscalefile);
        accScale = new double [linecnt];
        numScale = GetWaterAccScale(accscalefile,accScale);
    }

    Array1D <char> aaSeq_1darray(LONGEST_SEQ+1);
    Array2D <int> fij_2darray(LONGEST_SEQ, NUM_BLOSUM);
    Array2D <double> Qij_2darray(LONGEST_SEQ, NUM_BLOSUM);
    Array1D <double> score1_1darray(LONGEST_SEQ);
    Array1D <double> score2_1darray(LONGEST_SEQ);
    Array1D <double> parameter_1darray(NUM_PARAMETER);
    char *aaSeq = aaSeq_1darray.array1D;
    int **fij = fij_2darray.array2D;
    double **Qij = Qij_2darray.array2D;
    double *score1 = score1_1darray.array1D;
    double *score2 = score2_1darray.array1D;
    double *parameter = parameter_1darray.array1D;

    Array1D <double> alpha_1darray(LONGEST_SEQ);
    Array1D <double> beta_1darray(LONGEST_SEQ);
    double *alpha = alpha_1darray.array1D;       //number of independent residues in each column
    double *beta = beta_1darray.array1D;         //pseducounts in each column


    char seqmapfile[MAX_PATH+1] = "";
    char id[SIZE_CHAIN_ID+1] = "";

    Chain seqmap_chain;
    InitChain(&seqmap_chain);
    if(isOutputSAD)
    {
        seqmap_chain.aaSeq  = new char[LONGEST_SEQ+1] ;
        seqmap_chain.waterAcc  = new int[LONGEST_SEQ+1] ;
        seqmap_chain.shString  = new char[LONGEST_SEQ+1] ;
        seqmap_chain.secStruc  = new char[LONGEST_SEQ+1] ;
    }

    char pssmfile[MAX_PATH+1] = "";
    char Qijfile[MAX_PATH+1] = "";
    char rtname[MAX_PATH+1] = "";
    int length;
    int cntfile = 0;
    for(iss = filenamelist_set.begin() ; iss != filenamelist_set.end(); iss ++)
    {
        for(j = 0  ; j < NUM_PARAMETER; j ++) parameter[j] = 0.0;
        my_strcpy(pssmfile,(*iss).c_str(), MAX_PATH);
        if(GetPSSM(pssmfile, length, aaSeq, NULL, fij, score1, score2, parameter) <= 0)
        {
            fprintf(stderr, "Getting pssm file error!  pssm file %s", pssmfile);
            break;
        }

        if(isTreatAllZeroProfile)
        {
            TreatAllZeroFij(fij, length, score1, aaSeq, BLOSUM1D_alphabet);
#ifdef DEBUG
            fprintf(stdout, "TreatAllZeroFij is running\n");
#endif
        }/*Nanjiang, 2008-01-15*/

        if(parameter[1] != 0.00)
            lambda_std = parameter[1];

        if(!isNcSet) 
            GetNumIndependentRes(fij,length, alpha);
        else
        {
            for(j = 0 ; j < length; j ++) alpha[j] = Nc -1.0;
        }

        if(!isPseducountsSet)
            GetPseducount(score2, length, beta, max_pseducount, step_score2);
        else
        {
            for(j = 0 ; j < length; j ++) beta[j] = pseducounts;
        }

        fij2Qij(fij, length, alpha, beta, lambda_std, P,Sij, Qij);

        rootname(pssmfile, rtname);
        sprintf(Qijfile, "%s/%s.%s", outpath, rtname, qijfileext);
        int pssmIndex[NUM_BLOSUM];   /*pssmIndex is for reordering the profile in column according the given alphabet*/
        for(j = 0; j < NUM_BLOSUM; j ++) pssmIndex[j] = j;
        if(strcmp(outAlphabet,"") != 0)
        {
            for(i = 0 ; i < NUM_AA; i++)
            {
                pssmIndex[i] = Char2Digit(outAlphabet[i],  BLOSUM1D_alphabet);
            } 
        }
        else
        {
            my_strcpy(outAlphabet, BLOSUM1D_alphabet, NUM_AA);
        }

        if(isOutputSAD)
        {
			if (isNewSequence)
			{
				seqmap_chain.numRes = length ;
				for(i = 0 ; i < seqmap_chain.numRes; i ++)
				{
					seqmap_chain.secStruc[i] = '-';
					seqmap_chain.waterAcc[i] = 0;
					seqmap_chain.shString[i] = '-';
				}
				seqmap_chain.secStruc[seqmap_chain.numRes] = '\0';
				seqmap_chain.shString[seqmap_chain.numRes] = '\0';
			}
            else
            {
                if (idtype == 0)
                {
                    my_strcpy(id,rtname, SIZE_CHAIN_ID);
                    StdID(id);
                    GetSEQMAPFilePath(id, seqmapfile , seqmappath); /*supply seqmappath, 2009-07-15*/
                }
                else if (idtype == 1)
                {
                    my_strcpy(id,rtname, SIZE_CHAIN_ID);
                    sprintf(seqmapfile, "%s/%s.seqmap", seqmappath, id);
                }   /*changed 2007-10-22 by Nanjiang*/


                seqmap_chain.numRes = GetSEQMAP(seqmapfile, &seqmap_chain);

                if(seqmap_chain.numRes <= 0)
                { 
                    fprintf(stderr,"Error! seqmapfile %s can not open or is empty\n", seqmapfile);
                    continue; 
                }

                for(i = 0 ; i < seqmap_chain.numRes; i ++)
                {
                    if(seqmap_chain.secStruc[i] == ' ')
                    {
                        seqmap_chain.secStruc[i] = DSSP_SEC_RANDOM;
                    }
                    int daa = Char2Digit(seqmap_chain.aaSeq[i], BLOSUM1D_alphabet);
                    if(daa >= 0 && daa < NUM_AA)
                    {
                        //                        seqmap_chain.waterAcc[i] = max( Integer (seqmap_chain.waterAcc[i] * scale_wateracc / maxAcc[daa]+1), 0);
                        seqmap_chain.waterAcc[i] = ScaleWaterAcc(seqmap_chain.waterAcc[i], maxAcc[daa], accScale, numScale);
                        if(seqmap_chain.waterAcc[i] > 9 && seqmap_chain.waterAcc[i] < 0)
                        {
                            fprintf(stderr, "Error! seqmap_chain.WaterAcc[%d] = %d, not within [0-9]\n", i, seqmap_chain.waterAcc[i]);
                        }
                    }
                    else
                    {
                        seqmap_chain.waterAcc[i] = 0; // set waterAcc to 0 for non standard amino acid
                    }
                }
            }
        }

        if (isOutputBinaryFile)
        {
            sprintf(Qijfile, "%s%s", Qijfile, "bin");
            WriteBinaryPSSMFile(Qijfile, aaSeq, length, Qij, parameter, score1, score2, pssmIndex, outAlphabet, &seqmap_chain);
        }
        else
        {
            WritePSSM(Qijfile, aaSeq, length, Qij, parameter, score1, score2, pssmIndex, outAlphabet, isOutputInteger, alpha, beta, max_pseducount, step_score2, &seqmap_chain);
        }

        cntfile ++;
        if(!isQuietMode) 
        { fprintf(stdout,"%d \t %s output\n", cntfile, Qijfile); }
    }

    if(isOutputSAD)
    {
        DeleteChain(&seqmap_chain);
    }
    if (accScale != NULL){ delete [] accScale; }
    //
    return 0;
}/*}}}*/

