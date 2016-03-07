/*
 * =====================================================================================
 * 
 *        Filename:  pssm2modm.cpp
 *     Description:  convert get modmatrix from pssm file, modifed for tupings
 *     modmmatrix format
 *         Version:  1.0
 *         Created:  09/01/2007 05:01:43 PM CET
 *        Revision:  none
 *        Compiler:  g++
 *          Author:  Nanjiang Shu (Shu), nanjiang@struc.su.se
 *         Company:  Structural Chemistry, Stockholm Univesity
 * =====================================================================================
 */

#include <cstdio>
#include <cmath>
#include <set>
#include <cstring>
#include <string>
#include "mytemplate.h"
#include "array.h"
#include "myfunc.h"
#include "mypro.h"

// const char AAAlphabet_Tuping[] = "AVLIPFMKRHGSTCYNEWDQ" ;
//const char *rescodes = "ARNDCQEGHILKMFPSTWYVBZX";
//const char *ncbicodes= "XAXCDEFGHIKLMNPQRSTVWXYXXX";

#undef NUM_PARAMETER
#define NUM_PARAMETER 8

#undef NUM_FILE
#define NUM_FILE 10000

#undef NUM_AA
#define NUM_AA 20

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

bool isQuietMode      = false;
bool isOutputSAD      = false; //whether to output shape string, water accessibility (0-9) and dssp secondary structure at the third column
bool isTreatAllZeroProfile = true; /*by default treat the all zero weighted matrix profile as Hot-spot, that is set "100" on this amino acid, or all 5 in case the residue is 'X', 2009-07-27*/
bool isNewSequence    = false; //whether the supplied pssm is derived from a new sequence whose structural information is not available

int8 typeProfile = 1; /*profile type, default = 1
                        0 - float ProfileSAD
                        1 - int8, ProfileSADByte
                        2 - float, Profile
                        3 - int8, ProfileByte*/

//ChangeLog/*{{{*/
/*****************************************************************************
 * ChangeLog 2007-06-18
 *   using different scale to map raw water accessibility from DSSP to 0~9 grade
 *   ChangeLog 2007-10-22
 *   --idtype option added
 *
 * ChangeLog 2007-11-01
 *   output binary file also, reading binary file is about 100 times faster then
 *   reading the ascii file
 * ChangeLog 2009-07-15 
 *      option: --seqmap added
 *      option: -list added
 *      the argument parsing changed 
 * ChangeLog 2009-07-27
 *      TreatAllZeroFij added for weighted percentages
 ****************************************************************************/
/*}}}*/
void PrintHelp()
{
    fprintf(stdout,"Usage: pssm2modm [options] pssmfile\n");
    fprintf(stdout," Note that if --sad option is enabled, seqmap file should be available\n");
    fprintf(stdout,"Options: \n");
    fprintf(stdout,"  -list pssmfilelist    : set the list of pssm files\n");
    fprintf(stdout,"  -t log|per            : log or percentages, default = log\n");
    fprintf(stdout,"  -d outpath            : set the output path, default = ./\n");
    fprintf(stdout,"  -b                    : output modm matrix file in binary format\n");
    fprintf(stdout,"  -a  <string>          : supply the alphabet other than BLUSOM alphabet for outputing the pssm\n");
    fprintf(stdout,"  --sad                 : output shape shtring, wateracc, dssp secondary structure\n");
    fprintf(stdout,"                        : then the alphabet = AAAlphabet_Tuping\n");
    fprintf(stdout,"  --maxacc file         : file for maximal wateracc for each amino acid\n");
    fprintf(stdout,"  --newseq              : if new sequence, structural information is not available\n");
    fprintf(stdout,"  --accscale file       : file for scalue of wateracc\n");
    fprintf(stdout,"  --ext <string>        : extension for modmfile, default = modm\n");
    fprintf(stdout,"  --idtype 0 | 1        : 0 -- standardized id (5 chars), 1 -- exact id, default = 1\n");
    fprintf(stdout,"  --seqmap path         : set the seqmap path, default = $DATADIR/seqmap\n");
    fprintf(stdout,"  --typeprofile 0|1|2|3 : set the data type of profile, default = 1\n");
    fprintf(stdout,"                        : 0 and 1, with SAD,    0 -- for float type, 1 -- for short type\n");
    fprintf(stdout,"                        : 2 and 3, without SAD, 2 -- for float type, 3 -- for short type\n");
    fprintf(stdout,"  -q                    : quiet mode, do not output report\n");
    fprintf(stdout,"  -h|--help       : print this help message and exit\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"Created on 2007-01-09, updated 2009-07-27, Nanjiang Shu\n");
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
void WriteMODM(const char *outfile, char *aaSeq, int length, int **pssm, double* parameter, double *score1, double *score2, int* pssmIndex, char *outAlphabet,  Chain *pSeqmapChain, int type_modm)/*{{{*/
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
    if(type_modm == MODM_LOG)
    { fprintf(fpout,"# pssm is derived from log odd score in PSSM file, \n"); }
    else 
    { fprintf(fpout,"# pssm is derived from weighted percentages in PSSM file, \n"); }

    if(isOutputSAD)
    {
        fprintf(fpout,"# SAD : shape string, water accessibility, secondary structure in dssp\n");
    }

    fprintf(fpout,"# \n");

    fprintf(fpout,"%-5s %-2s","Num","AA");
    if(isOutputSAD)
    {
        fprintf(fpout," %3s", "SAD");
    }

    for(j = 0 ; j < NUM_AA ; j++) fprintf(fpout,"%4c",outAlphabet[j]);

    fprintf(fpout, "\n");
    for(i = 0 ; i < length ; i ++)
    {
        fprintf(fpout,"%5d %-2c",i+1,aaSeq[i]);
        if(isOutputSAD)
        {
            fprintf(fpout, " %c%1d%c", pSeqmapChain->shString[i], pSeqmapChain->waterAcc[i], pSeqmapChain->secStruc[i]);
        }
        for(j = 0 ; j < NUM_AA ; j++) fprintf(fpout,"%4d",pssm[i][pssmIndex[j]]);
        if(score1 != NULL)
            fprintf(fpout,"  %4.2f", score1[i]);
        if(score2!= NULL)
            fprintf(fpout," %4.2f", score2[i]);
        fprintf(fpout,"\n");
    }
    fprintf(fpout,"                      K         Lambda\n"   );
    fprintf(fpout,"Standard Ungapped    %.4lf      %.4lf \n", parameter[0], parameter[1]   );
    fprintf(fpout,"Standard Gapped      %.4lf      %.4lf \n", parameter[2], parameter[3]   );
    fprintf(fpout,"PSI Ungapped         %.4lf      %.4lf \n", parameter[4], parameter[5]   );
    fprintf(fpout,"PSI Gapped           %.4lf      %.4lf \n", parameter[6], parameter[7]   );

    fclose(fpout);
}
/*}}}*/
void WriteBinaryMODMFile(const char *outfile, char *aaSeq, int length, int **pssm, double* parameter, double *score1, double *score2, int* pssmIndex, char *outAlphabet, Chain *pSeqmapChain)/*{{{*/
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
            WriteBinaryMODM(outfile, outAlphabet, length, profile, parameter,typeProfile);
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
    int j = 0;
    char outpath[MAX_PATH+1] = "./";
    char maxaccfile[MAX_PATH+1] = ""; // file listint the maximal water accessibility for each amino acid
    char accscalefile[MAX_PATH+1] = "";// file list the scale for mapping water accessibility to 0~9
    char modmfileext[MAX_PATH+1] = "modm"; //default qijfile extension is Qij
    char outAlphabet[NUM_BLOSUM+1] = ""; //if outAlphabet not set, using the default BLOSUM alphabet
    int type_modm = MODM_LOG;
    int idtype = 0; /*default idtype = 0, meaning standardized 5 character chain identifier*/
    bool isOutputBinaryFile = false;

    char pssmFileList[MAX_PATH+1] =""; /*added 2009-07-15*/

    set <string> filenamelist_set;
    set <string> ::iterator iss;
    //char **filenamelist = NULL;
    //filenamelist = new char*[NUM_FILE];
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
            else if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--outpath") == 0)
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
            else if (strcmp(argv[i], "--wb") == 0)
            {
                isOutputBinaryFile = true;
                i ++;
            }
            else if (strcmp(argv[i], "--typeprofile") == 0)
            {
                if( ( i = option_parser_numeric(argc, argv, i, typeProfile, true, int8(0), int8(3))) == -1)
                    return -1;
            }
            else if (strcmp(argv[i], "--newseq") == 0)
            {
                isNewSequence = true;
                i ++;
            }
            else if (strcmp(argv[i], "-t") == 0)
            {
                if(strncasecmp(argv[i+1], "log", 1) == 0)
                {
                    type_modm = MODM_LOG;
                }
                else 
                {
                    type_modm = MODM_PER;
                }
                i += 2;
            }
            else if (strcmp(argv[i], "--idtype") == 0)
            {
                if( ( i = option_parser_numeric(argc, argv, i, idtype, true, 0, 1)) == -1)
                    return -1;
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
                if( ( i = option_parser_filename(argc, argv, i, modmfileext)) == -1)
                    return -1;
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
        fprintf(stderr,"Error! the typeProfile should be set to 0 or 1 when SAD is enabled.\n");
        return -1;
    }
    if ((! isOutputSAD)&& typeProfile < 2)
    {
        fprintf(stderr,"Error! the typeProfile should be set to 2 or 3 when SAD is disabled.\n");
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

    if(filenamelist_set.size() == 0)
    {
        fprintf(stderr,"Error! no pssm file set\n");
        return -1;
    }

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

    //int numfile = filenamelist_set.size();

    //char alphabet[NUM_BLOSUM+1] = "";

    Array1D <char> aaSeq_1darray(LONGEST_SEQ+1);
    Array2D <int> fij_2darray(LONGEST_SEQ, NUM_BLOSUM);
    Array1D <double> score1_1darray(LONGEST_SEQ);
    Array1D <double> score2_1darray(LONGEST_SEQ);
    Array1D <double> parameter_1darray(NUM_PARAMETER);
    char *aaSeq = aaSeq_1darray.array1D;
    int **fij = fij_2darray.array2D;
    double *score1 = score1_1darray.array1D;
    double *score2 = score2_1darray.array1D;
    double *parameter = parameter_1darray.array1D;


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
    char seqmapfile[MAX_PATH+1] = "";
    char modmfile[MAX_PATH+1] = "";
    char rtname[MAX_PATH+1] = "";
    int length;
    int cntfile = 0;
    for(iss = filenamelist_set.begin() ; iss != filenamelist_set.end(); iss ++)
    {
        for(j = 0  ; j < NUM_PARAMETER; j ++) parameter[j] = 0.0;
        my_strcpy(pssmfile,(*iss).c_str(), MAX_PATH);
        if(type_modm == MODM_LOG)
        {
            if(GetPSSM(pssmfile, length, aaSeq, fij, NULL, score1, score2, parameter) <= 0)
            {
                fprintf(stderr, "getting pssm file error!, pssm file %s", pssmfile);
                break;
            }
        }
        else 
        {
            if(GetPSSM(pssmfile, length, aaSeq, NULL, fij, score1, score2, parameter) <= 0)
            {
                fprintf(stderr, "getting pssm file error!, pssm file %s", pssmfile);
                break;
            }
            if (isTreatAllZeroProfile) /*this works only when using weighted percentages*/
            {
                TreatAllZeroFij(fij, length, score1, aaSeq, BLOSUM1D_alphabet);
            }
        }



        rootname(pssmfile, rtname);
        sprintf(modmfile, "%s/%s.%s", outpath, rtname, modmfileext);
        int pssmIndex[NUM_BLOSUM];
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
                if(idtype == 0)
                {
                    my_strcpy(id,rtname, SIZE_CHAIN_ID);
                    StdID(id);
                    GetSEQMAPFilePath(id, seqmapfile , seqmappath); /*supply seqmappath, 2009-07-15*/
                }
                else if (idtype == 1)
                {
                    my_strcpy(id,rtname, SIZE_CHAIN_ID);
                    sprintf(seqmapfile, "%s/%s.seqmap", seqmappath, id);

                }

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
            sprintf(modmfile, "%s%s", modmfile, "bin");
            WriteBinaryMODMFile(modmfile, aaSeq, length, fij, parameter, score1, score2, pssmIndex, outAlphabet, &seqmap_chain );
        }
        else
        {
            WriteMODM(modmfile, aaSeq, length, fij, parameter, score1, score2, pssmIndex, outAlphabet, &seqmap_chain, type_modm);
        }

        cntfile ++;
        if(!isQuietMode) 
        { fprintf(stdout,"%d \t %s output\n", cntfile, modmfile); }
    }

    if(isOutputSAD)
    {
        DeleteChain(&seqmap_chain);
    }
    if (accScale != NULL){ delete [] accScale; }
    //
    return 0;
}
/*}}}*/
