/*
 * =====================================================================================
 * 
 *        Filename:  mtx2modm.cpp
 *     Description:  convert mtx of the makemat output to pssm matrix file
 *                   the score in pssm matrix file is log odds score
 *         Version:  1.0
 *         Created:  12/06/2006 05:01:43 PM CET
 *        Compiler:  g++
 *          Author:  Nanjiang Shu (Shu), nanjiang@struc.su.se
 *         Company:  Structural Chemistry, Stockholm Univesity
 * 
 * =====================================================================================
 */

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "mytemplate.h"
#include "array.h"
#include "myfunc.h"
#include "mypro.h"

const char *rescodes = "ARNDCQEGHILKMFPSTWYVBZX";
const char *ncbicodes= "XAXCDEFGHIKLMNPQRSTVWXYXXX";

#undef NUM_PARAMETER
#define NUM_PARAMETER 12

#undef NUM_FILE
#define NUM_FILE 20000

#undef NUM_AA
#define NUM_AA 20

bool isOutputInteger=false;

void PrintHelp()
{
    fprintf(stdout,"Usage: mtx2modm [options] mtxfile [more mtxfiles ...]\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"Options: \n");
    fprintf(stdout,"  -d outpath    : set the output path, default = ./\n");
    fprintf(stdout,"  -s scale      : scalue for the profiles in mtx file, default = 100.0\n");
    fprintf(stdout,"  -i            : output the integer rounded profile, default is real number\n");
    fprintf(stdout,"  -l listfile   : file list for mtx file, one file per line\n");
    fprintf(stdout,"  --pssm        : output score1 and score2 output, pssm files should be in the same\n");
    fprintf(stdout,"                : folder as mtx file, with the same rootname\n");
    fprintf(stdout,"  -a <string>   : alphabet for profile, default using BLOSUM_alphabet\n");
    fprintf(stdout,"  -q            : quiet mode, do not output report\n");
    fprintf(stdout,"  -h|--help     : print this help message and exit\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"Created 2007-01-17, last modifed 2007-10-10, Nanjiang Shu\n");

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
int aanum(int ch)/*{{{*/
/* Convert AA letter to numeric code (0-22) */
    // convert the aa in ncbicodes order to blosum order
    // say, 'A' in ncbicodes is 1, aacvs[1] = 0, which is the indexing of 'A' 
    // in blosum alphabet
{
    static int aacvs[] =
    {
        999, 0, 20, 4, 3, 6, 13, 7, 8, 9, 22, 11, 10, 12, 2, 22, 14, 5, 1, 15, 16, 22, 19, 17, 22, 18, 21
    };
    return (isalpha(ch) ? aacvs[ch & 31] : 22);
}
/*}}}*/
int getmtx(const char* mtxfile, int &length, char *aaSeq, double* parameter, double **profile, double scale)/*{{{*/
/* Read PSI AA frequency data from .mtx file, parameter is for storing twelve
 * parameters*/

{
    int i,j;
    FILE *fpmtx;
    fpmtx = fopen(mtxfile, "r");
    checkfilestream(fpmtx, mtxfile, "r");

    int linesize;
    int maxline = LONGEST_SEQ;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D; 


    fgetline(fpmtx, line, maxline);
    if (sscanf(line, "%d", &length) != 1)
    {    
        fprintf(stderr, "Bad mtx file - no sequence length!\n");
        return (-1);
    }
    if (length > LONGEST_SEQ)
    {
        fprintf(stderr, "Error! Input sequence longer than LONGEST_SEQ=%d!\n", LONGEST_SEQ);
        return (-1);
    }

    fgetline(fpmtx, line, maxline);
    if (sscanf(line, "%s", aaSeq) != 1)
    {
        fprintf(stderr, "Bad mtx file - no sequence!\n");
        return (-1);
    }

    i = 0; 
    while((linesize = fgetline(fpmtx, line, maxline)) != EOF)
    {
        if(sscanf(line, "%lf", &parameter[i++]) != 1)
        {
            fprintf(stderr, "Bad mtx file - no parameter %d!\n", i-1);
            return (-1);
        }
        if(i >= NUM_PARAMETER) break;
    }

    int tmp;
    Array1D <int> pi_1darray(NUM_BLOSUM);
    int *pi = pi_1darray.array1D;
    i = 0;
    while((linesize = fgetline(fpmtx, line, maxline)) != EOF)
    {
        if (sscanf(line, "%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d", &tmp, &pi[ALA],&tmp, &pi[CYS], &pi[ASP], &pi[GLU], &pi[PHE], &pi[GLY], &pi[HIS],  &pi[ILE],  &pi[LYS],  &pi[LEU],  &pi[MET],  &pi[ASN], &pi[PRO], &pi[GLN], &pi[ARG],  &pi[SER],  &pi[THR],  &pi[VAL],  &pi[TRP],&tmp,&pi[TYR]) != 23)
        {
            fprintf(stderr, "Bad mtx format!\n");
            return (-1);
        }
        for(j = 0 ; j < 20; j ++) profile[i][j] = double(pi[j]) /scale;
        i ++;
        if (i >= length) break;
    }
    fclose(fpmtx);

    if(i < length)
    {
        fprintf(stderr, "Bad mtx format! sequence length not equals number of profile vectors!\n");
        return (-1);
    }
    
    return length;
}

/*}}}*/
void WriteMODMFile(const char *modmfile, char *aaSeq, int length, double **profile, int *profileIndex, char *outAlphabet, double *score1, double *score2, double* parameter, bool isOutputInteger = false)/*{{{*/
{
    FILE *fpout;
    fpout = fopen(modmfile, "w");
    checkfilestream(fpout, modmfile, "w");
    int i,j;
    fprintf(fpout,"#Num: amino acid sequence number\n");
    fprintf(fpout,"#AA: 1 letter amino acid\n");
    fprintf(fpout,"# pssm is converted from mtx file of makemat output\n");
    
    //int width = 5;

    fprintf(fpout,"%4s %-2s","Num","AA");

    if(isOutputInteger)
        for(j = 0 ; j < NUM_AA ; j++) fprintf(fpout,"%4c",outAlphabet[j]);
    else
        for(j = 0 ; j < NUM_AA ; j++) fprintf(fpout," %5c",outAlphabet[j]);

    fprintf(fpout, "\n");
    for(i = 0 ; i < length ; i ++)
    {

        fprintf(fpout,"%4d %-2c",i+1,aaSeq[i]);
        if(isOutputInteger)
            for(j = 0 ; j < NUM_AA ; j++) fprintf(fpout,"%4d",Integer(profile[i][profileIndex[j]]));
        else
            for(j = 0 ; j < NUM_AA ; j++) fprintf(fpout," %5.2lf",profile[i][profileIndex[j]]);

        if(score1 != NULL)
            fprintf(fpout,"  %4.2f", score1[i]);
        if(score2!= NULL)
            fprintf(fpout," %4.2f", score2[i]);

        fprintf(fpout,"\n");
    }
    fprintf(fpout,"# gapped PSI-BLAST, lambda = %lf, K = %lf\n", parameter[4], parameter[5]);
    
    fclose(fpout);
}
/*}}}*/
int main(int argc, char** argv)/*{{{*/
{

    if (argc < 2)
    {
        PrintHelp();
        return -1;
    }
    int i = 0;
    int j = 0;
    char outpath[MAX_PATH+1] = "./";
    char listfile[MAX_PATH+1] = "";
    char outAlphabet[NUM_BLOSUM+1] = ""; //if outAlphabet not set, using the default BLOSUM alphabet
    double scale = 100.0;
    //char **filenamelist = NULL;
    bool isOutputLast2Col = false;
    bool isQuietMode = false;

    set <string> filenamelist_set;
    set <string> ::iterator iss;

    int numfile = 0;
    
    int errmsg = 0;
    
    i = 1;
    while(i < argc)/*{{{*/
    {
        if(strcmp(argv[i], "-h") == 0|| strcmp(argv[i], "--help") == 0)
        {
            PrintHelp();
            errmsg = 1;
            break;
        }
        else if (strcmp(argv[i], "-d") == 0)
        {
            my_strcpy(outpath,argv[i+1], MAX_PATH);
            i += 2;
        }
        else if (strcmp(argv[i], "-s") == 0)
        {
            scale = atof(argv[i+1]);
            i += 2;
        }
        else if (strcmp(argv[i], "-i") == 0)
        {
            isOutputInteger = true;
            i ++;
        }
        else if (strcmp(argv[i], "-l") == 0 || strcmp(argv[i], "--list") == 0)
        {
            my_strcpy(listfile,argv[i+1], MAX_PATH);
            i += 2;
        }
        else if (strcmp(argv[i], "--alphabet") == 0 || strcmp(argv[i], "-a") == 0)
        {
            my_strcpy(outAlphabet,argv[i+1], NUM_BLOSUM);
            if(!IsValidAlphabet(outAlphabet))
            {
                errmsg  = -1;
                break;
            }
            i += 2;
        }
        else if (strcmp(argv[i], "-q") == 0)
        {
            isQuietMode = true;
            i ++;
        }
        else if (strcmp(argv[i], "--pssm") == 0)
        {
            isOutputLast2Col = true;
            i ++;
        }
        else
        {
            filenamelist_set.insert(argv[i]);
            i ++;
        }
    }/*}}}*/
    numfile = filenamelist_set.size();

    if(errmsg != 0)
    {
        return errmsg;
    }
    if(strcmp(listfile,"") == 0 && filenamelist_set.size() == 0)
    {
        fprintf(stderr,"Error! neither listfile nor mtx file in the argument list are set\n");
        return -1;
    }
    else if(strcmp(listfile,"") != 0)
    {
        FILE *fpin;
        fpin = fopen(listfile,"r");
        checkfilestream(fpin,listfile,"r");
        int linesize;
        int maxline = 300;
        Array1D <char> line_1darray(maxline+1);
        char *line = line_1darray.array1D;
        while((linesize = fgetline(fpin, line , maxline)) != EOF)
        {
            if(linesize > 0) filenamelist_set.insert(line);
        }
        fclose(fpin);
    }
        
    VerifyFolder(outpath);

    Array1D <char> aaSeq_1darray(LONGEST_SEQ+1);
    Array2D <double> profile_2darray(LONGEST_SEQ, NUM_BLOSUM);
    Array1D <double> parameter_1darray(NUM_PARAMETER);
    char *aaSeq = aaSeq_1darray.array1D;
    double **profile = profile_2darray.array2D;
    double *parameter = parameter_1darray.array1D;
    double *score1 = NULL;
    double *score2 = NULL;
    char modmfile[MAX_PATH+1] = "";
    char pssmfile[MAX_PATH+1] = "";
    char mtxfile[MAX_PATH+1] = "";
    char rtname[MAX_PATH+1] = "";
    char filepath[MAX_PATH+1] = "";
    int length;
    int lengthpssm;

    int cntfile = 0;
    for(iss = filenamelist_set.begin() ; iss !=  filenamelist_set.end(); iss ++)
    {
        score1 = NULL;
        score2 = NULL;

        my_strcpy(mtxfile, (*iss).c_str(), MAX_PATH);
        rootname(mtxfile, rtname);
        getfilepath(mtxfile, filepath);
        if(getmtx(mtxfile, length, aaSeq, parameter, profile, scale) <= 0)
        {
            fprintf(stderr, "%d \t %s getting mtx file error!",cntfile+1, mtxfile);
            continue;
        }
        lengthpssm = length ;

        int profileIndex[NUM_BLOSUM]; //order the profile according to the outAlphabet
        for(j = 0; j < NUM_BLOSUM; j ++) profileIndex[j] = j;
        if(strcmp(outAlphabet,"") != 0)
        {
            for(i = 0 ; i < NUM_AA; i++)
            {
                profileIndex[i] = Char2Digit(outAlphabet[i],  BLOSUM1D_alphabet);
            } 
        }
        else
        {
            my_strcpy(outAlphabet, BLOSUM1D_alphabet, NUM_BLOSUM);
        }


        if(isOutputLast2Col)
        {
            score1 = new double [LONGEST_SEQ];
            score2 = new double [LONGEST_SEQ];
            sprintf(pssmfile, "%s/%s.pssm", filepath, rtname);
            lengthpssm = GetPSSM(pssmfile, lengthpssm, NULL, NULL, NULL, score1, score2);
        }
        if ( length != lengthpssm)
        {
            fprintf(stderr,"Error! length of mtx file and pssm not the same for %s", rtname);
        }
        else
        {
            sprintf(modmfile, "%s/%s.modm", outpath, rtname);
            WriteMODMFile(modmfile, aaSeq, length, profile, profileIndex, outAlphabet, score1, score2, parameter, isOutputInteger);
            if(!isQuietMode) fprintf(stdout,"%d  %s output\n", cntfile+1, modmfile);
        }
        if (score1 != NULL) delete [] score1;
        if (score2 != NULL) delete [] score2;
    }
    //
    return 0;
}
/*}}}*/
