/*
 * =====================================================================================
 *
 *       Filename:  test-readmodm.cpp
 *    Description:  read modm file, in both text and binary mode
 *        Version:  1.0
 *        Created:  10/31/2007 02:26:42 PM CET
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Nanjiang Shu (Shu), nanjiang@struc.su.se
 *        Company:  Structural Chemistry, Stockholm Univesity
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <string.h>
#include "myfunc.h"
#include "mypro.h"
#include "array.h"
#include <stdlib.h>

#define BLOCK_SIZE 1000
bool isWriteBinaryFile = false;
bool isPrintOut = false;
bool isHaveSAD = false;
/* the structure of the binary file
 * alphabet
 * the number of positions
 * profiles for each position
 * psi-blast parameter*/

#undef SIZE_ALPHABET
#define SIZE_ALPHABET 100

void PrintHelp()
{
    fprintf(stdout,"Usage: test-readmodm [options] modmfile | -l modmfilelist\n");
    fprintf(stdout,"options:\n");
    fprintf(stdout,"  --mode 0|1     : 0 - text, 1 - binary\n");
    fprintf(stdout,"  -o|--out file  : output the result to outfile, default = stdout\n");
    fprintf(stdout,"  --outpath path : output the result to outfile, default = stdout\n");
    fprintf(stdout,"  -a  str        : supply alphabet, force the alphabet output to the binary file as supplied\n");
    fprintf(stdout,"  -b             : write binary file to basename.bin\n");
    fprintf(stdout,"  -p             : print the modm file in ascii file, for checking\n");
    fprintf(stdout,"  -h|--help      : print this help message and exit\n");
    fprintf(stdout,"  --sad          : indicates that the modm file include SAD(shape, waterAcc and dsspSec) information\n");
    fprintf(stdout,"  --notime       : do not output time consuming information\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"Created on 2007-10-31, Last modified on 2007-11-03, Nanjiang Shu\n");
}
void PrintVerboseHelp() { }
void WriteMODM(char *alphabet, Profile* profile, int numRes, double *parameter, FILE *fpout)/*{{{*/
{
    int i,j;
    fprintf(fpout,"#Num: amino acid sequence number\n");
    fprintf(fpout,"#AA: 1 letter amino acid\n");
    fprintf(fpout,"# profile is derived from log odd score in PSSM file, \n");
    fprintf(fpout,"# \n");

    int width = 5;

    fprintf(fpout,"%-5s %-2s","Num","AA");

    for(j = 0 ; j < 20 ; j++) fprintf(fpout,"%4c",alphabet[j]);

    fprintf(fpout, "\n");
    for(i = 0 ; i < numRes ; i ++)
    {
        fprintf(fpout,"%5d %-2c",i+1,profile[i].aa);
        for(j = 0 ; j < 20 ; j++) { fprintf(fpout,"%4d",Integer(profile[i].p[j])); }
        fprintf(fpout,"  %4.2f", profile[i].score1);
        fprintf(fpout," %4.2f", profile[i].score2);
        fprintf(fpout,"\n");
    }
    fprintf(fpout,"                      K         Lambda\n"   );
    fprintf(fpout,"Standard Ungapped    %.4lf      %.4lf \n", parameter[0], parameter[1]   );
    fprintf(fpout,"Standard Gapped      %.4lf      %.4lf \n", parameter[2], parameter[3]   );
    fprintf(fpout,"PSI Ungapped         %.4lf      %.4lf \n", parameter[4], parameter[5]   );
    fprintf(fpout,"PSI Gapped           %.4lf      %.4lf \n", parameter[6], parameter[7]   );
}
/*}}}*/
void WriteMODM(char *alphabet, ProfileSAD* profile, int numRes, double *parameter, FILE *fpout)/*{{{*/
{
    int i,j;
    fprintf(fpout,"#Num: amino acid sequence number\n");
    fprintf(fpout,"#AA: 1 letter amino acid\n");
    fprintf(fpout,"# profile is derived from log odd score in PSSM file, \n");
    fprintf(fpout,"# \n");

    int width = 5;

    fprintf(fpout,"%-5s %-2s","Num","AA");
    fprintf(fpout," %3s","SAD");

    for(j = 0 ; j < 20 ; j++) fprintf(fpout,"%4c",alphabet[j]);

    fprintf(fpout, "\n");
    for(i = 0 ; i < numRes ; i ++)
    {
        fprintf(fpout,"%5d %-2c",i+1,profile[i].aa);
        fprintf(fpout," %1c%1d%1c",profile[i].shape, profile[i].waterAcc, profile[i].dsspSec);
        for(j = 0 ; j < 20 ; j++) { fprintf(fpout,"%4d",Integer(profile[i].p[j])); }
        fprintf(fpout,"  %4.2f", profile[i].score1);
        fprintf(fpout," %4.2f", profile[i].score2);
        fprintf(fpout,"\n");
    }
    fprintf(fpout,"                      K         Lambda\n"   );
    fprintf(fpout,"Standard Ungapped    %.4lf      %.4lf \n", parameter[0], parameter[1]   );
    fprintf(fpout,"Standard Gapped      %.4lf      %.4lf \n", parameter[2], parameter[3]   );
    fprintf(fpout,"PSI Ungapped         %.4lf      %.4lf \n", parameter[4], parameter[5]   );
    fprintf(fpout,"PSI Gapped           %.4lf      %.4lf \n", parameter[6], parameter[7]   );
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
    char modmfile[MAX_PATH+1] = "";
    char modmfilelist[MAX_PATH+1] = "";
    char outfile[MAX_PATH+1] = "";
    char suppliedAlphabet[SIZE_ALPHABET+1] = ""; /*force supplied alphabet to the alphabet in the output binary file*/
    char outpath[MAX_PATH+1] = ""; /*default output path for the binary file is the same as the ascii file*/
    double value = 0.0;
    int mode = 0;
    const char control_option[] = "bp"; //options which control the program, and does not take parameters

    char tmpfile[MAX_PATH+1] = "/tmp/fileXXXXXX";
    if (mkstemp(tmpfile) == -1)
    {
        fprintf(stderr, "Error! No temporary file can be created\n");
        return -1;
    }
    FILE *fptmp = fopen(tmpfile,"w");
    bool isIDsSet = false;
    bool isOutputTimeInfo = true;

    i = 1;
    while(i < argc)/*{{{*/
    {
        if(argv[i][0] == '-' && !isNonOptionArg) //options
        {
            isNonOptionArg = false;
            if(IsInCharSet(argv[i][1], control_option))//if argv[i][1] is in control_option, it might be used as -aqs
            {
                for(j = 1 ; j < strlen(argv[i]); j++)
                {
                    switch (argv[i][j])
                    {
                        case 'b': isWriteBinaryFile = true; break;
                        case 'p': isPrintOut = true; break;
                        default : fprintf(stderr,"Invalid option, non-control option '%c' can be used together with contorl-option\n", argv[j]); return -1;
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
            else if(strcmp(argv[i],"-a") == 0 )
            {
                my_strcpy(suppliedAlphabet, argv[i+1], SIZE_ALPHABET);
                i += 2;

            }
            else if( (strcmp(argv[i],"-o") == 0) || (strcmp(argv[i], "--out") == 0))  
            {
                if( ( i = option_parser_filename(argc, argv, i, outfile)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--outpath") == 0))  
            {
                if( ( i = option_parser_filename(argc, argv, i, outpath)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"-l") == 0) )  
            {
                if( ( i = option_parser_filename(argc, argv, i, modmfilelist)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--mode") == 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, mode, true, 0, 1)) == -1)
                    return -1;
            }
            else if (strcmp(argv[i], "--sad") == 0)/*indicating the modm file contains SAD information*/
            {
                isHaveSAD = true;
                i ++;
            }
            else if (strcmp(argv[i], "--notime") == 0)/*do not output time info*/
            {
                isOutputTimeInfo = false;
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
        else //non-option argument
        {
            if(!isIDsSet) 
            {
                fptmp = fopen(tmpfile,"w");
                if(fptmp == NULL)
                {
                    fprintf(stderr,"can not open file '%s' for write\n", tmpfile);
                    assert( fptmp != NULL);
                }
            }
            fprintf(fptmp,"%s\n",argv[i]);
            isIDsSet = true;
            i ++;
        }
    }/*}}}*/

    if(isIDsSet) fclose(fptmp);


    if(!isIDsSet && strcmp(modmfilelist, "") == 0)
    {
        fprintf(stderr,"Error! neither ids and idListFile set\n");
        PrintHelp();
        return -1;
    }
    else if(isIDsSet && strcmp(modmfilelist,"") == 0)
    {
        my_strcpy(modmfilelist,tmpfile,MAX_PATH);
    }

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

    clock_t start, finish;
    double duration;
    int tmp;
    start = clock();

    FILE *fpFileList = fopen(modmfilelist, "r");
    checkfilestream(fpFileList, modmfilelist, "r");
    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;

    int sizeAlphabet = 0;
    char alphabet[SIZE_ALPHABET+1] = "";
    int numRes = 0;
    Array1D <Profile> profile_1darray(2000);
    Profile *profile = profile_1darray.array1D;
    Array1D <ProfileSAD> profileSAD_1darray(2000);
    ProfileSAD *profileSAD = profileSAD_1darray.array1D;
    double parameter[8];
    int status_sscanf = 0;

    int cntFile  = 0;
    while((linesize = fgetline(fpFileList,line, maxline)) != EOF)
    {
        my_strcpy(modmfile, line, MAX_PATH);

        for(j = 0; j < 8; j ++)
        {
            parameter[j] = 0.0;
        }
        strcpy(alphabet,"");


        if(mode == 0) /*read text file*/
        {
            FILE *fpin = fopen(modmfile,"r");
            checkfilestream(fpin,modmfile,"r");
            char str[100] = "";
            char first_non_blank_char = ' ';
            int cntRes = 0;
            while((linesize = fgetline(fpin, line, maxline)) != EOF)
            {
                if(line[0] == '#' || linesize <=0 ) 
                {
                    continue;
                }
                sscanf(line, " %c", &first_non_blank_char);
                if (first_non_blank_char == 'N')/*alphabet line*/
                {
                    SpanExcluding(line, str);
                    if(!isHaveSAD)
                    {
                        my_strcpy(alphabet, str+5, SIZE_ALPHABET);
                    }
                    else
                    {
                        my_strcpy(alphabet, str+8, SIZE_ALPHABET);
                    }
                }
                else if (isdigit(first_non_blank_char))
                {
                    if(!isHaveSAD)
                    {
                        profile[cntRes].score1 = 0.0;
                        profile[cntRes].score2 = 0.0;
                        status_sscanf = sscanf(line,"%d %c %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",/*{{{*/
                                &profile[cntRes].aaSeqIndex, 
                                &profile[cntRes].aa,
                                &profile[cntRes].p[0],
                                &profile[cntRes].p[1],
                                &profile[cntRes].p[2],
                                &profile[cntRes].p[3],
                                &profile[cntRes].p[4],
                                &profile[cntRes].p[5],
                                &profile[cntRes].p[6],
                                &profile[cntRes].p[7],
                                &profile[cntRes].p[8],
                                &profile[cntRes].p[9],
                                &profile[cntRes].p[10],
                                &profile[cntRes].p[11],
                                &profile[cntRes].p[12],
                                &profile[cntRes].p[13],
                                &profile[cntRes].p[14],
                                &profile[cntRes].p[15],
                                &profile[cntRes].p[16],
                                &profile[cntRes].p[17],
                                &profile[cntRes].p[18],
                                &profile[cntRes].p[19],
                                &profile[cntRes].score1,
                                &profile[cntRes].score2
                                    );/*}}}*/
                        if(status_sscanf < 22)
                        {
                            fprintf(stderr,"current line=%s\n!", line );
                            fprintf(stderr,"The profile maybe incomplete, status_sscanf = %d, number %d, file %s\n!", status_sscanf, cntRes+1, modmfile );
                            assert (status_sscanf >= 22);
                        }
                    }
                    else
                    {
                        profileSAD[cntRes].score1 = 0.0;
                        profileSAD[cntRes].score2 = 0.0;
                        status_sscanf = sscanf(line,"%d %c %1c%1d%1c %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",/*{{{*/
                                &profileSAD[cntRes].aaSeqIndex, 
                                &profileSAD[cntRes].aa,
                                &profileSAD[cntRes].shape,
                                &profileSAD[cntRes].waterAcc,
                                &profileSAD[cntRes].dsspSec,
                                &profileSAD[cntRes].p[0],
                                &profileSAD[cntRes].p[1],
                                &profileSAD[cntRes].p[2],
                                &profileSAD[cntRes].p[3],
                                &profileSAD[cntRes].p[4],
                                &profileSAD[cntRes].p[5],
                                &profileSAD[cntRes].p[6],
                                &profileSAD[cntRes].p[7],
                                &profileSAD[cntRes].p[8],
                                &profileSAD[cntRes].p[9],
                                &profileSAD[cntRes].p[10],
                                &profileSAD[cntRes].p[11],
                                &profileSAD[cntRes].p[12],
                                &profileSAD[cntRes].p[13],
                                &profileSAD[cntRes].p[14],
                                &profileSAD[cntRes].p[15],
                                &profileSAD[cntRes].p[16],
                                &profileSAD[cntRes].p[17],
                                &profileSAD[cntRes].p[18],
                                &profileSAD[cntRes].p[19],
                                &profileSAD[cntRes].score1,
                                &profileSAD[cntRes].score2
                                    );/*}}}*/
                        if(status_sscanf < 25)
                        {
                            fprintf(stderr,"current line=%s\n!", line );
                            fprintf(stderr,"The profile maybe incomplete, status_sscanf = %d, number %d, file %s\n!", status_sscanf, cntRes+1, modmfile );
                            assert (status_sscanf >= 25);
                        }
                    }
                    cntRes ++;
                }
                else if (first_non_blank_char == 'K')
                {
                    int j = 0 ;
                    while((linesize = fgetline(fpin, line, maxline)) != EOF)
                    {
                        status_sscanf = sscanf(line,"%s %s %lf %lf", str, str, &parameter[j], &parameter[j+1]);
                        if(status_sscanf < 4)
                        {
                            fprintf(stderr,"Warning, file %s may not have psi-blast parameters\n", modmfile);
                        }
                        j += 2;
                        if (j >= 8 ) break;
                    }
                }
            }
            fclose(fpin);
            numRes = cntRes;
            if (isWriteBinaryFile)
            {
                char binaryfile[MAX_PATH+1] = "";
                if(strcmp(outpath, "") == 0)  /*if outpath is not set, output to the same folder as the ascii file*/
                {
                    sprintf(binaryfile, "%sbin", modmfile);   /*if outpath is set*/
                }
                else
                {

                    char rtname[MAX_PATH+1] = "";
                    char fileext[MAX_PATH+1] = "";
                    rootname(modmfile, rtname);
                    getfileext(modmfile, fileext);
                    sprintf(binaryfile, "%s/%s.%sbin", outpath, rtname, fileext);
                }
//                fprintf(stdout,"%d \t write out binary file to %s\n", cntFile, binaryfile) ;
                if(strcmp(suppliedAlphabet,"") != 0)
                {
                    my_strcpy(alphabet, suppliedAlphabet, SIZE_ALPHABET);
                }
                if(!isHaveSAD)
                {
                    WriteBinaryMODM(binaryfile, alphabet, numRes, profile, parameter);
                }
                else
                {
                    WriteBinaryMODM(binaryfile, alphabet, numRes, profileSAD, parameter);
                }
            }
        }
        else if(mode == 1)
        {

            if(!isHaveSAD)
            {
                GetBinaryMODM(modmfile, alphabet, numRes, profile, parameter);
                if (isPrintOut)
                {
                    WriteMODM(alphabet, profile, numRes, parameter, fpout);
                }
            }
            else
            {
                GetBinaryMODM(modmfile, alphabet, numRes, profileSAD, parameter);
                if (isPrintOut)
                {
                    WriteMODM(alphabet, profileSAD, numRes, parameter, fpout);
                }
            }
        }
        cntFile ++;
    }

    finish = clock();
    duration = double(finish-start)  /double(CLOCKS_PER_SEC);
    //printf("CLOCKS_PER_SEC=%d, start=%lf, finish=%lf\n", CLOCKS_PER_SEC, start, finish);
    if (isOutputTimeInfo)
    {
        if(mode == 0)
        {
            fprintf(stdout,"reading text file cost %lf seconds\n", duration);
        }
        else if(mode == 1)
        {
            fprintf(stdout,"reading bindary file cost %lf seconds\n", duration);
        }

    }

    if(fpout != NULL && fpout != stdout) fclose(fpout);

    return 0;
}
/*}}}*/
