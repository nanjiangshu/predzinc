/*
 * =====================================================================================
 * 
 *        Filename:  gating_gistPred.cpp
 *     Description:  gating network for combining the predicitons using
 *     				 single-site vectors and pair-based vectors
 *         Version:  1.0
 *         Created:  07/28/2006 11:34:36 AM CEST
 *        Revision:  none
 *        Compiler:  gcc
 *          Author:  Nanjiang Shu (Shu), nanjiang@struc.su.se
 *         Company:  Structural Chemistry, Stockholm Univesity
 * 
 * =====================================================================================
 */
#include <algorithm>
#include <cmath>
#include <cassert>
#include <cstring>
#include "array.h"
#include "mytemplate.h"
#include "myfunc.h"
#include "mypro.h"

#define USE_AVERAGE 0
#define USE_MAXIMUM 1

void PrintHelp()
{
    fprintf(stdout,"Usage: gating_gistPred [options] -i gistPredFile1 gistPredFile2\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"options\n");
    fprintf(stdout,"  -ns site1 site2     : number of sites used for gistPredFile1 and  gistPredFile2, default = 1 2\n");
    fprintf(stdout,"  -p1 A B             : slope and shift for gistPredFile1  default = 2 0.5\n");
    fprintf(stdout,"  -p2 A B             : slope and shift for gistPredFile2  default = 4 0.5\n");
    fprintf(stdout,"  --operation avg|max : setting operation for merging multi residues, avg, average, default = avg\n");
    fprintf(stdout,"  -o outfile          : outfile,  default = stdout\n");
    fprintf(stdout,"  -h | --help         : print this help message and exit\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"created on 2006-08-20, last modifed 2007-06-20, Nanjiang Shu\n");
}
double SigmoidScore(double a, double b, double x)/*{{{*/
{
    return 1.0/(1.0+exp(-a*x-b));
}
/*}}}*/
double GatingScore(double x, double y)/*{{{*/
{
    return x + (1-x)*y;
}
/*}}}*/
int GatingGistPred(const char* gistPredictFile1, const char* gistPredictFile2, int numSite1,int numSite2,double a1,double b1,double a2,double b2, int operation, FILE *fpout)/*{{{*/
{

    FILE *fpGistPred1;
    FILE *fpGistPred2;
    //double p1;
    //double p2;
     
    GistPredChain *pChain = NULL;

    int i,j;
    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D; 


    set <string> idList1_set;
    set <string> idList2_set;
    set <string> idList_union_set;
    set <string> ::iterator iss;
    char id[SIZE_CHAIN_ID+1] = "";

    fpGistPred1 = fopen(gistPredictFile1,"r");
    fpGistPred2 = fopen(gistPredictFile2,"r");
    checkfilestream(fpGistPred1, gistPredictFile1, "r");
    checkfilestream(fpGistPred2, gistPredictFile2, "r");
    // get the number of unique ids in file1 and file2/*{{{*/

    while((linesize = fgetline(fpGistPred1, line, maxline)) != EOF)
    {
        if(linesize > 0 && line[0] != ' '&& line[0] != '#' && strncmp(line,"item",4) != 0)
        {
            ssubstitute(line,CHAR_VECTOR_ID_SEPRATOR,'\0');
            my_strcpy(id,line,SIZE_CHAIN_ID);
            idList1_set.insert(id);
        }
    }

    while((linesize = fgetline(fpGistPred2, line, maxline)) != EOF)
    {
        if(linesize > 0 && line[0] != ' '&& line[0] != '#' && strncmp(line,"item",4) != 0)
        {
            ssubstitute(line,CHAR_VECTOR_ID_SEPRATOR,'\0');
            my_strcpy(id,line,SIZE_CHAIN_ID);
            idList2_set.insert(id);
        }
    }
    /*}}}*/
    fclose(fpGistPred1);
    fclose(fpGistPred2);

    int numChain1 = idList1_set.size();
    int numChain2 = idList2_set.size();
    set_union(idList1_set.begin(), idList1_set.end(),idList2_set.begin(), idList2_set.end(), inserter(idList_union_set, idList_union_set.begin()));
    int numChainUnion = idList_union_set.size();
    Array2D <char> idList_union_2darray(numChainUnion, SIZE_CHAIN_ID+1);
    char ** idList_union = idList_union_2darray.array2D;
    i = 0 ;
    for ( iss = idList_union_set.begin(); iss != idList_union_set.end(); iss ++)
    {
        my_strcpy(idList_union[i], (*iss).c_str(), SIZE_CHAIN_ID);
        i ++;
    }
        
    idList1_set.clear();
    idList2_set.clear();
    //int  length;
    //int  numSite = max(numSite1, numSite2);
    //Array1D <int> aaSeqIndex_1darray(numSite+1);
    //Array1D <char> res_1char_list_1darray(numSite+2);
    //int idx;
    //int  *aaSeqIndex = aaSeqIndex_1darray.array1D;
    //char *res_1char_list = res_1char_list_1darray.array1D;

    Array1D <GistPredChain> chains1_1darray(numChain1+1);
    Array1D <GistPredChain> chains2_1darray(numChain2+1);
    GistPredChain *chains1 = chains1_1darray.array1D;
    GistPredChain *chains2 = chains2_1darray.array1D;
    for(i = 0 ; i < numChain1 ; i++) InitGistPredChain(&chains1[i]);
    for(i = 0 ; i < numChain2 ; i++) InitGistPredChain(&chains2[i]);

    
    Array2D <char> idList1_2darray(numChain1, SIZE_CHAIN_ID+1);
    Array2D <char> idList2_2darray(numChain2, SIZE_CHAIN_ID+1);
    char **idList1 = idList1_2darray.array2D;
    char **idList2 = idList2_2darray.array2D;
    //----------------------------------------------------------------------------
    /*get residues in each gist prediction file */
    //----------------------------------------------------------------------------
    //operation is when there are multiple discriminants for the same residue,
    //which rule to use to integrate them. average or maximum
    GetGistPredResidue(gistPredictFile1, numSite1, chains1, operation);
    GetGistPredResidue(gistPredictFile2, numSite2, chains2, operation);
    
    // --- sort the chains according to id, and get the idx number
    Array1D <int> idx1_1darray(numChain1);
    Array1D <int> idx2_1darray(numChain2);
    int *idx1 = idx1_1darray.array1D;
    int *idx2 = idx2_1darray.array1D;
    for( i = 0 ; i < numChain1 ; i ++) idx1[i] = i;
    for( i = 0 ; i < numChain2 ; i ++) idx2[i] = i;

    Array2D <char> tmp_idList_2darray(max(numChain1, numChain2), SIZE_CHAIN_ID+1);
    char **tmp_idList = tmp_idList_2darray.array2D;

    for( i = 0 ; i < numChain1 ; i ++) 
        my_strcpy(tmp_idList[i], chains1[i].id, SIZE_CHAIN_ID);
    QuickSort_String(idx1, tmp_idList, 0, numChain1-1);
    for( i = 0 ; i < numChain1 ; i ++) 
        my_strcpy(idList1[i], tmp_idList[idx1[i]],  SIZE_CHAIN_ID);


    for( i = 0 ; i < numChain2 ; i ++) 
        my_strcpy(tmp_idList[i], chains2[i].id, SIZE_CHAIN_ID);
    QuickSort_String(idx2, tmp_idList, 0, numChain2-1);
    for( i = 0 ; i < numChain2 ; i ++) 
        my_strcpy(idList2[i], tmp_idList[idx2[i]],  SIZE_CHAIN_ID);

    int label;
    //double discriminant;
    double gating_score;
    char rmtID[SIZE_CHAIN_ID+1] = "";
    char vectorRecordID[300+1] = "";

    int proIndex1 = 0;
    int proIndex2 = 0;
    for( i = 0 ; i < numChainUnion ; i++)
    {
        proIndex1 = BinarySearch_String(idList_union[i], idList1, numChain1);
        proIndex2 = BinarySearch_String(idList_union[i], idList2, numChain2);
        my_strcpy(rmtID, idList_union[i], SIZE_CHAIN_ID);
        RemoveTID(rmtID);
        
        if(proIndex1 != -1 && proIndex2 != -1)
        {    // 
            set <int> aaSeqIndex_union_set;
            //set <int> aaSeqIndex1_set;
            //set <int> aaSeqIndex2_set;

            pChain = &(chains1[idx1[proIndex1]]);
            int numRes1 = pChain->numRes;
            Array1D <int> aaSeqIndex1_1darray(numRes1);
            int *aaSeqIndex1 = aaSeqIndex1_1darray.array1D;
            for(j = 0 ; j < pChain->numRes; j ++) 
            { 
                aaSeqIndex1[j] = pChain->aaSeqIndex[j];
                aaSeqIndex_union_set.insert(pChain->aaSeqIndex[j]);
            }

            pChain = &(chains2[idx2[proIndex2]]);
            int numRes2 = pChain->numRes;
            Array1D <int> aaSeqIndex2_1darray(numRes2);
            int *aaSeqIndex2 = aaSeqIndex2_1darray.array1D;
            for(j = 0 ; j < pChain->numRes; j ++) 
            { 
                aaSeqIndex2[j] = pChain->aaSeqIndex[j];
                aaSeqIndex_union_set.insert(pChain->aaSeqIndex[j]);
            }

            int numRes_union = aaSeqIndex_union_set.size();
            Array1D <int> aaSeqIndex_union_1darray(numRes_union);
            int *aaSeqIndex_union = aaSeqIndex_union_1darray.array1D;
            Set2Array(aaSeqIndex_union_set.begin(),aaSeqIndex_union_set.end(), aaSeqIndex_union);

            for(j = 0 ; j < numRes_union; j ++)
            {
                int idx_aaSeqIndex1 = binarysearch(aaSeqIndex_union[j], aaSeqIndex1, numRes1);
                int idx_aaSeqIndex2 = binarysearch(aaSeqIndex_union[j], aaSeqIndex2, numRes2);
                char *pRes_1char_list;
                int  *pAASeqIndex;

                if(idx_aaSeqIndex1 != -1 && idx_aaSeqIndex2 != -1)
                {
                    double sigmoid1;
                    double sigmoid2;
                    pChain = &(chains1[idx1[proIndex1]]);
                    sigmoid1 = SigmoidScore(a1,b1,pChain->discriminant[idx_aaSeqIndex1]);
                    pChain = &(chains2[idx2[proIndex2]]);
                    sigmoid2 = SigmoidScore(a2,b2,pChain->discriminant[idx_aaSeqIndex2]);
                           
                    gating_score = GatingScore(sigmoid1, sigmoid2);
                    pRes_1char_list = &(pChain->aaSeq[idx_aaSeqIndex2]);
                    pAASeqIndex     = &(pChain->aaSeqIndex[idx_aaSeqIndex2]);
                }
                else
                { 
                    double sigmoid;
                    if(idx_aaSeqIndex1 != -1)
                    {
                        pChain = &(chains1[idx1[proIndex1]]);
                        sigmoid = SigmoidScore(a1,b1,pChain->discriminant[idx_aaSeqIndex1]);
                        pRes_1char_list = &(pChain->aaSeq[idx_aaSeqIndex1]);
                        pAASeqIndex     = &(pChain->aaSeqIndex[idx_aaSeqIndex1]);
                    }
                    else if(idx_aaSeqIndex2 != -1)
                    {
                        pChain = &(chains2[idx2[proIndex2]]);
                        sigmoid = SigmoidScore(a2,b2,pChain->discriminant[idx_aaSeqIndex2]);
                        pRes_1char_list = &(pChain->aaSeq[idx_aaSeqIndex2]);
                        pAASeqIndex     = &(pChain->aaSeqIndex[idx_aaSeqIndex2]);
                    }
                    else
                    {
                        printf("Error! aaSeqIndex_union neither in aaSeqIndex1 or in aaSeqIndex2\n");
                        assert(idx_aaSeqIndex1 != -1 || idx_aaSeqIndex2 != -1);
                    }
                    gating_score = sigmoid;

                }
                label = (gating_score >= 0.5) ? (1) : (-1);
                WriteVectorRecordID(vectorRecordID, rmtID, pChain->length, j, pRes_1char_list, pAASeqIndex, 1);
                fprintf(fpout,"%s \t %2d \t %10.6lf\n",vectorRecordID, label, gating_score);
            }

        }
        else
        {
            double sigmoid;
			double *pA = NULL;
			double *pB = NULL;
            if(proIndex1 != -1)
            {
                pChain = &(chains1[idx1[proIndex1]]);
                pA = &a1; pB = &b1;
            }
            else if(proIndex2 != -1)
            {
                pChain = &(chains2[idx2[proIndex2]]);
                pA = &a2; pB = &b2;
            }
            else
            {
                fprintf(stderr, "Error! idListUnion neither in idList1 or in idList2\n");
            }

            for(j = 0 ; j < pChain->numRes ; j++)
            {
                sigmoid = SigmoidScore(*pA, *pB, pChain->discriminant[j]);
                gating_score = sigmoid;
                label = (gating_score >= 0.5) ? (1) : (-1);
                WriteVectorRecordID(vectorRecordID, rmtID, pChain->length, j, &(pChain->aaSeq[j]), &(pChain->aaSeqIndex[j]), 1);
                fprintf(fpout,"%s \t %2d \t %10.6lf\n",vectorRecordID, label, gating_score);
            }
        }
    }

    
    // free memory
    for(i = 0 ; i < numChain1 ; i++) DeleteGistPredChain(&chains1[i]);
    for(i = 0 ; i < numChain2 ; i++) DeleteGistPredChain(&chains2[i]);
	return 0;
}
/*}}}*/

int main(int argc, char** argv)/*{{{*/
{

    if( argc < 2 )
    {
        printf("too few arguments\n");
        PrintHelp();
        return -1;
    }

    int    i;
    int    numSite1                     = 1;
    int    numSite2                     = 2;
    char   outfile[MAX_PATH+1]          = "";
    char   gistPredictFile1[MAX_PATH+1] = "";
    char   gistPredictFile2[MAX_PATH+1] = "";
    double a1 = 2.0;
    double b1 = 0.5;
    double a2 = 4.0;
    double b2 = 0.5;
    int operation = USE_AVERAGE;

    i = 1;
    while(i < argc)
    {
        if( strcmp(argv[i],"-h") == 0 || strcmp(argv[i],"--help") == 0)
        {
            PrintHelp();
            return 0;
        }
        else if(strcmp(argv[i],"-ns") == 0)
        {
            if(IsNumeric(argv[i+1]))
                numSite1 = atoi(argv[i+1]);
            else
                fprintf(stderr,"Argument Error! two integer values should follow argument: -ns\n");

            if(IsNumeric(argv[i+2]))
                numSite2 = atoi(argv[i+2]);
            else
                fprintf(stderr,"Argument Error! two integer values should follow argument: -ns\n");

            i += 3;
        }
        else if(strcmp(argv[i],"--operation") == 0)
        {
            if(strcasecmp(argv[i+1], "avg") == 0 || strcasecmp(argv[i+1], "average") == 0)
                operation = USE_AVERAGE;
            else if(strcasecmp(argv[i+1], "max") == 0 || strcasecmp(argv[i+1], "maximum") == 0)
                operation = USE_MAXIMUM;
            else
                fprintf(stderr,"Argument Error! use avg or max for argument : '%s'\n", "--operation");

            i += 2;
        }
        else if(strcmp(argv[i],"-p1") == 0)
        {
            if(IsNumeric(argv[i+1]))
                a1 = atof(argv[i+1]);
            else
                fprintf(stderr,"Argument Error! two double values should follow argument: -p1\n");

            if(IsNumeric(argv[i+2]))
                b1 = atof(argv[i+2]);
            else
                fprintf(stderr,"Argument Error! two double values should follow argument: -p1\n");
            i += 3;
        }
        else if(strcmp(argv[i],"-p2") == 0)
        {
            if(IsNumeric(argv[i+1]))
                a2 = atof(argv[i+1]);
            else
                fprintf(stderr,"Argument Error! two double values should follow argument: -p2\n");

            if(IsNumeric(argv[i+2]))
                b2 = atof(argv[i+2]);
            else
                fprintf(stderr,"Argument Error! two double values should follow argument: -p2\n");
            i += 3;
        }
        else if(strcmp(argv[i],"-i") == 0)
        {
            if(argc-i < 3)
                fprintf(stderr,"Argument Error! two filenames should follow argument: -i\n");
            my_strcpy(gistPredictFile1, argv[i+1], MAX_PATH);
            my_strcpy(gistPredictFile2, argv[i+2], MAX_PATH);
            i += 3;
        }
        else if(strcmp(argv[i],"-o") == 0)
        {
            my_strcpy(outfile, argv[i+1], MAX_PATH);
            i += 2;
        }
        else
        { 
            fprintf(stderr, "Wrong argument: %s\n",argv[i] );
            PrintHelp();
            return -1;
        }
    }

    if(strcmp(gistPredictFile1,"") == 0 || strcmp(gistPredictFile2,"") == 0)
    {
        fprintf(stderr,"Error! Filename of gistPredictFile not set\n");
        return -1;
    }

    FILE *fpout;
    if(strcmp(outfile,"") == 0)
    {
        fpout = stdout;
    }
    else
    {
        fpout = fopen(outfile,"w");
        if(fpout == NULL)
        {
            fprintf(stderr,"can not open file '%s' for write\n", outfile);
            assert(fpout != NULL);
        }
    }
    GatingGistPred(gistPredictFile1, gistPredictFile2, numSite1,numSite2,a1,b1,a2,b2, operation, fpout);
}
/*}}}*/
