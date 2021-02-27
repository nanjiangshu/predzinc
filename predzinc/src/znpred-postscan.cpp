/*
 * =====================================================================================
 *        Filename:  znpred-postscan.cpp
 *     Description:  combine zinc-binding prediction from SVM predictors and
 *                   homology-based predictions 
 *         Version:  1.0
 *         Created:  12/06/2006 05:01:43 PM CET
 *        Compiler:  g++
 *          Author:  Nanjiang Shu (Shu), nanjiang@struc.su.se
 *         Company:  Structural Chemistry, Stockholm Univesity
 * =====================================================================================
 */
// znpred-postscan.cpp
// proof reading of the predzinc made by the gating network
// the general rule is: 
//  scan the predicted zinc binding residues from the gating network
//  if there are two or three residues within 150 are predicted as zinc binding
//  with discriminant > 0.6, then, there are likely to be zinc binding
//  an iosolated highly predicted residue is unlikely to be a zinc binding
//  residue
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "array.h"
#include "myfunc.h"
#include "mypro.h"
#include "mytemplate.h"
#include "subset.h"

#undef SIZE_VECTOR_RECORD
#define SIZE_VECTOR_RECORD 500

#define RULE_HOMO 0
#define RULE_SCORE 1

char CHDE_Alphabet[] = "CHDE";
double SubMatrixCHDE[4][4] =
{
    {0.8 , 0.0 , -0.6, -0.6 }   , //C
    {0.0 , 0.8 , -0.6, -0.6 }   , //H
    {-0.6, -0.6, 0.6 , -0.6 }   , //D
    {-0.6, -0.6, -0.6, 0.6 } //E
};

struct HomoPro
{
    char idTar[SIZE_CHAIN_ID+1];  /* chain identifier for the target sequence, main sequence */
    char idCan[SIZE_CHAIN_ID+1];  /* chain identifier for the candidate sequence, the one detected */
    int lengthTar;
    int lengthCan;
    int begPosTar; /* starting position for the target sequence alinged with the candidate seq, starting from 0 */
    int endPosTar; /* end position is not include in the sequence, thus [0,12], include 12 residues, if the aaSeqIndex is starting from 0*/
    int begPosCan;
    int endPosCan;
    double homoScore;
    bool isSegInfoSupplied; /*if the aligned part information for the homoPro supplied in the homoProFile*/
    bool isWholeChainAligned; /*if the whole protein should be considered as aligned*/
};

double proportion = 0.45;
double weight = 1.0;
int max_homo_pro = 3; // if there are more than max_homo_pro homoPros, take only the top max_homo_pro
bool isNotUsingSegInfo = false;

double cutoff_consv = 0.7;
double cutoff_score2 = 0.00;
int cutoff_homoScore = 6;
int cutoff_homoScore_2 = cutoff_homoScore ;  /*false positive rate for the number 1 and number 2 are the same, when homoScore is <cutoff_homoScore*/
int min_metalBoundRes = 3;
int max_metalBoundRes = 5;
bool isExcludeOtherChainRes = true;
bool isUsingTotalBoundRes = true;
int homofilename_format = 1; // homofilename format 0: Res_$ID.txt 1: $ID.fragpost

//ChangeLog
/*****************************************************************************
 * ChangeLog 2007-08-21
 *     pattern match score is also used in calculating score for homology-based
 *     prediction
 *     patternMatchScore, 
 *     mainly ranging from -2 to 3
 * ChangeLog 2011-10-10 
 *    Add homologues in the output started with annotation tag
 *    "#"
 ****************************************************************************/

/*****************************************************************************
 * formula 1, convert the gating score to discriminant
 * sigmoid_score = 1 / (1 + exp (-a*x-b))
 * thus
 *
 * x = -(( ln((1-y) /y)) + b)/ a 
 *              1-y
 *         ln (---- ) + b
 *               y
 *  x = - ________________      (formula 1)
 *             a
 ****************************************************************************/

void PrintHelp()
{
    fprintf(stdout,"usage: znpred-postscan [options] znpredfile\n");
    fprintf(stdout," post scan the zinc binding residue prediction\n");
    fprintf(stdout,"options:\n");
    fprintf(stdout," -n #site            : numSite for each sample\n");
    fprintf(stdout," --gate              : the score is gating score, in that case formula 1 is used\n");
    fprintf(stdout," --rule homo|score   : using homology or score for postscan\n");
    fprintf(stdout," --homo <string>     : path for homology matching files\n");
    fprintf(stdout," --modm <string>     : path for modm files\n");
    fprintf(stdout," -s1 score           : the cutoff score for selecting seed Zn-binding residue, default = 0.1\n");
    fprintf(stdout," -s2 score           : the cutoff score for selecting neighboring Zn-binding residue, default = 0.0\n");
    fprintf(stdout," -sa score           : average score for Zn-binding residue group\n");
    fprintf(stdout," -a <string>         : resList for reading MetalProtein, default = CH\n");
    fprintf(stdout," -num number         : minimal number of residues required\n");
    fprintf(stdout," -d dist             : maximum distance in sequence to be bound to the same zinc\n");
    fprintf(stdout," -o file             : output to file outfile, default=stdout\n");
    fprintf(stdout," -w|--weight <real>  : weight on score derived from homo protein versus discriminant, default=0.55\n");
    fprintf(stdout," --prop      <real>  : proportion on score derived from homo protein versus discriminant, default=0.55\n");
    fprintf(stdout," --c-homo    <int>   : cutoff score for percentage of frag match of a homo pro, default=33\n");
    fprintf(stdout," --max-homopro <int> : maximal number of homomogous proteins to be used, default = 3\n");
    fprintf(stdout," --not-use-seg <int> : not using the matched segment information in the detected homologous, proteins, that is, consider the whole protein as homologous\n");
    fprintf(stdout," --homo-format<int>  : format for the homofile name, 0--Res_$ID.txt, 1--$ID.fragpost, default = 0\n");
    fprintf(stdout," --ssbond <string>   : ssbond file, default=$DATADIR/passe.ssbond\n");
    fprintf(stdout," --metal  <string>   : closeMetalPro file, default=$DATADIR/closeMetalPro.dat\n");
    fprintf(stdout," --level  <int>      : level for metal proteins, default = 1\n");
    fprintf(stdout," -h|--help           : print this help message and exit\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"Created on 2007-02-22, updated 2011-10-10, Nanjiang Shu\n");
    //fprintf(stdout,"\n");
    //fprintf(stdout,"\n");
    //fprintf(stdout,"\n");
    //fprintf(stdout,"\n");
}

void PrintVerboseHelp()/*{{{*/
{

}
/*}}}*/
void InitHomoPro(HomoPro *p)/*{{{*/
{
    strcpy(p->idTar, "");
    strcpy(p->idCan, "");
    p->lengthTar = 0;
    p->lengthCan = 0;
    p->begPosTar = 0;
    p->endPosTar = 0;
    p->begPosCan = 0;
    p->endPosCan = 0;
    p->homoScore = 0.0;
    p->isSegInfoSupplied = false;
    p->isWholeChainAligned = false;
};/*}}}*/
void CopyHomoPro(HomoPro *to, HomoPro *from)/*{{{*/
{
    my_strcpy(to->idTar, from->idTar, SIZE_CHAIN_ID);
    my_strcpy(to->idCan, from->idCan, SIZE_CHAIN_ID);
    to->lengthTar = from -> lengthTar;
    to->lengthCan = from -> lengthCan;
    to->begPosTar = from -> begPosTar;
    to->endPosTar = from -> endPosTar;
    to->begPosCan = from -> begPosCan;
    to->endPosCan = from -> endPosCan;
    to->homoScore = from -> homoScore;
    to->isSegInfoSupplied = from -> isSegInfoSupplied;
    to->isWholeChainAligned = from -> isWholeChainAligned;
};/*}}}*/
int GetHomoPro(const char *homofile, HomoPro *homoPro, int max_homo_pro)/*{{{*/
{
    FILE *fpHomo = NULL;
    fpHomo = fopen(homofile, "r");
    int numHomoPro = 0;
    int idxHomoPro = 0; /*index used when searching in homoProIDList*/
    
    Array2D <char> homoProIDList_array2d(max_homo_pro, SIZE_CHAIN_ID+1);/*added 2008-12-08 */
    char ** homoProIDList = homoProIDList_array2d.array2D; /*this list is stored as the predicted homopro ids are recorded, it is neither unique nor sorted*/

    if(fpHomo != NULL)
    {
        int linesize;
        int maxline = 300;
        Array1D <char> line_1darray(maxline+1);
        char *line = line_1darray.array1D;
        //read in homoprotiens
        int cntHomoPro = 0;
        int num;
        char idTar[SIZE_CHAIN_ID+1] = "";  /* chain identifier for the target sequence, main sequence */
        char idCan[SIZE_CHAIN_ID+1] = "";  /* chain identifier for the candidate sequence, the one detected */
        int lengthTar = 0;
        int lengthCan = 0;
        int begPosTar = 0; /* starting position for the target sequence alinged with the candidate seq, starting from 0 */
        int endPosTar = 0; /* end position is not include in the sequence, thus [0,12], include 12 residues, if the aaSeqIndex is starting from 0*/
        int begPosCan = 0;
        int endPosCan = 0;
        double homoScore = 0.0;
        double maxHomoScore = 0.0;
        int status_sscanf = 0;

        while((linesize= fgetline(fpHomo, line, maxline)) != EOF)
        {
            if(strncmp(line,"Num",3) == 0) break;
        }

        while((linesize= fgetline(fpHomo, line, maxline)) != EOF)
        {
            if(linesize <= 0 || line[0] == '#') continue;
            status_sscanf = sscanf(line,"%d %s %lf %s %d %d %d %s %d %d %d", 
                    &num, idCan, &homoScore,
                    idTar, &lengthTar, &begPosTar, &endPosTar,
                    idCan, &lengthCan, &begPosCan, &endPosCan);
            StdID(idTar);
            StdID(idCan);
            if(cntHomoPro == 0)
            {
                maxHomoScore = homoScore; /*since homoScore ordered descendingly, the first homoScore is the maxHomoScore*/
            }

            if(homoScore >= cutoff_homoScore
                    || (maxHomoScore < cutoff_homoScore && homoScore >= cutoff_homoScore_2)) 
            {
                if(status_sscanf >=2)
                {
                    my_strcpy(homoPro[cntHomoPro].idCan, idCan, SIZE_CHAIN_ID);
                }
                if(status_sscanf >=3)
                {
                    homoPro[cntHomoPro].homoScore = homoScore;
                }
                if(status_sscanf >=4)
                {
                    my_strcpy(homoPro[cntHomoPro].idTar, idTar, SIZE_CHAIN_ID);
                    homoPro[cntHomoPro].lengthTar = lengthTar;
                    homoPro[cntHomoPro].begPosTar = begPosTar;
                    homoPro[cntHomoPro].endPosTar = endPosTar;
                    my_strcpy(homoPro[cntHomoPro].idCan, idCan, SIZE_CHAIN_ID);
                    homoPro[cntHomoPro].lengthCan = lengthCan;
                    homoPro[cntHomoPro].begPosCan = begPosCan;
                    homoPro[cntHomoPro].endPosCan = endPosCan;
                    homoPro[cntHomoPro].isSegInfoSupplied = true;

                    double avgLen = double (lengthTar+lengthCan) /2.0;
                    double lengthProportion = double (min(lengthCan, lengthTar)) / max(lengthTar, lengthCan);
                    if( lengthCan < 200 && lengthTar < 200
                            && lengthProportion < 0.7
                            && abs(begPosTar-begPosCan) / avgLen < 0.1
                            && Coverage(begPosTar, endPosTar, begPosCan, endPosCan) / avgLen >= 0.5)
                    {
                        homoPro[cntHomoPro].isWholeChainAligned = true;
                    }
                    else
                    {
                        homoPro[cntHomoPro].isWholeChainAligned = false;
                    }

                    if(isNotUsingSegInfo)
                    {
                        homoPro[cntHomoPro].isWholeChainAligned = true;
                        homoPro[cntHomoPro].isSegInfoSupplied = false;
                    }
                }
                else
                {
                    homoPro[cntHomoPro].isSegInfoSupplied = false;
                    homoPro[cntHomoPro].isWholeChainAligned = true;
                }

                if ((idxHomoPro = LinearSearch_String(idCan, homoProIDList, cntHomoPro)) != -1 ) /*use only one copy of the homoPro*/
                {
                }
                else
                {
                    my_strcpy(homoProIDList[cntHomoPro], idCan, SIZE_CHAIN_ID);
                    cntHomoPro ++;
                }
            }
            if(cntHomoPro >= max_homo_pro)
            {
                break;
            }
        }
        fclose(fpHomo);
        numHomoPro = cntHomoPro ; 

    }
    else
    {
        numHomoPro = -1;
    }
    return numHomoPro;
}/*}}}*/
double SigmoidScore(double a, double b, double x)/*{{{*/
{
    return 1.0/(1.0+exp(-a*x-b));
}
/*}}}*/
int	NormalizeZnPredScore(double *znscore, int num)/*{{{*/
{
	int i;
	for( i = 0 ; i < num ; i ++)
	{
		if(znscore[i] < -0.2)
		{
			znscore[i] = SigmoidScore(10, 0.2, znscore[i]);
		}
		else if(znscore[i] < 0.0)
		{
			znscore[i] = SigmoidScore(10, 0.1, znscore[i])- 0.02;
		}
		else if(znscore[i] < 0.5)
		{
			znscore[i] = SigmoidScore(6, 0.2, znscore[i]) - 0.08;
		}
		else
		{
			znscore[i] = SigmoidScore(3.5, 0.2, znscore[i]);
		}
	}
	return num;
}/*}}}*/

char* GetHOMOFilePath(const char *rmtID, char *homofile, const char *homopath, int homofilename_format)/*{{{*/
{
    if(homofilename_format == 0)
    {
        sprintf(homofile, "%s/Res_%s.txt", homopath, rmtID);
    }
    else// if(homofilename_format == 1)
    {
        sprintf(homofile, "%s/%s.fragpost", homopath, rmtID);
    }

    return homofile;

}
/*}}}*/
double AntiSigmoidScore(double y, double a , double b)/*{{{*/
{
    return ( -(log ((1-y)/y)+b)/a );
}/*}}}*/
double GetScore(GistPredChain *pChain, int idx_seed, int *idx_neigh, int num_neigh)/*{{{*/
    // get the score of the residue group based on the number distance to the
    // seed residue and the discriminant values.
{
    double a = 0.1;
    double b = 2.0;
    double sigmoid_score = 0.0;
    double gating_score = 0.0;
    int dist = 0;
    double avg_score = 0.0;
    int i ;
    for(i = 0 ; i < num_neigh; i ++)
    {
        dist = abs(pChain->aaSeqIndex[idx_neigh[i]] - pChain->aaSeqIndex[idx_seed]);
        sigmoid_score = SigmoidScore(a,b,double(dist));
        gating_score = GatingScore(sigmoid_score, pChain->discriminant[idx_neigh[i]]);
        avg_score += gating_score;
    }
    avg_score /= num_neigh;
    return avg_score;
}
/*}}}*/
int select_best_binding_group(GistPredChain *pChain, int *idx, int num_local_HRes, int num_binding_group)/*{{{*/
/*****************************************************************************
 * giving a number of potential zinc binding residues, selec the best zinc
 * binding residue group
 * based on the discriminant value and the distance to the seed residue
 ****************************************************************************/
{
    int i, jj;
    int countSubset = 0;
    //int numSubset = 0;
    int nSet = num_local_HRes - 1;// because idx[0] is for the seed residue, and must be included, then the size of the set should be minus by 1
    int nSub = 0;
    Array1D<int> idxtmp_1darray(nSet);
    int *idxtmp = idxtmp_1darray.array1D;
    int idx_max_score;
    double max_score;
    double score;


    if( num_local_HRes >= 3 )
    {
        if(num_local_HRes >= 4) nSub = 3;
        else if(num_local_HRes >= 3) nSub = 2;

        countSubset = combin2(nSet,nSub);
        Array2D <int> subset_2darray(countSubset,nSub);
        int **subset = subset_2darray.array2D;
        int numSubset = GetAllSubset(subset, nSet, nSub);
        max_score = -9999.0;
        idx_max_score = -1;
        for ( i = 0 ; i < numSubset; i ++)
        {
            for(jj = 0 ; jj < nSub; jj ++) idxtmp[jj] = idx[subset[i][jj] +1];
            score = GetScore(pChain, idx[0], idxtmp, nSub);
            if(score > max_score)
            {
                max_score = score;
                idx_max_score = i;
            }
        }
        //assign the label for the group with the max score
        pChain->label[idx[0]] = num_binding_group;
        for( jj = 0 ; jj < nSub; jj ++)
            pChain->label[ idx[subset[idx_max_score][jj]+1] ] = num_binding_group;
        return (nSub + 1);
    }
    else if (num_local_HRes ==2 )
    {
        if( abs (pChain->aaSeqIndex[idx[0]] - pChain->aaSeqIndex[idx[1]]) <= 10)
        {
            for( jj = 0 ; jj < nSub; jj ++)
                pChain->label[ idx[jj] ] = num_binding_group;
            return num_local_HRes;
        }
    }

    return 0;
}
/*}}}*/
int ScanZnPred(GistPredChain *pChain, double cutoff_discriminant_1, double cutoff_discriminant_2, int cutoff_window, int min_local_numHRes, double cutoff_avg_score)/*{{{*/
/*****************************************************************************
 * cutoff_discriminant_1 : cutoff of the discriminant value for selecting the seed of ZB residues
 * cutoff_discriminant_2 : cutoff of the discriminant when extending the zinc binding residue group
 * cutoff_window         : the maximum distance for residues in sequence to be bound by the same zinc atom
 * min_local_numHRes     : minimal number of residues with discriminant >=
 *                         cutoff_discriminant needed within the cutoff_window in order to be
 *                         considered as zinc binding
 * cutoff_avg_score      : minimal score for the average discriminant of a group of residues 
 *                         within cutoff_window in order to be considered as zinc binding
 *
 * the mechanism for predicting a residue with low discriminant value as zinc
 * binding has not been incorporated yet. When predicting proteins of Jonas, I
 * found it might be valuable to predict one residue near two highly predicted
 * zinc binding residues as zinc binding
 ****************************************************************************/

{
    // first find all residues with discriminant above cutoff_discriminant
    int i,j;
    Array1D <int> idx_seed_1darray(pChain->numRes);
    Array1D <int> idx_local_1darray(pChain->numRes);
    int *idx_seed = idx_seed_1darray.array1D; // index of the residues in pChain [0 .. pChain->numRes] with discriminant >= cutoff_discriminant_1
    int *idx_local = idx_local_1darray.array1D; // index of residues in pChain [0 .. pChain->numRes] for possible zinc binding residue group near seed residue
    int numHRes = 0;
    for( i = 0 ; i < pChain->numRes; i++)
    {
        pChain->label[i] = -1; // initialize the label to -1
        if(pChain->discriminant[i] >= cutoff_discriminant_1)
        {
            idx_seed[numHRes] = i;
            numHRes ++;
        }
    }
    // for each seed residue , check if there are >= min_local_numHRes within cutoff_window
    // satisfying the rule 
    // The procedure is heuristic, that is, the residues has already been
    // marked as zinc binding will not be check any more
    int num_binding_group = 1;
    for ( i = 0 ; i < numHRes; i ++)
    {
        if(pChain->label[idx_seed[i]] > 0) // if the residue has already been marked as zinc binding, ignore checking
            continue;

        int num_local_HRes = 0; //local HRes within 2*cutoff_window
        idx_local[num_local_HRes] = idx_seed[i];
        num_local_HRes ++;
        for ( j = idx_seed[i]-1 ; j >= 0; j --) // searching backward
        {
            if((pChain->aaSeqIndex[idx_seed[i]] - pChain->aaSeqIndex[j]) > cutoff_window)
                break;

            if(pChain->label[j] > 0) // if the residue has already been marked as zinc binding, ignore checking 
                continue;

            if(pChain->discriminant[j] >= cutoff_discriminant_2)
            {
                idx_local[num_local_HRes] = j;
                num_local_HRes ++;
            }
        }

        for ( j = idx_seed[i]+1 ; j < pChain->numRes; j ++) // searching forward
        {
            if((pChain->aaSeqIndex[j] - pChain->aaSeqIndex[idx_seed[i]]) > cutoff_window)
                break;

            if(pChain->label[j] > 0) // if the residue has already been marked as zinc binding, ignore checking 
                continue;

            if(pChain->discriminant[j] >= cutoff_discriminant_2)
            {
                idx_local[num_local_HRes] = j;
                num_local_HRes ++;
            }
        }

        if(num_local_HRes < min_local_numHRes)
            continue;
        else
        {
            if(select_best_binding_group(pChain, idx_local, num_local_HRes, num_binding_group) > 0 )// if there is a zinc binding group predicted, increment num_binding_group
                num_binding_group ++;
        }
    }

    int num_zb_res = 0;
    for( i = 0 ; i < pChain->numRes ; i++)
    {
        if ( pChain->label [i] > 0) num_zb_res ++;
    }
    return num_zb_res;
}
/*}}}*/
int SearchPattern(GistPredChain *pChain, MetalPro *pMetalPro, int idxAtomEnv, HomoPro *pHomoPro, int *idxPatternRes, int &numPatternRes, double &patternMatchScore)/*{{{*/
/*****************************************************************************
 * for a given atomEnv, search the pattern of atomEnv in predChain, output the
 * index of residues with matched pattern in idxPatternRes 2007-06-11
 * adding using the aligned part of the homoPro, 2007-07-19
 ****************************************************************************/
{
    int i,j;
    //double lengthdiff = abs(pMetalPro->length-pChain->length) / double(pMetalPro->length+pChain->length) *2  ;
    double avgLen = double (pMetalPro->length + pChain->length) /2.0;


    int max_dest_num_res = pMetalPro->atomEnv[idxAtomEnv].numRes ;
    Array1D <char> destAA_1darray(max_dest_num_res+1);
    Array1D <int> destAASeqIndex_1darray(max_dest_num_res);
    char *destAA = destAA_1darray.array1D;
    int *destAASeqIndex = destAASeqIndex_1darray.array1D;
    int aaSeqIndex = 0;
    char aa = ' ';

    int cntDestNumRes = 0;

    for(i = 0; i < max_dest_num_res; i ++)
    {
        aa = pMetalPro->res[pMetalPro->atomEnv[idxAtomEnv].parentResIndex[i]].aa; 
        aaSeqIndex =  pMetalPro->res[pMetalPro->atomEnv[idxAtomEnv].parentResIndex[i]].aaSeqIndex; 
        if( pHomoPro->isWholeChainAligned ||
               (aaSeqIndex >= pHomoPro->begPosCan && aaSeqIndex < pHomoPro->endPosCan))
        {
            destAA[cntDestNumRes] = aa;
            destAASeqIndex[cntDestNumRes] = aaSeqIndex;
            cntDestNumRes ++;
        }
    }
    destAA[cntDestNumRes] = '\0';
    int destNumRes = cntDestNumRes;



    if( destNumRes < 1) /*if destNumRes == 0, no predicted residue should be used*/
    {
        return 0;
    }                   
    /*pattern searching is not suitable for <2 residue pattern */ 

    /*selecting the residues within the aligned region for the target sequence*/ 
    Array1D <char> aaTar_1darray(pChain->numRes);
    Array1D <int>  aaSeqIndexTar_1darray(pChain->numRes);
    Array1D <int>  idxPChainTar_1darray(pChain->numRes);
    char *aaTar = aaTar_1darray.array1D; /*aa for residues within alinged region of resides in the pChain*/
    int *aaSeqIndexTar = aaSeqIndexTar_1darray.array1D;/*aaSeqIndex for residues within alinged region of residues in the pChain */
    int *idxPChainTar = idxPChainTar_1darray.array1D;/*index to the above residue in pChain*/
    /* for example, 
     * aaTar[0] = 'H'
     * aaSeqIndexTar[0] = 79
     * idxPChainTar[0] = 0
     * means residue H is the 79th residue in the chain, 
     * but the 0th residue for the list of predicted zn-binding residues
     * */

    int cntResTar = 0;
    for(j = 0 ; j < pChain->numRes; j ++)
    {
        if(pHomoPro->isWholeChainAligned ||
                (pChain->aaSeqIndex[j] >= pHomoPro->begPosTar && pChain->aaSeqIndex[j] < pHomoPro->endPosTar))
        {
            aaTar[cntResTar] = pChain->aaSeq[j];
            aaSeqIndexTar[cntResTar] = pChain->aaSeqIndex[j];
            idxPChainTar[cntResTar] = j;
            cntResTar ++;
        }
    }
    int numResTar = cntResTar;/*number of predicted residues within the aligned region*/

    if(numResTar <= destNumRes) /*if number of residues within the aligned region is less than destNumRes, select all of them in the pattern*/
    {
        for(j = 0; j < numResTar; j ++)
        {
            idxPatternRes[j] = idxPChainTar[j];
        }
        numPatternRes = numResTar;
    }
    else
    {
        
        int countSubset = combin2(numResTar, destNumRes);
        Array2D <int> subset_2darray(countSubset,destNumRes);
        int **subset = subset_2darray.array2D;
        Array1D <double> matchScore_1darray(countSubset);
        double *matchScore = matchScore_1darray.array1D;

        GetAllSubset(subset, numResTar, destNumRes);
        
        for(i = 0 ;  i < countSubset; i ++)
        {
            double score = 0.0;
            int daa1 = 0;
            int daa2 = 0;
            for(j = 0 ; j < destNumRes; j++)
            {
                daa1 = Char2Digit(aaTar[subset[i][j]], CHDE_Alphabet);
                daa2 = Char2Digit(destAA[j], CHDE_Alphabet);
                score += SubMatrixCHDE[daa1][daa2];

                /*assuming if the shift of residue position < 1/4 avgLen, the
                 * score is still positive */
                score += ( - (abs((aaSeqIndexTar[subset[i][j]] - pHomoPro->begPosTar)- (destAASeqIndex[j]-pHomoPro->begPosCan))/avgLen)); 

                if(j > 0)
                {
                    score +=(- (abs( 
                            (aaSeqIndexTar[subset[i][j]] - aaSeqIndexTar[subset[i][j-1]]) -
                            (destAASeqIndex[j] - destAASeqIndex[j-1]) 
                            ) /( double (destAASeqIndex[j] - destAASeqIndex[j-1]) + 1.0)));
                    /*2 tourate, if the destpattern is 
                     * 12 14 25 79
                     * and target is 12 18 25 91
                     *
                     * 14-12 = 2, then 2*2 = 4, so if the first residue and
                     * second residue in the target pattern is less than 4, the
                     * score is still positive*/
                }
            }
            matchScore[i] =(score) ; /*divided by 10 to reduced the score to around -1~1*/
        }
        int index_max_score = max_element_index (matchScore, 0, countSubset-1);
        patternMatchScore = matchScore[index_max_score]; /*2007-08-21, return patternMatchScore as well*/

        for(j = 0 ; j < destNumRes; j ++)
        {
            idxPatternRes[j] = idxPChainTar[subset[index_max_score][j]];
        }
        numPatternRes = destNumRes;
    }

//      if target sequence and candidate are well aligned,
//      C10 C14 C54                     #target sequence,
//                C45   C63  C89 C102   #ZN binding, candidate sequence
//
//      C10, C14, C54 in the target sequence is unlikely to be matched with
//      zinc-binding residues in the candidate sequence          
    int cntResWithinRange = 0;
    for(j = 0 ; j < numPatternRes; j ++)
    {
        int k ;
        bool isWithinRange = false;
        for(k = 0 ; k < destNumRes; k ++)
        {
            if( abs(
                        (pChain->aaSeqIndex[idxPatternRes[j]]-pHomoPro->begPosTar) -
                        (destAASeqIndex[k] - pHomoPro->begPosCan)
                   ) /double(pChain->length) <0.1)
            {
                isWithinRange = true;
                break;
            }
        }
        cntResWithinRange += isWithinRange;
    }
    if(cntResWithinRange < 2 || patternMatchScore < 0.0) /*use only when the patternMatchScore is positive, 2008-12-08*/
    {
        numPatternRes = 0;
    }

    return numPatternRes;
}/*}}}*/
int ScanZnPred_HOMO(GistPredChain *pChain, MetalPro* metalPro, char **metalProIDList, int numMetalPro, SSBondPro *ssbondPro, char **ssbondProIDList, int numSSBondPro, const char *homopath, double cutoff_discriminant, double *znScoreFromHomo , FILE* fpout)/*{{{*/
{
    int i,j,k ;
    int numZnRes = 0;
    char homofile[MAX_PATH+1] = "";
    char rmtID[SIZE_CHAIN_ID+1] = "";
    my_strcpy(rmtID, pChain->id, SIZE_CHAIN_ID);
    RemoveTID(rmtID);
//    sprintf(homofile, "%s/Res_%s.txt", homopath, rmtID);
    GetHOMOFilePath(rmtID, homofile, homopath, homofilename_format);
    int numHomoPro = 0;

    /*the metalProIDList should be sorted first in order to use
     * BinarySearch_String, this bug has been fixed on 2008-12-07*/
    Array2D <char> metalProIDList_sorted_2darray(numMetalPro,SIZE_CHAIN_ID+1);
    char **metalProIDList_sorted = metalProIDList_sorted_2darray.array2D;
    Array1D <int> idxSortMetalPro_1darray(numMetalPro);
    int *idxSortMetalPro = idxSortMetalPro_1darray.array1D;
    for(i = 0 ; i < numMetalPro ; i++) idxSortMetalPro[i] = i;
    for(i = 0 ; i < numMetalPro; i++) { my_strcpy(metalProIDList[i], metalPro[i].id, SIZE_CHAIN_ID); }
    QuickSort_String(idxSortMetalPro, metalProIDList, 0, numMetalPro-1);
    for(i = 0 ; i < numMetalPro; i ++) { my_strcpy(metalProIDList_sorted[i], metalProIDList[idxSortMetalPro[i]], SIZE_CHAIN_ID);}

    Array2D <char> ssbondProIDList_sorted_2darray(numSSBondPro,SIZE_CHAIN_ID+1);
    char **ssbondProIDList_sorted = ssbondProIDList_sorted_2darray.array2D;
    Array1D <int> idxSortSSBondPro_1darray(numSSBondPro);
    int *idxSortSSBondPro= idxSortSSBondPro_1darray.array1D;
    for(i = 0 ; i < numSSBondPro ; i++) idxSortSSBondPro[i] = i;
    for(i = 0 ; i < numSSBondPro; i++) { my_strcpy(ssbondProIDList[i], ssbondPro[i].id, SIZE_CHAIN_ID); }
    QuickSort_String(idxSortSSBondPro, ssbondProIDList, 0, numSSBondPro-1);
    for(i = 0 ; i < numSSBondPro; i ++) { my_strcpy(ssbondProIDList_sorted[i], ssbondProIDList[idxSortSSBondPro[i]], SIZE_CHAIN_ID);}
    /*---------------------*/



    Array1D <HomoPro> homoPro_1darray(max_homo_pro+1);
    HomoPro *homoPro = homoPro_1darray.array1D;
    for(i = 0 ; i < max_homo_pro+1; i ++) { InitHomoPro(&homoPro[i]); }

    numHomoPro = GetHomoPro(homofile, homoPro, max_homo_pro);
    fprintf(stdout,"homofile=%s",homofile);

    if(numHomoPro > 0)
    {
        fprintf(fpout,"#Homologues used in prediction of %s:\n", pChain->id);
        for (i=0;i<numHomoPro;i++)
        {
            fprintf(fpout,"#Homolog %d for %s: %s %.3lf\n", i+1, pChain->id, homoPro[i].idCan, homoPro[i].homoScore);
        }
        fprintf(fpout,"#\n");

        double scale_znScore = 1.0/pow(double(numHomoPro), 0.6);

        Array2D <double> znScore_2darray(max_homo_pro, pChain->numRes);
        double **znScore = znScore_2darray.array2D;
        //initialize znScore
        for(i = 0; i < max_homo_pro; i++)
        { for(j = 0 ; j < pChain->numRes; j ++)
            { znScore[i][j] = 0.0; }
        }
        
        char keyMetal[] = "ZN";
        const char *alterMetal[] = {"CD", "NI"};
        int numAlterMetal = sizeof(alterMetal) / sizeof(char*);

        const char *alterMetal2[] = {"FE", "CU"};
        int numAlterMetal2 = sizeof(alterMetal2) / sizeof(char*);

        Array1D <int> idxPatternRes_1darray(pChain->numRes);
        int *idxPatternRes = idxPatternRes_1darray.array1D; //index array for pattern searching
        int numPatternRes = 0 ;
        double patternMatchScore = 0.0;
        int iHomo = 0;

        bool isAllAlterMetalBinding  = false;
        bool isAllAlterMetalBinding2 = false;

        int cntAlterMetalBindingPro = 0;
        int cntAlterMetalBindingPro2 = 0;
        for(iHomo = 0 ; iHomo < numHomoPro; iHomo++)
        {   //global check if all homo protein is alter2 metal binding protein
            int indexMetalPro = 0;
            bool isAlterMetalBinding = false;
            bool isAlterMetalBinding2 = false;

            MetalPro *pMetalPro;
            if ( (indexMetalPro = BinarySearch_String( homoPro[iHomo].idCan, metalProIDList_sorted, numMetalPro)) != -1) //if the homoPro a metal binding protein
            {
                indexMetalPro = idxSortMetalPro[indexMetalPro]; //map the index to the unsorted array
                pMetalPro = &(metalPro[indexMetalPro]);
                for(i = 0 ; i < pMetalPro->numMetalAtom; i++)
                {
                    if (LinearSearch_String((const char*)pMetalPro->atomEnv[i].metalAtomName, alterMetal, numAlterMetal) != -1)// if one protein is not all alterMetal binding, set isAllAlterMetalBinding to false
                    {
                        isAlterMetalBinding = true;
                    }
                    else if (LinearSearch_String((const char*)pMetalPro->atomEnv[i].metalAtomName, alterMetal2, numAlterMetal2) != -1)
                    {
                        isAlterMetalBinding2 = true;
                    }

                    cntAlterMetalBindingPro += isAlterMetalBinding;
                    cntAlterMetalBindingPro2 += isAlterMetalBinding2;
                }
            }
        }
        if (cntAlterMetalBindingPro / double(numHomoPro) > 0.5)
        {
            isAllAlterMetalBinding = true;
        }
        if (cntAlterMetalBindingPro2 / double(numHomoPro) > 0.5)
        {
            isAllAlterMetalBinding2 = true;
        }


        /*calculate the znScore of each residue according to the information
         * of the predicted homoPro*/

        for(iHomo = 0 ; iHomo < numHomoPro; iHomo++)
        {
            int indexMetalPro = 0;
            int indexSSBondPro = 0;
            HomoPro *pHomoPro = &(homoPro[iHomo]);
            
            MetalPro *pMetalPro;
            SSBondPro *pSSBondPro; 
            bool isZnBound = false;
            bool isAlterMetalBound = false;
            bool isAlterMetalBound2 = false;


            if((indexMetalPro = BinarySearch_String( homoPro[iHomo].idCan, metalProIDList_sorted, numMetalPro))!= -1) //if the homoPro a metal binding protein
            {
                indexMetalPro = idxSortMetalPro[indexMetalPro]; //map the index to the unsorted array
                pMetalPro = &(metalPro[indexMetalPro]);
                if(pChain->length >= 20)
                {
                    bool isAllZn = true;
                    for(i = 0 ; i < pMetalPro->numMetalAtom; i++)
                    {
                        if(strcasecmp(pMetalPro->atomEnv[i].metalAtomName, keyMetal) != 0 ||
                                LinearSearch_String((const char*)pMetalPro->atomEnv[i].metalAtomName, alterMetal, numAlterMetal) != -1)
                        {
                            isAllZn = false;
                        }
                    }

                    //if(pMetalPro->numBoundRes == pChain->numRes && homoScore[iHomo] >= 40 && isAllZn)
                    //{
                    //    SearchPattern(pChain, pMetalPro, i, idxPatternRes, numPatternRes);
                    //    for(j = 0 ; j < numPatternRes; j ++)
                    //    {
                    //        znScore[iHomo][idxPatternRes[j]] = min(homoScore[iHomo] / 40.0 /numHomoPro, 1.0 / numHomoPro);
                    //    }

                    //    //for(j = 0 ; j < pChain->numRes ; j ++)
                    //    //{
                    //    //    znScore[iHomo][j] =  min(homoScore[iHomo] / 40.0 / numHomoPro, 1.5/ numHomoPro); 
                    //    //}
                    //    isZnBound = true;
                    //}
                    //else
                    {
                        for(i = 0 ; i < pMetalPro->numMetalAtom; i++)
                        {
                            if(strcasecmp(pMetalPro->atomEnv[i].metalAtomName, keyMetal) == 0)
                            {
                                //searching for pattern
                                SearchPattern(pChain, pMetalPro, i, &homoPro[iHomo], idxPatternRes, numPatternRes, patternMatchScore);
#ifdef DEBUG_SEARCHPATTERN
                                printf("SearchPattern Debug, search pattern for \n");
                                printf("ID = %s\n", pChain->id);
                                printf("Metal AtomEnv for the Homo Protein %s is \n", pMetalPro->id);
                                for(int jj = 0 ; jj < pMetalPro->atomEnv[i].numRes; jj ++)
                                {
                                    printf("%c %d\n", pMetalPro->res[pMetalPro->atomEnv[i].parentResIndex[jj]].aa, pMetalPro->res[pMetalPro->atomEnv[i].parentResIndex[jj]].aaSeqIndex+1);
                                }
                                printf("Matched pattern is \n");
                                for(int jj = 0; jj < numPatternRes; jj ++)
                                {
                                    printf("%c %d\n", pChain->aaSeq[idxPatternRes[jj]], pChain->aaSeqIndex[idxPatternRes[jj]]+1);
                                }
#endif
                                for(j = 0 ; j < numPatternRes; j ++)
                                {
                                    znScore[iHomo][idxPatternRes[j]] = min(pHomoPro->homoScore / 50.0*scale_znScore * pow(1.35, patternMatchScore), 1.5*scale_znScore);
                                }
                                isZnBound = true;
                            }
                            else if (LinearSearch_String((const char*)pMetalPro->atomEnv[i].metalAtomName, alterMetal, numAlterMetal) != -1)
                            {
                                SearchPattern(pChain, pMetalPro, i, pHomoPro,idxPatternRes, numPatternRes, patternMatchScore);
                                for(j = 0 ; j < numPatternRes; j ++)
                                {
                                    znScore[iHomo][idxPatternRes[j]] = min(pHomoPro->homoScore/ 50.0 *scale_znScore* pow(1.35, patternMatchScore), 1.2 * scale_znScore);
                                }
                                isAlterMetalBound = true;
                            }
                            else if (LinearSearch_String((const char*)pMetalPro->atomEnv[i].metalAtomName, alterMetal2, numAlterMetal2) != -1 &&
                                    pMetalPro->atomEnv[i].totalBoundRes  == 4
                                    && !isAllAlterMetalBinding2)
                            {
                                SearchPattern(pChain, pMetalPro, i, pHomoPro, idxPatternRes, numPatternRes, patternMatchScore);
                                for(j = 0 ; j < numPatternRes; j ++)
                                {
                                    znScore[iHomo][idxPatternRes[j]] = min(pHomoPro->homoScore/ 50.0*scale_znScore* pow(1.35, patternMatchScore), 1.0 * scale_znScore);
                                }
                                isAlterMetalBound2 = true;
                            }

                            if(numPatternRes > 0)
                            {
                                fprintf(stdout,"patternMatchScore = %.4lf\n", patternMatchScore);
                            }

                        }
                    }

                    if(!isZnBound && !isAlterMetalBound &&  !isAlterMetalBound2) // 2007-07-13, bug fixed, isAllZn should be !isAllZn
                    {
                        /*give penalty to the residues within the aligned region,
                         * but not zinc binding*/
                        for(j = 0 ; j < pChain->numRes ; j ++)
                        {
                            if(znScore[iHomo][j] <= 0.0) /*znScore <= 0.0, means that residue has not been predicted as zinc-binding*/
                            {
                                if(pHomoPro->isWholeChainAligned || 
                                        (pChain->aaSeqIndex[j] >= pHomoPro->begPosTar && pChain->aaSeqIndex[j] < pHomoPro->endPosTar))
                                {

                                    znScore[iHomo][j] = max (-pHomoPro->homoScore/ 50.0*scale_znScore, -1.2  * scale_znScore);
                                }
                            }
                        }
                    }
                    
//                     if(homoScore[iHomo] >= 50)//highly confident homology/*{{{*/
//                     {
//                         if(isZnBound)
//                         {
//                             int cntZnBoundRes = 0 ;
//                             int j;
//                             for(i = 0 ; i < pMetalPro->numBoundRes; i++)
//                             {
//                                 for(j = 0 ; j < pMetalPro->res[i].numMetalBound; j++)
//                                 {
//                                     if(strcmp (pMetalPro->atomEnv[pMetalPro->res[i].parentAtomEnvIndex[j]].metalAtomName, keyMetal) == 0)
//                                     {
//                                         cntZnBoundRes ++; break;
//                                     }
//                                 }
//                             }
// 
//                             //sort pChain->discriminant by DESCENDING order
//                             int idx[20];
//                             for(j = 0 ; j < 20; j ++) idx[j] = j;
//                             QuickSort_index(idx, pChain->discriminant, pChain->numRes, DESCENDING);
//                             for(i = 0 ; i < cntZnBoundRes; i ++)
//                             {
//                                 if(pChain->discriminant[idx[i]] > -0.3)
//                                     pChain->label[idx[i]] = 1;
//                             }
//                         }
//                         else
//                         {
//                             for(i = 0 ; i < pChain->numRes; i++)
//                             {
//                                 pChain->label[i] = -1;
//                             }
//                         }
//                     }/*}}}*/
//                     else //if only moderate confident homology/*{{{*/
//                     {
//                         if(isZnBound)
//                         {
//                             for(i = 0 ; i < pChain->numRes; i++)
//                             {
//                                 if(pChain->discriminant[i] >= cutoff_discriminant)
//                                     pChain->label[i] = 1;
//                             }
//                         }
//                         else
//                         {
//                             for(i = 0 ; i < pChain->numRes; i++)
//                             {
//                                 if(pChain->discriminant[i] <= cutoff_discriminant + 0.05)
//                                     pChain->label[i] = -1;
//                             }
//                         }
// 
//                     }/*}}}*/
                }
            }

            if((indexSSBondPro = BinarySearch_String( homoPro[iHomo].idCan, ssbondProIDList_sorted, numSSBondPro))!= -1 )//if the homoPro a SSBondProtein
            {
                indexSSBondPro = idxSortSSBondPro[indexSSBondPro];
                pSSBondPro = &(ssbondPro[indexSSBondPro]);
                Array1D <int> ssbondResSeqIndex_1darray(pSSBondPro->numSSBondRes+10);
                int *ssbondResSeqIndex = ssbondResSeqIndex_1darray.array1D;
                int cntSSBondRes = 0;
                for(j = 0; j < pSSBondPro->numSSBond; j ++)
                {
                    for(k = 0 ; k < pSSBondPro->ssbond[j].numRes; k ++)
                    {
                        ssbondResSeqIndex[cntSSBondRes] = pSSBondPro->ssbond[j].res[k].aaSeqIndex;
                        cntSSBondRes ++;
                    }
                }
                int numSSBondRes = cntSSBondRes;
                
                for(j = 0 ; j < pChain->numRes; j++)
                {
//                    if(pHomoPro->isWholeChainAligned || (pChain->aaSeqIndex[j] >= pHomoPro->begPosTar && pChain->aaSeqIndex[j] < pHomoPro->endPosTar))
                    {
                        bool isWithinRange = false;
                        for(k = 0 ; k < numSSBondRes; k ++)
                        {
                            if( abs(
                                        (pChain->aaSeqIndex[j]-pHomoPro->begPosTar) -
                                        (ssbondResSeqIndex[k] - pHomoPro->begPosCan)
                                   ) /double(pChain->length) < 0.05)
                            {
                                isWithinRange = true;
                                break;
                            }
                        }
                        if(isWithinRange && pChain->aaSeq[j] == 'C' && znScore[iHomo][j] <= 0.0)
                        {
                            znScore [iHomo][j] += max (-pHomoPro->homoScore/ 50.0 *scale_znScore, -1.5*scale_znScore); 
                        }
                        else if(znScore[iHomo][j] <= 0.0) /*if this residue has been matched with zinc binding residues*/
                        {
                            znScore [iHomo][j] += max (-pHomoPro->homoScore/ 50.0 *scale_znScore, -1.0 *scale_znScore); 
                        }
                    }
                }
            }

            if(indexSSBondPro == -1 && indexMetalPro == -1)
            {
                for(j = 0 ; j < pChain->numRes; j++)
                {
                    if(pHomoPro->isWholeChainAligned || 
                            (pChain->aaSeqIndex[j] >= pHomoPro->begPosTar && pChain->aaSeqIndex[j] < pHomoPro->endPosTar))
                    {

                    znScore [iHomo][j] = max (-pHomoPro->homoScore/ 50.0 *scale_znScore, -1.5*scale_znScore);
                    }
                }
            }
        }

        // network evaluation of the zn-binding state based on znScore[][]
        for(j = 0 ; j < pChain->numRes; j ++)
        {
            double score = 0.0;
            //calculate znScore for each candidate zinc-binding residue
            for(i = 0 ; i < numHomoPro; i ++)
            {
                score += znScore[i][j];
            }
            if(znScoreFromHomo != NULL)
            {
                znScoreFromHomo[j] = score;
            }

            if (pChain->aaSeq[j] == 'D' || pChain->aaSeq[j] == 'E')/*2007-08-21, for DE, homology based prediction is not good, set the proportion to a lower value*/
            {
                score = (score * 0.1 * weight+ (1.00-0.1) * pChain->discriminant[j]);
            }
            else
            {
                score = (score * proportion * weight+ (1.00-proportion) * pChain->discriminant[j]);
            }

            pChain->discriminant[j] = score;

//             if(score >= 0.35)
//                 pChain->label[j] = 1;
//             else
//                 pChain->label[j] = -1;
        }
    }
	else
	{
		if(numHomoPro >= 0)
		{
			fprintf(stdout,"%s numHomoPro = 0\n", pChain->id);//debug
		}
        else
        {
            fprintf(stdout,"%s no HOMO file\n", pChain->id);
        }

		//if there's no homology detected, applying post scan according to the
		//discriminant
		// ******* to be implemented
        //

        for(j = 0 ; j < pChain->numRes; j ++)
        {
            double score = 0.0; /*2007-07-17, if no homo-file exist, set the score from homopro as 0.0*/
            if(znScoreFromHomo != NULL)
            {
                znScoreFromHomo[j] = score;
            }
            score = (score * proportion * weight+ (1.00-proportion) * pChain->discriminant[j]);
            pChain->discriminant[j] = score;
        }
	}

    int cntZnRes = 0 ;

    for(i = 0 ; i < pChain->numRes; i++)
    {
        if(pChain->discriminant[i] >= 0.0) 
        {
            pChain->label[i] = 1;
            cntZnRes ++;
        }
        else
        {
            pChain->label[i] = -1;
        }
    }

    if(pChain->numRes <= 1)/*this solve 1NJQA*/
    {
        for(i = 0 ; i < pChain->numRes; i++)
        {
            if(pChain->discriminant[i] <= 0.6)
            {
                pChain->label[i] = -1;
                pChain->discriminant[i] -= 0.4;
                cntZnRes = 0;
            }
        }
    }
    /*if there is only one residue highly predicted, while others are all low, decrease the score for that residue*/
    if(pChain->numRes > 1)
    {
        Array1D <int> idxDiscriminant_1darray(pChain->numRes);
        int *idxDiscriminant = idxDiscriminant_1darray.array1D;
        for(i = 0 ; i < pChain->numRes; i ++)
        {
            idxDiscriminant[i] = i;
        }
        QuickSort_index(idxDiscriminant, pChain->discriminant, 0, pChain->numRes-1, DESCENDING);
        if(pChain->discriminant[idxDiscriminant[0]] - pChain->discriminant[idxDiscriminant[1]] > 0.4 && pChain->discriminant[idxDiscriminant[1]] < 0.4)
        {
            pChain->discriminant[idxDiscriminant[0]] = pChain->discriminant[idxDiscriminant[1]];
        }
    }

	NormalizeZnPredScore(pChain->discriminant, pChain->numRes);
    
    numZnRes = cntZnRes;
    return numZnRes;
}
/*}}}*/
void ReportZnPred(GistPredChain *pChain, FILE *fpout, double *exScore1 = NULL, double *exScore2 = NULL)/*{{{*/
{
    int i;
    char vectorRecordID[SIZE_VECTOR_RECORD+1] = "";
    char rmtID[SIZE_CHAIN_ID+1] = "";
    my_strcpy(rmtID, pChain->id, SIZE_CHAIN_ID);
    RemoveTID(rmtID);
    for(i = 0 ; i < pChain->numRes; i ++)
    {
        WriteVectorRecordID(vectorRecordID, rmtID, pChain->length, i, &(pChain->aaSeq[i]), &(pChain->aaSeqIndex[i]), 1);
        fprintf(fpout,"%s\t%2d\t%6.4lf",vectorRecordID, pChain->label[i], pChain->discriminant[i]);
        if(exScore1 != NULL)
        {
            fprintf(fpout, "\t%6.4lf", exScore1[i]);
        }
        if(exScore2 != NULL)
        {
            fprintf(fpout, "\t%6.4lf", exScore2[i]);
        }
        fprintf(fpout,"\n");
    }
}
/*}}}*/

int main(int argc, char** argv)/*{{{*/
{


    if(argc < 2) 
    {
        PrintHelp();
        return 0;
    }

    char outfile[MAX_PATH+1] = "";
    char znpredfile[MAX_PATH+1] = "";
    int level = 1;
    double cutoff_discriminant_1 = 0.70;
    double cutoff_discriminant_2 = 0.60;
    double cutoff_avg_score = 0.45;
    int    min_local_numHRes = 2;
    int    cutoff_window = 50;
    bool isNonOptionArg = false;
    int numSite = 1;
    int operation = USE_AVERAGE;
    int i;

    bool isGatingScore = false;
    double a = 3.0;
    double b = 0.5;

    char modmpath[MAX_PATH+1] = "/misc/casiodata3/wk/passe/modm-mtx";
    char homopath[MAX_PATH+1] = "/misc/casiodata3/wk/passe/Result_nanjiang";
    char metalProFile[MAX_PATH+1] = "";
    char ssbondProFile[MAX_PATH+1] = "";
    char resList[NUM_BLOSUM+1] = "CHDE";
    double cutoff_discriminant = 0.0;
    char **keyMetalList = NULL;
    int numKeyMetal = 0;
    
    int rule = RULE_HOMO ;

    i = 1;
    //int cntID = 0;
    while(i < argc)/*{{{*/
    {
        if(argv[i][0] == '-' && !isNonOptionArg) //options
        {
            isNonOptionArg = false;
            if(strcmp(argv[i],"-h") == 0 ||strcmp(argv[i],"--help")==0 )
            {
                PrintHelp(); 
                return 0;
            }
            if(strcmp(argv[i],"-H") == 0 )
            {
                //fprintf(stdout,"%s\n", HELPMSG);
                PrintVerboseHelp();
                return 0;
            }
            else if(strcmp(argv[i],"-n") == 0)  
            {
                numSite = atoi(argv[i+1]);
                i += 2;
            }
            else if(strcmp(argv[i],"--gate") == 0)  
            {
                isGatingScore = true;
                i ++;
            }
            else if(strcmp(argv[i],"-s1") == 0)  
            {
                cutoff_discriminant_1 = atof(argv[i+1]);
                i += 2;
            }
            else if (strcmp(argv[i],"-s2") == 0)  
            {
                cutoff_discriminant_2 = atof(argv[i+1]);
                i += 2;
            }
            else if (strcmp(argv[i],"-num") == 0)  
            {
                min_local_numHRes = atoi(argv[i+1]);
                i += 2;
            }
            else if(strcmp(argv[i],"-d") == 0)  
            {
                cutoff_window = atoi(argv[i+1]);
                i += 2;
            }
            else if(strcmp(argv[i],"-w") == 0 || strcmp(argv[i], "--weight") == 0)  
            {
                weight = atof(argv[i+1]);
                i += 2;
            }
            else if(strcmp(argv[i],"--c-homo") == 0)  
            {
                cutoff_homoScore = atoi(argv[i+1]);
                i += 2;
            }
            else if(strcmp(argv[i],"--max-homopro") == 0)  
            {
                max_homo_pro = atoi(argv[i+1]);
                i += 2;
            }
            else if(strcmp(argv[i],"--not-use-seg") == 0)  
            {
                isNotUsingSegInfo = false;
                i += 1;
            }
            else if(strcmp(argv[i],"--homo-format") == 0)  
            {
                homofilename_format = atoi(argv[i+1]);
                i += 2;
            }
            else if(strcmp(argv[i],"--prop") == 0 )  
            {
                proportion = atof(argv[i+1]);
                i += 2;
            }
            else if(strcmp(argv[i],"-sa") == 0)  
            {
                cutoff_avg_score = atof(argv[i+1]);
                i += 2;
            }
            else if(strcmp(argv[i],"-o") == 0)  
            {
                my_strcpy(outfile, argv[i+1], MAX_PATH);
                i += 2;
            }
            else if(strcmp(argv[i],"--homo") == 0)  
            {
                my_strcpy(homopath, argv[i+1], MAX_PATH);
                i += 2;
            }
            else if(strcmp(argv[i],"-a") == 0)  
            {
                my_strcpy(resList, argv[i+1], NUM_BLOSUM);
                i += 2;
            }
            else if(strcmp(argv[i],"--level") == 0)  
            {
                level = atoi(argv[i+1]);
                i += 2;
            }
            else if(strcmp(argv[i],"--modm") == 0)  
            {
                my_strcpy(modmpath, argv[i+1], MAX_PATH);
                i += 2;
            }
            else if(strcmp(argv[i],"--metal") == 0)  
            {
                my_strcpy(metalProFile, argv[i+1], MAX_PATH);
                i += 2;
            }
            else if(strcmp(argv[i],"--ssbond") == 0)  
            {
                my_strcpy(ssbondProFile, argv[i+1], MAX_PATH);
                i += 2;
            }
            else if(strcmp(argv[i],"--rule") == 0)  
            {
                if(strncasecmp(argv[i+1], "homo",1) == 0)
                    rule = RULE_HOMO;
                else if(strncasecmp(argv[i+1], "score",1) == 0)
                    rule = RULE_SCORE;

                i += 2;
            }
            else if (strcmp(argv[i], "--") == 0)//next item is non option argument
            {
                isNonOptionArg = true;
                i ++;
                continue;
            }
            else
            {
                fprintf(stderr,"Error! Wrong argument '%s'\n", argv[i]);
                return -1;
            }
        }
        else //non-option argument
        {
            my_strcpy(znpredfile, argv[i], MAX_PATH);
            i ++;
        }
    }
    /*}}}*/
    FILE *fpGistPred = NULL;
    fpGistPred = fopen(znpredfile, "r");
    checkfilestream(fpGistPred, znpredfile, "r");

    FILE *fpout = NULL;
    if(strcmp(outfile, "") == 0)
        fpout = stdout;
    else
    {
        fpout = fopen(outfile,"w");
        checkfilestream(fpout, outfile,"w");
    }

    GistPredChain chain;
    InitGistPredChain(&chain);
    AllocGistPredChain(&chain, LONGEST_SEQ);
    int numRes = 0;
    int numZnRes = 0;
    if(rule == RULE_SCORE)
    {
        while((numRes = GetGistPredResidue(fpGistPred, numSite, &chain, operation))!= EOF)
        {
            numZnRes = ScanZnPred(&chain, cutoff_discriminant_1,cutoff_discriminant_2, cutoff_window, min_local_numHRes, cutoff_avg_score);
            ReportZnPred(&chain, fpout);
        }
    }
    else if(rule == RULE_HOMO)
    {
        // get the metal binding residues for each metal binding protein, metal  /*{{{*/
        // atoms are limited in keyMetalList, if numKeyMetal == 0, all metals are
        // used
        //*   metalPro   -- store metal binding residues bound by metal atoms which bind
        //*                 to >= min_metalBoundRes and <= max_metalBoundRes residues, within resList
        //*                 and filtered by cutoff_score2
        //*   metalPro1  -- same as metalPro, but the cutoff_score2 restriction is removed, the resList rule keeps
        //*   metalPro2  -- same as metalPro, but both resList and cutoff_score2 restriction are removed, that is 
        //*                 store all metal binding residues bound by metal atoms which
        //*                 bind to >= min_metalBoundRes and <= max_metalBoundRes (in
        //*                 the whole protein) residues. isExcludeOtherChainRes = true,
        //*                 that is only residues on the nrPDB chain is included, 
        using namespace MetalBindingProtein;
        Array1D <MetalPro> metalPro_1darray(MAX_NUM_METALPRO);
        Array1D <MetalPro> metalPro1_1darray(MAX_NUM_METALPRO);
        Array1D <MetalPro> metalPro2_1darray(MAX_NUM_METALPRO);
        MetalPro *metalPro  = metalPro_1darray.array1D;
        MetalPro *metalPro1 = metalPro1_1darray.array1D;
        MetalPro *metalPro2 = metalPro2_1darray.array1D;
        for(i = 0; i < MAX_NUM_METALPRO ; i++)
        {
            InitMetalPro(&metalPro[i]);
            InitMetalPro(&metalPro1[i]);
            InitMetalPro(&metalPro2[i]);
        }
        int numMetalPro;
        int numMetalPro1;
        int numMetalPro2;
        int numMetalRes;
        int numMetalRes1;
        int numMetalRes2;
        if(level == 0)
        {
            GetMetalBoundRes3(metalProFile, metalPro, metalPro1, metalPro2, numMetalPro, numMetalPro1, numMetalPro2, numMetalRes, numMetalRes1, numMetalRes2, cutoff_score2, resList, modmpath, isExcludeOtherChainRes, keyMetalList, numKeyMetal, min_metalBoundRes, max_metalBoundRes, isUsingTotalBoundRes);
        }
        else
        {
            GetMetalBoundRes2(metalProFile, metalPro1, metalPro2, numMetalPro1, numMetalPro2, numMetalRes1, numMetalRes2, resList, modmpath, isExcludeOtherChainRes, keyMetalList, numKeyMetal, min_metalBoundRes, max_metalBoundRes, isUsingTotalBoundRes);
        }
        // *}}}*/
        // Get MetalProLevel according to the level/*{{{*/
        MetalPro *metalProLevel = NULL;
        int numMetalProLevel = 0;
        if(level == 0)
        {
            metalProLevel = metalPro;
            numMetalProLevel = numMetalPro;
        }
        else if(level == 1)
        {
            metalProLevel = metalPro1;
            numMetalProLevel = numMetalPro1;
        }
        else if(level == 2)
        {
            metalProLevel = metalPro2;
            numMetalProLevel = numMetalPro2;
        }/*}}}*/

        //get SSBond /*{{{*/
        //ssbondPro1:  is including all SSBonded residues on the chain
        //ssbondPro :  is with those cysteins with score2 < 0.1 removed
        Array1D <SSBondPro> ssbondPro_1darray(MAX_SSBOND_PRO);
        Array1D <SSBondPro> ssbondPro1_1darray(MAX_SSBOND_PRO);
        SSBondPro *ssbondPro  = ssbondPro_1darray.array1D;
        SSBondPro *ssbondPro1 = ssbondPro1_1darray.array1D;
        for(i = 0 ; i < MAX_SSBOND_PRO; i++)
        {
            InitSSBondPro(&ssbondPro[i]);
            InitSSBondPro(&ssbondPro1[i]);
        }
        int numSSBondPro;
        int numSSBondPro1;
        int numSSBondRes;
        int numSSBondRes1;
        if(level == 0)
        {
            GetSSBondRes2(ssbondProFile, ssbondPro,ssbondPro1,numSSBondPro, numSSBondPro1, numSSBondRes, numSSBondRes1, cutoff_score2, modmpath);
        }
        else
        {
            GetSSBondRes(ssbondProFile, ssbondPro1, numSSBondPro1, numSSBondRes1);
        }
        /*}}}*/
        // Get ssbondProLevel according to the level/*{{{*/
        SSBondPro *ssbondProLevel = NULL;
        int numSSBondProLevel = 0;
        if(level == 0)
        {
            ssbondProLevel = ssbondPro;
            numSSBondProLevel = numSSBondPro;
        }
        else if(level == 1)
        {
            ssbondProLevel = ssbondPro1;
            numSSBondProLevel = numSSBondPro1;
        }
        /*}}}*/


        Array2D <char> metalProIDList_2darray(numMetalProLevel, SIZE_CHAIN_ID+1);
        Array2D <char> ssbondProIDList_2darray(numSSBondProLevel, SIZE_CHAIN_ID+1);
        char **metalProIDList = metalProIDList_2darray.array2D;
        char **ssbondProIDList = ssbondProIDList_2darray.array2D;
        for(i = 0 ; i < numMetalProLevel; i ++) { my_strcpy(metalProIDList[i], metalProLevel[i].id , SIZE_CHAIN_ID); }
        for(i = 0 ; i < numSSBondProLevel; i ++) { my_strcpy(ssbondProIDList[i], ssbondProLevel[i].id , SIZE_CHAIN_ID); }

        while((numRes = GetGistPredResidue(fpGistPred, numSite, &chain, operation))!= EOF)
        {
            Array1D <double> antiSigScore_1darray(numRes);
            double *antiSigScore = antiSigScore_1darray.array1D;
            Array1D <double> znScoreFromHomo_1darray(numRes);
            double *znScoreFromHomo = znScoreFromHomo_1darray.array1D;
            for(i = 0 ; i  < numRes; i ++)
            {
                antiSigScore[i] = INIT_DOUBLE;
                znScoreFromHomo[i] = INIT_DOUBLE;
            }
            if(isGatingScore)
            {
                for(i = 0 ; i < chain.numRes ; i ++)
                {
                    chain.discriminant[i] = AntiSigmoidScore(chain.discriminant[i], a , b);
                    antiSigScore[i] = chain.discriminant[i];
                }
            }
            numZnRes = ScanZnPred_HOMO(&chain, metalProLevel, metalProIDList, numMetalProLevel, ssbondProLevel, ssbondProIDList, numSSBondProLevel, homopath, cutoff_discriminant, znScoreFromHomo, fpout);
#ifdef DEBUG
            for(i = 0 ; i < numRes ; i ++)
            {
                if(antiSigScore[i] == INIT_DOUBLE)
                {
                    fprintf(stderr,"antiSigScore[%d] = %lf\n", i, antiSigScore[i]);
                }
                if(znScoreFromHomo[i] == INIT_DOUBLE)
                {
                    fprintf(stderr,"znScoreFromHomo[%d] = %lf\n",i, znScoreFromHomo[i]);
                }
            }
#endif
            ReportZnPred(&chain, fpout, antiSigScore, znScoreFromHomo);
        }

        for(i = 0; i < MAX_NUM_METALPRO ; i++)
        {
            DeleteMetalPro(&metalPro[i]);
            DeleteMetalPro(&metalPro1[i]);
            DeleteMetalPro(&metalPro2[i]);
        }
        for(i = 0 ; i < MAX_SSBOND_PRO; i++)
        {
            DeleteSSBondPro(&ssbondPro[i]);
            DeleteSSBondPro(&ssbondPro1[i]);
        }

    }


    fclose(fpGistPred);
    if(fpout != stdout && fpout != NULL) fclose(fpout);
    

    DeleteGistPredChain(&chain);

    return EXIT_SUCCESS;
}
/*}}}*/

