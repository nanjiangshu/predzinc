/*
 * =====================================================================================
 * 
 *        Filename:  create_svm_vector.cpp
 *     Description:  Create feature vectors for gist SVM 
 *         Version:  1.0
 *         Created:  10/31/2006 06:08:07 PM CEST
 *        Revision:  none
 *        Compiler:  g++
 *          Author:  Nanjiang Shu (Shu), nanjiang@struc.su.se
 *         Company:  Structural Chemistry, Stockholm Univesity
 * 
 * =====================================================================================
 */
#undef DATATYPE_SMATRIX
#define DATATYPE_SMATRIX double

#define _ASSERT_
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include "array.h"
#include "mytemplate.h"
#include "myfunc.h"
#include "mypro.h"
#include "subset.h"

#undef SIZE_VECTOR_RECORD
#define SIZE_VECTOR_RECORD 500

#undef PRED
#undef TRAIN
#define PRED 0
#define TRAIN 1

#undef METAL
#undef SS
#define METAL 0
#define SS 1

/*ChangeLog {{{
 * ChangeLog 2011-10-11
 *     modmfilelistfile also supported in predict mode
 * }}}*/


//bool isWriteVectorVectorFile = true;
//bool isWriteVectorLableFile = false;
bool isMaskPolyHis = false;
bool isUseConsvI = false;
bool isAddPairDist = true;
bool isUseSHCRRule = true;
char metalElementListFile [MAX_PATH+1] = "";
char metalTranMatrixDir[MAX_PATH+1] = "";


int profile_encoding_type = PSSM_PROFILE_ENCODE;

void PrintHelp()
{
    fprintf(stdout,"Usage: create_svm_vector modmfile OR -l modmfilelistfile\n");
    fprintf(stdout,"options\n");
    //fprintf(stdout," -l filelist  : batch mode, create svm vectors for a list of sequences\n");
    fprintf(stdout," --mode pred|train   : create vectors for prediction or build the training set\n");
    fprintf(stdout,"                     : if the mode is train, metalProFile should be set for metal and \n");
    fprintf(stdout,"                     : ssbondProFile should be set for SS, default = pred\n");
    fprintf(stdout," -t Metal|SS         : create vector for metal-binding or ssbond residues,  default=Metal\n");
    fprintf(stdout,"  --keymetal  strs   : list keyMetals, if keyMetal does not set, all metals are included\n");
    fprintf(stdout," -l modmfilelistfile     : filelist for modmfiles\n");
    fprintf(stdout," -ns value           : numSite, 1 for single site vector and 2 for pair-based vector, default=1\n");
    fprintf(stdout," --profile pssm|seq  : encode for profile, using pssm profle or sequences binary encoding, default = pssm\n");
    fprintf(stdout," --vector vectorFile : set name for vector file, default is rootname(seq-file).vector\n");
    fprintf(stdout," --label  labelFile  : output label file also, in that case, metal-binding file should be supplied\n");
    fprintf(stdout," --modm-type LOG|PER : modm file type, log odd or percentage, default=log\n");
    fprintf(stdout," -s value            : cutoff_score2, default = 0.01\n");
    fprintf(stdout," -c value            : cutoff_consv,  default = 0.70\n");
    fprintf(stdout," -a resList          : residues type to be selected, default='CHDE' for zinc-binding and 'C' for ssbond\n");
    fprintf(stdout," -K                  : parameter for vector , default = 10 for zinc-binding and 10 for ssbond\n");
    fprintf(stdout," -W                  : parameter for vector , default = 5\n");
    fprintf(stdout," -P                  : parameter for vector , default = 24\n");
    fprintf(stdout," --encode 0|1        : parameter for encoding scheme, 0: old scheme, 1: new scheme, default = 1,\n");
    fprintf(stdout," --min-hcres value   : minimal number of HCRes required,  default = 3\n");
    fprintf(stdout," --win-pair value    : window for residue pair   , default = 7\n");
    fprintf(stdout," --win-min-pair value: window for the minimal residue pair , default = 5\n");
    fprintf(stdout," --win3 value        : window for 3-residue group, default = 150\n");
    fprintf(stdout," --win4 vlaue        : window for 4-residue group, default = 150\n");
    fprintf(stdout," --ts3 vlaue         : minimal transtion score for Metal3, used when numSite >=2, default = 0.65\n");
    fprintf(stdout," --ts4 vlaue         : minimal transtion score for Metal4, used when numSite >=2, default = 1.0\n");
    fprintf(stdout," --bound-res i j     : number of binding residues to Metal, from i to j,  default = 3 4\n");
    fprintf(stdout," --level 0|1|2       : statistical level, default = 1, for fast loading metal pro\n");

    fprintf(stdout," --metal metalfile   : close metal list file, default = $PREDZINC/data/closeMetalPro.dat\n");
    fprintf(stdout," --metaltrans DIR    : Path for transmatrix\n");
    fprintf(stdout," --ssbond ssbondFile : ssbond list file, default = $PREDZINC/data/ssbond.dat\n");
    fprintf(stdout," --modmpath path     : path for modm files,   default = $DATADIR/modm\n");
    //fprintf(stdout," --outpath path    : path for output files, default = ./\n");
    fprintf(stdout," --not-add-pair-dist : add distance separating pair in the vector, default is using\n");
    fprintf(stdout," --not-mask-polyhis  : do not mask poly histag , default is to mask\n");
    fprintf(stdout," --not-use-ic        : do not use integrated conservation level , default is to use \n");
    fprintf(stdout," --not-use-shcr      : do not use shcr selection rule by limiting the window size, default is to use \n");
    fprintf(stdout," --use-sc            : use only bound residues from the same chain, default is to use residues from all chains\n");
    fprintf(stdout," --neg-filter int    : filter negtive samples in case of the file is too big, default = 1, not to filter\n");
    fprintf(stdout," -h | --help         : print this help message and exit\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"When create training database, the chain id should be in the format according to the readme file\n");
    fprintf(stdout,"Created 2006-10-31, update 2011-10-11, Nanjiang Shu\n");
}

double NormalizeConsv(double consv, char aa)/*{{{*/
{
    double norm_consv = 0.0;
    if (aa == 'C')
    {
        norm_consv = (consv - (-5.0)) / ( 12.0 - (-5.0));
    }
    else if (aa == 'H')
    {
        norm_consv = (consv - (-5.0)) / ( 11.0 - (-5.0));
    }
    else if (aa == 'D')
    {
        norm_consv = (consv - (-5.0)) / ( 8.0 - (-5.0));
    }
    else if (aa == 'E')
    {
        norm_consv = (consv - (-5.0)) / ( 8.0 - (-5.0));
    }

    if (norm_consv > 1.0)
    {
        norm_consv = 1.0;
    }
    else if (norm_consv < 0.0)
    {
        norm_consv = 0.0;
    }
    return norm_consv;

}/*}}}*/
int GetVectorDimension(int K, int W, int numSite, int P, int encoding_scheme, bool isAddPairDist)/*{{{*/
/*****************************************************************************
 * 2007-06-28
 * calculate the size of feature vector for each sample
 ****************************************************************************/
{
    int vectorDim  = 0;
    if(encoding_scheme == 0) // added 2007-06-28
    {
        vectorDim = (2*K + numSite + numSite-1)*P;
        if(isAddPairDist && numSite > 1)
        {
            vectorDim += numSite -1 ;
        }
    }
    else // using GetSVMVector_2 to encoding the feature vector
    {
        int dimVectorPerSite_Center = GetDimensionVectorPerSite(P, true);
        int dimVectorPerSite_Non_Center = GetDimensionVectorPerSite(P, false);
        vectorDim = (2*K +  2*W* (numSite-1) )*dimVectorPerSite_Non_Center+ numSite * dimVectorPerSite_Center;
        if(isAddPairDist && numSite > 1)
        {
            vectorDim += SIZE_ENCODE_DISTANCE * (numSite-1);
        }

    }
    return vectorDim;
}/*}}}*/
int GetMetalResPair(Residue *HCRes, int numHCRes, int numSite, int sizeGroup, int cutoff_window_pair, int cutoff_window, int cutoff_min_window_pair, double cutoff_transcore, double ***tranM, char *alphabetTranM, set <string> &resPair)/*{{{*/
{
    using namespace MetalBindingProtein;
    int i,j;

    Array1D <int> aaSeqIndex_1darray(numSite+1);
    Array1D <char> res_1char_list_1darray(numSite+2);
    int  *aaSeqIndex = aaSeqIndex_1darray.array1D;
    char *res_1char_list = res_1char_list_1darray.array1D;


    Array1D <int> a_1darray(MAX_BOUND_SITE);
    Array1D <int> idx_1darray(MAX_BOUND_SITE);
    int *a = a_1darray.array1D;
    int *idx = idx_1darray.array1D;

    int nSub;
    int nSet;
    //int nSub2;
    //int nSet2;

    int countSubset; /* number of subset, which is the combination number of C(nSet, nSub)*/ 
    countSubset = combin2(MAX_BOUND_SITE,2);
    Array2D <int> subset_2darray(countSubset,3);
    int **subset = subset_2darray.array2D;
    //double score;

    Array1D <int> win_pair_1darray(MAX_BOUND_SITE+1);
    int *win_pair = win_pair_1darray.array1D;

    Array1D <Residue> tmpResidues_1darray(numHCRes);
    Residue *tmpResidues = tmpResidues_1darray.array1D;
    for(i = 0 ; i < numHCRes; i++) InitResidue(&tmpResidues[i]);

    char str1[SIZE_VECTOR_RECORD+1] = "";
    char str2[SIZE_VECTOR_RECORD+1] = "";

    if(numHCRes == 2)
    {
        for(i = 0 ; i < numSite ; i++)
        {
            aaSeqIndex[i] = HCRes[i].aaSeqIndex;
            res_1char_list[i] = HCRes[i].aa;
        }
        res_1char_list[numSite] = '\0';

        sprintf(str1,"%c_%d",res_1char_list[0],aaSeqIndex[0]+1);
        for(i = 1 ; i < numSite ; i++)
        {
            sprintf(str2,"_%c_%d", res_1char_list[i], aaSeqIndex[i]+1);
            strcat(str1,str2);
        }
        resPair.insert(str1);
    }

    if(numHCRes >= sizeGroup)
    {
        nSub = sizeGroup;
        nSet = numHCRes;
        for(i = 0 ; i < nSub ; i++) a[i] = 0;
        bool done = true;
        for( ; ; )
        {
            ksub_next4(nSet,nSub,a,&done);
            if(done) break;

            for(i = 0; i < nSub; i ++) idx[i] = a[i]-1;

            int win = HCRes[idx[nSub-1]].aaSeqIndex - HCRes[idx[0]].aaSeqIndex; 

            for(i = 0 ; i < nSub - 1; i ++) 
                win_pair[i] = HCRes[idx[i+1]].aaSeqIndex - HCRes[idx[i]].aaSeqIndex;
            int min_win_pair =  *min_element(win_pair, win_pair+nSub-1) ;

            if(win <= cutoff_window && min_win_pair <= cutoff_min_window_pair)
            {
                for(i = 0 ; i < nSub ; i ++)
                {
                    CopyResidue(&tmpResidues[i],&HCRes[idx[i]]);
                }
                //score = TranScore(tranM, tmpResidues,alphabetTranM,nSub);

#ifdef DEBUG
                fprintf(stdout, "trans-score %6.3lf\n", score);
#endif

                //if(score >= cutoff_transcore)
                {
                    int nSet2 = nSub;
                    int nSub2  = numSite;
                    int numSubset = GetAllSubset(subset, nSet2, nSub2);

                    for(i = 0 ; i < numSubset ;i ++)
                    {
                        for(j = 0 ; j < numSite ; j++)
                        {
                            aaSeqIndex[j]     = tmpResidues[subset[i][j]].aaSeqIndex;
                            res_1char_list[j] = tmpResidues[subset[i][j]].aa;
                        }
                        res_1char_list[numSite] = '\0';


                        // debug start
                        //if(aaSeqIndex[0] > pMODM->length -1 || aaSeqIndex[nSub2-1] > pMODM->length -1)
                        //{
                            //printf("aaSeqIndex1 = %d, aaSeqIndex2 = %d, length = %d", aaSeqIndex[0], aaSeqIndex[nSub2-1], pMODM->length);
                            //assert( aaSeqIndex[0] < pMODM->length && aaSeqIndex[nSub2-1] < pMODM->length);
                        //}
                        // debug end
                        if(aaSeqIndex[nSub2-1] - aaSeqIndex[0] <= cutoff_window_pair)
                        {
                            sprintf(str1,"%c_%d",res_1char_list[0],aaSeqIndex[0]+1);
                            for(j = 1 ; j < numSite ; j++)
                            {
                                sprintf(str2,"_%c_%d", res_1char_list[j], aaSeqIndex[j]+1);
                                strcat(str1,str2);
                            }
                            resPair.insert(str1);
                            // printf("resPair=%s\n",str);// debug
                        }
                    }
                }
            }
        }
    }

    return resPair.size();
}
/*}}}*/
int GetMetalResPairAll(Residue *HCRes, int numHCRes, int cutoff_window_pair, int cutoff_window3, int cutoff_window4, int cutoff_min_window_pair, double cutoff_transcore_Metal3, double cutoff_transcore_Metal4, int **resGroupIndex, char **resGroupAA)/*{{{*/
/*****************************************************************************
 * Given zinc-binding residue candidates selected based on highly conserved
 * CHDEs, 
 * return the number of residue groups, residue groups are in two dimensional array
 ****************************************************************************/
{
    int i, j ;

    using namespace MetalBindingProtein;
    // get transmission matrix/*{{{*/
    char datadir[MAX_PATH+1]  ="";
    GetDataDir(datadir);
    Array3D <DATATYPE_SMATRIX>  tranM_Metal3_3darray(2, MAX_BOUND_SITE,MAX_BOUND_SITE);
    Array3D <DATATYPE_SMATRIX>  tranM_Metal4_3darray(3,MAX_BOUND_SITE,MAX_BOUND_SITE);
    DATATYPE_SMATRIX ***tranM_Metal3 = tranM_Metal3_3darray.array3D;
    DATATYPE_SMATRIX ***tranM_Metal4 = tranM_Metal4_3darray.array3D;
    char alphabetTranM[NUM_BLOSUM+1];
    char tranMatrixFile[MAX_PATH+1];
    for( i = 3 ; i <= 4 ; i ++)  // Metal3 and Metal4
    {
        for(j = 0 ; j < i-1; j++)
        {
            sprintf(tranMatrixFile,"%s/%s%d%d.mtx",metalTranMatrixDir,"tran",i,j);
            if (i == 3) GetPatternTransMatrix(tranMatrixFile,tranM_Metal3[j], alphabetTranM);
            if (i == 4) GetPatternTransMatrix(tranMatrixFile,tranM_Metal4[j], alphabetTranM);
        }
    }
    /*}}}*/

    int numSite = 2;
    set <string> resPair;
    resPair.clear();
    set <string>::iterator iss;  /* iterator for set <string>*/ 
    char str[SIZE_VECTOR_RECORD+1] = "";

    GetMetalResPair(HCRes, numHCRes, numSite, 3, cutoff_window_pair, cutoff_window3, cutoff_min_window_pair, cutoff_transcore_Metal3,  tranM_Metal3, alphabetTranM, resPair);
    GetMetalResPair(HCRes, numHCRes, numSite, 4, cutoff_window_pair, cutoff_window4, cutoff_min_window_pair, cutoff_transcore_Metal3,  tranM_Metal4, alphabetTranM, resPair);

    i = 0;
    for( iss = resPair.begin() ; iss != resPair.end(); iss++ )
    {
        my_strcpy(str, (*iss).c_str(), SIZE_VECTOR_RECORD);
        char *pch;
        char delim[] = "_";
        pch = strtok (str,delim);
        j = 0 ; 
        while (pch != NULL)
        {
            if(j%2 == 0)
                sscanf(pch,"%c", &resGroupAA[i][j/2]);
            else
            {
                sscanf(pch,"%d",&resGroupIndex[i][j/2]);
                resGroupIndex[i][j/2] --;
            }
            pch = strtok (NULL, delim);
            j ++;
        }
        i ++;
    }
    return resPair.size();
}
/*}}}*/
void WriteVectorHeader(int vectorDim, FILE *fpVector)/*{{{*/
{
    fprintf(fpVector,"%s", "item");
    int j;
    for(j = 0; j < vectorDim; j++) 
    {
        fprintf(fpVector, "\tv%d", j);
    }
    fprintf(fpVector,"\n"); 
}
/*}}}*/
void WriteLabelHeader(FILE *fpLabel)/*{{{*/
{
    fprintf(fpLabel,"%s", "item");
    fprintf(fpLabel,"\t%s", "label");
    fprintf(fpLabel,"\n");
    //int j;
}
/*}}}*/

int GetIsBoundArray(char **metalProIDList, int *idxSortMetalPro, MetalPro *metalPro, int numMetalPro, MODM *pMODM, bool *isBoundArray)/*{{{*/
{
    int indexMetalPro;
    int i;
    for(i = 0 ; i < pMODM->length; i++) isBoundArray[i] = false;
    if( (indexMetalPro = BinarySearch_String(pMODM->id,metalProIDList, numMetalPro)) != -1 )   //make sure metalProIDList is sorted ascendingly, 2007-06-21
    {
        indexMetalPro = idxSortMetalPro[indexMetalPro]; //map the index to the unsorted array
        for(i = 0 ; i < metalPro[indexMetalPro].numBoundRes; i++)
        {
            isBoundArray[metalPro[indexMetalPro].resSeqIndex[i]] = true;
        }
        return  metalPro[indexMetalPro].numBoundRes;
    }
    else 
        return 0;
}
/*}}}*/
int GetIsBoundArray(char **ssbondProIDList, int *idxSortSSBondPro, SSBondPro *ssbondPro, int numSSBondPro, MODM *pMODM, bool *isBoundArray)/*{{{*/
{
    int indexSSBondPro;
    int i, j;
    for(i = 0 ; i < pMODM->length; i++) isBoundArray[i] = false;
    if( (indexSSBondPro = BinarySearch_String(pMODM->id,ssbondProIDList, numSSBondPro)) != -1 ) // make sure the ssbondProIDList is sorted
    {
        indexSSBondPro = idxSortSSBondPro[indexSSBondPro];//map the index to the unsorted array
        for(i = 0 ; i < ssbondPro[indexSSBondPro].numSSBond; i++)
        {
            for(j = 0 ; j < ssbondPro[indexSSBondPro].ssbond[i].numRes; j++)
            isBoundArray[ssbondPro[indexSSBondPro].ssbond[i].res[j].aaSeqIndex] = true;
        }
        return  ssbondPro[indexSSBondPro].numSSBondRes;
    }
    else 
        return 0;
}
/*}}}*/

int GetBondLabel(int *resIndex, int numSite, bool *isBoundArray)/*{{{*/
/*****************************************************************************
 * BondLable = 1 only when all residues in the resIndex are bounded
 ****************************************************************************/
{
    int i;
    for(i = 0 ; i < numSite ; i++)
    {
        if(isBoundArray[resIndex[i]] == false)
            return -1;
    }
    return 1;
}
/*}}}*/

int WriteVector(int **resGroupIndex, char **resGroupAA, int numResGroup, int K, int W, int P, MODM *pMODM, int numSite, int vectorDim, bool *isBoundArray, int neg_filter, FILE *fpVector, FILE *fpLabel, double cutoff_score2, int encoding_scheme)/*{{{*/
/*****************************************************************************
 *  changelog 2007-06-28
 *      adding encoding_scheme, 
 *      encoding_scheme == 0, using the old GetSVMVector
 *      encoding_scheme == 1, using GetSVMVector_2
 ****************************************************************************/
{
    int i, j;
    int cntVector = 0;
    int cnt_neg = 0;
    int label;
    char vectorRecordID[SIZE_VECTOR_RECORD+1] = "";
    char rmtID[SIZE_CHAIN_ID+1] = "";
    SpanExcluding(pMODM->id, rmtID);


    Array1D <double> vector_1darray(vectorDim);
    double *vector = vector_1darray.array1D;
    for(i = 0 ; i < vectorDim; i ++) { vector[i] = INIT_DOUBLE; }

    for(i = 0 ; i < numResGroup ; i ++)
    {
        if( isBoundArray != NULL)
        { 
            label = GetBondLabel(resGroupIndex[i], numSite, isBoundArray); 
        }
        else
        { 
            label = 1; 
        }// originally label == 1, bug fixed, 2007-06-20, if isBoundArray == NULL, means it is for prediction, then set all label to 1 first
            
#ifdef DEBUG
        fprintf(stdout, "lable=%d\n", label);
#endif
        if (label == 1 || ((cnt_neg % neg_filter)== 0))
        {

            int returnSizeVector = 0;
            if(encoding_scheme == 0)
            {

                returnSizeVector = GetSVMVector(resGroupIndex[i], numSite, pMODM, K, W,P, vector, vectorDim, cutoff_score2, profile_encoding_type);  /* bug fixed, 2007-04-13, input the parameter cutoff_score2 as well, since in the function GetSVMVector, the default cutoff_score2 is 0.1 */
            }
            else /*(encoding_scheme == 1)*/
            {
                returnSizeVector = GetSVMVector_2(resGroupIndex[i],numSite, pMODM, K, W, P, vector, vectorDim, cutoff_score2, profile_encoding_type);
            }

            if(returnSizeVector > 0)
            {
#ifdef DEBUG
                for(j = 0; j < vectorDim; j ++)
                {
                    if(vector[j] == INIT_DOUBLE)
                    {
                        fprintf(stdout,"vectorDim=%d, vector[%d]=%lf\n", vectorDim, j, vector[j]);
                    }
                }
#endif
                WriteVectorRecordID(vectorRecordID, rmtID, pMODM->length,cntVector, resGroupAA[i],resGroupIndex[i], numSite);
                fprintf(fpVector,"%s", vectorRecordID);
                for(j = 0; j < vectorDim ; j++) { fprintf(fpVector,"\t%g",vector[j]);}
                fprintf(fpVector,"\n");

                if(fpLabel != NULL)
                {
                    fprintf(fpLabel,"%s\t%+d\n", vectorRecordID, label);
                }
                cntVector ++;
            }
        }
        if(label != 1) 
        {
            cnt_neg ++;
        }
    }
    return cntVector;
}
/*}}}*/

int CreateVector_Metal(int K, int W, int P, int numSite, int vectorDim, double cutoff_score2, double cutoff_consv, int cutoff_window_pair, int cutoff_window3, int cutoff_window4, int cutoff_min_window_pair, double cutoff_transcore_Metal3, double cutoff_transcore_Metal4, int min_numHCRes, char* resList, MODM *pMODM, FILE *fpVector, FILE *fpLabel, int neg_filter, bool *isBoundArray, int encoding_scheme)/*{{{*/
/*****************************************************************************
 * Get feature vectors for all preliminary selected zinc-binding residues and
 * write only the vectors, not the header line.
 * If isBoundArray = NULL, create vectors for testing set,
 * in that case, all labels will be set to 1
 ****************************************************************************/
{
    int i;

    //int  cnt = 0;
    set <string> resPair;
    Array1D <bool>   isPolyHis_1darray(pMODM->length);
    bool * isPolyHis = isPolyHis_1darray.array1D;

    Array1D <int> HCResAASeqIndex_1darray(pMODM->length);
    Array1D <int> shcr_res_idx_1darray(pMODM->length);
    int *HCResAASeqIndex = HCResAASeqIndex_1darray.array1D;
    int *shcr_res_idx = shcr_res_idx_1darray.array1D;
    int numHCRes;
    int numSHCRRes;
    Array1D <Residue> HCRes_1darray(pMODM->length);
    Residue *HCRes =  HCRes_1darray.array1D;

    for ( i = 0 ; i < pMODM->length ; i ++)
    {
        InitResidue(&HCRes[i]);
        isPolyHis[i] = false;
    }

    MaskPolyHis(pMODM->aaSeq, isPolyHis, pMODM->length);

    numHCRes = GetHCRes(pMODM, isPolyHis, HCRes, cutoff_score2, cutoff_consv, resList, isUseConsvI, isMaskPolyHis);

    for(i = 0 ; i < numHCRes; i ++) HCResAASeqIndex[i] = HCRes[i].aaSeqIndex;

    if(isUseSHCRRule) // add 2007-04-24
    {
        if(!IsSatisfySHCRRule(HCResAASeqIndex,numHCRes,min_numHCRes, cutoff_window_pair,cutoff_window3, cutoff_window4,  shcr_res_idx, numSHCRRes )) //bug fixed 2007-04-13, in IsSatisfySHCRRule(), when min_numHCRes is set to < 2
        {
            return 0;
        }
        //CC//: update HCRes by shcr_res_idx
        numHCRes = UpdateHCRes(HCRes, numHCRes, shcr_res_idx, numSHCRRes);
    }

    int  **resGroupIndex = NULL;
    char **resGroupAA = NULL;

    int numResGroup = 0;
    if(numSite == 1)
    {
        resGroupIndex = Create2DArray(resGroupIndex, numHCRes, numSite);
        resGroupAA = Create2DArray(resGroupAA, numHCRes, numSite+1);

        for( i = 0 ; i < numHCRes ; i++)
        {
            resGroupIndex[i][0] = HCRes[i].aaSeqIndex;
            resGroupAA[i][0]    = HCRes[i].aa; resGroupAA[i][1] = '\0';
        }
        numResGroup = numHCRes;
    }
    else if(numSite == 2)
    {
        resGroupIndex = Create2DArray(resGroupIndex, numHCRes*(numHCRes-1)/2+1, numSite);
        resGroupAA = Create2DArray(resGroupAA, numHCRes*(numHCRes-1)/2+1, numSite+1);

        numResGroup = GetMetalResPairAll(HCRes, numHCRes,  cutoff_window_pair, cutoff_window3, cutoff_window4, cutoff_min_window_pair, cutoff_transcore_Metal3, cutoff_transcore_Metal4, resGroupIndex, resGroupAA);
        
    }
    
    int numVector;
    numVector = WriteVector(resGroupIndex, resGroupAA, numResGroup,  K, W, P, pMODM, numSite, vectorDim, isBoundArray, neg_filter, fpVector, fpLabel, cutoff_score2, encoding_scheme);
    // close file streams
    // delete allocated memory
    if(numSite == 1)
    {
        Delete2DArray(resGroupIndex, numHCRes);
        Delete2DArray(resGroupAA, numHCRes);
    }
    else if(numSite == 2)
    {
        Delete2DArray(resGroupIndex, numHCRes*(numHCRes-1)/2+1);
        Delete2DArray(resGroupAA, numHCRes*(numHCRes-1)/2+1);
    }

    return numVector;
}
/*}}}*/

int CreateVectorPredict_Metal(int K, int W, int P, int numSite, double cutoff_score2, double cutoff_consv, int cutoff_window_pair, int cutoff_window3, int cutoff_window4, int cutoff_min_window_pair, double cutoff_transcore_Metal3, double cutoff_transcore_Metal4, int min_numHCRes, char* resList, set<string> modmfileset, int type_modm, const char *vectorFile, int encoding_scheme)/*{{{*/
/*****************************************************************************
 * Create feature vectors for the protein sequence going to predict
 * return the number of vectors
 * updated 2011-10-11, support multiple modmfiles
 ****************************************************************************/
{
    int neg_filter = 1;
    bool *isBoundArray = NULL;
    FILE *fpVector ;
    FILE *fpLabel = NULL;
    fpVector = fopen(vectorFile,"w");
    checkfilestream(fpVector,vectorFile,"w");

    char id[SIZE_CHAIN_ID+1] = "";
    int vectorDim = 0;
    vectorDim = GetVectorDimension(K, W, numSite, P, encoding_scheme, isAddPairDist);

    // each residue pair is represented by a vector of size (2k+2+1)*p+1, where
    // k and p are the same as described in the single site vector, the second
    // 2 is for residues consisting the pair, the first 1 is for the unified
    // residues and the second 1 indicates the number of residues separating
    // the residues of the pair. 
    WriteVectorHeader(vectorDim, fpVector);

    char str1[2*MAX_PATH+1] = "";
    char str2[2*MAX_PATH+1] = "";

    MODM modm;
    MODM *pMODM;
    InitMODM(&modm);
    AllocMODM(&modm, LONGEST_SEQ, NUM_BLOSUM);
    modm.type_modm = type_modm;

    //char command[600+1] = "";
    char modmfilepath[500]="";
    int numVector = 0;
    set <string> ::iterator iss;
    for (iss = modmfileset.begin();iss != modmfileset.end();iss++){
        strcpy(modmfilepath, (*iss).c_str());
        if (strcmp(modmfilepath, "") == 0){
            continue;
        }
        rootname(modmfilepath,str1); // get ID for the record
        SpanExcluding(str1,str2," ");
        my_strcpy(id,str2,SIZE_CHAIN_ID);

        //fprintf(stderr,"modmfilepath=%s\n", modmfilepath); [>DEBUG<]

        modm.length = GetMODM(modmfilepath, modm.M, modm.alphabetMODM, modm.aaSeq, modm.score1, modm.score2);
        if (modm.length <= 0){
            fprintf(stderr,"length of modmfile %s <=0. Ingore. \n", modmfilepath);
            continue;
        }

        int j;
        /*added 2011-10-11, reset bad score2 values to 0.0, e.g. inf*/
        for ( j=0;j<modm.length;j++){
            if (modm.score2[j] < 0.0 || modm.score2[j] > 10.0) {
                modm.score2[j] = 0.0;
            }
        }

        my_strcpy(modm.id,id,SIZE_CHAIN_ID);
        pMODM = &modm;

        numVector += CreateVector_Metal( K,  W,  P,  numSite,  vectorDim,  cutoff_score2,  cutoff_consv,  cutoff_window_pair,  cutoff_window3,  cutoff_window4,  cutoff_min_window_pair,  cutoff_transcore_Metal3,  cutoff_transcore_Metal4,  min_numHCRes, resList, pMODM,  fpVector,  fpLabel,  neg_filter, isBoundArray, encoding_scheme);

    }
    fclose(fpVector);
    DeleteMODM(&modm,LONGEST_SEQ);
    return numVector;
}
/*}}}*/

int CreateVectorTrain_Metal(int K, int W, int P, int numSite, double cutoff_score2, double cutoff_consv, int cutoff_window_pair, int cutoff_window3, int cutoff_window4, int cutoff_min_window_pair, double cutoff_transcore_Metal3, double cutoff_transcore_Metal4, bool isExcludeOtherChainRes, char **keyMetalList, int numKeyMetal, int min_metalBoundRes, int max_metalBoundRes, int min_numHCRes, bool isUsingTotalBoundRes, char* resList, const char *modmfilelistfile, int type_modm, const char* modmpath, const char *metalProFile,  int neg_filter, const char *vectorFile, const char *labelFile, int level, int encoding_scheme)/*{{{*/
/*****************************************************************************
 * encoding_scheme: added 2007-06-28
 * return the number of vectors
 ****************************************************************************/
{
    //#define K          5       // extension beyond binding pair
    //#define W          5       // extension within binding pair
    //#define P          24      // modm profile(20)+statistical value(2)+indicator outof seq(1)+hydrophobicity(1)

    int i;
    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D; 

    char id[SIZE_CHAIN_ID+1] = "";
    
    //  get zinc binding residues for each zinc binding protein, Metal3 or Metal4/*{{{*/
    // metalPro   -- store zinc binding residues in (Metal3 and Metal4) within resList and filtered by cutoff_score2
    // metalPro1  -- store zinc binding residues in (Metal3 and Metal4) within resList
    // metalPro2  -- store all zinc binding residues in (Metal3 and Metal4)
    using namespace MetalBindingProtein;
    /*this is the single site predictor, zinc protein using strucutre MetalPro*/
    int numMetalPro;
    int numMetalPro1;
    int numMetalPro2;
    int numMetalRes;
    int numMetalRes1;
    int numMetalRes2;
    Array1D <MetalPro> metalPro_1darray (MAX_NUM_METALPRO);
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
    if(level == 0)
    {
        GetMetalBoundRes3(metalProFile, metalPro, metalPro1, metalPro2, numMetalPro, numMetalPro1, numMetalPro2, numMetalRes, numMetalRes1, numMetalRes2, cutoff_score2, resList, modmpath, isExcludeOtherChainRes, keyMetalList, numKeyMetal, min_metalBoundRes, max_metalBoundRes, isUsingTotalBoundRes);
    }
    else
    {
        GetMetalBoundRes2(metalProFile, metalPro1, metalPro2, numMetalPro1, numMetalPro2, numMetalRes1, numMetalRes2, resList, modmpath, isExcludeOtherChainRes, keyMetalList, numKeyMetal, min_metalBoundRes, max_metalBoundRes, isUsingTotalBoundRes); 
    }
/*}}}*/
    MetalPro *metalProLevel = NULL;
    int numMetalProLevel = 0;
    // Get metalProLevel according to the level/*{{{*/
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
    }
    else
    { }/*}}}*/

    Array2D <char > metalProIDList_2darray(numMetalProLevel, SIZE_CHAIN_ID);
    char ** metalProIDList = metalProIDList_2darray.array2D;
    Array2D <char> metalProIDList_sorted_2darray(numMetalProLevel,SIZE_CHAIN_ID+1);
    char **metalProIDList_sorted = metalProIDList_sorted_2darray.array2D;
    Array1D <int> idxSortMetalPro_1darray(numMetalProLevel);
    int *idxSortMetalPro = idxSortMetalPro_1darray.array1D;

    for(i = 0 ; i < numMetalProLevel ; i++) idxSortMetalPro[i] = i;
   
    for(i = 0 ; i < numMetalProLevel; i++) { my_strcpy(metalProIDList[i], metalProLevel[i].id, SIZE_CHAIN_ID); }
    QuickSort_String(idxSortMetalPro, metalProIDList, 0, numMetalProLevel-1);
    for(i = 0 ; i < numMetalProLevel; i ++) { my_strcpy(metalProIDList_sorted[i], metalProIDList[idxSortMetalPro[i]], SIZE_CHAIN_ID);}

    bool *isBoundArray = NULL;
    FILE *fpVector ;
    FILE *fpLabel;
    fpVector = fopen(vectorFile,"w");
    fpLabel = fopen(labelFile,"w");
    checkfilestream(fpVector,vectorFile,"w");
    checkfilestream(fpLabel,labelFile,"w");

    
    int vectorDim  = 0;
    vectorDim = GetVectorDimension(K, W, numSite, P, encoding_scheme, isAddPairDist);

    char str1[2*MAX_PATH+1] = "";
    char str2[2*MAX_PATH+1] = "";

    MODM modm;
    MODM *pMODM;
    InitMODM(&modm);
    AllocMODM(&modm, LONGEST_SEQ, NUM_BLOSUM);
    modm.type_modm = type_modm;

    //char command[600+1] = "";

    
    Array1D <bool> isBoundArray_1darray(LONGEST_SEQ);
    isBoundArray = isBoundArray_1darray.array1D;
    // each residue pair is represented by a vector of size (2k+2+1)*p+1, where
    // k and p are the same as described in the single site vector, the second
    // 2 is for residues consisting the pair, the first 1 is for the unified
    // residues and the second 1 indicates the number of residues separating
    // the residues of the pair. 
    WriteVectorHeader(vectorDim, fpVector);
    WriteLabelHeader(fpLabel);

    FILE *fpFileList = fopen(modmfilelistfile, "r");
    checkfilestream(fpFileList, modmfilelistfile, "r");
    char modmfilepath[MAX_PATH+1] = "";
    int numVector = 0;

    
    while((linesize = fgetline(fpFileList, line , maxline)) != EOF)
    {
        if(linesize <= 0) continue;
        my_strcpy(modmfilepath, line, MAX_PATH);

        rootname(modmfilepath,str1); // get ID for the record
        SpanExcluding(str1,str2," ");
        my_strcpy(id,str2,SIZE_CHAIN_ID);
        StdID(id);


        modm.length = GetMODM(modmfilepath, modm.M, modm.alphabetMODM, modm.aaSeq, modm.score1, modm.score2);
        my_strcpy(modm.id,id,SIZE_CHAIN_ID);
        pMODM = &modm;

#ifdef DEBUG
        //debug start
		int totalCH = 0;
        for(int kk = 0 ; kk < pMODM->length; kk ++)
        {
            if(pMODM->aaSeq[kk] == 'C' || pMODM->aaSeq[kk] == 'H')
                totalCH ++;
        }
        //debug end
#endif
        GetIsBoundArray(metalProIDList_sorted, idxSortMetalPro, metalProLevel, numMetalProLevel, pMODM, isBoundArray);

        numVector += CreateVector_Metal( K,  W,  P,  numSite,  vectorDim,  cutoff_score2,  cutoff_consv,  cutoff_window_pair,  cutoff_window3,  cutoff_window4,  cutoff_min_window_pair,  cutoff_transcore_Metal3,  cutoff_transcore_Metal4,  min_numHCRes, resList, pMODM,  fpVector,  fpLabel,  neg_filter, isBoundArray, encoding_scheme);

    }

#ifdef DEBUG
    printf("totalCH = %d, at line %d of file %s\n", totalCH, __LINE__, __FILE__);
#endif
    fclose(fpFileList);
    fclose(fpVector);
    fclose(fpLabel);


    DeleteMODM(&modm,LONGEST_SEQ);
    for(i = 0; i < MAX_NUM_METALPRO ; i++)
    {
        DeleteMetalPro(&metalPro[i]);
        DeleteMetalPro(&metalPro1[i]);
        DeleteMetalPro(&metalPro2[i]);
    }
    return numVector;
}
/*}}}*/

int main(int argc, char** argv)/*{{{*/
{
    if( argc < 2 ) {
        PrintHelp();
        return -1;
    }

    int errsig = 0 ; /*error signal, 0 --> normal, no error*/

    int    i,j;
    int    mode                     = PRED;
    int    type_bond                = METAL;
    int    numSite                  = 1;
    int    type_modm                = MODM_LOG;
    int    level                    = 1 ; // added 2007-06-21, for fast loading MetalBoundRes

    double cutoff_score2            = 0.01;
    double cutoff_consv             = 0.70;
    int    cutoff_window_pair       = 150;
    int    cutoff_window3           = 300;
    int    cutoff_window4           = 400;
    int    cutoff_min_window_pair   = 20;
    double cutoff_transcore_Metal3     = 0.65;
    double cutoff_transcore_Metal4     = 1.0;
    int    neg_filter               = 1;
    bool   isExcludeOtherChainRes   = true;
    bool   isUsingTotalBoundRes     = true;

    char   resList[NUM_BLOSUM+1]    = "CHDE";
    int    max_metalBoundRes            = 5;
    int    min_metalBoundRes            = 3;
    int    min_numHCRes             = 2;
    int    K                        = 12;
    int    W                        = 5;
    int    P                        = 6;
    int    encoding_scheme          = 1;  // default encoding_scheme = 1, that is using the new encoding scheme now, 2007-06-28
    char   modmfile[MAX_PATH+1]     = "";
    char   modmfilelistfile[MAX_PATH+1] = "";
    char   vectorFile[MAX_PATH+1]   = "";
    char   labelFile[MAX_PATH+1]    = "";

    char   metalProFile[MAX_PATH+1] = "";
    char   ssbondProFile[MAX_PATH+1] = "";
    char   modmpath[MAX_PATH+1]    = "";

    char datadir[MAX_PATH+1] = "";
    GetDataDir(datadir);
    sprintf(metalProFile,"%s/revised-closeMetalPro.dat",datadir);
    sprintf(ssbondProFile,"%s/passe.ssbond",datadir);
    sprintf(modmpath,"%s/modm",datadir);
    

    char str[100+1] = "";

    int numKeyMetal = 0;
    char **keyMetalList = NULL;
    set <string> keyMetalList_set;
    set <string> ::iterator iss;
    // read in metal element list/*{{{*/
    int numMetalEle = 0 ;
    Array1D <Element> metalEle_1darray(NUM_METAL_ELEMENT);
    Element * metalEle = metalEle_1darray.array1D;
    Array2D <char> metalEleList_2darray(NUM_METAL_ELEMENT, SIZE_ATOM_ELEMENT+1);
    char **metalEleList = metalEleList_2darray.array2D;
    numMetalEle = GetMetalElementList( metalEle, metalElementListFile);
    for(i = 0 ; i < numMetalEle; i++)
        my_strcpy(metalEleList[i], metalEle[i].name, SIZE_ATOM_ELEMENT);
    /*}}}*/

    bool isNonOptionArg = false;
    i = 1;
    while(i < argc)/*{{{*/
    {
        if(argv[i][0] == '-' && !isNonOptionArg)
        {
            isNonOptionArg = false;
            if( strcmp(argv[i],"-h") == 0 || strcmp(argv[i],"--help") == 0) {
                PrintHelp();
                return 0;
            } else if(strcmp(argv[i],"--mode") == 0) {
                if(strcasecmp(argv[i+1],"pred") == 0 ){
                    mode = PRED;
                } else if(strcasecmp(argv[i+1],"train") == 0 ){
                    mode = TRAIN;
                } else {
                    fprintf(stderr,"Error! type pred or train for argument %s\n", "--mode");
                    errsig = -1; break;
                }
                i += 2;
            } else if(strcmp(argv[i],"-t") == 0) {
                if(strncasecmp(argv[i+1],"METAL", 1) == 0 ){
                    type_bond = METAL;
                } else if(strncasecmp(argv[i+1],"SS",1) == 0 ){
                    type_bond = SS;
                } else {
                    fprintf(stderr,"Error! type 'ss' or 'metal' for argument %s\n", "-t");
                    errsig = -1; break;
                }
                i += 2;
            } else if(strcmp(argv[i], "--keymetal") == 0) {
                for(j = i +1 ; j < argc ; j++) {
                    if(argv[j][0] != '-') {
                        my_strcpy(str, argv[j], SIZE_ATOM_ELEMENT); my_strupr(str);
                        if(BinarySearch_String(str, metalEleList, numMetalEle) != -1) {
                            keyMetalList_set.insert(str);
                        } else {
                            fprintf(stderr,"Error! '%s' is not metal element\n", argv[j]);
                            return -1;
                        }
                    } else {
                        break;
                    }
                }
                i = j;
            } else if(strcmp(argv[i],"-ns") == 0) {
                numSite = atoi(argv[i+1]);
                i += 2;
            } else if(strcmp(argv[i],"--vector") == 0) {
                my_strcpy(vectorFile, argv[i+1], MAX_PATH);
                i += 2;
            } else if(strcmp(argv[i],"--label") == 0) {
                my_strcpy(labelFile, argv[i+1], MAX_PATH);
                i += 2;
            } else if(strcmp(argv[i],"--level") == 0) {
                level = atoi(argv[i+1]);
                if(type_bond == METAL) {
                    if(level < 0 || level > 2) {
                        fprintf(stdout,"warning! level %d out of range [0-2]\n", level);
                        return -1;
                    }
                } else if( type_bond == SS) {
                    if(level < 0 || level > 1) {
                        fprintf(stdout,"warning! level %d out of range [0-1]\n", level);
                        return -1;
                    }
                }
                i += 2;
            } else if(strcmp(argv[i],"--modm-type") == 0) {
                if(strncasecmp(argv[i+1],"LOG", 1) == 0 ){
                    type_modm = MODM_LOG;
                } else if(strncasecmp(argv[i+1],"PER", 1) == 0 ){
                    type_modm = MODM_PER;
                } else {
                    fprintf(stderr,"Error! type LOG or PER for argument %s\n",
                            "--modm-type");
                    errsig = -1; break;
                }
                i += 2;
            } else if(strcmp(argv[i],"--profile") == 0) {
                if(strncasecmp(argv[i+1],"pssm", 1) == 0 ) {
                    profile_encoding_type = PSSM_PROFILE_ENCODE;
                } else if(strncasecmp(argv[i+1],"seq", 1) == 0 ){
                    profile_encoding_type = SEQUENCE_BINARY_ENCODE;
                } else {
                    fprintf(stderr,"Error! type PSSM or SEQ for argument %s\n", 
                            "--encoding");
                    errsig = -1; break;
                }
                i += 2;
            } else if(strcmp(argv[i],"-s") == 0) {
                cutoff_score2 = atof(argv[i+1]);
                i += 2;
            } else if(strcmp(argv[i],"-c") == 0) {
                cutoff_consv = atof(argv[i+1]);
                i += 2;
            } else if(strcmp(argv[i],"-a") == 0) {
                my_strcpy(resList,argv[i+1],NUM_BLOSUM);
                i += 2;
            } else if(strcmp(argv[i],"--min-hcres") == 0) {
                min_numHCRes = atoi(argv[i+1]);
                i += 2;
            } else if(strcmp(argv[i],"--modmpath") == 0) {
                my_strcpy(modmpath, argv[i+1], MAX_PATH);
                i += 2;
            } else if(strcmp(argv[i],"--win-pair") == 0) {
                cutoff_window_pair = atoi(argv[i+1]);
                i += 2;
            } else if(strcmp(argv[i],"--win-min-pair") == 0) {
                cutoff_min_window_pair = atoi(argv[i+1]);
                i += 2;
            } else if(strcmp(argv[i],"--win3") == 0) {
                cutoff_window3 = atoi(argv[i+1]);
                i += 2;
            } else if(strcmp(argv[i],"--win4") == 0) {
                cutoff_window4 = atoi(argv[i+1]);
                i += 2;
            } else if(strcmp(argv[i],"--ts3") == 0) {
                cutoff_transcore_Metal3 = atof(argv[i+1]);
                i += 2;
            } else if(strcmp(argv[i],"--ts4") == 0) {
                cutoff_transcore_Metal4 = atof(argv[i+1]);
                i += 2;
            }
            //else if(strcmp(argv[i],"-m") == 0)
            //{
            //my_strcpy(metalProFile, argv[i+1], MAX_PATH);
            //i += 2;
            //}
            else if(strcmp(argv[i],"-K") == 0) {/* parameters for vector window*/
                K = atoi(argv[i+1]);
                i += 2;
            } else if(strcmp(argv[i],"-W") == 0){/* parameters for vector window*/
                W = atoi(argv[i+1]);
                i += 2;
            } else if(strcmp(argv[i],"-P") == 0){/* parameters for vector window*/
                P = atoi(argv[i+1]);
                i += 2;
            } else if(strcmp(argv[i],"--encode") == 0){/* parameters encoding_scheme*/ 
                encoding_scheme = atoi(argv[i+1]);
                if(encoding_scheme < 0 || encoding_scheme > 1) {
                    fprintf(stderr, "Error! encoding_scheme,can be 0 or 1 only\n");
                    return -1;
                }
                i += 2;
            } else if(strcmp(argv[i],"--bound-res") == 0) {
                if(IsNumeric(argv[i+1])){
                    min_metalBoundRes = atoi(argv[i+1]);
                } else {
                    fprintf(stderr,"Argument Error! two integer values should follow argument: %s\n", "--bound-res");
                    errsig = -1; break;
                }

                if(IsNumeric(argv[i+2])){
                    max_metalBoundRes = atoi(argv[i+2]);
                } else {
                    fprintf(stderr,"Argument Error! two integer values should follow argument: %s\n", "--bound-res");
                    errsig = -1; break;
                }
                i += 3;
            }
            //else if(strcmp(argv[i],"-i") == 0)
            //{
            //my_strcpy(seqfile, argv[i+1], MAX_PATH);
            //my_strcpy(modmfile, argv[i+1], MAX_PATH);
            //i += 3;
            //}
            else if(strcmp(argv[i],"--metal") == 0) {
                my_strcpy(metalProFile, argv[i+1], MAX_PATH);
                i += 2;
            } else if(strcmp(argv[i],"-metaltrans")==0 ||
                    strcmp(argv[i],"--metaltrans") == 0) {
                my_strcpy(metalTranMatrixDir, argv[i+1], MAX_PATH);
                i += 2;
            } else if(strcmp(argv[i],"--ssbond") == 0) {
                my_strcpy(ssbondProFile, argv[i+1], MAX_PATH);
                i += 2;
            }
            //else if(strcmp(argv[i],"--write-vector") == 0)
            //{
            //isWriteVectorVectorFile = true;
            //i += 1;
            //}
            else if(strcmp(argv[i],"--not-mask-polyhis") == 0) {
                isMaskPolyHis = false;
                i += 1;
            } else if(strcmp(argv[i],"--not-use-ic") == 0) {
                isUseConsvI = false;
                i += 1;
            } else if(strcmp(argv[i],"--not-use-shcr") == 0) {
                isUseSHCRRule = false;
                i += 1;
            } else if(strcmp(argv[i],"--use--sc") == 0) {
                isUsingTotalBoundRes = false;
                i += 1;
            } else if(strcmp(argv[i],"--not-add-pair-dist") == 0) {
                isAddPairDist = false;
                i += 1;
            } else if(strcmp(argv[i],"--neg-filter") == 0) {
                neg_filter = atoi(argv[i+1]);
                i += 2;
            } else if(strcmp(argv[i],"-l") == 0 || strcmp(argv[i],"--l") == 0
                    || strcmp(argv[i], "-list") == 0) {
                my_strcpy(modmfilelistfile, argv[i+1], MAX_PATH);
                i += 2;
            } else if(strcmp(argv[i], "--") == 0) {
                isNonOptionArg = true;
                i ++;
            } else {
                fprintf(stderr,"\nError! argument '%s' is not in option list, type --help to check\n", argv[i]);
                return -1;
            }
        } else { 
            my_strcpy(modmfile, argv[i], MAX_PATH);
            i += 1;
            //fprintf(fptmp,"%s\n", argv[i]);
            //PrintHelp();
            //errsig = -1; break;
        }
    }/*}}}*/
    if(errsig != 0) {
        return errsig;
    }

    if(encoding_scheme == 1){ /* check P if the using the new encoding_scheme 2007-06-28*/
        if(P > 6) {
            fprintf(stderr,"Warning!, for the new encoding_scheme, P should be in [1-6]\n");
            return -1;
        }
    }
    numKeyMetal = keyMetalList_set.size();
    Array2D <char> keyMetalList_2darray(numKeyMetal, SIZE_ATOM_ELEMENT+1);
    keyMetalList = keyMetalList_2darray.array2D;
    i = 0 ;
    for(iss = keyMetalList_set.begin(); iss != keyMetalList_set.end(); iss ++) {
        my_strcpy(keyMetalList[i], (*iss).c_str(), SIZE_ATOM_ELEMENT);
        i ++;
    }

    if(mode == PRED) {
        if(strcmp(vectorFile,"") == 0) {
            fprintf(stderr,"vector file not set. Please set by the option -vector. Exit...");
            return -1 ;
        }
        set <string> modmfileset;
        modmfileset.clear();
        
        if(strcmp(modmfile,"") != 0) {
            modmfileset.insert(modmfile);
        }
        if (strcmp(modmfilelistfile, "") != 0) {
            FILE *fp = NULL;
            fp = fopen(modmfilelistfile, "r");
            if (fp!= NULL){
                int linesize;
                int maxline = 300;
                Array1D <char> line_1darray(maxline+1);
                char *line = line_1darray.array1D; 
                while((linesize=fgetdelim(fp, line, WHITE_SPACE, maxline)) != EOF) {
                    if (linesize > 0){
                        //fprintf(stderr,"line=%s\n", line); [>debug<]
                        modmfileset.insert(line);
                    }
                }
                fclose(fp);
            }else{
                fprintf(stderr,"Can not open modmfilelistfile %s for read. Exit...\n", modmfilelistfile);
            }
        }

        if (modmfileset.size() <= 0){
            fprintf(stderr,"No modmfile is set. Exit...\n");
            return 1;
        } else{
            if(type_bond == METAL) {
                //fprintf(stderr,"nummodmfile=%d\n", modmfileset.size());
                //set <string> ::iterator iss;
                //for (iss = modmfileset.begin();iss != modmfileset.end();iss++){
                    //fprintf(stderr,"iss=<%s>\n", (*iss).c_str());
                //}
                CreateVectorPredict_Metal(K, W, P, numSite, cutoff_score2, cutoff_consv, cutoff_window_pair, cutoff_window3, cutoff_window4, cutoff_min_window_pair, cutoff_transcore_Metal3, cutoff_transcore_Metal4, min_numHCRes, resList, modmfileset, type_modm, vectorFile, encoding_scheme);
            } else if(type_bond == SS) {
                fprintf(stderr,"create_svm_vector for SS-bond has not been implemented yet. Exit...\n");
                return 1;
                //CreateVectorPredict_SSBond(K, W, P, numSite, cutoff_score2, cutoff_consv, cutoff_window_pair, min_numHCRes, resList, modmfile, type_modm, vectorFile, isAddPairDist, encoding_scheme);
            }
        }
    } else if(mode == TRAIN) {
        if(strcmp(modmfilelistfile,"") == 0) {
            fprintf(stderr,"Error! modmfilelistfile not set\n");
            return -1;
        }

        char rtname[MAX_PATH+1] = "";
        rootname(modmfilelistfile,rtname);
        if(strcmp(vectorFile,"") == 0) {
            sprintf(vectorFile,"%s-ns%d.vector",rtname, numSite);
        }
        if(strcmp(labelFile,"") == 0) {
            sprintf(labelFile,"%s-ns%d.label",rtname, numSite);
        }
            
        if(type_bond == METAL) {
            CreateVectorTrain_Metal(K, W, P, numSite, cutoff_score2, cutoff_consv, cutoff_window_pair, cutoff_window3, cutoff_window4, cutoff_min_window_pair, cutoff_transcore_Metal3, cutoff_transcore_Metal4, isExcludeOtherChainRes, keyMetalList, numKeyMetal, min_metalBoundRes, max_metalBoundRes, isUsingTotalBoundRes, min_numHCRes, resList, modmfilelistfile, type_modm, modmpath, metalProFile,  neg_filter, vectorFile, labelFile, level, encoding_scheme);

        } else if(type_bond == SS) {
            fprintf(stderr,"create_svm_vector for SS-bond has not been implemented yet. Exit...\n");
            return 1;
            //CreateVectorTrain_SSBond(K, W, P, numSite, cutoff_score2, cutoff_consv, cutoff_window_pair, min_numHCRes, resList, modmfilelistfile, type_modm, modmpath,ssbondProFile,  neg_filter, vectorFile, labelFile, level, isAddPairDist, encoding_scheme);
        }
    }
    return EXIT_SUCCESS;
}
/*}}}*/
