#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cstring>
#include <iostream>
#include <cmath>
#include <set>
#include <algorithm>
#include "array.h"
#include "mytemplate.h"
#include "myfunc.h"
#include "mypro.h"
#include "subset.h"

using namespace std;

int DIGIT_INDEL = -9;
int GAP_DIGIT = DIGIT_INDEL;

char CHAR_INDEL = '-' ;// character representation of insertion and deletion
char GAP_CHAR = CHAR_INDEL;
char CHAR_VECTOR_ID_SEPRATOR =  ';'; //the separator for items in the record id of svm vectors

/*!! IMPORTANT!! extern defined in the header file, and the acutally instantiation is carried
 * out in here in the source file, 2010-04-20 */
const char *STD3CharAA_alphabet[] = /*{{{*/
{
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "UNK"
};
/*}}}*/
const char UNKNOWN_3CharAA[]= "UNK";
const char STD1CharAA_alphabet[] = "ARNDCQEGHILKMFPSTWYVX"; 

const char *BLOSUM3D_alphabet[] = /*{{{*/
// not ordered
{
    "ALA",//A
    "ARG",//R
    "ASN",//N
    "ASP",//D
    "CYS",//C
    "GLN",//Q
    "GLU",//E
    "GLY",//G
    "HIS",//H
    "ILE",//I
    "LEU",//L
    "LYS",//K
    "MET",//M
    "PHE",//F
    "PRO",//P
    "SER",//S
    "THR",//T
    "TRP",//W
    "TYR",//Y
    "VAL",//V
    "ASX",//X
    "GLX",
    "UNK",
    "UNK"
};
/*}}}*/
const char BLOSUM1D_alphabet[] = "ARNDCQEGHILKMFPSTWYVBZX*";
const char AAAlphabet_Tuping[] = "AVLIPFMKRHGSTCYNEWDQ" ;
const int MAPINDEX_BLOSUM_Tuping[] = {0, 19, 10, 9, 14, 13, 12, 11, 1, 8, 7, 15, 16, 4, 18, 2, 6, 17, 3, 5 };//index that map AAAlphabet_Tuping to BLOSUM1D_alphabet

const int AAS_Code [] = /*{{{*/
{
    0,//A--ALA
    20,//B--Non
    13,//C--CYS
    18,//D--ASP
    16,//E--GLU
    5,//F--PHE
    10,//G--GLY
    9,//H--HIS
    3,//I--ILE
    20,//J--Non
    7,//K--LYS
    2,//L--LEU
    6,//M--MET
    15,//N--ASN
    20,//O--Non
    4,//P--PRO
    19,//Q--GLN
    8,//R--ARG
    11,//S--SER
    12,//T--THR
    20,//U--Non
    1,//V--VAL
    17,//W--TRP
    20,//X--Non
    14,//Y--TYR
    20//Z--Non
};/*}}}*/

const char *PDB3Char_AllRes [] = // all recognizable 3-letter residues in PDB, ordered alphabetically/*{{{*/
//  from Bio/SCOP/Raf.py, and scop 
{
        "143", //C
        "1LU", //L
        "2AS", //D
        "2LU", //L
        "2MR", //R
        "3AH", //H
        "5HP", //E
        "AAR", //R
        "ABA", //C
        "ABA", //N
        "ACA", //A
        "ACL", //R
        "ACY", //G
        "AEI", //T
        "AGM", //R
        "AIB", //A
        "ALA", //A
        "ALM", //A
        "ALO", //T
        "ALS", //C
        "ALY", //K
        "ans", //la
        "ARG", //R
        "ARM", //R
        "ARO", //R
        "ASA", //D
        "ASB", //D
        "ASI", //N
        "ASK", //D
        "ASL", //D
        "ASN", //N
        "ASP", //D
        "ASQ", //D
        "ASX", //B
        "AYA", //A
        "BCS", //C
        "BFD", //D
        "BHD", //D
        "BMT", //T
        "BNN", //A
        "BTA", //E
        "BTC", //G
        "BUC", //C
        "BUG", //L
        "C5C", //C
        "C6C", //C
        "CAF", //C
        "CAS", //C
        "CAY", //C
        "CCS", //C
        "CCY", //S
        "CEA", //C
        "CGU", //E
        "CHG", //A
        "CLB", //S
        "CLD", //S
        "CLE", //L
        "CME", //C
        "CMT", //C
        "CR2", //G
        "CRG", //S
        "CRO", //G
        "CRO", //S
        "CRQ", //Q
        "CSA", //C
        "CSB", //C
        "CSD", //A
        "CSE", //C
        "CSH", //S
        "CSO", //C
        "CSP", //C
        "CSR", //C
        "CSS", //C
        "CSW", //C
        "CSX", //C
        "CSY", //S
        "CSZ", //C
        "CXM", //M
        "CY1", //C
        "CY3", //C
        "CY4", //C
        "CYF", //C
        "CYG", //C
        "CYM", //C
        "CYQ", //C
        "CYS", //C
        "CZZ", //C
        "DAH", //F
        "DAL", //A
        "DAR", //R
        "DAS", //D
        "DCY", //C
        "DGL", //E
        "DGN", //Q
        "DHA", //A
        "DHI", //H
        "DIL", //I
        "DIV", //V
        "DLE", //L
        "DLY", //K
        "DNP", //A
        "DOH", //D
        "DPN", //F
        "DPR", //P
        "DSN", //S
        "DSP", //D
        "DTH", //T
        "DTR", //W
        "DTY", //Y
        "DVA", //V
        "EFC", //C
        "EHP", //F
        "ESD", //S
        "FGL", //C
        "FLA", //A
        "FME", //M
        "FTR", //W
        "GGL", //E
        "GL3", //G
        "GLH", //E
        "GLN", //Q
        "GLU", //E
        "GLX", //Z
        "GLY", //G
        "GLZ", //G
        "GMA", //E
        "GPL", //K
        "GSC", //G
        "H2P", //H
        "HAC", //A
        "HAR", //R
        "HIC", //H
        "HIP", //H
        "HIS", //H
        "HMR", //R
        "HPH", //F
        "HPQ", //F
        "HTR", //W
        "HYP", //P
        "IAS", //D
        "IAS", //N
        "IIC", //S
        "IIL", //I
        "ILE", //I
        "INI", //K
        "IYR", //Y
        "KCX", //K
        "LEU", //L
        "LLP", //K
        "LLY", //K
        "LTR", //W
        "LYM", //K
        "LYN", //K
        "LYS", //K
        "LYX", //K
        "LYZ", //K
        "M3L", //K
        "MAA", //A
        "MEN", //N
        "MET", //M
        "MGN", //Q
        "MHO", //M
        "MHS", //H
        "MIS", //S
        "MLE", //L
        "MLY", //E
        "MLY", //K
        "MLZ", //K
        "MPQ", //G
        "MSA", //G
        "MSE", //M
        "MSO", //M
        "MTY", //F
        "MVA", //V
        "NC1", //S
        "NEM", //H
        "NEP", //H
        "NIY", //Y
        "NLE", //L
        "NLN", //L
        "NLP", //L
        "NMC", //G
        "NPH", //C
        "OAS", //S
        "OCS", //C
        "OCY", //C
        "OMT", //M
        "PAQ", //Y
        "PBI", //S
        "PCA", //E
        "PEC", //C
        "PHD", //D
        "PHE", //F
        "PHI", //F
        "PHL", //F
        "PIA", //A
        "PR3", //C
        "PRO", //P
        "PRR", //A
        "PRS", //P
        "PTR", //Y
        "PVL", //S
        "PYA", //A
        "PYX", //C
        "S1H", //S
        "SAC", //S
        "SAR", //G
        "SBD", //S
        "SBL", //S
        "SCH", //C
        "SCS", //C
        "SCY", //C
        "SEB", //S
        "SEC", //C
        "SEG", //S
        "SEL", //S
        "SEP", //S
        "SER", //S
        "SET", //S
        "SHC", //C
        "SHR", //K
        "SMC", //C
        "SME", //M
        "SNC", //C
        "SNN", //D
        "SOC", //C
        "STY", //Y
        "SUI", //D
        "SVA", //S
        "THR", //T
        "TIH", //A
        "TPL", //W
        "TPO", //T
        "TPQ", //A
        "TRF", //W
        "TRG", //K
        "TRN", //W
        "TRO", //W
        "TRP", //W
        "TRQ", //W
        "TRW", //W
        "TYB", //Y
        "TYI", //Y
        "TYN", //Y
        "TYQ", //Y
        "TYR", //Y
        "TYS", //Y
        "TYY", //Y
        "UNK", //X
        "VAL", //V
        "YCM", //C
        "YOF"  //Y
};/*}}}*/
const char PDB1Char_AllRes[] = "CLDLRHERCNARGTRAAATCKlRRRDDNDDNDDBACDDTAEGCLCCCCCCSCEASSLCCGSGSQCCACSCCCCCCSCMCCCCCCCCCFARDCEQAHIVLKADFPSDTWYVCFSCAMWEGEQEZGGEKGHARHHHRFFWPDNSIIKYKLKKWKKKKKKANMQMHSLEKKGGMMFVSHHYLLLGCSCCMYSECDFFFACPAPYSACSSGSSCCCSCSSSSSCKCMCDCYDSTAWTAWKWWWWWYYYYYYYXVCY";

const char *PDB3CharAA_alphabet[] = //recognized amino acid residue code in pdb file, ordered /*{{{*/
{ 
    "ACE",
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "ASX",
    "CYS",
    "FOR",
    "GLN",
    "GLU",
    "GLX",
    "GLY",
    "HIS",
    "HYP",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PCA",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "UNK", 
    "VAL"
}; /*}}}*/
const char PDB1CharAA_alphabet[] = "XARNDBCXQEZGHXILKMFPXSTWYXV";

const double Hydrophobicity[] = { /*{{{*/
    0.616, /* A */
    0.000, /* R */
    0.236, /* N */
    0.028, /* D */
    0.680, /* C */
    0.251, /* Q */
    0.043, /* E */
    0.501, /* G */
    0.165, /* H */
    0.943, /* I */
    0.943, /* L */
    0.283, /* K */
    0.738, /* M */
    1.00,  /* F */
    0.711, /* P */
    0.359, /* S */
    0.450, /* T */
    0.878, /* W */
    0.880, /* Y */
    0.825, /* V */
    0.000  /* X */
}; /*  The hydrophobicities given are the "Scaled" values from
# computational log(P) determinations by the "Small Fragment Approach" (see,
# "Development of Hydrophobicity Parameters to Analyze Proteins Which Bear
# Post- or Cotranslational Modifications" Black, S.D. and Mould, D.R. (1991)
# Anal. Biochem. 193, 72-82). The equation used to scale raw log(P) values to
# the scaled values given is as follows: Scaled Parameters = (Raw Parameters +
# 2.061)/4.484 .
ordered by STD1CharAA_alphabet, the last one is for X
*//*}}}*/
const double Background_AA_Freq[] = { 0.078, 0.051, 0.045, 0.054, 0.019, 0.043, 0.063, 0.074, 0.022, 0.052, 0.09,0.057, 0.022, 0.039, 0.052, 0.071, 0.059, 0.013, 0.032, 0.064 } ; // ordered by BLOSUM1D_alphabet, from sdsc
const double Ln_Background_AA_Freq[] = { -2.551, -2.976, -3.101, -2.919, -3.963, -3.147, -2.765, -2.604, -3.817, -2.957, -2.408, -2.865, -3.817, -3.244, -2.957, -2.645, -2.830, -4.343, -3.442, -2.749 };//ln of Background_AA_Freq, ordered by BLOSUM1D_alphabet
#define NUM_SHAPES 9   //scalar of shag_matrix
const char SHAPE_alphabet[] = "AKSRUVTG-";

int numNonMetalHetGroup = 60;
const char *nonMetalHetGroup[] =/*{{{*/
// ordered
{
    "101",
    "12A",
    "1AR",
    "1GL",
    "2AS",
    "2GL",
    "3AA",
    "3AT",
    "3DR",
    "3PO",
    "6HA",
    "6HC",
    "6HG",
    "6HT",
    "A26",
    "AA6",
    "ABD",
    "AC1",
    "ACO",
    "AIR",
    "AMU",
    "AMX",
    "AP5",
    "AMG",
    "APU",
    "B9A",
    "BCA",
    "BNA",
    "CAA",
    "CBS",
    "CGS",
    "CMC",
    "CND",
    "CO8",
    "COA",
    "COF",
    "COS",
    "DCA",
    "DGD",
    "FAB",
    "FAD",
    "FAG",
    "FAM",
    "FDA",
    "GPC",
    "IB2",
    "NAD",
    "NAH",
    "NAI",
    "NAL",
    "NAP",
    "NBD",
    "NDP",
    "PAD",
    "SAD",
    "SAE",
    "T5A",
    "TRE",
    "UP5",
    "ZID"
};/*}}}*/

int EncodeMatrixScore1[][SIZE_ENCODE_SCORE1] =
/* consider score1 is independt, that is residues with high score1 is not
 * related to the residue with low score1, so the encoding vector are
 * perpendicular*/
{
    { 4,0,0,0,0},//<0.2    
    { 0,4,0,0,0},//0.2~0.5  
    { 0,0,4,0,0},//0.5~1.0  
    { 0,0,0,4,0},//1.0~2.0  
    { 0,0,0,0,4} //2.0~+    
};
double EncodeScaleScore1[] =
{
    0.2,
    0.5,
    1.0,
    2.0 
};

#define SIZE_ENCODE_SCORE2 5
int EncodeMatrixScore2[][SIZE_ENCODE_SCORE2] = 
{
    { 4,0,0,0,0},//0~0.05    
    { 0,4,0,0,0},//0.05~0.25  
    { 0,0,4,0,0},//0.25~0.8   
    { 0,0,0,4,0},//0.8~1.5     
    { 0,0,0,0,4} //1.5~+     
};
double EncodeScaleScore2[] =
{
    0.05,
    0.25,
    0.8,
    1.5 
};


#define SIZE_ENCODE_AA_CENTER 5 /*this is for residues of interest centered in a window*/
char EncodeMatrixAA_Center_alphabet[] = "CHDE";
int EncodeMatrixAA_Center[][SIZE_ENCODE_AA_CENTER] =
{
    {4 ,1 ,0 ,0,0},//C     
    {1 ,4 ,0 ,0,0},//H     
    {0 ,0 ,4 ,0,0},//D     
    {0 ,0 ,0 ,4,0},//E     
    {0 ,0 ,0 ,0,4} //others
};

#define SIZE_ENCODE_AA_NON_CENTER 5 /*this is for residues of interest centered in a window*/
char EncodeMatrixAA_Non_Center_alphabet[] = "CHDE";
int EncodeMatrixAA_Non_Center[][SIZE_ENCODE_AA_NON_CENTER] =
{
    {4 ,1 ,0 ,0,0},//C     
    {1 ,4 ,0 ,0,0},//H     
    {0 ,0 ,4 ,0,0},//D     
    {0 ,0 ,0 ,4,0},//E     
    {0 ,0 ,0 ,0,4} //others
};
/*int EncodeMatrixAA_Non_Center[][SIZE_ENCODE_AA_NON_CENTER] =[>{{{<]            */
/*{                                                                              */
/*    {3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, [>A<]     */
/*    {0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, [>R<]     */
/*    {0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, [>N<]     */
/*    {0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, [>D<]     */
/*    {0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, [>C<]     */
/*    {0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, [>Q<]     */
/*    {0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, [>E<]     */
/*    {0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, [>G<]     */
/*    {0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, [>H<]     */
/*    {0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, [>I<]     */
/*    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, [>L<]     */
/*    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0}, [>K<]     */
/*    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0}, [>M<]     */
/*    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0}, [>F<]     */
/*    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0}, [>P<]     */
/*    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0}, [>S<]     */
/*    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0}, [>T<]     */
/*    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0}, [>W<]     */
/*    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0}, [>Y<]     */
/*    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0}, [>V<]     */
/*    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3}  [>others<]*/
/*};                                                                             */
/*[>}}}<]                                                                        */

#define SIZE_ENCODE_WATERACC 5 /*wateracc is not allowed to be used in zn prediction, because this is derived from 3D structure, 2007-07-24*/
int EncodeMatrixWaterAcc[][SIZE_ENCODE_WATERACC] =
{
    {3 ,0 ,0 ,0,0},//0     
    {0 ,3 ,0 ,0,0},//1,2,3
    {0 ,0 ,3 ,0,0},//4,5     
    {0 ,0 ,0 ,3,0},//6,7,8     
    {0 ,0 ,0 ,0,3} //9
};
int EncodeScaleWaterAcc[] =
{
    0,
    3,
    5,
    8 
};


#define SIZE_ENCODE_DISTANCE 5
int EncodeMatrixDistance[][SIZE_ENCODE_DISTANCE] =
{
    {4 ,0 ,0 ,0,0},//1,2     
    {0 ,4 ,0 ,0,0},//3
    {0 ,0 ,4 ,0,0},//4~5
    {0 ,0 ,0 ,4,0},//6~20
    {0 ,0 ,0 ,0,4} //>20
};
int EncodeScaleDistance[] =
{
    2,
    3,
    5,
    20
};

int EncodeMatrixHydrophobicity[][SIZE_ENCODE_HYDROPHOBICITY] =
{                                                             
   { 4,0,0},//0~0.3    R, D, E, H, N, Q, K                            
   { 0,4,0},//0.3~0.65  S, T, G, A                 
   { 0,0,4} //0.7~0.9  C, P,M,V,W,Y,I,L,F                        
};                                                            
double EncodeScaleHydrophobicity[] =                          
{                                                             
   0.3,                                                     
   0.65,                                                      
};                                                            

namespace ZnBindingProtein
{
    int MAX_NUM_ZNPRO = 400;
    int MAX_NUMRES    = 30;  // maximal number of residues bind to zinc atom per chain
    int MAX_ATOMENV   = 30;
    int MAX_BOUND_SITE = 10;   // maximal number of coordinate bound
    int MAX_BOUND_METAL_PER_RES = 10; // maximal number of metal atoms bound to a residue, set to a large value, e.g 1AOO, Cys9 bind to 4 Ag(I) atoms, if set the cutoff_distance to a larger value, then, this value can be more than 5
}

namespace MetalBindingProtein
{
    int MAX_NUM_METALPRO = 1500;
    int MAX_NUMRES       = 200;  // maximal number of residues bind to zinc atom per chain
    int MAX_ATOMENV      = 50;
    int MAX_BOUND_SITE    = 10;   // maximal number of coordinate bound
    int MAX_BOUND_METAL_PER_RES = 10; // maximal number of metal atoms bound to a residue
}

/*functions for reading in structured file*/
/*functions for reading PDB ATOM record line*/
string int2string(int number)/*{{{*/
{
    stringstream ss;//create a stringstream
    ss << number;//add number to the stream
    return ss.str();//return a string with the contents of the stream
}/*}}}*/

void ScanfAtomSerial(const char *PDBAtomRecordLine, int &atomSerial)/*{{{*/
{
    char str[SIZE_LINE_PDB+1] = ""; int i;
    for ( i = 0 ; i < 5 ; i++) str[i] = PDBAtomRecordLine[i+6]; str[i] = '\0';
    sscanf(str,"%d",&atomSerial);//serial
}/*}}}*/
void ScanfAtomResName(const char *PDBAtomRecordLine, char *atomResName)/*{{{*/
{
    char str[SIZE_LINE_PDB+1] = "";  int i;
    for ( i=0 ;i < 3 ; i++) str[i] = PDBAtomRecordLine[i+17]; str[i] = '\0';
    sscanf(str,"%s",atomResName);//residue name
}/*}}}*/
void ScanfAtomOrigName(const char *PDBAtomRecordLine, char * atomOrigName)/*{{{*/
{
    char str[SIZE_LINE_PDB+1] = ""; int i = 0;
    for ( i=0 ;i<4;i++) str[i] = PDBAtomRecordLine[i+12]; str[i] = '\0';
    my_strcpy(atomOrigName,str,SIZE_ATOM_ORIGNAME);
}/*}}}*/

void ScanfCoorRecord_Atom(const char *line, Atom *pAtom)/*{{{*/
{
    InitAtom(pAtom);
    char str[30+1] = "";
    int i ;
    for ( i=0 ;i<6;i++) str[i] = line[i+0]; str[i] = '\0';
    sscanf(str,"%s",pAtom->recordID);

    for ( i=0 ;i<5;i++) str[i] = line[i+6]; str[i] = '\0';
    pAtom->serial = atoi(str);   // serial
    for ( i=0 ;i<4;i++) str[i] = line[i+12]; str[i] = '\0';
    sscanf(str,"%s",pAtom->name);// atom name
    my_strcpy(pAtom->origName,str,SIZE_ATOM_ORIGNAME); // record four colums for atom name record
    // ignore the end whitespaces

    pAtom->altLoc = line[16];

    for ( i=0 ;i<3;i++) str[i] = line[i+17]; str[i] = '\0';
    sscanf(str,"%s",pAtom->resName);//residue name

    pAtom->chainID = line[21];//chainID = ' ', when chainID is empty

    for ( i=0 ;i<4;i++) str[i] = line[i+22]; str[i] = '\0';
    pAtom->resSeq  = atoi(str); //residue seq, it's not necessary corresponding that  in SEQRES
    pAtom->iCode = line[26];

    for ( i=0 ;i<8;i++) str[i] = line[i+30]; str[i] = '\0';
    pAtom->x = atof(str);

    for ( i=0 ;i<8;i++) str[i] = line[i+38]; str[i] = '\0';
    pAtom->y = atof(str);

    for ( i=0 ;i<8;i++) str[i] = line[i+46]; str[i] = '\0';
    pAtom->z = atof(str);

    for ( i=0 ;i<6;i++) str[i] = line[i+54]; str[i] = '\0';
    pAtom->occupancy = atof(str);

    for ( i=0 ;i<6;i++) str[i] = line[i+60]; str[i] = '\0';
    pAtom->tempFactor = atof(str); // temperature factor

    for ( i=0 ;i<4;i++) str[i] = line[i+72]; str[i] = '\0';
    sscanf(str,"%s",pAtom->segID); // segment identifier

    for ( i=0 ;i<2;i++) str[i] = line[i+76]; str[i] = '\0';
    sscanf(str,"%s",pAtom->element); // element name

    for ( i=0 ;i<2;i++) str[i] = line[i+78]; str[i] = '\0';
    sscanf(str,"%s",pAtom->charge); // charge of atoms
}/*}}}*/
void ScanfCoorRecord_Atom_Simp1(const char* line, Atom *pAtom)/*{{{*/
    /******************************************************
     * read in part of the information of the ATOM record *
     *      serial                                        *
     *      altLoc                                        *
     *      name                                          *
     *      resName                                       *
     *      chainID                                       *
     *      resSeq                                        *
     *      iCode                                         *
     *      x                                             *
     *      y                                             *
     *      z                                             *
     *      occupancy                                     *
     ******************************************************/
{
    InitAtom(pAtom);
    char str[30+1] = "";
    int i ;
    for ( i=0 ;i<6;i++) str[i] = line[i+0]; str[i] = '\0';
    sscanf(str,"%s",pAtom->recordID);

    for ( i=0 ;i<5;i++) str[i] = line[i+6]; str[i] = '\0';
    sscanf(str,"%d",&pAtom->serial);//serial
    for ( i=0 ;i<4;i++) str[i] = line[i+12]; str[i] = '\0';
    sscanf(str,"%s",pAtom->name);// atom name
    my_strcpy(pAtom->origName,str, SIZE_ATOM_ORIGNAME); // record four colums for atom name record
    pAtom->altLoc = line[16];

    for ( i=0 ;i<3;i++) str[i] = line[i+17]; str[i] = '\0';
    sscanf(str,"%s",pAtom->resName);//residue name

    pAtom->chainID = line[21];//if chainID is empty, set chainID as NULLCHAR, ' '

    for ( i=0 ;i<4;i++) str[i] = line[i+22]; str[i] = '\0';
    sscanf(str,"%d",&pAtom->resSeq);//residue seq, it's not necessary corresponding that in SEQRES records
    pAtom->iCode = line[26];

    for ( i=0 ;i<8;i++) str[i] = line[i+30]; str[i] = '\0';
    pAtom->x = atof(str);    // x coordinate

    for ( i=0 ;i<8;i++) str[i] = line[i+38]; str[i] = '\0';
    pAtom->y = atof(str);         // y coordinate

    for ( i=0 ;i<8;i++) str[i] = line[i+46]; str[i] = '\0';
    pAtom->z = atof(str);         // z coordinate

    for ( i=0 ;i<6;i++) str[i] = line[i+54]; str[i] = '\0';
    pAtom->occupancy = atof(str);  // occupancy

}/*}}}*/
int Scanf_SEQRES_Record(const char* line,int resBegin, char* title, int& serNum,char& chainID,int& numRes,char **resName)  /*{{{*/
{
    int resEnd = resBegin;
    char  str[20];
    int i,j;
    sscanf(line,"%6s",title);

    for(i=0 ;i<2; i++) str[i] = line[i+8]; str[i]='\0';
    serNum = atoi(str);

    chainID = line[11];

    for(i=0 ;i<4; i++) str[i] = line[i+13]; str[i]='\0';
    numRes = atoi(str);

    for(i = 0 ; i < 13 ; i ++)
    {
        for(j = 0 ; j < 3 ; j ++)
            str[j] = line[i*4+19+j];
        str[j] = '\0';
        SpanExcluding(str,resName[i]);
    }

    for(i = 0 ; i < 13 ; i ++)
    {
        if(strcmp(resName[i],"") != 0)
            resEnd ++;
    }

    return resEnd;
}/*}}}*/
void Scanf_SEQRES_Para(const char* line,char* title, int& serNum,char& chainID,int& numRes, char* resName)/*{{{*/ 
{
    char  str[20+1] = "";
    int i;
    serNum  = 0;
    chainID = '\0';
    numRes  = 0;
    strcpy(title,"");
    strcpy(resName, "");

    sscanf(line,"%6s",title);

    for(i=0 ;i<2; i++) str[i] = line[i+8]; str[i]='\0';
    serNum = atoi(str);

    chainID = line[11]; //chainID=' ', when chainID is emtpy

    for(i=0 ;i<4; i++) str[i] = line[i+13]; str[i]='\0';
    numRes = atoi(str);

    for(i=0 ;i<3; i++) str[i] = line[i+19]; str[i]='\0';
    sscanf(str,"%3s",resName);
}/*}}}*/
int Scanf_SEQRES_Seq(const char* line,int resBegin, char **resName)/*{{{*/
    /****************************************************************************
     * read in the protein sequence in SEQRES record, 
     * give one line of SEQRES record,
     * and the begining index of chain add the residues in this SEQRES
     * record to the sequence and return the ending index
     ***************************************************************************/
{
    char str[20] = "";
    int resEnd = resBegin;
    int i , j;
    for(i = 0 ; i < 13 ; i ++)
        strcpy(resName[i],"");

    for(i = 0 ; i < 13 ; i ++)
    {
        for(j = 0 ; j < 3 ; j ++)
            str[j] = line[i*4+19+j];
        str[j] = '\0';
        sscanf(str,"%3s",resName[i]);
    }

    for(i = 0 ; i < 13 ; i ++)
    {
        if(strcmp(resName[i],"") != 0)
            resEnd ++;
    }

    return resEnd;
}/*}}}*/

void ScanfSSBondRecord(FILE* fp, SSBondPro* pSSBondPro, int numSSBond)/*{{{*/
{
    char   *pch;
    int     i, j;
    char    delim[]     = ",";
    int     maxline     = 200;
    char    line[200+1] = "";
    SSBond *pSSBond;

    pSSBondPro->numSSBondRes = 0;
    for(i = 0 ; i < numSSBond ; i++)
    {
        pSSBond = &(pSSBondPro->ssbond[i]);
        //line 1, 
        fgetline( fp, line, maxline);
        j = 0 ;
        pch = strtok (line,delim);
        while (pch != NULL)
        {
            if( j == 0) 
                sscanf(pch,"%d", &pSSBond->numRes);
            else
                ScanfResRecord3(pch, &pSSBond->res[j-1], RESSEQ );
            j++;
            pch = strtok (NULL, delim);
        }
        // line2 
        fgetline( fp, line, maxline);
        j = 0 ;
        pch = strtok (line,delim);
        while (pch != NULL)
        {
            if( j == 0) 
                sscanf(pch,"%d", &pSSBond->numRes);
            else
                ScanfResRecord3(pch, &pSSBond->res[j-1], AASEQINDEX );
            j++;
            pch = strtok (NULL, delim);
        }
        pSSBondPro->numSSBondRes += (pSSBond->numRes);
#ifdef _ASSERT_
        assert( pSSBond->numRes <= 2 && pSSBond->numRes >= 1 );
#endif
    }
}
/*}}}*/
void ScanfDSSPResRecord(const char* line, DSSP_Residue* pDSSPRes)/*{{{*/
{
    InitDSSPRes(pDSSPRes);
    int i , j;
    char str[100] = "";
    for(i=0 ;i<5;i++)	str[i] = line[i]; str[i] = '\0';
    sscanf(str,"%d",&pDSSPRes->seqResSer);

    for(i=0 ;i<4;i++)	str[i] = line[i+6]; str[i] = '\0';
    sscanf(str,"%d",&pDSSPRes->resSeq);

    pDSSPRes->iCode =		line[10];
    pDSSPRes->chainID =		line[11];
    pDSSPRes->aa =			line[13];
    pDSSPRes->chainBreak =	line[14];
    pDSSPRes->ss =			line[16];
    pDSSPRes->turn3 =		line[18];
    pDSSPRes->turn4 =		line[19];
    pDSSPRes->turn5 =		line[20];
    pDSSPRes->geoBend =		line[21];
    pDSSPRes->chira =		line[22];
    pDSSPRes->bBridge1 =	line[23];
    pDSSPRes->bBridge2 =	line[24];

    for(i=0 ;i<4;i++)	str[i] = line[i+25]; str[i] = '\0';
    sscanf(str,"%d",&pDSSPRes->bp1);

    for(i=0 ;i<4;i++)	str[i] = line[i+29]; str[i] = '\0';
    sscanf(str,"%d",&pDSSPRes->bp2);

    pDSSPRes->bSheet =	line[33];

    for(i=0 ;i<4;i++)	str[i] = line[i+34]; str[i] = '\0';
    sscanf(str,"%d",&pDSSPRes->acc);

    for(j = 0 ; j < 4  ; j ++)
    {
        for(i = 0 ; i < 11 ; i++) str[i] = line[j*11 + i + 39]; str[i] = '\0';
        sscanf(str,"%d,%lf",&pDSSPRes->hbond[j].pos,&pDSSPRes->hbond[j].e);
    }

    for(i=0 ;i<8;i++)	str[i] = line[i+83]; str[i] = '\0';
    sscanf(str,"%f",&pDSSPRes->tco);

    for(i=0 ;i<6;i++)	str[i] = line[i+91]; str[i] = '\0';
    sscanf(str,"%f",&pDSSPRes->kappa);

    for(i=0 ;i<6;i++)	str[i] = line[i+97]; str[i] = '\0';
    sscanf(str,"%f",&pDSSPRes->alpha);

    for(i=0 ;i<6;i++)	str[i] = line[i+103]; str[i] = '\0';
    sscanf(str,"%f",&pDSSPRes->phi);

    for(i=0 ;i<6;i++)	str[i] = line[i+109]; str[i] = '\0';
    sscanf(str,"%f",&pDSSPRes->psi);

    for(i=0 ;i<7;i++)	str[i] = line[i+115]; str[i] = '\0';
    sscanf(str,"%f",&pDSSPRes->x);

    for(i=0 ;i<7;i++)	str[i] = line[i+122]; str[i] = '\0';
    sscanf(str,"%f",&pDSSPRes->y);

    for(i=0 ;i<7;i++)	str[i] = line[i+129]; str[i] = '\0';
    sscanf(str,"%f",&pDSSPRes->z);
}/*}}}*/

int ScanfKernerlRecordID(const char *vectorRecordID, char *id, int &length, int &idx, char *res_1char_list, int *aaSeqIndex, int numSite)/*{{{*/

    /*****************************************************************************
     * Scanf the vectorRecordID
     * return the number of items scaned
     ****************************************************************************/
{   
    char *pch;
    int  status_sscanf = 0;
    char delim[] = WHITE_SPACE;
    int linesize = strlen(vectorRecordID);

    assert( linesize > 0 && linesize < 9999999);
    Array1D <char> str_1darray(linesize+10) ;
    char *str = str_1darray.array1D;
    my_strcpy(str,vectorRecordID, linesize);

    int i;
    ssubstitute(str,  CHAR_VECTOR_ID_SEPRATOR, ' ');
    pch = strtok (str,delim);
    i = 0;
    while (pch != NULL) // 0    1     2      3     4            5     6
    {                   // id length  idx   res1  aaSeqIndex1  res2  aaSeqIndex2
        if( i == 0) 
        {
            status_sscanf = 1;
            my_strcpy(id, pch, SIZE_CHAIN_ID);
        }
        else if( i == 1)
            status_sscanf = sscanf(pch,"%d", &length);
        else if( i == 2)
            status_sscanf = sscanf(pch,"%d", &idx);
        else if(i >= 3 && i < 3 + numSite*2)
        {
            if(i%2 == 1)
                status_sscanf = sscanf(pch,"%c",&res_1char_list[(i-3)/2]);
            else
            {
                status_sscanf = sscanf(pch,"%d", &aaSeqIndex[(i-3)/2]);
                aaSeqIndex[(i-3)/2] --;
            }
        }
        else
        { }
        assert( status_sscanf == 1);
        i++;
        pch = strtok (NULL, delim);
    }
    res_1char_list[numSite] = '\0';
    return i;
}
/*}}}*/
int  ScanfGistPredictCard(const char *line, char *id, int &length, int &idx, char *res_1char_list,int *aaSeqIndex, int &label, double &discriminant, int numSite)/*{{{*/
    /*****************************************************************************
     *  idx : indexing of vector vector for that chain
     ****************************************************************************/
{
    int  status_sscanf;
    int linesize = strlen(line);
    assert( linesize > 0 && linesize < 9999999);

    Array1D <char> vectorRecordID_1darray(linesize) ;
    char *vectorRecordID = vectorRecordID_1darray.array1D;

    status_sscanf = sscanf(line,"%s %d %lf",vectorRecordID, &label, &discriminant);
    assert(status_sscanf == 3);
    int cnt = ScanfKernerlRecordID(vectorRecordID, id, length, idx, res_1char_list,aaSeqIndex, numSite);

//    StdID(id);
    return (2 + cnt);
}
/*}}}*/
int  ScanfGistVectorLabelCard(const char *line, char *id, int &length, int &idx, char *res_1char_list,int *aaSeqIndex, int &label,  int numSite)/*{{{*/
{
    int  status_sscanf;
    int linesize = strlen(line);
    assert( linesize > 0 && linesize < 9999999);

    Array1D <char> vectorRecordID_1darray(linesize) ;
    char *vectorRecordID = vectorRecordID_1darray.array1D;

    status_sscanf = sscanf(line,"%s %d",vectorRecordID, &label);
    assert(status_sscanf == 2);
    int cnt = ScanfKernerlRecordID(vectorRecordID, id, length, idx, res_1char_list,aaSeqIndex, numSite);

    StdID(id);
    return (1 + cnt);
}
/*}}}*/
void ScanfResRecord3(const char *str, Residue *pRes, int tag)/*{{{*/
{
    int status_sscanf;
    if(tag == RESSEQ)
    {
        if((status_sscanf = sscanf(str,"%3s%d%c(%c)", pRes->resName, &pRes->resSeq, &pRes->resICode, &pRes->chainID)) < 4)
        {
            fprintf(stderr,"Res record=%s\n",str);
            assert( status_sscanf == 4 );
        }
        pRes->aa = AA3To1(pRes->resName);
    }
    else if(tag == AASEQINDEX)
    {
        char tmpc;
        if((status_sscanf = sscanf(str,"%3s%d%c(%c)", pRes->resName, &pRes->aaSeqIndex, &tmpc, &pRes->chainID)) < 4)
        {
            printf("Res record=%s\n", str);
            assert( status_sscanf == 4);
        }
        pRes->aaSeqIndex --; //because here when read in, aaSeqIndex starts from 1,while aaSeqIndex should start from 0
    }
    else
    { }
}
/*}}}*/
int ScanfCloseMetalRes(FILE* fpin, Residue* res,int numRes)/*{{{*/
    //*******************************************************************
    //ScanfCloseMetalRes()
    //*******************************************************************
    // read in the residues close to the metal atom in the following format
    // 		HIS37 (A) ,  HIS51 (A)    "--residue name and resSeq in pdb file
    // 		HIS37 (A) ,  HIS51 (A)    "--residue name and aaSeqIndex in sequence
    // 		NE2       ,  NE2          "--atom name 
    // 		2.004     ,  2.010        "--distance from atom to metal atom
    // 		A         ,  V             "--shape string simbol
    // 		100       ,  100           "--conservation level of the residue generated	psi-blast
    //*******************************************************************
{
    int  j;
    char delim[] = ",";
    int maxline = 200;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;

    char *pch;
    //scanf line1: list of "Res resSeq resIcode chainID"
    fgetline( fpin, line, maxline);
    strtrim(line);
    j = 0 ;
    pch = strtok (line,delim);
    while (pch != NULL)
    {
        ScanfResRecord3(pch, &res[j], RESSEQ );
        j++;
        pch = strtok (NULL, delim);
    }
    // line2   list of "Res aaSeqIndex ' ' chainID" 
    fgetline( fpin, line, maxline);
    strtrim(line);
    j = 0 ;
    pch = strtok (line,delim);
    while (pch != NULL)
    {
        ScanfResRecord3(pch, &res[j], AASEQINDEX );
        j++;
        pch = strtok (NULL, delim);
    }

    //line3 atoms bind to metal atoms, 2007-04-19
    fgetline( fpin, line, maxline); //skip this line, 2007-04-19, for future implementation


    //line4, distance from atoms to the metal atom, 2007-04-19
    fgetline( fpin, line, maxline);  //skip this line, 2007-04-19  for future implementation 


    //scanf line5: list of shape string
    fgetline( fpin, line, maxline);
    strtrim(line);
    j = 0 ;
    pch = strtok (line,delim);
    while (pch != NULL)
    {
        sscanf(pch, " %c ", &(res[j].shape)); /*bug fixed, 2010-04-12, pch added in front of " %c "*/
        j++;
        pch = strtok (NULL, delim);
    }
    //scanf line6: list of conservation value
    fgetline( fpin, line, maxline);
    strtrim(line);
    j = 0 ;
    pch = strtok (line,delim);
    while (pch != NULL)
    {
        res[j].consv = atoi(pch);
        j++;
        pch = strtok (NULL, delim);
    }

    return 0;
}/*}}}*/

void ScanfSCOPRecord(const char *line, SCOP *pSCOP)/*{{{*/
{
    char domstr[SIZE_DOMAIN_DEF_RECORD+1] = "";
    char clsstr[50] = "";
    sscanf(line,"%s %s %s %s %d cl=%d,cf=%d,sf=%d,fa=%d,dm=%d,sp=%d,px=%d",
            pSCOP->did,
            pSCOP->pdbid,
            domstr,
            clsstr,
            &pSCOP->idnum,
            &pSCOP->cl,
            &pSCOP->cf,
            &pSCOP->sf,
            &pSCOP->fa,
            &pSCOP->dm,
            &pSCOP->sp,
            &pSCOP->px);

    pSCOP->chainID = pSCOP->did[5];
    my_strupr(&pSCOP->chainID);
    sscanf(clsstr,"%c.%d.%d.%d", &pSCOP->cls, &pSCOP->fod, &pSCOP->sup, &pSCOP->fam);
    my_strcpy(pSCOP->domainDef, domstr, SIZE_DOMAIN_DEF_RECORD);

    char* pch;
    char delim[] = ",";
    int i = 0 ;
    /* there are three types of chainid in scop id, 
     * 1. '_' e.g. d1hcz_1, d1hcv__
     * 2. 0-9, a-z, e.g.  d1hd1a_ 
     * 3. '.'  e.g. d1hdt.1, which is a multi polypeptide domain
     * */

    i = 0;  /* i, iterator for pSCOP->domDef */
    pch = strtok( domstr, delim );
    while( pch != NULL )
    {
        if(pSCOP->did[5] == '_')
        {
            pSCOP->domDef.chainIDs[i] = ' ';
            if(strcmp(pch, "-") == 0)
                //d1lci__ 1lci    -   e.23.1.1    43349  
            {
                pSCOP->domDef.isWholeChain[i] = true;
                pSCOP->domDef.posF[i] = 0 ; // if it's a whole chain, initialize the position
                pSCOP->domDef.posT[i] = 0x7FFFFFFF;
                pSCOP->domDef.icodeF[i] = ' ';
                pSCOP->domDef.icodeT[i] = ' ';
            }
            else /* defined with range,*/
                //d1ldg_1   1ldg    18-163  c.2.1.5 30165  
                //d1ak2_1  1ak2  14-146,177-233
            {
                pSCOP->domDef.isWholeChain[i] = false;
                //sscanf(pch, "%d-%d", &(pSCOP->domDef.posF[i]), &(pSCOP->domDef.posT[i]));
                // 2007-07-12, there are icode in SCOP range definition, should
                // be inportant, e.g.
                // d1pca_2	1pca	4A-99A	d.58.3.1	39063	cl=53931,cf=54861,sf=54897,fa=54898,dm=54899,sp=54900,px=39063
                // d1jqga2  1jqg    A:4P-100P   d.58.3.1    71792 cl=53931,cf=54861,sf=54897,fa=54898,dm=54899,sp=75429,px=71792
                
                int size_str = strlen(pch);
                char *rangestr = new char [size_str+1];
                my_strcpy(rangestr, pch , size_str);
                ssubstitute(rangestr, '-', ' '); // change "4A-99A" to "4A 99A"
                char *str1 = new char [size_str +1];
                char *str2 = new char [size_str +1];
                sscanf(rangestr, "%s %s", str1,str2);
                int len1 = strlen(str1); 
                int len2 = strlen(str2);
                if( !isdigit (str1[len1-1])) {
                    pSCOP->domDef.icodeF[i] = str1[len1-1]; 
                    str1[len1-1] = '\0';
                } else { 
                    pSCOP->domDef.icodeF[i] =  ' '; 
                }

                if( !isdigit (str2[len2-1])) { 
                    pSCOP->domDef.icodeT[i] = str2[len2-1]; 
                    str2[len2-1] = '\0';
                } else { 
                    pSCOP->domDef.icodeT[i] =  ' '; 
                }

                pSCOP->domDef.posF[i] = atoi(str1);
                pSCOP->domDef.posT[i] = atoi(str2);
            
                delete [] rangestr;
                delete [] str1;
                delete [] str2;
            }
        } else  /*chainid is not null*/ {
            int size_str = strlen(pch);
            char *rangestr = new char [size_str+1];
            my_strcpy(rangestr, pch , size_str);
            ssubstitute(rangestr, ':', ' '); // change "A:4A-99A" to "A 4A-99A"
            ssubstitute(rangestr, '-', ' '); // change "A 4A-99A" to "A 4A 99A"
            char *str1 = new char [size_str +1];
            char *str2 = new char [size_str +1];
            char *str3 = new char [size_str +1];
            int numitem = sscanf(rangestr, "%s %s %s", str3,str1, str2);
            if(numitem < 3)
            {
                pSCOP->domDef.isWholeChain[i] = true;
                pSCOP->domDef.posF[i] = 0 ; // if it's a whole chain, initialize the position
                pSCOP->domDef.posT[i] = 0x7FFFFFFF;
                pSCOP->domDef.icodeF[i] = ' ';
                pSCOP->domDef.icodeT[i] = ' ';
                pSCOP->domDef.chainIDs[i] = str3[0];
            } else {
                pSCOP->domDef.isWholeChain[i] = false;
                int len1 = strlen(str1); 
                int len2 = strlen(str2);
                /*int len3 = strlen(str3);*/
                pSCOP->domDef.chainIDs[i] = str3[0];

                if( !isdigit (str1[len1-1]))
                {
                    pSCOP->domDef.icodeF[i] = str1[len1-1]; 
                    str1[len1-1] = '\0';
                } else { pSCOP->domDef.icodeF[i] =  ' '; }

                if( !isdigit (str2[len2-1])) { 
                    pSCOP->domDef.icodeT[i] = str2[len2-1]; 
                    str2[len2-1] = '\0';
                } else { pSCOP->domDef.icodeT[i] =  ' '; }

                pSCOP->domDef.posF[i] = atoi(str1);
                pSCOP->domDef.posT[i] = atoi(str2);

            }
            delete [] rangestr;
            delete [] str1;
            delete [] str2;
            delete [] str3;
        }
        //debug code
        if(pSCOP->did[6] == '_' && !pSCOP->domDef.isWholeChain[i])
        {
            fprintf(stderr,"did=%s  domainDef=%s",pSCOP->did, pSCOP->domainDef);
            assert(pSCOP->did[6] == '_' && pSCOP->domDef.isWholeChain[i]);
        }

        pch = strtok( NULL, delim );
        i ++;
    }
    pSCOP->domDef.numChain = i ;
    pSCOP->domDef.chainIDs[pSCOP->domDef.numChain] = '\0';
}
/*}}}*/

/*structure implementation, init, copy, delete, allocate*/
void InitChain(Chain * pChain)/*{{{*/
{
    pChain->seqtype  = UNKNOWN_SEQ_TYPE;
    pChain->numRes   = 0;
    pChain->chainID  = '\0'; // initialize chainID as '\0'
    pChain->title    = NULL; 
    pChain->aaSeq    = NULL;
    pChain->shString = NULL;
    pChain->secStruc = NULL;
    pChain->resSer   = NULL;
    pChain->resICode = NULL;
    pChain->consv    = NULL;
    pChain->waterAcc = NULL;
    pChain->phi      = NULL;
    pChain->psi      = NULL;
}/*}}}*/
void CopyChain(Chain * to, Chain * from)/*{{{*/
/*****************************************************************************
 * copy the chain "from" to "to", 
 * copy only when both from and to has allocated memory
 * mmeory overflow may happen
 ****************************************************************************/
{
    int i = 0;
    int numRes = from->numRes;
    to->seqtype  = from->seqtype;
    to->numRes   = from->numRes;
    to->chainID  = from->chainID; 
    if( from -> title != NULL && to -> title != NULL)
    {
        my_strcpy(from->title, to->title, strlen(to->title));
    }
    if( from -> aaSeq != NULL && to -> aaSeq != NULL)
    {
        for( i = 0 ; i < numRes ; i ++)
            to -> aaSeq[i] = from -> aaSeq[i];
        to -> aaSeq[i] = '\0';
    }
    if( from -> shString != NULL && to -> shString != NULL)
    {
        for( i = 0 ; i < numRes ; i ++)
            to -> shString[i] = from -> shString[i];
        to -> shString[i] = '\0';
    }
    if( from -> secStruc != NULL && to -> secStruc != NULL)
    {
        for( i = 0 ; i < numRes ; i ++)
            to -> secStruc[i] = from -> secStruc[i];
        to -> secStruc[i] = '\0';
    }
    if( from -> resSer != NULL && to -> resSer != NULL)
    {
        for( i = 0 ; i < numRes ; i ++)
            to -> resSer[i] = from -> resSer[i];
    }
    if( from -> resICode != NULL && to -> resICode != NULL)
    {
        for( i = 0 ; i < numRes ; i ++)
            to -> resICode[i] = from -> resICode[i];
    }
    if( from -> consv != NULL && to -> consv != NULL)
    {
        for( i = 0 ; i < numRes ; i ++)
            to -> consv[i] = from -> consv[i];
    }
    if( from -> waterAcc != NULL && to -> waterAcc != NULL)
    {
        for( i = 0 ; i < numRes ; i ++)
            to -> waterAcc[i] = from -> waterAcc[i];
    }
    if( from -> phi != NULL && to -> phi != NULL)
    {
        for( i = 0 ; i < numRes ; i ++)
            to -> phi[i] = from -> phi[i];
    }
    if( from -> psi != NULL && to -> psi != NULL)
    {
        for( i = 0 ; i < numRes ; i ++)
            to -> psi[i] = from -> psi[i];
    }
}/*}}}*/
void DeleteChain(Chain *pChain)/*{{{*/
{
    if(pChain->title    != NULL) { delete [] pChain->title ; }
    if(pChain->aaSeq    != NULL) { delete [] pChain->aaSeq ; }
    if(pChain->resICode != NULL) { delete [] pChain->resICode; }
    if(pChain->resSer   != NULL) { delete [] pChain->resSer; }
    if(pChain->shString != NULL) { delete [] pChain->shString; }
    if(pChain->secStruc != NULL) { delete [] pChain->secStruc; }
    if(pChain->consv    != NULL) { delete [] pChain->consv; }
    if(pChain->waterAcc != NULL) { delete [] pChain->waterAcc; }
    if(pChain->phi      != NULL) { delete [] pChain->phi; }
    if(pChain->psi      != NULL) { delete [] pChain->psi; }
    InitChain(pChain);
}/*}}}*/

void InitGistPredChain(GistPredChain * pChain)/*{{{*/
{
    pChain->numRes       = 0;
    pChain->length       = 0;
    pChain->aaSeq        = NULL;
    pChain->aaSeqIndex   = NULL;
    pChain->label        = NULL;
    pChain->discriminant = NULL;
    my_strcpy(pChain->id, "", SIZE_CHAIN_ID);
}/*}}}*/
void AllocGistPredChain(GistPredChain * pChain, int size)/*{{{*/
{
    pChain->aaSeq        = new char[size+1];
    pChain->aaSeqIndex   = new int[size];
    pChain->label        = new int[size];
    pChain->discriminant = new double[size];
}/*}}}*/
void CopyGistPredChain(GistPredChain *to, GistPredChain *from)/*{{{*/
{
    to -> length = from -> length;
    to -> numRes = from -> numRes;
    my_strcpy(to -> id, from -> id, SIZE_CHAIN_ID);

    int numRes = from -> numRes;
    int i;
    if( from -> aaSeqIndex != NULL && to -> aaSeqIndex != NULL)
    {
        for( i = 0 ; i < numRes ; i ++)
            to -> aaSeqIndex[i] = from -> aaSeqIndex[i];
    }
    if( from -> aaSeq != NULL && to -> aaSeq != NULL)
    {
        for( i = 0 ; i < numRes ; i ++)
            to -> aaSeq[i] = from -> aaSeq[i];
        to -> aaSeq[i] = '\0';
    }
    if( from -> label != NULL && to -> label != NULL)
    {
        for( i = 0 ; i < numRes ; i ++)
            to -> label[i] = from -> label[i];
    }
    if( from -> discriminant != NULL && to -> discriminant != NULL)
    {
        for( i = 0 ; i < numRes ; i ++)
            to -> discriminant[i] = from -> discriminant[i];
    }
}
/*}}}*/
void DeleteGistPredChain(GistPredChain *pChain)/*{{{*/
{
    if(pChain->aaSeq != NULL)
    {
        delete [] pChain->aaSeq ;
    }
    if(pChain->aaSeqIndex != NULL)
    {
        delete [] pChain->aaSeqIndex;
    }
    if(pChain->label != NULL)
    {
        delete [] pChain->label;
    }
    if(pChain->discriminant != NULL)
    {
        delete [] pChain->discriminant;
    }
    InitGistPredChain(pChain);
}/*}}}*/

void InitDomainDEF(DomainDEF *pDomainDEF)/*{{{*/
{
    pDomainDEF->numChain = 0;
    strcpy (pDomainDEF->chainIDs,"" );
    for(int i = 0 ; i < NUM_CHAIN_PER_DOMAIN; i++)
    {
        pDomainDEF->isWholeChain[i] = false;
        pDomainDEF->posF[i] = 0;
        pDomainDEF->posT[i] = 0;
        pDomainDEF->icodeF[i] = ' ';
        pDomainDEF->icodeT[i] = ' ';
    }
}
/*}}}*/
void InitSCOP(SCOP *pSCOP)/*{{{*/
{
    strcpy (pSCOP->did, "" );
    strcpy (pSCOP->pdbid, "" );
    pSCOP->chainID = '\0';
    strcpy (pSCOP->domainDef, "" );
    pSCOP->idnum = 0;
    pSCOP->cl = 0;
    pSCOP->cf = 0;
    pSCOP->sf = 0;
    pSCOP->fa = 0;
    pSCOP->dm = 0;
    pSCOP->sp = 0;
    pSCOP->px = 0;
    pSCOP->cls = ' ';
    pSCOP->fod = 0;
    pSCOP->sup = 0;
    pSCOP->fam = 0;
    InitDomainDEF(& (pSCOP->domDef));
}
/*}}}*/
void CopySCOP(SCOP* to, SCOP* from)/*{{{*/
{
	my_strcpy(to->did,from->did, SIZE_SCOP_ID);
	my_strcpy(to->pdbid,from->pdbid, SIZE_PDBID);
	my_strcpy(to->domainDef,from->domainDef, SIZE_DOMAIN_DEF_RECORD);
	to->idnum = from->idnum;
	to->chainID = from->chainID;
	to->cl = from->cl;
	to->cf = from->cf;
	to->sf = from->sf;
	to->fa = from->fa;
	to->dm = from->dm;
	to->sp = from->sp;
	to->px = from->px;

	to->cls = from->cls;
	to->fod = from->fod;
	to->sup = from->sup;
	to->fam = from->fam;

	to->domDef.numChain = from->domDef.numChain;

	my_strcpy(to->domDef.chainIDs,from->domDef.chainIDs, NUM_CHAIN_PER_DOMAIN);
	for(int i = 0 ; i < from->domDef.numChain; i++)
	{
		to->domDef.posF[i] = from->domDef.posF[i];
		to->domDef.posT[i] = from->domDef.posT[i];
		to->domDef.icodeF[i] = from->domDef.icodeF[i];
		to->domDef.icodeT[i] = from->domDef.icodeT[i];
		to->domDef.isWholeChain[i] = from->domDef.isWholeChain[i];
	}
}
/*}}}*/

void InitMODM(MODM *pMODM)/*{{{*/
{
    strcpy(pMODM -> id , "" );
    pMODM -> length       = 0;
    pMODM -> type_modm    = MODM_LOG;
    pMODM -> M            = NULL;
    pMODM -> log_M        = NULL;
    pMODM -> score1       = NULL;
    pMODM -> score2       = NULL;
    pMODM -> consv        = NULL;
    pMODM -> aaSeq        = NULL;
    pMODM -> alphabetMODM = NULL;
    pMODM -> waterAcc     = NULL;
    pMODM -> shString     = NULL;
    pMODM -> dsspSec      = NULL;

}/*}}}*/
void AllocMODM(MODM *pMODM, int length, bool isAllocLogM /*=false*/, int sizeAlphabet /*= NUM_BLOSUM*/)/*{{{*/
{
    pMODM -> M       = Create2DArray(pMODM->M,length,sizeAlphabet);
    if(isAllocLogM) 
        pMODM-> log_M = Create2DArray(pMODM->log_M, length, sizeAlphabet);
    pMODM -> score1  = new double[length];
    pMODM -> score2  = new double[length];
    pMODM -> consv   = new double[length];
    pMODM -> aaSeq   = new char[length+1];
    pMODM -> alphabetMODM = new char[sizeAlphabet+1];
}/*}}}*/
void CopyMODM(MODM *to, MODM *from)/*{{{*/
{
    to -> length = from -> length;
    to -> type_modm = from -> type_modm;
    int length = from -> length;
    int sizeAlphabet = strlen(from->alphabetMODM);
    int i,j;

    my_strcpy(to -> id, from -> id, SIZE_CHAIN_ID);
    if( from -> alphabetMODM != NULL && to -> alphabetMODM != NULL)
        my_strcpy(to->alphabetMODM, from->alphabetMODM, NUM_BLOSUM);
    else 
        assert ( from -> alphabetMODM != NULL && to -> alphabetMODM != NULL);

    if( from -> M != NULL && to -> M != NULL)
    {
        for ( i = 0 ; i < length ; i ++)
            for(j = 0 ; j < sizeAlphabet; j ++)
                to->M[i][j] = from->M[i][j];  
    }
    else
    {
        assert( !(from -> M != NULL && to -> M == NULL));// if from-M exist while to->M is NULL, assert
    }

    if( from -> log_M != NULL && to -> log_M != NULL)
    {
        for ( i = 0 ; i < length ; i ++)
            for(j = 0 ; j < sizeAlphabet; j ++)
                to->log_M[i][j] = from->log_M[i][j];  
    }
    else
    {
        assert( !(from -> log_M != NULL && to -> log_M == NULL));
    }

    if( from -> score1 != NULL && to -> score1 != NULL)
    {
        for ( i = 0 ; i < length ; i ++)
            to->score1[i] = from->score1[i];  
    }
    else
    {
        assert(!(from -> score1 != NULL && to -> score1 == NULL));
    }

    if( from -> score2 != NULL && to -> score2 != NULL)
    {
        for ( i = 0 ; i < length ; i ++)
            to->score2[i] = from->score2[i];  
    }
    else
    {
        assert(!(from -> score2 != NULL && to -> score2 == NULL));
    }

    if( from -> consv != NULL && to -> consv != NULL)
    {
        for ( i = 0 ; i < length ; i ++)
            to->consv[i] = from->consv[i];  
    }
    else
    {
        assert(!(from -> consv != NULL && to -> consv == NULL));
    }

    if( from -> aaSeq != NULL && to -> aaSeq != NULL)
    {
        my_strcpy(to->aaSeq, from->aaSeq, LONGEST_SEQ);
    }
    else
    {
        assert(! (from -> aaSeq != NULL && to -> aaSeq == NULL ));
    }

    if( from -> shString != NULL && to -> shString != NULL)
    {
        my_strcpy(to->shString, from->shString, LONGEST_SEQ);
    }
    if( from -> dsspSec != NULL && to -> dsspSec != NULL)
    {
        my_strcpy(to->dsspSec, from->dsspSec, LONGEST_SEQ);
    }
    if( from -> waterAcc != NULL && to -> waterAcc != NULL)
    {
        for ( i = 0 ; i < length ; i ++)
            to->waterAcc[i] = from->waterAcc[i];  
    }

}
/*}}}*/
void DeleteMODM(MODM *pMODM, int length)/*{{{*/
    // length is the x size of  M matrix, to avoid memory leak,
    // setting of length is compulsory
{
    if(pMODM->M != NULL)
    {
        Delete2DArray(pMODM->M, length);
    }
    if(pMODM->log_M != NULL)
    {
        Delete2DArray(pMODM->log_M, length);
    }
    if(pMODM->score1 != NULL)
    {
        delete [] pMODM->score1 ;
    }
    if(pMODM->score2 != NULL)
    {
        delete [] pMODM->score2 ;
    }
    if(pMODM->consv != NULL)
    {
        delete [] pMODM->consv ;
    }
    if(pMODM->aaSeq != NULL)
    {
        delete [] pMODM->aaSeq ;
    }
    if(pMODM->alphabetMODM != NULL)
    {
        delete [] pMODM->alphabetMODM ;
    }
    if(pMODM->waterAcc != NULL)
    {
        delete [] pMODM->waterAcc ;
    }
    if(pMODM->shString != NULL)
    {
        delete [] pMODM->shString ;
    }
    if(pMODM->dsspSec != NULL)
    {
        delete [] pMODM->dsspSec ;
    }
    InitMODM(pMODM);
}/*}}}*/


void InitResidue(Residue *pRes)/*{{{*/
{
    pRes -> resSeq     = 0;
    pRes -> aaSeqIndex = 0;
    pRes -> aa         = ' ';
    pRes -> chainID    = ' ';
    pRes -> resICode   = ' ';
    pRes -> shape      = INIT_SHAPE;
    pRes -> consv      = INIT_CONSV;
    pRes -> isBioMetalBound = false;
    pRes -> numAtom    = 0;
    pRes -> atom       = NULL;
    pRes -> numMetalBound = 0;
    pRes -> parentAtomEnvIndex = NULL;
    strcpy( pRes -> resName, "");
}/*}}}*/
void CopyResidue(Residue* to,Residue* from)/*{{{*/
{
    to->resSeq           = from->resSeq;
    to->resICode         = from->resICode;
    to->aaSeqIndex       = from->aaSeqIndex;
    to->chainID          = from->chainID;
    to->shape            = from->shape;
    to->consv            = from->consv;
    to->aa               = from->aa;
    to->isBioMetalBound  = from->isBioMetalBound;
    to->numMetalBound    = from->numMetalBound;
    if(to->parentAtomEnvIndex != NULL && from->parentAtomEnvIndex != NULL)
    {
        for(int i = 0 ; i < from -> numMetalBound; i++)
            to->parentAtomEnvIndex[i] = from->parentAtomEnvIndex[i];
    }
    if(to->atom != NULL && from -> atom != NULL)
    {
        for(int i = 0 ; i < from -> numAtom; i ++)
            CopyAtom (&(to->atom[i]), &(from->atom[i]));
    }
    my_strcpy(to->resName      , from->resName      , SIZE_RES_NAME) ;
}/*}}}*/
void DeleteResidue(Residue *pRes, int numAtom /*= 0 */)/*{{{*/
{
    if( pRes -> atom != NULL)
    {
        if(numAtom == 0) numAtom = pRes->numAtom;
        for(int i = 0 ; i < numAtom; i ++)
            DeleteAtom(&(pRes->atom[i]));
        delete [] pRes -> atom;
    }

    if( pRes -> parentAtomEnvIndex != NULL)
        delete [] pRes->parentAtomEnvIndex;
    InitResidue(pRes);
}
/*}}}*/

void InitAtom(Atom * pAtom)/*{{{*/
{
    pAtom->serial     = INIT_INT;
    pAtom->chainID    = ' ';
    pAtom->altLoc     = ' ';
    pAtom->resSeq     = INIT_INT;
    pAtom->iCode      = ' ';
    pAtom->occupancy  = INIT_FLOAT;
    pAtom->tempFactor = INIT_FLOAT;
    pAtom->x          = INIT_DOUBLE;
    pAtom->y          = INIT_DOUBLE;
    pAtom->z          = INIT_DOUBLE;
    strcpy(pAtom->recordID, "");
    strcpy(pAtom->name    , "");
    strcpy(pAtom->resName , "");
    strcpy(pAtom->origName, "");
    strcpy(pAtom->segID   , "");
    strcpy(pAtom->element , "");
    strcpy(pAtom->charge  , "");

    pAtom->dist = NULL;
    pAtom->numDist = 0;
}/*}}}*/
void CopyAtom(Atom* to, Atom* from)/*{{{*/
{
    //InitAtom(to);
    to->serial     = from->serial;
    to->chainID    = from->chainID;
    to->altLoc     = from->altLoc;
    to->resSeq     = from->resSeq;
    to->iCode      = from->iCode;
    to->occupancy  = from->occupancy;
    to->tempFactor = from->tempFactor;
    to->x          = from->x;
    to->y          = from->y;
    to->z          = from->z;
    my_strcpy(to->recordID, from->recordID, SIZE_RECORD_ID);    
    my_strcpy(to->name    , from->name    , SIZE_ATOM_NAME);    
    my_strcpy(to->resName , from->resName , SIZE_RES_NAME);
    my_strcpy(to->origName, from->origName, SIZE_ATOM_ORIGNAME);
    my_strcpy(to->segID   , from->segID   , SIZE_ATOM_SEGID);
    my_strcpy(to->element , from->element , SIZE_ATOM_ELEMENT);
    my_strcpy(to->charge  , from->charge  , SIZE_ATOM_CHARGE);

    to->numDist = from->numDist;
    if(to -> dist != NULL && from -> dist != NULL)
    {
        for(int i = 0 ; i < from -> numDist; i++)
            to -> dist[i] = from -> dist[i];
    }

}/*}}}*/
void DeleteAtom(Atom *pAtom)/*{{{*/
{
    if(pAtom->dist != NULL)
        delete [] pAtom->dist;
    InitAtom(pAtom);
}/*}}}*/

void InitDSSPRes(DSSP_Residue* dssp_res)/*{{{*/
{
    int i;
    dssp_res->seqResSer = INIT_INT;	
    dssp_res->resSeq = INIT_INT;	
    dssp_res->iCode = ' ';		
    dssp_res->chainID = ' ';
    dssp_res->aa = ' ';		
    dssp_res->chainBreak = ' ';
    dssp_res->ss = ' ';		
    dssp_res->turn3 = ' ';		
    dssp_res->turn4 = ' ';		
    dssp_res->turn5 = ' ';		
    dssp_res->geoBend = ' ';	
    dssp_res->chira = ' ';		
    dssp_res->bBridge1 = ' ';	
    dssp_res->bBridge2 = ' ';	
    dssp_res->bp1 = INIT_INT;		
    dssp_res->bp2 = INIT_INT;		
    dssp_res->bSheet = ' ';	
    dssp_res->acc = INIT_INT;
    for( i = 0 ; i < 4 ; i++)
    {	
        dssp_res->hbond[i].pos = INIT_INT;
        dssp_res->hbond[i].e = INIT_DOUBLE;	
    }
    dssp_res->tco = INIT_DOUBLE ;		
    dssp_res->kappa = INIT_DOUBLE ;	
    dssp_res->alpha = INIT_DOUBLE ;	
    dssp_res->phi = INIT_DOUBLE ;		
    dssp_res->psi = INIT_DOUBLE ;		
    dssp_res->x = INIT_DOUBLE ;		
    dssp_res->y = INIT_DOUBLE ;		
    dssp_res->z = INIT_DOUBLE ;
}/*}}}*/

void InitAtomEnv(AtomEnv *pAtomEnv)/*{{{*/
{
    pAtomEnv->numRes = 0;
    pAtomEnv->totalBoundRes = 0;
    pAtomEnv->metalAtomResSeq = INIT_RESSEQ;
    pAtomEnv->metalAtomChainID = ' ';
    pAtomEnv->seqLength= 0;
    pAtomEnv->res= NULL;
    pAtomEnv->parentResIndex = NULL;
    strcpy(pAtomEnv->metalAtomName   , "");
    strcpy(pAtomEnv->metalAtomResName, "");
    strcpy(pAtomEnv->id              , "");
    strcpy(pAtomEnv->metalAtomPDBID  , "");
}
/*}}}*/
void CopyAtomEnv(AtomEnv *pAtomEnv1, AtomEnv *pAtomEnv2, bool isCopyResidue /*= true*/)/*{{{*/
    // copy 2 --> 1
{
    int i;
    int numRes;
    pAtomEnv1->numRes           = pAtomEnv2->numRes;
    pAtomEnv1->totalBoundRes    = pAtomEnv2->totalBoundRes;
    pAtomEnv1->metalAtomResSeq  = pAtomEnv2->metalAtomResSeq;
    pAtomEnv1->seqLength        = pAtomEnv2->seqLength;
    pAtomEnv1->metalAtomChainID = pAtomEnv2->metalAtomChainID;
    my_strcpy(pAtomEnv1->metalAtomName,pAtomEnv2->metalAtomName, SIZE_METAL_ATOM_NAME);
    my_strcpy(pAtomEnv1->metalAtomResName,pAtomEnv2->metalAtomResName, SIZE_METAL_ATOM_RES_NAME);
    my_strcpy(pAtomEnv1->id,pAtomEnv2->id, SIZE_CHAIN_ID);
    numRes = pAtomEnv2->numRes;
    if(isCopyResidue && pAtomEnv1->res != NULL && pAtomEnv2->res != NULL)
    {
        for(i = 0; i < numRes ; i++)
            CopyResidue(&(pAtomEnv1->res[i]),&(pAtomEnv2->res[i]));
    }
    if(pAtomEnv1->parentResIndex != NULL && pAtomEnv2->parentResIndex != NULL)
    {
        for( i = 0 ; i < numRes ; i++)
            pAtomEnv1->parentResIndex[i] = pAtomEnv2->parentResIndex[i];
    }
}
/*}}}*/
void DeleteAtomEnv(AtomEnv *pAtomEnv, int numRes  /*= 0*/, int numAtom_Res /*= 0*/ )/*{{{*/
    //if the memory of res.metalAtomList is allocated by a constant, it can be
    //freed by specifying numMetalBound_Res in the second argument explicitly
{
    if(pAtomEnv->res != NULL)
    {
        if(numRes == 0) numRes = pAtomEnv->numRes;
        for(int i = 0 ; i < numRes; i++)
            DeleteResidue(&(pAtomEnv->res[i]), numAtom_Res);
        delete [] pAtomEnv->res;
    }
    if(pAtomEnv->parentResIndex != NULL)
    {
        delete [] pAtomEnv->parentResIndex;
    }
    InitAtomEnv(pAtomEnv);
}
/*}}}*/

void InitSSBond(SSBond *pSSBond )/*{{{*/
{
    pSSBond->isInterChain = false;
    pSSBond->numRes       = 0;
    for(int i = 0 ; i < 2 ; i++)
        InitResidue(&(pSSBond->res[i]));
}
/*}}}*/
void CopySSBond(SSBond *to, SSBond *from, bool isCopyResidue /*= true*/)/*{{{*/
{
    int i;
    to->isInterChain = from->isInterChain;
    to->numRes       = from->numRes;
    if(isCopyResidue)
    {
        for(i = 0 ; i < from->numRes; i++)
            CopyResidue(&(to->res[i]),&(from->res[i]));
    }
}
/*}}}*/

void InitSSBondPro(SSBondPro *pSSBondPro )/*{{{*/
{
    pSSBondPro->numSSBond    = 0;
    pSSBondPro->numSSBondRes = 0;
    pSSBondPro->length       = 0;
    pSSBondPro->ssbond        = NULL;
    strcpy(pSSBondPro->id,"");
}
/*}}}*/
void CopySSBondPro(SSBondPro *to, SSBondPro *from, bool isCopySSBond /*= true*/)/*{{{*/
{
    int i;
    int numSSBond;
    to->length = from->length;
    to->numSSBond = from->numSSBond;
    to->numSSBondRes = from-> numSSBondRes;
    my_strcpy(to->id , from->id, SIZE_CHAIN_ID);
    numSSBond = from->numSSBond;
    if(isCopySSBond)
    {
        for(i = 0 ; i < numSSBond; i++)
            CopySSBond(&(to->ssbond[i]),&(from->ssbond[i]));
    }
}
/*}}}*/
void DeleteSSBondPro(SSBondPro *pSSBondPro)/*{{{*/
{
    if(pSSBondPro->ssbond != NULL)
        delete [] pSSBondPro->ssbond;
    InitSSBondPro(pSSBondPro);
}
/*}}}*/


void InitMetalPro(MetalPro *pMetalPro )/*{{{*/
{
    strcpy(pMetalPro->id,"");
    pMetalPro->numBoundRes    = 0;
    pMetalPro->length        = 0;
    pMetalPro->numMetalAtom  = 0;
    pMetalPro->metalAtomList = NULL;
    pMetalPro->res           = NULL;
    pMetalPro->resSeqIndex   = NULL;
    pMetalPro->atomEnv       = NULL;
}
/*}}}*/
void InitMetalPro(MetalPro2 *pMetalPro)/*{{{*/
{
    pMetalPro->numMetalAtom = 0;
    pMetalPro->length   = 0;
    pMetalPro->atomEnv  = NULL;
    strcpy(pMetalPro->id,"");
}
/*}}}*/
void CopyMetalPro(MetalPro *pMetalPro1, MetalPro *pMetalPro2, bool isCopyResidue/* = true*/)/*{{{*/
{
    int i;
    int numBoundRes;
    pMetalPro1->numBoundRes   = pMetalPro2->numBoundRes;
    pMetalPro1->length       = pMetalPro2->length;
    pMetalPro1->numMetalAtom = pMetalPro2->numMetalAtom;
    my_strcpy(pMetalPro1->id , pMetalPro2->id, SIZE_CHAIN_ID);

    numBoundRes = pMetalPro2->numBoundRes;
    if(pMetalPro1->metalAtomList != NULL && pMetalPro2->metalAtomList != NULL)
    {
        for(i = 0 ; i < pMetalPro2->numMetalAtom; i++)
            my_strcpy(pMetalPro1->metalAtomList[i], pMetalPro2->metalAtomList[i], SIZE_ATOM_ELEMENT);
    }
    if(pMetalPro1->resSeqIndex != NULL && pMetalPro2->resSeqIndex != NULL)
    {
        for(i = 0 ; i < pMetalPro2->numBoundRes; i++)
            pMetalPro1->resSeqIndex[i] = pMetalPro2->resSeqIndex[i]; 
    }
    if(isCopyResidue)
    {
        for(i = 0 ; i < numBoundRes; i++)
            CopyResidue(&(pMetalPro1->res[i]),&(pMetalPro2->res[i]));
    }

    if(pMetalPro1->atomEnv != NULL && pMetalPro2->atomEnv != NULL)
    {
        for(i = 0 ; i < pMetalPro2->numMetalAtom; i++)
        {
            CopyAtomEnv(&(pMetalPro1->atomEnv[i]),&(pMetalPro2->atomEnv[i]), false);
        }
    }

}
/*}}}*/
void CopyMetalPro(MetalPro2 *pMetalPro1, MetalPro2 *pMetalPro2, bool isCopyAtomEnv /*= true*/)/*{{{*/
{
    int i;
    int numAtomEnv;
    pMetalPro1->numMetalAtom = pMetalPro2->numMetalAtom;
    pMetalPro1->length = pMetalPro2->length;
    my_strcpy(pMetalPro1->id , pMetalPro2->id, SIZE_CHAIN_ID);
    numAtomEnv = pMetalPro2->numMetalAtom;
    if(isCopyAtomEnv)
    {
        for(i = 0; i < numAtomEnv ; i++)
        {
            if(&(pMetalPro1->atomEnv[i]) != NULL)
            {
                CopyAtomEnv(&(pMetalPro1->atomEnv[i]),&(pMetalPro2->atomEnv[i]));
            }
        }
    }
}
/*}}}*/
void DeleteMetalPro(MetalPro *pMetalPro, int numMetalAtom /*= 0*/, int numBoundRes /* = 0*/)/*{{{*/
    /*****************************************************************************
     * when the memory of 
     *        pMetalPro->atomEnv
     *        pMetalPro->metalAtomList 
     * are  allocated by a constant, it can be freed by specifying numMetalAtom explicitly
     * 2007-04-11, Nanjiang Shu
     ****************************************************************************/
{   
    int i;
    if(numBoundRes == 0) numBoundRes = pMetalPro->numBoundRes;
    if(pMetalPro->res != NULL) 
    {
        for( i = 0 ; i < numBoundRes ; i++)
            DeleteResidue(&(pMetalPro->res[i]));
        delete [] pMetalPro->res;
    }

    if(numMetalAtom == 0) numMetalAtom = pMetalPro->numMetalAtom;
    if(pMetalPro->atomEnv != NULL)
    {
        for( i = 0 ; i < numMetalAtom ; i++)
            DeleteAtomEnv(&(pMetalPro->atomEnv[i]));
        delete [] pMetalPro->atomEnv ;
    }

    if(pMetalPro->metalAtomList != NULL)
    {
        Delete2DArray(pMetalPro->metalAtomList, numMetalAtom);
    }

    if(pMetalPro->resSeqIndex != NULL) 
    {
        delete [] pMetalPro->resSeqIndex;
    }
    InitMetalPro(pMetalPro);

}
/*}}}*/
void DeleteMetalPro(MetalPro2 *pMetalPro, int numMetalAtom /*= 0*/)/*{{{*/
{
    int i ;
    if(numMetalAtom == 0) numMetalAtom = pMetalPro->numMetalAtom;
    if(pMetalPro->atomEnv != NULL)
    {
        for( i = 0 ; i < numMetalAtom ; i++)
            DeleteAtomEnv(&(pMetalPro->atomEnv[i]));
        delete [] pMetalPro->atomEnv ;
    }
    InitMetalPro(pMetalPro);
}
/*}}}*/

void InitPredPro(PredPro *pPredPro)/*{{{*/
{
    strcpy(pPredPro->id,"");
    pPredPro->resSeqIndex = NULL;
    pPredPro->numRes = 0;
    pPredPro->speResSeqIndex = NULL;
    pPredPro->numSpeRes = 0;
}
/*}}}*/
void CopyPredPro(PredPro *pPredPro1, PredPro *pPredPro2)/*{{{*/
{
    int i;
    pPredPro1->numRes = pPredPro2-> numRes;
    pPredPro1->numSpeRes = pPredPro2-> numSpeRes;
    if(pPredPro1->resSeqIndex != NULL && pPredPro2->resSeqIndex != NULL)
    {
        for(i = 0 ; i < pPredPro2->numRes ; i++)
            pPredPro1->resSeqIndex[i] = pPredPro2->resSeqIndex[i];
    }
    if(pPredPro1->speResSeqIndex != NULL && pPredPro2->speResSeqIndex != NULL)
    {
        for(i = 0 ; i < pPredPro2->numSpeRes ; i++)
            pPredPro1->speResSeqIndex[i] = pPredPro2->speResSeqIndex[i];
    }
}
/*}}}*/
void DeletePredPro(PredPro *pPredPro)/*{{{*/
{
    if(pPredPro->resSeqIndex != NULL)
    {
        delete [] pPredPro->resSeqIndex;
    }
    if(pPredPro->speResSeqIndex != NULL)
    {
        delete [] pPredPro->speResSeqIndex;
    }
    InitPredPro(pPredPro);
}
/*}}}*/

void InitAlignFactor(AlignFactor *pAlignFactor)/*{{{*/
{
    pAlignFactor -> score            = 0.0;
    pAlignFactor -> zScore           = 0.0;
    pAlignFactor -> pozScore         = 0.0;
    pAlignFactor -> eValue           = 0.0;
    pAlignFactor -> idt_cnt          = 0;
    pAlignFactor -> identity         = 0.0;
    pAlignFactor -> identity_short   = 0.0;
    pAlignFactor -> sim_cnt          = 0;
    pAlignFactor -> similarity       = 0.0;
    pAlignFactor -> similarity_short = 0.0;
    pAlignFactor -> gap_cnt          = 0;
    pAlignFactor -> gapPercent       = 0.0;

}
/*}}}*/

void InitAtomFreq(AtomFreq *pAtomFreq)/*{{{*/
{
    pAtomFreq->atomCnt = 0;
    strcpy(pAtomFreq->atomName,"");
}
/*}}}*/
void CopyAtomFreq(AtomFreq *to, AtomFreq *from)/*{{{*/
{
    to -> atomCnt = from -> atomCnt;
    strcpy(to->atomName, from->atomName);
}
/*}}}*/

void InitMSA(MSA *pMSA)/*{{{*/
{
    pMSA->numSeq = 0;
    pMSA->alnSeqName  = NULL;
    pMSA->alnSeq  = NULL;
    pMSA->gaplessSeqIndex = NULL;
    pMSA->alnSeqLen = NULL;
    pMSA->seqLen = NULL;
    pMSA->pid =NULL;
}
/*}}}*/
void DeleteMSA(MSA *pMSA)/*{{{*/
{
    int i;
    for(i = 0 ; i < pMSA->numSeq; i++)
        delete [] pMSA->alnSeqName[i];
    delete [] pMSA->alnSeqName;

    for(i = 0 ; i < pMSA->numSeq; i++)
        delete [] pMSA->alnSeq[i];
    delete [] pMSA->alnSeq;

    for(i = 0 ; i < pMSA->numSeq; i++)
        delete [] pMSA->gaplessSeqIndex[i];
    delete [] pMSA->gaplessSeqIndex;

    delete [] pMSA->alnSeqLen;
    delete [] pMSA->seqLen;
    delete [] pMSA->pid;
}
/*}}}*/
void f_neglect_clustalw_header(FILE *fp)/*{{{*/
//neglect the header of clustalw output file
{
    int status_fpos;
    fpos_t pos;
    int linesize = 0;
    int maxline = 50;
    char line[50+1] = "";
    
    fgetline(fp, line ,maxline); //get rid of the first line

    status_fpos = fgetpos(fp,&pos);
    assert ( status_fpos == 0 );
    while((linesize = fgetline ( fp, line ,maxline)) != EOF)
    {
        if(linesize > 0 && line[0] != ' ')
        {
            status_fpos = fsetpos(fp,&pos);
            assert ( status_fpos == 0);
            break;
        }
        else
        {
            status_fpos = fgetpos(fp,&pos);
            assert ( status_fpos == 0 );
        }
    }
}
/*}}}*/
/*read in data from file*/
int GetSeq_ATOM(const char* pdbfilepath, const char *id, char *aaSeq, int *resSer, char *resICode, FILE *fpLog /*= stdout*/)/*{{{*/
/*****************************************************************************
 * get one chain sequence from the ATOM recorde in PDB files
 * return the number of residues in the chain
 ****************************************************************************/
{
    FILE *fpPDBFile = fopen(pdbfilepath,"r");
    checkfilestream(fpPDBFile, pdbfilepath, "r", true);

    int i;
    int maxline = 200;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;

    // get sequences of proteins in ATOM records
    bool isFindChain = false;
    int  atomSerial = -999;
    char atomResName[SIZE_RES_NAME+1] = "";
    int  atomResSeq;
    int  atomResICode;
    int  atomResSeqFormer     = INIT_RESSEQ;
    char atomResICodeFormer   = INIT_ICODE;

    char chainID;
    char title[SIZE_TITLE+1] = "";
    char str[200+1] = "" ;
    int cntRes = 0;

    while(fgetline(fpPDBFile, line, maxline) != EOF)
    {
        if(sscanf(line,"%6s",title) < 1) continue;

        if(strcmp(title,"ENDMDL") == 0)
        {
            fprintf(fpLog,"%s : End of model, take the first model\n",id);
            if(atomSerial == MAX_ATOM_SERIAL)
                fprintf(fpLog,"%s :  atom Serial is possibly overflow, CHECK!\n",id);
            break; // if it is a multi model structure, take the first model
        }
        if(strcmp(title,"ATOM") == 0 || strcmp(title,"HETATM")==0)
        {
            for ( i = 0 ; i < 5 ; i++) str[i] = line[i+6]; str[i] = '\0';
            sscanf(str,"%d",&atomSerial);//serial
        }

        if(strcmp(title,"ATOM")==0) 
        {
            chainID = line[21];
            if(chainID == id[4])
            {
                isFindChain = true;
                //ScanfCoorRecord_Atom(line, &atom);

                for ( i=0 ;i<3;i++) str[i] = line[i+17]; str[i] = '\0';
                sscanf(str,"%s",atomResName);//residue name

                for ( i=0 ;i<4;i++) str[i] = line[i+22]; str[i] = '\0';
                atomResSeq  = atoi(str); //residue seq, it's not necessary corresponding that  in SEQRES
                atomResICode = line[26];

                if(atomResSeq != atomResSeqFormer || atomResICode != atomResICodeFormer)
                {// new residue
                    aaSeq[cntRes]    = AA3To1(atomResName);
                    resSer[cntRes]   = atomResSeq;
                    resICode[cntRes] = atomResICode;
                    cntRes ++;

                    atomResSeqFormer   = atomResSeq;
                    atomResICodeFormer = atomResICode;
                }
            }								   
            else if(isFindChain) //if current chainID != the request chainID, and the chain of interest already found, break
                break;
        }
    }
    fclose(fpPDBFile);
    aaSeq[cntRes]    = '\0';
    resICode[cntRes] = '\0';

    return  cntRes;
}
/*}}}*/
int GetSeq_ATOM(const char* pdbfile, Chain *chain, char* chainIDList,  bool isGetAllChain /*= false*/,bool isIncludeHETATM /*= true*/, int max_seq_length /* = LONGEST_SEQ*/, FILE *fpLog /*= stdout*/)/*{{{*/
/*****************************************************************************
 * get a series of protein sequences with the chainID within chainIDList, from
 * the ATOM recorde in PDB files
 * 2007-04-18
 * ChangeLog: 2007-06-18
 *   mutated amino acid, for example, MSE151 in 1A8O is mutated from MET, this
 *   residue should be regarded as a residue in the sequence, 
 *   now converted to 'X',
 *   it can also be converted to the residue from which it is mutated
 *   return the number of chains read in
 * ChangeLog: 2008-07-07
 *   In some PDB entries (e.g. 2IU4), a chain is divided in two or more pieces
 *   and does not stored continuously, the old algorithm using the bool
 *   variables "isAllChainFound" will ignore the second part of the chain in
 *   the atom coordinate. This bug is fixed now.
 * ChangeLog 2008-09-18 
 *   add the variable isIncludeHETATM, so that whether HETATM record is
 *   included can be regulated
 ****************************************************************************/
{
    FILE *fpPDBFile = fopen(pdbfile,"r");
    checkfilestream(fpPDBFile,pdbfile,"r");

    int numAllRes_PDB = sizeof(PDB3Char_AllRes) / sizeof(char*);

    int i;
    int maxline = 200;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;

    // get sequences of proteins in ATOM records
    int  atomSerial = -999;
    char atomResName[SIZE_RES_NAME+1] = "";
    int  atomResSeq;
    int  atomResICode;
    int  atomResSeqFormer     = INIT_RESSEQ;
    char atomResICodeFormer   = INIT_ICODE;

    char chainID = '\0';
    char chainIDFormer    = '\0'; //initialize chainID as '\0'
    char recordid[SIZE_TITLE+1] = "";
    char str[200+1] = "" ;
    int cntRes = 0;

    int numChainID = strlen(chainIDList);
    Chain *pChain;
    Array1D <bool> isChainNFilled_1darray(numChainID);
    bool *isChainNFilled = isChainNFilled_1darray.array1D;
    for(i = 0 ; i < numChainID; i++) isChainNFilled[i] = false;

    char proteinTitle[SIZE_LINE_PDB+1] = ""; /*title of the protein in the pdb file*/
    char *pTitle = NULL;

    bool isAllChainFound = false;
    int indexChain = -1;

    if(isGetAllChain)
    {
        strcpy(chainIDList,"");
    }
    int cntChain = 0; /*cntChain is used only when isGetAllChain == true*/

    while(fgetline(fpPDBFile, line, maxline) != EOF)
    {
        if(sscanf(line,"%6s",recordid) < 1) continue;

        if(pTitle == NULL && strncmp(recordid, "TITLE", 5) == 0) /*retrieve the protein annotation*/
        {
            my_strcpy(proteinTitle, line+6, SIZE_LINE_PDB);
            pTitle = strtrim(proteinTitle);
        }
        else if (pTitle == NULL && strcmp(proteinTitle, "") == 0 && strncmp(recordid, "COMPND", 6) == 0)
        {
            my_strcpy(proteinTitle, line+6, SIZE_LINE_PDB);
            pTitle = strtrim(proteinTitle);
        }

        if(strcmp(recordid,"ENDMDL") == 0)
        {
            if (fpLog) {
                fprintf(fpLog,"%s : End of model, take the first model\n",pdbfile);
            }
            if(atomSerial == MAX_ATOM_SERIAL) {
                if(fpLog){
                    fprintf(fpLog,"%s :  atom Serial is possibly overflow, CHECK!\n",pdbfile);
                }
            }
            break; // if it is a multi model structure, take the first model
        }
        if(strcmp(recordid,"ATOM") == 0 ||( strcmp(recordid,"HETATM")==0 && isIncludeHETATM))
        {
            for ( i = 0 ; i < 5 ; i++) str[i] = line[i+6]; str[i] = '\0';
            sscanf(str,"%d",&atomSerial);//serial


            chainID = line[21];

            if(!isGetAllChain) 
            {
                if(find(isChainNFilled, isChainNFilled+numChainID, false) == (isChainNFilled + numChainID)) 
                    isAllChainFound = true;
            }

            indexChain = Char2Digit(chainID, chainIDList);

            if(isGetAllChain && chainID != chainIDFormer && indexChain == -1)
            {
                /*it is a new chain*/
                indexChain = cntChain; /*the indexChain of the first new chain is 0*/
                chainIDList[cntChain] = chainID;
                chainIDList[cntChain+1] = '\0';
                cntChain ++;

                chainIDFormer = chainID;
                /*if not a new chain, indexChain return the index to the
                 * current chain*/
            }

            for ( i=0 ;i<3;i++) str[i] = line[i+17]; str[i] = '\0';
            sscanf(str,"%s",atomResName);//residue name
            if(strcmp(recordid,"HETATM") == 0) // if is hetatm, check if it is a mutated amino acid residue
            {
                int idxres  =  0;
                if((idxres = BinarySearch_String((const char*)atomResName, PDB3Char_AllRes, numAllRes_PDB)) == -1)
                {
                    continue;
                }
            }

            if(indexChain >= 0)
            {
                pChain = &(chain[indexChain]);
                if(!isGetAllChain) isChainNFilled[indexChain] = true;

                for ( i=0 ;i<3;i++) str[i] = line[i+17]; str[i] = '\0';
                sscanf(str,"%s",atomResName);//residue name

                for ( i=0 ;i<4;i++) str[i] = line[i+22]; str[i] = '\0';
                atomResSeq  = atoi(str); //residue seq, it's not necessary corresponding that  in SEQRES
                atomResICode = line[26];


                if(atomResSeq != atomResSeqFormer || atomResICode != atomResICodeFormer)
                {// new residue
                    cntRes = pChain->numRes;

                    int idxres = 0;
                    if((idxres = BinarySearch_String((const char*)atomResName, PDB3Char_AllRes, numAllRes_PDB))!= -1)
                    { pChain->aaSeq[cntRes]    = PDB1Char_AllRes[idxres]; }
                    else 
                    { pChain->aaSeq[cntRes] = UNKNOWN_AA; }

                    pChain->resSer[cntRes]   = atomResSeq;
                    pChain->resICode[cntRes] = atomResICode;
                    pChain->chainID = chainID;
                    pChain->numRes ++;

                    if(pChain->numRes >= max_seq_length)
                    {
                        fprintf(stderr,"memory overflow, seqLength >= max_seq_length, file:%s, chainID:%c\n", pdbfile,chainID);
                        assert(pChain->numRes < max_seq_length);
                    }
                    atomResSeqFormer   = atomResSeq;
                    atomResICodeFormer = atomResICode;
                }
            }								   
            else if(isAllChainFound && !isGetAllChain) //if the current chainID is not in the chainIDList, and allChain in chainIDList has already been found, then break 
            {
                //break; /*2008-07-07, Nanjiang, have to scan the whole PDB
                //file to search for an entire chain, since in 2IU4, a single
                //chain is not stored continuously*/
            }
        }
    }
    fclose(fpPDBFile);

    cntChain = 0; // count the number of Chain actually find for resList
    if(!isGetAllChain)
    {
        for(i = 0 ; i < numChainID; i++)
            cntChain += isChainNFilled[i];
    }
    else
    {
        cntChain = strlen(chainIDList);  /*2008-07-07 */
        numChainID = cntChain;/*bug fixed 2008-04-01, the numChainID should be set when isGetAllChain == true, otherwise the seqtype of chains will not be set*/
    }

    for( i = 0 ; i < numChainID; i++)
    {
        chain[i].aaSeq[chain[i].numRes]    = '\0';
        chain[i].resICode[chain[i].numRes] = '\0';
    }

    int sizeTitle = strlen(pTitle);
    for(i = 0 ; i < numChainID; i ++)/*copy the title to each chain*/
    {
        pChain = &(chain[i]);
        pChain->title = new char[sizeTitle+1];
        my_strcpy(pChain->title, pTitle, sizeTitle);
        if (isGetAllChain)
        {
            /*if isGetAllChain, copy the actually retrieved chainIDList to
             * chainIDList*/
            chainIDList[i] = pChain->chainID;
        }
        if(IsDNASeq(pChain->aaSeq, pChain->numRes))
        {
            pChain->seqtype = DNA_SEQ;
        }
        else 
        {
            pChain->seqtype = AA_SEQ;
        }
    }
    if (isGetAllChain)
    {
        chainIDList[numChainID] = '\0';
    }

    return  cntChain;
}
/*}}}*/
int GetDSSPChain(const char* id, Chain* pChain, const char* dsspfilepath)/*{{{*/
/****************************************************************************
 * GetDSSPChain()
 *
 * id should be standardized chain identifier, 
 * Get the aaSeq chain in dssp, which contains the secondary structure and
 * SS-bond information
 * do not allocate memory in this procedure
 * return length of dssp chain, or -1 if dsspfile can not open
 ***************************************************************************/
{
    int cntRes;
    int linesize;
    int maxline = 500;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;
    DSSP_Residue dssp_res;
    char chainID = id[4];

    FILE* fpDSSPFile = fopen(dsspfilepath,"r");
    checkfilestream(fpDSSPFile, dsspfilepath, "r");

    while(fgetline(fpDSSPFile,line ,maxline)!=EOF)
    {
        if(line[2] == '#')
            break;
    }

    cntRes = 0 ; 
    while((linesize = fgetline(fpDSSPFile, line , maxline))!= EOF)
    {
        if(linesize < 1)	continue;
        ScanfDSSPResRecord(line,&dssp_res);
        if((dssp_res.chainID == chainID || chainID == ' ')&& dssp_res.aa != '!')
        {
            if(pChain->aaSeq != NULL)    pChain->aaSeq[cntRes]    = dssp_res.aa;
            if(pChain->resSer != NULL)   pChain->resSer[cntRes]   = dssp_res.resSeq;
            if(pChain->resICode != NULL) pChain->resICode[cntRes] = dssp_res.iCode;
            if(pChain->secStruc != NULL) pChain->secStruc[cntRes] = dssp_res.ss;
            if(pChain->waterAcc != NULL) pChain->waterAcc[cntRes] = dssp_res.acc;
            if(pChain->phi != NULL) pChain->phi[cntRes] = dssp_res.phi;
            if(pChain->psi != NULL) pChain->psi[cntRes] = dssp_res.psi;
            cntRes ++;
        }
    }
    fclose(fpDSSPFile);

    pChain->aaSeq[cntRes] = '\0';
    pChain->chainID       = id[4];
    pChain->numRes        = cntRes;
    my_strcpy(pChain->pdbid, id, SIZE_PDBID);

    if( cntRes >= LONGEST_SEQ)
    {
        fprintf(stderr, "Warning!, memory overflow for pChain-> aaSeq\n");
        fprintf(stderr, "id=%s, cntRes=%d\n", id, cntRes);
        assert ( cntRes < LONGEST_SEQ);
    }
    return cntRes;
}/*}}}*/
int GetSEQMAP(const char* seqmapfilepath, Chain *pChain)/*{{{*/
/*****************************************************************************
 * GetSEQMAP: seqmap is the mapping between aaSeqIndex and resSeq in PDB ATOM
 * record
 * get chain of seqmap from seqmap file
 * do not allocate memory in this procedure
 * return length of seqmap chain, or -1 if dsspfile can not open
 * the order of seqmap column is 
 * "AA", "seqidx", "seqres","icode", "shape", "waterAcc", "dsspSec"
 * AA, seqidx, seqres, icode, four columns are necessary
 * others are optional
 ****************************************************************************/
{
    int status_sscanf ;
    int maxline = 200;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;

    FILE *fp;
    fp = fopen (seqmapfilepath,"r");
    checkfilestream(fp, seqmapfilepath, "r");
    f_neglect_comment(fp);

    char aa;
    int  aaSeqIndex;
    int  resSer;
    char resICode;
    char shape;
    int waterAcc;
    char dsspSec;
    fgetline(fp,line,maxline);
    int cntRes = 0;
    while(fgetline(fp,line,maxline) != EOF)
    {
        status_sscanf = sscanf(line, "%c %d %d %c %c %d %c", &aa, &aaSeqIndex, &resSer, &resICode, &shape, &waterAcc, &dsspSec);
        //assert ( status_sscanf == 4);

        if( resICode == NULL_ICODE)
            resICode = ' ';
        aaSeqIndex --;

        if(pChain->aaSeq != NULL && status_sscanf >=1) pChain->aaSeq[cntRes]    = aa;
        if(pChain->resSer != NULL && status_sscanf >= 3) pChain->resSer[cntRes]   = resSer;
        if(pChain->resICode != NULL && status_sscanf >= 4) pChain->resICode[cntRes] = resICode;
        if(pChain-> shString != NULL && status_sscanf >= 5) pChain->shString[cntRes] = shape;
        if(pChain->waterAcc != NULL && status_sscanf >= 6) pChain->waterAcc[cntRes] = waterAcc;
        if(pChain->secStruc != NULL && status_sscanf >= 7) pChain->secStruc[cntRes] = dsspSec;
        cntRes ++;
    }
    fclose(fp);

    if(pChain->aaSeq != NULL)    pChain->aaSeq[cntRes]    = '\0';
    if(pChain->resICode != NULL) pChain->resICode[cntRes] = '\0';
    if(pChain->shString != NULL)    pChain->shString[cntRes] = '\0';
    if(pChain->secStruc != NULL)  pChain->secStruc[cntRes] = '\0';
    pChain->numRes = cntRes;

    if(cntRes > LONGEST_SEQ)
    {
        printf("Warning! memory overflow for pChain->aaSeq\n");
        assert( cntRes <= LONGEST_SEQ);
    }
    return cntRes ;
}
/*}}}*/
int GetMetalElementList(Element *metalEle, const char metalElementListFile[] /*= ""*/ /*"/data/metal_element_list.dat"*/ )/*{{{*/
    /*****************************************************************************
     * Get the metal element list from file metalElementListFile 
     * element name are converted to uppercase to meet the format of PDB
     * return: number of metal elements
     ****************************************************************************/
{
    char c_metalElementListFile[MAX_PATH+1] = "";
    if(strcmp(metalElementListFile,"") == 0)
    {
        char datadir[MAX_PATH+1] = "";
        GetDataDir(datadir);
        sprintf(c_metalElementListFile,"%s/%s",datadir,"metal_element_list.dat");
    }
    else
        my_strcpy(c_metalElementListFile,metalElementListFile,MAX_PATH);

    FILE* fp = fopen(c_metalElementListFile,"r");
    checkfilestream(fp, c_metalElementListFile, "r", true);
    int maxline = 100;
    char line[100+1] = "";
    int i = 0;
    int stat;
    int linesize;
    while((linesize = fgetline(fp, line ,maxline)) != EOF)
    {
        if(linesize <= 0) continue;
        stat = sscanf(line,"%d %2s",&metalEle[i].atomNum,metalEle[i].name) ;
        my_strupr(metalEle[i].name);
        if(stat == 2) i ++;
    }
    int numMetalEle = i;
    fclose(fp);
    return numMetalEle;
}/*}}}*/
int GetAtom_PDB(const char *pdbfile, Atom *atom, int *pNumAtom, bool isGetAllChain /*= true*/, int max_num_atom /*= MAX_ATOM_SERIAL*/, char* chainIDList /*= ""*/, bool isRestrictByContactAtom /* = true*/, Atom *metalAtom /*= NULL*/, int *pNumMetalAtom /*= NULL*/, int max_num_metal_atom /* = MAX_METAL_CHAIN*/)/*{{{*/
    /*****************************************************************************
     *  Get the coordinate data of all atoms in the pdbfile,
     *  if chainIDList is specified, only retrieve the atoms with chainID within
     *  chainIDList
     *  2007-04-19, Nanjiang Shu
     *  LOG: 2007-04-21 15:42:46 Saturday  Week 15 <nanjiang@casio.fos.su.se>
     *     bug with SelectAltLocAtom:
     *     the alternative location of atom positions have many different cases
     *     exist in PDB,
     *     the first is like in 1DOS, the alternative location of atoms are next to
     *     each other
     *  ATOM    957  N  1ASP A 109      19.045   6.098   9.618  0.60  6.73           N  
     *  ATOM    958  N  2ASP A 109      19.020   6.116   9.614  0.40  8.38           N  
     *  ATOM    959  CA 1ASP A 109      18.353   4.848   9.331  0.60  8.21           C  
     *  ATOM    960  CA 2ASP A 109      18.326   4.866   9.329  0.40 10.12           C  
     *  ATOM    961  C  1ASP A 109      19.277   3.621   9.284  0.60  8.36           C  
     *  ATOM    962  C  2ASP A 109      19.159   3.598   9.569  0.40 11.42           C  
     *     the secondd case is line in ????
     *     the location of all atoms in one residue appear first, and the
     *     alternative location of all atoms follows
     *  The original SelectAltLocAtom can pick out one position for alternative atoms successfully. However, since SelectAltLocAtom
     *  reset the file position to the line below the first reading line, the
     *  next alternative position will be considered as a unique atom position. In
     *  that case, all atoms, including the alternative atom position will be read
     *  in as individual atoms.
     *
     *  Solution: read in one Residue at a time, and then with the residue,
     *  determine which atom to be used.
     *
     *  supply another functionality of GetAtom_PDB, that is, MetalAtom can be read
     *  in from this procedure
     *  if metalAtom = NULL, metal atom record will not be read in
     *  if atom = NULL, atom record will not be read in.
     ****************************************************************************/
{
    // read the coordinate data of all ATOM records
    if(strcmp(chainIDList,"") != 0)
    {
        isGetAllChain = false;
    }

    int i;
    // read in metal element list/*{{{*/
    int numMetalEle = 0 ;
    Array1D <Element> metalEle_1darray(NUM_METAL_ELEMENT);
    Element * metalEle = metalEle_1darray.array1D;
    Array2D <char> metalEleList_2darray(numMetalEle, SIZE_ATOM_ELEMENT+1);
    char **metalEleList = metalEleList_2darray.array2D;

    if(metalAtom != NULL)
    {
        numMetalEle = GetMetalElementList( metalEle);
        for(i = 0 ; i < numMetalEle; i++)
            my_strcpy(metalEleList[i], metalEle[i].name, SIZE_ATOM_ELEMENT);
    }
    /*}}}*/

    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;

    FILE *fpPDBFile = fopen(pdbfile,"r");
    checkfilestream(fpPDBFile, pdbfile,"r");
    char record[SIZE_RECORD_ID+1] = "";
    char atomResName[SIZE_RES_NAME+1] = "";
    char atomOrigName[SIZE_ATOM_ORIGNAME+1] = "";
    char chainID =  ' ';

    Array1D <Atom> altAtoms_1darray(MAX_ALTATOM);
    Array2D <char> altAtomLines_1darray(MAX_ALTATOM,SIZE_LINE_PDB+1);

    int cntAtom = 0 ;
    int cntMetalAtom = 0;
    int atomSerial = 0;

    while((linesize = fgetline(fpPDBFile, line, maxline) ) != EOF)
    {
        if(linesize <= 0 || sscanf(line,"%6s",record) != 1) continue;
        else if(strcmp(record,"ANISOU") == 0) continue;
        else if(strcmp(record,"ENDMDL") == 0)
        {
            fprintf(stderr,"%s: end of model, take only the first model\n", pdbfile);
            if(atomSerial == MAX_ATOM_SERIAL)
                fprintf(stderr,"%s: atom serial number > MAX_ATOM_SERIAL (%d), possible overflow, CHECK!\n",pdbfile,MAX_ATOM_SERIAL);
            break; // if it is multimodel, take only the first model
        }
        else if(strcmp(record,"ATOM")==0 && atom != NULL) // if record is ATOM
        {
            ScanfAtomSerial(line, atomSerial); 
            chainID = line[POS_CHAIN_ID];
            if(isGetAllChain || IsInCharSet(chainID, chainIDList))
            {
                ScanfAtomResName(line, atomResName);
                if(BinarySearch_String((const char*)atomResName, PDB3CharAA_alphabet,strlen(PDB1CharAA_alphabet)) != -1) // is an amino acid residue residue
                {	
                    ScanfAtomOrigName(line, atomOrigName);
                    if(IsInContactAtomList(atomOrigName) || !isRestrictByContactAtom)
                    {
                        ScanfCoorRecord_Atom(line,&atom[cntAtom]);
                        if(atom[cntAtom].altLoc != ' ') // if alternative atoms exist
                            cntAtom = SelectAltLocAtom(atom, cntAtom, fpPDBFile);	
                        else 
                            cntAtom ++;

                        if(cntAtom > max_num_atom)
                        {
                            fprintf(stderr,"cntAtom > max_num_atom (= %d), memory overflow\n", max_num_atom);
                            assert (cntAtom <= max_num_atom);
                        }
                    }
                }
            }
        }
        else if(strcmp(record, "HETATM") == 0 && metalAtom != NULL)
        {
            ScanfAtomSerial(line, atomSerial); 
            ScanfCoorRecord_Atom(line,&metalAtom[cntMetalAtom]);
            if(IsMetalAtom(&metalAtom[cntMetalAtom],metalEleList,numMetalEle))
            {
                if(metalAtom[cntMetalAtom].altLoc != ' ')  // if has alternative atoms
                    cntMetalAtom = SelectAltLocAtom(metalAtom, cntMetalAtom,fpPDBFile);

                else
                    cntMetalAtom ++;

                if(cntMetalAtom > max_num_metal_atom)
                {
                    fprintf(stderr,"numMetal > max_num_metal_atom (= %d), memory overflow\n", max_num_metal_atom);
                    assert (cntMetalAtom <= max_num_metal_atom);
                }
            }
        }
    }
    fclose(fpPDBFile);

    *pNumAtom = cntAtom;
    if(metalAtom != NULL) *pNumMetalAtom = cntMetalAtom;

    return *pNumAtom;
}
/*}}}*/
int GetBkFreq(const char *infile, double *P, const char *alphabet, int num_aa /*= NUM_AA*/)/*{{{*/
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
    double f;
    int daa;
    f_neglect_comment(fpin);
    char *pstr = NULL;
    while((linesize = fgetline(fpin, line, maxline)) != EOF)
    {
        pstr = strtrim(line);
        if(sscanf(pstr,"%c %lf", &ch, &f) == 2)
        {
            daa = Char2Digit(ch, alphabet);
            if(daa >= 0)
            {
                P[daa] = f;
                i ++;
            }
        }
    }
    fclose(fpin);
    int numAA = i;
    if( numAA != num_aa)
    {
        fprintf(stderr, "Warning! number of frequencies given in background frequency file %s is not equals to %d\n", infile, num_aa);
    }
    
    //check if it is the percentage
    double sum = 0.0; for(i = 0; i < numAA ; i ++) sum+= P[i];
    if(fabs(sum - 100.0) < 1.0) // if is percentage, convert to frequency
    { for(i = 0 ; i < numAA; i ++) P[i]  /= 100.0; }
    else if (fabs(sum-1.0) > 0.1)
    {
        fprintf(stderr, "Warning! the given background frequency is not unified in file %s\n", infile);
    }
    return numAA;
}
/*}}}*/

int8 GetBinaryMODMPara(const char *file, int &sizeAlphabet, int &length, int8 &typeProfile)/*{{{*/
    /*get the parameters of the binary MODM file, return the typeProfile
     * type = 0: ProfileSAD
     * type = 1: ProfileSADByte , the default type
     *
     * type = 2: Profile
     * type = 3: ProfileByte, 2009-11-05*/
{
    FILE *fpin = fopen(file,"rb");
    size_t nread = 0;
    if (checkfilestream(fpin, file, "rb") == -1)
    { return -1; }
    nread = fread(&typeProfile, sizeof(int8), 1, fpin);
    if (nread != 1)
    { fprintf(stderr,"File read error! function: %s. File: %s.\n","GetBinaryMODMPara",file); }
    nread = fread(&sizeAlphabet, sizeof(int), 1, fpin);
    if (nread != 1)
    { fprintf(stderr,"File read error! function: %s. File: %s.\n","GetBinaryMODMPara",file); }
    nread = fread(&length, sizeof(int), 1, fpin);
    if (nread != 1)
    { fprintf(stderr,"File read error! function: %s. File: %s.\n","GetBinaryMODMPara",file); }
    fclose(fpin);
    return typeProfile;
}/*}}}*/
template <class T> int GetBinaryMODM(const char *file, char *alphabet, int &length, T *profile, double *parameter, int8 &typeProfile)/*{{{*/
/*Read Binary MODM by given the filename*/
{
    FILE *fpin = NULL;
    fpin = fopen(file, "rb");
    if (checkfilestream(fpin, file, "rb") == -1) {
        return -1;
    }else{
        unsigned int sizeAlphabet = 0;
        size_t nread =0;/*return value for fread, the number of bytes read in*/
        if ((nread = fread(&typeProfile, sizeof(int8), 1, fpin)) != 1){ /*read the typeProfile*/ 
            fprintf(stderr,"fread failed for typeProfile in GetBinaryMODM\n");
            return -1;
        }
        if ((nread = fread(&sizeAlphabet, sizeof(int), 1, fpin))!=1){ /*read the sizeAlphabet*/
            fprintf(stderr,"fread failed for sizeAlphabet in GetBinaryMODM\n");
            return -1;
        }
        if ((nread = fread(&length, sizeof(int),1,fpin))!=1){/*read the numRes of the profile*/
            fprintf(stderr,"fread failed for length in GetBinaryMODM\n");
            return -1;
        }else if (length < 0){
            fprintf(stderr,"Error! length (%d) of the profile for file %s <=0. Ignore.\n", length, file);
            return -1;
        }
        if ((nread = fread(alphabet, sizeof(char), sizeAlphabet, fpin))!=sizeAlphabet){  /*read the alphabet of the profile*/ 
            fprintf(stderr,"fread failed for alphabet in GetBinaryMODM\n");
            return -1;
        }
        if ((nread = fread(parameter, sizeof(double), 8, fpin))!=8){  /*read the parameter of the profile*/ 
            fprintf(stderr,"fread failed for parameter in GetBinaryMODM\n");
            return -1;
        }
        if ((nread = fread(profile, sizeof(T), length, fpin))!= (unsigned int) (length)){   /*read the profile*/ 
            fprintf(stderr,"fread failed for profile in GetBinaryMODM\n");
            return -1;
        }
        fclose(fpin);
        return 0;
    }
}
template int GetBinaryMODM <Profile> (const char *file, char *alphabet, int&length, Profile * profile, double *parameter, int8 &typeProfile);
template int GetBinaryMODM <ProfileByte> (const char *file, char *alphabet, int&length, ProfileByte * profile, double *parameter, int8 &typeProfile);
template int GetBinaryMODM <ProfileSAD> (const char *file, char *alphabet, int&length, ProfileSAD * profile, double *parameter, int8 &typeProfile);
template int GetBinaryMODM <ProfileSADByte> (const char *file, char *alphabet, int&length, ProfileSADByte * profile, double *parameter, int8 &typeProfile);
/*}}}*/
template <class T> int GetBinaryMODM_fp(FILE *fpin, long readsize, char *alphabet, int &length, T *profile, double *parameter, int8 &typeProfile)/*{{{*/
/* Read in binary modm by given the FILE handler
 * created 2011-10-13*/
{
    if (fpin != NULL){
        int sizeAlphabet = 0;
        long numByteReadIn=0;
        size_t nread =0;/*return value for fread, the number of bytes read in*/
        if ((nread = fread(&typeProfile, sizeof(int8), 1, fpin)) != 1){ /*read the typeProfile*/ 
            fprintf(stderr,"fread failed for typeProfile in GetBinaryMODM_fp\n");
            return -1;
        }
        numByteReadIn += nread*sizeof(int8);
        if ((nread = fread(&sizeAlphabet, sizeof(int), 1, fpin))!=1){ /*read the sizeAlphabet*/
            fprintf(stderr,"fread failed for sizeAlphabet in GetBinaryMODM_fp\n");
            return -1;
        } else if (sizeAlphabet < 0){
            fprintf(stderr,"Error. sizeAlphabet (%d) <=0 in GetBinaryMODM_fp\n", sizeAlphabet);
            return -1;
        }
        numByteReadIn += nread*sizeof(int);
        if ((nread = fread(&length, sizeof(int),1,fpin))!=1){/*read the numRes of the profile*/
            fprintf(stderr,"fread failed for length in GetBinaryMODM_fp\n");
            return -1;
        }
        numByteReadIn += nread*sizeof(int);
        if ((nread = fread(alphabet, sizeof(char), sizeAlphabet, fpin))!=(unsigned int) sizeAlphabet){  /*read the alphabet of the profile*/ 
            fprintf(stderr,"fread failed for alphabet in GetBinaryMODM_fp\n");
            return -1;
        }
        numByteReadIn += nread*sizeof(char);
        if ((nread = fread(parameter, sizeof(double), 8, fpin))!=8){  /*read the parameter of the profile*/ 
            fprintf(stderr,"fread failed for parameter in GetBinaryMODM_fp\n");
            return -1;
        }
        numByteReadIn += nread*sizeof(double);
        if ((nread = fread(profile, sizeof(T), length, fpin))!=(unsigned int)length){   /*read the profile*/ 
            fprintf(stderr,"fread failed for profile in GetBinaryMODM_fp, sizeof(T)=%lu, length=%d \n", sizeof(T), length);
            return -1;
        }
        numByteReadIn += nread*sizeof(T);
        if (numByteReadIn != readsize){
            fprintf(stderr,"numByteReadIn(%ld) != readsize(%ld) in GetBinaryMODM_fp\n", numByteReadIn, readsize);
            return -1;
        }
        return 0;
    } else{
        return -1;
    }
}
template int GetBinaryMODM_fp <Profile> (FILE *fpin, long readsize, char *alphabet, int&length, Profile * profile, double *parameter, int8 &typeProfile);
template int GetBinaryMODM_fp <ProfileByte> (FILE *fpin,long readsize, char *alphabet, int&length, ProfileByte * profile, double *parameter, int8 &typeProfile);
template int GetBinaryMODM_fp <ProfileSAD> (FILE *fpin, long readsize, char *alphabet, int&length, ProfileSAD * profile, double *parameter, int8 &typeProfile);
template int GetBinaryMODM_fp <ProfileSADByte> (FILE *fpin, long readsize, char *alphabet, int&length, ProfileSADByte * profile, double *parameter, int8 &typeProfile);
/*}}}*/

int GetBinaryFragShort5(const char *file, int &fragFileType, char **idList, int &numID, int &maxSizeID, int &length, short *posTar, short *numCan, int &totalFragCan, FragCanShort5 *fragCan)/*{{{*/
/* Read in the frag file in binary format, 2009-06-23 
 * the mmeory is allocated from the outer procedure 2009-06-25 
 * */
{
    FILE *fpin = fopen(file,"rb");
    if (checkfilestream(fpin, file, "rb") == -1)
    { return -1; }
    int i = 0;

    size_t nread =0;/*return value for fread, the number of bytes read in*/

    nread = fread(&fragFileType, sizeof(int), 1, fpin); /*read the fragFileType*/ 
    if (nread != 1)
    { fprintf(stderr,"File read error! function: %s. File: %s.\n","GetBinaryFragShort5",file); }

    nread = fread(&numID, sizeof(int), 1, fpin); /*read the number of unique ids*/
    if (nread != 1)
    { fprintf(stderr,"File read error! function: %s. File: %s.\n","GetBinaryFragShort5",file); }

    nread = fread(&maxSizeID, sizeof(int),1,fpin);/*read the maximal size of ids*/
    if (nread != 1)
    { fprintf(stderr,"File read error! function: %s. File: %s.\n","GetBinaryFragShort5",file); }

    nread = fread(&length, sizeof(int), 1, fpin);  /*read length of the target sequence*/ 
    if (nread != 1)
    { fprintf(stderr,"File read error! function: %s. File: %s.\n","GetBinaryFragShort5",file); }

    nread = fread(&totalFragCan, sizeof(int), 1, fpin);   /*read the total number of candidate fragments*/ 
    if (nread != 1)
    { fprintf(stderr,"File read error! function: %s. File: %s.\n","GetBinaryFragShort5",file); }

    /*read the idList*/
    for (i = 0; i < numID; i ++)
    {
        nread = fread(idList[i], sizeof(char), size_t(maxSizeID), fpin );
        if (nread != size_t(maxSizeID))
        { fprintf(stderr,"File read error! function: %s. File: %s.\n","GetBinaryFragShort5",file); }
    }
    /*read the candidate fragments for each target fragment*/
    int start = 0; /*the start index for fragCan*/
    for (i = 0; i < length; i ++)
    {
        nread = fread(&(posTar[i]), sizeof(short), 1, fpin);
        if (nread != 1)
        { fprintf(stderr,"File read error! function: %s. File: %s.\n","GetBinaryFragShort5",file); }
        nread = fread(&(numCan[i]), sizeof(short), 1, fpin);
        if (nread != 1)
        { fprintf(stderr,"File read error! function: %s. File: %s.\n","GetBinaryFragShort5",file); }
        FragCanShort5 *pFragCan = &(fragCan[start]);
        nread = fread(pFragCan, sizeof(FragCanShort5), size_t(numCan[i]), fpin);
        if (nread != size_t (numCan[i]))
        { fprintf(stderr,"File read error! function: %s. File: %s.\n","GetBinaryFragShort5",file); }
        start = start + numCan[i];
    }
    fclose(fpin);
    return 0;
}/*}}}*/
int GetBinaryFragShort6(const char *file, int &fragFileType, char **idList, int &numID, int &maxSizeID, int &length, short *posTar, short *numCan, int &totalFragCan, FragCanShort6 *fragCan)/*{{{*/
/* Read in the frag file in binary format, 2009-06-23 
 * the mmeory is allocated from the outer procedure 2009-06-25 
 * */
{
    FILE *fpin = fopen(file,"rb");
    if (checkfilestream(fpin, file, "rb") == -1)
    { return -1; }
    int i = 0;
    size_t nread =0;/*return value for fread, the number of bytes read in*/

    nread=fread(&fragFileType, sizeof(int), 1, fpin); /*read the fragFileType*/ 
    if (nread != 1)
    { fprintf(stderr,"File read error! function: %s. File: %s.\n","GetBinaryFragShort6",file); }
    nread=fread(&numID, sizeof(int), 1, fpin); /*read the number of unique ids*/
    if (nread != 1)
    { fprintf(stderr,"File read error! function: %s. File: %s.\n","GetBinaryFragShort6",file); }
    nread=fread(&maxSizeID, sizeof(int),1,fpin);/*read the maximal size of ids*/
    if (nread != 1)
    { fprintf(stderr,"File read error! function: %s. File: %s.\n","GetBinaryFragShort6",file); }
    nread=fread(&length, sizeof(int), 1, fpin);  /*read length of the target sequence*/ 
    if (nread != 1)
    { fprintf(stderr,"File read error! function: %s. File: %s.\n","GetBinaryFragShort6",file); }
    nread=fread(&totalFragCan, sizeof(int), 1, fpin);   /*read the total number of candidate fragments*/ 
    if (nread != 1)
    { fprintf(stderr,"File read error! function: %s. File: %s.\n","GetBinaryFragShort6",file); }

    /*read the idList*/
    for (i = 0; i < numID; i ++)
    {
        nread = fread(idList[i], sizeof(char), size_t(maxSizeID), fpin );
        if (nread != size_t(maxSizeID))
        { fprintf(stderr,"File read error! function: %s. File: %s.\n","GetBinaryFragShort6",file); }
    }
    /*read the candidate fragments for each target fragment*/
    int start = 0; /*the start index for fragCan*/
    for (i = 0; i < length; i ++)
    {
        nread = fread(&(posTar[i]), sizeof(short), 1, fpin);
        nread = fread(&(numCan[i]), sizeof(short), 1, fpin);
        FragCanShort6 *pFragCan = &(fragCan[start]);
        nread = fread(pFragCan, sizeof(FragCanShort6), size_t(numCan[i]), fpin);
        start = start + numCan[i];
    }
    fclose(fpin);
    return 0;
}/*}}}*/
int GetBinaryFragInt5(const char *file, int &fragFileType, char **idList, int &numID, int &maxSizeID, int &length, short *posTar, short *numCan, int &totalFragCan, FragCanInt5 *fragCan)/*{{{*/
/* Read in the frag file in binary format, 2009-06-23 
 * the mmeory is allocated from the outer procedure 2009-06-25 
 * */
{
    FILE *fpin = fopen(file,"rb");
    if (checkfilestream(fpin, file, "rb") == -1) { 
        return -1; 
    }else{
        int i = 0;
        size_t nread = 0;

        nread = fread(&fragFileType, sizeof(int), 1, fpin); /*read the fragFileType*/ 
        nread = fread(&numID, sizeof(int), 1, fpin); /*read the number of unique ids*/
        nread = fread(&maxSizeID, sizeof(int),1,fpin);/*read the maximal size of ids*/
        nread = fread(&length, sizeof(int), 1, fpin);  /*read length of the target sequence*/ 
        nread = fread(&totalFragCan, sizeof(int), 1, fpin);   /*read the total number of candidate fragments*/ 

        /*read the idList*/
        for (i = 0; i < numID; i ++)
        {
            nread = fread(idList[i], sizeof(char), size_t(maxSizeID), fpin );
        }
        /*read the candidate fragments for each target fragment*/
        int start = 0; /*the start index for fragCan*/
        for (i = 0; i < length; i ++)
        {
            nread = fread(&(posTar[i]), sizeof(short), 1, fpin);
            nread = fread(&(numCan[i]), sizeof(short), 1, fpin);
            FragCanInt5 *pFragCan = &(fragCan[start]);
            nread = fread(pFragCan, sizeof(FragCanInt5), size_t(numCan[i]), fpin);
            start = start + numCan[i];
        }
        fclose(fpin);
        return 0;
    }
}/*}}}*/
int GetBinaryFragInt6(const char *file, int &fragFileType, char **idList, int &numID, int &maxSizeID, int &length, short *posTar, short *numCan, int &totalFragCan, FragCanInt6 *fragCan)/*{{{*/
/* Read in the frag file in binary format, 2009-06-23 
 * the mmeory is allocated from the outer procedure 2009-06-25 
 * */
{
    FILE *fpin = fopen(file,"rb");
    if (checkfilestream(fpin, file, "rb") == -1)
    { return -1; }
    int i = 0;
    size_t nread =0;/*return value for fread, the number of bytes read in*/

    nread = fread(&fragFileType, sizeof(int), 1, fpin); /*read the fragFileType*/ 
    nread = fread(&numID, sizeof(int), 1, fpin); /*read the number of unique ids*/
    nread = fread(&maxSizeID, sizeof(int),1,fpin);/*read the maximal size of ids*/
    nread = fread(&length, sizeof(int), 1, fpin);  /*read length of the target sequence*/ 
    nread = fread(&totalFragCan, sizeof(int), 1, fpin);   /*read the total number of candidate fragments*/ 

    /*read the idList*/
    for (i = 0; i < numID; i ++)
    {
        nread = fread(idList[i], sizeof(char), size_t(maxSizeID), fpin );
    }
    /*read the candidate fragments for each target fragment*/
    int start = 0; /*the start index for fragCan*/
    for (i = 0; i < length; i ++)
    {
        nread = fread(&(posTar[i]), sizeof(short), 1, fpin);
        nread = fread(&(numCan[i]), sizeof(short), 1, fpin);
        FragCanInt6 *pFragCan = &(fragCan[start]);
        nread = fread(pFragCan, sizeof(FragCanInt6), size_t(numCan[i]), fpin);
        start = start + numCan[i];
    }
    fclose(fpin);
    return 0;
}/*}}}*/
int ReadMSA(const char* msafile, MSA *pMSA, int format /*= CLUSTALW*/)/*{{{*/
{
    int numSeq = 0;
    switch( format )
    {
        case CLUSTALW: numSeq = ReadMSA_clustalw(msafile, pMSA); 
                   break;

        default  : fprintf(stderr, "wrong msa format! see help for accepted file format" );
                   numSeq = -1;
                   break;
    }	
    return numSeq;
}
/*}}}*/
int ReadMSA_clustalw(const char* msafile, MSA *pMSA)/*{{{*/
/*****************************************************************************
 * Read in multiple sequence alignments in clustalw format, 
 * the multiple sequences alignments are arranged by structure MSA
 * 2008-01-28, Nanjiang
 ****************************************************************************/
{
    FILE *fpmsa;
    fpmsa = fopen(msafile, "r");
    if (checkfilestream(fpmsa, msafile, "r") == -1)
    {
        return -1;
    }
    int i,j;
    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D; 

    //get the number of alinged sequeces
    f_neglect_clustalw_header(fpmsa);
    int cntSeq = 0;
    int numSeq = 0;
    while((linesize = fgetline(fpmsa,line,maxline)) != EOF)
    {
        if(line[0] == ' ' || linesize <= 0 ) break;
        else cntSeq ++;
    }
    numSeq = cntSeq ;

    Array1D <char> name_1darray(maxline);
    Array1D <char> seq_1darray(maxline);
    char *name = name_1darray.array1D;
    char *seq = seq_1darray.array1D;

    Array1D <string> seqname_1darray(numSeq);
    Array1D <string> alnseq_1darray(numSeq);
    string* seqname = seqname_1darray.array1D;
    string* alnseq = alnseq_1darray.array1D;

    rewind(fpmsa);
    f_neglect_clustalw_header(fpmsa);

    cntSeq = 0;
    while((linesize = fgetline(fpmsa,line,maxline)) != EOF)
    {
        if(linesize > 0 && line[0] != ' ' )
        {
            sscanf(line,"%s %s", name, seq);
            seqname[cntSeq] = name;
            alnseq[cntSeq].append(seq);
            cntSeq ++;
        }
        else if (linesize == 0)
        {
            cntSeq = 0;
        }
    }
    fclose(fpmsa);
    
    //copy the aligned sequences from string to msa
    pMSA->numSeq = numSeq;
    pMSA->alnSeqName = new char*[numSeq];
    pMSA->alnSeq = new char*[numSeq];
    pMSA->gaplessSeqIndex = new int*[numSeq];
    pMSA->alnSeqLen = new int[numSeq];
    pMSA->seqLen = new int[numSeq];
    pMSA->pid = new float[numSeq];
    for(i = 0 ; i < numSeq ; i++)
    {
        pMSA->alnSeqName[i] = new char[seqname[i].length()+1];
        pMSA->alnSeq[i] = new char[alnseq[i].length()+1];
        pMSA->gaplessSeqIndex[i] = new int [alnseq[i].length()];
        my_strcpy(pMSA->alnSeqName[i],seqname[i].c_str(), seqname[i].length());
        my_strcpy(pMSA->alnSeq[i],alnseq[i].c_str(), alnseq[i].length());
        int cntAA = 0;
        int alnSeqLen = alnseq[i].length();
        for(j = 0 ; j < alnSeqLen ; j++)
        {
            if(pMSA->alnSeq[i][j] != CHAR_INDEL)
            {
                pMSA->gaplessSeqIndex[i][j] = cntAA;
                cntAA ++;
            }
            else
                pMSA->gaplessSeqIndex[i][j] = DIGIT_INDEL;
        }
        pMSA->seqLen[i] = cntAA;
        pMSA->alnSeqLen[i] = alnseq[i].length();
    }

    int indexQuery;
    int indexKeySeq;
    indexQuery = LinearSearch_String("Query", (const char **)pMSA->alnSeqName, pMSA->numSeq);
    if(indexQuery >= 0)
        indexKeySeq = indexQuery;
    else
        indexKeySeq = 0;

    int cntIDT;
    for(i = 0 ; i < numSeq ; i++)
    {
        cntIDT = 0;
        for(j = 0 ; j < pMSA->alnSeqLen[i]; j++)
        {
            if(pMSA->alnSeq[i][j] == pMSA->alnSeq[indexKeySeq][j] && pMSA->alnSeq[i][j] != CHAR_INDEL)
                cntIDT ++;
        }
        pMSA->pid[i] = float(cntIDT) / float(min(pMSA->seqLen[indexKeySeq], pMSA->seqLen[i]))*100.0;
    }

    return numSeq;
}
/*}}}*/

void GetDBFPList( vector <FILE*> &fpList, string dbname, int maxdbfileindex)/*{{{*/
{
    for (int i=0;i<=maxdbfileindex;i++){
        string filename=dbname + int2string(i) + string(".db");
        FILE *fp = fopen(filename.c_str(),"r");
        if (fp != NULL){
            fpList.push_back(fp);
        }else{
            fprintf(stderr,"Failed to open dbfile %s\n", filename.c_str());
            exit(1);
        }
    }
}/*}}}*/
int ReadDatabaseIndex(string dbname, map<string, dbindex> &dbindexmap,int &maxdbfileindex )/*{{{*/
{
    string indexfile=dbname+string(".index");
    ifstream ifp (indexfile.c_str());
    string line;
    dbindex tmpdbindex;
    if (ifp.is_open()){
        while (ifp.good()) {
            getline(ifp, line);
            if (line.size() > 0 && line.substr(0,3) != string("DEF")){
                char *id = new char [ line.size()] ; 
                if (sscanf( line.c_str(), "%s %d %ld %ld", id, &(tmpdbindex.dbfileindex)
                        ,&(tmpdbindex.offset),&(tmpdbindex.size)) != 4){
                    fprintf(stderr,"Read dbindex failed for %s\n", dbname.c_str());
                    return -1;
                } else {
                    dbindexmap.insert(pair<string, dbindex>(string(id), tmpdbindex));
                }
                delete [] id;
            }
        }
        ifp.close();
        return dbindexmap.size();
    } else {
        return -1;
    }
}/*}}}*/
int ReadInDatabase(const char *id, const char *frag_acc_path, const char* qijpath, const char* modmpath, int qijformat, int modmformat, int fragaccformat, int **matFrag, int **matQij, int **matMODM, int **matMerged, float *matFragScore1, float *matFragScore2,  float *matQijScore1, float *matQijScore2, float *matMODMScore1, float *matMODMScore2, char *aaSeq, int *digitAASeq, char *shapeSeq, char *dsspSecSeq, int type_dataset, bool isReadBinaryFile, int NPer_Frag_Database, int ratioScheme, int mergeSide, int dsspMapMethod /*= 0*/, int8 typeProfile /*=1*/ )/*{{{*/
/*****************************************************************************
 * read in the dataset for a chain
 * including FragAcc matrix, Qij matrix and modm matrix
 * the type_dataset and isReadBinaryFile should be given
 * return the length of the chain if successful
 * return -1 for failure
 ****************************************************************************/
{
    int i,j;
    Array1D <int> sumProfileFrag_1darray(LONGEST_SEQ);
    int *sumProfileFrag = sumProfileFrag_1darray.array1D; /*sum of the profile at each position for fragacc matrix*/
    int lengthSeqMODM = 0;
    int lengthSeqQij = 0;
    char fragMatFile[MAX_PATH+1] = ""; /* frag acc matrix file */
    char QijFile[MAX_PATH+1]     = ""; /* Qij matrix file      */
    char modmFile[MAX_PATH+1]    = ""; /* modm file            */


    if (!isReadBinaryFile)/*{{{*/
    {
        if (ratioScheme == 0 && NPer_Frag_Database == 0)
        {   /*if ratioScheme == 0 and NPer_Frag_Database == 0, the Qij file is used as fragacc file, so that the weight is 100% on Qij file*/
            GetQijFilePath(id, fragMatFile, qijpath, qijformat, isReadBinaryFile);
        }
        else 
        {
            GetFragMatFilePath(id, fragMatFile, frag_acc_path, fragaccformat, isReadBinaryFile);
        }
        for (i=0; i < LONGEST_SEQ; i++) { sumProfileFrag[i] = 0; }
        ReadInProfile(fragMatFile, matFrag, NULL, NULL, NULL, NULL, matFragScore1, matFragScore2, sumProfileFrag, LONGEST_SEQ);

        //open Qij profiles for reading rows whose sum of components are zero
        GetQijFilePath(id, QijFile, qijpath, qijformat, isReadBinaryFile);
        lengthSeqQij = ReadInProfile(QijFile, matQij, NULL, NULL, NULL, NULL, matQijScore1, matQijScore2, NULL, LONGEST_SEQ);
        if(lengthSeqQij < 0) { return -1 ; }

        //open normal psi-blast percentage profiles, modm matrix, in tuping's format
        GetMODMFilePath_Tuping(id, modmFile, modmpath, modmformat, isReadBinaryFile);
        lengthSeqMODM = ReadInProfile(modmFile, matMODM, aaSeq,shapeSeq, NULL, dsspSecSeq, matMODMScore1, matMODMScore2, NULL, LONGEST_SEQ);
        if(lengthSeqMODM < 0) { return -1 ; }
        if(lengthSeqQij != lengthSeqMODM)
        {
            fprintf(stderr, "%s, length not consistent! lengthSeqQij = %d, lengthSeqMODM = %d\n", id, lengthSeqQij, lengthSeqMODM);
            /*assert(lengthSeqMODM == lengthSeqQij); */ /*just print the warning
                message, it is not necessary to abort the program, since i used the asSeqIndex to actually read the profile to M, 2009-06-09*/
        }
        aaSeq[lengthSeqMODM] = '\0';
        shapeSeq[lengthSeqMODM] = '\0';
        dsspSecSeq[lengthSeqMODM] = '\0';
        for(i = 0; i < lengthSeqMODM; i ++)
        { dsspSecSeq[i] = DSSP8to3(dsspSecSeq[i], dsspMapMethod); /*map 8-state dsspsec to 3-state, H,G-->Helix, E-->strand, others-->Random coil*/ }

        TreatAllZeroFij(matMODM, lengthSeqMODM, matMODMScore1, aaSeq, AAAlphabet_Tuping);

#ifdef DEBUG_READ_MODM
        fprintf(stdout,"modm from %s\n", modmFile);
        for(i = 0; i < lengthSeqMODM; i ++)
        {
            WriteMODMProfile(i+1, aaSeq[i], matMODM[i], 0,0, AAAlphabet_Tuping, stdout);
        }
        fprintf(stdout,"\n");
#endif
    }/*}}}*/
    else /*isReadBinaryFile == true*//*{{{*/
    {
        Array1D <ProfileSADByte> profile_1darray(LONGEST_SEQ);
        ProfileSADByte *profile = profile_1darray.array1D;
        char alphabet[50+1] = "";
        int length = 0; /*number of profiles read in*/
        double parameter[8];
        int aaSeqIndex = 0;

        /*2007-11-11, Nanjiang, changed to read in Qij first and then if fragacc exist, replace it with fragacc*/
        /*read in Qij matrix file, for those residue positions does not have fragacc profile, using the Qij profile instead*/
        GetQijFilePath(id, QijFile, qijpath, qijformat, isReadBinaryFile);    /*{{{*/
        if (GetBinaryMODM(QijFile, alphabet, length, profile, parameter, typeProfile) == -1)
        {
            fprintf(stderr, "can not open QijFile %s\n", QijFile);
            return -1;
        }
        lengthSeqQij = length ;
        for(i = 0 ; i < length; i ++)
        {
            aaSeqIndex = profile[i].aaSeqIndex-1;
            if (aaSeqIndex < 0)
            {
                fprintf(stderr,"Error aaSeqIndex (=%d) < 0, %s : %d  \n", aaSeqIndex, QijFile, i);
                assert(aaSeqIndex >= 0);
            }
            for (j=0; j<20; j++) {matQij[aaSeqIndex][j] = int(profile[i].p[j]); }
            matQijScore1[aaSeqIndex] = profile[i].score1; /*change i to aaSeqIndex, 2008-02-18, Nanjiang*/
            matQijScore2[aaSeqIndex] = profile[i].score2;
        }/*}}}*/

        /*read in frag acc matrix file, replace the Psi_blosum_Frag matrix with fragacc profiles */
        if (ratioScheme == 0 && NPer_Frag_Database == 0)
        {
            GetQijFilePath(id, fragMatFile, qijpath, qijformat, isReadBinaryFile);
        }
        else 
        {
            GetFragMatFilePath(id, fragMatFile, frag_acc_path, fragaccformat, isReadBinaryFile);
        }
        if (GetBinaryMODM(fragMatFile, alphabet, length, profile, parameter , typeProfile) == -1)
        {
            fprintf(stderr, "can not open fragMatFile %s\n", fragMatFile);
            //continue;
        }
        for(i=0; i<lengthSeqQij; i++) { sumProfileFrag[i] = 0; } /*init sumProfileFrag*/
        for(i = 0 ; i < length; i ++)
        {
            aaSeqIndex = profile[i].aaSeqIndex-1;
            if (aaSeqIndex < 0)
            {
                fprintf(stderr,"Error aaSeqIndex (=%d) < 0, %s : %d  \n", aaSeqIndex, fragMatFile, i);
                assert(aaSeqIndex >= 0);
            }
            for (j=0; j<20; j++)
            {    
                matFrag[aaSeqIndex][j] =  int (profile[i].p[j]);
                sumProfileFrag[aaSeqIndex] +=  matFrag[aaSeqIndex][j];
            }
            matFragScore1[aaSeqIndex] = profile[i].score1; /*bug fixed, i changed to aaSeqIndex, 2008-02-18, Nanjiang*/
            matFragScore2[aaSeqIndex] = profile[i].score2;
        }

        /*read in normal MODM matrix file*/
        GetMODMFilePath_Tuping(id, modmFile, modmpath, modmformat, isReadBinaryFile);
        if (GetBinaryMODM(modmFile, alphabet, length, profile, parameter, typeProfile) == -1)
        {
            fprintf(stderr, "can not open modmFile %s\n", modmFile);
            return -1;
        }
        lengthSeqMODM = length;
        if (lengthSeqMODM != lengthSeqQij)
        {
            fprintf(stderr, "Warning! lengthSeqMODM(%d) not equal to lengthSeqQij(%d)\n", lengthSeqMODM , lengthSeqQij);
            fprintf(stderr, "%s: lengthSeqMODM = %d\n", modmFile , lengthSeqMODM);
            fprintf(stderr, "%s: lengthSeqQij = %d\n", QijFile, lengthSeqQij);
            /*assert(lengthSeqMODM == lengthSeqQij); */ /*2009-06-10, this will
            not affect the result, so no need to assertion, just give the warning sign*/
        }
        for(i = 0 ; i < lengthSeqMODM; i ++)
        {
            aaSeqIndex = profile[i].aaSeqIndex-1;
            if (aaSeqIndex < 0)
            {
                fprintf(stderr,"Error aaSeqIndex (=%d) < 0, %s : %d  \n", aaSeqIndex, modmFile, i);
                assert(aaSeqIndex >= 0);
            }
            for (j=0; j<20; j++) {matMODM[aaSeqIndex][j] = int(profile[i].p[j]); }
            matMODMScore1[aaSeqIndex] = profile[i].score1;
            matMODMScore2[aaSeqIndex] = profile[i].score2;

            aaSeq[aaSeqIndex]   = profile[i].aa;   /*amino acids*/
            shapeSeq[aaSeqIndex] = profile[i].shape; /*shape strings*/
            dsspSecSeq[aaSeqIndex] = DSSP8to3(profile[i].dsspSec, dsspMapMethod);
        }
        aaSeq[lengthSeqMODM] = '\0';
        shapeSeq[lengthSeqMODM] = '\0';
        dsspSecSeq[lengthSeqMODM] = '\0';
        TreatAllZeroFij(matMODM, lengthSeqMODM, matMODMScore1, aaSeq, AAAlphabet_Tuping);

    }/*}}}*/

    int daa = 0;
    for(i = 0 ; i < lengthSeqMODM; i++)
    {
        /*get digit amino acid 0-19*/
        daa = aaSeq[i] - 'A';
        if (  (daa>=0) && (daa<=25)  ) { daa = AAS_Code[daa]; }
        else { daa = 20; }
        digitAASeq[i] = daa;
    }
    if(mergeSide == type_dataset)/*{{{*/
    {
        /*profile combination of Frag- and psi-blast*/
        for(i = 0; i < lengthSeqMODM; i ++)
        {
            double ratio = 0.0;
            double profileMerged[NUM_20_AA+1];
            if ( sumProfileFrag[i] >= 10  ) /*if the frag profile is valid, isAllZeroProfile = false*/
            {
                if(ratioScheme == 0)/*ratioScheme = 0, merge fragacc and modm by specified ratio, NPer_Frag_Database*/
                {
                    ratio = NPer_Frag_Database/100.0;
                }
                else if(ratioScheme == 1) /*ratioScheme = 1, merge fragacc and modm according to the score2, 2008-01-18, Nanjiang*/
                {
                    ratio = GetMergeRatio1(i, matFragScore2, matMODMScore2, lengthSeqMODM);  /*bug fixed, matFragScore*/
                }
                else if (ratioScheme == 2)
                {
                    ratio = GetMergeRatio2(i, matFragScore2, matMODMScore2, lengthSeqMODM); 
                }
                else if (ratioScheme == 3)
                {
                    ratio = GetMergeRatio3(i, matFragScore2, matMODMScore2, lengthSeqMODM); 
                }

#ifdef DEBUG_MERGESIDE
                //if (LinearSearch_String (id,  testIDList_DEBUG_MERGESIDE, numID_DEBUG_MERGESIDE) != -1 )
                //{

                    fprintf(stdout,"id=%s, aa = %c, No = %d, mergeSide=%d, ratioScheme=%d, type_dataset=%d, ratio = %.3lf\n", id, aaSeq[i],i+1, mergeSide, ratioScheme, type_dataset, ratio);
                //}
#endif
                double sum = 0.0;
                for (j=0; j<20; j++)
                {
                    profileMerged[j] = ratio*matFrag[i][j]+(1.0-ratio)*matMODM[i][j];
                    sum += profileMerged[j];
                }
                for (j=0; j<20; j++)
                {
                    matMerged[i][j] = Integer(profileMerged[j]*100.0/sum);/*normalize the sum of the profile to 100, Nanjiang, 2008-01-15*/
                }
#ifdef DEBUG_MERGESIDE /*test the mergeside and read in matrix*/
                if (LinearSearch_String (id,  testIDList_DEBUG_MERGESIDE, numID_DEBUG_MERGESIDE) != -1 )
                {
                    fprintf(stdout,"Profile for chain %s, three lines, matFrag, matMODM, matMerged\n", id);
                    WriteMODMProfile(i+1, aaSeq[i], matFrag[i], matFragScore1[i],matFragScore2[i], AAAlphabet_Tuping, stdout);
                    WriteMODMProfile(i+1, aaSeq[i], matMODM[i], matMODMScore1[i],matMODMScore2[i], AAAlphabet_Tuping, stdout);
                    //                WriteMODMProfile(i+1, aaSeq[i], profileMerged, 0.11,0.11, AAAlphabet_Tuping, stdout);
                    WriteMODMProfile(i+1, aaSeq[i], matMerged[i], 0,0, AAAlphabet_Tuping, stdout);
                }
#endif
            }
            else //if all components are zero in frag profiles, use the Qij profiles
            {
                for (j=0; j<20; j++) { matMerged[i][j] = matQij[i][j]; }
            }
        }
    }/*}}}*/
    else /*do not merge the fragacc and qij on the training set, take the qij profile directly, 2008-01-18, Nanjiang*//*{{{*/
    {
#ifdef DEBUG_MERGESIDE
        if (LinearSearch_String (id,  testIDList_DEBUG_MERGESIDE, numID_DEBUG_MERGESIDE) != -1 )
        {
            fprintf(stdout,"id=%s, mergeSide=%d, type_dataset=%d, do not merge\n", id, mergeSide, type_dataset );
            fprintf(stdout,"Profile for chain %s, matQij\n", id);
            WriteMODMProfile(i+1, aaSeq[i], matQij[i], matQijScore1[i],matQijScore2[i], AAAlphabet_Tuping, stdout);
        }
#endif
        for(i = 0; i < lengthSeqMODM; i ++)
        {
            for(j = 0; j < NUM_20_AA; j ++) { matMerged[i][j] = matQij[i][j];}
        }
    }/*}}}*/
    return lengthSeqMODM;
}/*}}}*/
int ReadInDatabase_dumpedfile(const char *id, map<string,dbindex>&dbindexfragacc, map<string,dbindex>&dbindexqij, map<string,dbindex>&dbindexmodm, vector <FILE*>& fpList_fragacc, vector<FILE*>&fpList_qij, vector<FILE*>&fpList_modm, int **matFrag, int **matQij, int **matMODM, int **matMerged, float *matFragScore1, float *matFragScore2,  float *matQijScore1, float *matQijScore2, float *matMODMScore1, float *matMODMScore2, char *aaSeq, int *digitAASeq, char *shapeSeq, char *dsspSecSeq, int type_dataset, bool isReadBinaryFile, int NPer_Frag_Database, int ratioScheme, int mergeSide, int dsspMapMethod, int8 typeProfile)/*{{{*/
/*****************************************************************************
 * read in the dataset for a chain
 * including FragAcc matrix, Qij matrix and modm matrix
 * the type_dataset and isReadBinaryFile should be given
 * return the length of the chain if successful
 * return -1 for failure
 ****************************************************************************/
{
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

    int i,j;
    Array1D <int> sumProfileFrag_1darray(LONGEST_SEQ);
    int *sumProfileFrag = sumProfileFrag_1darray.array1D; /*sum of the profile at each position for fragacc matrix*/
    int lengthSeqMODM = 0;
    int lengthSeqQij = 0;

    FILE *fp_fragacc=NULL;
    FILE *fp_qij=NULL;
    FILE *fp_modm = NULL;

    int status_fseek = 0;
    fp_fragacc = fpList_fragacc[dbindexfragacc[id].dbfileindex];
    fp_qij = fpList_qij[dbindexqij[id].dbfileindex];
    fp_modm = fpList_modm[dbindexmodm[id].dbfileindex];

    if (!isReadBinaryFile) {/*{{{*/
        for (i=0; i < LONGEST_SEQ; i++) { sumProfileFrag[i] = 0; }
        if((ratioScheme == 0 && NPer_Frag_Database == 0) || mergeSide != type_dataset) {
            /*changed 2009-10-22, if mergeSide != type_dataset, also do not
             * need to read in FragMatFile, then replace it by qijfile if
             * ratioScheme == 0 and NPer_Frag_Database == 0, the Qij file is
             * used as fragacc file, so that the weight is 100% on Qij file*/
            if ((status_fseek = fseek(fp_qij, dbindexqij[id].offset, SEEK_SET)) != 0){
                fprintf(stderr,"Fatal! fseek failed for id %s in db qij. Exit.\n",id);
                exit(1);
            }
            ReadInProfile_fp(fp_qij, dbindexqij[id].size, matFrag, NULL, NULL, NULL, NULL, matFragScore1, matFragScore2, sumProfileFrag, LONGEST_SEQ);
        } else {
            if ((status_fseek = fseek(fp_fragacc, dbindexfragacc[id].offset, SEEK_SET)) != 0){
                fprintf(stderr,"Fatal! fseek failed for id %s in db fragacc. Exit.\n",id);
                exit(1);
            }
            ReadInProfile_fp(fp_fragacc, dbindexfragacc[id].size, matFrag, NULL, NULL, NULL, NULL, matFragScore1, matFragScore2, sumProfileFrag, LONGEST_SEQ);
        }

        //open Qij profiles for reading rows whose sum of components are zero
        if ((status_fseek = fseek(fp_qij, dbindexqij[id].offset, SEEK_SET)) != 0){
            fprintf(stderr,"Fatal! fseek failed for id %s in db qij. Exit.\n",id);
            exit(1);
        }
        if ((lengthSeqQij = ReadInProfile_fp(fp_qij, dbindexqij[id].size, matQij, NULL, NULL, NULL, NULL, matQijScore1, matQijScore2, NULL, LONGEST_SEQ)) < 0) {
            return -1 ; 
        }

        //open normal psi-blast percentage profiles, modm matrix, in tuping's format
        if ((status_fseek = fseek(fp_modm, dbindexmodm[id].offset, SEEK_SET)) != 0){
            fprintf(stderr,"Fatal! fseek failed for id %s in db modm. Exit.\n",id);
            exit(1);
        }
        if ((lengthSeqMODM = ReadInProfile_fp(fp_modm, dbindexmodm[id].size, matMODM, aaSeq,shapeSeq, NULL, dsspSecSeq, matMODMScore1, matMODMScore2, NULL, LONGEST_SEQ)) < 0){
             return -1 ;
        } /*}}}*/
    } else { /*isReadBinaryFile == true*/ /*{{{*/
        Array1D <ProfileSADByte> profile_1darray(LONGEST_SEQ);
        ProfileSADByte *profile = profile_1darray.array1D;
        char alphabet[50+1] = "";
        int length = 0; /*number of profiles read in*/
        double parameter[8];
        int aaSeqIndex = 0;

        /*2007-11-11, Nanjiang, changed to read in Qij first and then if fragacc exist, replace it with fragacc*/
        /*read in Qij matrix file, for those residue positions does not have fragacc profile, using the Qij profile instead*/
        if ((status_fseek = fseek(fp_qij, dbindexqij[id].offset, SEEK_SET)) != 0){
            fprintf(stderr,"Fatal! fseek failed for id %s in db modm. Exit.\n",id);
            exit(1);
        }
        if (GetBinaryMODM_fp(fp_qij, dbindexqij[id].size, alphabet, length, profile, parameter, typeProfile) < 0 ) {
            fprintf(stderr, "Reading qij matrix failed for %s. offset=%ld, size=%ld. __LINE__=%d\n", id, dbindexqij[id].offset, dbindexqij[id].size, __LINE__);
            return -1;
        }
        lengthSeqQij = length ;
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

        /*read in frag acc matrix file, replace the Psi_blosum_Frag matrix with fragacc profiles */
        if ((ratioScheme == 0 && NPer_Frag_Database == 0)|| mergeSide != type_dataset )
            /*changed 2009-10-22, if mergeSide != type_dataset, also do not
             * need to read in FragMatFile, then replace it by qijfile*/
        {
            if ((status_fseek = fseek(fp_qij, dbindexqij[id].offset, SEEK_SET)) != 0){
                fprintf(stderr,"Fatal! fseek failed for id %s in db qij. Exit.\n",id);
                exit(1);
            }
            if (GetBinaryMODM_fp(fp_qij, dbindexqij[id].size, alphabet, length, profile, parameter, typeProfile ) < 0) {
                fprintf(stderr, "Reading qij matrix failed for %s. offset=%ld, size=%ld. __LINE__=%d\n", id, dbindexqij[id].offset, dbindexqij[id].size, __LINE__);
                return -1;
            }
        } else {
            if ((status_fseek = fseek(fp_fragacc, dbindexfragacc[id].offset, SEEK_SET)) != 0){
                fprintf(stderr,"Fatal! fseek failed for id %s in db fragacc. Exit.\n",id);
                exit(1);
            }
            if (GetBinaryMODM_fp(fp_fragacc, dbindexfragacc[id].size, alphabet, length, profile, parameter, typeProfile ) < 0 ) {
                fprintf(stderr, "Reading fragacc matrix failed for %s. offset=%ld, size=%ld\n", id, dbindexfragacc[id].offset, dbindexfragacc[id].size);
                return -1;
            }
        }
        for(i=0; i<lengthSeqQij; i++) { sumProfileFrag[i] = 0; } /*init sumProfileFrag*/
        for(i = 0 ; i < length; i ++) {
            aaSeqIndex = profile[i].aaSeqIndex-1;
            if (aaSeqIndex < 0) {
                fprintf(stderr,"Error aaSeqIndex (=%d) < 0, fragacc %s : %d\n", aaSeqIndex, id, i);
                assert(aaSeqIndex >= 0);
            }
            for (j=0; j<20; j++) {    
                matFrag[aaSeqIndex][j] =  int (profile[i].p[j]);
                sumProfileFrag[aaSeqIndex] +=  matFrag[aaSeqIndex][j];
            }
            matFragScore1[aaSeqIndex] = profile[i].score1; /*bug fixed, i changed to aaSeqIndex, 2008-02-18, Nanjiang*/
            matFragScore2[aaSeqIndex] = profile[i].score2;
        }

        /*read in normal MODM matrix file*/
        if ((status_fseek = fseek(fp_modm, dbindexmodm[id].offset, SEEK_SET)) != 0){
            fprintf(stderr,"Fatal! fseek failed for id %s in db modm. Exit.\n",id);
            exit(1);
        }
        if (GetBinaryMODM_fp(fp_modm, dbindexmodm[id].size, alphabet, length, profile, parameter, typeProfile) < 0 ) {
            fprintf(stderr, "Reading modm matrix failed for %s. offset=%ld, size=%ld\n", id, dbindexmodm[id].offset, dbindexmodm[id].size);
            return -1;
        }

        lengthSeqMODM = length;
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
    }/*}}}*/
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
                if (LinearSearch_String (id,  testIDList_DEBUG_MERGESIDE, numID_DEBUG_MERGESIDE) != -1 ) {

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
                if (LinearSearch_String (id,  testIDList_DEBUG_MERGESIDE, numID_DEBUG_MERGESIDE) != -1 ) {
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
        if (LinearSearch_String (id,  testIDList_DEBUG_MERGESIDE, numID_DEBUG_MERGESIDE) != -1 ) {
            fprintf(stdout,"id=%s, mergeSide=%d, type_dataset=%d, do not merge\n", id, mergeSide, type_dataset );
            fprintf(stdout,"Profile for chain %s, matQij\n", id);
            WriteMODMProfile(i+1, aaSeq[i], matQij[i], matQijScore1[i],matQijScore2[i], AAAlphabet_Tuping, stdout);
        }
#endif
        for(i = 0; i < lengthSeqMODM; i ++) {
            for(j = 0; j < NUM_20_AA; j ++) { matMerged[i][j] = matQij[i][j];}
        }
    }/*}}}*/
    return lengthSeqMODM;
}/*}}}*/

template <class T> int ReadSMatrix(const char *filename, T **S, char *alphabet)/*{{{*/
    /*****************************************************************************
     * read in substitution matrix for protein sequence alignment, 
     * e.g., BLOSUM62, which is also the default substitution matrix
     ****************************************************************************/
{
    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D; 

    char *pch;
    int i,j;
    FILE* fpout = NULL;
    fpout = fopen(filename,"r");
    checkfilestream(fpout,filename, "r");

    while(fgetline(fpout, line,maxline) != EOF)
    {
        if(line[0] != '#') break ;
    }
    SpanExcluding(line,alphabet," ");
    int n = strlen(alphabet);
    char delim[] = WHITE_SPACE;
    i = 0;
    while((linesize = fgetline(fpout, line,maxline))!= EOF)
    {
        pch = strtok(line, delim);
        j = 0;
        while(pch != NULL)
        {
            if(j >= 1 && j < n+1) S[i][j-1] = T (atof(pch));
            j ++;
            pch = strtok(NULL,delim);
        }
        i ++;
    }
    if (i != n)
    {
        fprintf(stderr,"Error! The number of rows and columns in the matrix %s are not equal.\n",filename);
        assert(i == n);
    }
    fclose(fpout);
    return n;
}
template int ReadSMatrix<int>   (const char* filename, int **S   , char *alphabet);
template int ReadSMatrix<float> (const char* filename, float **S , char *alphabet);
template int ReadSMatrix<double>(const char* filename, double **S, char *alphabet);
/*}}}*/


int GetAltAtomIndex(Atom *atoms, int n)/*{{{*/
    // return the index of selected Alternative atom according to the maximal
    // occupancy
{
    if(n == 1 ) return 0;
    int i ;
    int idx = 0;
    float maxOcc;
    maxOcc = MIN_FLOAT;
    for(i = 0 ; i < n ; i++)
    {
        if(atoms[i].occupancy > maxOcc)
        {
            maxOcc = atoms[i].occupancy;
            idx = i ;
        }
    }
    return idx;
}/*}}}*/
int GetAASeqIndex(const char* id, int resSeq, char resICode, Chain* pChain, int startIndex /*= 0*/)/*{{{*/
{
    int i;
    if( strncmp( id, pChain->pdbid, 4) != 0 || pChain-> chainID != id[4])
        return -1;
    else
    {
        for(i = startIndex ; i < pChain->numRes; i ++)
        {
            if( pChain->resSer[i] == resSeq && pChain->resICode[i] == resICode)
            {
                return i;
            }
        }
    }
    return (-1); //if not found in pChain
}/*}}}*/
int GetResidueIndex(Residue *pRes, Residue *resGroup, int numRes)/*{{{*/
    /*****************************************************************************
     * get the residue index from a group of residues, 
     * return the index of residue if found 
     * and return -1 if the residue does not exist in the resGroup
     * 2007-04-13,Nanjiang Shu
     ****************************************************************************/
{
    int i;
    if(numRes == 0 ) return -1;
    for ( i = 0 ; i < numRes; i++)
    {
        if(   pRes->chainID  == resGroup[i].chainID
                && pRes->resSeq   == resGroup[i].resSeq
                && pRes->resICode == resGroup[i].resICode )
            return i;
    }
    return -1;
}/*}}}*/
int GetResidueIndex(Atom *pAtom, Residue *resGroup, int numRes)/*{{{*/
    /*******************************************************************
     * Get the residue index for the residue the atom belongs, from a group of residues
     * return the index of the residue if found 
     * return -1 if the residue does not exist in the resGroup
     *******************************************************************/
{
    int i;
    if(numRes == 0 ) return -1;
    for ( i = 0 ; i < numRes; i++)
    {
        if(   pAtom->chainID == resGroup[i].chainID
                && pAtom->resSeq  == resGroup[i].resSeq 
                && pAtom->iCode   == resGroup[i].resICode )
            return i;
    }
    return -1;
}/*}}}*/

double GetContactDist(const char* eleName)/*{{{*/
{
    const char *metalNameList1[] = {"CD","CU", "CO", "FE", "MG", "NI","MN","ZN"};
    const char *metalNameList2[] = {"BA","K","NA"};
    const char *metalNameList3[] = {"CS"};
    int numList1 = sizeof(metalNameList1) / sizeof(char*);
    int numList2 = sizeof(metalNameList2) / sizeof(char*);
    int numList3 = sizeof(metalNameList3) / sizeof(char*);

    if(BinarySearch_String(eleName,metalNameList1,numList1) != -1 )
        return 2.9;
    else if( BinarySearch_String(eleName,metalNameList2, numList2) != -1)
        return 3.25;
    else if( BinarySearch_String(eleName,metalNameList3, numList3) != -1)
        return 3.5;
    else
        return 3.0;
}
/*}}}*/

int GetZnBoundRes(const char* metalProFile, MetalPro* znPro, int min_numRes /*= 3*/ , int max_numRes /*= 4*/, bool isUsingTotalBoundRes /*= true*/)/*{{{*/
    /*****************************************************************************
     *  get zinc binding residues for each zinc binding protein, the znPro is of
     *  structure MetalPro, that is all the zinc binding residues are stored
     *  directly under znPro, the redundent residues, i.e., those residues
     *  binds to more than one zinc ions, are removed 
     *  
     *  for GetZnBoundRes, res.parentAtomEnvIndex is not recored,  2007-04-11
     *
     *  isUsingTotalBoundRes, whether to use atomEnv.totalBoundRes or
     *  atomEnv.numRes for evaluating min_numRes and max_numRes rule, 2007-04-19
     *
     *  metalPro.atomEnv.totoalBoundRes is recorded
     *  but the allocated memory for metalPro.res is = metalPro.numBoundRes * sizeof(Residue)
     ****************************************************************************/
{

    char keyMetal[SIZE_ATOM_ELEMENT+1] = "ZN";
    // int  min_numRes  = 3; /* minimal binding resNum such that the binding site is counted*/ 
    // int  max_numRes  = 4;

    using namespace ZnBindingProtein;

    FILE* fpin; 
    fpin = fopen(metalProFile,"r");
    checkfilestream(fpin, metalProFile,"r");

    int  i, j;
    int  status_sscanf;
    int  index = 0;
    char id[SIZE_CHAIN_ID+1] = "";
    int  length = 0;
    int  numMetalAtom = 0;
    char metalAtomName[SIZE_METAL_ATOM_NAME+1] = "";
    int  metalAtomResSeq = 0;
    char metalAtomResName[SIZE_METAL_ATOM_RES_NAME+1] = "";
    char metalAtomChainID = 0;
    char chainID = ' ';
    int  numRes = 0;
    int  totalBoundRes = 0; // total number of residues bind to a metal atom, including residues not in the chain 
    int  maxline = 500;
    int  linesize;

    Array1D <char>    line_1darray(maxline+1);
    Array1D <Residue> tmpRes_1darray(MAX_BOUND_SITE);
    Array1D <int>     idx_1darray(MAX_BOUND_SITE);
    char    *line   = line_1darray.array1D;
    Residue *tmpRes = tmpRes_1darray.array1D;
    int     *idx    = idx_1darray.array1D;

    MetalPro tmpZnPro;
    InitMetalPro(&tmpZnPro);
    tmpZnPro.res = new Residue[MAX_NUMRES];
    for(i = 0 ; i < MAX_NUMRES; i++) 
    {
        InitResidue(&tmpZnPro.res[i]);
        //tmpZnPro.res[i].parentAtomEnvIndex = new int[MAX_BOUND_METAL_PER_RES];
    }

    tmpZnPro.resSeqIndex = new int[MAX_NUMRES];
    tmpZnPro.atomEnv = new AtomEnv[MAX_ATOMENV];
    for(i = 0 ; i < MAX_ATOMENV; i++)
    {
        InitAtomEnv(&tmpZnPro.atomEnv[i]);
        tmpZnPro.atomEnv[i].parentResIndex = new int[MAX_BOUND_SITE];
    }


    //set <int> idxPool;

    Residue *pRes;
    //AtomEnv *pAtomEnv;
    int cntZnPro = 0;
    bool isHaveZn;

    f_neglect_comment(fpin, '#');
    while((linesize = fgetline(fpin,line,maxline))!= EOF)
    {
        if(linesize <= 0) continue;
        if((status_sscanf = sscanf(line,"%d %5s %d %d",&index,id,&length,&numMetalAtom))< 4 ) 
        { 
            fprintf(stderr,"line =%s\n",line);
            assert(status_sscanf == 4);
        }

        isHaveZn = false;
        StdID(id);         //standardize nrPDB chain identifier
        chainID = id[4]; /* chainID always refers to the chainID of protein chain, which is equal to
                            id[4] for the standardized chain identifier*/ 

        //idxPool.clear();

        tmpZnPro.length     = length;
        tmpZnPro.numBoundRes = 0;
        my_strcpy(tmpZnPro.id,id,SIZE_CHAIN_ID);


        int cntNumRes = 0;
        int cntZnAtom = 0;
        if(numMetalAtom > 0 )
        {
            for(i = 0 ; i < numMetalAtom ; i++)
            {
                linesize = fgetline(fpin, line, maxline);
                if((status_sscanf = sscanf(line,"%2s %d :%1c %d %s",metalAtomName, &metalAtomResSeq, &metalAtomChainID,&totalBoundRes,metalAtomResName)) < 5 )
                {
                    fprintf(stderr,"line =%s\n", line);
                    assert(status_sscanf == 5);
                }
                if( totalBoundRes > 0)
                {
                    ScanfCloseMetalRes(fpin,tmpRes,totalBoundRes);//read in residues bound to this metal atom


                    // get the numRes for that chain
                    int cntRes = 0;
                    for(j = 0; j < totalBoundRes ; j++)
                    {
                        if(tmpRes[j].chainID == chainID)
                        {
                            idx[cntRes] = j; cntRes++;
                        }
                    }
                    numRes = cntRes; // numRes record the number of residues (binding to the meta atom) on a single chain,

                    if( strcasecmp(metalAtomName, keyMetal ) == 0 
                            && ( (!isUsingTotalBoundRes && numRes >= min_numRes && numRes <= max_numRes) 
                                ||(isUsingTotalBoundRes && totalBoundRes >= min_numRes && totalBoundRes <= max_numRes)))
                    {
                        isHaveZn = true;
                        tmpZnPro.atomEnv[cntZnAtom].totalBoundRes =totalBoundRes;
                        tmpZnPro.atomEnv[cntZnAtom].numRes = numRes;
                        my_strcpy(tmpZnPro.atomEnv[cntZnAtom].metalAtomName, metalAtomName, SIZE_METAL_ATOM_NAME);
                        my_strcpy(tmpZnPro.atomEnv[cntZnAtom].metalAtomResName, metalAtomResName, SIZE_METAL_ATOM_RES_NAME);
                        tmpZnPro.atomEnv[cntZnAtom].metalAtomChainID = metalAtomChainID;
                        my_strcpy(tmpZnPro.atomEnv[cntZnAtom].metalAtomPDBID,id,SIZE_PDBID);
                        my_strcpy(tmpZnPro.atomEnv[cntZnAtom].id,id,SIZE_CHAIN_ID);
                        tmpZnPro.atomEnv[cntZnAtom].seqLength = length;

                        for ( j = 0; j < numRes ; j++)
                        {
                            pRes = & (tmpRes[idx[j]]);
                            int resIndex = GetResidueIndex(pRes, tmpZnPro.res, cntNumRes);
                            if(resIndex == -1) //not found, new residue
                            {
                                CopyResidue( & (tmpZnPro.res[cntNumRes]), pRes);
                                tmpZnPro.resSeqIndex[cntNumRes] = pRes->aaSeqIndex;   //note that resSeqIndex might not be unique if isExcludeOtherChainRes = false
                                tmpZnPro.atomEnv[cntZnAtom].parentResIndex[j] = cntNumRes;
                                cntNumRes ++;
                            }
                            else
                            {
                                tmpZnPro.atomEnv[cntZnAtom].parentResIndex[j] = resIndex;
                            }

                            //int *pSeqIndex;
                            //pSeqIndex = find (tmpZnPro.resSeqIndex, tmpZnPro.resSeqIndex+cntNumRes, pRes->aaSeqIndex);
                            //if(idxPool.find(pRes -> aaSeqIndex) == idxPool.end()) [>not found,  new residue<] 
                            //if(pSeqIndex == tmpZnPro.resSeqIndex+cntNumRes)[> not found, new residue<] 
                            //{
                            //CopyResidue( & (tmpZnPro.res[cntNumRes]), pRes);
                            //idxPool.insert(pRes->aaSeqIndex);
                            //tmpZnPro.resSeqIndex[cntNumRes] = pRes->aaSeqIndex;
                            //tmpZnPro.atomEnv[cntZnAtom].parentResIndex[j] = cntNumRes;
                            //cntNumRes ++;
                            //}
                            //else
                            //{
                            //tmpZnPro.atomEnv[cntZnAtom].parentResIndex[j] = pSeqIndex - tmpZnPro.resSeqIndex;
                            //}
                        }
                        cntZnAtom ++;
                    }
                }
            }
            ////mapping atomEnv to the parentAtomEnvIndex under res
            //for ( i = 0 ; i < cntZnAtom; i++)
            //{
            //    pAtomEnv = &(tmpZnPro.atomEnv[i]);
            //    for(j = 0 ; j < pAtomEnv -> numRes; j ++)
            //    {
            //        int idxRes = pAtomEnv->parentResIndex[j];
            //        int cntMetalBound = tmpZnPro.res[idxRes].numMetalBound;
            //        tmpZnPro.res[idxRes].parentAtomEnvIndex[cntMetalBound] = i;
            //        tmpZnPro.res[idxRes].numMetalBound ++;
            //    }
            //}
        }
        tmpZnPro.numBoundRes = cntNumRes;
        tmpZnPro.numMetalAtom = cntZnAtom;
        if(isHaveZn)
        {
            znPro[cntZnPro].res = new Residue[tmpZnPro.numBoundRes];
            for( j = 0; j < tmpZnPro.numBoundRes; j++ )
            {
                InitResidue(&(znPro[cntZnPro].res[j]));
                //metalPro[cntZnPro].res[j].parentAtomEnvIndex = new int[tmpZnPro.res[j].numMetalBound];
            }
            znPro[cntZnPro].resSeqIndex = new int[tmpZnPro.numBoundRes];
            znPro[cntZnPro].atomEnv = new AtomEnv[tmpZnPro.numMetalAtom];
            for( j = 0; j < tmpZnPro.numMetalAtom; j++ ) InitAtomEnv(&(znPro[cntZnPro].atomEnv[j]));
            for( j = 0 ; j < tmpZnPro.numMetalAtom ; j ++)
            {
                znPro[cntZnPro].atomEnv[j].parentResIndex = new int[tmpZnPro.atomEnv[j].numRes];
            }
            CopyMetalPro(&znPro[cntZnPro], &tmpZnPro,true);
            cntZnPro ++;
        }
    }

    int numZnPro = cntZnPro;

    fclose(fpin);
    /*free memory*/
    DeleteMetalPro(&tmpZnPro, MAX_ATOMENV, MAX_NUMRES);

    return numZnPro;
}/*}}}*/
int GetZnBoundRes(const char* metalProFile, MetalPro2* znPro, int min_numRes /*= 3*/, int max_numRes /*= 4*/, bool isUsingTotalBoundRes /* = true*/)/*{{{*/
    /*****************************************************************************
     * override of GetZnBoundRes, 
     * here znPro use struc MetalPro2, in which residues are stored according to
     * the metal atom. In that case, residues binding to the same metal atom can be
     * identified
     ****************************************************************************/
{
    char keyMetal[SIZE_ATOM_ELEMENT+1] = "ZN";
    // int  min_numRes  = 3; /* minimal binding resNum such that the binding site is counted*/ 
    // int  max_numRes  = 4;

    // get the bindind residues for each zinc binding proteins
    using namespace ZnBindingProtein;
    char comment_char = '#';

    FILE* fpin; 
    fpin = fopen(metalProFile,"r");
    checkfilestream(fpin, metalProFile,"r");

    int  i, j;
    int status_sscanf;
    int  index = 0;
    char id[SIZE_CHAIN_ID+1] = "";
    char chainID;
    int  length;
    int  numMetalAtom;
    char metalAtomName[SIZE_METAL_ATOM_NAME+1] = "";
    int  metalAtomResSeq = 0;
    char metalAtomResName[SIZE_METAL_ATOM_RES_NAME+1] = "";
    char metalAtomChainID;
    int  totalBoundRes = 0;
    int  numRes = 0;

    int linesize;
    int maxline = 500;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;

    Array1D <Residue> tmpRes_array1D(LONGEST_SEQ);
    Residue *tmpRes = tmpRes_array1D.array1D;
    for(i = 0 ; i < LONGEST_SEQ ; i ++) InitResidue(&tmpRes[i]);
    Array1D <int> idx_1darray(LONGEST_SEQ);
    int  *idx       = idx_1darray.array1D;

    MetalPro2 tmpZnPro;
    InitMetalPro(&tmpZnPro);
    tmpZnPro.atomEnv = new AtomEnv[MAX_ATOMENV];
    for(i = 0 ; i < MAX_ATOMENV ; i++) InitAtomEnv(&(tmpZnPro.atomEnv[i]));
    for(i = 0; i < MAX_ATOMENV; i++)
    {
        tmpZnPro.atomEnv[i].res = new Residue[MAX_BOUND_SITE];
        for(j = 0 ; j < MAX_BOUND_SITE; j++) InitResidue(&(tmpZnPro.atomEnv[i].res[j]));
    }

    int numZnPro;
    int cntZnPro;
    bool isHaveZn;
    int numZnMetalAtom;

    f_neglect_comment(fpin, comment_char);
    cntZnPro = 0;
    while((linesize = fgetline(fpin,line,maxline))!= EOF)
    {

        index      = 0;
        isHaveZn   = false;
        numZnMetalAtom = 0;
        if(linesize <= 0 ) continue;
        if((status_sscanf = sscanf(line,"%d %5s %d %d",&index,id,&length,&numMetalAtom))< 4 ) 
        { 
            fprintf(stderr,"line =%s\n", line);
            assert(status_sscanf == 4);
        }
        if (index == 0) break;

        StdID(id);         //standardize nrPDB chain identifier
        chainID = id[4];

        tmpZnPro.length   = length;
        tmpZnPro.numMetalAtom = 0;
        my_strcpy(tmpZnPro.id,id,SIZE_CHAIN_ID);

        if( numMetalAtom <= 0) continue;

        for(i = 0 ; i < numMetalAtom ; i++)
        {
            //fscanf(fpin,"%2s:%1c %d %3s\n",metalAtomName,&metalAtomChainID,&numRes,metalAtomResName);
            linesize = fgetline(fpin, line, maxline);
            if((status_sscanf = sscanf(line,"%2s %d :%1c %d %s",metalAtomName,&metalAtomResSeq, &metalAtomChainID,&totalBoundRes,metalAtomResName)) < 4 )
            {
                fprintf(stderr, "line =%s\n", line);
                assert(status_sscanf == 5);
            }
            if( numRes > 0)
            {
                ScanfCloseMetalRes(fpin,tmpRes,totalBoundRes);
                // get the numRes for that chain
                int cntRes = 0;
                for(j = 0; j < totalBoundRes ; j++)
                {
                    if(tmpRes[j].chainID == chainID)
                    {
                        idx[cntRes] = j; cntRes++;
                    }
                }
                numRes = cntRes; // numRes record the number of residues (binding to the meta atom) on a single chain,

                if( strcasecmp(metalAtomName, keyMetal ) == 0 
                        &&((!isUsingTotalBoundRes && numRes >= min_numRes && numRes <= max_numRes)
                            ||(isUsingTotalBoundRes && totalBoundRes >= min_numRes && totalBoundRes <= max_numRes)))
                {
                    // printf("%d   id=%s\n",cntZnPro, id); //debug 
                    isHaveZn = true;
                    tmpZnPro.atomEnv[numZnMetalAtom].totalBoundRes = totalBoundRes;
                    tmpZnPro.atomEnv[numZnMetalAtom].numRes = numRes;
                    my_strcpy( tmpZnPro.atomEnv[numZnMetalAtom].metalAtomName, keyMetal,SIZE_METAL_ATOM_NAME);
                    for ( j = 0; j < numRes ; j++)
                    {
                        CopyResidue( & (tmpZnPro.atomEnv[numZnMetalAtom].res[j]), &tmpRes[idx[j]]);
                    }
                    numZnMetalAtom ++;
                }
            }
        }

        tmpZnPro.numMetalAtom = numZnMetalAtom;
        // allocate memory for znPro, 
        if( isHaveZn)
        {
#ifdef _ASSERT_
            //debug start
            if(znPro[cntZnPro].atomEnv != NULL)
            {
                printf("id=%s, cntZnPro = %d\n",id,cntZnPro);
                assert( znPro[cntZnPro].atomEnv == NULL);
            }
            //debug end
#endif            
            znPro[cntZnPro].atomEnv = new AtomEnv[numZnMetalAtom];
            for( i = 0 ; i < numZnMetalAtom ; i++)
            {    
                znPro[cntZnPro].atomEnv[i].res = new Residue[tmpZnPro.atomEnv[i].numRes];
                for(j = 0 ; j < tmpZnPro.atomEnv[i].numRes; j++) InitResidue(&(znPro[cntZnPro].atomEnv[i].res[j]));
            }


            // copy tmpZnPro to znPro[cntZnPro]
            CopyMetalPro(& znPro[cntZnPro], &tmpZnPro,true);
            cntZnPro ++;
        }
    }
    fclose(fpin);
    numZnPro = cntZnPro;

    /*free memory*/
    for(i = 0 ; i < LONGEST_SEQ; i++)
    {
        DeleteResidue(&tmpRes[i]);
    }

    DeleteMetalPro(&tmpZnPro,MAX_ATOMENV);

    return numZnPro;
}
/*}}}*/
int GetZnBoundRes3(const char* metalProFile, MetalPro* znPro, MetalPro* znPro1, MetalPro* znPro2, int &numZnPro, int &numZnPro1, int &numZnPro2, int &numZnRes, int &numZnRes1, int &numZnRes2, double cutoff_score2, char *resList, const char* modmpath, int min_ZnBoundRes /*= 3*/, int max_ZnBoundRes /*= 4*/, bool isUsingTotalBoundRes /*= true*/)/*{{{*/
    // get the zinc binding residues for each zinc binding protein, Zn3 or Zn4
    // znPro   -- store zinc binding residues in (Zn3 and Zn4) within resList and filtered by cutoff_score2
    // znPro1  -- store zinc binding residues in (Zn3 and Zn4) within resList
    // znPro2  -- store all zinc binding residues in (Zn3 and Zn4)
{
    using namespace ZnBindingProtein;
    int i, j, k;
    char modmfilepath[MAX_PATH+1] = "";

    MetalPro tmpZnPro;
    InitMetalPro(&tmpZnPro);
    tmpZnPro.res = new Residue[MAX_NUMRES];
    for( i = 0 ; i < MAX_NUMRES ; i++) 
    {    
        InitResidue(&(tmpZnPro.res[i]));
        //tmpZnPro.res[i].parentAtomEnvIndex = new int[MAX_BOUND_METAL_PER_RES];
    }

    tmpZnPro.resSeqIndex = new int[MAX_NUMRES];
    tmpZnPro.atomEnv = new AtomEnv[MAX_ATOMENV];
    for( i= 0 ; i < MAX_ATOMENV; i++)
    {
        InitAtomEnv(&(tmpZnPro.atomEnv[i]));
        tmpZnPro.atomEnv[i].parentResIndex = new int[MAX_BOUND_SITE];
    }


    AtomEnv *pAtomEnv;
    Residue *pRes;

    Residue *pRes1;
    Residue *pRes2;
    int cntZnRes;
    int cntZnPro;

    numZnPro2 = GetZnBoundRes(metalProFile, znPro2, min_ZnBoundRes, max_ZnBoundRes, isUsingTotalBoundRes); 
    /* numZnRes2 might not equal to the summation of all numRes over all metal
     * binding site, since there are cases one residue binding to more than one
     * metal ion */ 

    cntZnRes = 0;
    for(i = 0; i < numZnPro2; i++) 
        cntZnRes += znPro2[i].numBoundRes;
    numZnRes2 = cntZnRes;

    /*get the zinc binding protein information for znPro and znPro1*/
    cntZnPro = 0;
    cntZnRes = 0;
    for( i = 0; i < numZnPro2 ; i++)
    {
        CopyMetalPro( &tmpZnPro, &znPro2[i],false);
        int cntBoundRes = 0 ;

        for(j = 0 ; j < znPro2[i].numBoundRes; j++)
        {
            pRes1 = & (tmpZnPro.res[cntBoundRes]);
            pRes2 = & (znPro2[i].res[j]);
            if(IsInCharSet(pRes2->aa,resList))
            {
                CopyResidue(pRes1,pRes2);
                tmpZnPro.resSeqIndex[cntBoundRes] = pRes2->aaSeqIndex;
                cntBoundRes ++;
                cntZnRes ++;
            }
        }
        tmpZnPro.numBoundRes = cntBoundRes;

        int cntZnAtom = 0;
        for(j = 0 ; j < znPro2[i].numMetalAtom; j++)
        {
            pAtomEnv = &(znPro2[i].atomEnv[j]);
            int cntNumRes = 0;
            for(k = 0 ; k < znPro2[i].atomEnv[j].numRes; k ++)
            {
                pRes = &(znPro2[i].res[pAtomEnv->parentResIndex[k]]);
                int resIndex = GetResidueIndex(pRes, tmpZnPro.res, tmpZnPro.numBoundRes);
                if(resIndex != -1) //if this residue is in the updated tmpMetalPro.res, add it, 2007-05-04
                {
                    tmpZnPro.atomEnv[cntZnAtom].parentResIndex[cntNumRes] = resIndex;
                    cntNumRes ++;
                }
                //int* pSeqIndex = find (tmpZnPro.resSeqIndex, tmpZnPro.resSeqIndex+cntBoundRes, znPro2[i].resSeqIndex[znPro2[i].resSeqIndex[znPro2[i].atomEnv[j].parentResIndex[k]]]);
                //if(pSeqIndex != tmpZnPro.resSeqIndex+cntBoundRes) // if found
                //{    
                //tmpZnPro.atomEnv[cntZnAtom].parentResIndex[cntNumRes] = pSeqIndex-tmpZnPro.resSeqIndex;
                //cntNumRes ++;
                //}
            }
            tmpZnPro.atomEnv[cntZnAtom].numRes = cntNumRes;
            if (cntNumRes > 0)
                cntZnAtom ++;
        }
        tmpZnPro.numMetalAtom = cntZnAtom;

        if(cntBoundRes > 0)
        {
            znPro1[cntZnPro].res = new Residue[tmpZnPro.numBoundRes];
            for(j = 0 ; j < tmpZnPro.numBoundRes; j++) 
            {
                InitResidue(&(znPro1[cntZnPro].res[j]));
                //znPro1[cntZnPro].res[j].parentAtomEnvIndex = new int[tmpZnPro.res[j].numMetalBound];
            }

            znPro1[cntZnPro].resSeqIndex = new int[tmpZnPro.numBoundRes];
            znPro1[cntZnPro].atomEnv = new AtomEnv[tmpZnPro.numMetalAtom];
            for(j = 0 ; j < tmpZnPro.numMetalAtom; j++) InitAtomEnv(&(znPro1[cntZnPro].atomEnv[j]));
            for( j = 0 ; j < tmpZnPro.numMetalAtom ; j ++)
            {
                znPro1[cntZnPro].atomEnv[j].parentResIndex = new int[tmpZnPro.atomEnv[j].numRes];
            }
            CopyMetalPro(&znPro1[cntZnPro], &tmpZnPro);/* copy tmpZnPro to znPro[numZnPro]*/
            cntZnPro ++;
        }
    } 
    numZnPro1 = cntZnPro;
    numZnRes1 = cntZnRes;

    /*get znPro, which is restricted by both resList and cutoff_score2*/
    cntZnPro = 0;
    cntZnRes = 0;
    char id[SIZE_CHAIN_ID+1] = "";

    Array1D <double> score1_1darray(LONGEST_SEQ);
    Array1D <double> score2_1darray(LONGEST_SEQ);
    Array1D <char>   aaSeq_1darray(LONGEST_SEQ+1);
    Array2D <DATATYPE_MODM_MATRIX>    M_2darray(LONGEST_SEQ,        NUM_BLOSUM);
    double  *score1 = score1_1darray.array1D;
    double  *score2 = score2_1darray.array1D;
    char    *aaSeq  = aaSeq_1darray.array1D;
    DATATYPE_MODM_MATRIX    **M      = M_2darray.array2D;
    int length;
    char alphabetMODM[NUM_BLOSUM+1] = "";

    for( i = 0; i < numZnPro2 ; i++)
    {
        my_strcpy(id, znPro2[i].id, SIZE_CHAIN_ID); id[SIZE_CHAIN_ID] = '\0';
        if(GetMODMFilePath(id,modmfilepath, modmpath) != NULL)
            length = GetMODM(modmfilepath, M, alphabetMODM, aaSeq, score1, score2);
        else
            continue;


        CopyMetalPro( &tmpZnPro, &znPro2[i],false);
        int cntBoundRes = 0;
        for(j = 0 ; j < znPro2[i].numBoundRes; j++)
        {
            pRes1 = & (tmpZnPro.res[cntBoundRes]);
            pRes2 = & (znPro2[i].res[j]);
            if(IsInCharSet(pRes2->aa,resList) && score2[pRes2->aaSeqIndex] >= cutoff_score2)
            {
                CopyResidue(pRes1,pRes2);
                tmpZnPro.resSeqIndex[cntBoundRes] = pRes2->aaSeqIndex;
                cntBoundRes ++;
                cntZnRes ++;
            }
        }
        tmpZnPro.numBoundRes = cntBoundRes;

        int cntZnAtom = 0;
        for(j = 0 ; j < znPro2[i].numMetalAtom; j++)
        {
            pAtomEnv = &(znPro2[i].atomEnv[j]);
            int cntNumRes = 0;
            for(k = 0 ; k < znPro2[i].atomEnv[j].numRes; k ++)
            {
                pRes = & (znPro2[i].res[pAtomEnv->parentResIndex[k]]);
                int resIndex = GetResidueIndex(pRes, tmpZnPro.res, tmpZnPro.numBoundRes);
                if(resIndex != -1) //if this residue is in the updated tmpMetalPro.res, add it, 2007-05-04
                {
                    tmpZnPro.atomEnv[cntZnAtom].parentResIndex[cntNumRes] = resIndex;
                    cntNumRes ++;
                }

                //int* pSeqIndex = find (tmpZnPro.resSeqIndex, tmpZnPro.resSeqIndex+cntBoundRes, znPro2[i].resSeqIndex[znPro2[i].atomEnv[j].parentResIndex[k]]);
                //if(pSeqIndex != tmpZnPro.resSeqIndex+cntBoundRes) // if found
                //{    
                //tmpZnPro.atomEnv[cntZnAtom].parentResIndex[cntNumRes] = pSeqIndex-tmpZnPro.resSeqIndex;
                //cntNumRes ++;
                //}
            }
            tmpZnPro.atomEnv[cntZnAtom].numRes = cntNumRes;
            if (cntNumRes > 0)
                cntZnAtom ++;
        }
        tmpZnPro.numMetalAtom = cntZnAtom;

        if(cntBoundRes > 0)
        {
            znPro[cntZnPro].res = new Residue[tmpZnPro.numBoundRes];
            for(j = 0 ; j < tmpZnPro.numBoundRes; j++) InitResidue(&(znPro[cntZnPro].res[j]));
            znPro[cntZnPro].resSeqIndex = new int[tmpZnPro.numBoundRes];
            znPro[cntZnPro].atomEnv = new AtomEnv[tmpZnPro.numMetalAtom];
            for(j = 0 ; j < tmpZnPro.numMetalAtom; j++) InitAtomEnv(&(znPro[cntZnPro].atomEnv[j]));
            for( j = 0 ; j < tmpZnPro.numMetalAtom ; j ++)
            {
                znPro[cntZnPro].atomEnv[j].parentResIndex = new int[tmpZnPro.atomEnv[j].numRes];
            }
            CopyMetalPro(&znPro[cntZnPro], &tmpZnPro);/* copy tmpZnPro to znPro[numZnPro]*/
            cntZnPro ++;
        }
    }       /* znPro and znPro1 got*/ 
    numZnPro = cntZnPro;
    numZnRes = cntZnRes;

    DeleteMetalPro(&tmpZnPro, MAX_ATOMENV, MAX_NUMRES);

    return numZnPro2;
}
/*}}}*/
int GetZnBoundRes3(const char* metalProFile, MetalPro2* znPro, MetalPro2* znPro1, MetalPro2* znPro2, int &numZnPro, int &numZnPro1, int &numZnPro2,int &numZnRes, int &numZnRes1, int &numZnRes2, double cutoff_score2, char *resList, const char* modmpath, int min_ZnBoundRes /*= 3*/, int max_ZnBoundRes /*= 4*/, bool                    isUsingTotalBoundRes /*= true*/)/*{{{*/
{
    using namespace ZnBindingProtein;

    int i, j,k;
    int cntZnRes;
    int cntZnPro;

    Residue *pRes1;
    Residue *pRes2;
    char modmfilepath[MAX_PATH+1];
    int numResList = strlen(resList);

    MetalPro2 tmpZnPro;
    InitMetalPro(&tmpZnPro);
    tmpZnPro.atomEnv = new AtomEnv[MAX_ATOMENV];
    {
        for(i = 0 ; i < MAX_ATOMENV ; i++)
            InitAtomEnv(&(tmpZnPro.atomEnv[i]));
    }
    for(i = 0; i < MAX_ATOMENV; i++)
    {
        tmpZnPro.atomEnv[i].res = new Residue[MAX_BOUND_SITE];
        for(j = 0 ; j < MAX_BOUND_SITE; j++) InitResidue(&(tmpZnPro.atomEnv[i].res[j]));
    }

    numZnPro2 = GetZnBoundRes(metalProFile, znPro2, min_ZnBoundRes, max_ZnBoundRes); 
    /* numZnRes2 might not equal to the summation of all numRes over all metal
     * binding site, since there are cases one residue binding to more than one
     * metal ion */ 

    set <int> idxPool;
    cntZnRes = 0;
    for(i = 0; i < numZnPro2; i++)
    {
        idxPool.clear();
        for(j = 0; j < znPro2[i].numMetalAtom;j++)
        {
            for(k = 0; k < znPro2[i].atomEnv[j].numRes; k++)
                idxPool.insert(znPro2[i].atomEnv[j].res[k].aaSeqIndex);
        }
        cntZnRes += idxPool.size();
    }
    numZnRes2 = cntZnRes;

    /*get the zinc binding protein information for znPro and znPro1*/
    cntZnPro  = 0;
    cntZnRes  = 0;
    for( i = 0; i < numZnPro2 ; i++)
    {
        int numAtomEnv;
        int numRes;
        numAtomEnv  = 0;
        idxPool.clear();

        CopyMetalPro( &tmpZnPro, &znPro2[i],false);
        for(j = 0 ; j < znPro2[i].numMetalAtom; j++)
        {
            numRes = 0;
            for(k = 0 ; k < znPro2[i].atomEnv[j].numRes; k++)
            {
                pRes1 = &(tmpZnPro.atomEnv[numAtomEnv].res[numRes]);
                pRes2 = &(znPro2[i].atomEnv[j].res[k]);
                if(IsInCharSet(pRes2->aa,resList,numResList))
                {
                    CopyResidue(pRes1,pRes2);
                    idxPool.insert(pRes2->aaSeqIndex);
                    numRes ++;
                }
            }
            if(numRes > 0)
            {
#ifdef _ASSERT_
                if(&(tmpZnPro.atomEnv[numAtomEnv])== NULL)
                {
                    printf("ASSERT: numAtomEnv =%d\n",numAtomEnv);
                    assert( &(tmpZnPro.atomEnv[numAtomEnv]) != NULL);
                }
#endif
                CopyAtomEnv(&(tmpZnPro.atomEnv[numAtomEnv]),&(znPro2[i].atomEnv[j]),false); 
                tmpZnPro.atomEnv[numAtomEnv].numRes = numRes;
                // copy the znPro parameters without residue
                numAtomEnv ++;
            }
        }
        tmpZnPro.numMetalAtom = numAtomEnv;

        if(numAtomEnv > 0)
        {
            znPro1[cntZnPro].atomEnv = new AtomEnv[numAtomEnv];
            for( j = 0 ; j < numAtomEnv ; j++)
            {
                znPro1[cntZnPro].atomEnv[j].res = new Residue[tmpZnPro.atomEnv[j].numRes];
                for(k = 0 ; k < tmpZnPro.atomEnv[j].numRes ; k++) InitResidue(&(znPro1[cntZnPro].atomEnv[j].res[k]));
            }

            // copy tmpZnPro to znPro[numZnPro]
            CopyMetalPro(&znPro1[cntZnPro], &tmpZnPro);
            cntZnRes += idxPool.size();
            cntZnPro ++;
        }
    } 
    numZnPro1 = cntZnPro;
    numZnRes1 = cntZnRes;


    // getting znPro

    Array1D <double> score1_1darray(LONGEST_SEQ);
    Array1D <double> score2_1darray(LONGEST_SEQ);
    Array1D <char>   aaSeq_1darray(LONGEST_SEQ+1);
    Array2D <DATATYPE_MODM_MATRIX>    M_2darray(LONGEST_SEQ,        NUM_BLOSUM);
    double  *score1 = score1_1darray.array1D;
    double  *score2 = score2_1darray.array1D;
    char    *aaSeq  = aaSeq_1darray.array1D;
    DATATYPE_MODM_MATRIX    **M      = M_2darray.array2D;
    //char id[SIZE_CHAIN_ID+1];
    int length;
    char alphabetMODM[MAX_PATH+1];

    cntZnRes = 0;
    cntZnPro = 0;
    for( i = 0; i < numZnPro2 ; i++)
    {
        int numAtomEnv;
        int numRes;
        if(GetMODMFilePath(znPro2[i].id,modmfilepath, modmpath) != NULL)
        {
            length = GetMODM(modmfilepath, M, alphabetMODM, aaSeq, score1, score2);
        } else{
            continue;
        }

        numAtomEnv  = 0;
        idxPool.clear();

        CopyMetalPro( &tmpZnPro, &znPro2[i],false);
        for(j = 0 ; j < znPro2[i].numMetalAtom; j++)
        {
            numRes  = 0;
            for(k = 0 ; k < znPro2[i].atomEnv[j].numRes; k++)
            {
                pRes1 = &(tmpZnPro.atomEnv[numAtomEnv].res[numRes]);
                pRes2 = &(znPro2[i].atomEnv[j].res[k]);
                if(IsInCharSet(pRes2->aa,resList,numResList) && score2[pRes2->aaSeqIndex] >= cutoff_score2)
                {
                    CopyResidue(pRes1,pRes2);
                    idxPool.insert(pRes2->aaSeqIndex);
                    numRes ++;
                }
            }
            if(numRes > 0)
            {
#ifdef _ASSERT_
                if(&(tmpZnPro.atomEnv[numAtomEnv])== NULL)
                {
                    printf("ASSERT: numAtomEnv =%d\n",numAtomEnv);
                    assert( &(tmpZnPro.atomEnv[numAtomEnv]) != NULL);
                }
#endif                
                CopyAtomEnv(&(tmpZnPro.atomEnv[numAtomEnv]),&(znPro2[i].atomEnv[j]),false); 

                tmpZnPro.atomEnv[numAtomEnv].numRes = numRes;
                // copy the znPro parameters without residue
                numAtomEnv ++;
            }
        }
        tmpZnPro.numMetalAtom = numAtomEnv;

        if(numAtomEnv > 0)
        {
            znPro[cntZnPro].atomEnv = new AtomEnv[numAtomEnv];
            for(j = 0 ; j < numAtomEnv; j++) InitAtomEnv(&(znPro[cntZnPro].atomEnv[j]));
            for( j = 0 ; j < numAtomEnv ; j++)
            {
                znPro[cntZnPro].atomEnv[j].res = new Residue[tmpZnPro.atomEnv[j].numRes];
                for(k = 0 ; k < tmpZnPro.atomEnv[j].numRes; k++) InitResidue(&(znPro[cntZnPro].atomEnv[j].res[k]));
            }

            // copy tmpZnPro to znPro[numZnPro]
            CopyMetalPro(&znPro[cntZnPro], &tmpZnPro);
            cntZnRes += idxPool.size();
            cntZnPro ++;
        }
    } /* znPro, znPro1, znPro2 got*/ 
    numZnRes = cntZnRes;
    numZnPro = cntZnPro;

    DeleteMetalPro(&tmpZnPro,MAX_ATOMENV);

    return numZnPro2;
}
/*}}}*/
int GetSSBondRes(const char* ssbondProFile, SSBondPro *ssbondPro, int &numSSBondPro, int &numSSBondRes)/*{{{*/
{
    int i;
    int maxline = 200;
    char comment_char = '#';
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;

    FILE *fp = fopen(ssbondProFile,"r");
    checkfilestream(fp, ssbondProFile, "r");

    SSBondPro tmpssbondPro;
    InitSSBondPro(&tmpssbondPro);
    tmpssbondPro.ssbond = new SSBond[MAX_SSBOND_CHAIN];
    for(i = 0 ; i < MAX_SSBOND_CHAIN; i++)
        InitSSBond(&(tmpssbondPro.ssbond[i])) ;

    //SSBond *pssbond;
    int cntSSBondPro = 0;
    int cntSSBondRes = 0;
    //char delim[20+1] = ",";
    int index ;
    int length;
    int numSSBond;
    char id[SIZE_CHAIN_ID+1];

    f_neglect_comment(fp, comment_char);

    cntSSBondPro = 0;
    cntSSBondRes = 0;
    while(fgetline ( fp, line ,maxline) != EOF)
    {
        if( sscanf(line, "%d %5s %d %d", &index, id, &length, &numSSBond) < 4) continue;
        StdID(id);

        tmpssbondPro.length    = length;
        tmpssbondPro.numSSBond = numSSBond;
        my_strcpy(tmpssbondPro.id,id, SIZE_CHAIN_ID) ;

        ScanfSSBondRecord(fp, &tmpssbondPro, numSSBond);

        if(numSSBond > 0)
        {
            ssbondPro[cntSSBondPro].ssbond = new SSBond[numSSBond];
            for(i = 0 ; i < numSSBond; i++) InitSSBond(&(ssbondPro[cntSSBondPro].ssbond[i]));
            CopySSBondPro(&(ssbondPro[cntSSBondPro]), &tmpssbondPro, true);
            cntSSBondPro ++;
            cntSSBondRes += tmpssbondPro.numSSBondRes;
        }
    }
    fclose(fp);
    numSSBondPro = cntSSBondPro;
    numSSBondRes = cntSSBondRes;

    DeleteSSBondPro(&tmpssbondPro);

    return numSSBondPro;
}
/*}}}*/
int GetSSBondRes2(const char *ssbondProFile, SSBondPro *ssbondPro, SSBondPro *ssbondPro1,int &numSSBondPro, int &numSSBondPro1, int &numSSBondRes, int &numSSBondRes1, double cutoff_score2, const char* modmpath)/*{{{*/
{
    int i,j,k;

    char modmfilepath[MAX_PATH+1] = "";

    SSBondPro tmpssbondPro;
    InitSSBondPro(&tmpssbondPro);
    tmpssbondPro.ssbond = new SSBond[MAX_SSBOND_CHAIN];
    for(i = 0 ; i < MAX_SSBOND_CHAIN; i++)
        InitSSBond(&(tmpssbondPro.ssbond[i])) ;

    int cntSSBondPro = 0;
    int cntSSBondRes = 0;
    numSSBondPro = GetSSBondRes(ssbondProFile, ssbondPro1, numSSBondPro1, numSSBondRes1);
    //ssbondPro1 is including all ssbonded residues on the chain
    //ssbondPro is with those cysteins with score2 < 0.1 removed

    char id[SIZE_CHAIN_ID+1];

    Array1D <double> score1_1darray(LONGEST_SEQ);
    Array1D <double> score2_1darray(LONGEST_SEQ);
    Array1D <char>   aaSeq_1darray(LONGEST_SEQ+1);
    Array2D <DATATYPE_MODM_MATRIX>    M_2darray(LONGEST_SEQ,        NUM_BLOSUM);
    double  *score1 = score1_1darray.array1D;
    double  *score2 = score2_1darray.array1D;
    char    *aaSeq  = aaSeq_1darray.array1D;
    DATATYPE_MODM_MATRIX    **M      = M_2darray.array2D;
    int length;
    char alphabetMODM[NUM_BLOSUM+1] = "";
    Residue *pRes1;
    Residue *pRes2;
    int aaSeqIndex;


    cntSSBondPro = 0; 
    cntSSBondRes = 0 ;
    for ( i = 0 ; i < numSSBondPro1 ; i++)
    {
        my_strcpy(id, ssbondPro1[i].id, SIZE_CHAIN_ID); 
        if(GetMODMFilePath(id,modmfilepath, modmpath) != NULL)
            length = GetMODM(modmfilepath, M, alphabetMODM, aaSeq, score1, score2);
        else
            continue;
        CopySSBondPro( &tmpssbondPro, &ssbondPro1[i],false);
        int cntSSBond = 0;
        tmpssbondPro.numSSBondRes = 0;
        for(j = 0 ; j < ssbondPro1[i].numSSBond; j++)
        {
            CopySSBond(&(tmpssbondPro.ssbond[cntSSBond]), &(ssbondPro1[i].ssbond[j]), false);
            int cntNumRes = 0;
            for(k = 0 ; k < ssbondPro1[i].ssbond[j].numRes; k++)
            {
                pRes1 = &(tmpssbondPro.ssbond[cntSSBond].res[cntNumRes]);
                pRes2 = &(ssbondPro1[i].ssbond[j].res[k]);
                aaSeqIndex = pRes2->aaSeqIndex;
                if(score2[aaSeqIndex] >= cutoff_score2)
                {
                    CopyResidue(pRes1,pRes2);
                    cntNumRes ++ ;
                }
            }
            if(cntNumRes > 0)
            {
                tmpssbondPro.ssbond[cntSSBond].numRes = cntNumRes;
                cntSSBond ++;
                tmpssbondPro.numSSBondRes += cntNumRes;
            }
        }
        tmpssbondPro.numSSBond = cntSSBond;

        if(cntSSBond > 0)
        {
            ssbondPro[cntSSBondPro].ssbond = new SSBond[cntSSBond];
            for(j = 0 ; j < cntSSBond; j ++) InitSSBond(&(ssbondPro[cntSSBondPro].ssbond[j]));
            CopySSBondPro(&ssbondPro[cntSSBondPro], &tmpssbondPro);/* copy tmpZnPro to znPro[numZnPro]*/
            cntSSBondPro ++;
            cntSSBondRes += tmpssbondPro.numSSBondRes;
        }
    }       /* znPro and znPro1 got*/ 

    numSSBondPro = cntSSBondPro;
    numSSBondRes = cntSSBondRes;

    // free memory
    DeleteSSBondPro(&tmpssbondPro);

    return numSSBondPro;
}
/*}}}*/

int GetSVMVector(int aaSeqIndex, MODM *pMODM, int K, int P,double* vector, double cutoff_score2 /*= 0.1*/)/*{{{*/
    /*****************************************************************************
     * GetSVMVector for single residue
     * 
     * return size of vector of successful
     * return -1 , if the failed, for example, if the last column of profile is of
     * too low value
     ****************************************************************************/
{
    // P = 24
    // they are 
    // profile 20
    // last two columns 2
    // flag indicating where overflow sequence 1 { -1, overflow, 1 not overflow }
    // hydrophobicity 1         
    double W1 = 1.0; /* weight for percentage */ 
    double W2 = 1.0; /* weight for score1*/ 
    double W3 = 1.0; /* weight for score2*/ 
    double W4 = 1.0; /* weight for hydrophobicity*/ 
    double W5 = 1.0;  /* weight for digit AA*/  

    DATATYPE_MODM_MATRIX    **M   = pMODM -> M;
    double  *score1 = pMODM -> score1;
    double  *score2 = pMODM -> score2;
    char    *aaSeq  = pMODM -> aaSeq;
    int      length = pMODM -> length;

    int i,j;
    // double hydro; // hydrophobicity
    int shift;
    int beg;

    // check the profile 
    for ( j = aaSeqIndex -K ; j <= aaSeqIndex + K ;j ++)
    {
        if(j >= 0 && j < length)
        {
            if( score2[j] < cutoff_score2)
                return -1;
        }
    }

    shift = aaSeqIndex-K;
    for ( j = aaSeqIndex -K ; j <= aaSeqIndex +K ;j ++)
    {
        beg = (j-shift)*P;
        if(j >= 0 && j < length)
        {
            for( i = beg ; i <beg +20; i ++)
            {
                // vector[i] = (double)M[j][i-beg] /100.0;
                vector[i] = (double)M[j][i-beg] * W1;
            }
            vector[i] = (double)score1[j] * W2;
            i++; vector[i] = (double)score2[j] * W3; 
            i++; vector[i] = 1.0;       
            i++; vector[i] = Hydrophobicity[Char2Digit(aaSeq[j],STD1CharAA_alphabet)] * W4;
            if(P > 24)
            {
                i++;
                vector[i] = (double)(Char2Digit(aaSeq[j],STD1CharAA_alphabet)+1)* W5;
            }
        }
        else
        {
            for( i = beg; i < beg +20 +2; i ++)
            {
                vector[i] = 0.0;
            }
            vector[i] = -1.0;
            i++; vector[i] = 0.0;
            if(P > 24)
            {
                i++; vector[i] = 0.0;
            }
        }

    }
    return (2*K+1)*P;
}/*}}}*/
int GetSVMVector(int aaSeqIndex1, int aaSeqIndex2, MODM *pMODM, int K, int W,int P, double *vector, double cutoff_score2 /*= 0.1*/)/*{{{*/
    /*****************************************************************************
     * GetSVMVector for resiude pair
     * 
     * return size of vector of successful
     * return -1 , if the failed, for example, if the last column of profile is of
     * too low value
     ****************************************************************************/
{

    // -----k----1----internal---1---k------
    // if internel > 2*W, take only the start and end W

    // P = 24
    // they are 
    // profile 20
    // last two columns 2
    // flag indicating where overflow sequence 1 { -1, overflow, 1 not overflow }
    // hydrophobicity 1         

    double W1 = 1.0;   // weight for percentage profile
    double W2 = 1.0;   // weight for score1
    double W3 = 1.0;   // weight for score2
    double W4 = 1.0;   // weight for hydrophobicity
    double W5 = 1.0;  /* weight for digit AA*/  

    DATATYPE_MODM_MATRIX    **M   = pMODM -> M;
    double  *score1 = pMODM -> score1;
    double  *score2 = pMODM -> score2;
    char    *aaSeq  = pMODM -> aaSeq;
    int      length = pMODM -> length;

    int i,j;
    // double hydro; // hydrophobicity
    int shift;
    int beg;
    int inter;
    Array1D <double> profile_1darray(P);
    double* profile = profile_1darray.array1D;

    // check the profile 
    for ( j = aaSeqIndex1 -K ; j <= aaSeqIndex2 +K ;j ++)
    {
        if(j >= 0 && j < length)
        {
            if( score2[j] < cutoff_score2)
                return -1 ;
        }
    }

    shift = aaSeqIndex1 -K;
    for ( j = aaSeqIndex1 -K ; j <= aaSeqIndex1 ;j ++)
    {
        beg = (j-shift)*P;
        if(j >= 0 && j < length)
        {
            for( i = beg ; i < beg +20; i ++)
            {
                // vector[i] = double(M[j][i-beg]) / 100.0;
                vector[i] = double(M[j][i-beg]) * W1;
            }
            vector[i] = score1[j] * W2 ;
            i++; vector[i] = score2[j] * W3; 
            i++; vector[i] = 1.0;      
            i++; vector[i] = Hydrophobicity[Char2Digit(aaSeq[j],STD1CharAA_alphabet)] * W4;
            if(P > 24)
            {
                i++;
                vector[i] = (double)(Char2Digit(aaSeq[j],STD1CharAA_alphabet)+1)* W5;
            }
        }
        else
        {
            for( i = beg ; i < beg +20 +2; i ++)
            {
                vector[i] = 0.0;
            }
            vector[i] = -1.0;
            i++; vector[i] = 0.0;
            if(P > 24)
            {
                i++; vector[i] = 0.0;
            }
        }
    }

    shift = aaSeqIndex2 -(K+2);
    for ( j = aaSeqIndex2; j <= aaSeqIndex2+K ;j ++)
    {
        beg = (j-shift)*P;
        if(j >= 0 && j < length)
        {
            for( i = beg; i < beg+20; i ++)
            {
                vector[i] = double(M[j][i-beg]) * W1;
            }
            vector[i] = score1[j];
            i++;vector[i] = score2[j];
            i++;vector[i] = 1.0;
            i++; vector[i] = Hydrophobicity[Char2Digit(aaSeq[j],STD1CharAA_alphabet)];
            if(P > 24)
            {
                i++;
                vector[i] = (double)(Char2Digit(aaSeq[j],STD1CharAA_alphabet)+1)* W5;
            }
        }
        else
        {
            for( i = beg ; i < beg+20 +2; i ++)
            {
                vector[i] = 0.0;
            }
            vector[i] = -1.0;
            i++;vector[i] = 0.0;
            if(P > 24)
            {
                i++; vector[i] = 0.0;
            }
        }
    }

    inter = aaSeqIndex2 - aaSeqIndex1 -1 ;
    beg =(K+1)*P;
    if( inter < 1)
    {
        for( i = beg ; i < beg+20 +2; i ++)
        {
            vector[i] = 0.0;
        }
        vector[i] = -1.0;
        i++;  vector[i] = 0.0;
        if(P > 24)
        {
            i++; vector[i] = 0.0;
        }
    }
    else
    {
        for(i = 0 ; i < P ; i++) profile[i] = 0.0;
        for(j = aaSeqIndex1+1; j < aaSeqIndex2; j++)
        {
            if(j - aaSeqIndex1 <=W || aaSeqIndex2 - j <=W)
            {
                for(i = 0 ; i < 20 ; i++) 
                    profile[i] += M[j][i] *W1;
                profile[i] += score1[j] * W2;
                i++; profile[i] += score2[j] * W3;
                i++;profile[i] += 1.0;      
                i++; profile[i] +=  Hydrophobicity[Char2Digit(aaSeq[j],STD1CharAA_alphabet)] *W4; 
                if(P > 24)
                {
                    i++;
                    profile[i] += (double)(Char2Digit(aaSeq[j],STD1CharAA_alphabet)+1)* W5;
                }
            }
        }
        //for(i = 0 ; i < P ; i++) profile[i] /= double(inter);
        for(i = beg ; i< beg +P ; i++) 
            vector[i] = profile[i-beg];
    }

    return (2*K+1+2)*P; 
}/*}}}*/
int GetSVMVectorPerSite(double *vector, int vector_beg, int seqPos, int P, MODM *pMODM, double W1, double W2, double W3, double W4, double W5, double cutoff_score2, int encoding_type /* = PSSM_PROFILE_ENCODE*/)/*{{{*/
/*****************************************************************************
 *  length : length of the protein sequence
 *  seqPos : seqPos is the aaSeqIndex for that position
 *  vector_beg: vector_beg is the beginning index for the vector
 *  changeLog 2007-06-26, within the sequence 1 , without the sequence 0
 ****************************************************************************/
{
    int i = 0;
    int j = seqPos;
    int beg = vector_beg;
    if(j >= 0 && j < pMODM->length && pMODM->score2[j] >= cutoff_score2)
    {
        for( i = beg ; i < beg +20; i ++)
        {
            // vector[i] = double(M[j][i-beg]) / 100.0;
            if(encoding_type == PSSM_PROFILE_ENCODE)
                vector[i] = double(pMODM->M[j][i-beg]) * W1;
            else if (encoding_type == SEQUENCE_BINARY_ENCODE)
                vector[i] = IS_EQUAL(pMODM->aaSeq[j], STD1CharAA_alphabet[i]);
            else 
                vector[i] = 0.0;
        }
        if(P >= 21) { vector[i] = pMODM->score1[j] * W2 ;}
        if(P >= 22) {i++; vector[i] = pMODM->score2[j] * W3;}
        if(P >= 23) {i++; vector[i] = 1.0; /*within the seqeunce*/    }
        if(P >= 24) { i++; vector[i] = Hydrophobicity[Char2Digit(pMODM->aaSeq[j],STD1CharAA_alphabet)] * W4; }
        if(P >= 25) { i++; vector[i] = (double)(Char2Digit(pMODM->aaSeq[j],STD1CharAA_alphabet)+1)* W5; }
    }
    else  /*out of the range of the sequence*/ 
    {
        for( i = beg ; i < beg +20 ; i ++)
        {
            vector[i] = 0.0;
        }
        if(P >= 21) { vector[i] = 0.0; }
        if(P >= 22) { i++; vector[i] = 0.0; }
        if(P >= 23) { i++; vector[i] = 0.0;  /*out of the range of the sequence, endcode as 0*/ }
        if(P >= 24) { i++; vector[i] = 0.0; }
        if(P >= 25) { i++; vector[i] = 0.0; }
    }

    if( i-beg != P-1)
    {

        fprintf(stderr,"vector size conflicts, i = %d, P = %d\n", i, P);
        assert (i-beg == (P-1));//check if not all planed positions of the vector are filled
    }
    return P;
}
/*}}}*/
int GetSVMVectorPerSite_2(double *vector, int vector_beg, int seqPos, int P, int dimVectorPerSite, MODM *pMODM, double w1, double w2, double w3, double w4, double w5, double cutoff_score2, bool isCenterResidue, int encoding_type /* = PSSM_PROFILE_ENCODE*/)/*{{{*/
/*****************************************************************************
 *  length : length of the protein sequence
 *  seqPos : seqPos is the aaSeqIndex for that position
 *  vector_beg: vector_beg is the beginning index for the vector
 *  return last position assigned in the vector, that is the vector_beg pos
 *  for the nect assignement
 ****************************************************************************/
{
    int i = 0;   // i is iterator for vector position
    int ik = 0;

    int j = seqPos;
    int beg = vector_beg;
    if(j >= 0 && j < pMODM->length && pMODM->score2[j] >= cutoff_score2)/*{{{*/
        // j is iterator for residue position
    {
        i = beg;
        if(P >= 1)
        {
            for( ik = 0 ; ik < SIZE_ENCODE_PSSM; ik ++)
            {
                if(encoding_type == PSSM_PROFILE_ENCODE)
                    vector[i+ik] = double(pMODM->M[j][ik]) * w1;
                else if (encoding_type == SEQUENCE_BINARY_ENCODE)
                    vector[i+ik] = IS_EQUAL(pMODM->aaSeq[j], STD1CharAA_alphabet[ik]);
                else 
                    vector[i+ik] = 0.0;
            }
            i += SIZE_ENCODE_PSSM;
        }
        if(P >= 2) 
        { 
            int idx = GetEncodeScore1(pMODM->score1[j]);
            for (ik = 0 ; ik < SIZE_ENCODE_SCORE1; ik ++)
            {
                vector[i+ik] = EncodeMatrixScore1[idx][ik]* w2 ;
            }
            i += SIZE_ENCODE_SCORE1;
        }
        if(P >= 3) 
        { 
            int idx = GetEncodeScore2(pMODM->score2[j]);
            for (ik = 0 ; ik < SIZE_ENCODE_SCORE2; ik ++)
            {
                vector[i+ik] = EncodeMatrixScore2[idx][ik]* w3 ;
            }
            i += SIZE_ENCODE_SCORE2;
        }
        if(P >= 4) 
        { 
            vector[i] = 1.0; /*within the seqeunce*/ 
            i ++;   
        }
        if(P >= 5) 
        { 
            if(isCenterResidue)
            {
                //int idx = GetEncodeAA(pMODM->aaSeq[j]);
                int daa = Char2Digit(pMODM->aaSeq[j], EncodeMatrixAA_Center_alphabet);
                if(daa < 0) { daa = SIZE_ENCODE_AA_CENTER-1;}
                for (ik = 0 ; ik < SIZE_ENCODE_AA_CENTER; ik ++)
                {
                    vector[i+ik] = EncodeMatrixAA_Center[daa][ik]* w4 ;
                }
                i += SIZE_ENCODE_AA_CENTER;
            }
            else
            {
                int daa = Char2Digit(pMODM->aaSeq[j], EncodeMatrixAA_Non_Center_alphabet);
                if(daa < 0) { daa = SIZE_ENCODE_AA_NON_CENTER-1;}
                for (ik = 0 ; ik < SIZE_ENCODE_AA_NON_CENTER; ik ++)
                {
                    vector[i+ik] = EncodeMatrixAA_Non_Center[daa][ik]* w4 ;
                }
                i += SIZE_ENCODE_AA_NON_CENTER;
            }
        }
        if(P >= 6) 
        { 
            int daa = Char2Digit(pMODM->aaSeq[j], STD1CharAA_alphabet);
            if  (daa < 0) { daa = SIZE_ENCODE_HYDROPHOBICITY-1;}
            int idx = GetEncodeHydrophobicity(Hydrophobicity[daa]);
            for (ik = 0 ; ik < SIZE_ENCODE_HYDROPHOBICITY; ik ++)
            {
                vector[i+ik] = EncodeMatrixHydrophobicity[idx][ik]* w5 ;
            }
            i += SIZE_ENCODE_HYDROPHOBICITY;
        }
    }/*}}}*/
    else  /*out of the range of the sequence*/ /*{{{*/
    {
        for( i = beg; i < beg + dimVectorPerSite; i ++)
        {
            vector[i] = 0.0;
        }
    }/*}}}*/

    if( i-beg != dimVectorPerSite)
    {
        fprintf(stderr,"vector size conflicts, i = %d, P = %d, dimVectorPerSite = %d\n", i, P, dimVectorPerSite);
        assert (i-beg == (dimVectorPerSite-1));//check if not all planed positions of the vector are filled
    }
    return (beg+dimVectorPerSite) ;
}
/*}}}*/
int GetDimensionVectorPerSite(int P, bool isCenterResidue)/*{{{*/
/*****************************************************************************
 * return the dimension of the vector per site
 ****************************************************************************/
{
    int dimVectorPerSite = 0;
    dimVectorPerSite = SIZE_ENCODE_PSSM     * bool (max(0, P-0)) +
                       SIZE_ENCODE_SCORE1   * bool (max(0, P-1)) +
                       SIZE_ENCODE_SCORE2   * bool (max(0, P-2)) +
                       1                    * bool (max(0, P-3)) +
                       isCenterResidue    * SIZE_ENCODE_AA_CENTER       * bool (max(0, P-4)) +
                       (!isCenterResidue) * SIZE_ENCODE_AA_NON_CENTER   * bool (max(0, P-4)) +
                       SIZE_ENCODE_HYDROPHOBICITY * bool (max(0, P-5));

    return dimVectorPerSite;
}/*}}}*/
int GetSVMVector(int *aaSeqIndex, int numSite, MODM *pMODM, int K, int W, int P, double *vector, int vectorDim, double cutoff_score2 /*= 0.1*/, int encoding_type /*= PSSM_PROFILE_ENCODE*/)/*{{{*/
/****************************************************************************
 * vector dimension can be 
 * (2*K + numSite + numSite-1)*P or
 * (2*K + numSite + numSite-1)*P + 1
 * the last one indicates the distance separating the residues for a
 * neighbouring pair
 * return the size of the vector
 * return value <= 0 means failure
 * 
 * Encoding_type: PSSM_PROFILE_ENCODE         PSSM based encoding , 2007-04-27
 *                SEQUENCE_BINARY_ENCODE      sequence based binary encoding
 **************************************************************************/
{
    // -----k----1--w--internal-w--1---w---internal--w-1-k
    // if internel > 2*W, take only the start and end W

    // P = 24
    // they are 
    // profile 20
    // last two columns 2
    // flag indicating where overflow sequence 1 { -1, overflow, 1 not overflow }
    // hydrophobicity 1         

    double W1 = 1.0;   // weight for percentage profile
    double W2 = 1.0;   // weight for score1
    double W3 = 1.0;   // weight for score2
    double W4 = 1.0;   // weight for hydrophobicity
    double W5 = 1.0;  /* weight for digit AA*/  

    DATATYPE_MODM_MATRIX    **M      = pMODM -> M;
    double  *score1 = pMODM -> score1;
    double  *score2 = pMODM -> score2;
    char    *aaSeq  = pMODM -> aaSeq;
    //int      length = pMODM -> length;

    int i,j;
    int s;
    // double hydro; // hydrophobicity
    int shift;
    int beg;
    //int end;
    int inter;
    Array1D <double> profile_1darray(P);
    double* profile = profile_1darray.array1D;

    //// check the profile 
    //beg = aaSeqIndex[0]-K;
    //end = aaSeqIndex[numSite-1]+K;
    //for ( j = beg ; j <= end ;j ++)
    //{
    //    if(j >= 0 && j < length)
    //    {
    //        if( score2[j] < cutoff_score2)
    //            return -1 ;
    //    }
    //}
    //modified 2007-08-08, do not neglect  the vector if there is only one
    //residue with score2 < cutoff_score2, but set the value of the vecotr to
    //zero


    shift = aaSeqIndex[0] -K;
    for ( j = aaSeqIndex[0] -K ; j <= aaSeqIndex[0] ;j ++)
    {
        beg = (j-shift)*P;
        GetSVMVectorPerSite(vector, beg, j, P, pMODM, W1, W2, W3, W4, W5, cutoff_score2, encoding_type);
    }

    shift = aaSeqIndex[numSite-1] -(K+2*(numSite-1));
    for ( j = aaSeqIndex[numSite-1]; j <= aaSeqIndex[numSite-1]+K ;j ++)
    {
        beg = (j-shift)*P;
        GetSVMVectorPerSite(vector, beg,j, P, pMODM, W1, W2, W3, W4, W5, cutoff_score2, encoding_type);
    }

    for(s = 0 ; s < numSite-1 ; s++)
    {
        inter = aaSeqIndex[s+1] - aaSeqIndex[s] -1 ;
        beg =(K + s*2+1)*P;
        if( inter < 1)
        {
            for( i = beg ; i < beg+20; i ++) { vector[i] = 0.0; }
            if(P >= 21) { vector[i] = 0.0;}
            if(P >= 22) { i++; vector[i] = 0.0;}
            if(P >= 23) { i++; vector[i] = 0.0;}
            if(P >= 24) { i++; vector[i] = 0.0; }
            if(P >= 25) { i++; vector[i] = 0.0; }
        }
        else
        {
            for(i = 0 ; i < P ; i++) profile[i] = 0.0;
            for(j = aaSeqIndex[s]+1; j < aaSeqIndex[s+1]; j++)
            {
                if(((j - aaSeqIndex[s] <=W)||(aaSeqIndex[s+1] - j <=W))&&(score2[j] >= cutoff_score2))
                {
                    for(i = 0 ; i < 20 ; i++) 
                    {
                        if(encoding_type == PSSM_PROFILE_ENCODE) 
                            profile[i] += M[j][i] *W1;
                        else if(encoding_type == SEQUENCE_BINARY_ENCODE)
                            profile[i] += IS_EQUAL(aaSeq[j], STD1CharAA_alphabet[i]);
                    }
                    if(P >= 21) {profile[i] += score1[j] * W2; }
                    if(P >= 22) {i++; profile[i] += score2[j] * W3; }
                    if(P >= 23) {i++; profile[i] += 1.0;     }
                    if(P >= 24) {i++; profile[i] +=  Hydrophobicity[Char2Digit(aaSeq[j],STD1CharAA_alphabet)] *W4; }
                    if(P >= 25) {i++; profile[i] += (double)(Char2Digit(aaSeq[j],STD1CharAA_alphabet)+1)* W5;
                    }
                }
            }
            //for(i = 0 ; i < P ; i++) profile[i] /= double(inter);
            for(i = beg ; i< beg +P ; i++) 
                vector[i] = profile[i-beg];
        }
    }

    beg = (2*K + numSite + numSite-1)*P;
    if(vectorDim > beg)
    {
        i = beg;
        for(s = 0 ; s < numSite-1 ; s++) 
        {
            vector[i] = aaSeqIndex[s+1] - aaSeqIndex[s] -1;
            i ++;
        }
    }
    return vectorDim;
}
/*}}}*/
int GetSVMVector_2(int *aaSeqIndex, int numSite, MODM *pMODM, int K, int W, int P, double *vector, int vectorDim, double cutoff_score2 /*= 0.1*/, int encoding_type /*= PSSM_PROFILE_ENCODE*/)/*{{{*/
/****************************************************************************
 * 2007-06-26
 * adding other features, 
 *  1. give higher weights for centering residues in a window
 *  2. extra features for residues type, 
 *  C        4  1  0  0 0 
 *  H        1  4  0  0 0
 *  D        0  0  2  0 0
 *  E        0  0  0  1 0
 *  others   0  0  0  0 1
 *  3. water accessibility
 *
 * encoding of residues in between residues when numSite > 1, do not
 * calculate the average, average number is not good in doing dot product
 * encode each residue instead, but the a specific number , W, if there are
 * less then W residues between two specRes, then set to zero.
 *
 * vector dimension can be 
 * P is not the dimension of vector per residue, but the number of type of
 * features used to encode each residue, the order is
 *      PSSM_profile  SIZE_ENCODE_PSSM
 *      score1        SIZE_ENCODE_SCORE1
 *      score2        SZIE_ENCODE_SCORE2
 *      rangeIndicator 1
 *      amino_acid    SIZE_ENCODE_AA_CENTER || SIZE_ENCODE_AA_NON_CENTER
 *      wateracc      SIZE_ENCODE_WATERACC
 *
 * so P can be [1,6]
 * and dimVectorPerSite can be [SZIE_ENCODE_PSSM, SIZE_ENCODE_AA_CENTER||SIZE_ENCODE_AA_NON_CENTER+SIZE_ENCODE_WATERACC+SZIE_ENCODE_PSSM+1+SIZE_ENCODE_SCORE1+SIZE_ENCODE_SCORE2]
 * 
 *  the dimension of feature vector per sample is 
 *      (2*K + numSite + 2*W* (numSite-1) )*dimVectorPerSite 
 *      or if distance between residues used
 *      (2*K + numSite + 2*W* (numSite-1) )*dimVectorPerSite  +  SIZE_ENCODE_INTERDIST * (numSite-1)
 *
 * encoding for distance separating residues of a neighbouring pair is
 *
 *
 * return the size of the vector
 * return value <= 0 means failure
 * 
 * Encoding_type: PSSM_PROFILE_ENCODE         PSSM based encoding , 2007-04-27
 *                SEQUENCE_BINARY_ENCODE      sequence based binary encoding
 *
 * 2007-07-24, WaterACC can not be used for zinc prediction, since wateracc is
 * derived from 3D structures
 **************************************************************************/
{
    // -----k----1--w--internal-w--1---w---internal--w-1-k
    // if internel > 2*W, take only the start and end W

    // 1. for each residue, first 20 features are 20 psiblast profile
    // 2. 21st is the information content, but encoded as 
    //    0~0.3
    // 3.     
    // last two columns 2
    // flag indicating where overflow sequence 1 { -1, overflow, 1 not overflow }
    // hydrophobicity 1         

    double w1 = 1.0;   // weight for percentage profile
    double w2 = 1.0;   // weight for score1
    double w3 = 1.0;   // weight for score2
    double w4 = 1.0;   // weight for digit AA
    double w5 = 1.0;   // weight for water accessibility
    double weight_scale = 2.5; // if weight_scale = 4.0, the weight of the centering residue is 1.0+1/2.5 = 1.4, 1.4^2 = 1.96~~2, means in the dot product, the centering are twice weight compared to the furthest residue

    int i; // iterator for vector
    int j; // iterator for residue position
    int ik; // inner layer iterator
    int s;
    // double hydro; // hydrophobicity
    // int shift;
    int beg;
    //int end;
    int inter;

    bool isCenterResidue = false;
    int dimVectorPerSite  =0;
    int dimVectorPerSite_Center = 0;//dimension of feature vector to encode the centering residue
    int dimVectorPerSite_Non_Center = 0;//dimension of feature vector used to encode the non centering residue

    // check if the profile statisfying cutoff_score2 rule
    // 2007-08-08, do not neglect the the vector if there is only one residue
    // with low score2, 
    // but set the value of this residue to zero, if the score2 is too low
    //
    //beg = aaSeqIndex[0]-K; [>here beg end is for residue position<]
    //end = aaSeqIndex[numSite-1]+K;
    //for ( j = beg ; j <= end ;j ++)
    //{
    //    if(j >= 0 && j < pMODM->length)
    //    {
    //        if( pMODM->score2[j] < cutoff_score2)
    //            return -1 ;
    //    }
    //}

    dimVectorPerSite_Center = GetDimensionVectorPerSite(P, true);
    dimVectorPerSite_Non_Center= GetDimensionVectorPerSite(P, false);

    beg = 0; /*here beg is for the position in the vector*/
    //shift = aaSeqIndex[0] -K;
    for ( j = aaSeqIndex[0] -K ; j <= aaSeqIndex[0] ;j ++) // the beginning K residues
        // j is iterator for residue
    {
//        beg = (j-shift)*dimVectorPerSite_Non_Center;
        isCenterResidue = (j == aaSeqIndex[0]) ? true : false;
        dimVectorPerSite = GetDimensionVectorPerSite(P, isCenterResidue);
        double w = 1.0 +(1.0 - (aaSeqIndex[0]-j)/double (K)) / weight_scale;
        w1 = w2 = w3 = w4 = w5 = w;
        GetSVMVectorPerSite_2(vector, beg, j, P, dimVectorPerSite, pMODM, w1, w2, w3, w4, w5,cutoff_score2, isCenterResidue, encoding_type);
        beg += dimVectorPerSite;
    }

    //shift = aaSeqIndex[numSite-1] -(K+(numSite-1)+2*W*(numSite-1));  // the end
    beg = (K+2*W*(numSite-1))*dimVectorPerSite_Non_Center + (numSite-1)*dimVectorPerSite_Center;
    for ( j = aaSeqIndex[numSite-1]; j <= aaSeqIndex[numSite-1]+K ;j ++)
    {
//        beg = (j-shift)*dimVectorPerSite;
        isCenterResidue = (j == aaSeqIndex[numSite-1]) ? true : false;
        dimVectorPerSite = GetDimensionVectorPerSite(P, isCenterResidue);
        double w = 1.0 +(1.0 - (j - aaSeqIndex[numSite-1])/double (K)) / weight_scale;
        w1 = w2 = w3 = w4 = w5 = w;
        GetSVMVectorPerSite_2(vector, beg,j, P, dimVectorPerSite, pMODM, w1, w2, w3, w4, w5,cutoff_score2, isCenterResidue, encoding_type);
        beg += dimVectorPerSite;
    }

    for(s = 0 ; s < numSite-1 ; s++)// residues in between special residues
    {
        inter = aaSeqIndex[s+1] - aaSeqIndex[s] -1 ;

        beg =(K+W*2*s)*dimVectorPerSite_Non_Center +(s+1)*dimVectorPerSite_Center;
//        shift = aaSeqIndex[s] - (K + s + 2*W*s);
        for (j = aaSeqIndex[s]+1; j <= aaSeqIndex[s]+ W; j ++)/*for W residues after the residue at aaSeqIndex[s]*/
        {
 //           beg = (j - shift) * dimVectorPerSite;
            dimVectorPerSite = dimVectorPerSite_Non_Center;
            isCenterResidue = false;
            w1 = w2 = w3 = w4 = w5 = 1.0;
            if(j < aaSeqIndex[s+1])
            {
                GetSVMVectorPerSite_2(vector, beg, j, P, dimVectorPerSite, pMODM, w1, w2, w3, w4, w5,cutoff_score2, isCenterResidue, encoding_type);
                if((W-(j-aaSeqIndex[s])) < (2*W-inter)) // overlap region            
                {
                    for(ik = beg ; ik < beg + dimVectorPerSite; ik ++)
                    {
                        vector[ik] /= M_SQRT2;
                    }
                }
                beg += dimVectorPerSite;
            }
            else
            {
                for(ik = beg; ik < beg + dimVectorPerSite; ik ++)
                {
                    vector[ik] = 0.0;
                }
                beg += dimVectorPerSite;
            }
        }

        beg =(K+W*(2*s+1))*dimVectorPerSite_Non_Center + (s+1)*dimVectorPerSite_Center;
//        shift = aaSeqIndex[s+1] - (K + s+1 +2*W*(s+1) );
        for(j = aaSeqIndex[s+1] - W ; j < aaSeqIndex[s+1]; j ++) /*for W residue before the residue at aaSeqIndex[s+1]*/
        {
 //           beg = (j- shift)*dimVectorPerSite;
            dimVectorPerSite = dimVectorPerSite_Non_Center;
            isCenterResidue = false;
            w1 = w2 = w3 = w4 = w5 = 1.0;
            if(j > aaSeqIndex[s])
            {
                GetSVMVectorPerSite_2(vector, beg, j, P, dimVectorPerSite, pMODM, w1, w2, w3, w4, w5, cutoff_score2,   isCenterResidue, encoding_type);
                if((W-(aaSeqIndex[s+1]-j)) < (2*W-inter)) // overlap region            
                {
                    for(ik = beg ; ik < beg + dimVectorPerSite; ik ++)
                    {
                        vector[ik] /= M_SQRT2;
                    }
                }
                beg += dimVectorPerSite;
            }
            else
            {
                for(ik = beg; ik < beg + dimVectorPerSite; ik ++)
                {
                    vector[ik] = 0.0;
                }
                beg += dimVectorPerSite;
            }
        }
    }

    beg = (2*K +  2*W * (numSite-1)) *dimVectorPerSite_Non_Center + numSite*dimVectorPerSite_Center;
    if(vectorDim > beg) // encode distance between residues in a residue group
    {
        i = beg;
        for(s = 0 ; s < numSite-1 ; s++) 
        {
            int idx = GetEncodeDistance(aaSeqIndex[s+1]-aaSeqIndex[s]);
            for (ik = 0 ; ik < SIZE_ENCODE_DISTANCE; ik ++)
            {
                vector[i+ik] = EncodeMatrixDistance[idx][ik] ;
            }
            i += SIZE_ENCODE_DISTANCE;
        }
    }
    return vectorDim;
}
/*}}}*/

int GetMetalBoundRes(const char* metalProFile, MetalPro* metalPro, int &numMetalPro, bool isExcludeOtherChainRes/* = true*/, char **keyMetalList/* = NULL*/, int numKeyMetal/* = 0*/, int min_metalBoundRes /*= 0 */, int max_metalBoundRes /*= 0*/, bool isUsingTotalBoundRes /* = true*/)/*{{{*/
    /*****************************************************************************
     *  read in metal binding residues into struct MetalPro, 
     *  from the file metalProFile with the format
     *  1     101M    154    1
     *        FE :      1   HEM  4.00
     *          HIS93 ( )    
     *          HIS94 ( )    
     *          *            
     *          -99     
     * get metal binding residues for each metal binding protein, the metalPro is of
     * structure MetalPro, that is all the metal binding residues are stored
     * directly under metalPro, the redundent residues, i.e., those residues
     * binds to more than one metal atoms, are removed. Each atomEnv can still be
     * recovered from metalPro.atomEnv[i].parentResIndex, which stores the index of
     * the residues in that metal binding protein.
     *
     * isExcludeOtherChainRes is set to true by default, means that the metal
     * binding residues on other chains will be removed.
     *
     * keyMetalList is set to NULL and numKeyMetal is set to 0 by default, means
     * that all metal binding proteins are included.
     *
     * the default min_numRes and max_numRes are set to 0, which means there is no
     * restriction for the number of residues bound to metal atoms
     *
     * the structure of metalPro is like is
     *
     * mtalPro |--res-| -indeces to atomEnv
     *         |      |
     *         |       
     *         |--atomEnv-|- indeces to residues 
     *         |
     *         |
     *
     * 2007-04-10
     ****************************************************************************/
{
    using namespace MetalBindingProtein;
    char comment_char = '#';

    bool isNumBoundResSet = false;
    bool isKeyMetalSet = false;

    if(max_metalBoundRes != 0 ||  min_metalBoundRes != 0)
    {
        isNumBoundResSet = true;
    }
    if(numKeyMetal != 0)
    {
        isKeyMetalSet = true;
    }

    FILE* fpin; 
    fpin = fopen(metalProFile,"r");
    checkfilestream(fpin, metalProFile,"r");

    int  i, j;
    int  status_sscanf;
    int  index = 0;
    char id[SIZE_CHAIN_ID+1] = "";
    int  length = 0;
    int  numMetalAtom = 0;
    char metalAtomName[SIZE_METAL_ATOM_NAME+1] = "";
    int  metalAtomResSeq = 0;
    char metalAtomResName[SIZE_METAL_ATOM_RES_NAME+1] = "";
    char metalAtomChainID = 0;
    char chainID = ' ';
    int  numRes = 0;
    int  totalBoundRes = 0;
    int  maxline = 500;
    int  linesize;

    Array1D <char>    line_1darray(maxline+1);
    Array1D <Residue> tmpRes_1darray(MAX_BOUND_SITE);
    Array1D <int>     idx_1darray(MAX_BOUND_SITE);
    char    *line   = line_1darray.array1D;
    Residue *tmpRes = tmpRes_1darray.array1D;
    int     *idx    = idx_1darray.array1D;

    MetalPro tmpMetalPro;
    InitMetalPro(&tmpMetalPro);
    tmpMetalPro.res = new Residue[MAX_NUMRES];
    tmpMetalPro.metalAtomList = Create2DArray(tmpMetalPro.metalAtomList , MAX_ATOMENV, SIZE_ATOM_ELEMENT+1); 
    for(i = 0 ; i < MAX_NUMRES; i++) 
    {    
        InitResidue(&tmpMetalPro.res[i]);
        tmpMetalPro.res[i].parentAtomEnvIndex = new int[MAX_BOUND_METAL_PER_RES];
    }

    tmpMetalPro.resSeqIndex = new int[MAX_NUMRES];
    tmpMetalPro.atomEnv = new AtomEnv[MAX_ATOMENV];
    for(i = 0 ; i < MAX_ATOMENV; i++)
    {
        InitAtomEnv(&tmpMetalPro.atomEnv[i]);
        tmpMetalPro.atomEnv[i].parentResIndex = new int[MAX_BOUND_SITE];
    }


    Residue *pRes;
    AtomEnv *pAtomEnv;
    int cntMetalPro = 0;
    bool isHaveMetal;
    f_neglect_comment(fpin, comment_char);

    while((linesize = fgetline(fpin,line,maxline))!= EOF)
    {
        if(linesize <= 0 ) continue;
        if((status_sscanf = sscanf(line,"%d %5s %d %d",&index,id,&length,&numMetalAtom))< 4 ) 
        { 
            fprintf(stderr,"line =%s\n",line);
            assert(status_sscanf == 4);
        }

        isHaveMetal = false;
        StdID(id);         //standardize nrPDB chain identifier
        chainID = id[4]; /* chainID always refers to the chainID of protein chain, which is equal to
                            id[4] for the standardized chain identifier*/ 

        //idxPool.clear();

        tmpMetalPro.length     = length;
        tmpMetalPro.numBoundRes = 0;
        my_strcpy(tmpMetalPro.id,id,SIZE_CHAIN_ID);


        int cntNumRes = 0;
        int cntMetalAtom = 0;
        if(numMetalAtom > 0 )
        {
            for(i = 0 ; i < numMetalAtom ; i++)
            {
                for(j = 0 ; j < MAX_BOUND_SITE; j++ ) { InitResidue(&tmpRes[j]); }

                linesize = fgetline(fpin, line, maxline);
                if((status_sscanf = sscanf(line,"%2s %d :%1c %d %s",metalAtomName,&metalAtomResSeq, &metalAtomChainID,&totalBoundRes,metalAtomResName)) < 5 )
                {
                    fprintf(stderr,"line =%s\n", line);
                    assert(status_sscanf == 5);
                }
                my_strupr(metalAtomName); //notice that the metal element name in output is like Na, while in processing are all upper case.

                if( totalBoundRes > 0)
                {
                    ScanfCloseMetalRes(fpin,tmpRes,totalBoundRes);//read in residues bound to this metal atom
                    //for(j = 0 ; j < numRes; j++) my_strcpy(tmpRes[j].metalAtomName, metalAtomName, SIZE_METAL_ATOM_NAME);

                    // get the numRes for that chain
                    if(isExcludeOtherChainRes) // remove residues bound to that metal atom but not on the current chain
                    {
                        int cntRes = 0;
                        for(j = 0; j < totalBoundRes ; j++)
                        {
                            if(tmpRes[j].chainID == chainID) { idx[cntRes] = j; cntRes++; }
                        }
                        numRes = cntRes; // numRes record the number of residues (binding to the meta atom) on a single chain,
                    }
                    else
                    {
                        numRes = totalBoundRes;
                        for(j = 0 ; j < numRes; j ++) idx[j] = j;
                    }

                    if((!isKeyMetalSet || BinarySearch_String(metalAtomName, keyMetalList, numKeyMetal) != -1)
                            && (!isNumBoundResSet 
                                || (!isUsingTotalBoundRes && numRes >= min_metalBoundRes && numRes <= max_metalBoundRes)
                                || (isUsingTotalBoundRes && totalBoundRes >= min_metalBoundRes && totalBoundRes <= max_metalBoundRes)
                               )
                      )
                    {
                        isHaveMetal = true;
                        my_strcpy(tmpMetalPro.metalAtomList[cntMetalAtom], metalAtomName, SIZE_ATOM_ELEMENT);
                        tmpMetalPro.atomEnv[cntMetalAtom].numRes = numRes;
                        tmpMetalPro.atomEnv[cntMetalAtom].totalBoundRes = totalBoundRes;
                        my_strcpy(tmpMetalPro.atomEnv[cntMetalAtom].metalAtomName, metalAtomName, SIZE_METAL_ATOM_NAME);
                        my_strcpy(tmpMetalPro.atomEnv[cntMetalAtom].metalAtomResName, metalAtomResName, SIZE_METAL_ATOM_RES_NAME);
                        tmpMetalPro.atomEnv[cntMetalAtom].metalAtomResSeq = metalAtomResSeq;
                        tmpMetalPro.atomEnv[cntMetalAtom].metalAtomChainID = metalAtomChainID;
                        my_strcpy(tmpMetalPro.atomEnv[cntMetalAtom].metalAtomPDBID,id,SIZE_PDBID);
                        my_strcpy(tmpMetalPro.atomEnv[cntMetalAtom].id,id,SIZE_CHAIN_ID);
                        tmpMetalPro.atomEnv[cntMetalAtom].seqLength = length;

                        for ( j = 0; j < numRes ; j++)
                        {
                            pRes = & (tmpRes[idx[j]]);
                            //int *pSeqIndex;
                            //pSeqIndex = find (tmpMetalPro.resSeqIndex, tmpMetalPro.resSeqIndex+cntNumRes, pRes->aaSeqIndex);
                            //LOG: 2007-04-13 16:42:24 Friday  Week 14 <nanjiang@casio.fos.su.se>
                            // when isExcludeOtherChainRes is false, only pSeqIndex can not distinguish
                            // whether a residue in unique. chain id shold be included for identification
                            // of a new residue
                            //if(pSeqIndex == tmpMetalPro.resSeqIndex+cntNumRes)[> not found, new residue<] 

                            int resIndex = GetResidueIndex(pRes, tmpMetalPro.res, cntNumRes);
                            if(resIndex == -1) //not found, new residue
                            {
#ifdef _ASSERT_
                                if(cntNumRes >= MAX_NUMRES )
                                {
                                    fprintf(stderr,"tmpMetalPro.id=%s, cntNumRes=%d\n", tmpMetalPro.id, cntNumRes);
                                    assert(cntNumRes < MAX_NUMRES );
                                }
                                if(cntMetalAtom >= MAX_ATOMENV )
                                {
                                    fprintf(stderr,"tmpMetalPro.id=%s, cntMetalAtom = %d\n", tmpMetalPro.id, cntMetalAtom);
                                    assert(cntMetalAtom < MAX_ATOMENV);
                                }
#endif
                                CopyResidue( & (tmpMetalPro.res[cntNumRes]), pRes);
                                tmpMetalPro.resSeqIndex[cntNumRes] = pRes->aaSeqIndex;   //note that resSeqIndex might not be unique if isExcludeOtherChainRes = false
                                tmpMetalPro.atomEnv[cntMetalAtom].parentResIndex[j] = cntNumRes;
                                cntNumRes ++;
                            }
                            else
                            {
                                tmpMetalPro.atomEnv[cntMetalAtom].parentResIndex[j] = resIndex;
                            }
                        }
                        cntMetalAtom ++;
                    }
                }
            }
            //mapping atomEnv to the parentAtomEnvIndex under res
            for ( i = 0 ; i < cntMetalAtom; i++)
            {
                pAtomEnv = &(tmpMetalPro.atomEnv[i]);
                for(j = 0 ; j < pAtomEnv -> numRes; j ++)
                {
                    int idxRes = pAtomEnv->parentResIndex[j];
                    int cntMetalBound = tmpMetalPro.res[idxRes].numMetalBound;
#ifdef _ASSERT_
                    if(idxRes >= MAX_NUMRES )
                    {
                        fprintf(stderr,"tmpMetalPro.id=%s, idxRes=%d\n", tmpMetalPro.id, idxRes);
                        assert(idxRes < MAX_NUMRES );
                    }
                    if(cntMetalBound >= MAX_BOUND_METAL_PER_RES )
                    {
                        fprintf(stderr,"tmpMetalPro.id=%s, cntMetalBound = %d\n", tmpMetalPro.id, cntMetalBound);
                        assert(cntMetalBound < MAX_BOUND_METAL_PER_RES);
                    }
#endif
                    tmpMetalPro.res[idxRes].parentAtomEnvIndex[cntMetalBound] = i;
                    tmpMetalPro.res[idxRes].numMetalBound ++;
                }
            }
        }
        tmpMetalPro.numBoundRes = cntNumRes;
        tmpMetalPro.numMetalAtom = cntMetalAtom;
        if(isHaveMetal)
        {
            metalPro[cntMetalPro].res = new Residue[tmpMetalPro.numBoundRes];
            for( j = 0; j < tmpMetalPro.numBoundRes; j++ )
            {
                InitResidue(&(metalPro[cntMetalPro].res[j]));
                metalPro[cntMetalPro].res[j].parentAtomEnvIndex = new int[tmpMetalPro.res[j].numMetalBound];
            }

            metalPro[cntMetalPro].resSeqIndex = new int[tmpMetalPro.numBoundRes];
            metalPro[cntMetalPro].atomEnv = new AtomEnv[tmpMetalPro.numMetalAtom];
            for( j = 0; j < tmpMetalPro.numMetalAtom; j++ ) InitAtomEnv(&(metalPro[cntMetalPro].atomEnv[j]));
            for( j = 0 ; j < tmpMetalPro.numMetalAtom ; j ++)
            {
                metalPro[cntMetalPro].atomEnv[j].parentResIndex = new int[tmpMetalPro.atomEnv[j].numRes];
            }
            metalPro[cntMetalPro].metalAtomList = Create2DArray(metalPro[cntMetalPro].metalAtomList,cntMetalAtom, SIZE_ATOM_ELEMENT+1);
            CopyMetalPro(&metalPro[cntMetalPro], &tmpMetalPro,true);
            cntMetalPro ++;
        }
    }

    numMetalPro = cntMetalPro;

    fclose(fpin);
    /*free memory*/
    DeleteMetalPro(&tmpMetalPro, MAX_ATOMENV, MAX_NUMRES);

    return numMetalPro;
}/*}}}*/
int RestrictMetalBoundRes(MetalPro *metalPro1, int &numMetalPro1, int &numMetalRes1, MetalPro* metalPro2, int numMetalPro2, char *resList, int level, const char *modmpath /* = ""*/, double cutoff_score2 /* = 0.0*/)/*{{{*/
    /*****************************************************************************
     * restrict the metalPro1 by different restrict rules, 
     * level = 1, restricted by resList
     * level = 0, restricted by resList and cutoff_score2
     * input: metalPro2, 
     * copy metalPro2 to metalPro1 with the restriction rule
     ****************************************************************************/
{
    using namespace MetalBindingProtein;
    int i, j, k;
    char modmfilepath[MAX_PATH+1] = "";

    Residue *pRes1;
    Residue *pRes2;
    AtomEnv *pAtomEnv;
    Residue *pRes;
    int cntMetalRes;//count all metal binding residues in metalPro
    int cntMetalPro;

    MetalPro tmpMetalPro;
    InitMetalPro(&tmpMetalPro);
    tmpMetalPro.res = new Residue[MAX_NUMRES];
    tmpMetalPro.metalAtomList = Create2DArray(tmpMetalPro.metalAtomList , MAX_ATOMENV, SIZE_ATOM_ELEMENT+1); 
    for(i = 0 ; i < MAX_NUMRES; i++) 
    {    
        InitResidue(&tmpMetalPro.res[i]);
        tmpMetalPro.res[i].parentAtomEnvIndex = new int[MAX_BOUND_METAL_PER_RES];
    }

    tmpMetalPro.resSeqIndex = new int[MAX_NUMRES];
    tmpMetalPro.atomEnv = new AtomEnv[MAX_ATOMENV];
    for(i = 0 ; i < MAX_ATOMENV; i++)
    {

        InitAtomEnv(&tmpMetalPro.atomEnv[i]);
        tmpMetalPro.atomEnv[i].parentResIndex = new int[MAX_BOUND_SITE];
    }

    double *score1 = NULL;
    double *score2 = NULL;
    char *aaSeq = NULL;
    DATATYPE_MODM_MATRIX    **M = NULL;
    int length;
    char alphabetMODM[NUM_BLOSUM+1] = "";
    char id[SIZE_CHAIN_ID+1] = "";

    //Array1D <double> score1_1darray(LONGEST_SEQ);
    Array1D <double> score2_1darray(LONGEST_SEQ);
    //Array1D <char>   aaSeq_1darray(LONGEST_SEQ+1);
    //Array2D <DATATYPE_MODM_MATRIX>    M_2darray(LONGEST_SEQ, NUM_BLOSUM);
    //score1 = score1_1darray.array1D;
    score2 = score2_1darray.array1D;
    //aaSeq  = aaSeq_1darray.array1D;
    //M      = M_2darray.array2D;


    /*get the metal binding protein information for metalPro and metalPro1*/
    cntMetalPro = 0;
    cntMetalRes = 0;
    for( i = 0; i < numMetalPro2 ; i++)
    {
        CopyMetalPro( &tmpMetalPro, &metalPro2[i],false); // copy metalPro2[i] to tmpMetalPro, but does not copy metalPro2[i].res
        int cntBoundRes = 0 ; //count metal binding residues bound to metalPro2[i]

        if(level == 0)
        {
            my_strcpy(id, metalPro2[i].id, SIZE_CHAIN_ID); id[SIZE_CHAIN_ID] = '\0';
            if(GetMODMFilePath(id,modmfilepath, modmpath) != NULL)
                length = GetMODM(modmfilepath, M, alphabetMODM, aaSeq, score1, score2);
            else
            {
                assert(GetMODMFilePath(id,modmfilepath,modmpath) != NULL);
            }
        }

        //get the residues statisfying restriction rule
        for(j = 0 ; j < metalPro2[i].numBoundRes; j++)
        {
            assert(cntBoundRes >= 0 && cntBoundRes < MAX_NUMRES);
            pRes1 = & (tmpMetalPro.res[cntBoundRes]);
            pRes2 = & (metalPro2[i].res[j]);
            if(IsInCharSet(pRes2->aa,resList) && 
                    (level == 1  // metalPro1, restrict by resList
                     || (score2[pRes2->aaSeqIndex] >= cutoff_score2 && level == 0))  // metalPro restricted by both resList and cutoff_score2
              )
            {
                CopyResidue(pRes1,pRes2);
                tmpMetalPro.resSeqIndex[cntBoundRes] = pRes2->aaSeqIndex;
                cntBoundRes ++;
                cntMetalRes ++;
            }
        }
        tmpMetalPro.numBoundRes = cntBoundRes;

        int cntMetalAtom = 0;
        for(j = 0 ; j < metalPro2[i].numMetalAtom; j++)
        {
            pAtomEnv = &(metalPro2[i].atomEnv[j]);
            int cntNumRes = 0;
            for(k = 0 ; k < metalPro2[i].atomEnv[j].numRes; k ++)
            {
                pRes = &(metalPro2[i].res[pAtomEnv->parentResIndex[k]]);
                int resIndex = GetResidueIndex(pRes, tmpMetalPro.res, tmpMetalPro.numBoundRes);//get the resIndex for in the updated tmpMetalPro.res, 2007-05-04
                // bug fixed, cntNumRes change to tmpMetalPro.numBoundRes,  2007-05-25
                if(resIndex != -1) //if this residue is in the updated tmpMetalPro.res, add it
                {
                    tmpMetalPro.atomEnv[cntMetalAtom].parentResIndex[cntNumRes] = resIndex;
                    cntNumRes ++;
                }

                //int* pSeqIndex = find (tmpMetalPro.resSeqIndex, tmpMetalPro.resSeqIndex+cntBoundRes, metalPro2[i].resSeqIndex[metalPro2[i].resSeqIndex[metalPro2[i].atomEnv[j].parentResIndex[k]]]);
                //if(pSeqIndex != tmpMetalPro.resSeqIndex+cntBoundRes) // if found
                //{    
                //tmpMetalPro.atomEnv[cntMetalAtom].parentResIndex[cntNumRes] = pSeqIndex-tmpMetalPro.resSeqIndex;
                //cntNumRes ++;
                //}
            }
            tmpMetalPro.atomEnv[cntMetalAtom].numRes = cntNumRes;
            if (cntNumRes > 0)
                cntMetalAtom ++;
        }
        tmpMetalPro.numMetalAtom = cntMetalAtom;

        //re-generate tmpMetalPro.res[j].parentAtomEnvIndex, 2007-04-25
        // Initialize tmpMetalPro.res[j].numMetalBound, parentAtomEnvIndex
        for(j = 0 ; j < tmpMetalPro.numBoundRes; j++)
        {
            tmpMetalPro.res[j].numMetalBound = 0;
        }
        //mapping atomEnv to the parentAtomEnvIndex under res
        for ( j = 0 ; j < tmpMetalPro.numMetalAtom; j++)
        {
            pAtomEnv = &(tmpMetalPro.atomEnv[j]);
            for(k = 0 ; k < pAtomEnv -> numRes; k ++)
            {
                int idxRes = pAtomEnv->parentResIndex[k];
                int cntMetalBound = tmpMetalPro.res[idxRes].numMetalBound;
                assert(idxRes >= 0 && idxRes < MAX_NUMRES);
#ifdef DEBUG
                if( cntMetalBound < 0 || cntMetalBound >= MAX_BOUND_METAL_PER_RES)
                {
                    DEBUG_STDERR(cntMetalBound,"%d");
                }
#endif /* !DEBUG */
                assert(cntMetalBound >= 0 && cntMetalBound < MAX_BOUND_METAL_PER_RES);
                tmpMetalPro.res[idxRes].parentAtomEnvIndex[cntMetalBound] = j;
                tmpMetalPro.res[idxRes].numMetalBound ++;
            }
        }

        if(tmpMetalPro.numBoundRes > 0)
        {
            metalPro1[cntMetalPro].res = new Residue[tmpMetalPro.numBoundRes];
            for(j = 0 ; j < tmpMetalPro.numBoundRes; j++) 
            {
                InitResidue(&(metalPro1[cntMetalPro].res[j]));
                metalPro1[cntMetalPro].res[j].parentAtomEnvIndex = new int[tmpMetalPro.res[j].numMetalBound];
            }

            metalPro1[cntMetalPro].resSeqIndex = new int[tmpMetalPro.numBoundRes];
            metalPro1[cntMetalPro].atomEnv = new AtomEnv[tmpMetalPro.numMetalAtom];
            for(j = 0 ; j < tmpMetalPro.numMetalAtom; j++) InitAtomEnv(&(metalPro1[cntMetalPro].atomEnv[j]));
            for( j = 0 ; j < tmpMetalPro.numMetalAtom ; j ++)
            {
                metalPro1[cntMetalPro].atomEnv[j].parentResIndex = new int[tmpMetalPro.atomEnv[j].numRes];
            }
            CopyMetalPro(&metalPro1[cntMetalPro], &tmpMetalPro);/* copy tmpMetalPro to metalPro[numMetalPro]*/
            cntMetalPro ++;
        }
    } 
    numMetalPro1 = cntMetalPro;
    numMetalRes1 = cntMetalRes;

    DeleteMetalPro(&tmpMetalPro,MAX_ATOMENV,  MAX_NUMRES);

    return numMetalPro1;
}
/*}}}*/
int GetMetalBoundRes3(const char* metalProFile, MetalPro* metalPro, MetalPro* metalPro1, MetalPro* metalPro2, int &numMetalPro, int &numMetalPro1, int &numMetalPro2, int &numMetalRes, int &numMetalRes1, int &numMetalRes2, double cutoff_score2, char *resList, const char* modmpath, bool isExcludeOtherChainRes/* = true*/, char **keyMetalList/* = NULL*/, int numKeyMetal/* = 0*/, int min_metalBoundRes /*= 0 */, int max_metalBoundRes /*= 0*/, bool isUsingTotalBoundRes /* = true*/)/*{{{*/
    /*****************************************************************************
     * get the metal binding residues for each metal binding protein 
     *   metalPro   -- store metal binding residues bound by metal atoms which bind
     *                 to >= min_metalBoundRes and <= max_metalBoundRes residues, within resList
     *                 and filtered by cutoff_score2
     *   metalPro1  -- same as metalPro, but the cutoff_score2 restriction is removed, the resList rule keeps
     *   metalPro2  -- same as metalPro, but both resList and cutoff_score2 restriction are removed, that is 
     *                 store all metal binding residues bound by metal atoms which
     *                 bind to >= min_metalBoundRes and <= max_metalBoundRes (in
     *                 the whole protein) residues. isExcludeOtherChainRes = true,
     *                 that is only residues on the nrPDB chain is included, 
     *                 
     ****************************************************************************/
{

    numMetalPro2 = GetMetalBoundRes(metalProFile, metalPro2, numMetalPro2, isExcludeOtherChainRes, keyMetalList, numKeyMetal, min_metalBoundRes, max_metalBoundRes, isUsingTotalBoundRes); 
    /* numMetalRes2 might not equal to the summation of all numRes over all metal
     * binding site, since there are cases one residue binding to more than one
     * metal ion */ 

    int i;
    int cntMetalRes = 0;
    for(i = 0; i < numMetalPro2; i++) 
        cntMetalRes += metalPro2[i].numBoundRes;
    numMetalRes2 = cntMetalRes;

    int level;
    /*get metalPro1, which is restricted resList */
    level = 1;
    numMetalPro1 = RestrictMetalBoundRes(metalPro1, numMetalPro1, numMetalRes1, metalPro2, numMetalPro2, resList, level);

    /*get metalPro, which is restricted by both resList and cutoff_score2*/
    level = 0;
    numMetalPro = RestrictMetalBoundRes(metalPro, numMetalPro, numMetalRes, metalPro2, numMetalPro2, resList, level, modmpath, cutoff_score2);

    return numMetalPro2;
}
/*}}}*/
int GetMetalBoundRes2(const char* metalProFile, MetalPro* metalPro1, MetalPro* metalPro2, int &numMetalPro1, int &numMetalPro2, int &numMetalRes1, int &numMetalRes2, char *resList, const char* modmpath, bool isExcludeOtherChainRes/* = true*/, char **keyMetalList/* = NULL*/, int numKeyMetal/* = 0*/, int min_metalBoundRes /*= 0 */, int max_metalBoundRes /*= 0*/, bool isUsingTotalBoundRes /* = true*/)/*{{{*/
    /*****************************************************************************
     * get the metal binding residues for each metal binding protein 
     *   metalPro1   -- store metal binding residues bound by metal atoms which bind
     *                 to >= min_metalBoundRes and <= max_metalBoundRes residues, within resList
     *                 no cutoff_score2 restriction, thus modm files do not
     *                 need to be read
     *   metalPro2  -- same as metalPro1, but both resList and cutoff_score2 restriction are removed, that is 
     *                 store all metal binding residues bound by metal atoms which
     *                 bind to >= min_metalBoundRes and <= max_metalBoundRes (in
     *                 the whole protein) residues. isExcludeOtherChainRes = true,
     *                 that is only residues on the nrPDB chain is included, 
     *                 
     ****************************************************************************/
{

    numMetalPro2 = GetMetalBoundRes(metalProFile, metalPro2, numMetalPro2, isExcludeOtherChainRes, keyMetalList, numKeyMetal, min_metalBoundRes, max_metalBoundRes, isUsingTotalBoundRes); 
    /* numMetalRes2 might not equal to the summation of all numRes over all metal
     * binding site, since there are cases one residue binding to more than one
     * metal ion */ 

    int i;
    int cntMetalRes = 0;
    for(i = 0; i < numMetalPro2; i++) 
        cntMetalRes += metalPro2[i].numBoundRes;
    numMetalRes2 = cntMetalRes;

    int level;
    /*get metalPro1, which is restricted resList */
    level = 1;
    numMetalPro1 = RestrictMetalBoundRes(metalPro1, numMetalPro1, numMetalRes1, metalPro2, numMetalPro2, resList, level);

    return numMetalPro2;
}
/*}}}*/

int GetNumBoundRes(set<int> const& PrP, set<int> const& PrN, MetalPro* pZnPro, set <int> &TP, set <int> &FP, set<int> &TN, set<int> &FN) /*{{{*/
{
    int i;
    set <int> idxZnRes;
    for (i = 0; i< pZnPro->numBoundRes; i++)
        idxZnRes.insert(pZnPro->res[i].aaSeqIndex);

    set_intersection(PrP.begin(),PrP.end(), idxZnRes.begin(), idxZnRes.end(), inserter(TP,TP.begin()));
    set_difference(  PrP.begin(),PrP.end(), idxZnRes.begin(), idxZnRes.end(), inserter(FP,FP.begin()));

    set_intersection(PrN.begin(),PrN.end(), idxZnRes.begin(), idxZnRes.end(), inserter(FN,FN.begin()));
    set_difference(  PrN.begin(),PrN.end(), idxZnRes.begin(), idxZnRes.end(), inserter(TN,TN.begin()));

    return TP.size();
}/*}}}*/
int GetNumBoundRes(set<int> const& idx, SSBondPro* pPro, set <int> &idxTrue, set <int> &idxFalse) /*{{{*/
{
    int i,j ;
    set <int> idxRes;
    for (i = 0; i< pPro->numSSBond; i++)
    {    
        for(j = 0 ; j < pPro->ssbond[i].numRes; j++)
            idxRes.insert(pPro->ssbond[i].res[j].aaSeqIndex);
    }

    set_intersection(idx.begin(),idx.end(), idxRes.begin(), idxRes.end(), inserter(idxTrue,idxTrue.begin()));
    set_difference(idx.begin(),idx.end(), idxRes.begin(), idxRes.end(), inserter(idxFalse,idxFalse.begin()));
    return idxTrue.size();
}/*}}}*/
int GetNumBoundRes(set<int> const& idx, MetalPro* pZnPro, set <int> &idxTrue, set <int> &idxFalse) /*{{{*/
{
    int i;
    set <int> idxZnRes;
    for (i = 0; i< pZnPro->numBoundRes; i++)
        idxZnRes.insert(pZnPro->res[i].aaSeqIndex);
    set_intersection(idx.begin(),idx.end(), idxZnRes.begin(), idxZnRes.end(), inserter(idxTrue,idxTrue.begin()));
    set_difference(idx.begin(),idx.end(), idxZnRes.begin(), idxZnRes.end(), inserter(idxFalse,idxFalse.begin()));
    return idxTrue.size();
}/*}}}*/
int GetNumBoundRes(set<int> const& idx, MetalPro2* pZnPro, set <int> &idxTrue, set<int> &idxFalse) /*{{{*/
{
    int i,j;
    set <int> idxZnRes;
    for (i = 0; i< pZnPro->numMetalAtom; i++)
    {
        for(j = 0; j < pZnPro->atomEnv[i].numRes; j++)
            idxZnRes.insert(pZnPro->atomEnv[i].res[j].aaSeqIndex);
    }
    set_intersection(idx.begin(),idx.end(), idxZnRes.begin(), idxZnRes.end(), inserter(idxTrue,idxTrue.begin()));
    set_difference(idx.begin(),idx.end(), idxZnRes.begin(), idxZnRes.end(), inserter(idxFalse,idxFalse.begin()));
    return idxTrue.size();
}/*}}}*/
int GetHCRes(MODM *pMODM, bool *isPolyHis, Residue* HCRes, double cutoff_score2, double cutoff_consv, char* resList, bool isUseConsvI, bool isMaskPolyHis)/*{{{*/
    // return the number of HCRes
{
    int i;
    int daa;
    DATATYPE_CONSV consv;
    int      type_modm    = pMODM -> type_modm;
    int      length       = pMODM -> length;
    DATATYPE_MODM_MATRIX    **M         = pMODM -> M;
    //double  *score1       = pMODM -> score1;
    double  *score2       = pMODM -> score2;
    char    *aaSeq        = pMODM -> aaSeq;
    char    *alphabetMODM = pMODM -> alphabetMODM;
    int sizeAlphabet = strlen(alphabetMODM);

    int cntRes  = 0;
    for (i = 0 ; i < length ; i++)
    {
        daa = Char2Digit(aaSeq[i],alphabetMODM,sizeAlphabet);
        if(daa < 0) continue;

        if(isUseConsvI)
            consv = GetIntegConsv(aaSeq[i], M[i], alphabetMODM, type_modm);
        else
            consv = M[i][daa];
        if(IsInCharSet(aaSeq[i], resList) 
                && consv >= cutoff_consv && score2[i] >= cutoff_score2
                && (!isPolyHis[i] || !isMaskPolyHis))
        {
            HCRes[cntRes].aa         = aaSeq[i];
            HCRes[cntRes].consv      = consv;
            HCRes[cntRes].aaSeqIndex = i;
            cntRes ++;
        }
    }
    return cntRes ;
}/*}}}*/

int ReadSeq_FASTA(const char *fileName, char* seq, int *pSeq_type /*= NULL*/, int maxlength /*= LONGEST_SEQ*/, char* annotationLine/*=NULL*/, int maxSizeAnnotationLine /*=50*/)/*{{{*/
    /****************************************************************************
     * ReadSeq_FASTA()
     * read in fasta format sequence file which contains single sequence
     * return the length of the sequence if successful
     * return -1 if can not open file

     * LOG: 2006-04-26 16:27:02 Wednesday  Week 17 <nanjiang@shu>
     * bug fixed for 0x0d character in unix system, in windows, '\n' is 0x0d0x0a,
     * while under unix, '\n' = 0x0a;
     * also the leading white spaces are neglected
     * LOG: 2006-06-14 16:36:56 Wednesday  Week 24 <nanjiang@casio>
     *   check the type (DNA or AA) of fasta sequence file, using the annotation line 
     * Updated 2010-04-22: the annotation line (without the leading ">") can be
     *   read in, note that the maxSizeAnnotationLine must be supplied
     ***************************************************************************/
{
    int   c;
    int   i;
    FILE *fp;
    if((fp = fopen(fileName,"r")) == NULL)
    {
        fprintf(stderr, "can not open file %s \n",fileName) ;
        return READ_FILE_ERROR;
    }

    do{  /* neglect the leading white spaces */ 
        c = getc(fp);
    }while(isspace(c));

    if(pSeq_type != NULL) 
        *pSeq_type = UNKNOWN_SEQ_TYPE; /* initializing sequence type*/

    if(c  == '>') 
    { 
        Array1D <char> line_1darray(maxSizeAnnotationLine +1);
        char *line = line_1darray.array1D;
        fgetline(fp,line,maxSizeAnnotationLine);
        if( pSeq_type != NULL)
        {
            if(strstr(line, "protein") != NULL)
                *pSeq_type = AA_SEQ;
            else if( strstr(line,"nucleic") != NULL)
                *pSeq_type = DNA_SEQ;
            else if( strstr(line,"shape") != NULL)
                *pSeq_type = SHAPE_SEQ;
            else
                *pSeq_type = UNKNOWN_SEQ_TYPE;
        }
        if (annotationLine != NULL)
        {
            my_strcpy(annotationLine,line, maxSizeAnnotationLine);  /*read in the annotation line*/ 
        }
    }
    else  
    {
        fseek(fp, -1, SEEK_CUR); /* backward 1 byte of file stream if there is no annotation line*/ 
    }

    i = 0 ;
    while((c = getc(fp)) != EOF) 
    {
        if(c == '>') break;  /* read in the first sequence if there are multiple sequences in the file*/ 
        if( ! isspace(c)) /* neglect white spaces and return characters in sequence region*/ 
        {
            seq[i] = c ; i ++ ;
            if(i >= maxlength)
            {
                fprintf(stderr,"Error, sequence longer then maxlength = %d\n", maxlength);
                return -1;
            }
        }
    }
    fclose(fp);

    seq[i] = '\0' ;
    return i; /* return the length of sequence*/ 
}/*}}}*/
int ReadNextSeq_FASTA(FILE *fp, char* seq, int *pSeq_type /*= NULL*/, int maxlength /*= LONGEST_SEQ*/, char* annotationLine/*=NULL*/, int maxSizeAnnotationLine /*=50*/)/*{{{*/
    /****************************************************************************
     * ReadNextSeq_FASTA()
     * read in the fasta format sequence from the file stream 
     * return the length of the sequence if successful
     * return 0 or minus value if no more sequence can be read from the file stream
     * The leading white spaces are ignored
     * check the type (DNA or AA) of fasta sequence file based the annotation line 
     * Last modified, 2007-02-12, Nanjiang Shu
     * Updated 2010-04-22: the annotation line (without the leading ">") can be
     * read in, note that the maxSizeAnnotationLine must be supplied
     ***************************************************************************/
{
    int   c;
    int   i;

    do{  /* ignore the leading white spaces */ 
        c = getc(fp);
    }while(isspace(c));

    if(pSeq_type != NULL) 
        *pSeq_type = UNKNOWN_SEQ_TYPE; /* initializing sequence type*/

    if(c  == '>') 
    { 
        Array1D <char> line_1darray(maxSizeAnnotationLine +1);
        char *line = line_1darray.array1D;
        fgetline(fp,line,maxSizeAnnotationLine);
        if( pSeq_type != NULL)
        {
            if(strstr(line, "protein") != NULL)
                *pSeq_type = AA_SEQ;
            else if( strstr(line,"nucleic") != NULL)
                *pSeq_type = DNA_SEQ;
            else if( strstr(line,"shape") != NULL)
                *pSeq_type = SHAPE_SEQ;
            else
                *pSeq_type = UNKNOWN_SEQ_TYPE;
        }
        if (annotationLine != NULL)
        {
            my_strcpy(annotationLine,line, maxSizeAnnotationLine);  /*read in the annotation line*/ 
        }
    }
    else  
    {
        fseek(fp, -1, SEEK_CUR); /* backward 1 byte of file stream if there is no annotation line*/ 
    }

    i = 0 ;
    while((c = getc(fp)) != EOF) 
    {
        if(c == '>') 
        {
            fseek(fp, -1, SEEK_CUR); /* backward 1 byte of file stream*/
            break;  /* read in the first sequence if there are multiple sequences in the file*/ 
        }
        if(! isspace(c)) /* neglect white spaces and return characters in sequence region*/ 
        {
            seq[i] = c ; i ++ ;
            if(i >= maxlength)
            {
                fprintf(stderr,"Error, sequence longer then maxlength = %d\n", maxlength);
                return 1;

            }
        }
    }
    seq[i] = '\0' ;
    if(c == EOF && i == 0)
        return EOF;
    else
        return i; /* return the length of sequence*/ 
}/*}}}*/
int GetNextResidue_PDB(Residue *pRes, FILE *fpPDBFile, bool isGetAllAtomLocation  /*=false*/)/*{{{*/
    /*****************************************************************************
     * Get all atom records of the current residue in pdb file, 
     * if isGetAllAtomLocation == false, for atoms with altLoc, only one of them will be chosen and stored in pRes
     * if isGetAllAtomLocation == true, all atom locations will be read in and stored in pRes
     * the memory pRes->atom will be allocated here, allocated memory = pRes->numAtom
     * 2004-04-23
     ****************************************************************************/
{
    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;

    int i,j;
    int atomSerial = 0;
    bool isEndMDLreached = false;
    Array1D <Atom> atom_1darray(MAX_ATOM_PER_RES_WITH_ALTATOM);
    Atom *atom = atom_1darray.array1D;
    for(i = 0 ; i < MAX_ATOM_PER_RES_WITH_ALTATOM; i++) { InitAtom(&atom[i]); }

    int status_fpos;
    fpos_t pos;
    int cntAtom = 0;
    char recordID[SIZE_RECORD_ID+1] = "";

    status_fpos = fgetpos(fpPDBFile,&pos);
    while((linesize = fgetline(fpPDBFile, line, maxline))!= EOF)
    {
        if( linesize <= 0) continue;
        else if( sscanf(line,"%6s",recordID) < 1 ) continue;
        else if( strcmp(recordID,"ANISOU")== 0 ) { continue; } // neglect the "ANISOU" record
        else if( strcmp(recordID,"ENDMDL") == 0)
        {
            isEndMDLreached = true;
            fprintf(stderr,"End of model, take the first model only\n");
            if(atomSerial == MAX_ATOM_SERIAL)
                fprintf(stderr,"Atom Serial is possibly overflow (> MAX_ATOM_SERIAL = %d), CHECK!\n", MAX_ATOM_SERIAL);
        }
        else if(strcmp(recordID,"ATOM") == 0 || strcasecmp(recordID,"HETATM") == 0) // if the line atom record,
            {
                ScanfCoorRecord_Atom(line,&atom[cntAtom]);
                atomSerial = atom[cntAtom].serial;
                if(atom[cntAtom].chainID == atom[0].chainID
                        && atom[cntAtom].resSeq == atom[0].resSeq
                        && atom[cntAtom].iCode  == atom[0].iCode ) // if still in the same residue

                {
                    cntAtom ++;
                    status_fpos = fgetpos(fpPDBFile, &pos); //retrieve the file position before the next getline
                    assert (status_fpos == 0);
                }
                else
                {
                    status_fpos = fsetpos(fpPDBFile,&pos); // set back the file pos to the previous position, if the current line is already in another residue
                    assert (status_fpos == 0);
                    break;
                }
            }
        else
        {
            status_fpos = fsetpos(fpPDBFile,&pos); // set back the file pos to the previous position, if the current line is another type of recordID, (of course not in the same residue )
            assert (status_fpos == 0);
            break;
        }
    }

    //dealing with this atoms
    if((linesize == EOF || isEndMDLreached) && cntAtom == 0) 
    {
        return EOF;
    }
    else if(cntAtom == 0)
    {
        return 0;
    }
    else
    {
        //index the atoms to be copied to pRes->atom
        Array1D <int> idx_1darray(cntAtom);
        int *idx = idx_1darray.array1D;
        Array1D <bool> isUsed_1darray(cntAtom);
        bool *isUsed = isUsed_1darray.array1D;
        int numAtom = 0;
        for(i = 0 ; i < cntAtom ; i++)
        {
            isUsed[i] = false;
            idx[i] = i;
        }
        if(!isGetAllAtomLocation)//if remove multiple alternative locations
        {
            Array1D <Atom> altAtom_1darray(cntAtom);
            Atom *altAtom = altAtom_1darray.array1D;
            Array1D <int> idxAltAtom_1darray(cntAtom);
            int *idxAltAtom = idxAltAtom_1darray.array1D;
            for(i = 0 ; i < cntAtom ; i++) 
            {
                InitAtom(&altAtom[i]);
                idxAltAtom[i] = i;
            }
            int cntAltAtom = 0 ;

            numAtom = 0;
            for(i = 0; i < cntAtom ; i++)
            {
                if(isUsed[i])
                    continue;
                else if(atom[i].altLoc == ' ')
                {
                    idx[numAtom] = i;
                    isUsed[i] = true;
                    numAtom ++;
                }
                else // if the alternative location for this atom exist
                {
                    cntAltAtom  = 0;
                    for(j = i; j < cntAtom; j++)
                    {
                        if(isUsed[j]) continue;// if this atom has already been used, continue;
                        if(strcasecmp(atom[j].name , atom[i].name) == 0
                                && atom[j].altLoc != ' '
                                && atom[j].altLoc != atom[i].altLoc) // if it is the alternative location of the same atom
                        {
                            CopyAtom(&altAtom[cntAltAtom], &atom[j] );
                            isUsed[j] = true;
                            idxAltAtom[cntAltAtom] = j;
                            cntAltAtom ++;
                        }
                    }
                    int indexSelAltAtom = GetAltAtomIndex(altAtom, cntAltAtom);
                    idx[numAtom] = idxAltAtom[indexSelAltAtom];
                    numAtom ++;
                }
            }
        }
        else
        {
            numAtom = cntAtom;
        }

        pRes->atom = new Atom[numAtom];
        for(i = 0 ; i < numAtom ; i++)
        {
            InitAtom( & (pRes->atom[i]));
            CopyAtom(&(pRes->atom[i]),&(atom[idx[i]]));
        }

        pRes->numAtom = numAtom;
        return numAtom;
    }
}
/*}}}*/

int SelectAltLocAtom(Atom* totalAtom, int &cntTotalAtom, FILE* fpPDBFile, char **metalEleList /*= NULL */, int numMetalEle /* = 0*/, bool isSelectMetal /* = false*/)     /*{{{*/
    //*************************************************************
    // SelectAltLocAtom()
    // if there are alternative locations for atom exist
    // read in the rest of atoms in this residue and select one atom position from
    // all alternative locations, return cntAtom, 
    // cntAtom will be updated after reading in atoms with alternative locations in
    // this residue, the file position of fpPDBFile will point to the next residue
    // 2007-04-21
    // when isSelectMetal is true, read in only the metal atom information to
    // totalAtom
    //**************************************************************
{
    char   recordID[SIZE_RECORD_ID+1] = "";
    fpos_t pos;
    int status_fpos = 0;

    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D; 

    int i,j;
    bool isEndMDLreached = false;
    int atomSerial = 0;

    int cntAtom = 0; // count of atoms within this residue
    Array1D <Atom> atom_1darray(MAX_ATOM_PER_RES_WITH_ALTATOM);
    Atom *atom = atom_1darray.array1D;
    for ( i = 0 ; i < MAX_ATOM_PER_RES_WITH_ALTATOM; i ++) InitAtom(&atom[i]);
    CopyAtom(&atom[0],&totalAtom[cntTotalAtom]); 
    cntAtom ++;


    fgetpos(fpPDBFile,&pos);
    // read in all the following atom with in residue
    while((linesize = fgetline(fpPDBFile, line, maxline))!= EOF)
    {
        if( linesize <= 0) continue;
        else if( sscanf(line,"%6s",recordID) < 1 ) continue;
        else if( strcmp(recordID,"ANISOU")== 0 ) { continue; } // neglect the "ANISOU" record
        else if( strcmp(recordID,"ENDMDL") == 0)
        {
            isEndMDLreached = true;
            fprintf(stderr,"End of model, take the first model only\n");
            if(atomSerial == MAX_ATOM_SERIAL)
                fprintf(stderr,"atom Serial is possibly overflow (> MAX_atom_SERIAL = %d), CHECK!\n", MAX_ATOM_SERIAL);
        }
        else if(strcmp(recordID,"ATOM") == 0 || strcasecmp(recordID,"HETATM") == 0) // if the line atom record
        {
            ScanfCoorRecord_Atom(line,&(atom[cntAtom]));
            atomSerial = atom[cntAtom].serial;
            if(!isSelectMetal || IsMetalAtom(&atom[cntAtom], metalEleList, numMetalEle))
            {
                if(atom[cntAtom].chainID == atom[0].chainID
                        && atom[cntAtom].resSeq == atom[0].resSeq
                        && atom[cntAtom].iCode  == atom[0].iCode ) // if still in the same residue

                {
                    cntAtom ++;
                    status_fpos = fgetpos(fpPDBFile, &pos); //retrieve the file position before the next getline
                    assert (status_fpos == 0);
                }
                else
                {
                    status_fpos = fsetpos(fpPDBFile,&pos); // set back the file pos to the previous position, if the current line is already in another residue
                    assert (status_fpos == 0);
                    break;
                }
            }
        }
        else
        {
            status_fpos = fsetpos(fpPDBFile,&pos); // set back the file pos to the previous position, if the current line is another type of recordID, (of course not in the same residue )
            assert (status_fpos == 0);
            break;
        }
    }

    Array1D <int> idx_1darray(cntAtom);
    int *idx = idx_1darray.array1D;
    Array1D <bool> isUsed_1darray(cntAtom);
    bool *isUsed = isUsed_1darray.array1D;
    int numAtom = 0;
    for(i = 0 ; i < cntAtom ; i++)
    {
        isUsed[i] = false;
        idx[i] = i;
    }
    Array1D <Atom> altAtom_1darray(cntAtom);
    Atom *altAtom = altAtom_1darray.array1D;
    Array1D <int> idxAltAtom_1darray(cntAtom);
    int *idxAltAtom = idxAltAtom_1darray.array1D;
    for(i = 0 ; i < cntAtom ; i++) 
    {
        InitAtom(&altAtom[i]);
        idxAltAtom[i] = i;
    }
    int cntAltAtom = 0 ;

    numAtom = 0;
    for(i = 0; i < cntAtom ; i++)
    {
        if(isUsed[i])
            continue;
        else if(atom[i].altLoc == ' ')
        {
            idx[numAtom] = i;
            isUsed[i] = true;
            numAtom ++;
        }
        else // if the alternative location for this atom exist
        {
            cntAltAtom  = 0;
            for(j = i; j < cntAtom; j++)
            {
                if(isUsed[j] == true) 
                    continue;// if this atom has already been used, continue;
                else if(strcasecmp(atom[j].name , atom[i].name) == 0
                        && atom[j].altLoc != ' '
                       ) // if it is the alternative location of the same atom
                {
                    CopyAtom(&altAtom[cntAltAtom], &atom[j] );
                    isUsed[j] = true;
                    idxAltAtom[cntAltAtom] = j;
                    cntAltAtom ++;
                }
            }
            int indexSelAltAtom = GetAltAtomIndex(altAtom, cntAltAtom);
            idx[numAtom] = idxAltAtom[indexSelAltAtom];
            numAtom ++;
        }
    }

    for(i = 0 ; i < numAtom ; i++)
    {
        CopyAtom(&(totalAtom[cntTotalAtom]),&(atom[idx[i]]));
        cntTotalAtom ++;
    }

    return cntTotalAtom;
}/*}}}*/

template <class T> int GetMODM(const char* modmfilepath, T ** M /*= NULL*/, char* alphabetMODM /*= NULL*/, char* aaSeq /*= NULL*/,double* score1 /*= NULL*/,double* score2 /*= NULL*/, double *parameter /*=NULL*/, int *seqIndex/*=NULL*/, int startIndex/*=1*/)/*{{{*/
/*****************************************************************************
 * GetMODM 	 read in aaSeq also                                      
 * get modification matrix which is derived from PSI-BLAST PSSM         
 * for an menu line in modm file like this                              
 * Num AA     A   R   N   D   C   Q   E   G   H   I   L   K   W   Y   V 
 * numTag = 2, then followed by the matrix value                        
 *
 * seqIndex is used when the sequence in the MODM file is only partly recorded
 * startIndex is used for specifying which is the start index for the sequence, default is 1
 ****************************************************************************/
{
    int   i;
    int   cnt;
    int linesize;
    int maxline = 500;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D; 
    FILE *fpMODM;
    char *pch;
    int numTag; // number of coloums before the matrix values
    int sizeAlphabet;

    fpMODM = fopen(modmfilepath,"r");
    checkfilestream(fpMODM,modmfilepath,"r");

    f_neglect_comment(fpMODM,'#');

    numTag = 2;
    i = 0;
    linesize = fgetline(fpMODM,line,maxline);
    pch = strtok (line,WHITE_SPACE);
    while (pch != NULL)
    {
        if(i >= numTag) alphabetMODM[i-numTag] = *pch;
        i++;
        pch = strtok (NULL, WHITE_SPACE);
    }
    alphabetMODM[i-numTag] = '\0';
    sizeAlphabet = strlen(alphabetMODM);
    char delim[] = WHITE_SPACE;
    char first_non_blank_char = ' ';
    char tmpstr[100+1] = "";
    //LOG: 2006-02-03 12:56:10 Friday  Week 05 <nanjiang@shu>
    // expand the alphabetMODM for some unexpected characters, say "BZX*"
    //     ExpandAlphabet_AA(alphabetMODM);
    cnt = 0;
    while((linesize = fgetline(fpMODM, line ,maxline))!= EOF)
    {
        if(sscanf(line," %c", &first_non_blank_char) != 1) continue;
        if(isdigit(first_non_blank_char))
        {
            pch = strtok (line,delim);
            int j = 0 ;
            while (pch != NULL)
            {
                if(j == 0 && seqIndex != NULL)       {seqIndex[cnt] = atoi(pch); seqIndex[cnt] -= startIndex;}
                else if(j == 1 && aaSeq != NULL)      aaSeq[cnt] = *pch;
                else if(j >= numTag && j < numTag+sizeAlphabet && M != NULL)  M[cnt][j-numTag] = T(atof(pch));
                else if(j == numTag +sizeAlphabet   && score1 != NULL)        score1[cnt] = atof(pch);
                else if(j == numTag +sizeAlphabet+1 && score2 != NULL)        score2[cnt] = atof(pch);
                pch = strtok (NULL, delim);
                j++;
            }
            cnt ++;
            if( cnt >= LONGEST_SEQ) break;
        }
        else if (first_non_blank_char == 'K')
        {
            if(parameter == NULL) break;
            int j = 0;
            while((linesize = fgetline(fpMODM, line, maxline)) != EOF)
            {
                sscanf(line, "%s %s %lf %lf", tmpstr, tmpstr, &parameter[j], &parameter[j+1]);
                j += 2;
                if(j >= 8) break;
            }
            break;
        }
    }
    if(aaSeq != NULL) aaSeq[cnt] = '\0';
    fclose(fpMODM);	

    if(cnt >= LONGEST_SEQ ) return (-1);
    else return cnt;
}
template  int GetMODM <int>    (const char* modmfilepath, int **M /*= NULL*/   , char* alphabetMODM /*= NULL*/, char* aaSeq /*= NULL*/, double* score1 /*= NULL*/, double* score2 /*= NULL*/, double *parameter /* = NULL*/, int *seqIndex /*=NULL*/, int startIndex/*=1*/);
template  int GetMODM <float>  (const char* modmfilepath, float **M /*= NULL*/ , char* alphabetMODM /*= NULL*/, char* aaSeq /*= NULL*/, double* score1 /*= NULL*/, double* score2 /*= NULL*/, double *parameter /* = NULL*/, int *seqIndex /*=NULL*/, int startIndex/*=1*/);
template  int GetMODM <double> (const char* modmfilepath, double **M /*= NULL*/, char* alphabetMODM /*= NULL*/, char* aaSeq /*= NULL*/, double* score1 /*= NULL*/, double* score2 /*= NULL*/, double *parameter /* = NULL*/, int *seqIndex /*=NULL*/, int startIndex/*=1*/);

/*}}}*/
template <class T> void CalLogM(T **M, T **log_M, int xSize, int ySize,  double roundoff_scale /* = 1.0*/)/*{{{*/
/*****************************************************************************
 * calculate log(Mij) from Mij
 * roundoff_scale should set to a larger number, say, 100.0, if T is int
 ****************************************************************************/
{
    int i,j;
    assert(M != NULL);
    assert(log_M != NULL);
    double number_close_zero = 1e-6;

    for(i = 0 ; i < xSize; i++)
    {
        for(j = 0; j < ySize; j++)
        {
#ifdef ASSERT
            assert(M[i][j] < 0)
#endif
            if(double(M[i][j]) < number_close_zero)
                log_M[i][j] = T (log(number_close_zero)*roundoff_scale);
            else
                log_M[i][j] = T (log(double(M[i][j])) * roundoff_scale);
        }
    }
}

template  void CalLogM <int>    (int   **M , int    **log_M, int xSize, int ySize, double roundoff_scale /*= 1.0*/);
template  void CalLogM <float>  (float **M , float  **log_M, int xSize, int ySize, double roundoff_scale /*= 1.0*/);
template  void CalLogM <double> (double **M, double **log_M, int xSize, int ySize, double roundoff_scale /*= 1.0*/);
/*}}}*/
int GetTupingQij(const char *infile, int **M, char *aaSeq, char *shapeSeq, char *dsspSeq, int *waterAcc, double *score1, double *score2, int *seqIndex, int startIndex/*=1*/, const int *mapIndex/*=MAPINDEX_BLOSUM_Tuping*/)/*{{{*/
/*****************************************************************************
 * 2007-05-31
 * Read in the matrix information for the modm file output with tuping's format
 ****************************************************************************/
{
    FILE *fpin = fopen(infile, "r");
    checkfilestream(fpin, infile, "r");
    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;
    char first_non_black_char = ' ';
    int cnt = 0;
    char tmpstr[100] = "";
    
    char delim[] = WHITE_SPACE;

    while((linesize = fgetline(fpin, line, maxline)) != EOF)
    {
       assert( sscanf(line," %c", &first_non_black_char) == 1);
       if(isdigit(first_non_black_char))
       {
           char *pch;
           pch = strtok(line, delim);
           int j = 0 ;
           while(pch != NULL)
           {
               if (j == 0) 
               {
                   if(seqIndex != NULL)
                   {
                       seqIndex[cnt] = atoi(pch);
                       seqIndex[cnt] -= startIndex;
                   }
               }
               else if(j == 1) 
               {
                   if(aaSeq != NULL) aaSeq[cnt] = *pch;
               }
               else if(j == 2) 
               {
                   if(shapeSeq != NULL) shapeSeq[cnt] = pch[0];
                   if(waterAcc != NULL) 
                   {
                       tmpstr[0] = *(pch+1);
                       tmpstr[1] = '\0';
                       waterAcc[cnt] = atoi(tmpstr);
                   }
                   if(dsspSeq != NULL) dsspSeq[cnt] = pch[2];
               }
               else if (j >= 3 && j < 23)
               {
                   if(M != NULL) M[cnt][mapIndex[j-3]]=atoi(pch);
               }
               else if( j == 23)
               {
                   if(score1 != NULL) score1[cnt] = atof(pch);
               }
               else if( j == 24)
               {
                   if(score2 != NULL) score2 [cnt] = atof(pch);
               }

               pch = strtok(NULL, delim);
               j ++;
           }

           cnt ++;
       }
    }
    fclose(fpin);
    int length = cnt;
    if(aaSeq != NULL) aaSeq[length] = '\0';
    if(shapeSeq != NULL) shapeSeq[length] = '\0';
    if(dsspSeq != NULL) dsspSeq[length] = '\0';

    return length;
}
/*}}}*/
//int GetMODM(const char* modmfilepath, int** M, char* alphabetMODM,char* aaSeq = NULL,double* score1 = NULL,double* score2 = NULL)[>{{{<]
//{
//    int   i;
//    int   cnt;
//    char  tmpstr[100];
//    int   maxline    = 200;
//    char *line        = new char[maxline+1];
//    FILE *fpMODM;
//    int numTag; // number of coloums before the matrix value 
//    // stringstream ss(stringstream::in | stringstream::out);

//    if((fpMODM = fopen(modmfilepath,"r")) == 0)
//    {
//        printf("Can not open MODM file %s\n", modmfilepath);
//        assert( fpMODM != NULL);
//    }

//    while(fgetline(fpMODM,line, maxline ) != EOF)
//    {
//        if(line[0] != '#')
//            break;
//    }
//    numTag = 2;


//    // ss.str(line);
//    for(i = 0 ; i < numTag; i++) ss >> tmpstr;
//    ss.getline(line,linesize);
//    if(ss.bad() &&  ! ss.eof())
//    {
//        printf("stringstream fatal error\n");
//        assert(!(ss.bad() && !ss.eof()));
//    }
//    ss.clear(); 
//    SpanExcluding(line,alphabetMODM);
//    //LOG: 2006-02-03 12:56:10 Friday  Week 05 <nanjiang@shu>
//    // expand the alphabetMODM for some unexpected characters, say "BZX*"
//    //     ExpandAlphabet_AA(alphabetMODM);
//    cnt = 0;
//    while(fgetline(fpMODM, line ,linesize) != EOF)
//    {
//        ss.str(line);
//        for(i = 0 ; i < numTag  ; i++) 
//        {
//            ss >> tmpstr;
//            if(i == 1 && aaSeq != NULL)
//                aaSeq[cnt] = tmpstr[0];
//        }
//        for(i = 0 ; i < 20 ; i++) { ss >> M[cnt][i]; }
//        if(score1 != NULL) ss >> score1[cnt];
//        if(score2 != NULL) ss >> score2[cnt];

//        if(ss.bad() &&  ! ss.eof())
//        {
//            printf("stringstream fatal error\n");
//            assert(!(ss.bad() && !ss.eof()));
//        }
//        ss.clear(); 

//        cnt ++;
//    }
//    if(aaSeq != NULL) aaSeq[cnt] = '\0';
//    fclose(fpMODM);	
//    delete [] line;
//    return cnt;
//}[>}}}<]
int GetPSSM(const char *infile, int &length, char *aaSeq /*=NULL*/, int **Mij/*=NULL*/, int **fij/*=NULL*/, double *score1/*=NULL*/, double *score2/*=NULL*/, double *parameter/*=NULL*/, int *seqIndex/*=NULL*/)/*{{{*/
    // set variable to NULL for not reading in the data for that variable, e.g.
    // GetPSSM(pssmfile, length, aaSeq, NULL, fij, NULL, NULL, NULL);
    // reads in only fij matrix, if the last several items will be set to NULL, it
    // can be ignored when calling
    //
    // seqIndex is the first column of pssm file, which should be always from 1 to length
    // only for checking purpose
{
    FILE *fpin;
    fpin = fopen(infile, "r");
    checkfilestream(fpin, infile, "r");

    int linesize;
    int maxline = 500;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D; 

    int i;
    char first_non_blank_char;
    char *pch;
    char delim[] = WHITE_SPACE;
    char tmpstr[300];
    i =0 ;
    while((linesize = fgetline(fpin, line, maxline)) != EOF)
    {
        if(sscanf(line," %c", &first_non_blank_char) != 1) continue;
        if(isdigit(first_non_blank_char))
        {
            //int tmpnum;
            pch = strtok(line,delim);
            int j = 0;
            while(pch != NULL)
            {
                if(j == 0 && seqIndex != NULL){
                    seqIndex[i] = atoi(pch);
                } else if(j == 1 && aaSeq != NULL) {
                    aaSeq[i]    = pch[0];
                } else if(j >= 2  && j < 22 && Mij != NULL) {
                    Mij[i][j-2] = atoi(pch);
                } else if(j >= 22 && j < 42 && fij != NULL) {
                    fij[i][j-22]= atoi(pch);
                } else if(j == 42 && score1 != NULL){
                    score1[i]   = atof(pch);
                } else if(j == 43 && score2 != NULL){
                    if(strcmp(pch, "inf") == 0){
                        score2[i]=0.0;
                    }else{
                        score2[i]   = atof(pch);
                    }
                }
                pch = strtok(NULL, delim);
                j ++;
            }
            i ++;
            if(i >= LONGEST_SEQ) break;
        }
        else if (first_non_blank_char == 'K')
        {
            if(parameter == NULL) break;
            int j = 0;
            while((linesize = fgetline(fpin, line, maxline)) != EOF)
            {
                sscanf(line, "%s %s %lf %lf", tmpstr, tmpstr, &parameter[j], &parameter[j+1]);
                j += 2;
                if(j >= 8) break;
            }
            break;
        }
    }
    fclose(fpin);

    length = i;
    if (length >= LONGEST_SEQ) return -1;
    else return length;
}
/*}}}*/
int ReadInProfile(const char *file, int **M, char *aaSeq /*= NULL*/, char *shapeSeq /*= NULL*/, int *waterAcc /*= NULL*/, char *dsspSec /*= NULL*/, float *score1 /*= NULL*/, float *score2 /*= NULL*/, int *sumProfile /*= NULL*/, int maxLength /*= LONGEST_SEQ*/)/*{{{*/
/*****************************************************************************
 * Read in profile matrix. 
 *
 * This is suitable for read in data from Qij, MODM and fragAcc matrices with
 * SAD (shape-waterAcc-dsspSec) column
 *
 * Note that: The aaSeqIndex is always read from the first column in the file, and thus
 * not all items in the aaSeq are filled. 
 *
 * This function is different from GetTupingQij in that 
 * 1. in GetTupingQij, aaSeq and other arrays are always stored from 0 and fully filled
 *    an extra array seqIndex is kept to record the aaSeqIndex of each row
 * 2. mapIndex can be supplied in case that the alphabet order in the matrix
 *    file is different from the required. However, in this function, the
 *    alphabet AAAlphabet_Tuping is assumed 
 *
 * return the number of residles read in 
 * return -1 if the file fails to open
 * 2008-01-18, Nanjiang
 ****************************************************************************/
{
    int j = 0;
    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;
    char first_non_black_char = ' ';
    int cnt = 0;
    int aaSeqIndex = 0;
    float s1;
    float s2;
    char aa;
    char strSAD[10] = ""; /*string for storing shape acc and dsspSec*/
    char tmpstr[10] = "";
    int profile[NUM_AA];

    FILE *fpin = NULL;
    fpin = fopen(file,"r");
    if (  fpin != NULL  )/*{{{*/
    {
        while((linesize = fgetline(fpin, line ,maxline)) != EOF) //read in matrix/*{{{*/
        {
            assert( sscanf(line," %c", &first_non_black_char) == 1);
            if(isdigit(first_non_black_char))
            {
                sscanf(line,"%d %c %s %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %f %f",/*{{{*/
                        &aaSeqIndex,&aa,strSAD,
                        &profile[0],
                        &profile[1],
                        &profile[2],
                        &profile[3],
                        &profile[4],
                        &profile[5],
                        &profile[6],
                        &profile[7],
                        &profile[8],
                        &profile[9],
                        &profile[10],
                        &profile[11],
                        &profile[12],
                        &profile[13],
                        &profile[14],
                        &profile[15],
                        &profile[16],
                        &profile[17],
                        &profile[18],
                        &profile[19],
                        &s1,&s2);/*}}}*/
                aaSeqIndex -- ; /*aaSeqIndex in the file starts from 1*/
                if (   (aaSeqIndex<0) || (aaSeqIndex>=LONGEST_SEQ)   )//for the first two lines
                { continue; }

                if (aaSeq != NULL) { aaSeq[aaSeqIndex] = aa;}
                if (shapeSeq != NULL) { shapeSeq[aaSeqIndex] = strSAD[0];}
                if (waterAcc != NULL) 
                { 
                    tmpstr[0] = strSAD[1]; tmpstr[1] = '\0';
                    waterAcc[aaSeqIndex] = atoi(tmpstr);
                } /*bug fixed, should be converted into integer*/
                if (dsspSec != NULL) { dsspSec[aaSeqIndex] = strSAD[2];}
                for (j=0; j<NUM_20_AA; j++)
                {    
                    M[aaSeqIndex][j] =  profile[j];
                    if (sumProfile != NULL ){sumProfile[aaSeqIndex] += M[aaSeqIndex][j];}
                }
                if (score1 != NULL) {score1[aaSeqIndex] = s1;}
                if (score2 != NULL) {score2[aaSeqIndex] = s2;}
                cnt ++;
            }
        }/*}}}*/
        fclose(fpin);
        return cnt;
    }/*}}}*/
    else
    {
        return -1;
    }
}
/*}}}*/
int ReadInProfile_fp(FILE* fpin, unsigned long readsize, int **M, char *aaSeq /*= NULL*/, char *shapeSeq /*= NULL*/, int *waterAcc /*= NULL*/, char *dsspSec /*= NULL*/, float *score1 /*= NULL*/, float *score2 /*= NULL*/, int *sumProfile /*= NULL*/, int maxLength /*= LONGEST_SEQ*/) /*{{{*/
/*****************************************************************************
 * Read in profile matrix by given FILE handler, 
 * created 2011-10-13. 
 *
 * This is suitable for read in data from Qij, MODM and fragAcc matrices with
 * SAD (shape-waterAcc-dsspSec) column
 *
 * Note that: The aaSeqIndex is always read from the first column in the file, and thus
 * not all items in the aaSeq are filled. 
 *
 * This function is different from GetTupingQij in that 
 * 1. in GetTupingQij, aaSeq and other arrays are always stored from 0 and fully filled
 *    an extra array seqIndex is kept to record the aaSeqIndex of each row
 * 2. mapIndex can be supplied in case that the alphabet order in the matrix
 *    file is different from the required. However, in this function, the
 *    alphabet AAAlphabet_Tuping is assumed 
 *
 * return the number of residles read in 
 * return -1 if errors are encountered
 ****************************************************************************/
{
    if ( fpin != NULL  ) {
        Array1D <char> buffer_1darray(readsize+1);
        char *buffer = buffer_1darray.array1D;
        size_t status_fread =  fread (buffer,1,readsize,fpin);
        buffer[readsize] = '\0';
        if (status_fread != readsize){
            fprintf(stderr,"fread failed in the function ReadInProfile_fp\n");
            return -1;
        }else{
            int j = 0;
            char first_non_black_char = ' ';
            int cnt = 0;
            int aaSeqIndex = 0;
            float s1;
            float s2;
            char aa;
            char strSAD[10] = ""; /*string for storing shape acc and dsspSec*/
            char tmpstr[10] = "";
            int profile[NUM_AA];

            char *pch=NULL;
            char delim[]     = "\n";
            pch = strtok (buffer,delim);
            while (pch != NULL) {
                int status_sscanf = sscanf(pch," %c", &first_non_black_char);
                if(status_sscanf == 1 && isdigit(first_non_black_char)) {
                    sscanf(buffer,"%d %c %s %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %f %f",/*{{{*/
                            &aaSeqIndex,&aa,strSAD,
                            &profile[0],
                            &profile[1],
                            &profile[2],
                            &profile[3],
                            &profile[4],
                            &profile[5],
                            &profile[6],
                            &profile[7],
                            &profile[8],
                            &profile[9],
                            &profile[10],
                            &profile[11],
                            &profile[12],
                            &profile[13],
                            &profile[14],
                            &profile[15],
                            &profile[16],
                            &profile[17],
                            &profile[18],
                            &profile[19],
                            &s1,&s2);/*}}}*/
                    aaSeqIndex -- ; /*aaSeqIndex in the file starts from 1*/
                    if (   (aaSeqIndex<0) || (aaSeqIndex>=LONGEST_SEQ)   )//for the first two lines
                    { continue; }

                    if (aaSeq != NULL) { aaSeq[aaSeqIndex] = aa;}
                    if (shapeSeq != NULL) { shapeSeq[aaSeqIndex] = strSAD[0];}
                    if (waterAcc != NULL) 
                    { 
                        tmpstr[0] = strSAD[1]; tmpstr[1] = '\0';
                        waterAcc[aaSeqIndex] = atoi(tmpstr);
                    } /*bug fixed, should be converted into integer*/
                    if (dsspSec != NULL) { dsspSec[aaSeqIndex] = strSAD[2];}
                    for (j=0; j<NUM_20_AA; j++)
                    {    
                        M[aaSeqIndex][j] =  profile[j];
                        if (sumProfile != NULL ){sumProfile[aaSeqIndex] += M[aaSeqIndex][j];}
                    }
                    if (score1 != NULL) {score1[aaSeqIndex] = s1;}
                    if (score2 != NULL) {score2[aaSeqIndex] = s2;}
                    cnt ++;
                }
                pch = strtok (NULL, delim);
            }
            return cnt;
        }
    } else {
        return -1;
    }
}
/*}}}*/

template <class T> T GetIntegConsv(char aa, T *V, char* alphabetMODM, int type_modm /*= MODM_PER*/)/*{{{*/
{
    T consv;
    int sizeAlphabetMODM = strlen(alphabetMODM);
    int daa;
    int daaH;
    int daaC;

    daaH = Char2Digit('H', alphabetMODM,sizeAlphabetMODM);
    daaC = Char2Digit('C', alphabetMODM,sizeAlphabetMODM);
    assert( daaH >= 0);
    assert( daaC >= 0);

    if((daa = Char2Digit(aa,alphabetMODM,sizeAlphabetMODM)) >= 0)
        consv = V[daa];
    else 
        return T(INIT_CONSV);

    if( type_modm == MODM_PER)
    {
        if(aa == 'H')
        {
            consv= consv + V[daaC];
        }
        else if(aa == 'C')
        {
            consv = consv + V[daaH];
        }
        else if(aa == 'D' || aa == 'E')
        {
            consv = consv + V[daaC] + V[daaH];
        }
        else
        { }
    }
    else if ( type_modm == MODM_LOG)
    {
        if(aa == 'H')
        {
            consv = T(log(exp(double(consv))+ exp(double(V[daaC]))));
        }
        else if(aa == 'C')
        {
            consv = T(log(exp(double(consv))+ exp(double(V[daaH]))));
        }
        else if(aa == 'D' || aa == 'E')
        {
            consv = T(log(exp(double(consv))+  exp(double(V[daaC])) + exp(double(V[daaH]))));
        }
        else
        { }
    }

    return consv;
}
template int    GetIntegConsv  <int>   (char aa, int *V   , char *alphabetMODM, int type_modm /* = MODM_PER */);
template float  GetIntegConsv  <float> (char aa, float *V , char *alphabetMODM, int type_modm /* = MODM_PER */);
template double GetIntegConsv  <double>(char aa, double *V, char *alphabetMODM, int type_modm /* = MODM_PER */);
/*}}}*/

//int GetEncodeAA(char aa)[>{{{<]
//{
//    if(aa == 'C')
//    { return 0; }
//    else if (aa == 'H')
//    { return 1; }
//    else if (aa == 'D')
//    { return 2; }
//    else if (aa == 'E')
//    { return 3; }
//    else
//    { return 4; }
//}
//[>}}}<]
int GetEncodeScore1(double score1)/*{{{*/
{
    for(int i = 0 ; i < SIZE_ENCODE_SCORE1-1; i ++)
    {
        if(score1 < EncodeScaleScore1[i])
        {
            return i;
        }
    }
    return SIZE_ENCODE_SCORE1-1;
}
/*}}}*/
int GetEncodeScore2(double score2)/*{{{*/
{
    for(int i = 0 ; i < SIZE_ENCODE_SCORE2-1; i ++)
    {
        if(score2 < EncodeScaleScore2[i])
        {
            return i;
        }
    }
    return SIZE_ENCODE_SCORE2-1;
}
/*}}}*/
int GetEncodeWaterAcc(int waterAcc)/*{{{*/
{
    for(int i = 0 ; i < SIZE_ENCODE_WATERACC-1; i ++)
    {
        if(waterAcc <= EncodeScaleWaterAcc[i])
        {
            return i;
        }
    }
    return SIZE_ENCODE_WATERACC-1;
}
/*}}}*/
int GetEncodeHydrophobicity(double hydrophobicity)/*{{{*/
{
    for(int i = 0 ; i < SIZE_ENCODE_HYDROPHOBICITY-1; i ++)
    {
        if(hydrophobicity < EncodeScaleHydrophobicity[i])
        {
            return i;
        }
    }
    return SIZE_ENCODE_HYDROPHOBICITY-1;
}
/*}}}*/
int GetEncodeDistance(int dist)/*{{{*/
{
    for(int i = 0 ; i < SIZE_ENCODE_DISTANCE-1; i ++)
    {
        if(dist <= EncodeScaleDistance[i])
        {
            return i;
        }
    }
    return SIZE_ENCODE_DISTANCE-1;
}
/*}}}*/

/* boolean functions*/
bool IsNonMetalHetGroup(const char *resName)/*{{{*/
{
    int i = BinarySearch_String(resName, nonMetalHetGroup, numNonMetalHetGroup);
    if (i != -1)
        return true;
    else return 
        false;
}
/*}}}*/
bool IsDNASeq(const char *aaSeq, int n /*= 0*/)/*{{{*/
    //check if the amino acid sequence read from SEQRES record is DNA sequence
{
    if(n == 0) n = strlen(aaSeq);
    int cnt = 0 ;
    int i ;
    for(i = 0 ; i < n  ; i ++)
    {
        if(aaSeq[i] == CHAR_NON_RESIDUE)
            cnt ++;
    }
    if(float(cnt)/n > 0.5)
        return true;
    else
        return false;
}/*}}}*/
bool IsMetalAtom(Atom *pAtom, char **metalEleList, int numMetalEle)/*{{{*/
    /*****************************************************************************
     * IsMetalAtom()
     * Check if the ATOM record in PDB file is a metal atom, notice that the
     * atomName must be the origAtomName, that is, having all four characters.
     *
     * copy the metal element name to pAtom->element if it is confirmed as a metal
     * atom
     ****************************************************************************/
{ 
    if(strcmp(pAtom->element, "") != 0 && !IsDigit(pAtom->element))
    {
        if(strcmp(pAtom->element,"G") == 0) // wrong annotation of HG to G, for example in pdbcode 1E8B
        {
            fprintf(stderr, "Warning! atom.element name invalid, 'G', possible 'Hg'\n");
        }

        if(BinarySearch_String(pAtom->element, metalEleList, numMetalEle) != -1)
            return true;
        else
            return false;
    }
    else// pAtom->element does not have element information
    {
        if(pAtom->origName[0] == ' ') // the metal atom is always of the format 'FE1 '
            return false ;     // the first position is NOT a whitespace

        if(IsNonMetalHetGroup(pAtom->resName) )
            return false;

        int indexEle;
        char atomName[SIZE_ATOM_ELEMENT+1] = "";
        my_strcpy(atomName,pAtom->origName,SIZE_ATOM_ELEMENT);
        char *pch = &atomName[0];
        while(*pch)
        {
            if(*pch < 'A' || *pch > 'Z')
            {
                *pch = '\0';
                break;
            }
            pch ++;
        }
        if((indexEle = BinarySearch_String(atomName, metalEleList, numMetalEle))!= -1)
        {
            my_strcpy(pAtom->element, metalEleList[indexEle], SIZE_ATOM_ELEMENT); // copy the element to pAtom-> element
            return true;
        }
        else
            return false;
    }
}/*}}}*/
bool IsInContactAtomList(const char *atomOrigName)/*{{{*/
    /*****************************************************************************
     * atoms which can bind to metal
     * +--------------------------------------------------------------------+
     * |LOG: 2006-03-17 11:33:30 Friday  Week 11                            |
     * |     excluding carbon atoms and atoms in the backbone might be bette|
     * |     return false if carbon atoms or atoms in the backbone          |
     * |     else return true                                               |
     * +--------------------------------------------------------------------+
     * exclude carbon atoms, hydrogen atoms and backbone atoms
     * that is, assume the above atoms could not bond to metal atoms
     * there are some proteins structure including the coordinate of hydrogen
     * atoms, for example, "1A5T"
     * 
    //LOG: 2006-06-16 11:13:38 Friday  Week 24 <nanjiang@casio>
    //  backbone oxygen and nitrogen can bind to metal ions, for example
    //  1JIWI, ZN bind to backbone oxygen and nitrogen of SER1I
    //  1KHOA, ZN bind to backbone oxygen and nitrogen of TRP1A
    //  1LM4A, FE bind to backbone nitrogen of GLY-10A
    //
    //  so the contact Atom List is 'N', 'O', 'S' for standard amino acid
     ****************************************************************************/
{
    const char *atomList[] = {"N","O","S"};
    int numAtomList = 3;
    char eleName[SIZE_ATOM_ELEMENT+1] = "";
    sprintf(eleName,"%c", atomOrigName[1]);
    if(BinarySearch_String((const char*)eleName, atomList, numAtomList) != -1)
        return true;
    else
        return false;
    /*{{{*/
    /*
     *excluding backbone atoms, carbon and hydrogen atoms
     *char str[100] = "";
     *sscanf(atomOrigName,"%s",str);
     *if( (atomOrigName[1] == 'C' )      // carbon
     *       || atomOrigName[1] == 'H'  // hydrogen 
     *       || (strcmp(str,"O") == 0)  // backbone oxygen
     *       || (strcmp(str,"N") == 0)  // backbone nitrogen
     *       )
     *   return false;
     *else 
     *   return true;
     */
    /*}}}*/
}/*}}}*/
bool IsSSBonded(int aaSeqIndex, SSBondPro *ssbondPro)/*{{{*/
{
    int  i,j ; 
    for (i = 0; i < ssbondPro->numSSBond ; i++)
    {
        for(j = 0 ; j < ssbondPro->ssbond[i].numRes; j ++)
        {
            if(ssbondPro->ssbond[i].res[j].aaSeqIndex == aaSeqIndex) 
                return true;
        }
    }
    return false;
}/*}}}*/
bool IsZnBound(int aaSeqIndex, MetalPro *znPro)/*{{{*/
{
    int  i; 
    for (i = 0; i < znPro->numBoundRes ; i++)
    {
        if(aaSeqIndex == znPro->res[i].aaSeqIndex)
            return true;
    }
    return false;
}/*}}}*/
bool IsZnBound(int aaSeqIndex1, int aaSeqIndex2, MetalPro *znPro, bool operation)/*{{{*/
{
    int  i; 
    bool isResBound1 = false;
    bool isResBound2 = false;
    for (i = 0; i < znPro->numBoundRes ; i++)
    {
        if(aaSeqIndex1 == znPro->res[i].aaSeqIndex)
            isResBound1 = true;
        if(aaSeqIndex2 == znPro->res[i].aaSeqIndex)
            isResBound2 = true;
    }
    if(operation == AND)
        return (isResBound1 && isResBound2);
    else 
        return (isResBound1 || isResBound2);
}/*}}}*/
bool IsZnBound(int aaSeqIndex1, int aaSeqIndex2, MetalPro2 *znPro, bool operation)/*{{{*/
{
    int i, j;
    bool isResBound1 = false;
    bool isResBound2 = false;
    int aaSeqIndex;
    if(znPro->numMetalAtom < 1)
        return false;

    for (i = 0; i < znPro->numMetalAtom ; i++)
    {
        isResBound1 = isResBound2 = false;// only when both residues bond to the same atoms are considerd both bound
        for(j = 0; j < znPro->atomEnv[i].numRes ; j++)
        {
            aaSeqIndex =  znPro->atomEnv[i].res[j].aaSeqIndex;
            if(aaSeqIndex1 == aaSeqIndex)
            {
                isResBound1 = true;
            }
            else if(aaSeqIndex2 == aaSeqIndex)
            {
                isResBound2 = true;
            }
        }

        if(operation == AND)
        {
            if (isResBound1 && isResBound2)
                return true;
        }
        else 
        {
            if (isResBound1 || isResBound2)
                return true;
        }
    }

    return false;    
}/*}}}*/
bool IsNewResidue(Atom *pAtom, Residue *resGroup, int numRes)/*{{{*/
    /*******************************************************************
     * Check if the residue to which the atom belongs is a new residue *
     * for a group of residues                                         *
     *******************************************************************/
{
    int i;
    if(numRes == 0 ) return true;
    for ( i = 0 ; i < numRes; i++)
    {
        if(   pAtom->chainID == resGroup[i].chainID
                && pAtom->resSeq  == resGroup[i].resSeq 
                && pAtom->iCode   == resGroup[i].resICode )
            return false;
    }
    return true;
}/*}}}*/
bool IsNewResidue(Residue *pRes, Residue *resGroup, int numRes)/*{{{*/
    /*****************************************************************************
     * check if the residue is a new residue compared to a group of residues, or
     * check if the residue does not exist in the residue group
     ****************************************************************************/
{
    int i;
    if(numRes == 0 ) return true;
    for ( i = 0 ; i < numRes; i++)
    {
        if(   pRes->chainID  == resGroup[i].chainID
                && pRes->resSeq   == resGroup[i].resSeq
                && pRes->resICode == resGroup[i].resICode )
            return false;
    }
    return true;
}/*}}}*/
bool IsHEMPro(MODM *pMODM, bool *isPolyHis)/*{{{*/
{
    int i ;
    Array1D <Residue> HCRes_1darray(pMODM->length);
    Residue *HCRes = HCRes_1darray.array1D;
    char resList[] = "CHDE";
    double cutoff_score2 = 0.15;
    double cutoff_consv  = 9.0;
    char pattern[] = "HCCHHCCH";
    bool isUseConsvI = false;
    bool isMaskPolyHis = true;
    int numHCRes = 0;
    numHCRes = GetHCRes(pMODM, isPolyHis, HCRes, cutoff_score2, cutoff_consv, resList, isUseConsvI, isMaskPolyHis);
    Array1D <char> HCResSeq_1darray(numHCRes+2);
    char *HCResSeq = HCResSeq_1darray.array1D;
    for( i = 0 ; i < numHCRes ; i++)
    {
        HCResSeq[i] = HCRes[i].aa;
    }
    HCResSeq[i] = '\0';
    if(strstr(HCResSeq, pattern) != NULL)
        return true;
    else 
        return false;
}
/*}}}*/
bool IsSatisfyHCResRule(char* id, MODM *pMODM,  bool *isPolyHis,  /*{{{*/ 
        char* resList, double cutoff_consv, double cutoff_score2, int window, int min_numHCRes, 
        Residue* HCRes, int& numHCRes, /* return parameters*/ 
        Residue* LCRes, int& numLCRes, /* return parameters*/ 
        bool isMaskPolyHis, bool isUseConsvI
        )
/*****************************************************************************
 * HCRes -- high conservation level 
 *       (with consv >= cutoff_consv and score2 >= cutoff_score2) residues in resList
 * LCRes -- non high conservation level residues in resList and with 
 *          score2 >=  cutoff_score2
 ****************************************************************************/
{

    int i;
    DATATYPE_CONSV consv;
    int daa;         // digital amino acid
    int   numResList  = strlen(resList);
    bool isSatisfied1 = false;

    int      length       = pMODM->length;
    int      type_modm    = pMODM->type_modm;
    DATATYPE_MODM_MATRIX    **M         = pMODM->M;
    //double  *score1       = pMODM->score1;
    double  *score2       = pMODM->score2;
    char    *aaSeq        = pMODM->aaSeq;
    char    *alphabetMODM = pMODM->alphabetMODM;

    // bool isSatisfied2 = false;

    // int window1 = 100; //window for all residues
    // int window2 = 30;  // window for 2 residues

    int cntHCRes = 0;
    int cntLCRes = 0;
    for(i = 0 ; i < length ; i++)
    {
        if(IsInCharSet(aaSeq[i], resList, numResList) && score2[i] >= cutoff_score2 && 
                (!isPolyHis[i] || !isMaskPolyHis))
        {
            daa = Char2Digit(aaSeq[i],alphabetMODM);
            if(daa < 0) continue;
            if(!isUseConsvI)
                consv = M[i][daa];
            else
                consv = GetIntegConsv(aaSeq[i], M[i], alphabetMODM, type_modm);
            if(consv >= cutoff_consv )
            {
                HCRes[cntHCRes].consv      = consv;
                HCRes[cntHCRes].aaSeqIndex = i;
                HCRes[cntHCRes].aa         = aaSeq[i];
                cntHCRes ++;
            }
            else
            {
                LCRes[cntLCRes].consv      = consv;
                LCRes[cntLCRes].aaSeqIndex = i;
                LCRes[cntLCRes].aa         = aaSeq[i];
                cntLCRes ++;
            }
        }
    }
    numHCRes = cntHCRes;
    numLCRes = cntLCRes;

    int diffi = 0 ;
    if(cntHCRes < min_numHCRes)
    {
        isSatisfied1 = false;
        // isSatisfied2 = false;
    }
    else
    {
        for(i = 0; i < numHCRes - min_numHCRes; i++)
        {
            diffi =  HCRes[i+min_numHCRes-1].aaSeqIndex - HCRes[i].aaSeqIndex;
            if(diffi <= window)
            {
                isSatisfied1 = true;
                break;
            }
        }
        // for(i = 0; i < cntResSatisfied -2; i++)/*{{{*/
        // {
        // diffi = HCRes[i+2-1].aaSeqIndex - HCRes[i].aaSeqIndex;
        // if(diffi <= window2)
        // {
        // isSatisfied2 = true;
        // break;
        // }
        // }/*}}}*/
    }
    return (isSatisfied1 );
}/*}}}*/
int  IsSpecChain(const char* id, const char chainid /* = '0'*/)/*{{{*/
    /*****************************************************************************
     * IsSpecChain()
     * Check if the special chainID, e.g. '0', ':' and '_'  indicates a single chain
     * protein or the real 
     * chainID in a multi-chain protein, such as 1GAV
     *
     * return 1 : if chainid really exist in the pdb file
     * return 0 : if chainid does not exist in the pdbfile
     * return -1: if pdb file can not open
     ****************************************************************************/
{
    char path[MAX_PATH+1] = "";
    char pdbid[SIZE_PDBID+1] = "";
    FILE *fp = NULL;

    my_strcpy(pdbid,id,SIZE_PDBID); 
    if (GetPDBFilePath(pdbid,path) != NULL)
        fp = fopen(path,"r");
    else
        return -1;

    int maxline = SIZE_LINE_PDB;
    char line[SIZE_LINE_PDB+1] = "";
    char title[SIZE_RECORD_ID+1] = "";
    bool isSEQRES = false;
    while(fgetline(fp,line,maxline) != EOF)
    {
        sscanf(line,"%6s",title);
        if(strcmp(title,"SEQRES") == 0)
        {
            isSEQRES = true;
            if(line[11] == chainid)
                return 1;
        }
        else if(isSEQRES)
            break;

        if(strcmp(title,"ATOM") == 0 )
            break;
    }
    fclose(fp);
    return 0 ;
}/*}}}*/


/*converter*/
void StdRegExp(char* pattern)/*{{{*/
    /*****************************************************************************
     * regular expression transform,
     * x --> .
     * ( --> {
     * ) --> }
     ****************************************************************************/
{
    int n = strlen(pattern);
    int i;
    for(i = 0 ; i < n ; i++)
    {
        if( pattern[i] == 'x')
            pattern[i] = '.';
        else if( pattern[i] == '(')
            pattern[i] = '{';
        else if( pattern[i] == ')')
            pattern[i] = '}';
    }
}/*}}}*/
char AA3To1(const char* aa3, int mode /*= 0*/)/*{{{*/
/* convert 3char amino acid residue name to 1char res name  */
{
    int i;
    if (mode == 0) /* limited, using the STD3CharAA_alphabet, including 20 amino acids, non standard 3-char residues converted as 'X'*/
    {
        // check for non-residue molecule
        if (BinarySearch_String(aa3,PDB3CharAA_alphabet,strlen(PDB1CharAA_alphabet)) == -1)
            return CHAR_NON_RESIDUE;

        for(i = 0 ; i < NUM_STDAA ; i ++)
        {
            if(strcmp(aa3,STD3CharAA_alphabet[i]) == 0)
                return STD1CharAA_alphabet[i];
        }
        return UNKNOWN_AA; // if not in the alphabet, treate as UNKNOWN
        //LOG: 2006-02-03 14:35:39 Friday  Week 05 <nanjiang@shu>
        // STD1CharAA_alphabet has been changed
        // unknown amico acids treat as 'X' instead of '*'
    }
    else /*full mode, search the 3-char residue name in the PDB3Char_AllRes*/
    {
        int idxres = 0;
        if((idxres = BinarySearch_String(aa3, PDB3Char_AllRes, strlen(PDB1Char_AllRes)))!= -1)
        { return PDB1Char_AllRes[idxres]; }
        else 
        { return UNKNOWN_AA; }
    }
}/*}}}*/
const char* AA1To3(char ch)/*{{{*/
{
    const char *pch;
    if((pch = strchr(STD1CharAA_alphabet,ch ))!= NULL)
    {
        int index = pch - STD1CharAA_alphabet;
        return STD3CharAA_alphabet[index] ;
    }
    else 
        return UNKNOWN_3CharAA;
}
/*}}}*/
void CharToDigit_Protein(const char* str, int8* a,int n)/*{{{*/
{
    int i;
    const char *pch;
    for( i = 0 ; i < n ; i ++ )
    {
        if((pch = strchr(STD1CharAA_alphabet,toupper(str[i]) ))!= NULL)
        {
            a[i] = pch - STD1CharAA_alphabet;
        }
        else 
        {
            a[i] = Char2Digit(UNKNOWN_AA, STD1CharAA_alphabet);
            fprintf(stdout,"Warning, amino acid '%c' no in STD1CharAA_alphabet\n",str[i]);
        }
    }
}/*}}}*/
void DigitToChar_Protein(const int8 *a, int n, char* str)/*{{{*/
{
    int i;
    for( i = 0 ; i < n ; i ++ )
    {
        if(a[i] != DIGIT_INDEL)
            str[i] = STD1CharAA_alphabet[a[i]];
        else
            str[i] = CHAR_INDEL;
    }
    str[i] = '\0';
}/*}}}*/
int DigitShape(char shape)/*{{{*/
{
    const char *pch;
    if((pch = strchr(SHAPE_alphabet,shape)) != NULL)
    {
        int index = pch - SHAPE_alphabet;
        return index;
    }
    else
    {
        fprintf(stderr,"Error! shape '%c' not in SHAPE_alphabet '%s'\n", shape, SHAPE_alphabet);
        assert(strchr(SHAPE_alphabet,shape) != NULL);
    }
    return -1;
}/*}}}*/
int DigitAA(char* resName)/*{{{*/
    // convert 3char residue name to digit-index in STD1CharAA_alphabet
{
    char ch = AA3To1(resName);
    const char *pch;
    if((pch = strchr(STD1CharAA_alphabet,ch)) != NULL)
    {
        int index = pch - STD1CharAA_alphabet;
        return index;
    }
    else
    {
        printf("Error! residue '%s' not in STD3CharAA_alphabet \n", resName);
        assert(strchr(STD1CharAA_alphabet,ch) != NULL);
    }
    return -1;
}/*}}}*/
char DSSP8to3(char dsspSec, int method /*=0*/)/*{{{*/
/*****************************************************************************
 * Map 8-state DSSP secondary structure definition to 3-state secondary
 * structure definition
 * H, G   --> H (helix)
 * E      --> S (strand)
 * others and not '-' --> R (random coil)
 * '-' to 'R' changed 2009-07-11, Nanjiang
 ****************************************************************************/
{
    if (  (dsspSec=='H') || (dsspSec=='G')  ) {
        return 'H';
    } else if (  dsspSec == 'E'  ) {
        return 'S';
    } else if (dsspSec != '-') {
        return 'R';
    } else {
        return dsspSec;
    }
}
/*}}}*/

char* StdID(char* id)/*{{{*/
{
    id  = my_strupr(id);
    if(id[4] == '\0')
    { id[4] = ' '; id[5] = '\0'; } 
    else if(id[4] == '_')
    {
        if(id[5] != '\0') { id[4] = id[5]; id[5] = '\0'; }
        else { id[4] = ' '; id[5] = '\0'; }
    }
    return id;
}
/*}}}*/
char* StdElementName(const char* elementName, char *stdElementName)/*{{{*/
    //elementName will not be changed
{
    my_strcpy(stdElementName,elementName, SIZE_ATOM_ELEMENT);
    if(elementName[0]!= '\0' && islower(elementName[0]))
        stdElementName[0] = toupper(elementName[0]);
    if(elementName[1] != '\0' && isupper(elementName[1]))
        stdElementName[1] = tolower(elementName[1]);

    return stdElementName;
}
/*}}}*/
char* Stdid2Seqresid(char* id)/*{{{*/
    /*****************************************************************************
     * convert standardized chain identifier to the id format in file
     * pdb_seqres.txt, that is 
     *      $PDBID_$chainID
     * in lower case
     ****************************************************************************/
{
    char chainID;
    char pdbid[SIZE_PDBID+1];
    chainID = id[4];
    chainID = toupper(chainID);
    if(chainID == ' ') chainID = '\0';

    my_strcpy(pdbid,id,SIZE_PDBID); 
    my_strlwr(pdbid);
    sprintf(id,"%s_%c",pdbid,chainID);
    return id;
}
/*}}}*/
void RemoveTID(char* id)/*{{{*/
    /*****************************************************************************
     * remove the trailing white space in standardized chain identifier
     ****************************************************************************/
{
    if(id[4] == ' ') { id[4] = '\0'; }
}
/*}}}*/
char *StdPDBFileName2PDBID(const char* rtname_pdbfilename, char *pdbid)/*{{{*/
/*****************************************************************************
 * get the pdbid from the standard pdbfilename, i.e. pdb$pdbid or $pdbid
 * 2008-02-06, Nanjiang
 ****************************************************************************/
{
    if(strncasecmp(rtname_pdbfilename, "pdb", 3) == 0)
    {
        my_strcpy(pdbid, rtname_pdbfilename + 3, SIZE_PDBID);
        return pdbid;
    }
    else if (strlen(rtname_pdbfilename ) == 4)
    {
        my_strcpy(pdbid, rtname_pdbfilename, SIZE_PDBID);
        return pdbid;
    }
    else
    {
        fprintf(stderr,"Error! PDB filename is not standard! rtname = %s\n", rtname_pdbfilename);
        return NULL;
    }
}/*}}}*/


// write formatted text functions
void WriteSubMatrix(FILE* fpSubMatrix, double** subM, int dim,  char* alphabet,int *cnt, char formatValue[] /*= "%4.0lf"*/, char formatFreq[] /*= "%7d"*/)/*{{{*/
    /*********************************************************************
     * WriteSubMatrix()
     * write out the substitution matrix
     * dim: is the dimension of the matrix
     * formatValue is the output format for the matrix value, default = "%4.0f",
     * formatFreq  is the output format for frequency value of each row, default = "%7d"
     * when assigning formatValue, the size should be specified directly after '%'
     *
     * cnt[i] is the effective count of profile positions use to calculate the matrix,
     * set cnt arguments to "NULL" when do not want to print out frequency column
     * 
     * only the trailing arguments can be defaulted
     ********************************************************************/
{
    int i, j;
    char formatAlphabet[100+1]; /* printf format for alphebet like A   B   C   D*/ 
    char formatTextCount[100+1]; /* printf format for training "Count", have the same length as formatFreq */ 
    char formatRowHead[]  = "%2c";


    char *pch;
    char str[10];
    pch = strchr(formatValue,'%');
    if(pch != NULL && IsDigit(*(pch+1)))
    {    
        i = 0;
        while(IsDigit(*++pch)) 
        {
            str[i] = *pch; i++;
        }
        str[i] = '\0';
        sprintf(formatAlphabet,"%c%s%c",'%',str,'c');
    }
    else
        strcpy(formatAlphabet,"%4c");

    pch = strchr(formatFreq,'%');
    if(pch != NULL && IsDigit(*(pch+1)))
    {
        i = 0;
        while(IsDigit(*++pch)) 
        {
            str[i] = *pch; i++;
        }
        str[i] = '\0';
        sprintf(formatTextCount,"%c%s%c",'%',str,'s');
    }
    else
        strcpy(formatTextCount,"%7s");


    fprintf(fpSubMatrix,formatRowHead,' ');
    for(j = 0 ; j < dim ; j ++) 
        fprintf(fpSubMatrix, formatAlphabet,alphabet[j]);
    if(cnt != NULL) fprintf(fpSubMatrix,formatTextCount,"Count");
    fprintf(fpSubMatrix,"\n");

    for ( i = 0 ; i < dim ; i++)
    {
        fprintf(fpSubMatrix,formatRowHead,alphabet[i]);
        for ( j = 0 ; j < dim ; j++)
            fprintf(fpSubMatrix,formatValue, subM[i][j]);
        if(cnt != NULL) fprintf(fpSubMatrix,formatFreq,cnt[i]);
        fprintf(fpSubMatrix,"\n");
    }
    fprintf(fpSubMatrix,"\n");
}/*}}}*/
void WriteResRecord3(char *str, Residue *pRes, int tag)/*{{{*/
{
    if(tag == RESSEQ)
    {
        sprintf(str,"%s%d%c(%c) ",pRes->resName, pRes->resSeq,pRes->resICode,pRes->chainID);
    }
    else if(tag == AASEQINDEX)
    {
        sprintf(str,"%s%d%c(%c) ",pRes->resName, pRes->aaSeqIndex+1,' ',pRes->chainID);
    }
    else
    { }

}/*}}}*/
void WriteCloseMetalRes(FILE* fpout, Residue* res, int numRes)/*{{{*/
    //*******************************************************************
    //WriteCloseMetalRes()
    //*******************************************************************
    // write the residues close to the metal atom in the following format
    // 		HIS37 (A)   HIS51 (A)    "--residue name and resSeq in pdb file
    // 		HIS37 (A)   HIS51 (A)    "--residue name and aaSeqIndex in sequence
    // 		A           V            "--shape string symbol
    // 		100         100          "--conservation level of the residue generated	psi-blast
    //*******************************************************************
{
    int j,k ;
    char str[100+1] = "";
    char separator[] = " , ";

    int fieldwidth_item = 12;
    if(numRes <= 0)
        return; //if there is no residues, return immediately

    //line1: list of "Res resSeq resIcode chainID"
    fprintf(fpout,"\t\t");
    for ( j = 0 ; j < numRes ; j ++)
    {
        WriteResRecord3(str, &res[j], RESSEQ);
        fprintf(fpout,"%-*s ",fieldwidth_item, str);
        if(j < numRes -1) fprintf(fpout,"%s", separator);
    }
    fprintf(fpout,"\n");

    //line2: line of "Res aaSeqIndex ' ' chainID"
    fprintf(fpout,"\t\t");
    for ( j = 0 ; j < numRes ; j ++)
    {
        WriteResRecord3(str, &res[j], AASEQINDEX);
        fprintf(fpout,"%-*s ",fieldwidth_item, str);
        if(j < numRes -1) fprintf(fpout,"%s", separator);
    }
    fprintf(fpout,"\n");

    //line3: line of atomName
    fprintf(fpout,"\t\t");
    for ( j = 0 ; j < numRes ; j ++)
    {
        strcpy(str,"");
        for(k = 0 ; k < res[j].numAtom; k ++)
        {
            strcat(str,res[j].atom[k].name); strcat(str," ");
        }
        fprintf(fpout,"%-*s ",fieldwidth_item, str);
        if(j < numRes -1) fprintf(fpout,"%s", separator);
    }
    fprintf(fpout,"\n");


    //line4: line of distance from atom to the metal atom
    fprintf(fpout,"\t\t");
    char str1[100] = "";
    for ( j = 0 ; j < numRes ; j ++)
    {
        strcpy(str,"");
        for(k = 0 ; k < res[j].numAtom; k ++)
        {
            sprintf(str1,"%.3lf",res[j].atom[k].dist[0]);
            strcat(str,str1); strcat(str," ");
        }
        fprintf(fpout,"%-*s ",fieldwidth_item, str);
        if(j < numRes -1) fprintf(fpout,"%s", separator);
    }
    fprintf(fpout,"\n");

    //line5: list of shape string
    fprintf(fpout,"\t\t"); // head line third
    for ( j = 0 ; j < numRes ; j ++)
    {
        fprintf(fpout,"%-*c ",fieldwidth_item, res[j].shape);
        if(j < numRes -1) fprintf(fpout,"%s", separator);
    }
    fprintf(fpout,"\n");


    //line6: list of conservation value;
    fprintf(fpout,"\t\t"); // head line third
    for ( j = 0 ; j < numRes ; j ++)
    {
        fprintf(fpout,"%-*d ",fieldwidth_item, Integer(res[j].consv));
        if(j < numRes -1) fprintf(fpout,"%s", separator);
    }
    fprintf(fpout,"\n");

}/*}}}*/
void WriteMetalEnvRes(FILE* fpMetalEnvShape, AtomEnv* pAtomEnv, bool isPrintDebugInfo /* = false*/)/*{{{*/
    /*****************************************************************************
     * Write metalEnv for each binding metal ion in the following format
     *
     * 1A5T  |0| 334|4|CYS50 ( )    CYS59 ( )    CYS62 ( )    CYS65 ( )   // resSerin PDB file 
     * 1A5T  |1| 334|4|CYS50 ( )    CYS59 ( )    CYS62 ( )    CYS65 ( )   // aaSeqIndex 
     * 1A5T  |2| 334|4|R            A            R            A           // shape string
     * 1A5T  |3| 334|4|61           59           61           62          // conservation level 
     *
     ****************************************************************************/
{
    int  i;
    char str[100+1] = "";
    char shape;
    DATATYPE_CONSV  consv;

    //line1:"res resSeq resIcode (chainID)"
    fprintf(fpMetalEnvShape,"%s |0|%4d|%d|",pAtomEnv->id,pAtomEnv->seqLength,pAtomEnv->numRes);
    for(i = 0 ; i <pAtomEnv->numRes ; i ++)
    {
        sprintf(str,"%s%d%c(%c)",pAtomEnv->res[i].resName,pAtomEnv->res[i].resSeq,pAtomEnv->res[i].resICode,pAtomEnv->res[i].chainID);
        fprintf(fpMetalEnvShape,"%-12s ",str);
    }
    fprintf(fpMetalEnvShape,"\n");

    //line2:"res aaSeqIndex+1 ' ' (chainID)"
    fprintf(fpMetalEnvShape,"%s |1|%4d|%d|",pAtomEnv->id,pAtomEnv->seqLength,pAtomEnv->numRes);
    for(i = 0 ; i < pAtomEnv->numRes ; i ++)
    {
        sprintf(str,"%s%d%c(%c)",pAtomEnv->res[i].resName,pAtomEnv->res[i].aaSeqIndex+1,' ',pAtomEnv->res[i].chainID);
        fprintf(fpMetalEnvShape,"%-12s ",str);
    }
    fprintf(fpMetalEnvShape,"\n");

    //line3: "shapeString"
    fprintf(fpMetalEnvShape,"%s |2|%4d|%d|",pAtomEnv->id,pAtomEnv->seqLength,pAtomEnv->numRes);
    for(i = 0 ; i < pAtomEnv->numRes ; i ++)
    {
        shape = pAtomEnv->res[i].shape;
        fprintf(fpMetalEnvShape,"%-12c ",pAtomEnv->res[i].shape);
        // shapeFreqs[i][DigitShape(shape)] ++ ;
        // AAFreqs[i][DigitAA(pAtomEnv->res[i].resName)] ++;
    }
    fprintf(fpMetalEnvShape,"\n");	

    //line4: "conservation level"
    fprintf(fpMetalEnvShape,"%s |3|%4d|%d|",pAtomEnv->id,pAtomEnv->seqLength,pAtomEnv->numRes);
    for(i = 0 ; i < pAtomEnv->numRes ; i ++)
    {
        consv = pAtomEnv->res[i].consv;
        fprintf(fpMetalEnvShape,"%-12.2lf ",pAtomEnv->res[i].consv);
        if(isPrintDebugInfo)
        {
            if(pAtomEnv->res[i].consv == INIT_CONSV)
            {
                fprintf(stderr,"conservation value in closeMetalFile\t%s%c\n",pAtomEnv->metalAtomPDBID,pAtomEnv->res[i].chainID);

                if(strlen(pAtomEnv->metalAtomPDBID) < 4)
                {
                    fprintf(stderr, "invalid PDBID: %s\n",pAtomEnv->metalAtomPDBID);
                }
            }
        }
    }
    fprintf(fpMetalEnvShape,"\n");	
}
/*}}}*/
void WriteVectorRecordID(char *vectorRecordID, char *rmtID, int length, int idx, char *res_1char_list,int* aaSeqIndex, int numSite)/*{{{*/
    /*****************************************************************************
     * rmtID          : chain identifier with the last white space removed
     * length         : sequence length
     * idx            : vector vector index for that chain ID, 0, 1, 2 ...
     * res_1char_list : 1 char residue array
     * aaSeqIndex     : aaSeqIndex array
     * numSite        : number of site used for that vector, that is, number of residues used in each vector
     *
     * vectorRecordID is of format
     * ID_length_idx_Res1_aaSeqIndex1_Res2_aaSeqIndex2_...
     ****************************************************************************/
{
    int i;
    char str[20] = "";
    strcpy(vectorRecordID,"");
    sprintf(vectorRecordID,"%s%c%d%c%d",rmtID,CHAR_VECTOR_ID_SEPRATOR, length,CHAR_VECTOR_ID_SEPRATOR, idx);
    for(i = 0 ; i < numSite ; i++)
    {
        sprintf(str,"%c%c%c%d",CHAR_VECTOR_ID_SEPRATOR, res_1char_list[i], CHAR_VECTOR_ID_SEPRATOR, aaSeqIndex[i]+1);
        strcat(vectorRecordID,str);
    }
}
/*}}}*/
void WriteCHDEComment(FILE *fpout)/*{{{*/
{
    fprintf(fpout,"# Tag    -- mark for the residue\n");
    fprintf(fpout,"# Index  -- serial number for the amino acid in sequence\n");
    fprintf(fpout,"# AA     -- amino acid\n");
    fprintf(fpout,"# A R... -- Ala, Arg ...\n");
    fprintf(fpout,"# s1     -- information per position of the profile generated by PSI-BLAST\n");
    fprintf(fpout,"# s2     -- relative weight of gapless real matches to pseudocounts generated by PSI-BLAST\n");
}
/*}}}*/
void WriteCHDE(FILE *fpout, char *id, int *speResIndex, int numSpeRes, int *resIndex, int numRes, MODM *pMODM, bool isWriteHeader /*= true*/, bool isWriteCounter /*= false*/)/*{{{*/
    //void WriteCHDE(FILE *fpout, char *id, MetalPro *pZnPro, Residue *HCRes, int numHCRes, MODM *pMODM)
    // for one chain,
    // given an array of residue aaSeqIndex and another array of aaSeqIndex for
    // special residues, 
{
    int i,j;
    int idx;
    char tag[20] = "";

    int *indexArray = NULL;
    if(numSpeRes > 0)
    {
        indexArray = new int [numSpeRes+1];
        for( i = 0 ; i < numSpeRes ; i++)
            indexArray[i] = speResIndex[i];
        sort(indexArray, indexArray+numSpeRes);
    }

    if(isWriteHeader)
    {
        fprintf(fpout,">%s, sequence length = %4d", id, pMODM->length);
        if(isWriteCounter)
            fprintf(fpout,", numer of residues = %4d, number of '%c' residues = %4d", numRes, '*', numSpeRes);
        fprintf(fpout,"\n");
        // Tag: indicating whether zinc-binding or not, '*' : zinc binding, ' ' : non zinc-binding
        // s1: information of the profile
        // s2: last column in PSI-BLAST profile
        // Tag Index AA A R ... s1 s2     --header line
        fprintf(fpout,"%-5s %5s %2s", "Tag", "Index", "AA");
        for(j = 0 ; j < 20 ; j++)
            fprintf(fpout, " %3c", STD1CharAA_alphabet[j]);
        fprintf(fpout,"  %4s  %4s", "s1","s2");
        fprintf(fpout,"\n");
    }
    for( i = 0 ; i < numRes ; i++)
    {
        idx = resIndex[i];
        if(numSpeRes > 0)
        {
            if(binary_search(indexArray,indexArray+numSpeRes,idx ))
                strcpy(tag,"*");
            else
                strcpy(tag,"");
        }
        else 
            strcpy(tag,"");

        fprintf(fpout,"%-5s %5d %2c", tag, idx+1, pMODM->aaSeq[idx]);
        for(j = 0 ; j < 20 ; j ++)
            fprintf(fpout," %3d" , Integer(pMODM->M[idx][j]));
        fprintf(fpout,"  %4.2f  %4.2f", pMODM->score1[idx], pMODM->score2[idx]);
        fprintf(fpout,"\n");
    }
    if(indexArray != NULL)
        delete [] indexArray;
}
/*}}}*/
void WritePredictStat(FILE* fp,  int ReP, int ReN, int PrP, int PrN, int TP, int FP, int TN, int FN, double Sn, double Sp,int ReP1,int ReN1, int ReP2,int  ReN2,double cutoff_consv,int cutoff_window,char* resList, bool isWriteHeader /*= false*/)/*{{{*/
{
    if(isWriteHeader)
    {
        fprintf(fp,"# zinc binding residue prediction based on HCRes only\n");
        fprintf(fp,"# zinc ions with coordination number >= 3 and <= 4 are considered as biological binding\n");
        fprintf(fp,"# on residue level:\n");
        fprintf(fp,"#       real positives are those zinc binding residues with coordination between [3,4],\n");
        fprintf(fp,"#       within resList, say \"CHDE\", and with score2 >= cutoff_score2 (default 0.1)\n");
        fprintf(fp,"#       - ReP1 are those biologically zinc binding residues within resList, but without\n");
        fprintf(fp,"#         cutoff_score2 restrict, ReN1 also need to be in resList\n");
        fprintf(fp,"#       - ReP2 are all biologically zinc binding residues\n");
        fprintf(fp,"#         ReP2 + ReN2 is equal to the total number of residues in nrPDB\n");
        fprintf(fp,"#\n");
        fprintf(fp,"# on protein level:\n");
        fprintf(fp,"#       real positives are those proteins with at least one residue denoted as real positive\n");
        fprintf(fp,"#       - proteins in ReP1 should have at least one residue in resList biologically binding to zinc \n");
        fprintf(fp,"#         but without cutoff_score2 restrict\n");
        fprintf(fp,"#\n");
        fprintf(fp,"# consv-- conservation level in percentage, derived from PSI-BLAST weighted percentages\n");
        fprintf(fp,"# ReP  -- number of real positive samples\n");
        fprintf(fp,"# ReN  -- number of real negtive samples\n");
        fprintf(fp,"# PrP  -- number of predicted positive samples, \n");
        fprintf(fp,"# PrN  -- number of predicted negtive samples, PrN = (ReP+ReN) - PrP\n");
        fprintf(fp,"# TP   -- number of true positives,  TP = intersect(ReP, PrP)\n");
        fprintf(fp,"# FP   -- number of false positives, FP = PrP - TP\n");
        fprintf(fp,"# TN   -- number of true negtives ,  TN = PrN - FN\n");
        fprintf(fp,"# FN   -- number of false positives ,FN = intersect(ReP, PrN)\n");
        fprintf(fp,"# Sn   -- sensitivity, (or recall)   = TP / (TP + FN)\n");
        fprintf(fp,"# Sp   -- specificity, (or accuracy) = TP / (TP + FP)\n");
        fprintf(fp,"# ReP1 -- real positives,   without cutoff_score2 filted\n");
        fprintf(fp,"# ReN1 -- real negtives,  without cutoff_score2 filtered \n");
        fprintf(fp,"# ReP2 -- real positives for all zinc binding proteins\n");
        fprintf(fp,"# ReN2 -- all non zinc binding proteins, ReP2+ReN2 is equal to the number of nrPDB\n");
        fprintf(fp,"#\n");
        fprintf(fp,"# resList : %s\n",resList);
        fprintf(fp,"%8s", "window");
        fprintf(fp,"%8s", "consv");
        fprintf(fp,"%8s", "ReP");
        fprintf(fp,"%8s", "ReN");
        fprintf(fp,"%8s", "PrP");
        fprintf(fp,"%8s", "PrN");
        fprintf(fp,"%8s", "TP");
        fprintf(fp,"%8s", "FP");
        fprintf(fp,"%8s", "TN");
        fprintf(fp,"%8s", "FN");
        fprintf(fp,"%8s", "Sn");
        fprintf(fp,"%8s", "Sp");
        fprintf(fp,"%8s", "ReP1");
        fprintf(fp,"%8s", "ReN1");
        fprintf(fp,"%8s", "ReP2");
        fprintf(fp,"%8s", "ReN2");
        fprintf(fp,"\n");
    }
    // consv  ReP  ReN  PrP  PrN  TP  FP  TN   FN   Sn  Sp    ReP1   ReN1   ReP2  ReN2
    fprintf(fp,"%8d", cutoff_window);
    fprintf(fp,"%8.2lf", cutoff_consv);
    fprintf(fp,"%8d", ReP);
    fprintf(fp,"%8d", ReN);
    fprintf(fp,"%8d", PrP);
    fprintf(fp,"%8d", PrN);
    fprintf(fp,"%8d", TP);
    fprintf(fp,"%8d", FP);
    fprintf(fp,"%8d", TN);
    fprintf(fp,"%8d", FN);
    fprintf(fp,"%8.3lf", Sn);
    fprintf(fp,"%8.3lf", Sp);
    fprintf(fp,"%8d", ReP1);
    fprintf(fp,"%8d", ReN1);
    fprintf(fp,"%8d", ReP2);
    fprintf(fp,"%8d", ReN2);
    fprintf(fp,"\n");
}/*}}}*/
void WriteMODMTitle(const char* alphabet, FILE* fpout /*= stdout*/ )/*{{{*/
{
    int n = strlen(alphabet);
    int i;
    fprintf(fpout,"%4s %-2s","Num","AA");
    for(i = 0 ; i < n ; i++) fprintf(fpout,"%4c",alphabet[i]);
    fprintf(fpout,"\n");
}/*}}}*/
void WriteMODMProfile(int index, char aa, int* V, double score1, double score2,const char* alphabet, FILE* fpout /*= stdout*/ )/*{{{*/
{
    int n = strlen(alphabet);
    int i;
    fprintf(fpout,"%4d %-2c",index,aa);
    for(i = 0 ; i < n ; i++) fprintf(fpout,"%4d",V[i]);
    fprintf(fpout,"  %3.2lf  %3.2lf",score1,score2);
    fprintf(fpout,"\n");
}/*}}}*/
void WriteFastaSeq(char *str, FILE *fpout /*=stdout*/,int beg/* = 0*/ , int end /*= 0x7FFFFFFF*/ , int linelength  /*= 70*/)/*{{{*/
{
    int size = 0;
    size = strlen(str);
    if(end >= size) 
    {  end = size -1;  }
    int idx = max(beg,0);
    str += idx;
    int i = 0; 
    while (idx <= end)
    {
        putc(*str,fpout); i++; idx++; str++;
        if(i >= linelength)
        {
            putc('\n',fpout); i = 0;
        }
    }
    if(i != 0) putc('\n',fpout);
}
/*}}}*/
void WritePDBAtomRecord(Atom *pAtom, FILE *fpout/* = stdout*/ )/*{{{*/
{
    fprintf(fpout,"%-6s%5d%1s%4s%1c%3s%1s%1c%4d%1c%3s%8.3lf%8.3lf%8.3lf%6.2f%6.2f%6s%-4.4s%-2.2s%-2.2s\n",
            pAtom -> recordID,
            pAtom -> serial,
            "",
            pAtom -> origName,
            pAtom -> altLoc,
            pAtom -> resName,
            "",
            pAtom -> chainID,
            pAtom -> resSeq,
            pAtom -> iCode,
            "",
            pAtom -> x,
            pAtom -> y,
            pAtom -> z,
            pAtom -> occupancy,
            pAtom -> tempFactor,
            "",
            pAtom -> segID,
            pAtom -> element,
            pAtom -> charge);
}
/*}}}*/

template <class T> int WriteBinaryMODM(const char *outfile, const char *alphabet, int length, T * profile, double *parameter, int8 typeProfile)/*{{{*/
{
    FILE *fpout;
    fpout = fopen(outfile, "wb");
    if (checkfilestream(fpout, outfile, "wb") == -1)
    {
        return -1;
    }
    int sizeAlphabet = strlen(alphabet);
    
    fwrite(&typeProfile, sizeof(int8), 1, fpout); /*write the typeProfile*/ 
    fwrite(&sizeAlphabet, sizeof(int), 1, fpout); /*write the sizeAlphabet*/
    fwrite(&length, sizeof(int),1,fpout);/*write the numRes of the profile*/
    fwrite(alphabet, sizeof(char), sizeAlphabet, fpout);  /*write the alphabet of the profile*/ 
    fwrite(parameter, sizeof(double), 8, fpout);  /*write the parameter of the profile*/ 
    fwrite(profile, sizeof(T), length, fpout);   /*write the profile*/ 
    fclose(fpout);
    return 0;
}
template int WriteBinaryMODM <Profile> (const char *outfile, const char *alphabet, int length, Profile * profile, double *parameter, int8 typeProfile);
template int WriteBinaryMODM <ProfileByte> (const char *outfile, const char *alphabet, int length, ProfileByte * profile, double *parameter, int8 typeProfile);
template int WriteBinaryMODM <ProfileSAD> (const char *outfile, const char *alphabet, int length, ProfileSAD * profile, double *parameter, int8 typeProfile);
template int WriteBinaryMODM <ProfileSADByte> (const char *outfile, const char *alphabet, int length, ProfileSADByte * profile, double *parameter, int8 typeProfile);
/*}}}*/

int WriteBinaryFragShort5(const char *outfile, int fragFileType, char **idList, int numID, int maxSizeID, int length, short *posTar, short *numCan, int totalFragCan, FragCanShort5 *fragCan)/*{{{*/
/* Wrtie the frag file in binary format, 2009-06-23 
 * totalFragCan: the total number of candidate fragments
 * All size related parameters are written at the beginning of the file, so that they
 * can be read in the beginning so that memory can be allocated. 2009-06-25
 * */
{
    FILE *fpout;
    fpout = fopen(outfile, "wb");
    if (checkfilestream(fpout, outfile, "wb") == -1)
    {
        return -1;
    }
    int i = 0;
    
    fwrite(&fragFileType, sizeof(int), 1, fpout); /*write the fragFileType*/ 
    fwrite(&numID, sizeof(int), 1, fpout); /*write the number of unique ids*/
    fwrite(&maxSizeID, sizeof(int),1,fpout);/*write the maximal size of ids*/
    fwrite(&length, sizeof(int), 1, fpout);  /*write length of the target sequence*/ 
    fwrite(&totalFragCan, sizeof(int), 1, fpout);   /*write the total number of candidate fragments*/ 

    /*write the idList*/
    for (i = 0; i < numID; i ++)
    {
        fwrite(idList[i], sizeof(char), maxSizeID, fpout );  /*idList is allocated as (numID, maxSizeID+1)*/
    }
    /*write the candidate fragments for each target fragment*/
    int start = 0; /*the start index for fragCan*/
    for (i = 0; i < length; i ++)
    {
        fwrite(&(posTar[i]), sizeof(short), 1, fpout);
        fwrite(&(numCan[i]), sizeof(short), 1, fpout);
        FragCanShort5 *pFragCan = &(fragCan[start]);
        fwrite(pFragCan, sizeof(FragCanShort5), numCan[i], fpout);
        start = start + numCan[i];
    }
    fclose(fpout);
    return 0;
}/*}}}*/
int WriteBinaryFragShort6(const char *outfile, int fragFileType, char **idList, int numID, int maxSizeID, int length, short *posTar, short *numCan, int totalFragCan, FragCanShort6 *fragCan)/*{{{*/
/* Wrtie the frag file in binary format, 2009-06-23 
 * totalFragCan: the total number of candidate fragments
 * All size related parameters are written at the beginning of the file, so that they
 * can be read in the beginning so that memory can be allocated. 2009-06-25
 * */
{
    FILE *fpout;
    fpout = fopen(outfile, "wb");
    if (checkfilestream(fpout, outfile, "wb") == -1)
    {
        return -1;
    }
    int i = 0;
    
    fwrite(&fragFileType, sizeof(int), 1, fpout); /*write the fragFileType*/ 
    fwrite(&numID, sizeof(int), 1, fpout); /*write the number of unique ids*/
    fwrite(&maxSizeID, sizeof(int),1,fpout);/*write the maximal size of ids*/
    fwrite(&length, sizeof(int), 1, fpout);  /*write length of the target sequence*/ 
    fwrite(&totalFragCan, sizeof(int), 1, fpout);   /*write the total number of candidate fragments*/ 

    /*write the idList*/
    for (i = 0; i < numID; i ++)
    {
        fwrite(idList[i], sizeof(char), maxSizeID, fpout );  /*idList is allocated as (numID, maxSizeID+1)*/
    }
    /*write the candidate fragments for each target fragment*/
    int start = 0; /*the start index for fragCan*/
    for (i = 0; i < length; i ++)
    {
        fwrite(&(posTar[i]), sizeof(short), 1, fpout);
        fwrite(&(numCan[i]), sizeof(short), 1, fpout);
        FragCanShort6 *pFragCan = &(fragCan[start]);
        fwrite(pFragCan, sizeof(FragCanShort6), numCan[i], fpout);
        start = start + numCan[i];
    }
    fclose(fpout);
    return 0;
}/*}}}*/
int WriteBinaryFragInt5(const char *outfile, int fragFileType, char **idList, int numID, int maxSizeID, int length, short *posTar, short *numCan, int totalFragCan, FragCanInt5 *fragCan)/*{{{*/
/* Wrtie the frag file in binary format, 2009-06-23 
 * totalFragCan: the total number of candidate fragments
 * All size related parameters are written at the beginning of the file, so that they
 * can be read in the beginning so that memory can be allocated. 2009-06-25
 * */
{
    FILE *fpout;
    fpout = fopen(outfile, "wb");
    if (checkfilestream(fpout, outfile, "wb") == -1)
    {
        return -1;
    }
    int i = 0;
    
    fwrite(&fragFileType, sizeof(int), 1, fpout); /*write the fragFileType*/ 
    fwrite(&numID, sizeof(int), 1, fpout); /*write the number of unique ids*/
    fwrite(&maxSizeID, sizeof(int),1,fpout);/*write the maximal size of ids*/
    fwrite(&length, sizeof(int), 1, fpout);  /*write length of the target sequence*/ 
    fwrite(&totalFragCan, sizeof(int), 1, fpout);   /*write the total number of candidate fragments*/ 

    /*write the idList*/
    for (i = 0; i < numID; i ++)
    {
        fwrite(idList[i], sizeof(char), maxSizeID, fpout );  /*idList is allocated as (numID, maxSizeID+1)*/
    }
    /*write the candidate fragments for each target fragment*/
    int start = 0; /*the start index for fragCan*/
    for (i = 0; i < length; i ++)
    {
        fwrite(&(posTar[i]), sizeof(short), 1, fpout);
        fwrite(&(numCan[i]), sizeof(short), 1, fpout);
        FragCanInt5 *pFragCan = &(fragCan[start]);
        fwrite(pFragCan, sizeof(FragCanInt5), numCan[i], fpout);
        start = start + numCan[i];
    }
    fclose(fpout);
    return 0;
}/*}}}*/
int WriteBinaryFragInt6(const char *outfile, int fragFileType, char **idList, int numID, int maxSizeID, int length, short *posTar, short *numCan, int totalFragCan, FragCanInt6 *fragCan)/*{{{*/
/* Wrtie the frag file in binary format, 2009-06-23 
 * totalFragCan: the total number of candidate fragments
 * All size related parameters are written at the beginning of the file, so that they
 * can be read in the beginning so that memory can be allocated. 2009-06-25
 * */
{
    FILE *fpout;
    fpout = fopen(outfile, "wb");
    if (checkfilestream(fpout, outfile, "wb") == -1)
    {
        return -1;
    }
    int i = 0;
    
    fwrite(&fragFileType, sizeof(int), 1, fpout); /*write the fragFileType*/ 
    fwrite(&numID, sizeof(int), 1, fpout); /*write the number of unique ids*/
    fwrite(&maxSizeID, sizeof(int),1,fpout);/*write the maximal size of ids*/
    fwrite(&length, sizeof(int), 1, fpout);  /*write length of the target sequence*/ 
    fwrite(&totalFragCan, sizeof(int), 1, fpout);   /*write the total number of candidate fragments*/ 

    /*write the idList*/
    for (i = 0; i < numID; i ++)
    {
        fwrite(idList[i], sizeof(char), maxSizeID, fpout );  /*idList is allocated as (numID, maxSizeID+1)*/
    }
    /*write the candidate fragments for each target fragment*/
    int start = 0; /*the start index for fragCan*/
    for (i = 0; i < length; i ++)
    {
        fwrite(&(posTar[i]), sizeof(short), 1, fpout);
        fwrite(&(numCan[i]), sizeof(short), 1, fpout);
        FragCanInt6 *pFragCan = &(fragCan[start]);
        fwrite(pFragCan, sizeof(FragCanInt6), numCan[i], fpout);
        start = start + numCan[i];
    }
    fclose(fpout);
    return 0;
}/*}}}*/
// programs for zn protein detection
double TransMatrixWeight(Residue * pRes1, Residue * pRes2)/*{{{*/
{
    double w                   = 1.0;
    // double wres[4]             = { 1.0, 1.0, 1.0 ,1.0};         //for CHDE
    // char   resList[NUM_BLOSUM] = "HCDE";
    // int    daa1                = Char2Digit(pRes1->aa,resList);
    // int    daa2                = Char2Digit(pRes2->aa,resList);
    // w *= wres[daa1] ;
    // w *= wres[daa2] ;

    //w *= (pRes1 -> consv / 100.0);
    //w *= (pRes2 -> consv / 100.0) ;
    w *= (pRes1 -> consv / 10.0 );
    w *= (pRes2 -> consv / 10.0 ) ;

    return w;  
}/*}}}*/
int MaskPolyHis(const char *aaSeq, bool *isPolyHis, int length /*= 0*/)/*{{{*/
{ /*return number of poly his tag, return 0 if the length is < sizePolyHis */
    int sizePolyHis = 6;
    if(length == 0) length = strlen(aaSeq);
    if(length < sizePolyHis)
        return 0;

    int i,j;
    //int termiLen = 100;
    char pattern[100+1] = "H{6,}";
    sprintf(pattern,"H{%d,}",sizePolyHis);
    Array1D <regmatch_t> allmatch_1darray(length/sizePolyHis);
    regmatch_t *allmatch = allmatch_1darray.array1D;
    int numMatch;
    numMatch = reg_findall( aaSeq, pattern,allmatch, false);

    for (i = 0 ; i < numMatch ; i++)
    {
        for( j = allmatch[i].rm_so; j < allmatch[i].rm_eo ; j ++)
        {
            isPolyHis[j] = true;
        }
    }
    return numMatch;
}
/*}}}*/
void FilterSSBondPred(const char* gistPredictFile_ss, const char* gistPredictFile_zn, const char* gistPredictFile_zn_2, int numSite_ss, int numSite_zn, double cutoff_discriminant_ss)/*{{{*/
{
    FILE *fpin1;
    FILE *fpin2;
    FILE *fpout;
    fpin1 = fopen(gistPredictFile_ss,"r");
    fpin2 = fopen(gistPredictFile_zn,"r");
    fpout = fopen(gistPredictFile_zn_2,"w");
    checkfilestream(fpin1, gistPredictFile_ss, "r", true);
    checkfilestream(fpin2, gistPredictFile_zn, "r", true);
    checkfilestream(fpout, gistPredictFile_zn_2, "w", true);

    int maxline = 400;
    int linesize = 0;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;
    set <string> strs;
    char id[SIZE_CHAIN_ID+1] = "";
    int tint;
    int i;
    int max_ss_res;
    int cnt_ss_res;
    /*determining the number of proteins in gistPredictFile_ss file    */

    f_neglect_comment(fpin1);
    fgetline(fpin1,line,maxline);
    max_ss_res = MIN_INT;
    while((linesize =  fgetline(fpin1, line, maxline))!= EOF)
    {
        ssubstitute(line,  CHAR_VECTOR_ID_SEPRATOR, ' ');
        sscanf(line,"%5s %d %d", id, &tint, &cnt_ss_res);
        StdID(id);
        cnt_ss_res += 1;
        strs.insert(id);
        if(cnt_ss_res > max_ss_res)
            max_ss_res = cnt_ss_res;
    }

    int numSSID = strs.size();
    Array2D <char> ssIDList_2darray( numSSID+1, SIZE_CHAIN_ID+1);
    Array2D <int> ssResIdx_2darray( numSSID+1, max_ss_res+1);
    Array1D <int> numSSRes_1darray( numSSID+1);
    char **ssIDList = ssIDList_2darray.array2D;
    int  **ssResIdx = ssResIdx_2darray.array2D;
    int   *numSSRes = numSSRes_1darray.array1D;

    /*get ss-bonded protein and residues*/

    int length;
    int idx;
    Array1D <int> aaSeqIndex_1darray(numSite_ss+1);
    Array1D <char> res_1char_list_1darray(numSite_ss+2);
    int  *aaSeqIndex = aaSeqIndex_1darray.array1D;
    char *res_1char_list = res_1char_list_1darray.array1D;
    int label;
    double discriminant;

    char idFormer[SIZE_CHAIN_ID+1] = "";
    int cntID  = 0;
    int cntRes = 0;

    rewind(fpin1);
    f_neglect_comment(fpin1);
    fgetline(fpin1,line,maxline);
    while(1)
    {
        if((linesize = fgetline(fpin1,line,maxline)) != EOF)
        {
            ScanfGistPredictCard(line, id, length, idx, res_1char_list,aaSeqIndex, label, discriminant, numSite_ss);
        }
        if(strcmp(id,idFormer) != 0 || linesize == EOF)
        {
            if(strcmp(idFormer,"") != 0)
            {
                numSSRes[cntID] = cntRes;
                my_strcpy(ssIDList[cntID], idFormer, SIZE_CHAIN_ID);
                cntID ++;
            }
            cntRes = 0;
            my_strcpy(idFormer,id, SIZE_CHAIN_ID);
        }
        if(linesize == EOF) break;

        if(discriminant > cutoff_discriminant_ss)
        {
            for(i = 0 ; i < numSite_ss; i++) 
            {
                ssResIdx[cntID][cntRes] = aaSeqIndex[i];
            }
            cntRes += numSite_ss;
        }
    }
    numSSID = cntID;

    /*filter predicted ss-bonded proteins from zinc binding proteins*/
    fprintf(fpout,"# filter predicted ss-bonded residues\n");
    int ssIDIndex;
    bool isSSBonded = false;
    char rmtID[SIZE_CHAIN_ID+1] = "";
    char vectorRecordID[400+1] = "";
    while ((linesize = fgetline(fpin2, line, maxline)) != EOF)
    {
        if(line[0] == '#') 
            fprintf(fpout,"%s\n",line);
        else if (strncmp(line, "item", 4) == 0)
        {
            fprintf(fpout,"%s\n",line);
        }
        else
        {
            ScanfGistPredictCard(line, id, length, idx, res_1char_list,aaSeqIndex, label, discriminant, numSite_zn);
            isSSBonded = false;
            if( (ssIDIndex = BinarySearch_String(id, ssIDList, numSSID)) != -1)
            {
                for(i = 0 ; i < numSite_zn ; i++)
                {
                    if(binary_search(ssResIdx[ssIDIndex], ssResIdx[ssIDIndex]+numSSRes[ssIDIndex],aaSeqIndex[i]))
                    {
                        isSSBonded = true;
                        break;
                    }
                    if(numSSRes[ssIDIndex] >= 4 && res_1char_list[i] == 'C')/*{{{*/
                        // if there are a lot of cysteins are predicted as ss
                        // bonded, it is likely that all highly conserved cys
                        // are ss-bonded
                    {
                        isSSBonded = true;
                        break;
                    }/*}}}*/
                }
            }
            if(isSSBonded)
            {
                my_strcpy(rmtID,id,SIZE_CHAIN_ID);
                RemoveTID(rmtID);
                sscanf(line,"%s", vectorRecordID);
                fprintf(fpout,"%s   %d  %10.6lf\n",vectorRecordID, -1, -1.0);
            }
            else
                fprintf(fpout,"%s\n",line);
        }
    }
    fclose(fpin1);
    fclose(fpin2);
    fclose(fpout);
}
/*}}}*/
void AtomFrequencyAna_AtomEnv(AtomEnv *atomEnvs, int numAtomEnv,Element *metalEle, int numMetalEle, char *str) /*{{{*/ 
{	
    int i ,j;
    int freq[100];
    for(i=0 ;i < numMetalEle ; i ++) freq[i] = 0;

    for(i = 0 ; i<numAtomEnv ;i ++)
    {
        for(j= 0 ; j<numMetalEle;j++)
        {
            if(strcasecmp(atomEnvs[i].metalAtomName,metalEle[j].name) == 0)
            {
                freq[j] ++ ;
                break;
            }
        }	
    }
    char tmpstr[20] = "";
    char stdElementName[SIZE_ATOM_ELEMENT+1] = "";
    for(i = 0 ; i < numMetalEle ; i++)
    {
        if(freq[i] != 0)
        {
            strcpy(tmpstr,"");
            sprintf(tmpstr,"\t%s\t%d",StdElementName(metalEle[i].name,stdElementName), freq[i]);
            strcat(str,tmpstr);
        }
    }
}/*}}}*/
double TranScore(double ***M,Residue* res,char *alphabet, int numRes)/*{{{*/
{
    double score = 0.0;
    int daa[2];
    int i ;
    int sizeAlphabet = strlen(alphabet);
    for(i = 0 ; i < numRes-1 ;i++)
    {
        daa[0] = Char2Digit(res[i].aa,alphabet,sizeAlphabet);
        daa[1] = Char2Digit(res[i+1].aa,alphabet,sizeAlphabet);
        // debug
        if (daa[0] < 0 || daa[1] < 0)
        {
            printf("alphabet=%s, aa=%c\n",alphabet, res[i].aa);
        }
        // debug end
        assert( daa[0] >= 0 );
        assert( daa[1] >= 0 );
        score += M[i][daa[0]][daa[1]] * TransMatrixWeight(&res[i], &res[i+1]);
    }
    return score;
}/*}}}*/
void MergeGistResidue(GistPredChain *pChain, int *idx, int start, int end,  int &label, double &discriminant, int operation)/*{{{*/
{
    //double sum = 0.0;
    //double max = MIN_DOUBLE;
    int i;
    int n = end - start;
    Array1D <double> a_1darray(n);
    double *a = a_1darray.array1D;
    for(i = 0 ; i < n ; i ++) a[i] = pChain->discriminant[idx[i+start]]; 
    if(operation == USE_MAXIMUM) 
    {
        discriminant = *max_element(a, a + n);
    }
    else
    {
        discriminant = Average(a, 0, n - 1);
    }
    if(discriminant >= 0) label = 1;
    else label = -1;
}
/*}}}*/
void SortUniqueResidue(GistPredChain *pChain, int operation)/*{{{*/
    /*****************************************************************************
     * sort the predicted zinc-binding residues according to the sequence series
     * number, merge the redundant residues according to the specified method
     ****************************************************************************/
{
    int i,j;
    //int numRes;
    GistPredChain tmpChain;
    InitGistPredChain(&tmpChain);
    AllocGistPredChain(&tmpChain, pChain->numRes);
    CopyGistPredChain(&tmpChain, pChain);
    //--- sort the pChain according to the index
    Array1D <int> idx_1darray(tmpChain.numRes);
    int *idx = idx_1darray.array1D;
    for(i = 0 ; i < tmpChain.numRes; i++) idx[i] = i;
    QuickSort_index(idx, tmpChain.aaSeqIndex, 0, tmpChain.numRes-1);

    int label;
    double discriminant;

    //--- Unique the redundant residues by method defined with operation
    i = 0 ;
    int cntRes_unique = 0;
    while( i < tmpChain.numRes)
    {
        j = 0 ;  // cnt for the same residues
        while(tmpChain.aaSeqIndex[idx[i+j]] == tmpChain.aaSeqIndex[idx[i]] )
        {
            j++ ;
            if(i+j >= tmpChain.numRes) break;
        }

        MergeGistResidue(&tmpChain, idx, i, i+j, label, discriminant, operation);
        pChain->aaSeqIndex[cntRes_unique]   = tmpChain.aaSeqIndex[idx[i]];
        pChain->aaSeq[cntRes_unique]        = tmpChain.aaSeq[idx[i]];
        pChain->label[cntRes_unique]        = label;
        pChain->discriminant[cntRes_unique] = discriminant;
        i += j;
        cntRes_unique ++;
    }
    pChain->numRes = cntRes_unique;
    pChain->aaSeq[cntRes_unique] = '\0';

    //free memory
    DeleteGistPredChain(&tmpChain);
}
/*}}}*/
int GetGistPredResidue(const char* gistPredFile, int numSite, GistPredChain *chains, int operation)/*{{{*/
/*****************************************************************************
 * GetGistPredResidue for all chains, enough chains should be allocated for the
 * pointer *chains
 * Bug fixed 2010-05-19, AllocGistPredChain(&tmpChain, LONGEST_SEQ);, when
 * numSite > 1, the total records for one chain can be more than the number of
 * amino acids in the chain. make it the length of the file, which will be
 * safely longer than the number of amino acids.
 ****************************************************************************/
{
    int i = 0;
    int linesize = 0;
    int maxline = 500;

    int maxRecord = 0; /*maxRecord is obtained as the number of lines of the file, which is safely more than cntRes, 2010-05-19*/
    maxRecord = fgetlinecnt(gistPredFile, maxline) * numSite;

    Array1D <char> line_1darray(maxline+10);
    char *line = line_1darray.array1D; 

    char id[SIZE_CHAIN_ID+1] = "";
    char idFormer[SIZE_CHAIN_ID+1] = "";
    int  length = 0;
    Array1D <int>  aaSeqIndex_1darray(numSite+1);
    Array1D <char> res_1char_list_1darray(numSite+2);
    int  *aaSeqIndex = aaSeqIndex_1darray.array1D;
    char *res_1char_list = res_1char_list_1darray.array1D;
    int idx = 0;
    int    label = 0;
    double discriminant = 0.0;

    GistPredChain tmpChain;
    InitGistPredChain(&tmpChain);
    AllocGistPredChain(&tmpChain, maxRecord);


    FILE *fpGistPred = fopen(gistPredFile, "r");
    checkfilestream(fpGistPred, gistPredFile, "r");

    int cntPro = 0 ;
    int cntRes = 0 ;

    while(1)/*{{{*/
    {
        if( (linesize = fgetline(fpGistPred,line,maxline))!= EOF )
        {
            if(linesize > 0 && line[0] != '#' && strncmp(line, "item",4) != 0)
            {
                ScanfGistPredictCard(line, id, length, idx, res_1char_list,aaSeqIndex, label, discriminant, numSite);
                /*cnt1 ++;*/
            }
            else 
                continue;
        }

        if((strcmp (id, idFormer) != 0  || linesize == EOF)) // end of one protein record
        {
            if(strcmp(idFormer,"") != 0 ) // if not the first record
            {
                tmpChain.numRes = cntRes;
                tmpChain.aaSeq[cntRes] = '\0';
                my_strcpy(tmpChain.id , idFormer, SIZE_CHAIN_ID);
                SortUniqueResidue(&tmpChain, operation);
                AllocGistPredChain(&chains[cntPro], tmpChain.numRes);
                CopyGistPredChain(&chains[cntPro],&tmpChain);
                cntPro ++;
            }
            cntRes = 0 ;
            my_strcpy(idFormer,id, SIZE_CHAIN_ID);
        }

        if(linesize == EOF) break;

        for(i = 0 ; i < numSite ; i++)
        {
            tmpChain.length = length ;
            tmpChain.aaSeqIndex[cntRes]   = aaSeqIndex[i];
            tmpChain.aaSeq[cntRes]        = res_1char_list[i];
            tmpChain.label[cntRes]        = label;
            tmpChain.discriminant[cntRes] = discriminant;
            cntRes ++;
            if (cntRes > maxRecord)
            {
                fprintf(stderr,"cntRes (%d) > maxRecord (%d), reading error\n", cntRes, maxRecord);
                assert(cntRes > maxRecord);
            }
        }
    }/*}}}*/
    fclose(fpGistPred);

    DeleteGistPredChain(&tmpChain);
    return cntRes;
}
/*}}}*/
int GetGistPredResidue(FILE* fpGistPred, int numSite, GistPredChain *pChain, int operation)/*{{{*/
    /*****************************************************************************
     * Read in one chain each time
     ****************************************************************************/
{
    int i;
    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D; 

    char id[SIZE_CHAIN_ID+1] = "";
    char idFormer[SIZE_CHAIN_ID+1] = "";
    int  length;
    Array1D <int>  aaSeqIndex_1darray(numSite+1);
    Array1D <char> res_1char_list_1darray(numSite+2);
    int  *aaSeqIndex = aaSeqIndex_1darray.array1D;
    char *res_1char_list = res_1char_list_1darray.array1D;
    int idx;
    int    label;
    double discriminant;


    int cntRes = 0 ;
    fpos_t pos ;
    fpos_t posTwo[2];

    int ww=0;
    assert(fgetpos(fpGistPred, &posTwo[0]) == 0);
    assert(fgetpos(fpGistPred, &posTwo[1]) == 0);

    while(1)/*{{{*/
    {
        assert(fgetpos(fpGistPred, &pos) == 0);
        if( (linesize = fgetline(fpGistPred,line,maxline))!= EOF )
        {
            if(linesize > 0 && line[0] != '#' && strncmp(line, "item",4) != 0)
            {
                ScanfGistPredictCard(line, id, length, idx, res_1char_list,aaSeqIndex, label, discriminant, numSite);
                assert(fgetpos(fpGistPred, &(posTwo[ww])) == 0);
                ww=(ww+1)%2;
                /*cnt1 ++;*/
            } else {
                continue;
            }
        }

        if((strcmp (id, idFormer) != 0  || linesize == EOF)) // end of one protein record
        {
            if(strcmp(idFormer,"") != 0 ) // if not the first record
            {
                pChain->numRes = cntRes;
                pChain->aaSeq[cntRes] = '\0';
                my_strcpy(pChain->id , idFormer, SIZE_CHAIN_ID);
                SortUniqueResidue(pChain, operation);
                if(strcmp(id, idFormer) != 0)
                {
                    assert(fsetpos(fpGistPred, &posTwo[ww]) == 0);
                }
                break;
            }
            //cntRes = 0 ;
            my_strcpy(idFormer,id, SIZE_CHAIN_ID);
        }
        if( linesize == EOF) break;

        for(i = 0 ; i < numSite ; i++)
        {
            pChain->length = length ;
            pChain->aaSeqIndex[cntRes]   = aaSeqIndex[i];
            pChain->aaSeq[cntRes]        = res_1char_list[i];
            pChain->label[cntRes]        = label;
            pChain->discriminant[cntRes] = discriminant;
            cntRes ++;
        }
    }/*}}}*/
    if(linesize == EOF && cntRes == 0) return EOF;
    else return cntRes;
}
/*}}}*/

int UpdateHCRes(Residue *HCRes, int &numRes, int *idx, int numSHCRRes)/*{{{*/
{
    int i;
    Array1D <Residue> tmpres_1darray(numRes);
    Residue *tmpres = tmpres_1darray.array1D;
    for(i = 0; i < numRes; i ++) { InitResidue(&tmpres[i]); } //uninitialized data detected, fixed on 2007-06-20
    for(i = 0 ; i < numRes ; i ++)
    {
        CopyResidue(&tmpres[i], &HCRes[i]);
    }
    for(i = 0 ; i < numSHCRRes ; i ++)
    {
        CopyResidue(&HCRes[i], &tmpres[idx[i]]);
    }
    numRes = numSHCRRes;
    return numRes;
}
/*}}}*/
int HCRes_window(Residue* HCRes, int& numRes, int sizeGroup, int cutoff_window)/*{{{*/
    // filter the HCRes by rule the restrict rule: all sizegroups residues
    // should be within cutoff_window. 
    // for example, if sizeGroups = 2, it means all selected residues should be
    // in a pair whose distance is within cutoff_window.
    //
    //
    // HCRes and numRes will be changed
    // return updated numRes;
    // sizeGroup : number of residues for each group
{
    if (numRes < sizeGroup)
        return -1;

    int i;

    Array1D <int> a_1darray(numRes);
    Array1D <int> idx_1darray(numRes);
    Array1D <Residue> tmpres_1darray(numRes);

    int *a = a_1darray.array1D;
    int *idx = idx_1darray.array1D;
    Residue *tmpres = tmpres_1darray.array1D;
    for(i = 0 ; i < numRes ; i ++)
    {
        CopyResidue(&tmpres[i], &HCRes[i]);
    }

    int nSet = numRes;
    int nSub = sizeGroup;
    int win;
    set <int> shcr_res_idx_set;
    bool done = true;
    while(1)
    {
        ksub_next4(nSet,nSub,a,&done);
        if(done) break;
        for(i = 0; i < nSub; i++) idx[i] = a[i]-1; /* because aaSeqIndex starts from 0, while a[i] starts from 1 
                                                      not that the value of a[i] should not be changed, since it will affect the next subset */
        win = tmpres[idx[nSub-1]].aaSeqIndex - tmpres[idx[0]].aaSeqIndex; 
        if(win <= cutoff_window)
        {
            for(i = 0 ; i < nSub ; i++)
            {
                shcr_res_idx_set.insert(idx[i]);
            }
        }
    }
    numRes = shcr_res_idx_set.size();
    Array1D <int> idx_res_1darray(numRes);
    int *idx_res = idx_res_1darray.array1D;

    Set2Array(shcr_res_idx_set.begin(), shcr_res_idx_set.end(), idx_res);
    for(i = 0 ; i < numRes ; i++)
    {
        CopyResidue(&tmpres[idx_res[i]], &HCRes[i]);
    }
    return numRes;
}
/*}}}*/
int SHCR_window(int* HCResAASeqIndex, int numRes, int numSite, int cutoff_window, int cutoff_window_pair, int* shcr_res_idx, int &numSHCRRes)/*{{{*/
    // filter the HCRes by rule the restrict rule: all sizegroups residues
    // should be within cutoff_window. 
    // for example, if sizeGroups = 2, it means all selected residues should be
    // in a pair whose distance is within cutoff_window.
    //
    //
    // HCRes and numRes will be changed
    // return updated numRes;
    // sizeGroup : number of residues for each group
{
    if (numRes < numSite)
    {
        numSHCRRes = 0;
        return 0;
    }

    int i;

    Array1D <int> a_1darray(numRes);
    Array1D <int> idx_1darray(numRes);
    int *a = a_1darray.array1D;
    int *idx = idx_1darray.array1D;


    int nSet = numRes;
    int nSub = numSite;
    int win;
    int min_win_pair;
    Array1D <int> win_pair_1darray(numSite);
    int *win_pair = win_pair_1darray.array1D;
    set <int> shcr_res_idx_set;
    bool done = true;
    while(1)
    {
        ksub_next4(nSet,nSub,a,&done);
        if(done) break;
        for(i = 0; i < nSub; i++) idx[i] = a[i]-1; /* because aaSeqIndex starts from 0, while a[i] starts from 1 
                                                      not that the value of a[i] should not be changed, since it will affect the next subset */
        win = HCResAASeqIndex[idx[nSub-1]]- HCResAASeqIndex[idx[0]]; 
        if(win <= cutoff_window)
        {
            for(i = 0 ; i < nSub - 1; i ++) win_pair[i] = HCResAASeqIndex[idx[i+1]] - HCResAASeqIndex[idx[i]];
            min_win_pair =  *min_element(win_pair, win_pair+nSub-1) ;
            if(min_win_pair <= cutoff_window_pair)
            {
                for(i = 0 ; i < nSub ; i++) shcr_res_idx_set.insert(idx[i]);
            }
        }
    }
    numSHCRRes = shcr_res_idx_set.size();
    Set2Array(shcr_res_idx_set.begin(), shcr_res_idx_set.end(), shcr_res_idx);
    return numSHCRRes;
}
/*}}}*/
bool IsSatisfySHCRRule(int *HCResAASeqIndex, int numRes, int min_numHCRes, int cutoff_window_pair, int cutoff_window_zn3, int cutoff_window_zn4, int* shcr_res_idx, int &numSHCRRes)/*{{{*/
{
    int i;
    int numSite;


    if(numRes < min_numHCRes )
    {
        numSHCRRes = 0;
        return false;

    }

    if(numRes <= 1) //bug fixed 2007-04-13, for conditions when min_numHCRes is set to < 2
    {
        for(i = 0 ; i < numRes ; i++) shcr_res_idx[i] = i;
        numSHCRRes = numRes;
        return true;
    }
    else if(numRes == 2)
    {
        if(HCResAASeqIndex[1] - HCResAASeqIndex[0] <= cutoff_window_pair)
        {
            for(i = 0 ; i < numRes ; i++) shcr_res_idx[i] = i;
            numSHCRRes = numRes;
            return true;
        }
        else
        {
            numSHCRRes = 0;
            return false;
        }
    }
    Array1D <int> shcr_res_idx_zn3_1darray(numRes);
    Array1D <int> shcr_res_idx_zn4_1darray(numRes);
    int *shcr_res_idx_zn3 = shcr_res_idx_zn3_1darray.array1D;
    int *shcr_res_idx_zn4 = shcr_res_idx_zn4_1darray.array1D;
    int numSHCRRes_zn3;
    int numSHCRRes_zn4;


    numSite = 3;
    SHCR_window(HCResAASeqIndex, numRes, numSite, cutoff_window_zn3, cutoff_window_pair, shcr_res_idx_zn3, numSHCRRes_zn3);

    numSite = 4;
    SHCR_window(HCResAASeqIndex, numRes, numSite, cutoff_window_zn4, cutoff_window_pair, shcr_res_idx_zn4, numSHCRRes_zn4);

    set <int> shcr_res_idx_set;
    for(i = 0 ; i < numSHCRRes_zn3; i++) shcr_res_idx_set.insert(shcr_res_idx_zn3[i]);
    for(i = 0 ; i < numSHCRRes_zn4; i++) shcr_res_idx_set.insert(shcr_res_idx_zn4[i]);

    Set2Array(shcr_res_idx_set.begin(), shcr_res_idx_set.end(), shcr_res_idx);
    numSHCRRes = shcr_res_idx_set.size();

    return numSHCRRes;

}
/*}}}*/

// alignment related functions
template <class T> double GetScore_PPAlignment(T *va, T* vb, T *log_va, T* log_vb, int formula /*= 0 */, const double *P /*=NULL*/, const double *log_P /*= NULL */)/*{{{*/
    /*****************************************************************************
     * Get score for profile profile alignment,
     * ChangeLog: 2007-05-29
     * 1. log calculation is about 30 times time consuming than multiply calculation,
     * avoid it in the innerest of the calculation. 
     * precalculate log_va and log_vb.
     * 2. it is not very appropriate to set va[i]*( log(vb[i])/P[i] ) as 0 when
     * vb[0] is 0, it's better to set vb[i] to a small number, say, 1e-6 when vb[i]
     * is 0
     *
     *
     * formula: different formula for doing profile profile alignment
     * formula 0: 
     * 
     *              ----  
     *         s=   \                   b(i)                  a(i)   
     *              /      a(i) * ln (  ---- ) + b(i) * ln (  ---- ) 
     *              ----                P(i)                  P(i)   
     *             i=1..20 
     * formula 1:
     * 
     *              ----     
     *         s=   \     
     *              /     a(i) * b(i)   
     *              ----
     *             i=1..20 
     *
     ****************************************************************************/
{
    int i;
    int N = 20;
    double score = 0.0;
    if( P == NULL)     { P = Background_AA_Freq; }
    if( log_P == NULL) { log_P = Ln_Background_AA_Freq; }
    switch( formula )
    {
        case 0 : 
            double ta,tb;
            for(i = 0 ; i < N; i ++)
            {
                ta = double( va[i]*(log_vb[i]-log_P[i]));
                tb = double( vb[i]*(log_va[i]-log_P[i]));
                score += (ta+tb);
            }
            break;
        case 1 : 
            for(i = 0 ; i < N; i ++)
            {
                score += double(va[i]*vb[i]);
            }
            break;
        default  : break;
    }  
    return (score);
}
template double GetScore_PPAlignment <int>    (int *v1   , int *v2   , int *log_v1   , int *log_v2   , int formula /*=0*/, const double *P /*=NULL*/, const double *log_P /*=NULL*/);
template double GetScore_PPAlignment <float>  (float *v1 , float *v2 , float *log_v1 , float *log_v2 , int formula /*=0*/, const double *P /*=NULL*/, const double *log_P /*=NULL*/);
template double GetScore_PPAlignment <double> (double *v1, double *v2, double *log_v1, double *log_v2, int formula /*=0*/, const double *P /*=NULL*/, const double *log_P /*=NULL*/);
/*}}}*/

template <class T> int Alignment(char *Xstr, char *Ystr, char *alphabet, int m, int n,char *title1,char* title2,  /*{{{*/ 
        char* alignXstr0, char* alignYstr0, int* alignRel0,	
        AlignFactor *pAlignFactor, T gapOpen, T gapExt,T **SM, int seqtype /*= AA_SEQ */)
/*****************************************************************************
 * gapOpen and gapExt use float number, and thus matrix values are 
 * calculated in float
 ****************************************************************************/
{
    /* first we reserve memory for all the variables we will use */
    int i,j;
    int Imax,Jmax;
    int iAlign;
    float maxValue, tmp;
    int M = m+1;
    int N = n+1;
    int NALIGN = n+m+1; /* maximum size after alignment*/ 
    int sizealphabet = strlen(alphabet);

    Array2D <float> V_2darray(M,N);
    Array2D <int> trace_2darray(M,N);
    float **V   = V_2darray.array2D;     /* value label matrix */
    int **trace = trace_2darray.array2D; /* trace matrix       */

    Array1D <int8> X_1darray(M);
    Array1D <int8> Y_1darray(N);
    Array1D <int8> alignX_1darray(NALIGN);
    Array1D <int8> alignY_1darray(NALIGN);
    Array1D <char> alignXstr_1darray(NALIGN);
    Array1D <char> alignYstr_1darray(NALIGN);
    Array1D <int>  alignRel_1darray(NALIGN);

    int8 *X         = X_1darray.array1D;         /* digital X sequence         */
    int8 *Y         = Y_1darray.array1D;         /* digital Y sequence         */
    int8 *alignX    = alignX_1darray.array1D;    /* digital aligned X sequence */
    int8 *alignY    = alignY_1darray.array1D;    /* digital aligned Y sequence */
    char *alignXstr = alignXstr_1darray.array1D; /* aligned X sequence         */
    char *alignYstr = alignYstr_1darray.array1D; /* aligned Y sequence         */
    int  *alignRel  = alignRel_1darray.array1D;  /* alignment relation         */


    for(i = 0 ; i < m ; i ++) X[i] = Charcase2Digit(Xstr[i], alphabet, sizealphabet);
    for(i = 0 ; i < n ; i ++) Y[i] = Charcase2Digit(Ystr[i], alphabet, sizealphabet);

    // initializing
    for(i=0;i<=m;i++) V[i][0] = i*gapExt; /* do this for i=0,1,...,m */
    for(j=0;j<=n;j++) V[0][j] = j*gapExt;// for the first row and first column, set 0

    for(i=0;i<=m;i++) trace[i][0]=STOP;
    for(j=0;j<=n;j++) trace[0][j]=STOP;

    /* labeling of all nodes, this is the main loop of the algorithm  */
    float e,f,g; /* e --> up, f --> left, g --> diagnol*/ 
    float s ;
    for(i = 1 ; i <= m ; i ++)                 // matrix ------- seqY
    {    /* note: we begin at i=1 ! */         //        |
        for(j = 1 ; j <= n ; j ++)             //        |
        {                                      //        seqX
            s = T (SM[X[i-1]][Y[j-1]]);

            if(trace[i-1][j] != UP)
                e = V[i-1][j] + gapOpen + gapExt ;
            else
                e = V[i-1][j] + gapExt ;

            if(trace[i][j-1] != LEFT)
                f = V[i][j-1] + gapOpen + gapExt ;
            else
                f = V[i][j-1]  + gapExt ;			

            g = V[i-1][j-1] + s ;

            tmp = f ;
            trace[i][j] = LEFT ;

            if(e> tmp)
            {
                tmp = e ;
                trace[i][j] = UP ;
            }
            if(g > tmp)
            {
                tmp = g ;
                trace[i][j] = DIAG ;
            }
            V[i][j] = tmp ;
        }
    }

    // find the highest score in bottom row and rightmost column of the value
    // matrix
    maxValue = MIN_FLOAT ;
    Imax     = 0;
    Jmax     = 0;
    for ( i = 1; i <= m ; i ++ )  /* scan the rightmost column*/ 
    {
        if( V[i][n] > maxValue)
        {
            Imax = i ;
            Jmax = n ;
            maxValue = V[i][n]; /*2009-11-12, bug fixed, put the maxValue inside the comparison*/
        }
            //maxValue = V[i][n]; 
    }
    for ( j = 1; j <= n ; j ++ )  /* scan the bottom row*/ 
    {
        if( V[m][j] > maxValue)
        {
            Imax = m ;
            Jmax = j ;
            maxValue = V[m][j] ;
        }
    }

    //[>{{{<]
    //    FILE *fplog = fopen("log.dat","w");
    //    if(1)
    //    {
    //      fprintf(fplog, "alignment value matrix\n");;
    //      fprintf(fplog,"%4c %4c", ' ', ' '); for(j = 0 ; j < n ; j ++) fprintf(fplog, "%4c ", Ystr[j]);fprintf(fplog,"\n");
    //      for( i = 0 ; i <= m ; i ++)
    //      {   
    //          if(i == 0) fprintf(fplog,"%4c ", ' ');
    //          else       fprintf(fplog,"%4c ",Xstr[i-1]);
    //          for(j= 0 ; j <= n ; j ++)
    //          {
    //              fprintf(fplog,"%4.1f ",V[i][j]);
    //          }
    //          fprintf(fplog,"\n");		
    //      }
    //      fprintf(fplog,"\n\n");
    //      fprintf(fplog,"trace matrix \n");
    //      fprintf(fplog,"%4c %4c", ' ', ' '); for(j = 0 ; j < n ; j ++) fprintf(fplog, "%4c ", Ystr[j]);fprintf(fplog,"\n");
    //      for( i = 0 ; i <= m ; i ++)
    //      {
    //          if(i == 0) fprintf(fplog,"%4c ", ' ');
    //          else       fprintf(fplog,"%4c ",Xstr[i-1]);
    //          for(j= 0 ; j <= n ; j ++)
    //          {
    //              fprintf(fplog,"%4d ",trace[i][j]);
    //          }
    //          fprintf(fplog,"\n");
    //      }
    //      fprintf(fplog,"\n\n");
    //    }
    //    fprintf(fplog,"Imax =%d\njmax = %d\nmaxValue = %f\n",Imax,Jmax,maxValue);
    //    fclose(fplog);
    //[>}}}<]

    /* creating aligned sequences alignY and alignX */
    /* note: they are created in reverse order! */
    i = m;
    j = n;
    iAlign = 0; // the length of alignment

    //unaligned ends 

    while(i>Imax) 
    {      
        alignY[iAlign] = GAP_DIGIT; // -1 --> '-'      
        alignX[iAlign] = X[i-1];      
        i --;
        iAlign ++;
    }

    while(j > Jmax) 
    {
        alignY[iAlign] = Y[j-1];
        alignX[iAlign] = GAP_DIGIT;
        j --;
        iAlign ++;
    }

    /* when we come here we know that i==Imax and j==Jmax */
    /* it is from this position we make the jump to the virtual stop node */

    while(trace[i][j] != STOP) 
    {
        if(trace[i][j] == DIAG) 
        {
            alignY[iAlign]=Y[j-1];
            alignX[iAlign]=X[i-1];
            i--;
            j--;
            iAlign++;
        }
        else if(trace[i][j] == LEFT) 
        {
            alignY[iAlign]=Y[j-1];
            alignX[iAlign]= GAP_DIGIT;
            j--;
            iAlign++;
        }
        else if(trace[i][j] == UP) 
        {
            alignY[iAlign]= GAP_DIGIT;
            alignX[iAlign]=X[i-1];
            i--;
            iAlign++;
        }
    }
    /* unaligned beginning */
    while(i>0) 
    {  
        alignY[iAlign]= GAP_DIGIT;
        alignX[iAlign]=X[i-1];
        i--;
        iAlign++;
    }

    while(j>0) 
    {
        alignY[iAlign]=Y[j-1];
        alignX[iAlign]= GAP_DIGIT;
        j--;
        iAlign++;
    }

    AlignAna(alignX,alignY,iAlign,SM,pAlignFactor, alignRel);
    pAlignFactor->identity_short = pAlignFactor->identity * iAlign / min(m,n) ;
    pAlignFactor->similarity_short = pAlignFactor->similarity * iAlign / min(m,n) ;
	pAlignFactor->score = float(V[Imax][Jmax]);
//	alignFactor->pozScore = float(POZ_Score(SM,alignFactor->score,pairCnt));
    /*not finished, 2007-12-21*/

    // translate digital aligned sequences back into ascii characters
    for(i = 0 ; i < iAlign ; i ++) 
    {    
        if(alignX[i] == GAP_DIGIT) 
        {
            alignXstr[i] = GAP_CHAR;
        }
        else alignXstr[i] = Digit2Char(alignX[i], alphabet, sizealphabet);
    }
    for(i = 0 ; i < iAlign ; i ++) 
    {    
        if(alignY[i] == GAP_DIGIT) 
        {
            alignYstr[i] = GAP_CHAR;
        }
        else alignYstr[i] = Digit2Char(alignY[i], alphabet, sizealphabet);
    }

    int cnt = 0 ;
    for(i = iAlign-1 ; i >= 0 ; i -- )
    {
        alignXstr0[cnt] = alignXstr[i];
        alignYstr0[cnt] = alignYstr[i];
        alignRel0[cnt] = alignRel[i];
        cnt ++;
    }
    alignXstr0[cnt] = '\0';
    alignYstr0[cnt] = '\0';

    return iAlign;
}
template int Alignment<int> (char *Xstr, char *Ystr, char *alphabet, int m, int n,char *title1,char* title2, char* alignXstr0, char* alignYstr0, int* alignRel0, AlignFactor *pAlignFactor, int gapOpen, int gapExt,int **SM, int seqtype /*= AA_SEQ */);
template int Alignment<float> (char *Xstr, char *Ystr, char *alphabet, int m, int n,char *title1,char* title2, char* alignXstr0, char* alignYstr0, int* alignRel0, AlignFactor *pAlignFactor, float gapOpen, float gapExt, float **SM, int seqtype /*= AA_SEQ */);
template int Alignment<double> (char *Xstr, char *Ystr, char *alphabet, int m, int n,char *title1,char* title2, char* alignXstr0, char* alignYstr0, int* alignRel0, AlignFactor *pAlignFactor, double gapOpen, double gapExt, double **SM, int seqtype /*= AA_SEQ */);
/*}}}*/
template <class T> int Alignment_Profile(char *Xstr, char *Ystr, T **M1, T **M2, T **log_M1, T **log_M2, char *alphabet, int m, int n,char *title1,char* title2,  /*{{{*/ 
        char* alignXstr0, char* alignYstr0, int* alignRel0,	
        AlignFactor *pAlignFactor, float gapOpen, float gapExt,int seqtype /*= AA_SEQ */)
/*****************************************************************************
 * gapOpen and gapExt use float number, and thus matrix values are 
 * calculated in float
 ****************************************************************************/
{
    /* first we reserve memory for all the variables we will use */
    int i,j;
    int Imax,Jmax;
    int iAlign;
    float maxValue, tmp;
    int M = m+1;
    int N = n+1;
    int NALIGN = n+m+1; /* maximum size after alignment*/ 
    int sizealphabet = strlen(alphabet);

    Array2D <float> V_2darray(M,N);
    Array2D <int> trace_2darray(M,N);
    float **V   = V_2darray.array2D;     /* value label matrix */
    int **trace = trace_2darray.array2D; /* trace matrix       */

    //Initializing
    {
        for(i = 0 ; i < M ; i ++)
            for(j = 0 ; j < N ;j ++)
            {
                V[i][j] = 0.0;
                trace[i][j] = STOP;
            }
    }

    Array1D <int8> X_1darray(M);
    Array1D <int8> Y_1darray(N);
    Array1D <int8> alignX_1darray(NALIGN);
    Array1D <int8> alignY_1darray(NALIGN);
    Array1D <char> alignXstr_1darray(NALIGN);
    Array1D <char> alignYstr_1darray(NALIGN);
    Array1D <int>  alignRel_1darray(NALIGN);

    int8 *X         = X_1darray.array1D;         /* digital X sequence         */
    int8 *Y         = Y_1darray.array1D;         /* digital Y sequence         */
    int8 *alignX    = alignX_1darray.array1D;    /* digital aligned X sequence */
    int8 *alignY    = alignY_1darray.array1D;    /* digital aligned Y sequence */
    char *alignXstr = alignXstr_1darray.array1D; /* aligned X sequence         */
    char *alignYstr = alignYstr_1darray.array1D; /* aligned Y sequence         */
    int  *alignRel  = alignRel_1darray.array1D;  /* alignment relation         */


    for(i = 0 ; i < m ; i ++) X[i] = Charcase2Digit(Xstr[i], alphabet, sizealphabet);
    for(i = 0 ; i < n ; i ++) Y[i] = Charcase2Digit(Ystr[i], alphabet, sizealphabet);

    // initializing
    for(i=0;i<=m;i++) V[i][0] = i*gapExt; /* do this for i=0,1,...,m */
    for(j=0;j<=n;j++) V[0][j] = j*gapExt;// for the first row and first column, set 0

    for(i=0;i<=m;i++) trace[i][0]=STOP;
    for(j=0;j<=n;j++) trace[0][j]=STOP;

    /* labeling of all nodes, this is the main loop of the algorithm  */
    float e =0.0,f =0.0,g = 0.0; /* e --> up, f --> left, g --> diagnol*/ 
    float s = 0.0;
    for(i = 1 ; i <= m ; i ++)                 // matrix ------- seqY
    {    /* note: we begin at i=1 ! */         //        |
        for(j = 1 ; j <= n ; j ++)             //        |
        {                                      //        seqX
            //s = T (SM[X[i-1]][Y[j-1]]);
            s = float(GetScore_PPAlignment(M1[i-1], M2[j-1], log_M1[i-1], log_M2[j-1]));

            if(trace[i-1][j] != UP)
                e = V[i-1][j] + gapOpen + gapExt ;
            else
                e = V[i-1][j] + gapExt ;

            if(trace[i][j-1] != LEFT)
                f = V[i][j-1] + gapOpen + gapExt ;
            else
                f = V[i][j-1]  + gapExt ;			

            g = V[i-1][j-1] + s ;

            tmp = f ;
            trace[i][j] = LEFT ;

            if(e> tmp)
            {
                tmp = e ;
                trace[i][j] = UP ;
            }
            if(g > tmp)
            {
                tmp = g ;
                trace[i][j] = DIAG ;
            }
            V[i][j] = tmp ;
        }
    }

    // find the highest score in bottom row and rightmost column of the value
    // matrix
    maxValue = MIN_FLOAT ;
    Imax     = 0;
    Jmax     = 0;
    for ( i = 1; i <= m ; i ++ )  /* scan the rightmost column*/ 
    {
        if( V[i][n] > maxValue)
        {
            Imax = i ;
            Jmax = n ;
            maxValue = V[i][n];
        }
    }
    for ( j = 1; j <= n ; j ++ )  /* scan the bottom row*/ 
    {
        if( V[m][j] > maxValue)
        {
            Imax = m ;
            Jmax = j ;
            maxValue = V[m][j] ;
        }
    }

    //[>{{{<]
    //    FILE *fplog = fopen("log.dat","w");
    //    if(1)
    //    {
    //      fprintf(fplog, "alignment value matrix\n");;
    //      fprintf(fplog,"%4c %4c", ' ', ' '); for(j = 0 ; j < n ; j ++) fprintf(fplog, "%4c ", Ystr[j]);fprintf(fplog,"\n");
    //      for( i = 0 ; i <= m ; i ++)
    //      {   
    //          if(i == 0) fprintf(fplog,"%4c ", ' ');
    //          else       fprintf(fplog,"%4c ",Xstr[i-1]);
    //          for(j= 0 ; j <= n ; j ++)
    //          {
    //              fprintf(fplog,"%4.1f ",V[i][j]);
    //          }
    //          fprintf(fplog,"\n");		
    //      }
    //      fprintf(fplog,"\n\n");
    //      fprintf(fplog,"trace matrix \n");
    //      fprintf(fplog,"%4c %4c", ' ', ' '); for(j = 0 ; j < n ; j ++) fprintf(fplog, "%4c ", Ystr[j]);fprintf(fplog,"\n");
    //      for( i = 0 ; i <= m ; i ++)
    //      {
    //          if(i == 0) fprintf(fplog,"%4c ", ' ');
    //          else       fprintf(fplog,"%4c ",Xstr[i-1]);
    //          for(j= 0 ; j <= n ; j ++)
    //          {
    //              fprintf(fplog,"%4d ",trace[i][j]);
    //          }
    //          fprintf(fplog,"\n");
    //      }
    //      fprintf(fplog,"\n\n");
    //    }
    //    fprintf(fplog,"Imax =%d\njmax = %d\nmaxValue = %f\n",Imax,Jmax,maxValue);
    //    fclose(fplog);
    //[>}}}<]

    /* creating aligned sequences alignY and alignX */
    /* note: they are created in reverse order! */
    i = m;
    j = n;
    iAlign = 0; // the length of alignment

    //unaligned ends 

    while(i>Imax) 
    {      
        alignY[iAlign] = GAP_DIGIT; // -1 --> '-'      
        alignX[iAlign] = X[i-1];      
        i --;
        iAlign ++;
    }

    while(j > Jmax) 
    {
        alignY[iAlign] = Y[j-1];
        alignX[iAlign] = GAP_DIGIT;
        j --;
        iAlign ++;
    }

    /* when we come here we know that i==Imax and j==Jmax */
    /* it is from this position we make the jump to the virtual stop node */

    while(trace[i][j] != STOP) 
    {
        if(trace[i][j] == DIAG) 
        {
            alignY[iAlign]=Y[j-1];
            alignX[iAlign]=X[i-1];
            i--;
            j--;
            iAlign++;
        }
        else if(trace[i][j] == LEFT) 
        {
            alignY[iAlign]=Y[j-1];
            alignX[iAlign]= GAP_DIGIT;
            j--;
            iAlign++;
        }
        else if(trace[i][j] == UP) 
        {
            alignY[iAlign]= GAP_DIGIT;
            alignX[iAlign]=X[i-1];
            i--;
            iAlign++;
        }
    }
    /* unaligned beginning */
    while(i>0) 
    {  
        alignY[iAlign]= GAP_DIGIT;
        alignX[iAlign]=X[i-1];
        i--;
        iAlign++;
    }

    while(j>0) 
    {
        alignY[iAlign]=Y[j-1];
        alignX[iAlign]= GAP_DIGIT;
        j--;
        iAlign++;
    }

    AlignAna_Profile(alignX,alignY,iAlign,M1,M2, log_M1, log_M2, pAlignFactor, alignRel);
    pAlignFactor->identity_short = pAlignFactor->identity * iAlign / min(m,n) ;
    pAlignFactor->similarity_short = pAlignFactor->similarity * iAlign / min(m,n) ;
	pAlignFactor->score = float(V[Imax][Jmax]);

    // translate digital aligned sequences back into ascii characters
    for(i = 0 ; i < iAlign ; i ++) 
    {    
        if(alignX[i] == GAP_DIGIT) 
        {
            alignXstr[i] = GAP_CHAR;
        }
        else alignXstr[i] = Digit2Char(alignX[i], alphabet, sizealphabet);
    }
    for(i = 0 ; i < iAlign ; i ++) 
    {    
        if(alignY[i] == GAP_DIGIT) 
        {
            alignYstr[i] = GAP_CHAR;
        }
        else alignYstr[i] = Digit2Char(alignY[i], alphabet, sizealphabet);
    }

    int cnt = 0 ;
    for(i = iAlign-1 ; i >= 0 ; i -- )
    {
        alignXstr0[cnt] = alignXstr[i];
        alignYstr0[cnt] = alignYstr[i];
        alignRel0[cnt] = alignRel[i];
        cnt ++;
    }
    alignXstr0[cnt] = '\0';
    alignYstr0[cnt] = '\0';

    return iAlign;
}
template int Alignment_Profile<int>    (char *Xstr, char *Ystr, int ** M1  , int **M2   , int **log_M1   , int **log_M2   , char *alphabet, int m, int n, char *title1, char* title2, char* alignXstr0, char* alignYstr0, int* alignRel0, AlignFactor *pAlignFactor, float gapOpen, float gapExt, int seqtype /*= AA_SEQ */);
template int Alignment_Profile<float>  (char *Xstr, char *Ystr, float **M1 , float **M2 , float **log_M1 , float **log_M2 , char *alphabet, int m, int n, char *title1, char* title2, char* alignXstr0, char* alignYstr0, int* alignRel0, AlignFactor *pAlignFactor, float gapOpen, float gapExt, int seqtype /*= AA_SEQ */)  ;
template int Alignment_Profile<double> (char *Xstr, char *Ystr, double **M1, double **M2, double **log_M1, double **log_M2, char *alphabet, int m, int n, char *title1, char* title2, char* alignXstr0, char* alignYstr0, int* alignRel0, AlignFactor *pAlignFactor, float gapOpen, float gapExt, int seqtype /*= AA_SEQ */);
/*}}}*/
template <class T> void AlignAna(const int8* alignX, const int8* alignY, int length, T **subMatr, AlignFactor* pAlignFactor, int* alignRel)/*{{{*/
{
    int i ;
    int matchCnt = 0;
    int simiCnt  = 0;
    int gapCnt   = 0;
    for(i = length-1 ; i >= 0 ; i -- )
    {
        if(alignX[i] == DIGIT_INDEL || alignY[i] == DIGIT_INDEL)
        {
            gapCnt ++ ;
            alignRel[i] = GAP ;
        }
        else if(subMatr[alignX[i]][alignY[i]] > 0)
        {
            simiCnt ++ ;
            alignRel[i] = SIM ;
            if(alignX[i] == alignY[i])
            {
                matchCnt ++ ;
                alignRel[i] = IDT ;
            }
        }
        else
            alignRel[i] = MIS ;
    }
    pAlignFactor->idt_cnt    = matchCnt;
    pAlignFactor->sim_cnt    = simiCnt;
    pAlignFactor->gap_cnt    = gapCnt;
    pAlignFactor->identity   = float(matchCnt) / length;
    pAlignFactor->similarity = float((simiCnt - matchCnt)*1.0 + matchCnt) / length;
    pAlignFactor->gapPercent = float(gapCnt)   / length;
}
template  void AlignAna<int>   (const int8* alignX, const int8* alignY, int length, int **subMatr   , AlignFactor* pAlignFactor, int* alignRel);
template  void AlignAna<float> (const int8* alignX, const int8* alignY, int length, float **subMatr , AlignFactor* pAlignFactor, int* alignRel);
template  void AlignAna<double>(const int8* alignX, const int8* alignY, int length, double **subMatr, AlignFactor* pAlignFactor, int* alignRel);
/*}}}*/
template <class T> void AlignAna_Profile(const int8* alignX, const int8* alignY, int length, T **M1, T **M2, T **log_M1, T **log_M2, AlignFactor* pAlignFactor, int* alignRel)/*{{{*/
{
    int i ;
    int matchCnt = 0;
    int simiCnt  = 0;
    int gapCnt   = 0;
    int cntXSeq = 0; // count the amino acid seq in alignX, that is neglect the gap and return the real seqIndex, because M1 and log_M1 need the original sequence index
    int cntYSeq = 0;
    for(i = 0 ; i < length ; i ++ )
    {
        if(alignX[i] == DIGIT_INDEL || alignY[i] == DIGIT_INDEL)
        {
            gapCnt ++ ;
            alignRel[i] = GAP ;
        }
        else if(GetScore_PPAlignment(M1[cntXSeq],M2[cntYSeq], log_M1[cntXSeq], log_M2[cntYSeq]) > 0.0)
        {
            simiCnt ++ ;
            alignRel[i] = SIM ;
            if(alignX[i] == alignY[i])
            {
                matchCnt ++ ;
                alignRel[i] = IDT ;
            }
        }
        else
            alignRel[i] = MIS ;
        if(alignX[i] != DIGIT_INDEL) cntXSeq ++;
        if(alignY[i] != DIGIT_INDEL) cntYSeq ++;
    }
    pAlignFactor->idt_cnt    = matchCnt;
    pAlignFactor->sim_cnt    = simiCnt;
    pAlignFactor->gap_cnt    = gapCnt;
    pAlignFactor->identity   = float(matchCnt) / length;
    pAlignFactor->similarity = float((simiCnt - matchCnt)*1.0 + matchCnt) / length;
    pAlignFactor->gapPercent = float(gapCnt)   / length;
}
template void AlignAna_Profile<int>   (const int8* alignX, const int8* alignY, int length, int **M1   , int **M2   , int **log_M1   , int **log_M2   , AlignFactor* pAlignFactor, int* alignRel);
template void AlignAna_Profile<float> (const int8* alignX, const int8* alignY, int length, float **M1 , float **M2 , float **log_M1 , float **log_M2 , AlignFactor* pAlignFactor, int* alignRel);
template void AlignAna_Profile<double>(const int8* alignX, const int8* alignY, int length, double **M1, double **M2, double **log_M1, double **log_M2, AlignFactor* pAlignFactor, int* alignRel);
/*}}}*/

int SeqMatchAlign_Protein(const char *Xstr, const char* Ystr, char *title1,char* title2,char* alignXstr0, char* alignYstr0,int* alignRel0,FILE *fpout, bool isPrintToScreen )/*{{{*/
{
    /* first we reserve memory for all the variables we will use */
    int gapOpen = -7;
    int gapExt = -1 ;
    int match = 5;
    int misMatch = -10;
    int penalty_X_to_Ohter = -1; /*when align atomSeq to aaSeq, X residue might be mutated, the penalty should be low, 2007-11-09, Nanjiang*/

    int i,j;
    int Imax,Jmax;
    int iAlign;
    int maxValue, tmp;
    int m = strlen(Xstr);
    int n = strlen(Ystr);

    int8* X ;
    int8* Y ;

    int** V = NULL ;  /* value label matrix */

    int** trace = NULL;     /* trace matrix */
    char* alignYstr ;
    char* alignXstr ;
    int8* alignY; /* aligned Y sequence */
    int8* alignX; /* aligned X sequence */

    int **subMatr = NULL ; // substitution matrix
    int M = m+1;
    int N = n+1;
    int NALIGN = n+m+1;


    X = new int8[M] ;
    Y = new int8[N] ;
    CharToDigit_Protein(Xstr,X,m);
    CharToDigit_Protein(Ystr,Y,n);

    V     = Create2DArray(V, M , N);
    trace = Create2DArray(trace, M , N);

    alignYstr = new char[NALIGN+1] ;
    alignXstr = new char[NALIGN+1] ;
    alignY = new int8[NALIGN+1] ;
    alignX = new int8[NALIGN+1] ;

    subMatr = Create2DArray(subMatr , NUM_BLOSUM, NUM_BLOSUM);
    for (i = 0 ; i < NUM_BLOSUM; i++)
    {
        for ( j = 0 ; j < NUM_BLOSUM; j++)
        {
            if (i < 20 && j < 20)/*for 20 standard amino acids*/
            {
                if(i == j )
                    subMatr[i][j] = match ;
                else
                    subMatr[i][j] = misMatch;
            }
            else /*for non standard amino acid,e.g X, Z, B*/
            {
               subMatr[i][j] = penalty_X_to_Ohter;
            }
        }
    }
    /*{{{*/
    /*
       cout <<"blosum62 matrix "<< endl;
       cout << endl;
       for ( i = 0 ; i < NUM_BLOSUM ; i++)
       {
       for ( j = 0 ; j < NUM_BLOSUM ; j++)
       printf("%2d ", subMatr[i][j]);
       cout << endl;
       }
       cout << endl;
       */

    //here the real program begins


    //	printf("Ystr: %s\n",Ystr);
    //	printf("Ystr: ");
    //	for (i=0;i<n;i++) cout<<Y[i]<<' ';
    //	cout <<endl;  
    //	printf("Xstr: %s\n",Xstr);
    //	printf("Xstr: ");
    //	for (i=0;i<m;i++) cout<<X[i]<<' ';
    //	cout <<endl;
    //	exit(0);/*}}}*/

    //initialization

    for(i=0;i<=m;i++) V[i][0] = i*gapExt; /* do this for i=0,1,...,m */
    for(j=0;j<=n;j++) V[0][j] = j*gapExt;// for the first row and first column, set 0


    for(i=0;i<=m;i++) trace[i][0]=STOP;
    for(j=0;j<=n;j++) trace[0][j]=STOP;

    maxValue=0;
    Imax=0;
    Jmax=0;

    /* labeling of all nodes, this is the main loop of the algorithm  */

    int e,f,g;
    int s ;
    for(i=1;i<=m;i++) 
    {    /* note: we begin at i=1 ! */
        for(j=1;j<=n;j++)
        {  
            s = subMatr[X[i-1]][Y[j-1]];

            if(trace[i-1][j] != UP)
                e = V[i-1][j] + gapOpen + gapExt ;
            else
                e = V[i-1][j] + gapExt ;

            if(trace[i][j-1] != LEFT)
                f = V[i][j-1] + gapOpen + gapExt ;
            else
                f = V[i][j-1]  + gapExt ;			

            g = V[i-1][j-1] + s ;

            tmp = f ;
            trace[i][j] = LEFT ;

            if(e> tmp)
            {
                tmp = e ;
                trace[i][j] = UP ;
            }

            if(g > tmp)
            {
                tmp = g ;
                trace[i][j] = DIAG ;
            }

            V[i][j] = tmp ;

        }
    }


    // find the highest score
    maxValue = MIN_INT;
    for ( i = 1; i <= m ; i ++ )
    {
        if( V[i][n] > maxValue)
        {
            maxValue = V[i][n] ;
            Jmax = n ;
            Imax = i ;
        }
    }

    for ( j = 1; j <= n ; j ++ )
    {
        if( V[m][j] > maxValue)
        {
            maxValue = V[m][j] ;
            Imax = m ;
            Jmax = j ;
        }
    }

    // always tracing back from the bottom-right corner
    Imax = m ;
    Jmax = n;


    /* now create aligned sequences alignY and alignX */
    /* note: these are created in reverse order! */
    i=m;
    j=n;
    iAlign=0; // the length of alignment

    /* unaligned ends */

    while(i>Imax) 
    {      
        alignY[iAlign] = DIGIT_INDEL; // -1 --> '-'      
        alignX[iAlign] = X[i-1];      
        i --;
        iAlign ++;
    }

    while(j > Jmax) 
    {
        alignY[iAlign] = Y[j-1];
        alignX[iAlign] = DIGIT_INDEL;
        j --;
        iAlign ++;
    }

    /* when we come here we know that i==Imax and j==Jmax */
    /* it is from this position we make the jump to the virtual stop node */

    while(trace[i][j] != STOP) 
    {

        if(trace[i][j] == DIAG) 
        {
            alignY[iAlign]=Y[j-1];
            alignX[iAlign]=X[i-1];
            i--;
            j--;
            iAlign++;
        }
        else if(trace[i][j] == LEFT) 
        {
            alignY[iAlign]=Y[j-1];
            alignX[iAlign]= DIGIT_INDEL;
            j--;
            iAlign++;
        }
        else if(trace[i][j] == UP) 
        {
            alignY[iAlign]= DIGIT_INDEL;
            alignX[iAlign]=X[i-1];
            i--;
            iAlign++;
        }
    }

    /* unaligned beginning */

    while(i>0) 
    {  
        alignY[iAlign]= DIGIT_INDEL;
        alignX[iAlign]=X[i-1];
        i--;
        iAlign++;
    }

    while(j>0) 
    {
        alignY[iAlign]=Y[j-1];
        alignX[iAlign]= DIGIT_INDEL;
        j--;
        iAlign++;
    }


    /* print out the alignment result*/	
    AlignFactor alignFactor ;
    int *alignRel = new int[iAlign+1];
    AlignAna_Protein(alignX,alignY,iAlign,subMatr,&alignFactor, alignRel);

    if (isPrintToScreen) /*added 2007-11-07, Nanjiang*/
    {
        printf("GapOpen      = %d\n" , gapOpen);
        printf("GapExtension = %d\n",gapExt);
        printf("identity     = %3.2f%c\n", alignFactor.identity   * 100,'%');
        printf("similarity   = %3.2f%c\n", alignFactor.similarity * 100,'%');
        printf("gapPercent   = %3.2f%c\n", alignFactor.gapPercent * 100,'%');
        printf("identity of shorter seq= %3.2f%c\n", alignFactor.identity * iAlign /min(m,n) * 100,'%');
        printf("%s(%d) <--> %s(%d)\n",title1,m,title2,n);
        printf("\n");
    }

    if( fpout != NULL)  /*added 2007-11-07, Nanjiang*/ 
    {
        fprintf(fpout,"GapOpen      = %d\n" , gapOpen);
        fprintf(fpout,"GapExtension = %d\n",gapExt);
        fprintf(fpout,"identity     = %3.2f%c\n", alignFactor.identity * 100,'%');
        fprintf(fpout,"similarity   = %3.2f%c\n", alignFactor.similarity * 100,'%');
        fprintf(fpout,"gapPercent   = %3.2f%c\n", alignFactor.gapPercent * 100,'%');
        fprintf(fpout,"identity of shorter seq= %3.2f%c\n", alignFactor.identity * iAlign /min(m,n) * 100,'%');
        fprintf(fpout,"%s(%d) <--> %s(%d)\n",title1,m,title2,n);
        fprintf(fpout,"\n");
    }

    DigitToChar_Protein(alignY,iAlign,alignYstr);
    DigitToChar_Protein(alignX,iAlign,alignXstr);

    int cnt = 0 ;
    for(i = iAlign-1 ; i >= 0 ; i -- )
    {
        alignYstr0[cnt] = alignYstr[i];
        alignXstr0[cnt] = alignXstr[i];
        alignRel0[cnt] = alignRel[i];
        cnt ++;
    }

    int lineLength = 60 ;	
    if(fpout != NULL)   /*added 2007-11-07, Nanjiang*/ 
    {
        WriteAlignment(title1,title2,alignXstr0,alignYstr0,alignRel0,iAlign, lineLength,fpout);
    }

    //free memory
    delete [] Y ;
    delete [] X ;
    Delete2DArray(V, M);
    Delete2DArray(trace, M);
    delete [] alignYstr ;
    delete [] alignXstr ;
    delete [] alignY ;
    delete [] alignX ;
    Delete2DArray(subMatr , NUM_BLOSUM);
    delete [] alignRel ;

    return iAlign;
}
/*}}}*/
void WriteAlignmentHeader(float gapOpen, float gapExt, AlignFactor *pAlignFactor, char *title1, char *title2, int m, int n, FILE *fpout)/*{{{*/
{
    // write header
    fprintf(fpout,"# GapOpen      = %3.1f\n" , -gapOpen);
    fprintf(fpout,"# GapExtension = %3.1f\n" , -gapExt);
    fprintf(fpout,"# Identity     = %3.2f%c\n", pAlignFactor->identity * 100,'%');
    fprintf(fpout,"# Similarity   = %3.2f%c\n", pAlignFactor->similarity * 100,'%');
    fprintf(fpout,"# GapPercent   = %3.2f%c\n", pAlignFactor->gapPercent * 100,'%');
    fprintf(fpout,"# Identity of the shorter sequence = %3.2f%c\n", pAlignFactor->identity_short * 100,'%');
    fprintf(fpout,"# %s (%d) <--> %s (%d) rawScore = %g Evalue = %lg\n",title1,m,title2,n, pAlignFactor->score, pAlignFactor->eValue);
    fprintf(fpout,"# \n");
}
/*}}}*/
void WriteAlignment(const char* title1,const char* title2,const char* aXstr,const char* aYstr, int* aRel, int length, int lineLength,FILE* fpout, int outmode /*= HORIZONTAL */)/*{{{*/
    /*****************************************************************************
     * Write Alignment horizontally or vertically
     * if writing horizongtally, lineLength must be set, which is the number of
     * characters output for each line
     * format
     *
     * Horizongtal: 4 lines
     *           01234567890123456                  60     // ruler for aligned seq | accumulated length for aligned seq
     * X:     1  MST-HHSSGHHGG....                  59     // sequence X | accumulated length for sequence X
     *           |.| |  ||||. ....
     * Y:     1  MHTHH--SGHHKK.....                 58     // sequence Y | accumulated length for sequence Y
     *
     * vertical: 6 columns
     * # column1: indexing for alignment
     * # column2: indexing for sequence X
     * # column3: sequence X
     * # column4: similarity symbol, '|': identical, '.':similar, ' ':not-similar 
     * # column5: sequence Y
     * # column6: indexing for sequence Y
     *
     * 1   1 M = M   1
     * 2   2 S - H   2
     * 3   3 T = T   3
     * 4   . |   H   4
     * 5   4 H = H   5
     * 6   5 H   |   .
     * 7   6 S   |   .
     * 8   7 S = S   6
     * .
     ****************************************************************************/
{
    int i;
    int size_title = max(strlen (title1), strlen(title2));
    if(outmode == HORIZONTAL)
    {
        int col   = 0;
        int col1  = 0; // count for sequence 1, without gaps
        int col2  = 0; // count for sequence 2, without gaps
        int count = 0;
        int lineCnt = 0 ;
        int lineBegin = 0 , lineEnd = 0 ;
        while(lineEnd < length)
        {	
            lineBegin = lineCnt * lineLength ;
            lineEnd = min((lineCnt + 1) * lineLength,length);
            //-------------------------------------------  line 1
            fprintf(fpout,"%-*s  %4s ",size_title,"","");	
            //count = 0;
            for(i = lineBegin ;i < lineEnd ; i ++ ) 
            {
                col ++ ;
                count ++ ;
                if(count > 9) count = 0 ;
                fprintf(fpout,"%1d",count);
            }
            fprintf(fpout," %4d",col);		
            fprintf(fpout,"\n");	
            //-------------------------------------------- line 2
            col1 ++;
            fprintf(fpout,"%-*s: %4d ",size_title, title1,col1);
            for(i = lineBegin ;i < lineEnd ; i ++ ) 
            {
                fprintf(fpout,"%c",aXstr[i]);
                if(aXstr[i] != GAP_CHAR) col1 ++;
            }
            col1 --;
            fprintf(fpout," %4d",col1);
            fprintf(fpout,"\n");
            //-------------------------------------------- line 3
            fprintf(fpout,"%-*s  %4s ",size_title,"","");
            for(i = lineBegin ;i < lineEnd ; i ++ ) 
            {
                if(aRel[i] == IDT)
                {
                    fprintf(fpout,"%c",'|');
                }
                else if(aRel[i] == SIM)
                {
                    fprintf(fpout,"%c",'.');
                }
                else
                {
                    fprintf(fpout,"%c",' ');
                }
            }
            fprintf(fpout,"\n");
            //-------------------------------------------- line 4
            col2 ++;
            fprintf(fpout,"%-*s: %4d ",size_title,title2,col2);
            for(i = lineBegin ;i < lineEnd ; i ++ ) 
            {
                fprintf(fpout,"%c",aYstr[i]);
                if(aYstr[i] != GAP_CHAR) col2++;
            }
            col2 --;
            fprintf(fpout," %4d",col2);
            fprintf(fpout,"\n");
            //--------------------------------------------
            fprintf(fpout,"\n\n");
            lineCnt ++ ;
        }
        fprintf(fpout,"\n");

    }
    else//vertical 
    {
        char IDT_C = '=';
        char SIM_C = '-';
        char MIS_C = ' ';
        char GAP_CHAR_VER = '|';
        char GAP_CHAR_INDEX = '.';
        fprintf(fpout, "# column 1: indexing for the alignment.\n");
        fprintf(fpout, "# column 2: indexing for sequence \"%s\". '%c' means gap.\n",title1, GAP_CHAR_INDEX);
        fprintf(fpout, "# column 3: sequence \"%s\". '%c' means gap.\n", title1, GAP_CHAR_VER);
        fprintf(fpout, "# column 4: similarity symbol. '%c':identical, '%c':similar, '%c':not-similar.\n", IDT_C, SIM_C, MIS_C);
        fprintf(fpout, "# column 5: sequence \"%s\". '%c' means gap.\n",title2, GAP_CHAR_VER);
        fprintf(fpout, "# column 6: indexing for sequence \"%s\". '%c' means gap.\n", title2, GAP_CHAR_INDEX);
        int cnt1 = 0;
        int cnt2 = 0;
        char symbol = NULLCHAR;
        for( i = 0 ; i < length ; i++)
        {
            //column1, indexing of aligned sequence
            fprintf(fpout,"%4d ",i+1);
            //column2, indexing of sequece X
            if(aXstr[i] != GAP_CHAR)  
            {
                cnt1 ++;
                fprintf(fpout,"%4d ", cnt1);
            }
            else
                fprintf(fpout,"%4c ", GAP_CHAR_INDEX);
            //column3, sequence X
            if (aXstr[i] != GAP_CHAR) fprintf(fpout,"%1c ",aXstr[i]);
            else fprintf(fpout,"%1c ", GAP_CHAR_VER);
            //column4, alignment symbol
            if(aRel[i] == IDT) symbol = IDT_C;
            else if(aRel[i] == SIM) symbol = SIM_C; 
            else symbol = MIS_C;
            fprintf(fpout,"%1c ",symbol);
            //column5, sequence Y
            if (aYstr[i] != GAP_CHAR) fprintf(fpout,"%1c ",aYstr[i]);
            else fprintf(fpout,"%1c ", GAP_CHAR_VER);
            //column6, indexing of longer sequence
            if(aYstr[i] != GAP_CHAR)  
            {
                cnt2 ++;
                fprintf(fpout,"%4d", cnt2);
            }
            else
                fprintf(fpout,"%4c", GAP_CHAR_INDEX);
            //end of line
            fprintf(fpout,"\n");
        }
    }
}
/*}}}*/

/*functions for scop */
int FindSCOP(char* pdbid, char chainID, int seqF, int seqT, SCOP* pSCOP, FILE *fpSCOPspi, FILE *fpSCOPspd, char** pdbIDs, long* offsetspi, int numPDB) /*{{{*/ 
/*****************************************************************************
 * modified 2007-07-23 
 * seqF: the start position of domain                        
 * seqT: the end position of domain                          
 * pSCOP: pointer of the SCOP structure                      
 * fpSCOPspi: file pointer for the spi file                  
 * fpSCOPspd: file pointer for the spd file                  
 * pdbIDs: array for the pdbid                               
 * offsetspi: array of the offset for each pdbid in spi file 
 * numPDB: number of items in the spi file                   
 ****************************************************************************/
{
    int status_fseek = 0;
	int idx ;
	char str[100] = "";
	int numDomain;
	int numDomainMatch;
	int cnt;
	int i ,j;
    int linesize;
    int maxline = 500;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;

    int status_fscanf = 0;

    bool isUsingWholeChain = false;
    if(seqF == 0 && seqT == 0x7FFFFFFF) 
    {
        isUsingWholeChain = true;   
    }

	if((idx = BinarySearch_String(pdbid,pdbIDs,numPDB)) != -1)
	{
		if((status_fseek = fseek(fpSCOPspi,offsetspi[idx], SEEK_SET))!= 0)
			fprintf(stderr,"fseek faild \t stream = %s \t offset = %ld\n","fpSCOPspi",offsetspi[idx]);
#ifdef ASSERT
        assert(status_fseek == 0);
#endif
		status_fscanf = fscanf(fpSCOPspi,"%s %d",str,&numDomain);
        if (status_fscanf != 2)
        {
            fprintf(stderr,"fscanf error! prog:%s.\n","FindSCOP");
        }
		char domName[SIZE_SCOP_ID+1] = "";

		long offset = 0;
		Array2D <char> domIDs_2darray(numDomain,SIZE_SCOP_ID+1);
		Array1D <long> offsetspd_1darray(numDomain);
        char ** domIDs = domIDs_2darray.array2D;
        long * offsetspd = offsetspd_1darray.array1D;

		cnt = 0 ;
#ifdef DEBUG
        printf("seqF=%d\n", seqF);
        printf("seqT=%d\n", seqT);
        printf("numDomain=%d\n", numDomain);
        printf("isUsingWholeChain=%d\n", isUsingWholeChain);
        printf("chainID=%c,\n", chainID);
#endif
		for(i =0 ; i < numDomain ; i++)
		{
			status_fscanf = fscanf(fpSCOPspi,"%s %ld",domName,&offset);
			if(toupper(domName[5]) == toupper(chainID) || domName[5] == '.' || chainID == ' ')
			{
				my_strcpy(domIDs[cnt],domName, SIZE_SCOP_ID);
				offsetspd[cnt] = offset;
				cnt ++;
			}
		}

		SCOP scop;
        InitSCOP(&scop);
		int coverage; // the sequence coverage of the query domain and scop domain
		numDomainMatch = 0 ;
		for(i = 0 ; i < cnt ; i++)
		{
			if((status_fseek = fseek(fpSCOPspd,offsetspd[i],SEEK_SET) )!= 0)
				fprintf(stderr, "fseek faild: stream = %s \t offset = %ld\n","fpSCOPspd",offsetspd[i]);
#ifdef ASSERT
            assert(status_fseek == 0);
#endif
            linesize = fgetline(fpSCOPspd, line, maxline);
			ScanfSCOPRecord(line, &scop);
            bool isMatch = false;

            for(j = 0 ; j < scop.domDef.numChain ; j++)
            {
                if(chainID == scop.domDef.chainIDs[j] || chainID == ' ')
                {
                    coverage = Coverage(seqF,seqT,scop.domDef.posF[j],scop.domDef.posT[j])+1;
                    if(scop.domDef.isWholeChain[j] == true || isUsingWholeChain ||
                            (double(coverage) / (seqT-seqF) >= 0.5 /*midified 2007-07-23, new rule is: if >half of the sequence is within the SCOP domain, consider this sequence belonged to that domain*/
                             || double(coverage) / (scop.domDef.posT[j]-scop.domDef.posF[j]) >= 0.5
                             )
                           )
                    {
                        isMatch = true;
                        break;
                    }
#ifdef DEBUG
        printf("coverage=%d\n", coverage);
#endif
                }
            }
            if(isMatch)
            {
                CopySCOP(&pSCOP[numDomainMatch],&scop);
                numDomainMatch ++;
            }
		}

		return numDomainMatch;
	}
	else
		return -1;
}
/*}}}*/
int FindSCOP_scopid(char* scopid, SCOP* pSCOP, FILE *fpSCOPspi, FILE *fpSCOPspd, char** pdbIDs, long* offsetspi, int numPDB) /*{{{*/ 
/*****************************************************************************
 * 2007-07-10
 * pSCOP: pointer of the SCOP structure                      
 * fpSCOPspi: file pointer for the spi file                  
 * fpSCOPspd: file pointer for the spd file                  
 * pdbIDs: array for the pdbid                               
 * offsetspi: array of the offset for each pdbid in spi file 
 * numPDB: number of items in the spi file                   
 *
 * ChangeLog 2010-04-12: For scopid starts with g, e.g. "g1r8o.1", it is
 * converted to "d1r8o.1"
 ****************************************************************************/
{
    int status_fseek = 0;
    int status_fscanf = 0;
	int idx ;
	char str[100] = "";
	int numDomain;
	int numDomainMatch;
	int cnt;
	int i ;
    //bool isUsingWholeChain = false;
    char pdbid[SIZE_PDBID+1] = "";

    int linesize;
    int maxline = 500;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;

    int size_scopid = strlen(scopid);
    Array1D <char> tmpSCOPID_1darray(size_scopid + 1);
    char  *tmpSCOPID = tmpSCOPID_1darray.array1D;
    my_strcpy(tmpSCOPID, scopid, size_scopid);

    if (tmpSCOPID[0] == 'g') {
        tmpSCOPID[0] = 'd'; /*Changed 2010-04-12*/
    }

    my_strcpy(pdbid, tmpSCOPID+1, SIZE_PDBID);
    my_strupr(pdbid);
#ifdef DEBUG
    fprintf(stdout,"pdbid=%s, tmpSCOPID=%s\n", pdbid, tmpSCOPID);
#endif

	if((idx = BinarySearch_String(pdbid,pdbIDs,numPDB)) != -1)
	{
		if((status_fseek = fseek(fpSCOPspi,offsetspi[idx], SEEK_SET))!= 0)
			fprintf(stderr,"fseek faild \t stream = %s \t offset = %ld\n","fpSCOPspi",offsetspi[idx]);
#ifdef ASSERT
        assert(status_fseek == 0);
#endif
		status_fscanf = fscanf(fpSCOPspi,"%s %d",str,&numDomain);
        if (status_fscanf != 2)
        {
            fprintf(stderr,"fscanf faild! function:%s.\n","FindSCOP_scopid");
        }
		char domName[SIZE_SCOP_ID+1] = "";

		long offset = 0;
		Array2D <char> domIDs_2darray(numDomain,SIZE_SCOP_ID+1);
		Array1D <long> offsetspd_1darray(numDomain);
        char ** domIDs = domIDs_2darray.array2D;
        long * offsetspd = offsetspd_1darray.array1D;

		cnt = 0 ;
		for(i =0 ; i < numDomain ; i++)
		{
			status_fscanf = fscanf(fpSCOPspi,"%s %ld",domName,&offset);
            if (status_fscanf != 2)
            {
                fprintf(stderr,"fscanf faild! function:%s.\n","FindSCOP_scopid");
            }
			if(strcasecmp(domName, tmpSCOPID) == 0)
			{
				my_strcpy(domIDs[cnt],domName, SIZE_SCOP_ID);
				offsetspd[cnt] = offset;
				cnt ++;
			}
		}

		SCOP scop;
        InitSCOP(&scop);
		//int coverage; // the sequence coverage of the query domain and scop domain
		numDomainMatch = 0 ;
		for(i = 0 ; i < cnt ; i++)
		{
			if((status_fseek = fseek(fpSCOPspd,offsetspd[i],SEEK_SET) )!= 0)
				fprintf(stderr, "fseek faild: stream = %s \t offset = %ld\n","fpSCOPspd",offsetspd[i]);
#ifdef ASSERT
            assert(status_fseek == 0);
#endif
            linesize = fgetline(fpSCOPspd, line, maxline);
			ScanfSCOPRecord(line, &scop);
            CopySCOP(&pSCOP[numDomainMatch],&scop);
            numDomainMatch ++;
		}

		return numDomainMatch;
	}
	else
		return -1;
}
/*}}}*/

/*misc*/
double GetMergeRatio1(int aaSeqIndex, float *p1, float *p2, int length)/*{{{*/
/*****************************************************************************
 * get the ratio to merge matFrag and matMODM
 * p1: parameter 1
 * p2: parameter 2
 * the return ratio 0-1, is on parameter 1
 ****************************************************************************/
{
    int win = 3;/*window size*/
    int halfwin = win/2;
    double ratio = 0.0;
    double sum1, sum2;
    sum1=sum2=0.0;
    for(int i = aaSeqIndex-halfwin; i <= aaSeqIndex+halfwin; i ++)
    {
        if(i < 0 || aaSeqIndex >= length) continue;
        sum1 += p1[i];
        sum2 += p2[i];
    }
    if(sum1 <= 0.0 ) { ratio = 0.0; }
    else if (sum2 <=0.0) { ratio = 1.0;}
    else { ratio = sum1 /(sum1+sum2);}
    return ratio;
}
/*}}}*/
double GetMergeRatio2(int aaSeqIndex, float *p1, float *p2, int length)/*{{{*/
/*****************************************************************************
 * get the ratio to merge matFrag and matMODM
 * p1: parameter 1, for matFrag
 * p2: parameter 2, for matMODM
 * the return ratio 0-1, is on parameter 1
 * For GetMergeRatio, the most ratio is varying from 0.4-0.6 
 * new scheme: lower the ratio on fragacc profiles
 *
 * 2009-06-17
 ****************************************************************************/
{
    int win = 3;/*window size*/
    int halfwin = win/2;
    double ratio = 0.0;
    double sum1, sum2;
    sum1=sum2=0.0;
    for(int i = aaSeqIndex-halfwin; i <= aaSeqIndex+halfwin; i ++)
    {
        if(i < 0 || aaSeqIndex >= length) continue;
        sum1 += p1[i];
        sum2 += p2[i];
    }
    if(sum1 <= 0.0 ) { ratio = 0.0; }
    else if (sum2 <=0.0) { ratio = 0.5;}/*at most 0.5 on matFrag*/
    else 
    { 
        if (sum2 < 0.25*win)
        {
            sum2 *= 1.5;

        }
        else if (sum2 < 0.5 *win)
        {
            sum2 *= 2.0;

        }
        else if (sum2 < 1.0*win)
        {
            sum2 *= 3.0;
        }
        else
        {
            sum2 *= 4.0;
        }
        ratio = sum1 /(sum1+sum2);
    }
    return ratio;
}
/*}}}*/
double GetMergeRatio3(int aaSeqIndex, float *p1, float *p2, int length)/*{{{*/
/*****************************************************************************
 * get the ratio to merge matFrag and matMODM
 * p1: parameter 1, for matFrag
 * p2: parameter 2, for matMODM
 * the return ratio 0-1, is on parameter 1
 * For GetMergeRatio, the most ratio is varying from 0.4-0.6 
 * new scheme: lower the ratio on fragacc profiles
 *
 * 2009-06-17
 ****************************************************************************/
{
    int win = 3;/*window size*/
    int halfwin = win/2;
    double ratio = 0.0;
    double sum1, sum2;
    sum1=sum2=0.0;
    for(int i = aaSeqIndex-halfwin; i <= aaSeqIndex+halfwin; i ++)
    {
        if(i < 0 || aaSeqIndex >= length) continue;
        sum1 += p1[i];
        sum2 += p2[i];
    }
    if(sum1 <= 0.0 ) { ratio = 0.0; }
    else if (sum2 <=0.0) { ratio = 0.4;}/*at most 0.4 on matFrag*/
    else 
    { 
        if (sum2 < 0.25*win)
        {
            sum2 *= 1.5;

        }
        else if (sum2 < 0.5 *win)
        {
            sum2 *= 2.0;

        }
        else if (sum2 < 1.0*win)
        {
            sum2 *= 3.0;
        }
        else if (sum2 < 1.5*win)
        {
            sum2 *= 5.0;
        }
        else
        {
            sum2 *= 10.0; /*else use nper = 0*/
        }
        ratio = sum1 /(sum1+sum2+1e-6);
        if (ratio > 0.4)
        {
            ratio = 0.4; /*set the maximum of ratio to 0.4*/
        }
    }
    return ratio;
}
/*}}}*/
char *GetQijFilePath(const char *rmtID, char *filename, const char *path, int qijformat /*= QIJ_FORMAT_NANJIANG*/, bool isReadBinaryFile /*= true*/)/*{{{*/
{
    if (isReadBinaryFile)
    {
        if(qijformat == QIJ_FORMAT_TUPING)
            sprintf(filename, "%s/Qijmatrix_%s.txtbin", path, rmtID);
        else
            sprintf(filename,"%s/%s.Qijbin", path, rmtID);
    }
    else 
    {
        if(qijformat == QIJ_FORMAT_TUPING)
            sprintf(filename, "%s/Qijmatrix_%s.txt", path, rmtID);
        else
            sprintf(filename,"%s/%s.Qij", path, rmtID);
    }
    return filename;
}
/*}}}*/
char *GetFragMatFilePath(const char *rmtID, char *filename, const char *path, int fragaccformat /*= FRAGACC_FORMAT_NANJIANG*/, bool isReadBinaryFile /*= true*/)/*{{{*/
{
    if (isReadBinaryFile)
    {
        if (fragaccformat == FRAGACC_FORMAT_TUPING)
        { sprintf(filename, "%s/frag_%s.txtbin", path, rmtID); }
        else
        { sprintf(filename, "%s/%s.fragaccbin", path, rmtID); }
    }
    else
    { 
        if (fragaccformat == FRAGACC_FORMAT_TUPING)
        {   sprintf(filename, "%s/frag_%s.txt", path, rmtID);}
        else
        {   sprintf(filename, "%s/%s.fragacc", path, rmtID);}
    }
    return filename;
}
/*}}}*/
char *GetMODMFilePath_Tuping(const char *rmtID, char *filename, const char *path, int modmformat /*= MODM_FORMAT_NANJIANG*/, bool isReadBinaryFile /*= true*/)/*{{{*/
{
    if(isReadBinaryFile)
    {
        if(modmformat == MODM_FORMAT_TUPING)
            sprintf(filename, "%s/modmatrix_%s.txtbin", path, rmtID);
        else 
            sprintf(filename, "%s/%s.modmbin", path, rmtID);
    }
    else
    {
        if(modmformat == MODM_FORMAT_TUPING)
            sprintf(filename, "%s/modmatrix_%s.txt", path, rmtID);
        else 
            sprintf(filename, "%s/%s.modm", path, rmtID);
    }
    return filename;
}
/*}}}*/
void ReorderMatrix(int **M, const char *alphabet_from, const char *alphabet_to)/*{{{*/
{
    int sizeAlphabet = strlen(alphabet_from);
    if ( strlen(alphabet_from) != strlen(alphabet_to))
    {
        fprintf(stderr, "Error: strlen(alphabet_to(%s))=%lu Not-Equal-To strlen(alphabet_from(%s)=%lu)\n", alphabet_to, strlen(alphabet_to), alphabet_from, strlen(alphabet_from));
        assert(strlen(alphabet_from) == strlen(alphabet_to));
    }

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
template <class T> void TreatAllZeroFij(int **fij, int length, T *score1, char *aaSeq, const char* alphabet)/*{{{*/
/*****************************************************************************
 * For all zero profiles, set the fij on amino acid of the position to 100
 * 2008-01-15, Nanjiang
 * updated 2009-07-27
 ****************************************************************************/
{
    int i,j;
    int sumFij = 0;
    for(i = 0;i < length; i ++) {
        sumFij = Sum(fij[i], 0, 19);
        //fprintf (stdout,"%d, sumFij=%d\n", i, sumFij);
        if((sumFij < 50 ||score1[i] < 0.0) && ( sumFij < 80)) {/*in such conditions, modify the profile, 2009-07-27*/
            int daa = Char2Digit(aaSeq[i], alphabet);
            if (daa >= 0 && daa < 20) {/*if this amino acid is within the standard 20 amino acids*/
                for (j = 0;j<20;j++) {
                    if (j == daa) {
                        fij[i][j] = 100; /*if the amino acid name is one of the 20 standard amino acid and the MODM profile is all zero, set the freqency on that residue to 100 %*/
                    } else {
                        fij[i][j] = 0;
                    }
                }
            } else  { /*if non standard amino acids*/ 
                for(j = 0;j < 20;j++) {
                    fij[i][j] = Integer(Background_AA_Freq[j]*100.0) ;
                }
            }
        }
    }
}
template void TreatAllZeroFij<float> (int** fij, int length, float* score1 , char* aaSeq, const char* alphabet);
template void TreatAllZeroFij<double>(int** fij, int length, double* score1, char* aaSeq, const char* alphabet);
/*}}}*/

// functions for sequence handling

int GetAllSubset(int **subset, int n, int k)/*{{{*/
    /*****************************************************************************
     * GetAllSubset()
     * return number of subset
     * for example, give n = 4 ; k = 2
     * subset will store all the subset of (n,k), indexed from 0, that is
     * subset[0] = 0,1
     * subset[1] = 0,2
     * ...
     * subset[6] = 2,3
     * note: do not use this function when n and k are large, e.g. C(100,4) = 3921225
     * which need 60M memory to store subset, C(50,4) need 3.5M, acceptable
     ****************************************************************************/
{
    int *a = new int[k];
    int i ;
    int cnt;
    bool done = true;

    cnt = 0;
    for( ; ; )
    {
        ksub_next4(n,k,a,&done);
        if(done) break;
        for(i = 0 ; i < k ; i ++) subset[cnt][i] = a[i]-1;
        cnt ++;
    }
    delete [] a;
    return cnt;
}
/*}}}*/
