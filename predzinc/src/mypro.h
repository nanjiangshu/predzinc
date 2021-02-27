#ifndef HAS_MYPRO_H
#define HAS_MYPRO_H

#include <set>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <cassert>
#include "myfunc.h"

#include "Constant.h"
#include "DataType.h"

#define _ASSERT_

using namespace std;

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

// type of modm file
#ifndef DATATYPE_MODM_MATRIX
#define DATATYPE_MODM_MATRIX double
#endif

#ifndef DATATYPE_CONSV 
#define DATATYPE_CONSV double  // type of the conervation value, double or int
#endif

#ifndef DATATYPE_SMATRIX
#define DATATYPE_SMATRIX int // datatype for the substitution matrix
#endif


#ifndef MODM_FILE_TYPE
#define MODM_FILE_TYPE
#define MODM_LOG 0
#define MODM_PER 1
#endif

#ifndef ENCODING_TYPE
#define   ENCODING_TYPE
#define   PSSM_PROFILE_ENCODE    0
#define   SEQUENCE_BINARY_ENCODE 1
#endif

#ifndef MSA_FORMAT
#define MSA_FORMAT
#define CLUSTALW   0
#define MSF        1
#define GCG        2
#define FASTA      3
#define PHYLIP     4
#endif

//#ifndef SEQ_TYPE
//#define PROTEIN 0
//#define DNA     1
//#define OTHER   2
//#endif

#ifndef SVM_FILE_TYPE
#define FILE_PRED 0
#define FILE_LABEL 1
#endif

//for reading the gist predicted result
#define USE_AVERAGE 0
#define USE_MAXIMUM 1


#define NUM_DNA 4   // A, C, G, T, number of DNAs
#define NUM_AA 21  // A,R,N... , number of amino acids
#define NUM_BLOSUM 24   //scalar of BLOSUM matrix

#ifndef NUM_NUC
#define NUM_NUC 12
#endif

#ifndef NUM_SHAPE
#define NUM_SHAPE 8
#endif

#define NUM_CHAINS 70   //maximal number of chains per protein

// for resdiue record
#define RESSEQ     0
#define AASEQINDEX 1

// alignment output mode, vertical or horizontal
#ifndef ALIGNMODE
#define ALIGNMODE
#define HORIZONTAL 0
#define VERTICAL 1
#endif

#ifndef SEQALIGN
#define SEQALIGN
//below are 4 tokens defined for dynamic programming
#define STOP 0
#define UP   1         //from up
#define LEFT 2         //from left
#define DIAG 3         //from diagonal 

// for the alignment result
#define IDT 3 // identity
#define SIM 2 // similar, positive
#define MIS 1 // mismatch
#define GAP 0 // gap
#endif

extern int DIGIT_INDEL ; // numeric representation of insertion and deletion
extern int GAP_DIGIT;
extern char CHAR_INDEL; // character representation of insertion and deletion
extern char GAP_CHAR ;
extern char CHAR_VECTOR_ID_SEPRATOR; //the separator for items in the record id of svm vectors


// logical operation
#define AND 1
#define OR  0

// define the constant for GetBegEnd, when doing the binding site analysis
#define FIRST 0
#define LAST  1
#define INTER 2


//#define INIT_CONSV       -99

#define UNKNOWN_AA       'X'  // unknown amino acid
#define UNKNOWN_SHAPE    '-'  // unknown shape symbol
#define UNKNOWN_DSSP_SEC '-'  // unknown dssp secondary structure symbol, for example, the missing residue in the ATOM record compared to SEQRES record 
#define DSSP_SEC_RANDOM  'R'  // symbol for the random coil, in dssp it is recorded as ' '
#define CHAR_DNA         '*'
#define CHAR_NON_RESIDUE 'X' // non residue in pdb seqres record, e.g, DNA or RNA
#define NULL_ICODE       '.'

#define NUM_20_AA         20
#define NUM_STDAA 21 //define the number of standard name of residues appear in PDB file
extern const char *STD3CharAA_alphabet[];
extern const char UNKNOWN_3CharAA[];
extern const char STD1CharAA_alphabet[];
//extern const char *EXT3CharAA_alphabet[] = {"ALX","GLX","UNK"};

extern const char *BLOSUM3D_alphabet[] ;
extern const char BLOSUM1D_alphabet[] ;
extern const char AAAlphabet_Tuping[] ;
extern const int MAPINDEX_BLOSUM_Tuping[] ;
extern const int AAS_Code[];

extern const char *PDB3Char_AllRes [] ;
extern const char PDB1Char_AllRes[] ;

extern const char *PDB3CharAA_alphabet[] ;
extern const char PDB1CharAA_alphabet[] ;

extern const double Hydrophobicity[] ;
extern const double Background_AA_Freq[] ;
extern const double Ln_Background_AA_Freq[] ;
#define NUM_SHAPES 9   //scalar of shag_matrix
extern const char SHAPE_alphabet[] ;

extern int numNonMetalHetGroup ;
extern const char *nonMetalHetGroup[] ;

//encoding for feature vectors
//encoding score1, information content
#define SIZE_ENCODE_PSSM   20
#define SIZE_ENCODE_SCORE1 5
extern int EncodeMatrixScore1[][SIZE_ENCODE_SCORE1] ;
extern double EncodeScaleScore1[] ;

#define SIZE_ENCODE_SCORE2 5
extern int EncodeMatrixScore2[][SIZE_ENCODE_SCORE2] ; 
extern double EncodeScaleScore2[] ;


#define SIZE_ENCODE_AA_CENTER 5 /*this is for residues of interest centered in a window*/
extern char EncodeMatrixAA_Center_alphabet[] ;
extern int EncodeMatrixAA_Center[][SIZE_ENCODE_AA_CENTER] ;

#define SIZE_ENCODE_AA_NON_CENTER 5 /*this is for residues of interest centered in a window*/
extern char EncodeMatrixAA_Non_Center_alphabet[] ;
extern int EncodeMatrixAA_Non_Center[][SIZE_ENCODE_AA_NON_CENTER] ;
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
extern int EncodeMatrixWaterAcc[][SIZE_ENCODE_WATERACC] ;
extern int EncodeScaleWaterAcc[] ;


#define SIZE_ENCODE_DISTANCE 5
extern int EncodeMatrixDistance[][SIZE_ENCODE_DISTANCE] ;
extern int EncodeScaleDistance[] ;

/*#define SIZE_ENCODE_HYDROPHOBICITY 7*//*{{{*/
/*hydrophobicity [0.0, 1.0]*/
/*int EncodeMatrixHydrophobicity[][SIZE_ENCODE_HYDROPHOBICITY] = */
/*{                                                              */
/*    { 3,0,0,0,0,0,0},//0~0.02    R                             */
/*    { 0,3,0,0,0,0,0},//0.02~0.2  D, E, H                       */
/*    { 0,0,3,0,0,0,0},//0.2~0.5   N, Q, K, S, T                 */
/*    { 0,0,0,3,0,0,0},//0.5~0.65  G, A                          */
/*    { 0,0,0,0,3,0,0},//0.65~0.7  C,                            */
/*    { 0,0,0,0,0,3,0},//0.7~0.9   P,M,V,W,Y                     */
/*    { 0,0,0,0,0,0,3} //0.9~+     I,L,F                         */
/*};                                                             */
/*double EncodeScaleHydrophobicity[] =                           */
/*{                                                              */
/*    0.02,                                                      */
/*    0.2,                                                       */
/*    0.5,                                                       */
/*    0.65,                                                      */
/*    0.7,                                                       */
/*    0.9                                                        */
/*};                                                             */
/*}}}*/
#define SIZE_ENCODE_HYDROPHOBICITY 3
/*hydrophobicity [0.0, 1.0]*/
extern int EncodeMatrixHydrophobicity[][SIZE_ENCODE_HYDROPHOBICITY] ;
extern double EncodeScaleHydrophobicity[] ;                          


// define functions
#define GetPatternTransMatrix ReadSMatrix
#define AlignAna_Protein      AlignAna
#define Is0Chain(x) IsSpecChain(x, '0')

namespace ZnBindingProtein
{
    extern int MAX_NUM_ZNPRO ;
    extern int MAX_NUMRES    ;
    extern int MAX_ATOMENV   ;
    extern int MAX_BOUND_SITE ;
    extern int MAX_BOUND_METAL_PER_RES ;
}

namespace MetalBindingProtein
{
    extern int MAX_NUM_METALPRO ;
    extern int MAX_NUMRES       ;
    extern int MAX_ATOMENV      ;
    extern int MAX_BOUND_SITE    ; 
    extern int MAX_BOUND_METAL_PER_RES ;
}

/* to define a struct, there are two methods
 * 1).
 *   typedef struc
 *   {
 *       bla,bla
 *   }SNAME;
 *
 *   this method is valid in both C and C++
 * 2).
 *   struct SNAME2
 *   { 
 *       bla,bla
 *   };
 *   this method is valid only in C++
 */


struct Element //define the structure for chemical element/*{{{*/
{
	int atomNum;   // element atomic number
	char name[SIZE_ATOM_ELEMENT+1] ; // element name
};/*}}}*/
struct AtomFreq/*{{{*/
{
	char atomName[4];
	int atomCnt;
};/*}}}*/
struct PredPro/*{{{*/
{
    char id[SIZE_CHAIN_ID+1];
    int  *resSeqIndex;
    int  *speResSeqIndex;
    int numRes;
    int numSpeRes;
};/*}}}*/
struct Atom //define the structure for ATOM record in PDB file/*{{{*/
{ 
    char   recordID[SIZE_RECORD_ID+1];
	int    serial;        // atom serial number
	char   name[SIZE_ATOM_NAME+1];       // atom name
	char   origName[SIZE_ATOM_ORIGNAME+1];   // record the four columns in ATOM record
	char   altLoc;        // Alternate location indicator
	char   resName[SIZE_RES_NAME+1];    // 3-letter residue name
	char   chainID;       // chainID
	int    resSeq;        // residue sequence number ;
	char   iCode;         // Code for insertion of residues
	double x;             // x coordinate
	double y;             // y coordinate
	double z;             // z coordinate
	float  occupancy;     // occupancy
	float  tempFactor;    // temperature factor
	char   segID[SIZE_ATOM_SEGID+1];      // segment identifier
	char   element[SIZE_ATOM_ELEMENT+1];    // element symbol
	char   charge[SIZE_ATOM_CHARGE+1];     // charge on the atom

    double *dist;         // record the distance of this atom to other atoms
    int    numDist;       // number of distances
};/*}}}*/
struct Residue  //data structure for the residue/*{{{*/
{
	int  resSeq;             //residue sequential number in ATOM record
	char resICode;           //insertion code
	int  aaSeqIndex;         // residue sequence index number in SEQRES
	char resName[SIZE_RES_NAME+1];
	char aa;                 // one letter residue name
	char chainID;
	char shape;
	double  consv;           // conservation level of the residue according to PSSM
	bool isBioMetalBound;    // if the residue is biologic-ally metal bound
    int  numMetalBound;      // number of metals bound to the residue, for co-catalytic zinc, one Asp can bind to two metal atoms
	int *parentAtomEnvIndex; // used only when res is a mmeber of MetalPro, usually the size of parentAtomEnvIndex = 1
                             // but in some co-catalytic sites, one residue can bind to more than one metal atoms
	int  numAtom ;
	Atom *atom;
//	Atom atom[MAXATOM];
};/*}}}*/
struct AtomEnv  // environment of atoms, the residues within /*{{{*/
                // CUTOFF (in Angstroms) distance to an atom
{
	int      numRes;           // number of residues included in the atomEnv, for example including only residues on a single chain
    int      totalBoundRes;    // number of all residues in the protein bind to the metal, 2007-04-19
    int      metalAtomResSeq;  // residue serial number for metalAtom  in PDB, 2007-04-19
	char     metalAtomName[SIZE_METAL_ATOM_NAME+1];         // metal atom element symbol
	char     metalAtomResName[SIZE_METAL_ATOM_RES_NAME+1];  // residue name of metal atom, for example "HEM"
	char     metalAtomChainID;    // chainID annotation for the metal atom in PDB file
    char     metalAtomPDBID[SIZE_PDBID+1];
	char     id[SIZE_CHAIN_ID+1];// chain identifier of nrPDB we are looking for
	int      seqLength;          // length of the seqence of that chain
	Residue *res;
    int     *parentResIndex;    // used only when atomEnv is a member of MetalPro
                                // in that case, for example, if this metal binds to 3 residues, with
                                // parentResIndex=(0,2,4), then residues under pMetalPro->atomEnv[0] are
                                // int *pIndex = pMetalPro->atomEnv[0].parentResIndex
                                // pMetalPro->res[pIndex[0]]
                                // pMetalPro->res[pIndex[1]]
                                // pMetalPro->res[pIndex[2]]
};/*}}}*/
struct DomainDEF/*{{{*/ //// parsered domain definition, e.g. A:1-32,A:35-140 
{
	int  numChain;                             // number of chains, here chain might be the segment of a single chain, for example d1hcz_1 1-167,231-250, there are two segments;
	bool isWholeChain[NUM_CHAIN_PER_DOMAIN];   // whether the SCOP domain including a while chain
	char chainIDs[NUM_CHAIN_PER_DOMAIN+1];
	int  posF[NUM_CHAIN_PER_DOMAIN];           // start position in the sequence
	int  posT[NUM_CHAIN_PER_DOMAIN];           // end position in the sequence
    char icodeF[NUM_CHAIN_PER_DOMAIN+1]; // icode for the start position
    char icodeT[NUM_CHAIN_PER_DOMAIN+1]; // icode for the end position
};
/*}}}*/
struct SCOP/*{{{*/
{
	char      did[SIZE_SCOP_ID+1];                   // domain id, e.g d9icwa1
	char      pdbid[SIZE_PDBID+1];                   // pdbid, e.g 9icw
	char      chainID;                               //chain id, e.g A
	char      domainDef[SIZE_DOMAIN_DEF_RECORD+1];   // string for domain record definition, e.g. A:205-298
	int       idnum;    // domain id number, e.g. 17982
	int       cl;       // class id number
	int       cf;       // fold id number
	int       sf;       // superfamily id number
	int       fa;       // family id number
	int       dm;       // domain id number
	int       sp;       // species id number
	int       px;       //

	char      cls;      // class id, e.g. a, in a.60.6.1
	int       fod;      // fold id, e.g. 60 in a.60.6.1
	int       sup;      // superfamily id, e.g. 6 in a.60.6.1
	int       fam;      // family id, e.g. 1 in a.60.6.1
	DomainDEF domDef;   // parsered domain definition
};
/*}}}*/
struct MSA/*{{{*/  /*data structure for multiple sequence alignment*/
{
    int numSeq;
    char **alnSeqName; //record name for each sequence
    char **alnSeq;   // size= numSeq*seqlength , aligned sequence with gaps
    //char **seq;      // sequence without gaps
    int **gaplessSeqIndex; //size= numSeq *seqLength,  indexing for the alnSeq, only non gap residues are counted, gaps set to DIGIT_INDEL
    int *alnSeqLen;  // size=numSeq;
    int *seqLen; //size = numSeq, without gaps
    float *pid; //percent identity of the shorter sequence, to the keyseqence, key sequence is the sequence with id Query, otherwise the first sequence.
};/*}}}*/

//LOG: 2006-11-06 18:47:23 Monday  Week 45 <nanjiang@casio>
// suggestion, MetalPro can be in the structure
// struct MetalPro
// {
//    int numBoundRes;         // number of bound Residues
//    char id[SIZE_CHAIN_ID+1];
//    int length;
//    int numMetalAtom;
//    char **metalAtomList;   // binding metal atom list of the chain
//    int *resIndex;
//    Residue *res;
//    AtomEnv *atomEnv;// with metalPro->atomEnv->res = metalPro->res[index]
// }
struct MetalPro/*{{{*/
{
    char id[SIZE_CHAIN_ID+1];
    int numBoundRes;        // number of bound Residues
    int length;             // length of the protein sequence
    int numMetalAtom;       // number of metal atoms for the protein
    char **metalAtomList;   // list of all bound metal atoms in the chain, just listing, not unique
    Residue *res;           // residues binding to metals for this protein
	int *resSeqIndex;       // array of aaSeqIndex for residues binding to metals for this protein
	AtomEnv *atomEnv;       // atomEnv for this protein
};/*}}}*/

struct MetalPro2 //  in MetalPro2, residues are stored according to metal atom/*{{{*/
{
    int numMetalAtom;         // number of metal atoms bound to the polypeptide chain
    char id[SIZE_CHAIN_ID+1];
    int length;
    AtomEnv *atomEnv;
};/*}}}*/
struct Chain //data structure for one chain of a protein /*{{{*/
{
    int8  seqtype;    /*sequence type, AA_SEQ: amino acid seq, DNA_SEQ: nucleic acid seq, SHAPE_SEQ: shape strings*/
    char  *title;     /*title of the sequence, e.g. the annotations in the PDB file, 2008-02-06*/
    char  pdbid[SIZE_PDBID+1];
    char  chainID;    // chain ID
    int   numRes;     // residue number
    char *aaSeq;      // amino acid sequence
    int  *resSer;     // residue serial number (in ATOM record) array
    char *resICode;   // reside insertion code array
    char *shString;   // shapeString;
    char *secStruc;   // secondary structure
    int  *waterAcc;   // water accessibility
    float *phi;       // phi torsion angle
    float *psi;       // psi torsion angle
    DATATYPE_CONSV  *consv;      // conservation level
};/*}}}*/
struct GistPredChain/*{{{*/
{
    char    id[SIZE_CHAIN_ID+1];
    int     numRes;     // number of residues in GistPredChain, might not be equal to the length of the chain
    int     length;     // total number of residues of pdb chain
    char   *aaSeq;
    int    *aaSeqIndex;
    int    *label;
    double *discriminant;// discriminant for each residue
};/*}}}*/
struct DSSP_HBond/*{{{*/
{
	int    pos;   // positon of hydrogen bond
	double e;     // electric energy
};/*}}}*/
struct DSSP_Residue/*{{{*/
{
	int        seqResSer;    // sequential residue number
	int        resSeq;       // resSeq in PDB File
	char       iCode;        // iCode in PDB File
	char       chainID;
	char       aa;           // 1 letter amini acid name
	char       chainBreak;
	char       ss;           // secondary structure lable
	char       turn3;        // 3 turn helix
	char       turn4;        // 4 turn helix
	char       turn5;        // 5 turn helix
	char       geoBend;      // geometrical bend
	char       chira;        // chirality
	char       bBridge1;     // beta bridge lable
	char       bBridge2;     // beta bridge lable
	int        bp1;          // beta bridge partner resnum
	int        bp2;          // beta bridge partner resnum
	char       bSheet;       // beta sheet lable
	int        acc;          // solvent accessibility,ACC
	DSSP_HBond hbond[4];     // h-bond
	float      tco;          // cosine of angle between C=O of residue I and C=O of residue I-1
	float      kappa;        // virtual bond angle (bend angle) defined by the three C-alpha atoms of residues I-2,I,I+2.
	float      alpha;        // virtual torsion angle (dihedral angle) defined by the four C-alpha atoms of residues I-1,I,I+1,I+2
	float      phi;          // IUPAC peptide backbone torsion angles
	float      psi;          // IUPAC peptide backbone torsion angles
	float      x;            // Ca x coordinate
	float      y;            // Ca y coordinate
	float      z;            // Ca z coordinate
};/*}}}*/
struct SSBond/*{{{*/
{
    bool    isInterChain;   // interChain or intraChain
    int     numRes;         // in case of interchain ssbond, sometimes it will be 1
    Residue res[2];
};/*}}}*/
struct SSBondPro/*{{{*/
{
    char id[SIZE_CHAIN_ID+1];
    char length;
    int numSSBond;
    int numSSBondRes;
    SSBond *ssbond;
};/*}}}*/
struct AlignFactor/*{{{*/
{
	float score;
	float zScore;
	float pozScore;
	double eValue;
	int   idt_cnt;
	float identity;
	float identity_short;
	int   sim_cnt;
	float similarity;
	float similarity_short;
	int   gap_cnt;
	float gapPercent;
};/*}}}*/
struct MetalBindPro/*{{{*/
{
    char id[SIZE_CHAIN_ID+1];      //nrPDBb chain identifier
    bool isBound;
    DATATYPE_CONSV  maxConsv;   //maximal conservation for keyAA
};/*}}}*/
struct MODM/*{{{*/
{
    char     id[SIZE_CHAIN_ID+1];
    int      length;
    int      type_modm;
    DATATYPE_MODM_MATRIX    **M;     //the pssm matrix
    DATATYPE_MODM_MATRIX    **log_M; //log(Mij), this is used when profile-profile alignment is needed, the log(Mij) is precalculated to save the computational time, in that case, Mij should be weighted percentage
    double  *score1;
    double  *score2;
    double  *consv;  //conservation level for certain residues, say CHDE
    char    *aaSeq;
    char    *alphabetMODM;
    int     *waterAcc;
    char    *shString;
    char    *dsspSec;
};/*}}}*/
struct Profile /*{{{*/
{ /*profile, for each line of the pssm matrix*/
    int   aaSeqIndex;
    char  aa;
    float p[20];
    float score1;
    float score2;
};/*}}}*/
struct ProfileByte /*{{{*/
{ /*profile, for each line of the pssm matrix*/
    short   aaSeqIndex;
    char  aa;
    int8 p[20];
    float score1;
    float score2;
};/*}}}*/
struct ProfileSAD/*{{{*/
{ /*profile, for each line of the pssm matrix, including the Shape-Acc-Dsspsec information*/
    int   aaSeqIndex;
    char  aa;
    char  shape;
    int   waterAcc;
    char  dsspSec;
    float p[20];
    float score1;
    float score2;
};/*}}}*/
struct ProfileSADByte/*{{{*/
{ /*profile, for each line of the pssm matrix, including the Shape-Acc-Dsspsec information*/
    short   aaSeqIndex;
    char  aa;
    char  shape;
    int8   waterAcc;
    char  dsspSec;
    int8 p[20];
    float score1;
    float score2;
};/*}}}*/
/*structure for reading and writing binary Frag files generated by search_new
 * program 2009-06-23, Nanjiang*/
/*format of the binary frag file
 * short type (0 for six columns and 1 for five columns)
 * idList: the chain idlist sorted in accending order
 * posTar numCanFrag Nx|idxInner idxOuter posCan score| or
 * posTar numCanFrag Nx|idxInner posCan score|*/
struct FragCanShort5/*{{{*/
{
    short idxInner; /*id index within the frag file.*/
    short posCan;   /*position of the candidate fragment in sequence*/
    float score;    /*score of the candidate fragment to the target fragment*/
};/*}}}*/
struct FragCanShort6/*{{{*/
{
    short idxInner; /*id index within the frag file.*/
    short idxOuter; /*id index for the external supplied id list*/
    short posCan;   /*position of the candidate fragment in sequence*/
    float score;    /*score of the candidate fragment to the target fragment*/
};/*}}}*/
struct FragCanInt5/*{{{*/
{
    int idxInner; /*id index within the frag file.*/
    int posCan;   /*position of the candidate fragment in sequence*/
    float score;    /*score of the candidate fragment to the target fragment*/
};/*}}}*/
struct FragCanInt6/*{{{*/
{
    int idxInner; /*id index within the frag file.*/
    int idxOuter; /*id index for the external supplied id list*/
    int posCan;   /*position of the candidate fragment in sequence*/
    float score;    /*score of the candidate fragment to the target fragment*/
};/*}}}*/

struct dbindex{/*{{{*/
    int dbfileindex;
    long offset;
    unsigned long size;
};/*}}}*/

//#define MAX2(x,y) ((x) > (y) ? (x) : (y))
//#define MIN2(x,y) ((x) < (y) ? (x) : (y))

/*functions for reading in structured file*/
/* functions for reading PDB ATOM record line*/
void ScanfAtomSerial(const char *PDBAtomRecordLine, int &atomSerial);
void ScanfAtomResName(const char *PDBAtomRecordLine, char * atomResName);
void ScanfAtomOrigName(const char *PDBAtomRecordLine, char * atomOrigName);

void ScanfCoorRecord_Atom(const char *line, Atom *pAtom);
void ScanfCoorRecord_Atom_Simp1(const char *line, Atom *pAtom);
int  Scanf_SEQRES_Record(const char* line, int resBegin, char* title, int* serNum,char* chainID,int *numRes,char **resName);
void Scanf_SEQRES_Para(const char* line,char* title, int& serNum,char& chainID,int& numRes, char* resName);
int  Scanf_SEQRES_Seq(const char* line, int resBegin, char **resName);

void ScanfSSBondRecord(FILE* fp, SSBondPro* pSSBondPro, int numSSBond);
void ScanfDSSPResRecord(const char* line, DSSP_Residue* pDSSPRes);

int  ScanfGistVectorLabelCard(const char *line, char *id, int &length, int &idx, char *res_1char_list,int *aaSeqIndex, int &label,  int numSite);
int  ScanfGistPredictCard(const char *line, char *id, int &length, int &idx, char *res_1char_list,int *aaSeqIndex, int &label, double &discriminant, int numSite);
int  ScanfKernerlRecordID(const char *vectorRecordID, char *id, int &length, int &idx, char *res_1char_list, int *aaSeqIndex, int numSite);
void ScanfResRecord3(const char *str, Residue *pRes, int tag);
int  ScanfCloseMetalRes(FILE* fpNrCloseMetal, Residue* res,int numRes);

void ScanfSCOPRecord(const char *line, SCOP* pSCOP);

/*get data from file*/
int GetAtom_PDB(const char *pdbfile, Atom *atom, int *pNumAtom, bool isGetAllChain = true, int max_num_atom = MAX_ATOM_SERIAL, char* chainIDList = "", bool isRestrictByContactAtom = true, Atom *metalAtom = NULL, int *pNumMetalAtom = NULL, int max_num_metal_atom  = MAX_METAL_CHAIN);
int GetSeq_ATOM(const char* pdbfilepath, const char *id, char *aaSeq, int *resSer, char *resICode, FILE *fpLog = stdout);
int GetSeq_ATOM(const char* pdbfile, Chain *chain, char* chainIDList,  bool isGetAllChain = false, bool isIncludeHETATM = true, int max_seq_length = LONGEST_SEQ, FILE *fpLog = stdout);
int GetDSSPChain(const char* id, Chain* pChain, const char* dsspfilepath);
int GetSEQMAP(const char* seqmapfilepath, Chain *pChain);
int GetMetalElementList(Element *metalEle, const char metalElementListFile[] = "");
int GetBkFreq(const char *infile, double *P, const char *alphabet, int num_aa = NUM_AA);

int8 GetBinaryMODMPara(const char *file, int &sizeAlphabet, int &length, int8 &typeProfile);
template <class T> int GetBinaryMODM(const char *file, char *alphabet, int &length, T *profile, double *parameter, int8 &typeProfile);
template <class T> int GetBinaryMODM_fp(FILE *fpin, long readsize, char *alphabet, int &length, T *profile, double *parameter, int8 &typeProfile);

int GetBinaryFragShort5(const char *outfile, int &fragFileType, char **idList, int &numID, int &maxSizeID, int &length, short *posTar, short *numCan, int &totalFragCan, FragCanShort5 *fragCan);
int GetBinaryFragShort6(const char *outfile, int &fragFileType, char **idList, int &numID, int &maxSizeID, int &length, short *posTar, short *numCan, int &totalFragCan, FragCanShort6 *fragCan);
int GetBinaryFragInt5(const char *outfile, int &fragFileType, char **idList, int &numID, int &maxSizeID, int &length, int *posTar, int *numCan, int &totalFragCan, FragCanInt5 *fragCan);
int GetBinaryFragInt6(const char *outfile, int &fragFileType, char **idList, int &numID, int &maxSizeID, int &length, int *posTar, int *numCan, int &totalFragCan, FragCanInt6 *fragCan);
int ReadMSA(const char* msafile, MSA *pMSA, int format = CLUSTALW);
int ReadMSA_clustalw(const char* msafile, MSA *pMSA);

void GetDBFPList( vector <FILE*> &fpList, string dbname, int maxdbfileindex);
int ReadDatabaseIndex(string dbname, map<string, dbindex> &dbindexmap,int &maxdbfileindex );
int ReadInDatabase(const char *id, const char *frag_acc_path, const char* qijpath, const char* modmpath, int qijformat, int modmformat, int fragaccformat, int **matFrag, int **matQij, int **matMODM, int **matMerged, float *matFragScore1, float *matFragScore2,  float *matQijScore1, float *matQijScore2, float *matMODMScore1, float *matMODMScore2, char *aaSeq, int *digitAASeq, char *shapeSeq, char *dsspSecSeq, int type_dataset, bool isReadBinaryFile, int NPer_Frag_Database, int ratioScheme, int mergeSide, int dsspMapMethod = 0, int8 typeProfile =1 );
int ReadInDatabase_dumpedfile(const char *id, map<string,dbindex>&dbindexfragacc, map<string,dbindex>&dbindexqij, map<string,dbindex>&dbindexmodm, vector <FILE*> &fpList_fragacc, vector<FILE*>&fpList_qij, vector<FILE*>&fpList_modm, int **matFrag, int **matQij, int **matMODM, int **matMerged, float *matFragScore1, float *matFragScore2,  float *matQijScore1, float *matQijScore2, float *matMODMScore1, float *matMODMScore2, char *aaSeq, int *digitAASeq, char *shapeSeq, char *dsspSecSeq, int type_dataset, bool isReadBinaryFile, int NPer_Frag_Database, int ratioScheme, int mergeSide, int dsspMapMethod = 0, int8 typeProfile = 1);

//int GetPatternTransMatrix(const char* patterTransFile, double** M, char* alphabet) ;
//int ReadSMatrix(const char *filename, DATATYPE_SMATRIX **S, char *alphabet);
template <class T> int ReadSMatrix(const char *filename, T **S, char *alphabet);
int ReadSeq_FASTA(const char *fileName, char* seq, int *pSeq_type = NULL, int maxlength = LONGEST_SEQ, char* annotationLine = NULL, int maxSizeAnnotationLine = 50);
int ReadNextSeq_FASTA(FILE *fp, char* seq, int *pSeq_type = NULL, int maxlength = LONGEST_SEQ,char* annotationLine = NULL, int maxSizeAnnotationLine = 50);

template <class T> int GetMODM(const char* modmfilepath, T **M = NULL, char* alphabetMODM = NULL,char* aaSeq = NULL, double* score1 = NULL,double* score2 = NULL, double* parameter =NULL, int *seqIndex = NULL, int startIndex = 1);
int GetTupingQij(const char *infile, int **M = NULL, char *aaSeq = NULL, char *shapeSeq = NULL, char *dsspSeq = NULL, int *waterAcc = NULL, double *score1 = NULL, double *score2 = NULL, int *seqIndex = NULL, int startIndex = 1, const int *mapIndex=MAPINDEX_BLOSUM_Tuping);
int GetPSSM(const char *infile, int &length, char *aaSeq = NULL, int **Mij = NULL, int **fij = NULL, double *score1 = NULL, double *score2 = NULL, double *parameter = NULL, int *seqIndex = NULL);
template <class T> void CalLogM(T **M, T **log_M, int xSize, int ySize,  double roundoff_scale  = 1.0);
int ReadInProfile(const char *file, int **M, char *aaSeq = NULL, char *shapeSeq = NULL, int *waterAcc = NULL, char *dsspSec = NULL, float *score1 = NULL, float *score2 = NULL, int *sumProfile = NULL, int maxLength = LONGEST_SEQ);
int ReadInProfile_fp(FILE* fpin, unsigned long readsize, int **M, char *aaSeq = NULL, char *shapeSeq = NULL, int *waterAcc = NULL, char *dsspSec = NULL, float *score1 = NULL, float *score2 = NULL, int *sumProfile = NULL, int maxLength = LONGEST_SEQ) ;

int GetAltAtomIndex(Atom *atoms, int n);
int GetAASeqIndex(const char* id, int resSeq, char resICode, Chain* pChain, int startIndex = 0);
int GetResidueIndex(Atom *pAtom, Residue *resGroup, int numRes);
int GetResidueIndex(Residue *pRes, Residue *resGroup, int numRes);
template <class T> T  GetIntegConsv(char aa, T* V,  char* alphabetMODM, int type_modm = MODM_PER);

double GetContactDist(const char* eleName);

int GetZnBoundRes(const char* metalProFile, MetalPro* znPro, int min_numRes = 3 , int max_numRes = 4, bool isUsingTotalBoundRes = true);
int GetZnBoundRes(const char* metalProFile, MetalPro2* znPro, int min_numRes = 3, int max_numRes = 4, bool isUsingTotalBoundRes = true);
int GetZnBoundRes3(const char* metalProFile, MetalPro* znPro, MetalPro* znPro1, MetalPro* znPro2, int &numZnPro, int &numZnPro1, int &numZnPro2, int &numZnRes, int &numZnRes1, int &numZnRes2, double cutoff_score2, char *resList, const char* modmpath, int min_ZnBoundRes = 3, int max_ZnBoundRes = 4, bool isUsingTotalBoundRes = true);
int GetZnBoundRes3(const char* metalProFile, MetalPro2* znPro, MetalPro2* znPro1, MetalPro2* znPro2, int &numZnPro, int &numZnPro1, int &numZnPro2,int &numZnRes, int &numZnRes1, int &numZnRes2, double cutoff_score2, char *resList, const char* modmpath, int min_ZnBoundRes = 3, int max_ZnBoundRes = 4, bool isUsingTotalBoundRes = true);
int GetSSBondRes(const char* ssbondProFile, SSBondPro *ssbondPro, int &numSSBondPro, int &numSSBondRes);
int GetSSBondRes2(const char *ssbondProFile, SSBondPro *ssbondPro, SSBondPro *ssbondPro1,int &numSSBondPro, int &numSSBondPro1, int &numSSBondRes, int &numSSBondRes1, double cutoff_score2, const char* modmpath);

int GetSVMVector(int aaSeqIndex, MODM *pMODM, int K, int P,double* vector, double cutoff_score2 = 0.1);
int GetSVMVector(int aaSeqIndex1, int aaSeqIndex2, MODM *pMODM, int K, int W,int P, double *vector, double cutoff_score2 = 0.1);
int GetSVMVector(int *aaSeqIndex, int numSite, MODM *pMODM, int K, int W, int P, double *vector, int vectorDim, double cutoff_score2 = 0.1, int encoding_type = PSSM_PROFILE_ENCODE);
int GetDimensionVectorPerSite(int P, bool isCenterResidue);
int GetSVMVector_2(int *aaSeqIndex, int numSite, MODM *pMODM, int K, int W, int P, double *vector, int vectorDim, double cutoff_score2 = 0.1, int encoding_type = PSSM_PROFILE_ENCODE);
int GetSVMVectorPerSite(double *vector, int vector_beg, int seqPos, int P, MODM *pMODM, double W1, double W2, double W3, double W4, double W5, double cutoff_score2, int encoding_type  = PSSM_PROFILE_ENCODE);
int GetSVMVectorPerSite_2(double *vector, int vector_beg, int seqPos, int P, int dimVectorPerSite, MODM *pMODM, double w1, double w2, double w3, double w4, double w5, double cutoff_score2, bool isCenterResidue, int encoding_type = PSSM_PROFILE_ENCODE);

int RestrictMetalBoundRes(MetalPro *metalPro1, int &numMetalPro1, int &numMetalRes1, MetalPro* metalPro2, int numMetalPro2, char *resList, int level, const char *modmpath = "", double cutoff_score2 = 0.0);
int GetMetalBoundRes3(const char* metalProFile, MetalPro* metalPro, MetalPro* metalPro1, MetalPro* metalPro2, int &numMetalPro, int &numMetalPro1, int &numMetalPro2, int &numMetalRes, int &numMetalRes1, int &numMetalRes2, double cutoff_score2, char *resList, const char* modmpath, bool isExcludeOtherChainRes = true, char **keyMetalList = NULL, int numKeyMetal = 0, int min_metalBoundRes = 0 , int max_metalBoundRes = 0, bool isUsingTotalBoundRes  = true);
int GetMetalBoundRes2(const char* metalProFile, MetalPro* metalPro1, MetalPro* metalPro2, int &numMetalPro1, int &numMetalPro2, int &numMetalRes1, int &numMetalRes2, char *resList, const char* modmpath, bool isExcludeOtherChainRes = true, char **keyMetalList = NULL, int numKeyMetal = 0, int min_metalBoundRes = 0 , int max_metalBoundRes = 0, bool isUsingTotalBoundRes  = true);
int GetMetalBoundRes(const char* metalProFile, MetalPro* metalPro, int &numMetalPro, bool isExcludeOtherChainRes = true, char **keyMetalList = NULL, int numKeyMetal = 0, int min_numRes = 0 , int max_numRes = 0, bool isUsingTotalBoundRes = true);
int GetNumBoundRes(set<int> const& PrP, set<int> const& PrN, MetalPro* pZnPro, set <int> &TP, set <int> &FP, set<int> &TN, set<int> &FN) ;
int GetNumBoundRes(set<int> const& idx, SSBondPro* pPro, set <int> &idxTrue, set <int> &idxFalse);
int GetNumBoundRes(set<int> const& idx, MetalPro* pZnPro, set <int> &idxTrue, set <int> &idxFalse);
int GetNumBoundRes(set<int> const& idx, MetalPro2* pZnPro, set <int> &idxTrue, set<int> &idxFalse);
int GetHCRes(MODM *pMODM, bool *isPolyHis, Residue* HCRes, double cutoff_score2, double cutoff_consv, char* resList, bool isUseConsvI, bool isMaskPolyHis);

int GetNextResidue_PDB(Residue *pRes, FILE *fpPDBFile, bool isGetAllAtomLocation = false);
int SelectAltLocAtom(Atom* atom, int &cntTotalAtom,FILE* fpPDBFile, char **metalEleList = NULL, int numMetalEle = 0, bool isSelectMetal = false);

void GetBegEnd(int aaSeqIndex1, int aaSeqIndex2, int W, int length, int loc,int& beg, int& end);
void GetBegEnd(int aaSeqIndex1, int aaSeqIndex2, int aaSeqIndex3, int W, int& beg, int& end);

/*int GetEncodeAA(char aa);*/
int GetEncodeScore1(double score1);
int GetEncodeScore2(double score2);
int GetEncodeWaterAcc(int waterAcc);
int GetEncodeDistance(int dist);
int GetEncodeHydrophobicity(double hydrophobicity);

/* boolean functions*/
bool IsNonMetalHetGroup(const char *resName);
bool IsDNASeq(const char *aaSeq, int n = 0);
bool IsMetalAtom(Atom* pAtom,char** metalEleList,int numOfMetalEle);
bool IsInContactAtomList(const char *origName);
bool IsSSBonded(int aaSeqIndex, SSBondPro *ssbondPro);
bool IsZnBound(int aaSeqIndex, MetalPro *znPro);
bool IsZnBound(int aaSeqIndex1, int aaSeqIndex2, MetalPro *znPro, bool operation);
bool IsZnBound(int aaSeqIndex1, int aaSeqIndex2, MetalPro2 *znPro, bool operation);
bool IsNewResidue(Atom* pAtom, Residue *resGroup, int numRes);
bool IsNewResidue(Residue *pRes, Residue *resGroup, int numRes);
bool IsHEMPro(MODM *pMODM, bool *isPolyHis);
bool IsSatisfyHCResRule(char* id, MODM *pMODM,  bool *isPolyHis,  char* resList, double cutoff_consv, double cutoff_score2, int window, int min_numHCRes, Residue* HCRes, int& numHCRes, Residue* LCRes, int& numLCRes, bool isMaskPolyHis, bool isUseConsvI);
int IsSpecChain(const char* id, const char chainid = '0');

/*converter*/
string int2string(int number);
void StdRegExp(char* pattern);
char AA3To1(const char* aa3, int mode = 0);
const char* AA1To3(char ch);
void CharToDigit_Protein(const char* str, int8* a,int n);
void DigitToChar_Protein(const int8 *a, int n, char* str);
int  DigitShape(char shape);
int  DigitAA(char *resName);
char* StdID(char* id);  // standardize chain identifier, id will be changed
char* StdElementName(const char* elementName, char *stdElementName);//elementName will not be changed
char* Stdid2Seqresid(char* id);
void RemoveTID(char* id); // remove the trailing white space of chain identifier
char *StdPDBFileName2PDBID(const char* rtname_pdbfilename, char *pdbid);
char DSSP8to3(char dsspSec, int method = 0);


/*structure implementation, init, copy, delete, allocate*/
void InitChain(Chain *pChain);
void CopyChain(Chain *from, Chain *to);
void DeleteChain(Chain *pChain);

void InitGistPredChain(GistPredChain *pChain);
void AllocGistPredChain(GistPredChain *pChain, int size);
void CopyGistPredChain(GistPredChain *to, GistPredChain *from);
void DeleteGistPredChain(GistPredChain *pChain);

void InitMODM(MODM *pMODM);
void AllocMODM(MODM *pMODM, int length, bool isAllocLogM = false, int sizeAlphabet = NUM_BLOSUM);
void CopyMODM(MODM *to, MODM *from);
void DeleteMODM(MODM *pMODM, int length);

void InitDomainDEF(DomainDEF *pDomainDEF);
void InitSCOP(SCOP *pSCOP);
void CopySCOP(SCOP* to, SCOP* from);

void InitResidue(Residue *pRes);
void CopyResidue(Residue *pRes1, Residue *pRes2);
void DeleteResidue(Residue *pRes, int numAtom = 0) ;

void InitAtom(Atom * pAtom);
void CopyAtom(Atom *pAtom1, Atom *pAtom2);
void DeleteAtom(Atom *pAtom);

void InitDSSPRes(DSSP_Residue* dssp_res);

void InitAtomEnv(AtomEnv *pAtomEnv);
void CopyAtomEnv(AtomEnv *pAtomEnv1, AtomEnv * pAtomEnv2, bool isCopyResidue = true);
void DeleteAtomEnv(AtomEnv *pAtomEnv, int numRes = 0, int numAtom_Res = 0 );

void InitSSBond(SSBond* pSSBond);
void CopySSBond(SSBond *to, SSBond *from, bool isCopyResidue = true);

void InitSSBondPro(SSBondPro *pSSBondPro);
void CopySSBondPro(SSBondPro *to, SSBondPro *from, bool isCopySSBond = true);
void DeleteSSBondPro(SSBondPro *pSSBondPro);

void InitMetalPro(MetalPro   *pMetalPro);
void CopyMetalPro(MetalPro *pMetalPro1 , MetalPro *pMetalPro2 , bool isCopyAtomEnv = true);
void DeleteMetalPro(MetalPro *pMetalPro, int numMetalAtom = 0, int numBoundRes = 0);

void InitMetalPro(MetalPro2  *pMetalPro);
void CopyMetalPro(MetalPro2 *pMetalPro1, MetalPro2 *pMetalPro2, bool isCopyAtomEnv = true);
void DeleteMetalPro(MetalPro2 *pMetalPro, int numAtomEnv = 0);

void InitPredPro(PredPro *pPredPro);
void CopyPredPro(PredPro *pPredPro1, PredPro *pPredPro2);
void DeletePredPro(PredPro *pPredPro);

void InitAlignFactor(AlignFactor *pAlignFactor);

void InitAtomFreq(AtomFreq *pAtomFreq);
void CopyAtomFreq(AtomFreq *to, AtomFreq *from);

void InitMSA(MSA *pMSA);
void DeleteMSA(MSA *pMSA);
void f_neglect_clustalw_header(FILE *fp);

// write formatted text functions
void WriteSubMatrix(FILE* fpSubMatrix, double** subM, int dim,  char* alphabet,int *cnt, char formatValue[] = "%4.0lf", char formatFreq[] = "%7d");
void WriteResRecord3(char *str, Residue *pRes, int tag);
void WriteCloseMetalRes(FILE* fpCloseMetal, Residue* res, int numRes);
void WriteMetalEnvRes(FILE* fpMetalEnvShape, AtomEnv* pAtomEnv, bool isPrintDebugInfo = false);
void WriteVectorRecordID(char *vectorRecordID, char *rmtID, int length, int idx, char *res_1char_list,int* aaSeqIndex, int numSite);
void WriteCHDEComment(FILE *fpout);
void WriteCHDE(FILE *fpout, char *id, int *speResIndex, int numSpeRes, int *resIndex, int numRes, MODM *pMODM, bool isWriteHeader = true, bool isWriteCounter = true);
void WritePredictStat(FILE* fp,  int ReP, int ReN, int PrP, int PrN, int TP, int FP, int TN, int FN, double Sn, double Sp,int ReP1,int ReN1, int ReP2,int  ReN2, double cutoff_consv, int cutoff_window, char* resList, bool isWriteHeader = false);
void WriteMODMTitle(const char* alphabet, FILE* fpout = stdout );
void WriteMODMProfile(int index, char aa, int* V, double score1, double score2,const char* alphabet, FILE *fpout = stdout);
void WriteFastaSeq(char *str, FILE *fpout = stdout, int beg = 0, int end = 0x7FFFFFFF, int linelength = 70);
void WritePDBAtomRecord(Atom *pAtom, FILE *fpout = stdout);

template <class T> int WriteBinaryMODM(const char *outfile, const char *alphabet, int length, T * profile, double *parameter, int8 typeProfile);

int WriteBinaryFragShort5(const char *outfile, int fragFileType, char **idList, int numID, int maxSizeID, int length, short *posTar, short *numCan, int totalFragCan, FragCanShort5 *fragCan);
int WriteBinaryFragShort6(const char *outfile, int fragFileType, char **idList, int numID, int maxSizeID, int length, short *posTar, short *numCan, int totalFragCan, FragCanShort6 *fragCan);
int WriteBinaryFragInt5(const char *outfile, int fragFileType, char **idList, int numID, int maxSizeID, int length, int *posTar, int *numCan, int totalFragCan, FragCanInt5 *fragCan);
int WriteBinaryFragInt6(const char *outfile, int fragFileType, char **idList, int numID, int maxSizeID,  int length, int *posTar, int *numCan, int totalFragCan, FragCanInt6 *fragCan);

// programs for Zn protein detection
double TransMatrixWeight(Residue * pRes1, Residue * pRes2);
int MaskPolyHis(const char *aaSeq, bool *isPolyHis, int length = 0);
void FilterSSBondPred(const char* gistPredictFile_ss, const char* gistPredictFile_zn, const char* gistPredictFile_zn_2, int numSite_ss, int numSite_zn, double cutoff_discriminant_ss);
void AtomFrequencyAna_AtomEnv(AtomEnv * atomEnvs, int numOfAtomEnv, Element *metalEle, int numOfMetalEle, char * str);
double TranScore(double ***M,Residue* res,char *alphabet, int numRes);
void MergeGistResidue(GistPredChain *pChain, int *idx, int start, int end,  int &label, double &discriminant, int operation = USE_AVERAGE);
void SortUniqueResidue(GistPredChain *pChain, int operation = USE_AVERAGE);
int GetGistPredResidue(const char* gistPredFile, int numSite, GistPredChain *chains, int operation = USE_AVERAGE);
int GetGistPredResidue(FILE* fpGistPred, int numSite, GistPredChain *pChain, int operation = USE_AVERAGE);

int UpdateHCRes(Residue *HCRes, int &numRes, int *idx, int numSHCRRes);
int HCRes_window(Residue* HCRes, int& numRes, int sizeGroup, int cutoff_window);
int SHCR_window(int* HCResAASeqIndex, int numRes, int numSite, int cutoff_window, int cutoff_window_pair, int* shcr_res_idx, int &numSHCRRes);
bool IsSatisfySHCRRule(int *HCResAASeqIndex, int numRes, int min_numHCRes, int cutoff_window_pair, int cutoff_window_zn3, int cutoff_window_zn4, int* shcr_res_idx, int &numSHCRRes);


// alignment related functions
template <class T> double GetScore_PPAlignment(T *v1, T *v2, T *log_v1, T *log_v2, int formula = 0 , const double *P = NULL, const double *log_P = NULL);
template <class T> void AlignAna(const int8* alignX, const int8*alignY, int length, T **subMatr, AlignFactor *pAlignFactor, int* alignRel);
template <class T> void AlignAna_Profile(const int8* alignX, const int8*alignY, int length, T **M1, T **M2, T **log_M1, T **log_M2, AlignFactor *pAlignFactor, int* alignRel);
template <class T> int Alignment(char *Xstr, char *Ystr, char *alphabet, int m, int n,char *title1,char* title2, char* alignXstr0, char* alignYstr0, int* alignRel0,	AlignFactor *pAlignFactor,T gapOpen, T gapExt, T **SM, int seqtype = AA_SEQ);
template <class T> int Alignment_Profile(char *Xstr, char *Ystr, T **M1, T **M2, T **log_M1, T **log_M2, char *alphabet, int m, int n,char *title1,char* title2, char* alignXstr0, char* alignYstr0, int* alignRel0,	AlignFactor *pAlignFactor,float gapOpen, float gapExt,int seqtype = AA_SEQ);
int  SeqMatchAlign_Protein(const char *Xstr, const char* Ystr,char *title1,char* title2, char* alignXstr0, char* alignYstr0, int* alginRel0,FILE *fpout = stdout, bool isPrintToScreen = true);
void WriteAlignmentHeader(float gapOpen, float gapExt, AlignFactor *pAlignFactor, char *title1, char *title2, int m, int n, FILE *fpout);
void WriteAlignment(const char* title1,const char* title2,const char* aXstr,const char* aYstr, int* aRel, int length, int lineLength,FILE* fpout, int outmode = HORIZONTAL );

int GetAllSubset(int **subset, int n, int k);

/*functions for scop */
int FindSCOP(char* pdbid, char chainID, int seqF, int seqT, SCOP* pSCOP, FILE *fpSCOPspi, FILE *fpSCOPspd, char** pdbIDs, long* offsetspi, int numPDB);
int FindSCOP_scopid(char* scopid, SCOP* pSCOP, FILE *fpSCOPspi, FILE *fpSCOPspd, char** pdbIDs, long* offsetspi, int numPDB);

/*misc*/
double GetMergeRatio1(int aaSeqIndex, float *p1, float *p2, int length);
double GetMergeRatio2(int aaSeqIndex, float *p1, float *p2, int length);
double GetMergeRatio3(int aaSeqIndex, float *p1, float *p2, int length);
char *GetQijFilePath(const char *rmtID, char *filename, const char *path, int qijformat = QIJ_FORMAT_NANJIANG, bool isReadBinaryFile = true);
char *GetFragMatFilePath(const char *rmtID, char *filename, const char *path, int fragaccformat = FRAGACC_FORMAT_NANJIANG, bool isReadBinaryFile = true);
char *GetMODMFilePath_Tuping(const char *rmtID, char *filename, const char *path, int modmformat = MODM_FORMAT_NANJIANG, bool isReadBinaryFile = true);
void ReorderMatrix(int **M, const char *alphabet_from, const char *alphabet_to);
template <class T> void TreatAllZeroFij(int **fij, int length, T *score1, char *aaSeq, const char* alphabet);

#endif /*HAS_MYPRO_H*/
