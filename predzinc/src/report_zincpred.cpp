/*
 * =====================================================================================
 *       Filename:  report_zincpred.cpp
 *    Description:  reprot the result of zinc-binding site prediction
 *        Version:  1.0
 *        Created:  10/16/2007 11:53:02 AM CEST
 *       Revision:  none
 *       Compiler:  g++
 *         Author:  Nanjiang Shu (Shu), nanjiang@struc.su.se
 *        Company:  Structural Chemistry, Stockholm Univesity
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <string.h>
#include <map>
#include "array.h"
#include "mytemplate.h"
#include "myfunc.h"
#include "mypro.h"
/*ChangeLog {{{*/
/* ChangeLog 2011-10-10
 *      HTML output updated
 * ChangeLog 2011-10-11 
 *      Support multiple input files
 *}}}*/
double cutoff_znscore = 0.45;
bool isOutputHTML = true;
bool ishtml = true;
int lineLength=70;
char predzinc_version[100]="PredZinc version 1.3";
char layout_css[]=/*{{{*/"\n\
                  tr {font-family:Arial; font-size:14px; text-align:center;vertical-align:middle}\n\
                  td {font-family: Vernada; font-size: 12px; text-align:center;\n\
                      vertical-align:middle}\n\
                      td.aaseq {max-width:180px;word-wrap:break-word; font-family: Sanserif, Sanserif, mono; font-size: 11px; text-align: left; vertical-align:middle}\n\
                      tr.subtable {font-family:Arial; font-size:12px; text-align:center;vertical-align:top}\n\
                      td.subtable {font-family:Arial; font-size:12px; vertical-align:top}\n\
                      p.aaseq {\n\
                          font:11px/16px Sans-serif, mono;\n\
                              margin:0px 0px 0px 0px;\n\
                              padding:0px 0px 0px 0px;\n\
                      }\n\
                      body {\n\
                          margin:0px 0px 0px 0px;\n\
                              padding:10px 15px 10px 15px;\n\
                              font-family:verdana, arial, helvetica, sans-serif;\n\
                              color:#333;\n\
                              background-color:white;\n\
                      }\n\
                      h1 {\n\
                          margin:0px 0px 15px 0px;\n\
                              padding:0px;\n\
                              font-size:28px;\n\
                              line-height:28px;\n\
                              font-weight:900;\n\
                              color:#642c2c;\n\
                      }\n\
                      b.zinc{\n\
                          color:red;\n\
                              font-size:100%;\n\
                              font-family:monospace;\n\
                      }\n\
                      b.index{\n\
                          color:#8B008B;\n\
                              font-size:120%;\n\
                              font-family:monospace;\n\
                      }\n\
                      h2 {\n\
                          margin:6px 0px 12px 0px;\n\
                              padding:0px;\n\
                              font-size:14px;\n\
                              line-height:28px;\n\
                              font-weight:900;\n\
                              color:#000;\n\
                      }\n\
                      h2.headline {\n\
                          margin:12px 0px 6px 0px;\n\
                              padding:0px;\n\
                              background-color:white;\n\
                              color:#4767A6;\n\
                              font-size:20px;\n\
                              font-family:verdana,arial;\n\
                              line-height:200%;\n\
                              border: 0px solid #aaa;\n\
                      }\n\
                      h3.subtitle {\n\
                          margin:22px 0px 12px 0px;\n\
                              padding:0px;\n\
                              color:#4B0082;\n\
                              font-size:16px;\n\
                              font-family:verdana,arial;\n\
                              line-height:14px;\n\
                      }\n\
                      \n\
                      p {\n\
                          font:12px/21px verdana, arial, helvetica, sans-serif;\n\
                              margin:10px 0px 16px 0px;\n\
                              padding:0px 10px 0px 20px;\n\
                      }\n\
#Content>p {margin:0px,0px,0px,0px;}\n\
#Content>p+p {text-indent:0px;}\n\
                      \n\
                      pre {\n\
                          width: 100%;\n\
                          line-height:1.2em;\n\
                          border: 0px solid #aaa;\n\
                          padding: 0;\n\
                          font-size:13px;\n\
                          margin: 0px 0px 0px 0px;\n\
                      }\n\
                      \n\
                      a{\n\
                          color:blue;\n\
                          font-size:11px;\n\
                          text-decoration:none;\n\
                          font-weight:600;\n\
                          font-family:verdana, arial, helvetica, sans-serif;\n\
                      }\n\
a:link {color:blue;}\n\
                      a:visited {color:blue;}\n\
                      a:hover {background-color:#eee;}\n\
                      \n\
                      a.menulist {\n\
                          color:#09c;\n\
                          font-size:12px;\n\
                          text-decoration:none;\n\
                          font-weight:600;\n\
                          font-family:verdana, arial, helvetica, sans-serif;\n\
                      }\n\
                      a.menulist:link {color:#09c;}\n\
                      a.menulist:visited {color:#05e;}\n\
                      a.menulist:hover {background-color:#0174DF;}\n\
                      a.menulist:hover {color:white;}\n\
                      \n\
                      \n\
                      li{\n\
                          list-style:circle inside;\n\
                          font-size:14px;\n\
                          line-height:24px;\n\
                      }\n\
#Header {\n\
margin:10px 5px 10px 5px;\n\
    padding:17px 0px 20px 20px;\n\
    height:33px; \n\
    border:0px solid;\n\
    border-color:#365E2E;\n\
    border-bottom-width:15px;\n\
    line-height:11px;\n\
    background-color:#fff;\n\
    height:14px; \n\
}\n\
                      body>#Header {height:24px;}\n\
                      \n\
#Content {\n\
    margin:10px 50px 50px 180px;\n\
    padding:20px 5px 10px 15px;\n\
    border:0px groove #92ADF7;\n\
    border-left-width:8px;\n\
}\n\
                      \n\
#Menu {\n\
    position:absolute;\n\
    height:100%;\n\
    top:120px;\n\
    left:10px;\n\
    padding:10px 5px 10px 5px;\n\
    background-color:#fff;\n\
    border:0px solid #92ADF7;\n\
    border-right-width:0px;\n\
    line-height:17px;\n\
    width:150px;\n\
}\n\
                      body>#Menu {width:160px;}\n\
                      ";/*}}}*/

                      using namespace std;
void PrintHelp()
{
    fprintf(stdout,"usage: report_zincpred [options] zinc-pred-output-file\n");
    fprintf(stdout," Note: The input file is the output of znpred-postscan\n");
    fprintf(stdout,"options:\n");
    fprintf(stdout,"  --seqfile <file>    : supply the amino acid sequence file, can be multi-seq\n");
    fprintf(stdout,"  --parafile <file>   : supply file storing parameters for PREDZINC\n");
    fprintf(stdout,"  -o|--out <outfile>  : output the result to outfile, default = stdout\n");
    fprintf(stdout,"  --c-score  real     : threshold for zinc prediction score, default = 0.45\n"); 
    fprintf(stdout,"  --ishtml   y|n      : whether output HTML file, default = yes\n"); 
    fprintf(stdout,"  -h|--help           : print this help message and exit\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"Created on 2007-10-16, updated 2011-10-12, Nanjiang Shu\n");
}
void PrintVerboseHelp() { }
void GetPredZincVersion(string rundir)/*{{{*/
{
    FILE *fp = NULL;
    string versionFile = rundir + string("/../version");
    fp = fopen (versionFile.c_str(),"r");
    if (fp != NULL){
        fgetline(fp, predzinc_version, 100);
        fclose(fp);
    }

}/*}}}*/
void PrintReference(FILE *fpout) /*{{{*/
{
    fprintf(fpout,"Reference: \n");
    fprintf(fpout,"   Shu, N., Zhou, T. and Hovmoller, S. (2008) Prediction of zinc-binding\n");
    fprintf(fpout,"   sites in proteins from sequence, Bioinformatics, 24, 775-782.\n");
    fprintf(fpout,"\n");
}/*}}}*/
void PrintReferenceHTML(FILE *fpout) /*{{{*/
{
    char reference[]="Shu, N., Zhou, T. and Hovmoller, S. (2008) Prediction of zinc-binding sites in proteins from sequence, Bioinformatics, 24, 775-782.";
    fprintf(fpout,"<h3 class=\"subtitle\">Reference:</h3>\n");
    fprintf(fpout,"<p>\n");
    fprintf(fpout,"%s", reference);
    fprintf(fpout,"</p>\n");
}/*}}}*/
void WriteAASeq(const char *aaSeq, FILE *fpout) /*{{{*/
{
    fprintf(fpout,"Your input amino acid sequence (containing %d aa) is\n", int(strlen(aaSeq)));
    WriteFastaSeq((char*)aaSeq,fpout,0,LONGEST_SEQ*3,lineLength );
    fprintf(fpout,"\n");
}/*}}}*/
void WriteAASeqHTML(char *aaSeq, FILE *fpout) /*{{{*/
{
    fprintf(fpout,"<h3 class=\"subtitle\">Your input amino acid sequence (containing <font color=\"#FF00FF\">%d</font> aa) is\n", int(strlen(aaSeq)));
    fprintf(fpout,"<pre>\n");
    WriteFastaSeq(aaSeq,fpout,0,LONGEST_SEQ*3,lineLength );
    fprintf(fpout,"</pre>\n");
}/*}}}*/
void PrintParameter(const char *parafile, FILE *fpout) /*{{{*/
{
    fprintf(fpout,"Parameters:\n");
    fprintf(fpout,"\n");
    int buff;
    FILE *fp = fopen (parafile, "r");
    if (fp != NULL)
    {
        while ( (buff = fgetc(fp)) != EOF)
        {
            fputc(buff, fpout);
        }
        fprintf(fpout,"\n");
        fclose(fp);
    }
}/*}}}*/
void PrintParameterHTML(const char *parafile, FILE *fpout) /*{{{*/
{
    fprintf(fpout,"<h3 class=\"subtitle\">Parameters:</h3>\n");
    fprintf(fpout,"<p>\n");
    int buff;
    FILE *fp = fopen (parafile, "r");
    if (fp != NULL)
    {
        while ( (buff = fgetc(fp)) != EOF)
        {
            fputc(buff, fpout);
            if (buff=='\n'){
                fprintf(fpout,"<br>");
            }
        }
        fprintf(fpout,"\n");
        fclose(fp);
    }
    fprintf(fpout,"</p>\n");
    fprintf(fpout,"<br>\n");
}/*}}}*/

int ReportZnPred(const char *zincPredOutputFile, map <string,string> &aaSeqMap, const char* parafile, FILE *fpout = stdout)/*{{{*/
{
    GistPredChain chain;
    InitGistPredChain(&chain);
    AllocGistPredChain(&chain, LONGEST_SEQ);
    int numRes = 0;

    FILE *fpPred = NULL;
    fpPred = fopen(zincPredOutputFile,"r");
    checkfilestream(fpPred, zincPredOutputFile,"r");
    int numSite = 1;
    int operation = 0;
    operation = USE_AVERAGE;
    int i;


    fprintf(fpout,"Zinc-binding site prediction by %s (c) Shu.\n\n", predzinc_version);
    PrintReference(fpout);
    PrintParameter(parafile,fpout);
    fprintf(fpout,"\n");
    int cntChain = 0;
    while((numRes = GetGistPredResidue(fpPred, numSite, &chain, operation))!= EOF)
    {
        cntChain += 1;
        fprintf(fpout,"//BEGIN query %d\n", cntChain);

        if (aaSeqMap.find(chain.id) != aaSeqMap.end()){
            WriteAASeq(aaSeqMap[chain.id].c_str(),fpout);
        }

        Array1D <int> idx_1darray(chain.numRes+1);
        int *idx = idx_1darray.array1D;
        for (i = 0; i < chain.numRes+1; i ++) { idx[i] = i; }
        QuickSort_index(idx, chain.discriminant, 0, chain.numRes-1, DESCENDING);


        Array1D <int> idx_znres_1darray(chain.numRes);
        int *idx_znres = idx_znres_1darray.array1D;
        int cntZnRes = 0 ;

        for (i = 0 ; i < chain.numRes; i ++) {
            if (chain.discriminant[idx[i]] >= cutoff_znscore ) {
                idx_znres[cntZnRes] = idx[i];
                cntZnRes ++;
            }
        }

        if(cntZnRes <= 0) {
            fprintf(fpout,"No zinc-binding residues were predicted for protein \"%s\"\n", chain.id);
        } else {
            fprintf(fpout,"The following %d residues were predicted as zinc-binding for protein \"%s\" (with score >= %6.3lf), sorted by scores\n\n", cntZnRes, chain.id, cutoff_znscore);
            fprintf(fpout,"%3s %8s %6s\n", "Res", "SerialNo","Score");
            fprintf(fpout,"\n");
            for (i = 0 ; i < cntZnRes; i ++) {
                fprintf(fpout, "%3s %8d %6.3lf\n", AA1To3(chain.aaSeq[idx_znres[i]]), chain.aaSeqIndex[idx_znres[i]]+1, chain.discriminant[idx_znres[i]]);
            }
        }

        fprintf(fpout,"\n");
        fprintf(fpout,"Prediction scores for the rest %d selected residues, sorted by scores\n", chain.numRes-cntZnRes);
        fprintf(fpout,"\n");
        fprintf(fpout,"%3s %8s %6s\n", "Res", "SerialNo","Score");
        fprintf(fpout,"\n");

        for (i = 0 ; i < chain.numRes; i ++) {
            if (chain.discriminant[idx[i]] < cutoff_znscore){
                fprintf(fpout, "%3s %8d %6.3lf\n", AA1To3(chain.aaSeq[idx[i]]),
                        chain.aaSeqIndex[idx[i]]+1,
                        chain.discriminant[idx[i]]);
            }
        }
        fprintf(fpout, "//End\n");
    }
    if (fpPred != NULL ) {
        fclose(fpPred);
    }
    DeleteGistPredChain(&chain);
    return numRes;
}/*}}}*/
void WriteHTMLHeader(const char* title, FILE *fpout)/*{{{*/
{
    fprintf(fpout,"<html>\n");
    fprintf(fpout,"<head>\n");
    fprintf(fpout,"<title>%s</title>\n", title);
    fprintf(fpout,"<style type=\"text/css\" media=\"all\">\n");
    fprintf(fpout,"<!--\n");
    fprintf(fpout,"%s", layout_css);
    fprintf(fpout,"-->\n");
    fprintf(fpout,"</style>\n");
    fprintf(fpout,"</head>\n");
    fprintf(fpout,"<body>\n");
    fprintf(fpout,"<h2 class=\"headline\">Zinc-binding site prediction by %s (c) Shu.</h2>\n", predzinc_version);
}/*}}}*/
void WriteHTMLTail(FILE* fpout)/*{{{*/
{
    fprintf(fpout,"</head>\n");
    fprintf(fpout,"</html>\n");
}/*}}}*/
void WriteHTMLTableHeader(FILE *fpout)/*{{{*/
{
    char itemlist[][100]= {
        "Sequence with predicted ZB residues highlighted in <font color=\"red\">red</font>",
        "List of predicted ZB residues with scores",
        "Homologues used in the prediction"
    };
    fprintf(fpout,"<tr>\n");
    int i;
    for (i=0;i<3;i++){

        fprintf(fpout,"<th>\n");
        fprintf(fpout,"%s\n", itemlist[i]);
        fprintf(fpout,"</th>\n");
    }
    fprintf(fpout,"</tr>\n");
}/*}}}*/
void WriteHTMLTableContent(map <string,string> aaSeqMap, GistPredChain *pChain, int *idx_znres, int cntZnRes, char **homologIDList, double* homoscorelist, int numHomolog, FILE *fpout)/*{{{*/
{
    int i;

    set <int> aaSeqIndex_zinc;
    aaSeqIndex_zinc.clear();

    for (i=0;i<cntZnRes;i++){
        aaSeqIndex_zinc.insert(pChain->aaSeqIndex[idx_znres[i]]);
    }


    fprintf(fpout,"<tr>\n");


    /*  column 1 ===========================================*/
    fprintf(fpout,"<td  class=aaseq>\n");
    fprintf(fpout,"<p  class=aaseq>\n");
    if (aaSeqMap.find(pChain->id) != aaSeqMap.end()) {
        Array1D <char> aaSeq_1darray(aaSeqMap[pChain->id].size()+5);
        char *aaSeq=aaSeq_1darray.array1D;
        strcpy(aaSeq,aaSeqMap[pChain->id].c_str());

        int lengthseq=strlen(aaSeq);
        char font_color[100]="black";
        int font_size=1;
        bool isBolded=false;
        bool isZincRes=false;
        for (i=0;i<lengthseq;i++){
            if (strchr("CHDEchde", aaSeq[i])!= NULL){
                isBolded=true;
            }else{
                isBolded=false;
            }
            if (aaSeqIndex_zinc.find(i) != aaSeqIndex_zinc.end()){
                font_size=3;
                isZincRes=true;
                strcpy(font_color,"red");
            }else{
                font_size=1;
                isZincRes=false;
                strcpy(font_color,"black");
            }

            if (isBolded){
                fprintf(fpout,"<b>");
            }
            if (isZincRes){
                fprintf(fpout,"<font size=\"%d\" color=\"%s\">%c</font>", font_size, font_color, aaSeq[i]);
            } else {
                fprintf(fpout,"%c", aaSeq[i]);
            }
            if (isBolded){
                fprintf(fpout,"</b>");
            }
        }
    }
    fprintf(fpout,"\n");
    fprintf(fpout,"</p>\n");
    fprintf(fpout,"</td>\n");

    /*  column 2 ===========================================*/
    fprintf(fpout,"<td>\n");
    if(cntZnRes <= 0) {
        fprintf(fpout,"None");
    } else {
        fprintf(fpout,"<table>\n");
        fprintf(fpout,"<tr class=subtable>\n");
        fprintf(fpout,"<th>%s</th>\n", "No");
        fprintf(fpout,"<th>%s</th>\n", "AA");
        fprintf(fpout,"<th>%s</th>\n", "SeqIndex");
        fprintf(fpout,"<th>%s</th>\n", "ZnScore");
        fprintf(fpout,"</tr>\n");
        for (i = 0 ; i < cntZnRes; i ++) {
            fprintf(fpout,"<tr>\n");
            fprintf(fpout,"<td align=\"center\">%d</td>\n", i+1);
            fprintf(fpout,"<td align=\"center\">%s</td>\n", AA1To3(pChain->aaSeq[idx_znres[i]]));
            fprintf(fpout,"<td align=\"right\">%d</td>\n", pChain->aaSeqIndex[idx_znres[i]]+1);
            fprintf(fpout,"<td align=\"right\">%6.3lf</td>\n", pChain->discriminant[idx_znres[i]]);
            fprintf(fpout,"</tr>\n");
        }
        fprintf(fpout,"</table>\n");
    }
    fprintf(fpout,"</td>\n");

    /*  column 3 ===========================================*/
    char url[200]="";
    char pdbid[10]="";
    fprintf(fpout,"<td>\n");
    if (numHomolog <=0){
        fprintf(fpout,"None");
    }else{
        fprintf(fpout,"<table>\n");
        fprintf(fpout,"<tr class=subtable>\n");
        fprintf(fpout,"<th>%s</th>\n", "Homologues");
        fprintf(fpout,"<th>%s</th>\n", "Score");
        fprintf(fpout,"</tr>\n");

        for (i=0;i<numHomolog;i++){
            strncpy(pdbid,homologIDList[i],4);
            my_strupr(pdbid);
            fprintf(fpout,"<tr>\n");
            sprintf(url,"http://www.rcsb.org/pdb/explore/explore.do?pdbId=%s",pdbid);
            fprintf(fpout,"<td><a href=\"%s\" target=\"_blank\">%s</a></td>\n", url, homologIDList[i]);
            fprintf(fpout,"<td>%6.3lf</td>\n", homoscorelist[i]);
            fprintf(fpout,"</tr>\n");
        }
        fprintf(fpout,"</table>\n");
    }
    fprintf(fpout,"</td>\n");

    fprintf(fpout,"</tr>\n");
}
/*}}}*/
int GetHomologInfo(FILE* fp, char **idlist, double* scorelist)/*{{{*/
{
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;
    char idcan[100]="";
    double score=0.0;

    int status_pos = 0;
    fpos_t pos;
    int cnthomo=0;
    status_pos = fgetpos(fp, &pos);
    if (status_pos != 0){
        fprintf(stderr, "Warning! fgetpos() failed for GetHomologInfo\n");
    }
    char tmpstr1[100];
    char tmpstr2[100];
    int status_sscanf=0;
    while(fgetline(fp, line, maxline) != EOF){
        if (line[0]=='#'){
            if (strncmp("#Homolog",line, 8)==0){
                status_sscanf=sscanf(line,"%s %s %s %s %s %lf", tmpstr1,tmpstr2,tmpstr2,tmpstr2, idcan, &score);
                if(strcmp(tmpstr1,"#Homolog") == 0 && status_sscanf==6){
                    scorelist[cnthomo]=score;
                    strcpy(idlist[cnthomo],idcan);
                    cnthomo++;
                }
            }
        }else{
            status_pos=fsetpos(fp,&pos);
            if (status_pos != 0){
                fprintf(stderr, "Warning! fsetpos() failed for GetHomologInfo\n");
            }
            break;
        }
    }
    return cnthomo;
}/*}}}*/

int ReportZnPredHTML_new(const char *zincPredOutputFile, map <string,string> aaSeqMap, const char* parafile, FILE *fpout = stdout)/*{{{*/
{
    GistPredChain chain;
    InitGistPredChain(&chain);
    AllocGistPredChain(&chain, LONGEST_SEQ);
    int numRes = 0;
    //int numZnRes = 0;

    FILE *fpPred = NULL;
    fpPred = fopen(zincPredOutputFile,"r");
    checkfilestream(fpPred, zincPredOutputFile,"r");
    int numSite = 1;
    int operation = 0;
    operation = USE_AVERAGE;
    int i;

    int numHomolog=0;                    
    Array2D <char> homologIDList_2darray(100, 100);
    char ** homologIDList = homologIDList_2darray.array2D;
    Array1D <double> homoscorelist_1darray(100);
    double * homoscorelist= homoscorelist_1darray.array1D;

    numHomolog= GetHomologInfo(fpPred, homologIDList, homoscorelist);

    char title[]="Results of PredZinc";
    //WriteHTMLHeader(title, fpout);

    //PrintReferenceHTML(fpout);
    //PrintParameterHTML(parafile,fpout);

    fprintf(fpout,"<table border=\"1\" cellpadding=\"10\">\n");
    WriteHTMLTableHeader(fpout);
    while((numRes = GetGistPredResidue(fpPred, numSite, &chain, operation))!= EOF)
    {
        Array1D <int> idx_1darray(chain.numRes+1);
        int *idx = idx_1darray.array1D;
        for (i = 0; i < chain.numRes+1; i ++) { idx[i] = i; }
        QuickSort_index(idx, chain.discriminant, 0, chain.numRes-1, DESCENDING);

        Array1D <int> idx_znres_1darray(chain.numRes);
        int *idx_znres = idx_znres_1darray.array1D;
        int cntZnRes = 0 ;

        for (i = 0 ; i < chain.numRes; i ++) {
            if (chain.discriminant[idx[i]] >= cutoff_znscore ) {
                idx_znres[cntZnRes] = idx[i];
                cntZnRes ++;
            }
        }
        WriteHTMLTableContent(aaSeqMap, &chain, idx_znres, cntZnRes, homologIDList, homoscorelist, numHomolog, fpout);
        numHomolog= GetHomologInfo(fpPred, homologIDList, homoscorelist);
    }
    fprintf(fpout,"</table>\n");

    //fprintf(fpout, "<p>"); 
    //fprintf(fpout,"Predicted zinc-binding (ZB) residues are highlighted in red and with larger font size. Cys, His, Asp and Glu are bolded. Residues are predicted as zinc-binding if the score is >= %.3f<br>\n", cutoff_znscore);
    //fprintf(fpout,"CYS: cystein, HIS: histidine, ASP: aspartate, GLU: glutamate.<br>\n");
    //fprintf(fpout,"Clicking the ID of homologues links to the PDB website.<br>\n");
    //fprintf(fpout, "</p>"); 
    //WriteHTMLTail(fpout);

    if (fpPred != NULL )
    {
        fclose(fpPred);
    }
    DeleteGistPredChain(&chain);
    return numRes;
}/*}}}*/
int ReportZnPredHTML(const char *zincPredOutputFile, char* aaSeq, const char* parafile, FILE *fpout = stdout)/*{{{*/
{
    GistPredChain chain;
    InitGistPredChain(&chain);
    AllocGistPredChain(&chain, LONGEST_SEQ);
    int numRes = 0;
    //int numZnRes = 0;

    FILE *fpPred = NULL;
    fpPred = fopen(zincPredOutputFile,"r");
    checkfilestream(fpPred, zincPredOutputFile,"r");
    int numSite = 1;
    int operation = 0;
    operation = USE_AVERAGE;
    int i;

    fprintf(fpout,"<html>\n");
    fprintf(fpout,"<head>\n");
    fprintf(fpout,"<title>PREDZINC result in HTML format</title>\n");
    fprintf(fpout,"<style type=\"text/css\" media=\"all\" >@import \"../css/layoutSub.css\";</style>\n");
    fprintf(fpout,"</head>\n");
    fprintf(fpout,"<body>\n");

    fprintf(fpout,"<h2 class=\"headline\">Zinc-binding site prediction by %s (c) Shu.</h2>\n",predzinc_version);
    PrintReferenceHTML(fpout);
    PrintParameterHTML(parafile,fpout);
    fprintf(fpout,"\n");
    while((numRes = GetGistPredResidue(fpPred, numSite, &chain, operation))!= EOF)
    {
        WriteAASeqHTML(aaSeq,fpout);
        Array1D <int> idx_1darray(chain.numRes+1);
        int *idx = idx_1darray.array1D;
        for (i = 0; i < chain.numRes+1; i ++) { idx[i] = i; }
        QuickSort_index(idx, chain.discriminant, 0, chain.numRes-1, DESCENDING);


        Array1D <int> idx_znres_1darray(chain.numRes);
        int *idx_znres = idx_znres_1darray.array1D;
        int cntZnRes = 0 ;

        for (i = 0 ; i < chain.numRes; i ++) {
            if (chain.discriminant[idx[i]] >= cutoff_znscore ) {
                idx_znres[cntZnRes] = idx[i];
                cntZnRes ++;
            }
        }

        if(cntZnRes <= 0) {
            fprintf(fpout,"<h3 class=\"subtitle\">No zinc-binding residues were predicted for protein \"%s\"<h3>\n", chain.id);
        } else {
            fprintf(fpout,"<h3 class=\"subtitle\">The following <font color=\"#FF00FF\">%d</font> residues were predicted as zinc-binding for protein \"%s\" (with score >= <font color=\"#FF00FF\">%6.3lf</font>) </h3>\n", cntZnRes, chain.id, cutoff_znscore);
            fprintf(fpout,"<pre>\n");
            fprintf(fpout,"%3s %8s %6s\n", "Res", "SerialNo","Score");
            fprintf(fpout,"\n");
            for (i = 0 ; i < cntZnRes; i ++) {
                fprintf(fpout, "%3s %8d %6.3lf\n", AA1To3(chain.aaSeq[idx_znres[i]]), chain.aaSeqIndex[idx_znres[i]]+1, chain.discriminant[idx_znres[i]]);
            }
            fprintf(fpout,"</pre>\n");
        }

        fprintf(fpout,"\n");
        fprintf(fpout,"<h3 class=\"subtitle\">Prediction scores for the rest <font color=\"#FF00FF\">%d</font> selected residues, sorted by scores</h3>\n", chain.numRes-cntZnRes);
        fprintf(fpout,"<pre>\n");
        fprintf(fpout,"%3s %8s %6s\n", "Res", "SerialNo","Score");
        fprintf(fpout,"\n");

        for (i = 0 ; i < chain.numRes; i ++) {
            if (chain.discriminant[idx[i]] < cutoff_znscore){
                fprintf(fpout, "%3s %8d %6.3lf\n", AA1To3(chain.aaSeq[idx[i]]),
                        chain.aaSeqIndex[idx[i]]+1,
                        chain.discriminant[idx[i]]);
            }
        }
        fprintf(fpout, "//End\n");
    }
    fprintf(fpout,"</pre>\n");
    if (fpPred != NULL )
    {
        fclose(fpPred);
    }
    DeleteGistPredChain(&chain);
    return numRes;
}/*}}}*/
int main(int argc, char** argv)/*{{{*/
{
    bool isNonOptionArg = false;
    if(argc < 2) 
    {
        PrintHelp();
        return 0;
    }
    int i;
    char outfile[MAX_PATH+1] = "";
    char seqfile[MAX_PATH+1] = "";
    char parafile[MAX_PATH+1] = "";
    char zincPredOutputFile[MAX_PATH+1] = "";
    i = 1;
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
            else if( (strcmp(argv[i],"-seqfile") == 0) || (strcmp(argv[i], "--seqfile") == 0))  
            {
                if( ( i = option_parser_filename(argc, argv, i, seqfile)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"-parafile") == 0) || (strcmp(argv[i], "--parafile") == 0))  
            {
                if( ( i = option_parser_filename(argc, argv, i, parafile)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--ishtml") == 0))
            {
                ishtml = (strncasecmp(argv[i+1],"yes", 1) == 0 ) ? (true) : (false); 
                i = i + 2;
            }
            else if( (strcmp(argv[i], "--c-score") == 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, cutoff_znscore, true, 0.0, 1.0)) == -1)
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
            my_strcpy(zincPredOutputFile, argv[i], MAX_PATH);
            i ++;
        }
    }/*}}}*/

    if(strcmp(zincPredOutputFile, "") == 0) {
        fprintf(stderr, "Error! zinc-pred-output-file not supplied\n");
        return -1;
    }
    if(strcmp(seqfile, "") == 0) {
        fprintf(stderr, "Error! seqfile not supplied\n");
        return -1;
    }

    FILE *fpout = NULL;
    FILE *fpouthtml = NULL;
    char htmloutfile[MAX_PATH+1] = "";

    string rundir;
    Array1D <char> dir(500);
    getfilepath(argv[0], dir.array1D);
    rundir = dir.array1D;

    GetPredZincVersion(rundir);
    if(strcmp(outfile,"") == 0 || strcasecmp(outfile, "stdout") == 0) {
        fpout = stdout;
        fpouthtml = stdout;
    } else {
        fpout = fopen(outfile, "w");
        checkfilestream(fpout, outfile,"w");

        sprintf(htmloutfile,"%s.html", outfile);
        fpouthtml = fopen(htmloutfile, "w");
        checkfilestream(fpouthtml, htmloutfile,"w");
    }

    Array1D <char> aaSeq_1darray(LONGEST_SEQ*3);
    char* aaSeq=aaSeq_1darray.array1D;
    Array1D <char> annoLine_1darray(500);
    char* annoLine=annoLine_1darray.array1D;
    int seqType = 0;
    FILE *fpin = NULL;
    char seqid[500]="";
    fpin = fopen(seqfile, "r");
    map <string, string> aaSeqMap;
    if (fpin != NULL){
        while (ReadNextSeq_FASTA(fpin, aaSeq, &seqType, LONGEST_SEQ*3, annoLine, 500) != EOF)
        {   
            strcpy(seqid,"");
            sscanf(annoLine,"%s ", seqid);
            //fprintf(stderr,"\nseqid=<%s>\n", seqid);
            if (strcmp(seqid,"") != 0){
                aaSeqMap.insert(pair<string,string>(seqid, aaSeq));
            }
        }
        fclose(fpin);
    }

    //printf("numseq=%d\n", aaSeqMap.size());
    //printf("id=%s\n", "1AH7A");
    //printf("seq=%s\n", aaSeqMap["1AH7A"].c_str());

    ReportZnPred(zincPredOutputFile, aaSeqMap, parafile, fpout);
    if (ishtml) {
        ReportZnPredHTML_new(zincPredOutputFile, aaSeqMap, parafile, fpouthtml);
    }
    if(fpout != NULL && fpout != stdout) fclose(fpout);
    if(fpouthtml != NULL && fpouthtml != stdout) fclose(fpouthtml);

    return 0;
}
/*}}}*/
