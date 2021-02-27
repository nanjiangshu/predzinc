/*
 * =====================================================================================
 *       Filename:  myfunc.h
 *    Description:  source file for the library of common functions
 *        Version:  1.0
 *       Compiler:  g++
 *         Author:  Nanjiang Shu (Shu), nanjiang@struc.su.se
 *        Company:  Structural Chemistry, Stockholm Univesity
 * =====================================================================================
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <set>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <ctime>
#include <algorithm>

#include <regex.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "array.h"
#include "mytemplate.h"
#include "Constant.h"
#include "DataType.h"
#include "myfunc.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

using namespace std;

/*_io_stream*/
int checkfilestream(FILE *fp, const char* filename, const char *mode, bool isAssert /*=false*/)/*{{{*/
{
    if( fp == NULL)
    {
        fprintf(stderr,"Can not open file '%s' with mode '%s'\n", filename,mode);
        if(isAssert)
        {
            assert(fp != NULL);
        }
        return -1;
    }
    else
        return 0;
}
/*}}}*/

/*string operation*/
int my_strcpy(char* to, const char* from, int max, int sizefrom /*= 0*/)/*{{{*/
/******************************************************************************
 * my modification of strcpy
 * copy max characters from "from" to "to", add NULL terminator automatically
 * updated 2008-04-23, memcpy win in speed when the copying string is long
 * updated 2011-10-27:
 *   since calling strlen will be time consuming for very large from string,
 *   e.g. when "from" is the buffer of a whole trunk of file. Therefore, strlen
 *   part is removed.
 *****************************************************************************/
{
    if(sizefrom < 200) {
        strncpy(to,from,max);
    } else {
        memcpy(to,from,max);
    }
    to[max] = '\0';
    return max;
}/*}}}*/
char *my_strupr(char* str, int beg /* = 0*/, int end /* = 0x7FFFFFFF*/)/*{{{*/
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  my_strupr
 *  Description:  convert the string "from beg to end" to upper case, and return the string
 *  pointer. if called without beg and end, the whole string will be converted
 *  to uppercase
 * =====================================================================================
 */
{
	char *str0 = str;         /* preserve string pointer */
	int   size = strlen(str);
    if(end >= size) end = size -1;
    int diff = 'A' -'a';
    int idx = beg < 0 ? 0 : beg ;
    str += idx;
	while(idx <= end)
	{
		if(*str >= 'a' && *str <= 'z')
			*str += diff ;
		str ++; idx ++;
	}
	return str0;
}/*}}}*/
char *my_strlwr(char* str, int beg /* =0*/, int end /*= 0x7FFFFFFF*/)/*{{{*/
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  my_strlwr
 *  Description:  convert the string "from beg to end" to lower case, and return the string
 *  pointer. if called without beg and end, the whole string will be converted
 *  to lower
 * =====================================================================================
 */
{
	char *str0 = str;         /* preserve string pointer */
	int   size = strlen(str);
    if(end >= size) end = size -1;
    int diff = 'a' -'A';
    int idx = beg < 0 ? 0 : beg ;
    str += idx;
	while(idx <= end)
	{
		if(*str >= 'A' && *str <= 'Z')
			*str += diff ;
		str ++; idx ++;
	}
	return str0;
}/*}}}*/
char *strrtrim(char *str, const char *trim /*= " \t\n\r"*/)/*{{{*/
/*****************************************************************************
 * trim right (trailing) of the string, return the string after trimming
 * the str will be changed, however, use this function with the returned value,
 * e.g. 
 * char *newstr = strrtrim(origstr);
 ****************************************************************************/
{
    char  *end;
    if(!str) return NULL;

    end = str + strlen(str);
    while(end-- > str)
    {
        if(!strchr(trim, *end))
            return str;
        *end = 0;
    }
    return str;
}
/*}}}*/
char *strltrim(char *str, const char *trim /*= " \t\n\r"*/)/*{{{*/
/*****************************************************************************
 * trim left (leading) of the string, return the string after trimming
 * the str is not changed, use this function with the returned value,
 * e.g. 
 * char *newstr = strltrim(origstr);
 ****************************************************************************/
{
    if(!str) return NULL;
    /*if(!trim)*/
    /*    trim = " \t\r\n";*/
    while(*str)
    {
        if(!strchr(trim, *str))
            return str;
        ++str;
    }
    return str;
}
/*}}}*/
char *strtrim(char *str, const char *trim /*= " \t\n\r"*/)/*{{{*/
/*****************************************************************************
 * trim both the left and right part of the string, return the string after
 * trimming
 * only the end of the str is changed, use this function with the returned value,
 * e.g. 
 * char *newstr = strtrim(origstr);
 ****************************************************************************/
{
    return strltrim(strrtrim(str, trim), trim);
}
/*}}}*/
char *ssubstitute(char* str, char chForSub, char chAfterSub, int start /*= 0*/, int end /*= 0x7FFFFFFF */)/*{{{*/
/*****************************************************************************
 * substitute the character 'chForSub' in sub string str[start-end] by the
 * character 'chAfterSub'
 ****************************************************************************/
{
    int i;
	int   size = strlen(str);
    if(end >= size) end = size -1;
    for(i = start; i <= end; i++)
    {
        if(str[i] == chForSub) str[i] = chAfterSub;
    }
    return str;
}/*}}}*/
char *sreplace(char* to, char* from, int start /*= 0*/)/*{{{*/
/*****************************************************************************
 * replace characters in "to" by "from" from position "start", 
 * if strlen(from)+ start > strlen(to), stop at the end of string "to"
 ****************************************************************************/
{
    int  nTo = strlen(to);
    int  max;
    char *pStart;
    if( start >= 0 && start < nTo)
    {
        pStart = to + start;
        max = min(strlen(from), strlen(pStart));
        strncpy(pStart, from, max);
    }
    return to;
}/*}}}*/
char *SpanExcluding(const char* strForSpan,char* strAfterSpan, const char charSet[]/*=WHITE_SPACE*/)/*{{{*/
//**********************************************************************
// SpanExcluding()
// compress the string by removing the characters in the charSet
// and store the compressed string into strAfterSpan
// the search is case-sensitive
// default charSet is WHITE_SPACE
//**********************************************************************
{
	int l = strlen(strForSpan);
	int n = strlen(charSet);
	int i,j ; 
	j = 0 ;
	for(i = 0 ; i < l ;i ++)
	{
		if(!IsInCharSet(strForSpan[i],charSet, n))
		{
			strAfterSpan[j] = strForSpan[i] ;
			j ++ ;
		}
	}
	strAfterSpan[j] = '\0' ;
	return strAfterSpan;
}/*}}}*/
char *strchomp(char *str)/*{{{*/
// replace the \n with '\0'
{
    char *c = strchr(str, '\n');
    if(c) *c='\0';
    return str;
}
/*}}}*/

/*_io_ reading */
int fgetlinesize(FILE* fp)/*{{{*/
/*****************************************************************************
 * fgetlinesize()
 * get the number of characters from the current position to end_of_line
 ****************************************************************************/
{
    int nch = 0; /* record number of characters of the current line */
    int c;
    while((c = getc(fp)) != EOF)
    {
        if(c == 0x0d ) continue; /* in unix, '\n'= 0x0a, thus 0x0d will be cheated as another character*/ 
        if(c == '\n') break; /* in dos, '\n' is also equal to 0x0a, but the preceding 0x0d will not be read*/ 
        nch = nch + 1;
    }
    if(c == EOF && nch == 0) return EOF;
    else return nch;
}/*}}}*/
int fgetline(FILE* fp, char* line, int max/* = 0x7FFFFFFF*/)/*{{{*/
/*****************************************************************************
 * Read one line from fp, copying it to line array (but no more than max
 * chars). Does not place terminating \n in line array.  
 * it can be called without "max" flag, but should make sure the allocated
 * memory for "line" is larger than the longest line.
 *
 * Returns: line length, or 0 for empty line, or "EOF" for end-of-file.
 * 
 * since fgetc(), getc(), getchar(), returns int value,always use an int
 * variable to store the result of the fgetc(), getc() and getchar().  
 * getc() is faster than fgetc()
 *   LOG: 2006-04-26 16:31:29 Wednesday  Week 17 <nanjiang@shu>
 *   bug fixed for '\n' return keys, 
 *   now it is valid for both dos and unix
 ****************************************************************************/
{
    int nch = 0; /* record number of characters actually read */
    int c;
    max = max - 1;			/* leave room for '\0' */
    while((c = getc(fp)) != EOF)
    {
        if(c == 0x0d ) continue; /* in unix, '\n'= 0x0a, thus 0x0d will be cheated as another character*/ 
        if(c == '\n') break; /* in dos, '\n' is also equal to 0x0a, but the preceding 0x0d will not be read*/ 

        if(nch < max)
        {
            line[nch] = c;
            nch = nch + 1;
        }
    }
    line[nch] = '\0';

    if(c == EOF && nch == 0) return EOF;
    else return nch;
}/*}}}*/
int fgetdelim(FILE* fp, char* str, const char* delim /* = WHITE_SPACE*/, int max/* = 0x7FFFFFFF*/)/*{{{*/
/*****************************************************************************
 * read a string deliminated by any supplied character set
 ****************************************************************************/
{
    int nch = 0; /* record number of characters actually read */
    int c;
    max = max - 1;			/* leave room for '\0' */
    while((c = getc(fp)) != EOF)
    {
        if(IsInCharSet(c, delim))
            break;

        if(nch < max)
        {
            str[nch] = c;
            nch = nch + 1;
        }
    }
    str[nch] = '\0';

    if(c == EOF && nch == 0) return EOF;
    else return nch;
}/*}}}*/
int fgetlinecnt(FILE* fp, bool is_count_black_line /*= true*/)/*{{{*/
/*****************************************************************************
 * fgetlinecnt()
 * return the number of lines of the give file
 ****************************************************************************/
{
	int cnt = 0 ;
    int linesize;
    rewind(fp);
    while((linesize=fgetlinesize(fp)) != EOF) 
    {
        if(linesize >= 0 || is_count_black_line) cnt++;
    }
	return cnt ;
}/*}}}*/
int fgetlinecnt(const char* filename, bool is_count_black_line /*= true*/)/*{{{*/
{
	int cnt = 0 ;
    int linesize;
    FILE* fp = fopen(filename,"r");
    checkfilestream(fp, filename, "r", true);
    while((linesize=fgetlinesize(fp)) != EOF) 
    {
        if(linesize >= 0 || is_count_black_line) cnt++;
    }
    fclose(fp);
	return cnt ;
}/*}}}*/
int fgetlinecnt(const char* filename, int &maxline, bool is_count_blank_line /*= true*/)/*{{{*/
/*****************************************************************************
 * fgetlinecnt()
 * return the number of lines of file, also get the size of the longest line
 * updated 2009-07-15, print the filename when fp == 0
 ****************************************************************************/
{
    int   cnt;
    int   linesize;
    FILE *fp = fopen(filename,"r");
    checkfilestream(fp, filename, "r", true);
	cnt = 0 ;
    maxline = MIN_INT;
    while((linesize=fgetlinesize(fp)) != EOF) 
    {
        if(linesize != 0 || is_count_blank_line) cnt++;
        if(linesize > maxline) maxline = linesize;
    }
    fclose(fp);
	return cnt ;
}/*}}}*/
void f_neglect_comment(FILE* fp, const char comment_char /*= '#'*/)/*{{{*/
/*****************************************************************************
 * ignore / neglect the heading lines in a file which is started by a special
 * character, the default comment leading character is '#'
 ****************************************************************************/
{
    int status_fpos;
    fpos_t pos;
    int maxline = 50;
    char line[50+1] = "";
    
    status_fpos = fgetpos(fp,&pos);
    assert ( status_fpos == 0 );
    while(fgetline ( fp, line ,maxline) != EOF)
    {
        if(line[0] != comment_char)
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

/*file management*/
char *rootname(const char* filename, char* rtname, int max_rtname /*= MAX_PATH*/)/*{{{*/
/*****************************************************************************
 * rootname
 * given the file name, 
 * return the rootname of the filename
 ****************************************************************************/
{
    const char *pch;
    char *pstr;
    if((pch = strrchr(filename,'/')) != NULL)
        pstr = (char*) pch+1;
    else
        pstr = (char*) filename;

    if((pch = strrchr(pstr,'.')) != NULL)
        my_strcpy(rtname,pstr, min((int)(pch - pstr), max_rtname));
    else
        rtname = pstr;
    return rtname;
}
/*}}}*/
//char *rootname(const char *filename, char* rtname, int max_path = MAX_PATH )[>{{{<]
/*****************************************************************************
 * rootname
 * given the file name, 
 * return the rootname of the filename
 ****************************************************************************/
//{
//    char str[MAX_PATH+1] = "";
//    char* pch;
//    char* basename;
//    int n; 
//    int i;
//    my_strcpy(str, filename, MAX_PATH);
//    if( strcmp(str,"") == 0)
//        return NULL;
    
//    pch = strtok (str,"/");
//    while (pch != NULL)
//    {
//        basename = pch;
//        pch = strtok (NULL, "/");
//    }

//    if((pch = strrchr(basename, '.')) != NULL)
//    {
//        *pch = '\0';
//    }
        
//    my_strcpy(rtname,basename,max_path);

//    return rtname;
//}[>}}}<]
char *getfilepath(const char* filename, char *path, int max_path /*= MAX_PATH */)/*{{{*/
// return the file path
{
    const char *pch;
    if((pch = strrchr(filename,'/')) != NULL)
    {
        my_strcpy(path,filename, min((int)(pch - filename), max_path));
    }
    else
        my_strcpy(path,".",max_path);
    return path;
}
/*}}}*/
char *getfileext(const char* filename, char *ext, int max_ext /*= MAX_PATH */)/*{{{*/
// return the file extension
{
    const char *pch;
    char *pstr;
    if((pch = strrchr(filename,'.')) != NULL)
    {
        pstr = (char*) (pch + 1);
        my_strcpy(ext, pstr, max_ext);
    }
    else
        my_strcpy(ext,"",max_ext);
    return ext;
}
/*}}}*/

void VerifyFolder(const char* folder)/*{{{*/
/*****************************************************************************
 * VerifyFolder(const char* folder)
 * check the existance of folder, if the folder not exist, create it
 * else, do nothing
 ****************************************************************************/
{
    char command[MAX_COMMAND_LINE+1] = "";
    struct stat st;
    bool isFolderExist;
    int status_system;
    isFolderExist = ( stat(folder , &st ) == 0 );
    if (!isFolderExist )
    {
        sprintf(command,"%s %s", "mkdir", folder);
        status_system = system(command) ;
        assert(status_system != -1);
    }
}
/*}}}*/
int GetDataDir(char *datadir)/*{{{*/
/*****************************************************************************
 * Get DATADIR for the current host
 * return -1 if DATADIR not defined, and datadir set to /data as default
 * return 0 if succeed, set datadir as defined in DATADIR
 ****************************************************************************/
{
    char *pDATADIR = NULL;
    pDATADIR = getenv("DATADIR");
    if(pDATADIR == NULL || strcmp(pDATADIR,"") == 0)
    {
        fprintf(stderr,"Warning! Environmental variable DATADIR is NULL, set datadir to '/data'\n");
        my_strcpy(datadir,"/data",MAX_PATH);
        return -1;
    }
    else
    {
        my_strcpy(datadir,pDATADIR,MAX_PATH);
        return 0;
    }
}
/*}}}*/
int GetTMPDir(char *tmpdir)/*{{{*/
/*****************************************************************************
 * Get TMPDIR for the current host
 * return -1 if TMPDIR not defined, and tmpdir set to current dir
 * return 0 if succeed, set tmpdir as defined in TMPDIR
 ****************************************************************************/
{
    char *pTMPDIR = NULL;
    pTMPDIR = getenv("TMPDIR");
    if(pTMPDIR == NULL || strcmp(pTMPDIR,"") == 0)
    {
        char *pPWD = getenv("PWD");
        fprintf(stderr,"Warning! Environmental variable TMPDIR is NULL, set tmpdir to current dir %s\n", pPWD);
        my_strcpy(tmpdir,pPWD,MAX_PATH);
        return -1;
    }
    else
    {
        my_strcpy(tmpdir,pTMPDIR,MAX_PATH);
        return 0;
    }
}
/*}}}*/
int GetWorkDir(char *workdir)/*{{{*/
{
    char *pWORKDIR = NULL;
    pWORKDIR = getenv("WORKDIR");
    if(pWORKDIR == NULL || strcmp(pWORKDIR,"") == 0)
    {
        fprintf(stderr,"Warning! Environmental variable WORKDIR is NULL, set datadir to './'\n");
        my_strcpy(workdir,"./",MAX_PATH);
        return -1;
    }
    else
    {
        my_strcpy(workdir,pWORKDIR,MAX_PATH);
        return 0;
    }
}
/*}}}*/

char* GetPDBFilePath(const char* pdbid, char* pdbfilepath, const char pdbpath[] /*=""*/, const char pdbobsoletepath[] /*=""*/, const char pdbmodelspath[] /*=""*/, const char pdbmodelsobsoletepath[]/*=""*/)/*{{{*/
/*****************************************************************************
 * GetPDBFilePath() 
 * given the PDBID, return the filepath of pdb file 
 * return NULL, if the pdb file for the query pdbid does not exist, 
 ****************************************************************************/
{
    // char pdbpath[MAX_PATH] = "/data/pdb_dcp";
    // char pdbobsoletepath[MAX_PATH] = "/data/obsolete/pdb_dcp";
    char  subdir[4] = "";
    char  tmppdbid[SIZE_PDBID+1] = "";
    int   i;
    FILE *fp;
    char datadir[MAX_PATH+1]           = "";
    char c_pdbpath[MAX_PATH+1]         = "";
    char c_pdbobsoletepath[MAX_PATH+1] = "";
    char c_pdbmodelspath[MAX_PATH+1]         = "";
    char c_pdbmodelsobsoletepath[MAX_PATH+1] = "";

    if( strcmp( pdbpath,"") == 0 || strcmp(pdbobsoletepath,"") == 0 || strcmp( pdbmodelspath,"") == 0 || strcmp(pdbmodelsobsoletepath,"") == 0 ) {
        GetDataDir(datadir);
        if (strcmp(datadir,"")==0){
            fprintf(stderr,"DATADIR not set, please specify the location for PDB files\n");
            return NULL;
        }
    }


    if(strcmp(pdbpath,"") == 0) {
        sprintf(c_pdbpath,"%s/%s",datadir,"pdb/data/structures/divided/pdb_dcp");
    } else {
        my_strcpy(c_pdbpath, pdbpath,MAX_PATH);
    }

    if(strcmp(pdbobsoletepath,"") == 0) {
        sprintf(c_pdbobsoletepath,"%s/%s",datadir,"pdb/data/structures/obsolete/pdb_dcp");
    } else{
        my_strcpy(c_pdbobsoletepath, pdbobsoletepath,MAX_PATH);
    }
    if(strcmp(pdbmodelspath,"") == 0) {
        sprintf(c_pdbmodelspath,"%s/%s",datadir,"pdb/data/structures/models/current/pdb_dcp");
    } else {
        my_strcpy(c_pdbmodelspath, pdbmodelspath,MAX_PATH);
    }

    if(strcmp(pdbmodelsobsoletepath,"") == 0) {
        sprintf(c_pdbmodelsobsoletepath,"%s/%s",datadir,"pdb/data/structures/models/obsolete/pdb_dcp");
    } else{
        my_strcpy(c_pdbmodelsobsoletepath, pdbmodelsobsoletepath,MAX_PATH);
    }

    my_strcpy(tmppdbid,pdbid, SIZE_PDBID);
    my_strlwr(tmppdbid);
    for(i = 0 ; i < 2 ; i++) 
    { subdir[i] = tmppdbid[i+1]; }
    subdir[i] = '\0';

    sprintf(pdbfilepath,"%s/%s/pdb%s.ent",c_pdbpath,subdir,tmppdbid);	
    if((fp=fopen(pdbfilepath,"r"))==0) {
        sprintf(pdbfilepath,"%s/%s/pdb%s.ent",c_pdbobsoletepath,subdir, tmppdbid);	
        if((fp = fopen(pdbfilepath,"r"))==0) {
            sprintf(pdbfilepath,"%s/%s/pdb%s.ent",c_pdbmodelspath,subdir, tmppdbid);	
            if((fp = fopen(pdbfilepath,"r"))==0) {
                sprintf(pdbfilepath,"%s/%s/pdb%s.ent",c_pdbmodelsobsoletepath,subdir, tmppdbid);	
                if((fp = fopen(pdbfilepath,"r"))==0) {
                    strcpy(pdbfilepath,"");
                    fprintf(stderr, "no PDB file for PDBID %s\n",pdbid);
                    return (NULL);
                }
            }
        }
    }
    fclose(fp);
    return pdbfilepath;                            
}    /*}}}*/
char* GetMODMFilePath(const char* id, char* modmfilepath, const char modmpath[] /*= ""*/ )/*{{{*/
//***********************************************************************
//id: the unique chain identifier in nrPDB, for example "1D66 ", "1AC0A"
//id should be standardised, 
//***********************************************************************
{
	// char modmpath[MAX_PATH] = "/data/modm";
	char idtmp[SIZE_CHAIN_ID+1] = "";
    FILE *fp;

    char datadir[MAX_PATH+1]    = "";
    char c_modmpath[MAX_PATH+1] = "";

    if( strcmp( modmpath,"") == 0 )
    {
        GetDataDir(datadir);
        sprintf(c_modmpath,"%s/%s",datadir,"modm");
    }
    else
        my_strcpy(c_modmpath, modmpath,MAX_PATH);

	my_strcpy(idtmp,id,SIZE_CHAIN_ID);
    my_strupr(idtmp);
	if(idtmp[4] == ' ')  // because for the standardized id, chainID=' ' when
		idtmp[4] = '\0'; // chanID is empty, but for filename, it is '\0'
	sprintf(modmfilepath,"%s/%s.modm",c_modmpath,idtmp);
    if((fp=fopen(modmfilepath,"r"))==0)
    {
        fprintf(stderr, "no MODM file for id %s\n",id);
        return (NULL);
    }
    fclose(fp);
	return modmfilepath;
}/*}}}*/
char* GetPSSMFilePath(const char* id, char* pssmfilepath, const char pssmpath[] /*= ""*/)/*{{{*/
//***********************************************************************
// id: the unique chain identifier in nrPDB, for example "1D66 ", "1AC0A"
// return path of pssm file or
// return NULL if pssm not exist
//***********************************************************************
{
	char idtmp[SIZE_CHAIN_ID+1] = "";
    FILE *fp;

    char datadir[MAX_PATH+1]    = "";
    char c_pssmpath[MAX_PATH+1] = "";

    if( strcmp( pssmpath,"") == 0 )
    {
        GetDataDir(datadir);
        sprintf(c_pssmpath,"%s/%s",datadir,"pssm");
    }
    else
        my_strcpy(c_pssmpath, pssmpath,MAX_PATH);

	my_strcpy(idtmp,id,SIZE_CHAIN_ID);  
    my_strupr(idtmp);
	if(idtmp[4] == ' ')  // because for the standardized id, chainID=' ' when
		idtmp[4] = '\0'; // chanID is empty, but for filename, it is '\0'
	sprintf(pssmfilepath,"%s/%s.pssm",c_pssmpath,idtmp);
    if((fp=fopen(pssmfilepath,"r"))==0)
    {
        fprintf(stderr,"no PSSM file for id %s\n",id);
        return (NULL);
    }
    fclose(fp);
	return pssmfilepath;
}/*}}}*/
char* GetDSSPFilePath(const char* pdbid, char* dsspfilepath, const char dssppath[] /*= ""*/)/*{{{*/
//***********************************************************************
//input: pdbid
//return: dsspfilepath
//***********************************************************************
{
	char pdbidtmp[SIZE_PDBID+1] = "";
    FILE *fp;

    char datadir[MAX_PATH+1]    = "";
    char c_dssppath[MAX_PATH+1] = "";

    if( strcmp( dssppath,"") == 0 )
    {
        GetDataDir(datadir);
        sprintf(c_dssppath,"%s/%s",datadir,"dssp");
    }
    else
        my_strcpy(c_dssppath, dssppath,MAX_PATH);

	my_strcpy(pdbidtmp,pdbid,SIZE_PDBID); 
    my_strlwr(pdbidtmp);
	sprintf(dsspfilepath,"%s/%s.dssp",c_dssppath,pdbidtmp);
    if((fp=fopen(dsspfilepath,"r"))==0)
    {
        fprintf(stderr, "no DSSP file for PDBID %s\n",pdbid);
        return (NULL);
    }
	fclose(fp);
	return dsspfilepath;
}/*}}}*/
char* GetPDBAAFilePath(const char* id, char* pdbaafilepath, const char pdbaapath[] /*= ""*/)/*{{{*/
//***********************************************************************
//input : id       standardized chain identifier
//outout: pdbaafilepath
//***********************************************************************
{
	char pdbidtmp[SIZE_PDBID+1];
    char idtmp[SIZE_CHAIN_ID+1];
    char chainIDtmp;
    FILE *fp;

    char datadir[MAX_PATH+1]    = "";
    char c_pdbaapath[MAX_PATH+1] = "";

    if( strcmp( pdbaapath,"") == 0 )
    {
        GetDataDir(datadir);
        sprintf(c_pdbaapath,"%s/%s",datadir,"pdbaa");
    }
    else
        my_strcpy(c_pdbaapath, pdbaapath,MAX_PATH);

    chainIDtmp = toupper(id[4]);
    if(chainIDtmp == ' ') chainIDtmp = '\0';
        
	my_strcpy(pdbidtmp,id,SIZE_PDBID);
    my_strlwr(pdbidtmp);
    sprintf(idtmp,"%s_%c",pdbidtmp,chainIDtmp);
	sprintf(pdbaafilepath,"%s/%s.aa",c_pdbaapath,idtmp);
    if((fp=fopen(pdbaafilepath,"r"))==0)
    {
        fprintf(stderr, "no pdbaa file for id %s\n",id);
        return (NULL);
    }
	fclose(fp);
	return pdbaafilepath;
}/*}}}*/
char* GetSCOPAAFilePath(const char* id, char* scopaafilepath, const char scopaapath[] /*= ""*/)/*{{{*/
//***********************************************************************
//input : id       standardized chain identifier
//outout: pdbaafilepath
//***********************************************************************
{
    FILE *fp;
    char idtmp[SIZE_SCOP_ID+1] = "";

    char datadir[MAX_PATH+1]    = "";
    char c_scopaapath[MAX_PATH+1] = "";

    if( strcmp( scopaapath,"") == 0 )
    {
        GetDataDir(datadir);
        sprintf(c_scopaapath,"%s/%s",datadir,"scopaa");
    }
    else
        my_strcpy(c_scopaapath, scopaapath,MAX_PATH);
        
    my_strcpy(idtmp,id,SIZE_SCOP_ID);
    my_strlwr(idtmp);
	sprintf(scopaafilepath,"%s/%s.aa",c_scopaapath,idtmp);
    if((fp=fopen(scopaafilepath,"r"))==0)
    {
        fprintf(stderr, "no scopaa file for id %s\n",id);
        return (NULL);
    }
	fclose(fp);
	return scopaafilepath;
}/*}}}*/
char* GetSEQMAPFilePath(const char* id, char* seqmapfilepath, const char seqmappath[] /*= ""*/)/*{{{*/
{
    char idtmp[SIZE_CHAIN_ID+1] = "";
    FILE *fp;

    char datadir[MAX_PATH+1]      = "";
    char c_seqmappath[MAX_PATH+1] = "";

    if( strcmp( seqmappath,"") == 0 )
    {
        GetDataDir(datadir);
        sprintf(c_seqmappath,"%s/%s",datadir,"seqmap");
    }
    else
        my_strcpy(c_seqmappath, seqmappath,MAX_PATH);

	my_strcpy(idtmp,id,SIZE_CHAIN_ID);  
    my_strupr(idtmp);
	if(idtmp[4] == ' ')  // because for the standardized id, chainID=' ' when
		idtmp[4] = '\0'; // chanID is empty, but for filename, it is '\0'
	sprintf(seqmapfilepath,"%s/%s.seqmap",c_seqmappath,idtmp);
    if((fp=fopen(seqmapfilepath,"r"))==0)
    {
        fprintf(stderr,"no SEQMAP file for id %s\n",id);
        return (NULL);
    }
    fclose(fp);
	return seqmapfilepath; 
}/*}}}*/
char* GetShapeStringFilePath(const char* id, char* shapestringfilepath, const char shapestringpath[] /*= ""*/)/*{{{*/
{
    char idtmp[SIZE_CHAIN_ID+1] = "";
    FILE *fp;

    char datadir[MAX_PATH+1]      = "";
    char c_shapestringpath[MAX_PATH+1] = "";

    if( strcmp( shapestringpath,"") == 0 )
    {
        GetDataDir(datadir);
        sprintf(c_shapestringpath,"%s/%s/%s",datadir,"shapestring", "divided");
    }
    else
        my_strcpy(c_shapestringpath, shapestringpath,MAX_PATH);

	my_strcpy(idtmp,id,SIZE_CHAIN_ID);  
    //my_strupr(idtmp); /*to be case sensitive 2008-07-07*/
	if(idtmp[4] == ' ')  // because for the standardized id, chainID=' ' when
		idtmp[4] = '\0'; // chanID is empty, but for filename, it is '\0'
	sprintf(shapestringfilepath,"%s/%s.shapestring",c_shapestringpath,idtmp);
    if((fp=fopen(shapestringfilepath,"r"))==0)
    {
        fprintf(stderr,"no ShapeString file for id %s\n",id);
        return (NULL);
    }
    fclose(fp);
	return shapestringfilepath; 
}/*}}}*/




/*binary checking*/
bool IsBlankLine(const char *buf)/*{{{*/
{
    for ( ; *buf; buf++)
    {
        if (!isspace(*buf))
            return false;
    }
    return true;
}
/*}}}*/
bool IsInteger(double x)/*{{{*/
// whether the double value is quite close to the integer 
{
	int i ;
	int digit = 5 ;
	double temp ;
	x = fabs(x) ;
	for ( i = digit ; i >= 0 ; i--)
	{
		temp = pow(10.0 , double (i) );
		while ( x >= temp ) 
			x -= temp;
	}
	
	if( fabs(x) < 1e-5 || fabs(x - 1.0 ) < 1e-5 )
		return TRUE ;
	else 
		return FALSE ;
}
/*}}}*/
bool IsZero(Complex *cplx, int dim)/*{{{*/
//**********************************************************************
// check if the complex number is zero
//**********************************************************************
{
	int i ;
	for ( i = 0 ; i < dim ; i ++ )
	{
		if(fabs(cplx[i].x) > 0.0 || fabs(cplx[i].y) > 0.0) 
			return FALSE ;
	}
	return TRUE ;
}/*}}}*/
bool IsNumeric(const char* str)/*{{{*/
//**********************************************************************
//IsNumeric(char*)
//check if a string is a numeric number
{
	int i = 0;
    int cntdot = 0;
	while(str[i] != '\0')
	{
        if( i == 0)
        {
            if(!isdigit(str[i]) && str[i] != '+' && str[i] != '-' && str[i] != '.' )
                return false;
            if(str[i] == '.')
                cntdot ++;
        }
        else
        {
            if(! isdigit(str[i]) && str[i] != '.')
                return false;
            if(str[i] == '.')
                cntdot ++;
        }
        i ++;
	}
    if(cntdot >= 2)
        return false;
    else 
        return true;
}
/*}}}*/
bool IsDigit(const char* str)/*{{{*/
//**********************************************************************
//IsDigit(char*)
//check if a string contains only digit number [0-9]
//minus number or real number will return false
{
	int i = 0;
	while(str[i] != '\0')
	{
		if(! (str[i] >= '0' && str[i] <= '9'))
			return false;
		i ++;
	}
	return true;
}
/*}}}*/
bool IsDigit(const char ch)/*{{{*/
//check if a character is a digit number
{
	if( ch  >= '0' && ch  <= '9')
		return true;
	else
		return false;
}
/*}}}*/
bool IsLower(const char c)/*{{{*/
{
    if(c >= 'a' && c <= 'z')
        return true;
    else
        return false;
}
/*}}}*/
bool IsUpper(const char c)/*{{{*/
{
    if(c >= 'A' && c <= 'Z')
        return true;
    else
        return false;
}
/*}}}*/
bool IsInCharSet(const char ch, const char *charSet, int n /*= 0 */)/*{{{*/
/*****************************************************************************
 * check if the character "ch" is in charSet,
 ****************************************************************************/
{
	if(n == 0)
        n = strlen(charSet);
    int i;
	for(i = 0 ;i < n ; i ++)
	{
		if(ch == charSet[i])
			return true;
	}
	return false;
}/*}}}*/

/* mathematical functions */
void FFT(Complex *cplx, int dim, int pow2, bool isInv)/*{{{*/
//**********************************************************************
//FFT(), fast Fourier Transform
//**********************************************************************
// Fast Fourier Transform
// 2^pow2 = dim (number of compelx number)
// isInv : FALSE forward 
//		   TRUE  inverse
{
	int n2, j,  l, i, ib, k, k1, k2;
	int sgn;
	double tr, ti, arg, nu1;	// intermediate values in calcs. 
	double c, s;	       // cosine & sine components of Fourier trans. 


	n2 = dim / 2;
	nu1 = pow2 - 1.0;
	k = 0;
 
	sgn = isInv ? -1 : 1;	// sign change for inverse transform 	

	// Calculate the componets of the Fourier series of the function
	for( l = 0; l < pow2; l++ )
    {    
    	do    
	    {	
    		for( i = 0; i < n2; i++ )    
	        {	
    			j = int (k / ( pow( 2.0, nu1 ) ));    
				ib = BitSwap( j, pow2 );  //retern int ;	
	            arg = 2.0 * PI * ib / dim; 	
                c = cos( arg );     
                s = sgn * sin( arg );     // mutiply with sng ,sng == -1 for inverse. 
                k1 = k;     
                k2 = k1 + n2;     
                     
                tr = cplx[k2].x * c - cplx[k2].y * s ;    
    			ti = cplx[k2].x * s + cplx[k2].y * c ;    
	            cplx[k2].x = cplx[k1].x - tr; 	
                cplx[k2].y = cplx[k1].y - ti;     
                cplx[k1].x = cplx[k1].x + tr;     
                cplx[k1].y = cplx[k1].y + ti;            
                k++;             
    		}         
			k +=  n2;	
	     } while( k < dim - 1);	
         k = 0;    
         nu1 -= 1.0;    
         n2 /= 2;         
    }    
		
    for( k = 0; k < dim; k++ )   // return to original 
    {
        ib = BitSwap( k, pow2 );
        if( ib > k)
        {
			Swap( &cplx[k].x, &cplx[ib].x );
			Swap( &cplx[k].y, &cplx[ib].y );
        }    
	}

// If calculating the inverse transform, must divide the data by the number of
// data points.
    if( isInv )
	{
		double rn = 1.0 / dim ;
		for( k = 0; k < dim; k++)
		{
			cplx[k].x *= rn;
			cplx[k].y *= rn;
		}
	}

}/*}}}*/
int BitSwap(int i, int pow2)/*{{{*/
{
	int k ;
	int ib, i2;

	ib = 0;
    for( k = 0; k < pow2; k ++ )
    {
        i2  = i / 2;
        ib = ib * 2 + (i - 2 * i2);
        i = i2;
    }
    return( ib );
}/*}}}*/
int Integer(double x)/*{{{*/
//**********************************************************************
// Integer(double)
// contract the double value to integer according to ><0.5 law
//**********************************************************************
{
	double t ;
	t = x - int(x) ;
	if( x >= 0.0 )
	{
		if( t < 0.5 )
			return int(x);
		else
			return int(x) + 1 ;

	}
	else
	{
		if( t > -0.5 )
			return int(x) ;
		else
			return int(x) -1 ;
	}
}
/*}}}*/
void SmoothImage(unsigned short **image, int imageWidth, int imageHeight, int sizeN /*= 1*/)/*{{{*/
//****************************************************************************
// smoothing image by averaging neighbouring pixels
// sizeN -- size of neighbouring pixels, (2*size+1)x(2*size+1) pixels will be
// averaged, in default, sizeN = 1, thus, 3x3 pixels will be averaged
//*****************************************************************************
{
	int row, col;
	int i,j;
	float sum;
	//int size = 1;  // size of neighbour
	double weight = double (2*sizeN+1) * (2*sizeN+1) ;

	unsigned short ** imageTemp = NULL  ;
	imageTemp = Create2DArray(imageTemp, imageHeight + 1 , imageWidth + 1 );

	//remove the boundary
	for(row = sizeN + 1 ; row <= imageHeight - sizeN ; row ++)
		for(col = sizeN + 1 ; col <= imageWidth - sizeN ; col++ )
		{
			sum = 0.0;
			for( i = row - sizeN ; i <= row + sizeN ; i++)
				for( j = col - sizeN ; j<= col +sizeN ;j++)
					sum+=image[i][j];
			imageTemp[row][col] = (unsigned short) Integer(sum/weight);
		}
	// smoothing the image by averaging 3*3 neighbour
	for(row = sizeN+1 ; row <= imageHeight-sizeN ; row ++)
		for(col = sizeN+1 ; col <= imageWidth-sizeN ; col++ )
		{
			image[row][col] = imageTemp[row][col];
		}
	Delete2DArray(imageTemp , imageHeight +1 );
}/*}}}*/
void GetAmpPha(double re, double im, double& amp, double& pha)/*{{{*/
//**********************************************************************
//Get Phase from real and Imaginary , the phase is scaled in degree
//**********************************************************************
{
	amp = sqrt(re*re + im*im);
	if(amp < 1e-5) // if amp = 0 , pha is nonsense, so set to zero
		pha = 0.0;
	else
	{
		pha = acos(re / amp);
		pha *= (180.0 / PI);
		if(im < 0)
			pha =  -(pha) ;
	}
}/*}}}*/
double uniform_random(double min /*= 0.0*/, double max /*= 1.0*/)/*{{{*/
// return a uniform distributed random number in the range of (min,max)
{
    return ( (max-min) * double(rand())/( RAND_MAX+1.0) + min); 
}/*}}}*/
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
int GetNumDigit(int num)/*{{{*/
/*****************************************************************************
 * return the number of digits in num
 ****************************************************************************/
{
    if(num == 0) return 1;
    int i = 1;
    while((num/=10)!= 0)
    {
        i ++;
    }
    return i;
}
/*}}}*/


/*array operation, searching, sorting and shuffling*/
void Shuffle(int *array, int n, unsigned int rand_seed /*= time(NULL)*/)/*{{{*/
//**********************************************************************
//Shuffle an integer array with size n
//**********************************************************************
{
	int i ,j;
    srand(rand_seed);
	for(i = n-1 ; i > 0 ; i --)
	{
		j = rand() % i ;
		Swap(&array[i],&array[j]);
	}
}
/*}}}*/

template <class T> int BinarySearch_String(T keyStr, T* strs, int n) /*{{{*/
///*******************************************************
// BinarySearch()
// search a key string in a string database which is ordered by 
// alphabet
// n: the number of items of strs
// return the index of keyStr in strs if successful
// else return -1
///******************************************************
{
    int lo = 0;
    int hi = n - 1;
    int mid;
    int stat;
    /* Repeat while there are elements in range */
    while (lo <= hi) 
    {
        mid = (lo + hi) / 2; /* Compute midpoint */
        stat = strcmp(strs[mid],keyStr);
        if (stat == 0)            /* found element */
        { return(mid); } 
        else if (stat >0)      /* target is in first half */
        { hi = mid - 1; } 
        else /* array[mid] < target , target is in second half */
        { lo = mid + 1; }
    }
    return(-1); /* Nothing left in range; failure */
}
template int BinarySearch_String<char*> (char* keyStr, char** strs, int n);
template int BinarySearch_String<const char*> (const char* keyStr, const char** strs, int n);
 /*}}}*/
template <class T> int LinearSearch_String(T keyStr, T* strs, int n) /*{{{*/
///*******************************************************
// LinearSearch()
// search a key string in a string database linearly
// the order of the string database is not required
// n is the number of items in strs
// return the index of keyStr in strs if found
// return -1 if not found
///******************************************************
{
    int i;
	for(i = 0 ; i < n ; i ++)      
	{  
		if ( strcmp(keyStr,strs[i]) == 0 )
			return i ;
	} 
	return -1; 
}
template int LinearSearch_String<char*>(char* keyStr, char** strs, int n);
template int LinearSearch_String<const char*> (const char* keyStr, const char** strs, int n);
/*}}}*/

void BubbleSort(int* a, int n)/*{{{*/
///*****************************************************
// sort an integer array according to the ascending order
// a : array with integer data type
// n : size of array
{
	int i,j;
	for ( i = 0 ; i < n-1 ; i ++)
		for (j = n-1 ; j > i ; j --)
			if ( a[j-1] > a[j] )
			{
				Swap(&a[j-1],&a[j]);
			}
}/*}}}*/
void QuickSort_String(int *idx, char** strs, int low, int high)/*{{{*/
//*************************************************************
// QuickSort_String()
// Sort a string dataset alphabetically
// return the index array
// note that low and high refer the the index in a array
// so for example, for an array with n item, hight = n-1
//**************************************************************
{
    char* pivot;
    int m;
    int i;
    if(low < high)
    {
		Swap(&idx[low], &idx[(high+low)/2]);
	    pivot = strs[idx[low]];
		m = low;
	    for (i = low + 1; i <= high; i++)
		{
			if (strcmp(strs[idx[i]] , pivot) < 0)
			{
				m++;
                Swap(&idx[m], &idx[i]);
			}
	   }
	   Swap(&idx[low], &idx[m]);
	   QuickSort_String(idx, strs,low, m - 1);
	   QuickSort_String(idx, strs,m + 1, high);
	}
}/*}}}*/

int Grouping(char **strs, int numStrs, char **strGroup, int *subTotal,int SIZE_STRGROUP)/*{{{*/
/*****************************************************************************
 * grouping a list of strs, return the subtotal of each type of string, and
 * the number of groups
 * return numGroup
 ****************************************************************************/
{
    int i;
    int numStrGroup;
    int index;
    set <string> groups;
    set <string> :: iterator iss;
    for(i = 0 ; i < numStrs; i++)
        groups.insert(strs[i]);
    numStrGroup = groups.size();

    i = 0 ;
    for(iss = groups.begin(); iss != groups.end(); iss++)
    {
        my_strcpy(strGroup[i],(*iss).c_str(),SIZE_STRGROUP);
        i ++;
    }

    for( i = 0; i < numStrGroup ; i++) subTotal[i] = 0;

    for( i = 0 ; i < numStrs; i++)
    {
        if((index = BinarySearch_String(strs[i], strGroup,numStrGroup))!= -1)
            subTotal[index] ++;
    }
    
    return numStrGroup;
}
/*}}}*/


/*regular expression*/
int reg_findall(const char* string, const char* pattern,  regmatch_t * pmatch ,bool isOverlap /*= false*/)/*{{{*/
/*****************************************************
 *                                                   *
 * Description:  reguler expression, using gnu.regex *
 * string :     the string for searching
 * pattern:     regular expression 
 * pmatch :     pointer to the match structure, can be achieve (start,end)
 * isOverlap:   with overlapping or without overlapping ,default=false
 * return the number of matches
 * longest match
 *****************************************************/
{

    int i;
    regex_t regexp;
    int tag_regcomp = 0;
    int tag_regexec = 0;
    int cflags = REG_EXTENDED;
    int eflags = REG_NOTEOL;
    size_t nmatch = 1;
    regmatch_t match ;

    int n = strlen(string);
    // cout << "string =" << string << endl;
    // cout << "pattern= " << pattern << endl;
    char* substr = new char[n+1];
    strcpy(substr,string);
    if(n <= 0)
        return 0;

    tag_regcomp = regcomp(& regexp, pattern, cflags);
    if(tag_regcomp != 0)   // regexp compiling failed
    {
        // regfree(&regexp);
        cout << "pattern= " << pattern << endl;
        cout << "regexp compiling failed" << endl;
        return -1 ;
    }

    int numMatch = 0;
    int shift = 0;
    while(1)
    {
        tag_regexec = regexec(& regexp,  substr, nmatch, & match, eflags);
        if( tag_regexec == 0)
        {
            // cout << "numMatch= "<< numMatch << ":" << "shift=" << shift<< ":"<< substr<< endl;
            pmatch[numMatch].rm_so = match.rm_so + shift;
            pmatch[numMatch].rm_eo = match.rm_eo + shift;
            if(isOverlap)
            {
                shift += (match.rm_so +1 );
            }
            else
            {
                shift += match.rm_eo;
            }
            strcpy(substr,"");
            for(i = shift ; i < n ; i ++) { substr[i - shift] = string[i]; } substr[i-shift] = '\0';
            numMatch ++ ;
        }
        else
        {
            break;
        }
        if( strcmp (substr, "") == 0)
            break;
    }

    regfree(&regexp);
    delete [] substr;
    
    return numMatch;
}/*}}}*/


/*bioinformatics, programming on protein or nucleic database,*/

int Char2Digit(char aa, const char* alphabet, int n/* = 0*/)/*{{{*/
/****************************************************************************
 * Char2Digit()
 *
 * return the digital indexing of the character in alphabet, 
 * starting from 0
 * if the character is not in the alphabet, return -1
 ***************************************************************************/
{
    if( n == 0) 
        n = strlen(alphabet);
    const char *pch;
    if((pch = strchr(alphabet,aa))!= NULL)
    {
        int index = pch - alphabet;
        return index;
    }
    else
        return -1; // if not in alphabet
}/*}}}*/
int Charcase2Digit(char aa, const char* alphabet, int n /*= 0*/)/*{{{*/
/****************************************************************************
 * Charcase2Digit()
 *
 * return the digital indexing of the character in alphabet, case insensitive, 
 * starting from 0
 * if the character is not in the alphabet, return -1
 ***************************************************************************/
{
    if( n == 0) 
        n = strlen(alphabet);
    const char *pch;
    char aa_lower = tolower(aa);
    char aa_upper = toupper(aa);
    if((pch = strchr(alphabet,aa_lower))!= NULL 
            ||(pch = strchr(alphabet, aa_upper))!= NULL)
    {
        int index = pch - alphabet;
        return index;
    }
    else
        return -1; // if not in alphabet
}/*}}}*/
int Digit2Char(int dc, const char* alphabet, int n /*= 0*/)/*{{{*/
/****************************************************************************
 * Digit2Char()
 *
 * return the character in alphabet according to its indexing, 
 * starting from 0
 * if digit_char is out of range of the alphabet, return NULLCHAR
 ***************************************************************************/
{
    if( n == 0) n = strlen(alphabet);
    if( dc >= 0 && dc < n)
        return alphabet[dc];
    else
        return NULLCHAR; // if not in alphabet
}/*}}}*/


double Compute_ROC_score(int* label, int n)/*{{{*/
//*******************************************************************
// ROC analysis
// RFP analysis
//*******************************************************************
{
	int tp = 0;
	int fp = 0;
	__int64 iroc = 0;
	double roc = 0.0;

	int i ;
	for(i = 0 ; i < n ; i++)
	{
		if(label[i] ==1 )
			tp ++;
		else
		{
			fp ++;
			iroc += tp;
		}
	}

	if(tp == 0)
		roc = 0.0;
	else
	{
		if(fp == 0 )
			roc = 1.0;
		else
			roc = iroc /( double(tp)*fp);
	}

	return roc;
}/*}}}*/
double Compute_ROC50_score(int* label, int n)/*{{{*/
{	
	int tp = 0;
	int fp = 0;
	__int64 iroc50 = 0;
	double roc50 = 0.0;
	int i ;
	for(i = 0 ; i < n ; i++)
	{
		if(label[i] ==1 )
			tp ++;
		else
		{
			fp ++;
			iroc50 += tp;
			if(fp >= 50) 
				break;
		}
	}

	if(tp == 0)
		roc50 = 0.0;
	else
	{
		if(fp == 0 )
			roc50 = 1.0;
		else
			roc50 = iroc50 /( double(tp)*fp);
	}

	return roc50;
}/*}}}*/
double Compute_medianRFP_score(int* label, double* score, int n)/*{{{*/
{	
	int i;
	double median;
	double mRFP;
	
	int cnt = 0;
	Array1D <float> positiveScore(n);
	for(i = 0 ; i< n; i++)
	{
		if(label[i] == 1 ) 
		{
			positiveScore.array1D[cnt] = score[i];
			cnt ++;
		}
	}
	QuickSort(positiveScore.array1D,0,cnt-1);
	if((cnt%2) == 1)
		median = positiveScore.array1D[cnt/2];
	else
		median = (positiveScore.array1D[cnt/2-1] + positiveScore.array1D[cnt/2]) / 2.0 ;

	int negCnt = 0;
	int totCnt = 0;
	for(i = 0 ; i< n ; i++)
	{
		if(score[i] >= median)
		{
			if(label[i] != 1)
				negCnt ++;
			totCnt ++;
		}
	}
	mRFP = double(negCnt)/totCnt;
	return mRFP;
}/*}}}*/
double Compute_medianRFP50_score(int* label, double* score, int n)/*{{{*/
{
	int i;
	double median;
	double mRFP50;
	
	int cnt = 0;
	int n50;
	int fp = 0;
	Array1D <float> positiveScore(n);
	for(i = 0 ; i< n; i++)
	{
		if(label[i] == 1 ) 
		{
			positiveScore.array1D[cnt] = score[i];
			cnt ++;
		}
		else
		{
			fp ++;
			if(fp >= 50) break;
		}
	}
	n50 = i ;

	QuickSort(positiveScore.array1D,0,cnt-1);
	if((cnt%2) == 1)
		median = positiveScore.array1D[cnt/2];
	else
		median = (positiveScore.array1D[cnt/2-1] + positiveScore.array1D[cnt/2]) / 2.0 ;

	int negCnt = 0;
	int totCnt = 0;
	fp = 0;
	for(i = 0 ; i< n50 ; i++)
	{
		if(score[i] >= median)
		{
			if(label[i] != 1)
				negCnt ++;
			totCnt ++;
		}
	}
	mRFP50 = double(negCnt)/totCnt;
	return mRFP50;
}/*}}}*/

// functions for argument parser
int option_parser_filename(int argc, char **argv, int beg, char *filename)/*{{{*/
/*****************************************************************************
 * beg is the current index of argument list, e.g., argv[i] = "--out"
 ****************************************************************************/
{
    int i ; 
    bool isNonOptionArg = false;
    bool isFileNameSet = false;

    for(i = beg +1 ; i < argc ; i++)
    {
        if(argv[i][0] == '-' && strcmp(argv[i], "--") != 0 && !isNonOptionArg)
        {
            fprintf(stderr,"option '%s' must be followed by a filename, not option\n", argv[beg]);
            return -1;
        }
        else if(strcmp(argv[i], "--") == 0 && !isNonOptionArg)
        {
            isNonOptionArg = true;
        }
        else
        {
            my_strcpy(filename, argv[i], MAX_PATH);
            isFileNameSet = true;
            break;
        }
    }

    if(!isFileNameSet)
    {
        fprintf(stderr,"option '%s' must be followed by a filename\n", argv[beg]);
        return -1;
    }
    else 
    {
        return i+1;
    }
}
/*}}}*/
template <class T> int option_parser_numeric(int argc, char **argv, int beg, T &x, bool isRangeSet /*= false*/, T min /*= MIN_INT*/, T max /*= MAX_INT*/)/*{{{*/
/*****************************************************************************
 * beg is the current index of argument list, e.g., argv[i] = "--value"
 ****************************************************************************/
{
    int i ; 
    bool isValueSet = false;
    double tmp;
    i = beg +1;

    if (i < argc)
    {
        if(IsNumeric(argv[i]))
        {
            tmp = atof(argv[i]);
            if(isRangeSet)
            {
                if(tmp < min || tmp > max)
                {
                    fprintf(stderr,"Invalid value! Value after option '%s' must be in the range of [%g %g]\n", argv[beg], double(min), double(max));
                    return -1;
                }
            }
            x = T(tmp);
            isValueSet = true;
        }
    }

    if(!isValueSet)
    {
        fprintf(stderr,"option '%s' must be followed by a numerical value\n", argv[beg]);
        return -1;
    }
    else 
    {
        return i+1;
    }
}
template int option_parser_numeric<int>   (int argc, char **argv, int beg, int &x   , bool isRangeSet /*= false*/, int min /*= MIN_INT*/      , int max/* = MAX_INT*/);
template int option_parser_numeric<float> (int argc, char **argv, int beg, float &x , bool isRangeSet /*= false*/, float min /*= MIN_FLOAT*/  , float max /*= MAX_FLOAT*/);
template int option_parser_numeric<double>(int argc, char **argv, int beg, double &x, bool isRangeSet /*= false*/, double min/* = MIN_DOUBLE*/, double max /*= MAX_DOUBLE*/);
template int option_parser_numeric<int8>(int argc, char **argv, int beg, int8 &x, bool isRangeSet /*= false*/, int8 min/* = MIN_DOUBLE*/, int8 max /*= MAX_DOUBLE*/);
/*}}}*/

