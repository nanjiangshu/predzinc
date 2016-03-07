/*
 * =====================================================================================
 *        Filename:  getfilepath.cpp
 *     Description:  Return the path of the file
 *         Version:  1.0
 *         Created:  10/31/2006 06:08:07 PM CEST
 *        Compiler:  g++
 *          Author:  Nanjiang Shu (Shu), nanjiang@struc.su.se
 *         Company:  Structural Chemistry, Stockholm Univesity
 * =====================================================================================
 */
#include <cstdio>
#include <cstring>

#ifndef MAX_PATH
#define MAX_PATH 266
#endif

#ifndef HAVE_MIN_MAX
#define HAVE_MIN_MAX
#define MIN2(x,y)	((x) < (y) ? (x) : (y))
#define MIN3(x,y,z) (MIN2(x,y)<(z) ? (MIN2(x,y)) : (z))
#define MAX2(x,y)   ((x)<(y) ? (y) : (x))
#define MAX3(x,y,z) (MAX2(x,y)<(z) ? (z) : MAX2(x,y))
#endif
void PrintHelp()
{
    fprintf(stdout,"usage: getfilepath filenames [more filenames...]\n");
    fprintf(stdout,"  if there are white spaces in file name, use qoutes \" \"\n");
    fprintf(stdout,"  -h | --help:    print this message\n");
    fprintf(stdout,"Nanjiang Shu, last modified 2007-01-17\n");
}
int my_strcpy(char* to, const char* from, int max)/*{{{*/
/******************************************************************************
 * my modification of strcpy
 * copy max characters from "from" to "to", add NULL terminator automatically
 *****************************************************************************/
{
    int nch = 0;
    int n   = strlen(from);
    nch  = MIN2(n,max);
    memcpy(to,from,nch);
    to[nch] = '\0';
    return nch;
}
/*}}}*/
char *getfilepath(const char* filename, char *path, int max_path = MAX_PATH )/*{{{*/
// return the file path
{
    const char *pch;
    if((pch = strrchr(filename,'/')) != NULL)
    {
        my_strcpy(path,filename, MIN2(pch - filename, max_path));
    }
    else
        my_strcpy(path,".",max_path);
    return path;
}
/*}}}*/
int main(int argc, char** argv)/*{{{*/
{
    char path[MAX_PATH+1] = "";
    int i ; 
    if(argc <2) 
    {
        PrintHelp();
        return 1;
    }
    else
    {
        i = 1;
        while(i < argc)
        {
            if(strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
            {
                PrintHelp();
                break;
            }
            else
            {
                fprintf(stdout,"%s\n", getfilepath(argv[i],path));
                i++;
            }
        }
    }
    return 0;
}
/*}}}*/
