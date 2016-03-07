//rootname.c
//return the rootname of a filename
//standalone version, that is, do not need to include my own library

#include <cstdio>
#include <cstring>
#include "myfunc.h"

void PrintHelp()
{
    fprintf(stdout,"usage: rootname filenames [more filenames...]\n");
    fprintf(stdout," accepting pipe\n");
    fprintf(stdout,"  if there are white spaces in file name, use qoutes \" \"\n");
    fprintf(stdout,"  -h | --help:    print this message\n");
    fprintf(stdout,"Nanjiang Shu, last modified 2007-01-24\n");
}
int main(int argc, char** argv)/*{{{*/
{
    char rtname[MAX_PATH+1] = "";
    int i ; 
    if(argc <2) 
    {
        //PrintHelp();
        //return 1;
        char filename[MAX_PATH+1] = "";
        while ( fgetdelim(stdin, filename, WHITE_SPACE, MAX_PATH)!= EOF) 
        {
            fprintf(stdout,"%s\n", rootname(filename,rtname));
        }
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
                fprintf(stdout,"%s\n", rootname(argv[i],rtname));
                i++;
            }
        }
    }
    return 0;
}
/*}}}*/
