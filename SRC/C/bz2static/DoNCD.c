#define BZIP2_STATIC
#include "bzlib.h"
#include <stdio.h>
#include <string.h>

#define COMPRESSBUFSIZE (10*1024+1024)
char gCompressBuffer[COMPRESSBUFSIZE];

double ComputeNCD(char *pBuffer,int n1,int n2)
{
    int c1,c2,c12;
    int dlen;
    int rval; 
    double NCD;

    dlen=COMPRESSBUFSIZE;
    rval=BZ2_bzBuffToBuffCompress(gCompressBuffer,&dlen,pBuffer,n1,9,0,30); 
    c1=dlen;
    
    dlen=COMPRESSBUFSIZE;
    rval=BZ2_bzBuffToBuffCompress(gCompressBuffer,&dlen,pBuffer+n1,n2,9,0,30); 
    c2=dlen;
    
    dlen=COMPRESSBUFSIZE;    
    rval=BZ2_bzBuffToBuffCompress(gCompressBuffer,&dlen,pBuffer,n2+n1,9,0,30); 
    c12=dlen;
   
    if (c1<c2)
        NCD=((double)(c12-c1))/c2;
    else
        NCD=((double)(c12-c2))/c1;
    return NCD;    
}



