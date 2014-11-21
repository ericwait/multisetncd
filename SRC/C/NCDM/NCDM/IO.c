
/******************************************************************************

This program, "NCDM", the associated MATLAB scripts and all 
provided data, are copyright (C) 2013-2014 Andrew R. Cohen and Paul
M. B. Vitanyi.  All rights reserved.

This program uses bzip2 compressor as a static library.
See the file SRC\C\bz2static\LICENSE.txt for details on that software.

This software may be referenced as:

A.R.Cohen and P.M.B. Vitanyi, "Normalized Compression Distance of Multisets 
with Applications," IEEE Transactions on Pattern Analysis and Machine 
Intelligence. 2014. In Press. Also arXiv:1212.5711.  

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

2. The origin of this software must not be misrepresented; you must 
   not claim that you wrote the original software.  If you use this 
   software in a product, an acknowledgment in the product 
   documentation would be appreciated but is not required.

3. Altered source versions must be plainly marked as such, and must
   not be misrepresented as being the original software.

4. The name of the author may not be used to endorse or promote 
   products derived from this software without specific prior written 
   permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Andrew R. Cohen acohen@coe.drexel.edu
Paul M. B. Vitanyi Paul.Vitanyi@cwi.nl
NCDM  version 1.0 of 13 March 2013
NCDM  version 2.0 (release) November 2014

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NCDM.h"


int ReadTrellis(char *szName,int FS,int N)
{
	FILE *fp;
	int k;
	char *pBufferIndex;
	int nidxTrue,nbytes,nret;
	char szPath[256];

	//sprintf(szPath,szName,FS,N);
	//sprintf(szPath,szName,N);
	sprintf(szPath,szName,FS);

	if (0==gnMPI_ID)
		printf("reading time series data from: %s\n",szPath);

	fp=fopen(szPath,"r");       
	if (NULL==fp)
	{	
		printf("!!!!!!!!!! open Trellis file %s failed\n  bye now...\n",szPath);
		return -1;
	}   		
	
	// read header
	fscanf(fp,"%d,%d\n",&gm,&gnTrellisDataSize);
	gnTrellisDataSize++; // one extra ';'
	if (0==gnMPI_ID)
		printf("got gm=%d\n",gm);
	if (gpbTrellisData)
		free(gpbTrellisData);
	gpbTrellisData = calloc(gnTrellisDataSize,sizeof(char));		
	if (NULL==gpbTrellisData)
	{
		printf("gpbTrellisData calloc fails!@#\n");
		return -1;
	}
	if (gpbInputBuffer)
		free(gpbInputBuffer);
	gpbInputBuffer = calloc(gnTrellisDataSize,sizeof(char));		
	if (NULL==gpbInputBuffer)
	{
		printf("gpbTrellisData calloc fails!@#\n");
		return -1;
	}

	if (gpbCompressBuffer)
		free(gpbCompressBuffer);
	gnCompressBufferSize = (int) (1.1*(double)gnTrellisDataSize);
	gpbCompressBuffer = calloc(gnCompressBufferSize,sizeof(char));		
	if (NULL==gpbCompressBuffer)
	{
		printf("gpbCompressBuffer calloc fails!@#\n");
		return -1;
	}
	// alloc structs
	if (gpTrellis)
		free(gpTrellis);
	gpTrellis = calloc(gm,sizeof(TRELLIS));		
	if (NULL==gpTrellis)
	{
		printf("gpTrellis calloc fails!@#\n");
		return -1;
	}

	gnN_Class=0;

	pBufferIndex = gpbTrellisData;
	for (k=0;k<gm;k++)
	{
		nret=fscanf(fp,"%d\n",&nidxTrue);
		nret=fscanf(fp,"%d\n",&nbytes);
		fgets(pBufferIndex,nbytes,fp);
		gpTrellis[k].nIdxTrue=nidxTrue;
		if (nidxTrue > gnN_Class)
			gnN_Class = nidxTrue;
		gpTrellis[k].pbStart = pBufferIndex;
		gpTrellis[k].nlen = nbytes;
		pBufferIndex[nbytes-1]=';';
		pBufferIndex += nbytes;
	} // k
	
	fclose(fp);
	if (0==gnMPI_ID)
		printf("got %d classes\n",gnN_Class);
	if (gpnClassLUT)
		free(gpnClassLUT);
	gpnClassLUT = calloc(gm,sizeof(int));
	if (NULL==gpnClassLUT)
	{
		printf("gpnClassLUT calloc fails!@#\n");
		return -1;
	}
	// initialize
	memset(gpnClassLUT,-1,gm*sizeof(int));
	for (k=0;k<gnN_Class;k++)
		gpnClassLUT[k]=k+1;

	return 0;

} // ReadTrellis

