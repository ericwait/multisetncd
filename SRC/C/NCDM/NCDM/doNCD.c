
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



#define BZIP2_STATIC
#include "bzlib.h"
#include <stdio.h>
#include <string.h>
#include "NCDM.h"

int gnDenomXMax;

int Press(char *pBuffer, int nlen)
{
    
    unsigned int dlen;
    int rval; 

    dlen=gnCompressBufferSize;
    rval=BZ2_bzBuffToBuffCompress(gpbCompressBuffer,&dlen,pBuffer,nlen,9,0,30); 
    
	return dlen;

} // Press

// put a leading '[' and a trailing ']' on the buffer
void MarkBuffer(char * pBuffer,int nlen)
{
	*pBuffer='[';
	*(pBuffer+nlen-1)=']';
	return;
} // MarkBuffer

double Numerator(int * nrgInclude)
{
	int i;
	char *pBuffer;
	int nX,nx,nX_Size; // for numerator
	int nxmin,idx_xmin;
	int idx;
	double dNumerator;

	nxmin = MAXINT;
	gnDenomXMax=0;

	idx_xmin=-1;
	nX_Size=1;
	memset(gpbInputBuffer,0,gnTrellisDataSize);
	pBuffer = gpbInputBuffer+1;
	for (i=0;i<=gm;i++)
		if (nrgInclude[i]>=0)
		{
			idx = nrgInclude[i];
			memcpy(pBuffer,gpTrellis[idx].pbStart,gpTrellis[idx].nlen);
			MarkBuffer(gpbInputBuffer,gpTrellis[idx].nlen+1);
			nx=Press(gpbInputBuffer,gpTrellis[idx].nlen+1);
			if (nx<nxmin)
			{
				nxmin=nx;
				idx_xmin = i;
			}

			if (nx>gnDenomXMax)
				gnDenomXMax=nx;

			memset(gpbInputBuffer,0,gnTrellisDataSize);
		}

		for (i=0;i<=gm;i++)
			if (nrgInclude[i]>=0)
			{
				idx = nrgInclude[i];
				memcpy(pBuffer,gpTrellis[idx].pbStart,gpTrellis[idx].nlen);
				pBuffer += gpTrellis[idx].nlen;
				nX_Size += gpTrellis[idx].nlen;
			}
		MarkBuffer(gpbInputBuffer,nX_Size);
		nX = Press(gpbInputBuffer,nX_Size);	

		dNumerator = nX - nxmin;
		return dNumerator;
} // Numerator

double Denominator(int * nrgInclude)
{
	int i;
	char *pBuffer;
	int idx;
	int dxmax,idx_dxmax;
	int ndx; // element to exclude from denom
	int nXx,nX_Size; // size of denom with x excluded
	
	// ac 8 5 13
	//dxmax=(double)gnDenomXMax;
	//return dxmax;
 	
	// denominator - max(Z(X\x))
	ndx=0;
	dxmax=0;
	while (ndx<=gm)
	{
		while ((ndx<=gm) && (nrgInclude[ndx]<0))
			ndx++;
		if ((ndx<=gm) && (nrgInclude[ndx]>=0))
		{
				memset(gpbInputBuffer,0,gnTrellisDataSize);
				pBuffer = gpbInputBuffer+1;
				nX_Size = 1;
				// exclude ndx and compress
				for (i=0;i<=gm;i++)
					if ((nrgInclude[i]>=0) && (i!=ndx))
					{
						idx = nrgInclude[i];
						memcpy(pBuffer,gpTrellis[idx].pbStart,gpTrellis[idx].nlen);
						pBuffer += gpTrellis[idx].nlen;
						nX_Size += gpTrellis[idx].nlen;
					}
				
				MarkBuffer(gpbInputBuffer,nX_Size);
				nXx = Press(gpbInputBuffer,nX_Size);	
				if (nXx > dxmax)
				{
					idx_dxmax = ndx;
					dxmax = nXx;
				}
				ndx++;
		}
	}
	
	return dxmax;

} // Denominator
double  NCDM(int * nrgInclude,int iTest, double *pdDenominator)
{
	double n,n2;
	double d,d2;
	int idx=-1;
	int i,j;
	double dPairMin,dPair;

	n=Numerator(nrgInclude);
	d=Denominator(nrgInclude);
	if (NULL!=pdDenominator)
		*pdDenominator=d;

	if (0.==d)
		return 1.;

	if (iTest<0)
		return n/d;

	nrgInclude[0]=-1;
	nrgInclude[gm]=iTest;

	n2=Numerator(nrgInclude);
	d2=Denominator(nrgInclude);
	if ((NULL!=pdDenominator)&&(d2<d))
		*pdDenominator=d2;

	nrgInclude[gm]=-1;


	if (n2/d2<n/d)
		return n2/d2;
	else
		return n/d;
}
