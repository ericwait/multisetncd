
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


#define MESSAGE_DIST_I			-10
#define MESSAGE_NEWFS			-11
#define MESSAGE_NEWELEMENT		-12
#define MESSAGE_ALLDONE			-14
#define MESSAGE_CLUSTERS		-15
#define MESSAGE_CLUSTER_I		-16
#define MESSAGE_UPDATE_LUT		-17
#define MESSAGE_UPDATE_CLASSES	-18
#define MESSAGE_CLASS_DISTANCE	-19

#define CLUSTER_EXCLUDE		-1
#define CLUSTER_INCLUDE		-2

#define TRUE 1
#define FALSE 0

typedef int BOOL;
typedef struct _Trellis
{
	char *pbStart; // index into gpbTrellisData
	int nlen; // # bytes of trellis data
	int nIdxTrue;
} TRELLIS, *PTRELLIS;

#define ABS(x) ((x)>0.?x:-1.*(x))

extern char * gpbTrellisData; // raw data storage
extern int gm; // # of elements
extern PTRELLIS gpTrellis;
extern char * gpbCompressBuffer; // for compression into
extern int gnCompressBufferSize;
extern int gnTrellisDataSize;
extern char * gpbInputBuffer; // copy input to ncdm here
extern int gnMPI_ID;
extern int gnN_Class;
extern int gnWorkers;
extern int * gpnClassLUT;
extern FILE * gfpOut;
// IO.c
int ReadTrellis(char *szName,int FS,int N);

// doNCD.c
double  NCDM(int * nrgInclude,int iTest, double *pdDenominator);

int Partition(char * pszTrellisFile,int FS,int N,int nMinCluster);
int randn(int m);
