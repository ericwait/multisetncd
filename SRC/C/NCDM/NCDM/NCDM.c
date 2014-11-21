
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
#include <time.h>
#include <mpi.h>

#include "NCDM.h"

char * gpbTrellisData; // raw data storage
char * gpbCompressBuffer; // for compression into
char * gpbInputBuffer; // copy input to ncdm here
int gnTrellisDataSize;
int gnCompressBufferSize;
int gnN_Class;
int gnWorkers;

char gpszTrellisFile[256]; // file path to the input
char gpszOutputFile[256]; // where to write results
int gnMinCluster;

  
int gm; // # of elements
PTRELLIS gpTrellis;
int gnMPI_ID;

// Compute the NCD(class)
//	nClass 1..nClasses, the class we want to include
//	iExclude - the element that we are cross-validating - don't include with any class
//	iTest - the element we are cross-validating. If this is set, check e(xA) and e(Ax) and take
//			minimum. If this is not set (== -1) then we are computing e.g. e(A)
//	nrgInclude - scractch space to store list of elements to process
double DistanceToClass(int nClass,int iExclude,int iTest,int *nrgInclude)
{
	int idxInclude;
	int i,j;
	double d,dmaxSubset,dPair,dPairMin;

	// eA
	dmaxSubset=0.0;
	memset(nrgInclude,-1,(gm+1)*sizeof(int));
	idxInclude = 1;	
	nrgInclude[0]=iTest; // eAx or eBx vs. eA or eB
			
	for (i=0;i<gm;i++)
	{
		if ((nClass==gpTrellis[i].nIdxTrue) && (i!=iExclude))
		{	
			nrgInclude[idxInclude]=i;
			idxInclude++;
		}
	}

	// idxInclude starts at 1 and points to the 1st empty slot in the nrgInclude list
	// if there is a test, idxInclude must be at least 2 (test + 1 elements)
	// if there is not a test, idxInclude must be at least 3 (2 elements)
	if ( ((-1 == iTest) && (idxInclude<=2)) || ((iTest>=0)&&(idxInclude<=1)))
		d=0.;
	else
		d=NCDM(nrgInclude,iTest,NULL);

	return d;

} // DistanceToClass

double DistanceToCluster(int nClass,int iExclude,int iTest,int *nrgInclude,int *nrgCluster,double *pdDenominator)
{
	int idxInclude;
	int i,j;
	double d;

	// eA
	memset(nrgInclude,-1,(gm+1)*sizeof(int));
	idxInclude = 1;	
	nrgInclude[0]=iTest; // eAx or eBx vs. eA or eB
			
	for (i=0;i<gm;i++)
	{
		if ((nrgCluster[i]==nClass) && (i!=iExclude))
		{	
			nrgInclude[idxInclude]=i;
			idxInclude++;
		}
	}
	d=NCDM(nrgInclude,iTest,pdDenominator);
	if (d<-1.)
		d=1.; // for e.g. one element lists, so we dont get -Inf
	return d;

} // DistanceToCluster

void MPI_Worker()
{
	int nret;
	int *nrgInclude=NULL;
	int *nrgCluster=NULL;
	int i,iTest,nPred,nMessage;
	double *pe,*pex;
	MPI_Status stat; 
	clock_t t;
	int FS,N;
	double *pdClass=NULL;
	double dMinClass,dSum;
	int nidxClassMin,iClass,jClass;
	double eA,eB,eAx,eBx,eMin;
	double d1,d2;

	while (1)
	{					
		nret=MPI_Recv(&nMessage, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,&stat);
		if (MESSAGE_ALLDONE == nMessage)
			break;
		else if (MESSAGE_NEWFS == nMessage)
		{
			nret=MPI_Recv(&FS, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,&stat);
			nret=MPI_Recv(&N, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,&stat);
			// ACK DON'T USE GM ABOVE HERE!
			nret = ReadTrellis(gpszTrellisFile,FS,N);
			if (nret<0)
				return;		
			if (NULL==nrgInclude)
			{
				nrgInclude = calloc(gm+1,sizeof(int));
				if (NULL==nrgInclude)
				{
					printf("ACK - out of memory\n");
					return;
				}		
			}
			if (NULL==nrgCluster)
			{
				nrgCluster = calloc(gm,sizeof(int));
				if (NULL==nrgCluster)
				{
					printf("ACK - out of memory\n");
					return;
				}		
			}
			// # classes bounded by gm, grows in partitioning
			if (NULL==pdClass)
			{
				pdClass = calloc(gm,sizeof(double));
				if (NULL==pdClass)
				{
					printf("ACK - out of memory\n");
					return;
				}
				pe = calloc(gm,sizeof(double));
				if (NULL==pe)
				{
					printf("ACK - out of memory\n");
					return;
				}
				pex = calloc(gm,sizeof(double));
				if (NULL==pex)
				{
					printf("ACK - out of memory\n");
					return;
				}
			}
		} // MESSAGE_NEWFS
		else if (MESSAGE_DIST_I == nMessage)
		{
			MPI_Recv(&iTest, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,&stat);
			dMinClass=1.;
			nidxClassMin=-1;
			

			for (i=0;i<gnN_Class;i++)
			{										
				pe[i] = DistanceToClass(i+1,iTest,-1,nrgInclude);
				pex[i] = DistanceToClass(i+1,iTest,iTest,nrgInclude);
				
				//pdClass[i]=(pex[i]-pe[i])/pex[i];
				pdClass[i]=(pex[i]-pe[i]);
				if ((gnMinCluster==0)  && (0==pdClass[i]))
					continue;

				if ((pdClass[i]>-1) && (pdClass[i]<dMinClass))
				{
					dMinClass=pdClass[i];
					nidxClassMin=i;
				}
			}	
			if (-1==nidxClassMin)
			{
				printf(" ****** !!!!!!! ACK ACK no class found - nidxClassMin == 0 iTest = %d !!!!!!!!! \n",iTest);
				nidxClassMin=0; // ?
			}
			dSum=pex[nidxClassMin-1]-pe[nidxClassMin-1];
			MPI_Send(&iTest, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			MPI_Send(&nidxClassMin, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			if (1==gnMPI_ID)
				printf("1:gnN_Class = %d\n",gnN_Class);

			for (i=0;i<gnN_Class;i++)
			{	
				MPI_Send(&pdClass[i], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);		
			}
		} // MESSAGE_DIST_I
		else if (MESSAGE_CLUSTERS == nMessage)
		{
			nret=MPI_Recv(nrgCluster, gm, MPI_INT, 0, 0, MPI_COMM_WORLD,&stat);
			//if (1==gnMPI_ID)
			//	printf("1: got clusters\n");
		}
		else if (MESSAGE_CLUSTER_I == nMessage)
		{
			
			nret=MPI_Recv(&iTest, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,&stat);
			//if (1==gnMPI_ID)
			//	printf("1: got itest=%d\n",iTest);

			eA = DistanceToCluster(1,iTest,-1,nrgInclude,nrgCluster,NULL);
			eAx = DistanceToCluster(1,iTest,iTest,nrgInclude,nrgCluster,NULL);
			
			eB = DistanceToCluster(2,iTest,-1,nrgInclude,nrgCluster,NULL);
			eBx = DistanceToCluster(2,iTest,iTest,nrgInclude,nrgCluster,NULL);

			MPI_Send(&iTest, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			MPI_Send(&eAx, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			MPI_Send(&eA, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			MPI_Send(&eBx, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			MPI_Send(&eB, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);


		//if (1==gnMPI_ID)
		//		printf("1: sent back %d, %lf\n",nPred,eMin);
		} // MESSAGE_CLUSTER_I
		else if (MESSAGE_UPDATE_LUT == nMessage)
		{
			nret=MPI_Recv(gpnClassLUT, gm, MPI_INT, 0, 0, MPI_COMM_WORLD,&stat);
		}
		else if (MESSAGE_UPDATE_CLASSES == nMessage)
		{
			nret=MPI_Recv(nrgInclude, gm, MPI_INT, 0, 0, MPI_COMM_WORLD,&stat);
			gnN_Class=0;
			for (i=0;i<gm;i++)
			{
				gpTrellis[i].nIdxTrue=nrgInclude[i];
				if (gpTrellis[i].nIdxTrue > gnN_Class)
					gnN_Class = gpTrellis[i].nIdxTrue;
			}
		}
		else if (MESSAGE_CLASS_DISTANCE == nMessage)
		{
			MPI_Recv(&iClass, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,&stat);
			MPI_Recv(&jClass, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,&stat);
			
			memset(nrgInclude,-1,(gm+1)*sizeof(int));	
			for (i=0;i<gm;i++)
				if ((gpTrellis[i].nIdxTrue==iClass) ) 
					nrgInclude[i]=i;
			eA	= NCDM(nrgInclude,-1,NULL);

			memset(nrgInclude,-1,(gm+1)*sizeof(int));
			for (i=0;i<gm;i++)
				if ((gpTrellis[i].nIdxTrue==jClass) ) 
					nrgInclude[i]=i;
			eB = NCDM(nrgInclude,-1,NULL);
			 
			memset(nrgInclude,-1,(gm+1)*sizeof(int));
			for (i=0;i<gm;i++)
				if ((gpTrellis[i].nIdxTrue==iClass) || (gpTrellis[i].nIdxTrue==jClass)  ) 
					nrgInclude[i]=i;
			eAx = NCDM(nrgInclude,-1,NULL);
			dMinClass =  eAx-eA-eB;
			MPI_Send(&dMinClass, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		} // MESSAGE_CLASS_DISTANCE

	}	// while(1)				
} // ReceiveWorker

int MPI_Master()
{
	int nDest,iTest,nret,ii,nMessage,iClass;
	int *nrgResults;
	MPI_Status stat; 
	int nError=0,nMinError;
	int FS,N;
	clock_t t;
	double dSum,dOptimal,dSumXVal;
	int nTrue;
	FILE *fpOut,*fpIn,*fpTime;
	char fname[256];
	int fx,nx;
	time_t result ;

	fpTime=fopen("NIST timing.txt","a");
	t = clock();
	result = time(NULL);
    fprintf(fpTime,"%d,%s,",gnMinCluster, ctime(&result));
	fflush(fpTime);

	fpOut=fopen(gpszOutputFile,"w");
	if (NULL==fpOut)
	{
		printf("fopen file %s failed - nret = %d\n",errno,fname);
	}
	printf("reading init\n");
	nret = ReadTrellis(gpszTrellisFile,3,2); // just to get gm
	printf("done reading init\n");
	if (nret<0)
		return -1;
	nrgResults = calloc(gm,sizeof(int));
	if (NULL==nrgResults)
	{
		printf("ACK - out of memory\n");
		return -1;
	}
	nMinError=gm;
	//gfpOut=fopen("sntnq.txt","w");
	//if (NULL==gfpOut)
	//	printf("ACK BAD output file\n");

	for (FS=3;FS<=3;FS++)
	//for (FS=20;FS<21;FS++)
	{

		if ((FS&1 && !(FS&2)) || (!(FS&1) && FS&2))
				continue;		
		for (N=2;N<=2;N++)
		{
				
			//t = clock();
			
			Partition( gpszTrellisFile, FS, N, gnMinCluster);
			result = time(NULL);
			fprintf(fpTime,"%s,",ctime(&result));
			fflush(fpTime);

			printf("back from partition\n");
			for (iClass=0;iClass<gnN_Class;iClass++)
			{				
					fprintf(fpOut,"%d,",gpnClassLUT[iClass]);
			}
			fprintf(fpOut,"\n");
			nError = 0; 
			// don't send FS - partition does it and then updates classes!

			// send FS
			//nMessage = MESSAGE_NEWFS;
			//for (nDest=1;nDest<gnWorkers;nDest++)
			//{
			//	MPI_Send(&nMessage, 1, MPI_INT, nDest, 0, MPI_COMM_WORLD);   
			//	MPI_Send(&FS, 1, MPI_INT, nDest, 0, MPI_COMM_WORLD); 
			//	MPI_Send(&N, 1, MPI_INT, nDest, 0, MPI_COMM_WORLD); 
			//}
			// send out distance requests
			nDest = 1;
			nMessage = MESSAGE_DIST_I;
			// use ii since we may not always want to test all gm elements
			for (ii=0;ii<gm;ii++)
			{	
				//if (ii&1)
				//	iTest = gm-ii-1;
				//else
					iTest = ii;

				//iTest=ii;
				MPI_Send(&nMessage, 1, MPI_INT, nDest, 0, MPI_COMM_WORLD); 
				MPI_Send(&iTest, 1, MPI_INT, nDest, 0, MPI_COMM_WORLD);                               
				nDest++;
				if (nDest>=gnWorkers)
					nDest=1;					
			}
			//nret = ReadTrellis(gpszTrellisFile,FS); 
			// printf("FS=%d N=%d  ||  ",FS,N);
			//dOptimal=GetFullClassDistance();
			// wait for answers
			printf("0:gnN_Class = %d\n",gnN_Class);
			nDest = 1;
			dSumXVal=0.;
			for (ii=0;ii<gm;ii++)
			{
				MPI_Recv(&iTest, 1, MPI_INT, nDest, 0, MPI_COMM_WORLD,&stat);
				MPI_Recv(&nret, 1, MPI_INT, nDest, 0, MPI_COMM_WORLD,&stat);

				for (iClass=0;iClass<gnN_Class;iClass++)
				{	
					MPI_Recv(&dSum, 1, MPI_DOUBLE, nDest, 0, MPI_COMM_WORLD,&stat);
					fprintf(fpOut,"%f,",dSum);
				}
				fprintf(fpOut,"\n");
				dSumXVal+=dSum;
				if (ii != iTest)
				{
					// uh oh!
					printf("\n\n\n *** ACK OH NO *** \n\n\n");
				}
				// nret is on [0,gn_NClass-1]
				nrgResults[iTest]=gpnClassLUT[nret];
				//
				
				if (gpTrellis[iTest].nIdxTrue>=0)
					nTrue = gpnClassLUT[gpTrellis[iTest].nIdxTrue-1];
				else
					nTrue=-1;
				if (nTrue>0 && nrgResults[iTest]!=nTrue)
					nError++;

				// printf("%d:(%d,%d), ",iTest,nret,gpTrellis[iTest].nIdxTrue);
				//if (iTest<=1000)
				//	printf("%d: got %d (raw %d) - true %d (raw %d), nError = %d\n",iTest,nrgResults[iTest],nret,nTrue,gpTrellis[iTest].nIdxTrue,nError);
				nDest++;
				if (nDest>=gnWorkers)
					nDest=1;
			}
			//printf("\n");
			t = clock() - t;
			if (nError<nMinError)
				nMinError=nError;
			printf("  nerrors=%d, elapsed time : %f seconds, MinErrors=%d\n\n\n",nError,((float)t)/CLOCKS_PER_SEC,nMinError);
			fflush(fpOut);
		} // N
	} // FS
	// send workers shutdown messages
	nret=MESSAGE_ALLDONE;
    for (nDest=1;nDest<gnWorkers;nDest++)
	{
		MPI_Send(&nret, 1, MPI_INT, nDest, 0, MPI_COMM_WORLD); 
		MPI_Send(&nret, 1, MPI_INT, nDest, 0, MPI_COMM_WORLD); 
	}

	result = time(NULL);
	fprintf(fpTime,"%s\n",ctime(&result));
	fclose(fpTime); 
	return 0;

} // MPI_Master
int main(int argc, char * argv[])
{
	int nret;
	clock_t t;
	time_t result ;


	// parse command lie
	strcpy(gpszTrellisFile,argv[1]);
	strcpy(gpszOutputFile,argv[2]);
	gnMinCluster=atoi(argv[3]);	
	
	setvbuf( stdout, NULL, _IONBF, 0 ); // from MPI documents - no buffering - output now !

	MPI_Init(&argc,&argv); /* all MPI programs start with MPI_Init; all 'N' processes exist thereafter */
	MPI_Comm_size(MPI_COMM_WORLD,&gnWorkers); /* find out how big the SPMD world is */
	MPI_Comm_rank(MPI_COMM_WORLD,&gnMPI_ID); /* and this processes' rank is */

	if (0==gnMPI_ID)
	{

		printf("MPI - init! numprocs = %d, myid = %d\n",gnWorkers,gnMPI_ID);		
		nret = MPI_Master();
		t = clock() - t;
		printf("\n\n\n nErrors = %d ----> Total elapsed time : %f seconds\n",nret,((float)t)/CLOCKS_PER_SEC);
	}
	else
	{
		MPI_Worker();
	}
	MPI_Finalize(); /* MPI Programs end with MPI Finalize; this is a weak synchronization point */

	result = time(NULL);
    printf("%s", ctime(&result));

	return 0;
}