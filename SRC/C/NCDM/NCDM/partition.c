
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

#define MIN(a,b) (a<b?a:b)

int *gpnClassLUT;

//  http://members.cox.net/srice1/random/crandom.html
// random int on [0..m-1]
int randn(int m)
{
	 double r;  /* random value in range [0,1) */ 
	 double x;  /* random value in range [0,M) */ 
	 int y;  /* random integer in range [0,M) if M is an integer then range = [0,M-1] */  

	 /* r is a random floating point value in the range [0,1) {including 0, not including 1}. 
		 Note we must convert rand() and/or RAND_MAX+1 to floating point values to avoid integer division. 
		 In addition, Sean Scanlon pointed out the possibility that RAND_MAX may be the largest positive integer 
		 the architecture can represent, so (RAND_MAX+1) may result in an overflow, or more likely the value 
		 will end up being the largest negative integer the architecture can represent, 
		 so to avoid this we convert RAND_MAX and 1 to doubles before adding. */ 
		r = (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) ); 
    /* x is a random floating point value in the range [0,M) {including 0, not including M}. */                            
    x = (r * m); 
   /* y is a random integer in the range [0,M) {including 0, not including M}. If M is an integer 
		then the range is [0,M-1] {inclusive} */                             
    y = (int) x; 
	return y;

} // randn

void InitClusters(int k,int *nrgCluster,int nClass)
{
	int i,nSeed;

	for (i=0;i<gm;i++)
	{
		if (nClass==gpTrellis[i].nIdxTrue)
			nrgCluster[i]=CLUSTER_INCLUDE;
		else
			nrgCluster[i]=CLUSTER_EXCLUDE;
	}
	
	srand ( time(NULL) );
	// pick initial points
	for (i=0;i<k;i++)
	{
		do
			nSeed=randn(gm);
		while (nrgCluster[nSeed]==CLUSTER_EXCLUDE);
		nrgCluster[nSeed]=i+1; // 1 or 2 for now
		//printf("seed - %d\n",nSeed);
	}
	
	return;

} // InitClusters

void CreateNewClass(int nClass, int *nrgResults)
{
	int c1,c2;
	int i,nMessage,nDest;

	c1=0;c2=0;
	gnN_Class++;
	gpnClassLUT[gnN_Class-1]=gpnClassLUT[nClass-1];

	printf("splt %d : ");
	for (i=0;i<gm;i++)
	{
		if (2==nrgResults[i])
		{
			gpTrellis[i].nIdxTrue=gnN_Class;
			c2++;
		}
		else if (1==nrgResults[i])
		{
			c1++;
			printf("%d,",i);
		}
	}
	printf(" || ");
		for (i=0;i<gm;i++)
	{
		if (2==nrgResults[i])
		{
			printf("%d,",i);
		}
	}
		printf("\n");
}  // CreateNewClass


double GetPartitionClassDistance()
{
	int i,j;
	double dmin = 2.; // min distance between 2 classes
	double d;
	int nMessage;
	int nDest,nReceived;
	MPI_Status stat; 
	
	return -1.;
	printf("computing partitioned class distance = ");

	nDest = 1;
	nMessage = MESSAGE_CLASS_DISTANCE;
	for (i=1;i<=gnN_Class;i++)
		for (j=1;j<=gnN_Class;j++)
		{
			if (i==j)
				continue;
			if (gpnClassLUT[i-1]==gpnClassLUT[j-1])
				continue;
			MPI_Send(&nMessage, 1, MPI_INT, nDest, 0, MPI_COMM_WORLD); 
			MPI_Send(&i, 1, MPI_INT, nDest, 0, MPI_COMM_WORLD);
			MPI_Send(&j, 1, MPI_INT, nDest, 0, MPI_COMM_WORLD);
			nDest++;
			if (nDest>=gnWorkers)
				nDest=1;
		}
	
		nReceived=0;
	// wait for answers
	nDest = 1;
	for (i=1;i<=gnN_Class;i++)
		for (j=1;j<=gnN_Class;j++)
		{
			if (i==j)
				continue;
			if (gpnClassLUT[i-1]==gpnClassLUT[j-1])
				continue;
			MPI_Recv(&d, 1, MPI_DOUBLE, nDest, 0, MPI_COMM_WORLD,&stat);
			dmin+=d;
			nReceived++;
			nDest++;
			if (nDest>=gnWorkers)
				nDest=1;
		}

	dmin=dmin/(double)nReceived;
	printf("%lf\n",dmin);
	return dmin;
} // GetPartitionClassDistance


double GetFullClassDistance()
{
	int i,j;
	double dmin = 2.; // min distance between 2 classes
	double d;
	int nMessage;
	int nDest;
	MPI_Status stat; 
	
	printf("computing min class distance = ");

	nDest = 1;
	nMessage = MESSAGE_CLASS_DISTANCE;
	for (i=1;i<=gnN_Class;i++)
		for (j=i+1;j<=gnN_Class;j++)
		{
			MPI_Send(&nMessage, 1, MPI_INT, nDest, 0, MPI_COMM_WORLD); 
			MPI_Send(&i, 1, MPI_INT, nDest, 0, MPI_COMM_WORLD);
			MPI_Send(&j, 1, MPI_INT, nDest, 0, MPI_COMM_WORLD);
			nDest++;
			if (nDest>=gnWorkers)
				nDest=1;
		}
	
	// wait for answers
	nDest = 1;
	for (i=1;i<=gnN_Class;i++)
		for (j=i+1;j<=gnN_Class;j++)
	{
		MPI_Recv(&d, 1, MPI_DOUBLE, nDest, 0, MPI_COMM_WORLD,&stat);
		if (d<dmin)
			dmin=d;
		nDest++;
		if (nDest>=gnWorkers)
			nDest=1;
	}
	printf("%lf\n",dmin);
	return dmin;
} // GetFullClassDistance



int PartitionRep(int * pnc1,int * pnc2, double * pdeMin, double * pdeAB, int * nrgCluster,int *nrgResults, int *nrgInclude, int nClass,int nMinCluster)
{
	int nDest,iTest,ii,nret;
	int nMessage;
	MPI_Status stat; 
	int nError=0;
	int nIter=0;
	BOOL bDone=FALSE;
	clock_t t;
	int i;
	int idxMin;
	double eA,eB,eAx,eBx,eMin,eMinRep,ePred;
	double eMinA,eMinB;
	int nrep;
	int nc1,nc2;
	FILE * fpOut;

	fpOut=fopen("snt.txt","w");
	*pnc1=0;
	*pnc2=0;

	InitClusters(2,nrgCluster,nClass);
	nIter = 0;	
	bDone=FALSE;
	nc1=1;
	nc2=1; // start at 1 each - seeds!

	eMin = 3.; // eA+eB this iteration. minimize!
	eMinRep=eMin;eMinA=eMin;eMinB=eMin;
	while (!bDone && (nIter<20))
	{

		//t = clock();
		bDone = FALSE;
		// send clusters
		nMessage = MESSAGE_CLUSTERS;
		for (nDest=1;nDest<gnWorkers;nDest++)
		{
			MPI_Send(&nMessage, 1, MPI_INT, nDest, 0, MPI_COMM_WORLD);   
			MPI_Send(nrgCluster, gm, MPI_INT, nDest, 0, MPI_COMM_WORLD);
		}

		nDest = 1;
		nMessage=MESSAGE_CLUSTER_I;
		for (iTest=0;iTest<gm;iTest++)
		{
			if (gpTrellis[iTest].nIdxTrue != nClass)
				continue;
			MPI_Send(&nMessage, 1, MPI_INT, nDest, 0, MPI_COMM_WORLD); 
			MPI_Send(&iTest, 1, MPI_INT, nDest, 0, MPI_COMM_WORLD);                               
			nDest++;
			if (nDest>=gnWorkers)
				nDest=1;
		}

		if (0.==*pdeAB)
		{	// compute e(AB)
			printf("checking eab ...\n");
			memset(nrgInclude,-1,(gm+1)*sizeof(int));
			for (iTest=0;iTest<gm;iTest++)
			{					
				if (nClass==gpTrellis[iTest].nIdxTrue+1)
					nrgInclude[iTest]=iTest;
			}
			*pdeAB = NCDM(nrgInclude,-1,NULL); 
			//printf("class %d : eAB = %lf\n",nClass,*pdeAB);
		}
		
		// wait for answers		
		idxMin=-1;
		nDest = 1;

		for (ii=0;ii<gm;ii++)
		{
			if (gpTrellis[ii].nIdxTrue != nClass)
				continue;

			MPI_Recv(&iTest, 1, MPI_INT, nDest, 0, MPI_COMM_WORLD,&stat);
			if (ii != iTest)
			{
				// uh oh!
				printf("\n\n\n *** ACK OH NO *** \n\n\n");
			}
			if (gpTrellis[iTest].nIdxTrue!=nClass)
			{
				printf(" ACK ACK - PartitionRep :: got bad iTest - %d from nDest %d\n",iTest,nDest);
			}
					
			MPI_Recv(&eAx, 1, MPI_DOUBLE, nDest, 0, MPI_COMM_WORLD,&stat);	
			MPI_Recv(&eA, 1, MPI_DOUBLE, nDest, 0, MPI_COMM_WORLD,&stat);	
			MPI_Recv(&eBx, 1, MPI_DOUBLE, nDest, 0, MPI_COMM_WORLD,&stat);
			MPI_Recv(&eB, 1, MPI_DOUBLE, nDest, 0, MPI_COMM_WORLD,&stat);	
			
			//printf("%d,%d,%d,%f,%f,%f,%f\n",nIter,ii,nrgCluster[iTest],eA,eAx,eB,eBx);	
			if ((eAx-eA)<(eBx-eB))
			{
				ePred=eAx-eA;
				nret=1;
			}
			else
			{
				ePred=eBx-eB;
				nret=2;
			}

			if (nrgCluster[iTest]<1)
			{	// 1st pass - assign pairwise nearest neighbor to seeds				
				nrgCluster[iTest]=nret;	
				if (1==nret)
					nc1++;
				else
					nc2++;
			}
			else if ((nrgCluster[iTest]!=nret) && (nret>0))
			{	// wanna change - only the single best change gets made
			//if ((ePred>0.) && (ePred<eMin) && (eA>0.0) &&  (eB>0.) )
				if (ePred<eMin)
				{
					// printf("setting emin = %lf, nc1=%d,nc2=%d\n",eMin,nc1,nc2);
					eMin = ePred;
					eMinA = eA;
					eMinB = eB;
					idxMin = iTest;
				}
			}
			else if ((nrgCluster[iTest]==nret) && (nret>0))
			{	
				// no change
			}
			nDest++;
			if (nDest>=gnWorkers)
				nDest=1;
		} // for iTest

		if (nIter>0)
		{
			if ((3.==eMin) || (-1 == idxMin))
			{
				// printf("got back dMin=-2, idxMin=%d exiting\n",idxMin);
				bDone = 1;
				//printf("DONE : ");
			}
			else
			{
				// toggle the minimum
				//printf("toggling %d: from %d to %d || ",idxMin,nrgCluster[idxMin],2-nrgCluster[idxMin]+1);
				nrgCluster[idxMin]=2-nrgCluster[idxMin]+1;
				if (1==nrgCluster[idxMin])
				{
					nc1++;
					nc2--;
				}
				else
				{
					nc2++;
					nc1--;
				}
				eMinRep = eMin;
			}
		}
		//		
		printf(" eMinRep=%lf, nIter=%d,time=%f seconds, eMin=%lf\n",eMinRep,nIter,((float)t)/CLOCKS_PER_SEC,
			eMin);
				
		nIter++;
	} // while !bDone
	//printf("%lf: %lf (%d),%lf (%d) \n",eMinRep,eMinA,nc1,eMinB,nc2);
	if ( (eMinRep<*pdeMin) && (nc1>=nMinCluster) && (nc2>=nMinCluster) )
	{
		*pnc1=nc1;
		*pnc2=nc2;
		*pdeMin=eMinRep;
		//printf("setting *pdeMin to %lf\n",eMinRep);
		memcpy(nrgResults,nrgCluster,gm*sizeof(int));
	}

	return 1;
} // PartitionRep

int Partition(char * pszTrellisFile,int FS,int N,int nMinCluster)
{
	int nDest,iTest,ii,nret;
	int *nrgCluster;
	int *nrgResults; // current best
	int *nrgInclude;
	int nMessage;
	MPI_Status stat; 
	int nError=0;
	BOOL bDone=FALSE;
	clock_t t;
	int i;
	int idxMin;
	double eAB,eMinBest;
	int nrep;
	int nClass;
	int nc1,nc2;
	double dThresh;

	t = clock();
	nret = ReadTrellis(pszTrellisFile,FS,N); // just to get gm
	if (nret<0)
		return -1;
	nrgResults = calloc(gm,sizeof(int));
	if (NULL==nrgResults)
	{
		printf("ACK - out of memory\n");
		return -1;
	}
	nrgCluster = calloc(gm,sizeof(int));
	if (NULL==nrgCluster)
	{
		printf("ACK - out of memory\n");
		return -1;
	}
		
	nrgInclude = calloc(gm+1,sizeof(int));
	if (NULL==nrgInclude)
	{
		printf("ACK - out of memory\n");
		return;
	}	
	// send FS, N
	nMessage = MESSAGE_NEWFS;
	for (nDest=1;nDest<gnWorkers;nDest++)
	{
		MPI_Send(&nMessage, 1, MPI_INT, nDest, 0, MPI_COMM_WORLD);   
		MPI_Send(&FS, 1, MPI_INT, nDest, 0, MPI_COMM_WORLD); 
		MPI_Send(&N, 1, MPI_INT, nDest, 0, MPI_COMM_WORLD); 
	}
	if (nMinCluster>1)
	{
		dThresh = -1.; //GetFullClassDistance();
		eAB=0.;
		nClass = 1;
		while (nClass <= gnN_Class)
		{
			printf("partitioning class %d ... \n",nClass);
			// eMinBest is eA+eB. Since eAB is fixed and we want to maximize eAB-eA-eB, we minimize (eA+eB)
			eMinBest = 3.; 
			memset(nrgResults,-1,gm*sizeof(int));
			for (nrep=0;nrep<10;nrep++)
			{	
				PartitionRep(&nc1,&nc2, &eMinBest, &eAB, nrgCluster,nrgResults, nrgInclude,  nClass, nMinCluster);
			} // for nrep
			printf("\ndone class %d : eMinBest = %lf\n",nClass,eMinBest);
			if (eMinBest<3.) // (((eAB-eMinBest) > dThresh))
				CreateNewClass(nClass,nrgResults);
			else
			{
				nClass++;	
			}
		} // nClass
	} // nMinCluster>1
	else if (0== nMinCluster)
	{
		//// revert to single nearest neighbor
		for (i=0;i<gm;i++)
		{
			gpnClassLUT[i]=gpTrellis[i].nIdxTrue;
			gpTrellis[i].nIdxTrue=i+1;
		}
		gnN_Class=gm;
	}
	// else 1==nMinCluster, unpartitioned multiples distance

	t = clock() - t;
	nMessage=MESSAGE_UPDATE_LUT;
	for (nDest=1;nDest<gnWorkers;nDest++)
	{
		MPI_Send(&nMessage, 1, MPI_INT, nDest, 0, MPI_COMM_WORLD);   
		MPI_Send(gpnClassLUT, gm, MPI_INT, nDest, 0, MPI_COMM_WORLD);   
	}

	for (i=0;i<gm;i++)
		nrgResults[i]=gpTrellis[i].nIdxTrue;
	nMessage=MESSAGE_UPDATE_CLASSES;
	for (nDest=1;nDest<gnWorkers;nDest++)
	{
		MPI_Send(&nMessage, 1, MPI_INT, nDest, 0, MPI_COMM_WORLD);   
		MPI_Send(nrgResults, gm, MPI_INT, nDest, 0, MPI_COMM_WORLD);   
	}
	//dThresh=GetPartitionClassDistance();
	printf("done partitioning, time = %f...sent updated LUT and classes --> partition class distance = %lf\n",((float)t)/CLOCKS_PER_SEC,dThresh);

	return 0;


} // Partition
