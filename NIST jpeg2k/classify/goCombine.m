% 
% /******************************************************************************
% 
% This program, "NCDM", the associated MATLAB scripts and all 
% provided data, are copyright (C) 2013-2014 Andrew R. Cohen and Paul
% M. B. Vitanyi.  All rights reserved.
% 
% This program uses bzip2 compressor as a static library.
% See the file SRC\C\bz2static\LICENSE.txt for details on that software.
% 
% This software may be referenced as:
% 
% A.R.Cohen and P.M.B. Vitanyi, "Normalized Compression Distance of Multisets 
% with Applications," IEEE Transactions on Pattern Analysis and Machine 
% Intelligence. 2014. In Press. Also arXiv:1212.5711.  
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
% 
% 1. Redistributions of source code must retain the above copyright
%    notice, this list of conditions and the following disclaimer.
% 
% 2. The origin of this software must not be misrepresented; you must 
%    not claim that you wrote the original software.  If you use this 
%    software in a product, an acknowledgment in the product 
%    documentation would be appreciated but is not required.
% 
% 3. Altered source versions must be plainly marked as such, and must
%    not be misrepresented as being the original software.
% 
% 4. The name of the author may not be used to endorse or promote 
%    products derived from this software without specific prior written 
%    permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
% OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
% GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
% WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% 
% Andrew R. Cohen acohen@coe.drexel.edu
% Paul M. B. Vitanyi Paul.Vitanyi@cwi.nl
% NCDM  version 1.0 of 13 March 2013
% NCDM  version 2.0 (release) November 2014
% 
% ******************************************************************************/

% reads the results files, and builds the classifier input
% use before calling goNN (neural net) or goDT (ensemble of random subspace
% discriminant classifiers)
%
%
%
ROOT= '..\results\';
Dataset={...
    'trainNISTn5r100.mat','testNISTn5r100.mat';
    'trainNISTn5r100v2.mat','testNISTn5r100v2.mat';
    'trainNISTn5r100v3.mat','testNISTn5r100v3.mat';
    'trainNISTn5r100v4.mat','testNISTn5r100v4.mat';
    'trainNISTn5r100v5.mat', 'testNISTn5r100v5.mat';
    'trainNISTn5r100v6.mat', 'testNISTn5r100v6.mat';
    'trainNISTn5r100v7.mat', 'testNISTn5r100v7.mat';
    'trainNISTn5r100v8.mat', 'testNISTn5r100v8.mat'
    }
NSTRIDE=1000;    
rgTrain=zeros([60e3 NSTRIDE*size(Dataset,1)]);
combineNCDTrain=[];
for iData=1:size(Dataset,1)
    load([ROOT Dataset{iData,1}]);
    combineNCDTrain=[combineNCDTrain;NCDTraining];
%     Classify=Classify(:,1:503);
%      Classify=Classify(1:22000,:);
    size(Classify)
    for i=1:size(Classify,1)
        idxDest=4+(iData-1)*NSTRIDE;
        rgTrain(Classify(i,2),idxDest:idxDest+NSTRIDE-1)=Classify(i,4:end);
        rgTrain(Classify(i,2),3)=Classify(i,3);
    end
end
%remove elements that were used as part of the multisets
for iData=1:size(Dataset,1)
    for i=1:size(rgTrain,1)
        idxDest=4+(iData-1)*NSTRIDE;
        dx = rgTrain(i,idxDest:idxDest+NSTRIDE-1);
        if any(dx)
            continue;
        end
        rgTrain(i,3)=-1; % mark for deletion
    end
end
idxRemove=find(rgTrain(:,3)==-1);
rgTrain(idxRemove,:)=[];
rgTrain=rgTrain(:,3:end);

rgTest=zeros(10e3,NSTRIDE*size(Dataset,1));
for iData=1:size(Dataset,1)
    load([ROOT Dataset{iData,2}]);
%     Classify=Classify(:,1:503);
    size(Classify)
    for i=1:size(Classify,1)
        idxDest=4+(iData-1)*NSTRIDE;
        rgTest(Classify(i,2),idxDest:idxDest+NSTRIDE-1)=Classify(i,4:end);
        rgTest(Classify(i,2),3)=Classify(i,3);
    end
end
% remove missing
for iData=1:size(Dataset,1)
    for i=1:size(rgTest,1)
        idxDest=4+(iData-1)*NSTRIDE;
        dx = rgTest(i,idxDest:idxDest+NSTRIDE-1);
        if any(dx)
            continue;
        end
        rgTest(i,3)=-1; % mark for deletion
    end
end
idxRemove=find(rgTest(:,3)==-1);
rgTest(idxRemove,:)=[];
rgTest=rgTest(:,3:end);

save('combine.mat','rgTest','rgTrain');