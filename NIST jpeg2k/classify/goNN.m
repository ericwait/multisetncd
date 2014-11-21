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

% neural net classifier
% run goCombine to get rgTrain/rgTest
function res=goNN(rgTrain,rgTest)

trainData = rgTrain(:,2:end);
trainLabel = rgTrain(:,1);

testData = rgTest(:,2:end);
testLabel = rgTest(:,1);

[size(testData) size(trainData)]



trainLabel=zeros(size(rgTrain,1),10);
for i=1:size(rgTrain,1)
    trainLabel(i,rgTrain(i,1)+1)=1;
end

res=[];
for nn=15:20 %6000
    net=patternnet([ nn]);
    %     net.trainFcn='trainrp';
    %     net.performFcn='mse';
    net.trainParam.showWindow=false;
%     net.trainParam.showCommandLine=false;
    %      [net tr]=train(net,trainData',trainLabel','useParallel','yes');
    [net tr]=train(net,trainData',trainLabel');
    
    y=net(testData');
    y=y';
    pred=[];
    for i=1:size(y,1)
        [mm idx] = max(y(i,:));
        pred(i)=idx-1;
    end
    pred=pred';
    accuracy=length(find(pred==testLabel))/length(testLabel);
    
    y=net(trainData');
    y=y';
    pred=[];
    for i=1:size(y,1)
        [mm idx] = max(y(i,:));
        pred(i)=idx-1;
    end
    pred=pred';
    accuracyTrain=length(find(pred==rgTrain(:,1)))/length(rgTrain(:,1));
    
    
    res=[res;nn accuracy accuracyTrain tr.best_vperf accuracyTrain+tr.best_vperf];
    [mm idxBest]=max(res(:,3));
    fprintf('              %d: %.4f %.4f %.4f best = %d,%.4f,%.4f\n',nn,accuracy,accuracyTrain,tr.best_vperf,res(idxBest,1),res(idxBest,2),res(idxBest,3));
end
[mm idxBest]=max(res(:,3));
res=res(idxBest,:);

