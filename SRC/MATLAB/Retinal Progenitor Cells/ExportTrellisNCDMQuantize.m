% 
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
path(path,'..\')

% export a trellis for NCDM

% 'Trellis DCC59 54 43 dxy tot_dist phi_dist e s.mat' these are the 86
% cells that are used for the terminal differentiation prediction. idxCell
% is the cell fate.
% idxTrue is used for crossValidation and multiples classes
% load 'Trellis DCC59 54 43 dxy tot_dist phi_dist e s.mat'
%outcome=[1 ph,ph;2 bi,bi;3 bi,ph;4 ph,am;5 am,bi;6 am,am;7 pro,*]

% % this code does a 3 class problem, ph,ph vs. bi,ph vs. ph,am
% %     no bi,bi or am,bi or am,am
% for i=1:length(Trellis)
%     if 1==Trellis(i).idxCell
%         Trellis(i).idxTrue = 1;
%     elseif 3==Trellis(i).idxCell
%         Trellis(i).idxTrue = 2;
%     elseif 4==Trellis(i).idxCell
%         Trellis(i).idxTrue = 3;
%     else
%         Trellis(i).idxTrue = -1;
%     end
% end
% % remove the non-analyzed cell types
% nidxNoAnalyze=find(-1==[Trellis.idxTrue]);
% Trellis(nidxNoAnalyze)=[];

load 'Trellis SR2 progenitors 54 55 59 dxy tot_dist phi_dist e s.mat'

% this code is for the 2 class ph,ph vs. everyone else question
for i=1:length(Trellis)
    if 7==Trellis(i).idxCell
        Trellis(i).idxTrue = 1;
    else
        Trellis(i).idxTrue = 2;
    end
end


for nF=3:63 % start at 3 since [dx dy] always on together
    
    strBin=dec2base(nF,2);
    nBin=[];
    for i=1:length(strBin)
        nBin(length(strBin)-i+1)=str2num(strBin(i));
    end
    % dx,dy always off or on together 
    if xor(nBin(1),nBin(2)),continue, end
    
    FS = find(nBin);
    for N=2:26
        xy=[];
        for i=1:length(Trellis)
            xy=[xy;Trellis(i).Features(:,FS)];
        end
        
        mu=mean(xy);
        sig=cov(xy);
        s1=sig^-1;
        Breakpoints=[1/N:1/N:1-1/N];
        c2 = chi2inv(Breakpoints,size(Trellis(i).Features,2));
        for i= 1: length(Trellis)
            Trellis(i).sData=esym2d(N,c2,mu,s1,Trellis(i).Features(:,FS));
        end
        
        
        fname = ['SRQ\SR_nF_' num2str(nF) 'N_' num2str(N) '.txt'];
        fout = fopen(fname,'w');
        
        BufferSize = 0;
        for i=1:length(Trellis)
            V1 = Trellis(i).sData();
            if iscell(V1)
                V1=cell2mat(V1);
            end
            vv= mat2str(V1);
            vv=vv(2:end-1);
            sz=length(vv(:))+1;% +1 for \n
            Trellis(i).vv=vv;
            Trellis(i).sz=sz;
            BufferSize=BufferSize+sz;
        end
        fprintf(fout,'%d,%d\n',length(Trellis),BufferSize);
        for i=1:length(Trellis)
            fprintf(fout,'%d\n',Trellis(i).idxTrue);
            fprintf(fout,'%d\n',Trellis(i).sz);
            fprintf(fout,'%s\n',Trellis(i).vv);
        end
        fclose(fout);
        
    end
end % N