function [gnd_w] = EvidenceCode1(gnd, gnd_ec, selGOs, Ndata, parentG0s, ECw)
% evidence code加权
% 1	'TAS'									
% 0.8	'EXP'	'IDA'								
% 0.6	'IPI'	'IMP'	'IGI'	'IEP'						
% 0.4	'ISA'	'ISO'	'ISM'	'RCA'	'IGC' 'IBA'	'IBD'	'IKR'	'IRD'	'ISS'
% 0.2	'NAS'	'IC'								
% 0.1	'ND'	'IEA'
% Evidences={'TAS';'EXP';'IDA';'IPI';'IMP';'IGI';'IEP';'ISA';'ISO';'ISM';'RCA';'IGC';'IBA';'IBD';'IKR';'IRD';'ISS';'NAS';'IC';'IEA'};%GOA
%          ECw=[0.8    1      1     1    1     1    0.6   0.4   0.4    0.4  0.6   0.6   0.6   0.6    0.6   0.6  0.4   0.4   0.8  0.4];  % list1原始EC
%          ECw=[0.8    1      1     1    1     1    0.6   0.4   0.4    0.4  0.6   0.6   0.6   0.6    0.6   0.6  0.4   0.4   0.8  0.6];  % list2
%           ECw=[0.8    1      1     1    1     1    1   0.4   0.4    0.4  0.6   0.6   0.6   0.6    0.6   0.6  0.4   0.4   0.8  0.6];  % list3
% ECw=[1 0.8 0.8 0.8 0.8 0.8 0.8 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.5 0.5 0 0.4]; %intellgo2
% ECw=[1 0.8 0.8 0.8 0.8 0.8 0.8 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4];
%           ECw=[0.6  0.8    1   0.8    0.8   0.8   0.8   0.4    0.4   0.4  0.4   0.4   0.4    0.4   0.4   0.4  0.4   0.4  0.6   0.4];  % list6
%            ECw=[0.5    1    1      1     1     1     1    0.5   0.5    0.5  0.5   0.5   0.5    0.5   0.5   0.5  0.5   0.5   0.5   1];  % list11 
% ECw=[0.7	0.5	0.8	0.6	0.8	0.7	0.9	0.9	1	0.5	0.4	1	0.4	1	1	1	0.6	0.8	1	0.6];  % list12 
% ECw=[0.8	0.7	0.8	0.8	0.8	0.8	0.9	0.9	1	0.7	0.6	1	0.6	1	1	1	0.8	0.8	1	0.8];  % list13 
%  ECw=[0.024691358	0.5	0.030788177	0.689614593	0.021604938	0.051136364	0	0.012345679	0	0	0	0	0.156709109	0	0	0	0.029498525	0.041666667	0	0.033596838]; %list12mf
% A=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
% ECw=A-ECw;
gnd_w=gnd;
for ii=1:Ndata
    idxleaf=find(gnd_ec(ii,:)>0);
    evidences=gnd_ec(ii,idxleaf);
%     [value, eviidx]=sort(evidences,2,'descend');
%     idxleaf=idxleaf(eviidx);
    [ECwvalue, ECwidx]=sort(ECw(evidences),2,'ascend');
    idxleaf=idxleaf(ECwidx);
    parleaf=parentG0s(idxleaf);
%     [ECwvalue, ECwidx]=sort(ECw,2,'descend');
    for jj=1:length(idxleaf)
        parjjidx=getGOIdx(parleaf{jj},selGOs);  %parGOs
        gnd_w(ii,idxleaf(jj))=ECwvalue(jj);
        gnd_w(ii,parjjidx)=ECwvalue(jj);
    end
end
