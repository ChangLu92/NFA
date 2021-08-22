function [gnd_w] = EvidenceCode(gnd, gnd_ec, selGOs, Ndata, parentG0s)
% Giving weights to gene-term matrix using evidence codes
% 1	'TAS'									
% 0.8	'EXP'	'IDA'								
% 0.6	'IPI'	'IMP'	'IGI'	'IEP'						
% 0.4	'ISA'	'ISO'	'ISM'	'RCA'	'IGC' 'IBA'	'IBD'	'IKR'	'IRD'	'ISS'
% 0.2	'NAS'	'IC'								
% 0.1	'ND'	'IEA'
% % %  Coded by Chang Lu (lucifer1992@email.swu.edu.cn), College of Computer and Information Science,
% % %   Southwest University.
% % %   version 1.0 date:2016-08-07
%%%%%%%%%%
%   gnd:the protein-function association matrix
%   gnd_ec: evidence codes the protein-function association matrix
%   selGOs: the selected GOs vector
%   Ndata:  the number of genes
%   parentG0s: the parents GOs set of the selected GOs
%%%%%%%%%%
ECw=[1 0.8 0.8 0.6 0.6 0.6 0.6 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4];
gnd_w=gnd;
for ii=1:Ndata
    idxleaf=find(gnd_ec(ii,:)>0);
    evidences=gnd_ec(ii,idxleaf);
    [value, eviidx]=sort(evidences,2,'descend');
    idxleaf=idxleaf(eviidx);
    parleaf=parentG0s(idxleaf);
    for jj=1:length(idxleaf)
        parjjidx=getGOIdx(parleaf{jj},selGOs);  %parGOs
        gnd_w(ii,idxleaf(jj))=ECw(value(jj));
        gnd_w(ii,parjjidx)=ECw(value(jj));
    end
end
