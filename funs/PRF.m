function [tp,per_pre,per_re,per_f1,microprecision,microrecall,num_candidate]=PRF(gnd,Y,Z)
% % compute precision, recall and f1-measure.
% % %  Coded by Chang Lu (lucifer1992@email.swu.edu.cn), College of Computer and Information Science,
% % %   Southwest University.
% % %   version 1.0 date:2016-03-12
%%%%%%%%%%
%   gnd: the protein-function association matrix after adding noisy annotations 
%   Y: the original protein-function association matrix 
%   Z: the protein-function association matrix after dislodging noisy annotations
%%%%%%%%%%
[Ndata,~]=size(gnd); 
tp=zeros(Ndata,1);
per_pre=zeros(1,Ndata);
per_re=zeros(1,Ndata);
per_f1=zeros(1,Ndata);
num_noise=zeros(1,Ndata);
num_candidate=zeros(1,Ndata);
for ii=1:size(Y,1)
    aa=find(gnd(ii,:)-Y(ii,:)>0); %真实的噪声
    bb=find((gnd(ii,:)-Z(ii,:))>0); %被判定为的噪声
    if ~isempty(intersect(aa,bb))
       tp(ii)=length(intersect(aa,bb));
        per_pre(ii)=tp(ii)/length(bb);
        per_re(ii)=tp(ii)/length(aa);
        per_f1(ii)=2*per_pre(ii)*per_re(ii)/(per_pre(ii)+per_re(ii));
    end
     num_noise(ii)=length(aa);
    num_candidate(ii)=length(bb);
end
microprecision=sum(tp)/sum(num_candidate);
microrecall=sum(tp)/sum(num_noise);