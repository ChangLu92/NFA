function [ goIdx ] = getGOIdx(go, GOs)
%  get the GO's the index (according to the integer label) in the selGOs,
%  all input in the numeric
%  Coded by Guoxian Yu (guoxian85@gmail.com), College of Computer and Information Science,
%   Southwest University.
%   version 1.0 date:2014-02-26
%%%%%%%%%%
%   go: GO term id
%   GOs: the selected GOs vector
%%%%%%%%%%
goIdx=zeros(length(go),1);
for ii=1:length(go)
    goii=go(ii);
%     jj=1;
    [idx]=find(GOs==goii,1); %返回GOs中goii所在位置
    if(~isempty(idx))
        goIdx(ii)=idx;
    end
%     while jj<=length(GOs) && ~bFind
%         if GOs(jj)==goii
%             goIdx(ii)=jj;%the index of GONums for the GO
%             bFind=1;
%         end
%         jj=jj+1;
%     end
end
