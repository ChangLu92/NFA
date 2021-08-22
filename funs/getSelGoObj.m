function [ sub_goObj ] = getSelGoObj(selGOs,goObj) %从所有goobj中选取selGOs个项
% % % filter the goObj to speedup computation
% % %  Coded by Guoxian Yu (guoxian85@gmail.com), College of Computer and Information Science,
% % %   Southwest University.
% % %   version 1.0 date:2014-02-26
%%%%%%%%%%
%   goObj:  gene ontology read from GO file (obo format)
%   selGOs: the selected GOs vector
%%%%%%%%%%
terms_ids=get(goObj.terms,'id');
num_id=length(terms_ids);
ids=zeros(num_id,1);
for ii=1:num_id
    temp=terms_ids{ii};
    ids(ii)=temp;
end

mask=zeros(num_id,1);
for ii=1:length(selGOs)
    selGO=selGOs(ii);
    idx=find(ids==selGO);
    if(idx>0)
        mask(idx)=1;
    end
end
mask=logical(mask);
sel_terms=goObj.terms(mask);
sub_goObj=goObj(sel_terms);
sel_terms=get(sub_goObj.terms,'id');
fprintf('Selected %d GOs \n', length(sel_terms));
