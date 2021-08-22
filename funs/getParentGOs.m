function [ parGOs ] = getParentGOs(childGO, selGOs, goObj)
% % %get the parentGOs of childGO(without prefix 'GO:', all input in the numeric form).
% % %  Coded by Guoxian Yu (guoxian85@gmail.com), College of Computer and Information Science,
% % %   Southwest University.
% % %   version 1.0 date:2014-02-26
%%%%%%%%%%
%   childGO:  the children GOs set of the selected GOs
%   selGOs: the selected GOs vector
%   goObj:  gene ontology read from GO file (obo format)
%%%%%%%%%%

num_child=size(childGO,1);
parGOs=cell(size(childGO,1),1);
for ii=1:num_child
    child=childGO(ii);
    parGO=getancestors(goObj,child,'Exclude',true);
    parGOs{ii}=intersect(parGO,selGOs);
end

