 function [ childGOs ] = getChildGOs(parGO, selGOs, goObj)
% % %get the parentGOs of childGO(without prefix 'GO:', all input in the numeric form).
% % %  Coded by Guoxian Yu (guoxian85@gmail.com), College of Computer and Information Science,
% % %   Southwest University.
% % %   version 1.0 date:2014-02-26
%%%%%%%%%%
%   parGO:  the parents GOs set of the selected GOs
%   selGOs: the selected GOs vector
%   goObj:  gene ontology read from GO file (obo format)
%%%%%%%%%%
num_child=size(parGO,1);
childGOs=cell(num_child,1);
for ii=1:num_child
    go=parGO(ii);
    childGO=getdescendants(goObj,go,'Exclude',true);
    childGOs{ii}=intersect(childGO,selGOs);
end