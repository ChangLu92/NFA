 function [ selGOsDepth ] = getSelGOsDepth(selGOs, goObj, rootGO)
%% get the depth of selGOs, depth of root node is 0
%   written by Guoxian Yu (guoxian85@gmail.com), College of Computer and Information Sciences,
%   Southwest University.
%   version 2.0 date:2014-09-12 
tic;
maxDepth=100;%currently, the maximum depth is 15 only
descendants2=[];
%determine the maximum depth
for depth=1:maxDepth
     descendants=getdescendants(goObj,rootGO,'Depth',depth, 'Exclude',true);
     temp=setdiff(descendants,descendants2); 
     if isempty(temp)
        break;
     else
         descendants2=descendants;
     end
end
maxDepth=depth-1;

selGOsDepth=zeros(length(selGOs),1);%recorded the GO depth
root_Idx=getGOIdx(rootGO, selGOs);
selGOsDepth(root_Idx)=1;
for depth=maxDepth:-1:1
    descendants=getdescendants(goObj,rootGO,'Depth',depth, 'Exclude',true);
    desc_Idxs=getGOIdx(descendants, selGOs);
    desc_Idxs(desc_Idxs==0)=[];
    selGOsDepth(desc_Idxs)=depth+1;
end
 disp(['计算节点深度运行时间：',num2str(toc)]);