function [] = NoisyGOA_NtN( dataseth,datasetr,goObj,GOs )
%NOISYGOA_SP 此处显示有关此函数的摘要
%   此处显示详细说明
fprintf('start %s  at %s\n,==Method:%s==',dataseth, datestr(now),'NoisyGOA_NtN+EC');
load(dataseth);
load(datasetr);
% datapath=[pwd,filesep,'results3',filesep];%pwd is the current work directory
% addpath(datapath);

selGOs=GOs;
size_go=length(selGOs);
k=100;

% if size_go==3739
if size_go==3958
    gnd1=hGO.ccLabels;
    gnd2=rGO.ccLabels;
    rootGO=5575;   %ccroot
end

% if size_go==10031
if size_go==10217
    gnd1=hGO.mfLabels;
    gnd2=rGO.mfLabels;
    rootGO=3674;   %mfroot
end

% if size_go==25954
if size_go==26382
    gnd1=hGO.bpLabels;
    gnd2=rGO.bpLabels;
    rootGO=8150;  %bproot
end

gnd_h=full(gnd1);
gnd_r=full(gnd2);

%only test on annotated proteins
index=find(sum(gnd_h,2)==0);
gnd_h(index,:)=[];
gnd_r(index,:)=[];
[~,Nfun]=size(gnd_h);


minT=1;% the minimum size of member proteins
% maxT=300;
fun_stat_h=sum(gnd_h,1);
fun_stat_r=sum(gnd_r,1);
sel_funh_idx=find(fun_stat_h>=minT);
sel_funr_idx=find(fun_stat_r>=minT);
sel_fun_idx=union(sel_funh_idx,sel_funr_idx);
% sel_fun_idx=intersect(sel_funh_idx,sel_funr_idx);
selGOs=GOs(sel_fun_idx);
gndh=gnd_h(:,sel_fun_idx);
[Ndata, Nfun]=size(gndh);
rootidx=getGOIdx(rootGO,selGOs);

gnd_r=gnd_r(:,sel_fun_idx);
gnd_h=gnd_h(:,sel_fun_idx);

num_perprotein_noise=zeros(Ndata,1); %the number of noisy annotations of each protein
gnd=gnd_r-gnd_h;
sub_goObj=getSelGoObj(selGOs,goObj);%filter the goObj to speedup computation
%  DirectchildGOs=getDirectChildGOs(selGOs, selGOs,sub_goObj);
childGOs=getChildGOs(selGOs, selGOs,sub_goObj);
%  DirectparGOs=getDirectParentGOs(selGOs, selGOs,sub_goObj);
Depth  = getSelGOsDepth(selGOs,sub_goObj, rootGO);
% parGOs=getParentGOs(selGOs, selGOs,sub_goObj);

%%  count the number of noisy annotations
for ii=1:Ndata
    idx=find(gnd(ii,:)==-1);
    noiseidx{ii}=idx;
    num_perprotein_noise(ii)=length(idx);
end

%% Identifying Noisy Gene Ontology Annotations
[gnd_hw] = NtN(gnd_h, selGOs, childGOs, Depth ) ;  %加权
k=round(Nfun*0.5);
[U,S,V] = svd(gnd_hw); %返回一个与X 同大小的对角矩阵S，两个酉矩阵U 和V，
Sk=S(1:k,1:k);
Uk=U(:,1:k);
Vk=V(1:k,:);
newgnd=Uk*Sk*Vk;
  

%% 直接用稀疏表示

Idx=find(gnd_h==0);
newgnd(Idx)=0;

% dt=sum(newgnd,2);
% idx_dt = (dt>eps);
% dt(idx_dt) = 1./dt(idx_dt); 
% newgnd=diag(dt)*newgnd;

%  newgnd2=zeros(Ndata,Nfun);
%  for ii=1:Ndata
%      newgnd2(ii,:)= TPR_DAG(newgnd(ii,:),gnd_h(ii,:),DirectchildGOs,DirectparGOs,selGOs,rootidx);
%  end
% % % 
%   newgnd=newgnd2;
newgnd=newgnd+gnd_h;
newgnd=1./newgnd;
newgnd(newgnd>=inf&newgnd<=inf)=0;
 [value,ind]=sort(newgnd,2,'descend');

 
newgnd2=gnd_h;
for ii=1:Ndata
    Idx=find(newgnd2(ii,:)>0);
    noise=ind(ii,1:num_perprotein_noise(ii));
    childnoise=childGOs(noise);
    newgnd2(ii,noise)=0;
    for jj=1:length(childnoise)
        childidxp=getGOIdx(childnoise{jj},selGOs);
        childIdx=intersect(childidxp,Idx);
        newgnd2(ii,childIdx)=0;
    end
end

%% compute precision, recall and f1-measure
Y=gnd_r;
Z =newgnd2(1:Ndata,:);


% ratio=0.8;
% time=500;
% data=length(find(num_perprotein_noise>0));
% selNdata=round(Ndata*ratio);
% per_truePositive=zeros(time,selNdata);
% per_precision=zeros(time,selNdata);
% per_recall=zeros(time,selNdata);
% per_F1score=zeros(time,selNdata);
% for run=1:time
%     sortprotein=randperm(Ndata);
%     protein_idx=sortprotein(1:selNdata);
%     rangnd_h=gnd_h(protein_idx,:);
%     ranY=Y(protein_idx,:);
%     ranZ=Z(protein_idx,:);
%     [tp,per_pre,per_re,per_f1]=PRF(rangnd_h,ranY,ranZ);
%     per_truePositive(run,:)=tp;
%     per_precision(run,:)=per_pre;
%     per_recall(run,:)=per_re;  
%     per_F1score(run,:)=per_f1;
% end
% ave_truePositive=sum(per_precision,2)/data;  
% ave_precision=sum(per_precision,2)/data;  
% ave_recall=sum(per_precision,2)/data;
% ave_fvalue=sum(per_F1score,2)/data;
% 
% truePositives=sum(ave_truePositive)/time;  
% precisions=sum(ave_precision)/time;  
% recalls=sum(ave_recall)/time;
% fvalue=sum(ave_fvalue)/time;

[tp,per_pre,per_re,per_f1,Miprecisions,Mirecall,num_candidate]=PRF(gnd_h,Y,Z);
data=length(find(num_perprotein_noise>0));
[maprecisions,marecalls,mafvalue,miprecisions,mirecalls,mifvalue,ave_maprecision,ave_marecall,ave_mafvalue,ave_miprecision,ave_mirecall,ave_mifvalue]=bootstrapping(tp,per_pre,per_re,per_f1,num_candidate,num_perprotein_noise);
prec_seq='tp,per_pre,per_re,per_f1,Macropre,Macrore,Macrof1,Maprecisions,Marecall,Maf1,Miprecisions,Mirecall，Mif1,maprecisions,marecalls,mafvalue,miprecisions,mirecalls,mifvalue,ave_maprecision,ave_marecall,ave_mafvalue,ave_miprecision,ave_mirecall,ave_mifvalue,num_perprotein_noise,num_candidate';

precision=cell(30,1);
precision{1}=tp;
precision{2}=per_pre;
precision{3}=per_re;
precision{4}=per_f1;
precision{5}=sum(per_pre)/data;
precision{6}=sum(per_re)/data;
precision{7}=sum(per_f1)/data;
precision{8}=Miprecisions;
precision{9}=Mirecall;
precision{10}=2*Miprecisions*Mirecall/(Miprecisions+Mirecall);
precision{11}=maprecisions;
precision{12}=marecalls;
precision{13}=mafvalue;
precision{14}=miprecisions;
precision{15}=mirecalls;
precision{16}=mifvalue;
precision{17}=ave_maprecision;
precision{18}=ave_marecall;
precision{19}=ave_mafvalue;
precision{20}=ave_miprecision;
precision{21}=ave_mirecall;
precision{22}=ave_mifvalue;
precision{23}=num_perprotein_noise;
precision{24}=num_candidate;


% 
stds=cell(30,1);
stds{11}=std(ave_maprecision,0,1);
stds{12}=std(ave_marecall,0,1);
stds{13}=std(ave_mafvalue,0,1);
stds{14}=std(ave_miprecision,0,1);
stds{15}=std(ave_mirecall,0,1);
stds{16}=std(ave_mifvalue,0,1);

if rootGO==5575
    evalstr=['save results',filesep,dataseth, '_NoisyGOA_NtN_hr_cc.mat precision stds prec_seq'];
end
if rootGO==3674
    evalstr=['save results',filesep,dataseth, '_NoisyGOA_NtN_hr_mf.mat precision stds prec_seq'];
end
if rootGO==8150
    evalstr=['save results',filesep,dataseth, '_NoisyGOA_NtN_hr_bp.mat precision stds prec_seq'];
end
eval(evalstr);


fprintf('\n =====finish NoisyGOA_SP_hr time=%s\n',datestr(now));


end

