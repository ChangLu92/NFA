function [] = NoGOA_SP( dataseth,datasetr,goObj,GOs )
%   此处显示详细说明
fprintf('start %s  at %s\n,==Method:%s==',dataseth, datestr(now),'NoisyGOA_SP');
load(dataseth);
load(datasetr);

selGOs=GOs;
size_go=length(selGOs);
lambda=0.5;
k=5;

if size_go==3958
    gnd1=hGO.ccLabels;
    gnd2=rGO.ccLabels;
    rootGO=5575;   %ccroot
end

if size_go==10217
    gnd1=hGO.mfLabels;
    gnd2=rGO.mfLabels;
    rootGO=3674;   %mfroot
end

if size_go==26382
    gnd1=hGO.bpLabels;
    gnd2=rGO.bpLabels;
    rootGO=8150;  %bproot
end


gnd_h=gnd1;
gnd_r=gnd2;

%only test on annotated proteins
index=find(sum(gnd_h,2)==0);
gnd_h(index,:)=[];
gnd_r(index,:)=[];
proteins(index)=[];

minT=1;% the minimum size of member proteins
fun_stat_h=sum(gnd_h,1);
fun_stat_r=sum(gnd_r,1);
sel_funh_idx=find(fun_stat_h>=minT);
sel_funr_idx=find(fun_stat_r>=minT);
sel_fun_idx=union(sel_funh_idx,sel_funr_idx);
% sel_fun_idx=intersect(sel_funh_idx,sel_funr_idx);
selGOs=GOs(sel_fun_idx);
gndh=gnd_h(:,sel_fun_idx);
[Ndata, Nfun]=size(gndh);

gnd_r=gnd_r(:,sel_fun_idx);
gnd_h=gnd_h(:,sel_fun_idx);

rootidx=getGOIdx(rootGO,selGOs);
num_perprotein_noise=zeros(Ndata,1); %the number of noisy annotations of each protein
gnd=gnd_r-gnd_h;


childGOs=getChildGOs(selGOs, selGOs,goObj);
% DirectchildGOs=getDirectChildGOs(selGOs, selGOs,goObj);
% DirectparGOs=getDirectParentGOs(selGOs, selGOs,goObj);


% %% 50%训练，50%测试
% random=randperm(Ndata);
% Ndata_train=round(Ndata*0.5);
% % gnd_hec_train=gnd_hec(random(1:Ndata_train),:);
% gnd_train=gnd(random(1:Ndata_train),:);
% 
% gnd_r(random(1:Ndata_train),:)=[];
% gnd_h(random(1:Ndata_train),:)=[];
% % gnd_hec(random(1:Ndata_train),:)=[];
% gnd(random(1:Ndata_train),:)=[];
% 
% 
% 
% %% 测试
% num_perprotein_noise=zeros(Ndata-Ndata_train,1); 
% for ii=1:(Ndata-Ndata_train)
%     idx=find(gnd(ii,:)==-1);
%     noiseidx{ii}=idx;
% %     idxec=find(gnd_hec(ii,:)>0);
% %     idxinter=intersect(idx,idxec);
% %     gnd_ec(ii,idxinter)=gnd_hec(ii,idxinter);
%     num_perprotein_noise(ii)=length(idx);
% end
% Ndata=Ndata-Ndata_train;


%%  count the number of noisy annotations
for ii=1:Ndata
    idx=find(gnd(ii,:)==-1);
    noiseidx{ii}=idx;
    num_perprotein_noise(ii)=length(idx);
end

%% Identifying Noisy Gene Ontology Annotations
% [gnd_hw] = TFIDF(gnd_h, selGOs, Ndata,childGOs);  %加权
  S=SparsityW(gnd_h',lambda);%权重的求解
  S1=(S+S')/2;
 newgnd=S1*gnd_h;
%  newgnd=newgnd';
%  newgnd2=zeros(Ndata,Nfun);
%  for ii=1:Ndata
%      newgnd2(ii,:)= TPR_DAG(newgnd(ii,:),DirectchildGOs,DirectparGOs,selGOs,rootidx);
%  end
%% 直接用稀疏表示
%  newgnd=newgnd';
Idx=find(gnd_h==0);
newgnd(Idx)=0;
newgnd=newgnd+gnd_h;
newgnd=1./newgnd;
newgnd(newgnd>=inf&newgnd<=inf)=0;
 [~,ind]=sort(newgnd,2,'descend');
 
%   if rootGO==5575
%     evalstr=['save -v7.3 resultsyu',filesep,dataseth, '_SP_noise_cc.mat ind noiseidx selGOs S proteins'];
% end
% if rootGO==3674
%     evalstr=['save -v7.3 resultsyu',filesep,dataseth, '_SP_noise_mf.mat ind noiseidx selGOs S proteins'];
% end
% if rootGO==8150
%     evalstr=['save -v7.3 resultsyu',filesep,dataseth, '_SP_noise_bp.mat ind noiseidx selGOs S proteins'];
% end
% eval(evalstr);
 
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
    evalstr=['save results',filesep,dataseth, '_NoisyGOA_SP_hr_cc.mat precision prec_seq stds '];
end
if rootGO==3674
    evalstr=['save results',filesep,dataseth, '_NoisyGOA_SP_hr_mf.mat precision prec_seq stds '];
end
if rootGO==8150
    evalstr=['save results',filesep,dataseth, '_NoisyGOA_SP_hr_bp.mat precision prec_seq stds '];
end
eval(evalstr);

fprintf('\n =====finish NoisyGOA_SP_hr time=%s\n',datestr(now));


end

