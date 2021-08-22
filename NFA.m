function [] = NFA(dataseth,datasetr,goObj,GOs,ECw )
%predicting noisy GO annotations using evidences and sparse representation
% Chang Lu College of Computer and Information
% Science, Southwest University. Contact gxyu@swu.edu.cn, lucifer1992@email.swu.edu.cn
%%%%%%%%%%
%   goObj: gene ontology read from GO file (obo format)
%   GOs:  GO term id for each column of the protein-function association matrix
%   dataseth: Name of the dataset archived on 2015-11-09.
%   datasetr: Name of the dataset archived on 2016-04-11.
%%%%%%%%%%
fprintf('start %s  at %s\n,==Method:%s==',dataseth, datestr(now),'NFA');
load(dataseth);
load(datasetr);

selGOs=GOs;
size_go=length(selGOs);
lambda=0.5;

if size_go==3958
    gnd_h=hGO.ccLabels;
    gnd_r=rGO.ccLabels;
    gnd_hec=hGO.ccECs;
    rootGO=5575;   %ccroot
end

if size_go==10217
    gnd_h=hGO.mfLabels;
    gnd_r=rGO.mfLabels;
    gnd_hec=hGO.mfECs;
    rootGO=3674;   %mfroot
end

if size_go==26382
    gnd_h=hGO.bpLabels;
    gnd_r=rGO.bpLabels;
    gnd_hec=hGO.bpECs;
    rootGO=8150;  %bproot
end

%only test on annotated proteins
index=find(sum(gnd_h,2)==0);
gnd_h(index,:)=[];
gnd_r(index,:)=[];
gnd_hec(index,:)=[];

minT=1;% the minimum size of nember proteins
fun_stat_h=sum(gnd_h,1);
fun_stat_r=sum(gnd_r,1);
sel_funh_idx=find(fun_stat_h>=minT);
sel_funr_idx=find(fun_stat_r>=minT);
sel_fun_idx=union(sel_funh_idx,sel_funr_idx);

selGOs=GOs(sel_fun_idx);
gndh=gnd_h(:,sel_fun_idx);
[Ndata, Nfun]=size(gndh);

gnd_r=gnd_r(:,sel_fun_idx);
gnd_h=gnd_h(:,sel_fun_idx);
gnd_hec=gnd_hec(:,sel_fun_idx);

num_perprotein_noise=zeros(Ndata,1); %the number of noisy annotations of each protein
gnd=gnd_r-gnd_h;
sub_goObj=getSelGoObj(selGOs,goObj);%filter the goObj to speedup computation
childGOs=getChildGOs(selGOs, selGOs,sub_goObj);
parGOs=getParentGOs(selGOs, selGOs,sub_goObj);

%%  calculating the number of noisy annotations
for ii=1:Ndata
    idx=find(gnd(ii,:)==-1);
    num_perprotein_noise(ii)=length(idx);
end

%% Identifying Noisy Gene Ontology Annotations
% [gnd_hw] = EvidenceCode(gnd_h, gnd_hec, selGOs, Ndata, parGOs);
[gnd_hw] = EvidenceCode1(gnd_h, gnd_hec, selGOs, Ndata, parGOs, ECw);
S=SparsityW(gnd_hw',lambda);  %Sparse representation
S=(S+S')/2;  
V=S * gnd_hw;
idx_h=find(gnd_h==0);
V(idx_h)=0;
V=V+gnd_h;
V=1./V;
V(V>=inf&V<=inf)=0;
[val_V,idx_V]=sort(V,2,'descend'); 
newgnd=gnd_h;
for ii=1:Ndata
    noise_idx=idx_V(ii,1:num_perprotein_noise(ii));
    childnoise=childGOs(noise_idx);
    newgnd(ii,noise_idx)=0;
    for jj=1:length(childnoise)
        child_idx=getGOIdx(childnoise{jj},selGOs);
        newgnd(ii,child_idx)=0;
    end
end

%% compute MacroP, MacroR and MacroF; MicroP, MicroR and MicroF;
Y=gnd_r;
Z =newgnd(1:Ndata,:);
[tp,per_pre,per_re,per_f1,Miprecisions,Mirecall,num_candidate]=PRF(gnd_h,Y,Z);
num_noise=length(find(num_perprotein_noise>0));
[maprecisions,marecalls,mafvalue,miprecisions,mirecalls,mifvalue,ave_maprecision,ave_marecall,ave_mafvalue,ave_miprecision,ave_mirecall,ave_mifvalue]=bootstrapping(tp,per_pre,per_re,per_f1,num_candidate,num_perprotein_noise);
prec_seq='tp,per_pre,per_re,per_f1,Macropre,Macrore,Macrof1,Maprecisions,Marecall,Maf1,Miprecisions,Mirecall£¬Mif1,maprecisions,marecalls,mafvalue,miprecisions,mirecalls,mifvalue,ave_maprecision,ave_marecall,ave_mafvalue,ave_miprecision,ave_mirecall,ave_mifvalue,num_perprotein_noise,num_candidate';

precision=cell(30,1);
precision{1}=tp;
precision{2}=per_pre;
precision{3}=per_re;
precision{4}=per_f1;
precision{5}=sum(per_pre)/num_noise;
precision{6}=sum(per_re)/num_noise;
precision{7}=sum(per_f1)/num_noise;
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

stds=cell(30,1);
stds{11}=std(ave_maprecision,0,1);
stds{12}=std(ave_marecall,0,1);
stds{13}=std(ave_mafvalue,0,1);
stds{14}=std(ave_miprecision,0,1);
stds{15}=std(ave_mirecall,0,1);
stds{16}=std(ave_mifvalue,0,1);

if rootGO==5575
    evalstr=['save results',filesep,dataseth, '_NFA_cc.mat precision stds prec_seq'];
end
if rootGO==3674
    evalstr=['save results',filesep,dataseth, '_NFA_mf.mat precision stds prec_seq'];
end
if rootGO==8150
    evalstr=['save results',filesep,dataseth, '_NFA_bp.mat precision stds prec_seq'];
end
eval(evalstr);
fprintf('\n =====finish NFA time=%s\n',datestr(now));
end

