datapath=[pwd,filesep,'data',filesep];%pwd is the current work directory
addpath(datapath);                    
datapath=[pwd,filesep,'funs',filesep];%pwd is the current work directory
addpath(datapath);

load('GeneOntointer2016102920170331.mat'); %GO file (intersection of historical GO and recent GO)
datasets={'YeastGOAinternotH20161031';'YeastGOAinternotR20170410'}; %historical GOA, archieved date 2015-11-09 %recent GOA, archieved date 2016-04-11
datasets2={'ArabidopsisGOAinternotH20161031';'ArabidopsisGOAinternotR20170410'};
% Evidences={'TAS';'EXP';'IDA';'IPI';'IMP';'IGI';'IEP';'ISA';'ISO';'ISM';'RCA';'IGC';'IBA';'IBD';'IKR';'IRD';'ISS';'NAS';'IC';'IEA';'ND'};%GOA
         ECw=[1    1      1     1    1     1    1   0.6   0.6    0.6  0.6   0.6   0.6   0.6    0.6   0.6  0.6   0.8   0.6  0.6 0.4;
              0.8    1      1     1    1     1    0.6   0.4   0.4    0.4  0.6   0.6   0.6   0.6    0.6   0.6  0.4   0.4   0.6  0.4    0; 
              1     1      1     1     1     1     1    0.4   0.4    0.4  0.6   0.6   0.6   0.6    0.6   0.6  0.4   0.8   0.6  0.6  0.4];
GOs={ccGOs;mfGOs;bpGOs};
for ii=1:3
%    NoGOA(datasets{1},datasets{2},goObj_h,GOs{ii});
  NoisyGOA_NtN(datasets{1},datasets{2},goObj_h,GOs{ii});
  NoisyGOA_NtN(datasets2{1},datasets2{2},goObj_h,GOs{ii});
end


