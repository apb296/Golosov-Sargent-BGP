% This script runs a simulation on existing set of shocks with two
% different starting values btild =0; btild=1

close all
clear all
%% Set the Parallel Config
err=[];
try
    matlabpool
catch err
end
if isempty(err)
    
    
    if(matlabpool('size') > 0)
        matlabpool close
    end
    
    matlabpool open local;
    
end
  

% LOAD THE COEFF
LastIter=250;
load(['Data/c' num2str(LastIter) '.mat'])
SetParaStruc
%
% LOAD THE history of shocks
SimData=load(['Data/SimDataParallel.mat']);
sHist0=SimData.sHist;

K=5;

btild0grid=linspace(0,1.5,K);
NumSim=10000;

parfor ctrb=1:K
[sHist(:,ctrb),gHist(:,ctrb),u2btildHist(:,ctrb),RHist(:,ctrb),TauHist(:,ctrb),YHist(:,ctrb),TransHist(:,ctrb),btildHist(:,ctrb),c1Hist(:,ctrb),c2Hist(:,ctrb),l1Hist(:,ctrb),l2Hist(:,ctrb),IntHist(:,ctrb),IncomeFromAssets_Agent1Hist(:,ctrb),AfterTaxWageIncome_Agent1Hist(:,ctrb),AfterTaxWageIncome_Agent2Hist(:,ctrb)]=RunSimulations(LastIter,btild0grid(ctrb),NumSim,Para,sHist0);
end

save( [Para.datapath 'SimDataParallelCommonShocks.mat'],'sHist','gHist','u2btildHist','RHist','TauHist','YHist','TransHist','btildHist','btild0grid','c1Hist','c2Hist','l1Hist','l2Hist','Para','IntHist','AfterTaxWageIncome_Agent1Hist','AfterTaxWageIncome_Agent2Hist','IncomeFromAssets_Agent1Hist')


