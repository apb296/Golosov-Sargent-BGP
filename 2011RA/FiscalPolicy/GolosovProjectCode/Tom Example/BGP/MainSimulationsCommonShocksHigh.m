% This script runs a simulation on existing set of shocks with two
% different starting values btild =0; btild=1

close all
clear all
SetParaStruc
%% Set the Parallel Config
err=[];
try
    matlabpool('size')
catch err
end
if isempty(err)
    
    
    if(matlabpool('size') > 0)
        matlabpool close
    end
    
    matlabpool open local;
    
end
  

NumSim=10000;
sHist0=round(rand(10000,1))+1;
OldData=load( 'Data/SimDataParallelCommonShocks.mat');
sHist0=OldData.sHist;
sHist0 = [sHist0;2*ones(100,4)]; %4 is because sHist0 is 10000x4 matrix
NumSim = 10100;


K=4;

ex(1).casename='BM'; % benchmark
ex(2).casename='Pareto' ;% alternative pareto weights
ex(3).casename='GVol'; % mps in g
ex(4).casename='Ineq'; % higher diff in thetas
for ctrb=1:K
CoeffFileName=['Data/CompStats/c' ex(ctrb).casename '.mat'];
Sol=load(CoeffFileName);
Param(ctrb)=Sol.Para;
end

for ctrb=1:K
  CoeffFileName=['Data/CompStats/c' ex(ctrb).casename '.mat'];

[sHist(:,ctrb),gHist(:,ctrb),u2btildHist(:,ctrb),RHist(:,ctrb),...
TauHist(:,ctrb),YHist(:,ctrb),TransHist(:,ctrb),btildHist(:,ctrb),...
c1Hist(:,ctrb),c2Hist(:,ctrb),l1Hist(:,ctrb),l2Hist(:,ctrb),...
IntHist(:,ctrb),IncomeFromAssets_Agent1Hist(:,ctrb),...
AfterTaxWageIncome_Agent1Hist(:,ctrb),AfterTaxWageIncome_Agent2Hist(:,ctrb),...
GShockDiffHist(:,ctrb),TransDiffHist(:,ctrb),LaborTaxAgent1DiffHist(:,ctrb),...
LaborTaxAgent2DiffHist(:,ctrb),DebtDiffHist(:,ctrb),GiniCoeffHist(:,ctrb)]...
=RunSimulations(CoeffFileName,0,NumSim,Param(ctrb),sHist0);
end

%Added High to end of data path
save( [Para.datapath 'SimDataParallelCommonShocksHigh.mat'],'sHist',...
       'gHist','u2btildHist','RHist','TauHist','YHist','TransHist',...
       'btildHist','c1Hist','c2Hist','l1Hist','l2Hist','Para','IntHist',...
       'AfterTaxWageIncome_Agent1Hist','AfterTaxWageIncome_Agent2Hist',...
       'IncomeFromAssets_Agent1Hist','GShockDiffHist','TransDiffHist',...
       'LaborTaxAgent1DiffHist','LaborTaxAgent2DiffHist','DebtDiffHist',...
       'GiniCoeffHist')


