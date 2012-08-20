% This script runs K simulations of length NumSim with different choices of
% starting value of initial difference in assets. 

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
K=6;
btild0grid=linspace(-1.5,1.5,K);
NumSim=5000;

parfor ctrb=1:K
[sHist(:,ctrb),u2btildHist(:,ctrb),RHist(:,ctrb),TauHist(:,ctrb),YHist(:,ctrb),TransHist(:,ctrb),btildHist(:,ctrb)]=RunSimulations(2,LastIter,btild0grid(ctrb),NumSim,Para);
end
save( [Para.datapath 'SimDataParallel.mat'],'sHist','u2btildHist','RHist','TauHist','YHist','TransHist','btildHist')
% 
% figure()
% for ctrb=1:K
% subplot(3,2,ctrb)
% plot(TauHist(:,ctrb))
% end


