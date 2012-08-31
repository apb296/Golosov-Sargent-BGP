 clc
 clear all
 SetParaStruc;
 
 
% -------- EXPERIMENT  0 -----------------------------------------------
% This is the baseline settings
 Para.datapath=['Data/CompStats/'];
 mkdir(Para.datapath)
 casename='BM';
 Para.StoreFileName=['c' casename '.mat'];
 CoeffFileName=[Para.datapath Para.StoreFileName];

 % ---------------CHANGE THE DEAFULT PARAMTERS----------------------------
 % This is the BM case - no change
 
 %  --- SOLVE THE BELLMAN EQUATION --------------------------------------
 Main(Para)
 
 %  --- SIMULATE THE MODEL ----------------------------------------------
 
[sHist,gHist,u2btildHist,RHist,TauHist,YHist,TransHist,btildHist,...
c1Hist,c2Hist,l1Hist,l2Hist,IntHist,IncomeFromAssets_Agent1Hist,...
AfterTaxWageIncome_Agent1Hist,AfterTaxWageIncome_Agent2Hist]...
=RunSimulations(CoeffFileNane,Para.btild0,Para.NumSim,Para);
save( [Para.datapath 'SimData' casename '.mat'],'sHist','gHist','u2btildHist',...
'RHist','TauHist','YHist','TransHist','btildHist','btild0grid',...
'c1Hist','c2Hist','l1Hist','l2Hist','Para','IntHist',...
'AfterTaxWageIncome_Agent1Hist','AfterTaxWageIncome_Agent2Hist',...
'IncomeFromAssets_Agent1Hist');
 
 
% -------- EXPERIMENT  1 -----------------------------------------------
% This setting increases the Pareto weight of Agent 2

% This is the baseline settings
 Para.datapath=['Data/CompStats/'];
 mkdir(Para.datapath)
 casename='Pareto';
 Para.StoreFileName=['c' casename '.mat'];
 CoeffFileName=[Para.datapath Para.StoreFileName];
 
 
 % ---------------CHANGE THE DEAFULT PARAMTERS----------------------------
alpha_2=.75;
alpha_1=1-alpha_1;
Para.alpha_1=alpha_1*Para.n1;
Para.alpha_2=alpha_2*Para.n2;

 
 %  --- SOLVE THE BELLMAN EQUATION --------------------------------------
 Main(Para)
 
 
 % -------- EXPERIMENT  2 -----------------------------------------------
% This setting introduces a mean preserving spread in g
 
 clear all
% This is the baseline settings
 Para.datapath=['Data/CompStats/'];
 mkdir(Para.datapath)
 Para.StoreFileName=['cPareto.mat'];

 % ---------------CHANGE THE DEAFULT PARAMTERS----------------------------
 % MPS in g . We double the spread between the expenditure shocks
MeanExpenditureShocks=mean(Para.g);
ExpenditureShocksGap=Para.g(2)-Para.g(1);
NewExpenditureshocksGap= ExpenditureShocksGap*2;
Para.g(1)=MeanExpenditureShocks-NewExpenditureshocksGap/2;
Para.g(2)=MeanExpenditureShocks+NewExpenditureshocksGap/2;
 %  --- SOLVE THE BELLMAN EQUATION --------------------------------------
 Main(Para)

  %  --- SIMULATE THE MODEL ----------------------------------------------
 
[sHist,gHist,u2btildHist,RHist,TauHist,YHist,TransHist,btildHist,...
c1Hist,c2Hist,l1Hist,l2Hist,IntHist,IncomeFromAssets_Agent1Hist,...
AfterTaxWageIncome_Agent1Hist,AfterTaxWageIncome_Agent2Hist]=...
RunSimulations(CoeffFileNane,Para.btild0,Para.NumSim,Para);
save( [Para.datapath 'SimData' casename '.mat'],'sHist','gHist','u2btildHist',...
'RHist','TauHist','YHist','TransHist','btildHist','btild0grid',...
'c1Hist','c2Hist','l1Hist','l2Hist','Para','IntHist',...
'AfterTaxWageIncome_Agent1Hist','AfterTaxWageIncome_Agent2Hist',...
'IncomeFromAssets_Agent1Hist');
 
 
 
% -------- EXPERIMENT  3 -----------------------------------------------
% This setting introduces a mean preserving spread in theta
clear all
% This is the baseline settings
 Para.datapath=['Data/CompStats/'];
 mkdir(Para.datapath)
 Para.StoreFileName=['cIneq.mat'];

 % ---------------CHANGE THE DEAFULT PARAMTERS----------------------------
% increase difference between theta1 and theta 2 (were we again pick new 
%thetas so that "first best" undistorted output is the same, 
%but inequality is higher).
OldThetaSpread=Para.theta_1-Para.theta_2;
NewThetaSpread=OldThetaSpread*1.3;
[c1FB c2FB l1FB l2FB yFB g_yFB_h Agent1WageShareFB_h]=getFB(Para,2);
Output(1)=c1FB*Para.n1+c2FB*Para.n2+Para.g(1);
[c1FB c2FB l1FB l2FB yFB g_yFB_l Agent1WageShareFB_l]=getFB(Para,1);
Output(2)=c1FB*Para.n1+c2FB*Para.n2+Para.g(2);
AverageOutput=sum(Output)/2;
NewTheta=fsolve(@(ThetaGuess) ResNewTheta( ThetaGuess,NewThetaSpread, AverageOutput,Para),[Para.theta_1 Para.theta_2]);
 %  --- SOLVE THE BELLMAN EQUATION --------------------------------------
 Main(Para)
 
 
   %  --- SIMULATE THE MODEL ----------------------------------------------
 
[sHist,gHist,u2btildHist,RHist,TauHist,YHist,TransHist,btildHist,...
c1Hist,c2Hist,l1Hist,l2Hist,IntHist,IncomeFromAssets_Agent1Hist,...
AfterTaxWageIncome_Agent1Hist,AfterTaxWageIncome_Agent2Hist]...
=RunSimulations(CoeffFileNane,Para.btild0,Para.NumSim,Para);
save( [Para.datapath 'SimData' casename '.mat'],'sHist','gHist','u2btildHist',...
'RHist','TauHist','YHist','TransHist','btildHist','btild0grid',...
'c1Hist','c2Hist','l1Hist','l2Hist','Para','IntHist',...
'AfterTaxWageIncome_Agent1Hist','AfterTaxWageIncome_Agent2Hist',...
'IncomeFromAssets_Agent1Hist');
 

 

 