% This script plots the long simulation using data stored in
% Data/SimDataParallel.mat file. 

close all
clear all
clc

load( 'Data/SimDataParallelCommonShocks.mat')
mkdir ('Graphs/LongSimulations/RisingTaxes/CommonShocks');
mkdir ('Tex/LongSimulations/RisingTaxes/CommonShocks');

% ---- CHANGE THIS AS PER THE CASE ---
plotpath='Graphs/LongSimulations/RisingTaxes/CommonShocks/';
texpath='Tex/LongSimulations/RisingTaxes/CommonShocks/';
% ----


K=size(u2btildHist,2);
T=100;

% -- labor taxes ----------------------------------------------------------
X.data=TauHist;
X.sHist=sHist;
X.ylabel='tau';
X.name ='LaborTaxes';
PlotSimulationCommonshock( X,T,btild0grid,K,gHist,plotpath,texpath)
   

% -- btild ----------------------------------------------------------
X.data=btildHist;
X.sHist=sHist;
X.ylabel='b2~';
X.name ='RelativeAssetsAgent2';
PlotSimulationCommonshock( X,T,btild0grid,K,gHist,plotpath,texpath)


% -- Trans ----------------------------------------------------------
X.data=TransHist;
X.sHist=sHist;
X.ylabel='T';
X.name ='Transfers';
PlotSimulationCommonshock( X,T,btild0grid,K,gHist,plotpath,texpath)


% -- AfterTaxIncomeAgent1 ----------------------------------------------------------
X.data=AfterTaxWageIncome_Agent1Hist;
X.sHist=sHist;
X.ylabel='After-tax wage income';
X.name ='AfterTaxWageIncomeAgent1';
PlotSimulationCommonshock( X,T,btild0grid,K,gHist,plotpath,texpath)

% -- AfterTaxIncomeAgent2 ----------------------------------------------------------
X.data=AfterTaxWageIncome_Agent2Hist;
X.sHist=sHist;
X.ylabel='After-tax wage income';
X.name ='AfterTaxWageIncomeAgent2';
PlotSimulationCommonshock( X,T,btild0grid,K,gHist,plotpath,texpath)

% -- IncomeFromAssetsAgent2 ----------------------------------------------------------
X.data=IncomeFromAssets_Agent1Hist;
X.sHist=sHist;
X.ylabel='Asset Income';
X.name ='IncomeFromAssetsAgent1';
PlotSimulationCommonshock( X,T,btild0grid,K,gHist(1:end-1,:),plotpath,texpath)



