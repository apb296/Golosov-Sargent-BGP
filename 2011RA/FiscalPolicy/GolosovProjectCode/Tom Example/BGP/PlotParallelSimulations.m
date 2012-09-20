% This script plots the long simulation using data stored in
% Data/SimDataParallel.mat file. 

close all
clear all
clc

load( 'Data/SimDataParallel.mat')
mkdir ('Graphs/lowEpsilon/');
mkdir ('Tex/lowEpsilon/');
mkdir ('Graphs/lowEpsilon/');
mkdir ('Tex/lowEpsilon/');

% ---- CHANGE THIS AS PER THE CASE ---
plotpath='Graphs/lowEpsilon/';
texpath='Tex/lowEpsilon/';
% ----


K=size(u2btildHist,2);
IndxBenchMark=find(btild0grid==0); % find the index for the case u2btild=0, this is the benchmark case



BudgetCheck_Agent1=IncomeFromAssets_Agent1Hist+AfterTaxWageIncome_Agent1Hist(2:end,:)-c1Hist(2:end,:)+btildHist(2:end,:); 
BudgetCheck_Agent2=AfterTaxWageIncome_Agent2Hist(1:end-1,:)-c2Hist(1:end-1,:);
if( max(BudgetCheck_Agent1)+max(BudgetCheck_Agent2) < 1e-5)
    disp('No errors in budget constraints')
else
      disp('check errors in budget constraints')
end


% --------------- Moments -------------------------------------------------
BurnSampleRatio=.5;                                                         % Percentage of simulations to disregard

% structure of moments  
%--------------------------------------------------
% | Mean  |  Std Dev | Persistence | Corr with g
% -------------------------------------------------

% 1. g
Moments(1,1) =mean(gHist(end*BurnSampleRatio:end,IndxBenchMark));
Moments(1,2)=std(gHist(end*BurnSampleRatio:end,IndxBenchMark));
Moments(1,3)=corr(gHist(end*BurnSampleRatio:end,IndxBenchMark),gHist(end*BurnSampleRatio-1:end-1,IndxBenchMark));
Moments(1,4)=corr(gHist(end*BurnSampleRatio:end,IndxBenchMark),gHist(end*BurnSampleRatio:end,IndxBenchMark));

% 2. Tau
Moments(2,1) =mean(TauHist(end*BurnSampleRatio:end,IndxBenchMark));
Moments(2,2)=std(TauHist(end*BurnSampleRatio:end,IndxBenchMark));
Moments(2,3)=corr(TauHist(end*BurnSampleRatio:end,IndxBenchMark),TauHist(end*BurnSampleRatio-1:end-1,IndxBenchMark));
Moments(2,4)=corr(TauHist(end*BurnSampleRatio:end,IndxBenchMark),gHist(end*BurnSampleRatio:end,IndxBenchMark));

% 3. Trans
Moments(3,1) =mean(TransHist(end*BurnSampleRatio:end,IndxBenchMark));
Moments(3,2)=std(TransHist(end*BurnSampleRatio:end,IndxBenchMark));
Moments(3,3)=corr(TransHist(end*BurnSampleRatio:end,IndxBenchMark),TransHist(end*BurnSampleRatio-1:end-1,IndxBenchMark));
Moments(3,4)=corr(TransHist(end*BurnSampleRatio:end,IndxBenchMark),gHist(end*BurnSampleRatio:end,IndxBenchMark));

% 4. btild
Moments(4,1) =mean(btildHist(end*BurnSampleRatio:end,IndxBenchMark));
Moments(4,2)=std(btildHist(end*BurnSampleRatio:end,IndxBenchMark));
Moments(4,3)=corr(btildHist(end*BurnSampleRatio:end,IndxBenchMark),btildHist(end*BurnSampleRatio-1:end-1,IndxBenchMark));
Moments(4,4)=corr(btildHist(end*BurnSampleRatio:end,IndxBenchMark),gHist(end*BurnSampleRatio:end,IndxBenchMark));


% 5. AfterTaxIncome_Agent1
Moments(5,1) =mean(AfterTaxWageIncome_Agent1Hist(end*BurnSampleRatio:end,IndxBenchMark));
Moments(5,2)=std(AfterTaxWageIncome_Agent1Hist(end*BurnSampleRatio:end,IndxBenchMark));
Moments(5,3)=corr(AfterTaxWageIncome_Agent1Hist(end*BurnSampleRatio:end,IndxBenchMark),AfterTaxWageIncome_Agent1Hist(end*BurnSampleRatio-1:end-1,IndxBenchMark));
Moments(5,4)=corr(AfterTaxWageIncome_Agent1Hist(end*BurnSampleRatio:end,IndxBenchMark),gHist(end*BurnSampleRatio:end,IndxBenchMark));



% 6. AfterTaxIncome_Agent2
Moments(6,1) =mean(AfterTaxWageIncome_Agent2Hist(end*BurnSampleRatio:end,IndxBenchMark));
Moments(6,2)=std(AfterTaxWageIncome_Agent2Hist(end*BurnSampleRatio:end,IndxBenchMark));
Moments(6,3)=corr(AfterTaxWageIncome_Agent2Hist(end*BurnSampleRatio:end,IndxBenchMark),AfterTaxWageIncome_Agent2Hist(end*BurnSampleRatio-1:end-1,IndxBenchMark));
Moments(6,4)=corr(AfterTaxWageIncome_Agent2Hist(end*BurnSampleRatio:end,IndxBenchMark),gHist(end*BurnSampleRatio:end,IndxBenchMark));



% 7. Income from Assets Agent 1 
Moments(7,1) =mean(IncomeFromAssets_Agent1Hist(end*BurnSampleRatio:end,IndxBenchMark));
Moments(7,2)=std(IncomeFromAssets_Agent1Hist(end*BurnSampleRatio:end,IndxBenchMark));
Moments(7,3)=corr(IncomeFromAssets_Agent1Hist(end*BurnSampleRatio:end,IndxBenchMark),IncomeFromAssets_Agent1Hist(end*BurnSampleRatio-1:end-1,IndxBenchMark));
Moments(7,4)=corr(IncomeFromAssets_Agent1Hist(end*BurnSampleRatio:end,IndxBenchMark),gHist(end*BurnSampleRatio:end-1,IndxBenchMark));

% 8. Interest Rates
Moments(8,1) =mean(IntHist(end*BurnSampleRatio:end,IndxBenchMark));
Moments(8,2)=std(IntHist(end*BurnSampleRatio:end,IndxBenchMark));
Moments(8,3)=corr(IntHist(end*BurnSampleRatio:end,IndxBenchMark),IntHist(end*BurnSampleRatio-1:end-1,IndxBenchMark));
Moments(8,4)=corr(IntHist(end*BurnSampleRatio:end,IndxBenchMark),gHist(end*BurnSampleRatio:end-1,IndxBenchMark));


rowLabels = {'Shock : $g$','Labor Taxes: $\tau$', 'Transfers :$T$','Relative assets  :$\tilde{b}_2$', 'After-tax income (Agent 1)','After-tax income (Agent 2)','Income from assets (Agent 1)','Gross Int. rates'};
columnLabels = {'Mean','Std','AutoCorr','Corr with g'};
matrix2latex(Moments, [texpath 'SimulationMoments.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.4f', 'size', 'tiny');


figure()
subplot(3,1,1)
plot(AfterTaxWageIncome_Agent1Hist(2:end,IndxBenchMark))
xlabel('t')
title('After Tax income (Agent1)')
subplot(3,1,2)
plot(AfterTaxWageIncome_Agent2Hist(2:end,IndxBenchMark))
xlabel('t')
title('After Tax income (Agent2)')

subplot(3,1,3)
plot(IncomeFromAssets_Agent1Hist(1:end,IndxBenchMark))
xlabel('t')
title('Income from Assets (Agent1)')
print(gcf,'-dpng',[plotpath 'SimulationsIncome.png'])

figure()
hist(IntHist(1:end,IndxBenchMark))
title('Interest Rates')
print(gcf,'-dpng',[plotpath 'SimulationsInterest.png'])


figure()
subplot(2,1,1)
plot(TauHist(:,IndxBenchMark))
xlabel('t')
ylabel('$\tau$','Interpreter','Latex')
title('Labor taxes')
subplot(2,1,2)
plot(btildHist(:,IndxBenchMark))
xlabel('t')
ylabel('$\tilde{b}_2$','Interpreter','Latex')
title('Relative assets of Agent 2')
print(gcf,'-dpng',[plotpath 'SimulationsInitBtild_0.png'])

 figure()
 for ctrb=1:K
 subplot(3,2,ctrb)
 plot(TauHist(:,ctrb))
 xlabel('t')
ylabel('$\tau$','Interpreter','Latex')
title(['$\tilde{b}_0=$' num2str(btild0grid(ctrb)) ],'Interpreter','Latex')
 end
print(gcf,'-dpng',[plotpath 'SimulationsTau_all.png'])


 figure()
 for ctrb=1:K
 subplot(3,2,ctrb)
 plot(btildHist(:,ctrb))
 xlabel('t')
ylabel('$\tilde{b}_2$','Interpreter','Latex')
title(['$\tilde{b}_0=$' num2str(btild0grid(ctrb)) ],'Interpreter','Latex')
 end
 print(gcf,'-dpng',[plotpath 'SimulationsBtild_all.png'])

 u2bdiffFineGrid=linspace(min(u2btildHist(:,IndxBenchMark))*2.5,max(u2btildHist(:,IndxBenchMark))*1.5,35);
 RList=quantile(RHist(:,IndxBenchMark),[.5 .25 .5 .75]);


figure()
plot(btildHist(find(sHist(:,IndxBenchMark)==1),IndxBenchMark) ,TauHist(find(sHist(:,IndxBenchMark)==1),IndxBenchMark))
hold on
plot(btildHist(find(sHist(:,IndxBenchMark)==2),IndxBenchMark) ,TauHist(find(sHist(:,IndxBenchMark)==2),IndxBenchMark),'r')
T=100;

 figure()
subplot(2,1,1)
    X.Data=TauHist(end-T+1:end,IndxBenchMark);
    X.sHist=sHist(end-T+1:end,IndxBenchMark);
    X.name={'$\tau$'};  
    PlotSimul(X,1);
    title('Labor Taxes - Last 100 periods')
subplot(2,1,2)
 
    X.Data=TauHist(2:T+1,IndxBenchMark);
    X.sHist=sHist(2:T+1,IndxBenchMark);
    X.name={'$\tau$'};  
    PlotSimul(X,1);
        title('Labor Taxes - First 100 periods')
    print(gcf,'-dpng',[plotpath 'Simulation_TauTrunc.png'])

    
    figure()
subplot(2,1,1)
    X.Data=TransHist(end-T+1:end,IndxBenchMark);
    X.sHist=sHist(end-T+1:end,IndxBenchMark);
    X.name={'T'};  
    PlotSimul(X,1);
    title('Transfers - Last 100 periods')
subplot(2,1,2)
 
 
    X.Data=TransHist(2:T+1,IndxBenchMark);
    X.sHist=sHist(2:T+1,IndxBenchMark);
    X.name={'T'};  
    PlotSimul(X,1);
        title('Transfers - First100 periods')
    print(gcf,'-dpng',[plotpath 'Simulation_TransTrunc.png'])

    
    % Truncated Income plots 
    % A. First 100 periods
    figure()
    subplot(2,2,1)
    X.Data=AfterTaxWageIncome_Agent1Hist(2:T,IndxBenchMark);
    X.sHist=sHist(2:T,IndxBenchMark);
  X.name={'I_1'};  
    PlotSimul(X,1);
    title('After-tax income (Agent 1)')
    
     subplot(2,2,2)
    X.Data=AfterTaxWageIncome_Agent2Hist(2:T,IndxBenchMark);
    X.sHist=sHist(2:T,IndxBenchMark);
  X.name={'I_2'};  
    PlotSimul(X,1);
    title('After-tax income (Agent 2)')
    
    subplot(2,2,3)
    X.Data=IncomeFromAssets_Agent1Hist(1:T-1,IndxBenchMark);
    X.sHist=sHist(1:T-1,IndxBenchMark);
    X.name={'I_1_b'};  
    PlotSimul(X,1);
    title('Income from assets (Agent 1)')
   
     subplot(2,2,4)
    X.Data=IntHist(1:T-1,IndxBenchMark);
    X.sHist=sHist(1:T-1,IndxBenchMark);
      X.name={'R'};  
    PlotSimul(X,1);
    title('Int rates')
    
       print(gcf,'-dpng',[plotpath 'Simulation_IncomeTrunc.png'])
       
   % B. Last 100 periods
    figure()
    subplot(2,2,1)
    X.Data=AfterTaxWageIncome_Agent1Hist(end-T:end,IndxBenchMark);
    X.sHist=sHist(end-T:end,IndxBenchMark);
     X.name={'I_1'};  
    PlotSimul(X,1);
    title('After-tax income (Agent 1)')
    
     subplot(2,2,2)
    X.Data=AfterTaxWageIncome_Agent2Hist(end-T:end,IndxBenchMark);
    X.sHist=sHist(end-T:end,IndxBenchMark);
    X.name={'I_2'};  
    PlotSimul(X,1);
    title('After-tax income (Agent 2)')
     
    subplot(2,2,3) 
    X.Data=IncomeFromAssets_Agent1Hist(end-T:end,IndxBenchMark);
    X.sHist=sHist(end-T:end,IndxBenchMark);
   X.name={'I_1_b'};  
    PlotSimul(X,1);
    title('Income from assets (Agent 1)')
    
     subplot(2,2,4)
    X.Data=IntHist(end-T:end,IndxBenchMark);
    X.sHist=sHist(end-T:end,IndxBenchMark);
    X.name={'R'};  
    PlotSimul(X,1);
    title('Int rates')
    
       print(gcf,'-dpng',[plotpath 'Simulation_IncomeTrunc_b.png'])


  s_=1;
oldplotpath=plotpath;
LastIter=500;
load(['Data/c' num2str(LastIter) '.mat'])
SetParaStruc
GetTauPolicyPlots(u2bdiffFineGrid,mean(RHist(:,IndxBenchMark)),s_,LastIter,Para,oldplotpath)
   
    
   