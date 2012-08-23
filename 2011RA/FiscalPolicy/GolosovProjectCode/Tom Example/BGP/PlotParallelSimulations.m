% This script plots the long simulation using data stored in
% Data/SimDataParallel.mat file. 

close all
clear all

load( 'Data/SimDataParallel.mat')
mkdir ('Graphs/LongSimulations/FallingTaxes/');
mkdir ('Graphs/LongSimulations/RisingTaxes/');

% ---- CHANGE THIS AS PER THE CASE ---
plotpath='Graphs/LongSimulations/FallingTaxes';
% ----


K=size(u2btildHist,2);
IndxBenchMark=find(btild0grid==0); % find the index for the case u2btild=0, this is the benchmark case
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
 for ctrb=1:6
 subplot(3,2,ctrb)
 plot(TauHist(:,ctrb))
 xlabel('t')
ylabel('$\tau$','Interpreter','Latex')
title(['$\tilde{b}_0=$' num2str(btild0grid(ctrb)) ],'Interpreter','Latex')
 end
print(gcf,'-dpng',[plotpath 'SimulationsTau_all.png'])


 figure()
 for ctrb=1:6
 subplot(3,2,ctrb)
 plot(btildHist(:,ctrb))
 xlabel('t')
ylabel('$\tilde{b}_2$','Interpreter','Latex')
title(['$\tilde{b}_0=$' num2str(btild0grid(ctrb)) ],'Interpreter','Latex')
 end
 print(gcf,'-dpng',[plotpath 'SimulationsBtild_all.png'])

 u2bdiffFineGrid=linspace(min(u2btildHist(:,IndxBenchMark))*1.1,max(u2btildHist(:,IndxBenchMark))*1.1,35);
 RList=quantile(RHist(:,IndxBenchMark),[.5 .25 .5 .75]);
 s_=1;

 
LastIter=250;
load(['Data/c' num2str(LastIter) '.mat'])
SetParaStruc
GetTauPolicyPlots(u2bdiffFineGrid,mean(RHist(:,IndxBenchMark)),s_,LastIter,Para)


figure()
plot(btildHist(find(sHist(:,IndxBenchMark)==1),IndxBenchMark) ,TauHist(find(sHist(:,IndxBenchMark)==1),IndxBenchMark))
hold on
plot(btildHist(find(sHist(:,IndxBenchMark)==2),IndxBenchMark) ,TauHist(find(sHist(:,IndxBenchMark)==2),IndxBenchMark),'r')



 


%  cons(u1btildctr) = PolicyRules(1);
% end

figure(figFOCRes)
 subplot(2,2,Rctr)
 plot(u2bdiffFineGrid, FOCRes,'k')
 if Rctr==1
     legend('g_l','g_h')
 end
 xlabel('$x$','Interpreter','Latex')
 ylabel('FOCRes','Interpreter','Latex')
 title(['$R=$' num2str(RList(Rctr))])



