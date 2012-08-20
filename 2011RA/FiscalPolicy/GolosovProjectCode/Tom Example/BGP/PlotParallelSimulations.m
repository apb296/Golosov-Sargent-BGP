% This script plots the long simulation using data stored in
% Data/SimDataParallel.mat file. 

close all
clear all

load( 'Data/SimDataParallel.mat')
mkdir ('Graphs/LongSimulations/FallingTaxes/');
mkdir ('Graphs/LongSimulations/RisingTaxes/');

% ---- CHANGE THIS AS PER THE CASE ---
plotpath='Graphs/LongSimulations/RisingTaxes';
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
 
 