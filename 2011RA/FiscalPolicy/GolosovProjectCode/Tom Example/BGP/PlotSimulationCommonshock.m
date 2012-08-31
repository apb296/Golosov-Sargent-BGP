function  PlotSimulationCommonshock( X,T,btild0grid,K,gHist,plotpath,texpath)
BurnSampleRatio=.5;                                                         % Percentage of simulations to disregard
figure()
subplot(2,1,1)
plot(X.data(:,1))
xlabel('t')
ylabel(X.ylabel,'Interpreter','Latex')
title('btild0=0')
subplot(2,1,2)
plot(X.data(:,K))
xlabel('t')
ylabel(X.ylabel,'Interpreter','Latex')
title(['btild0=' num2str(btild0grid(K))])
print(gcf,'-dpng',[plotpath 'LongSimulations' X.name '.png'])


 figure()
    subplot(2,1,1)
    XX.Data=X.data(end-T+1:end,1);
    XX.sHist=X.sHist(end-T+1:end,1);
    XX.name=X.ylabel;  
    PlotSimul(XX,1);
    title([X.name ' - Last 100 periods with btild0=0'])
    subplot(2,1,2)
     XX.Data=X.data(end-T+1:end,K);
    XX.sHist=X.sHist(end-T+1:end,1);
    XX.name=X.ylabel;  
    PlotSimul(XX,1);
    title([X.name ' - Last 100 periods with btild0=' num2str(btild0grid(K))])
print(gcf,'-dpng',[plotpath 'TruncSimulations' X.name 'Last100.png'])

% -- moments -------------------------------------------------------------
Moments(1,1) =mean(X.data(end*BurnSampleRatio:end,1));
Moments(1,2)=std(X.data(end*BurnSampleRatio:end,1));
Moments(1,3)=corr(X.data(end*BurnSampleRatio:end,1),X.data(end*BurnSampleRatio-1:end-1,1));
Moments(1,4)=corr(X.data(end*BurnSampleRatio:end,1),gHist(end*BurnSampleRatio:end,1));


Moments(2,1) =mean(X.data(end*BurnSampleRatio:end,K));
Moments(2,2)=std(X.data(end*BurnSampleRatio:end,K));
Moments(2,3)=corr(X.data(end*BurnSampleRatio:end,K),X.data(end*BurnSampleRatio-1:end-1,1));
Moments(2,4)=corr(X.data(end*BurnSampleRatio:end,K),gHist(end*BurnSampleRatio:end,1));


rowLabels = {'btild0=0',['btild0=' num2str(btild0grid(K))]};
columnLabels = {'Mean','Std','AutoCorr','Corr with g'};
matrix2latex(Moments, [texpath X.name 'Moments.tex'] , 'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', 'c', 'format', '%-6.4f', 'size', 'tiny');



end

