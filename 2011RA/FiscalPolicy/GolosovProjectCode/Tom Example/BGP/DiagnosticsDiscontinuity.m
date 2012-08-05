close all
clear all
SetParaStruc

clear numsolved
% CAPTION : fig:flagPoints - This plots the sucess of the optimizer to
% solve the FOC at the points selected in the state space for the final set of coeffecients. The red points
% denote failure.   
startIter=2;
endIter=15;
for iter=startIter:endIter
    
    
    load(['Data/c' num2str(iter) '.mat'])
    numsolved(iter)=length(IndxSolved);
%[Tau0,Rprime0,u2btildprime0]=SolveTime0(c,V,1,Para)
end
plotpath=[];
    xSolved=x_state(IndxSolved,:);
xUnSolved=x_state(IndxUnSolved,:);

    figure()

    scatter(squeeze(xSolved(:,1)),squeeze(xSolved(:,2)),'b','filled')
    hold on
    scatter(squeeze(xUnSolved(:,1)),squeeze(xUnSolved(:,2)),'r','filled')
    hold on
    xlabel('$x$','Interpreter','Latex')
    ylabel('$R$','Interpreter','Latex')
    print(gcf,'-dpng',[plotpath 'flagPoints' num2str(iter) '.png'])

% CAPTION : fig:NumFOCSolved - This plots shows the number of points that
% the FOC had a solution across iterations
 figu2BtildePrime =figure('Name','x');
 figBtildePrime =figure('Name','btild');
 figRprime=figure('Name','R');

figure()
plot(numsolved)
xlabel('Iter')
ylabel('Number of points FOC were saisfied')


figure()
plot((cdiff(max(round(iter*.1),1):iter)'))
xlabel('Iteration');
ylabel('Max of Coefficient Difference');
print(gcf,'-dpng',[plotpath 'CoeffConvergence.png'])



u2btildLL=Para.u2btildMin;
u2btildUL=Para.u2btildMax;
ucbtild_bounds = [u2btildLL,u2btildUL];
Rbounds=[min(Para.RGrid),max(Para.RGrid)];

 
 u2bdiffFineGrid=linspace(min(Para.u2bdiffGrid),max(Para.u2bdiffGrid),35);
 %RList=linspace(min(Para.RGrid),max(Para.RGrid),4);
 RList=linspace(2.5,2.51,4);
 s_=1;
for Rctr=1:4 
 for u2btildctr=1:length(u2bdiffFineGrid)
    R=RList(Rctr);
     u2btild=u2bdiffFineGrid(u2btildctr);
   [PolicyRulesInit]=GetInitialApproxPolicy([u2btild R s_] ,x_state,PolicyRulesStore);
    [PolicyRules, V_new,exitflag]=CheckGradNAG(u2btild,R,s_,c,V,PolicyRulesInit,Para)
    if exitflag==1
        IndxPrint(u2btildctr)=1;
    else
        IndxPrint(u2btildctr)=0;
        [u2btild R s_]
    end
    
% 
  u2BtildePrime(u2btildctr,:)=PolicyRules(end-1:end);
  BtildePrime(u2btildctr,:)=PolicyRules(end-5:end-4);
  Rprime(u2btildctr,:)=PolicyRules(end-3:end-2);
 end

%  cons(u1btildctr) = PolicyRules(1);
% end

figure(figBtildePrime)
 subplot(2,2,Rctr)
 plot(u2bdiffFineGrid(logical(IndxPrint)), BtildePrime(logical(IndxPrint),1),'k')
 hold on
 plot(u2bdiffFineGrid(logical(IndxPrint)), BtildePrime(logical(IndxPrint),2),':k')
 hold on
 if Rctr==1
     legend('g_l','g_h')
 end
 %plot(u2bdiffFineGrid, 0*u2bdiffFineGrid,':k');
 hold on
 %plot(u2bdiffFineGrid,repmat([u2btildLL u2btildUL],length(u2bdiffFineGrid),1)-[u2bdiffFineGrid' u2bdiffFineGrid'] ,':r')
% 
 xlabel('$x$','Interpreter','Latex')
 ylabel('$\tilde{b}_2$','Interpreter','Latex')
 title(['$R=$' num2str(RList(Rctr))])

 figure(figu2BtildePrime)
 subplot(2,2,Rctr)
 plot(u2bdiffFineGrid(logical(IndxPrint)), u2BtildePrime(logical(IndxPrint),1)- u2bdiffFineGrid(logical(IndxPrint))','k')
 hold on
 plot(u2bdiffFineGrid(logical(IndxPrint)), u2BtildePrime(logical(IndxPrint),2)-u2bdiffFineGrid(logical(IndxPrint))',':k')
 hold on
 if Rctr==1
     legend('g_l','g_h')
 end
 %plot(u2bdiffFineGrid, 0*u2bdiffFineGrid,':k');
 hold on
 %plot(u2bdiffFineGrid,repmat([u2btildLL u2btildUL],length(u2bdiffFineGrid),1)-[u2bdiffFineGrid' u2bdiffFineGrid'] ,':r')
% 
 xlabel('$x$','Interpreter','Latex')
 ylabel('$x*-x$','Interpreter','Latex')
 title(['$R=$' num2str(RList(Rctr))])

 
% 
 figure(figRprime)
 subplot(2,2,Rctr)
 plot(u2bdiffFineGrid(logical(IndxPrint)), Rprime(logical(IndxPrint),1),'k');
 hold on
 plot(u2bdiffFineGrid(logical(IndxPrint)), Rprime(logical(IndxPrint),2),':k');
 xlabel('$x$','Interpreter','Latex')
 ylabel('$R^{*}$','Interpreter','Latex')
  title(['$R=$' num2str(RList(Rctr))])
end
% % Caption : fig:ValueFunction - This plot depicts the value function 
figure()
NumIter=max(round((iter-1)/10),1);
MaxIter=iter;
ListIterations=(startIter:NumIter:endIter);
 for l=1:length(ListIterations)
     load(['Data/c' num2str(ListIterations(l)) '.mat'])

for Rctr=1:4
    subplot(2,2,Rctr)
    fplot(@(u2btild) funeval(c(s_,:)',V(s_),[u2btild RList(Rctr)]),[ucbtild_bounds(1) ucbtild_bounds(2)],'-k');
xlabel('$x$','Interpreter','Latex')
title(['$R=$' num2str(RList(Rctr))],'Interpreter','Latex')
hold on
end
end
% 
% 
% fplot(@(u1btild) funeval(c,V,u1btild),ucbtild_bounds,'-k');
% xlabel('$x$','Interpreter','Latex')
% ylabel('$V(x)$','Interpreter','Latex')
%  

% % Caption : fig:ValueFunction - This plot depicts the value function 
figure()
NumIter=max(round((iter-1)/10),1);
MaxIter=iter;
xlist=linspace(min(Para.u2bdiffGrid),max(Para.u2bdiffGrid),4)
ListIterations=(startIter:NumIter:endIter);
 for l=1:length(ListIterations)
     load(['Data/c' num2str(ListIterations(l)) '.mat'])

for xctr=1:4
    subplot(2,2,xctr)
    fplot(@(R) funeval(c(s_,:)',V(s_),[xlist(xctr) R]),[Rbounds(1) Rbounds(2)],'-k');
xlabel('$R$','Interpreter','Latex')
title(['$x=$' num2str(xlist(xctr))],'Interpreter','Latex')
hold on
end
end
% 
% 
% fplot(@(u1btild) funeval(c,V,u1btild),ucbtild_bounds,'-k');
% xlabel('$x$','Interpreter','Latex')
% ylabel('$V(x)$','Interpreter','Latex')
%  
clc
close all
x_state(1:length(Para.RGrid):Para.GridSize/2,1);
R=min(Para.RGrid);
u2bdiffFineGrid=linspace(1.8,1.9,5)
PolPlot=[];
for u2btildctr=1:length(u2bdiffFineGrid)
u2btild= u2bdiffFineGrid(u2btildctr);
s_=1;
  [PolicyRulesInit]=GetInitialApproxPolicy([u2btild R s_] ,x_state,PolicyRulesStore);
    [PolicyRules, V_new,exitflag]=CheckGradNAG(u2btild,R,s_,c,V,PolicyRulesInit,Para);
    PolPlot(u2btildctr,:)=PolicyRules(:,1);
     DiagonsticCheckConstraints(u2btild,R,s_,c,V,PolicyRulesInit*.96,Para)
end
figure()
plot(u2bdiffFineGrid,PolPlot)


figure()
plot(x_state(1:length(Para.RGrid):Para.GridSize/2,1),(PolicyRulesStore(1:length(Para.RGrid):Para.GridSize/2,end-1)-x_state(1:length(Para.RGrid):Para.GridSize/2,1)))