close all
clear all
SetParaStruc

clear numsolved
% CAPTION : fig:flagPoints - This plots the sucess of the optimizer to
% solve the FOC at the points selected in the state space for the final set of coeffecients. The red points
% denote failure.   
startIter=2;
endIter=25;
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

figure()
plot(numsolved)
xlabel('Iter')
ylabel('Number of points FOC were saisfied')


figure()
plot((cdiff(max(round(iter*.1),1):iter)'))
xlabel('Iteration');
ylabel('Max of Coefficient Difference');
print(gcf,'-dpng',[plotpath 'CoeffConvergence.png'])

load(['Data/c'  num2str(iter-1)  '.mat'])

c_old=c;

load(['Data/c'  num2str(iter)  '.mat'])
c=c_old;


u2btildLL=Para.u2btildMin;
u2btildUL=Para.u2btildMax;
ucbtild_bounds = [u2btildLL,u2btildUL];
Rbounds=[min(Para.RGrid),max(Para.RGrid)];

%Caption : fig:FunctionalConvergence - This figure plots the value function
% with respect to $\tilde{b}_2$ across iterations. The four panels refer to
% vaules of R. The red line is the first iteration

NumIter=max(round((iter-startIter)/10),1);
MaxIter=iter;
ListIterations=(startIter:NumIter:endIter);
figure()
% Fix s_
s_=1;
RList=linspace(Rbounds(1),Rbounds(2),4);
for l=1:length(ListIterations)
     load(['Data/c' num2str(ListIterations(l)) '.mat'])

for Rctr=1:4
    subplot(2,2,Rctr)
if l==1
fplot(@(u2btild) funeval(c(s_,:)',V(s_),[u2btild RList(Rctr)]),[ucbtild_bounds(1) ucbtild_bounds(2)],'-r');
else
    fplot(@(u2btild) funeval(c(s_,:)',V(s_),[u2btild RList(Rctr)]),[ucbtild_bounds(1) ucbtild_bounds(2)],'-k');
end
xlabel('$x$','Interpreter','Latex')
title(['$R=$' num2str(RList(Rctr))],'Interpreter','Latex')
hold on
end
end

print(gcf,'-dpng',[plotpath 'FunctionalConvergence.png'])

% 
% % ChebError




 numtest=50;
 for n=1:numtest
%     
%     
     u2btild=ucbtild_bounds(1)+(ucbtild_bounds(2)-ucbtild_bounds(1))*rand;
     R=Rbounds(1)+(Rbounds(2)-Rbounds(1))*rand;
     
xTarget(n,:)=[u2btild R s_];
 [PolicyRulesInit]=GetInitialApproxPolicy(xTarget(n,:),x_state,PolicyRulesStore);
 [PolicyRules, V_new,exitflag]=CheckGradNAG(u2btild,R,s_,c,V,PolicyRulesInit,Para) ;
 VDirect=funeval(c(s_,:)',V(s_),xTarget(n,1:2));
 Check(n)=(VDirect-V_new)/V_new;
 
 %Do optimization
 %Vopt = CheckOpt(u2btild,R,s_,c,V,PolicyRulesInit,Para);
 %Check2(n) = (Vopt-V_new)/V_new;
 %[Vopt1,Vopt2] = CheckOpt(u2btild,R,s_,c,V,PolicyRules(1:3),Para);
 %Check3(n) = (Vopt1-V_new)/V_new;
 %Check4(n) = (Vopt2-V_new)/V_new;
if ~(exitflag==1)
 colFOC(n,:)=[1 0 0];
 else
     colFOC(n,:)=[0 0 1];
 end
 

 end
 figure()
 plot(Check)
 xlabel('Number of Test Points')
 ylabel('Percentage Error')
 print(gcf,'-dpng',[plotpath 'ChebError.png'])
 
%  figure()
%  plot(Check2);
%  figure()
%  plot(Check3);
%  figure()
%  plot(Check4);
 
% 
% % Caption : fig:ValueFunction - This plot depicts the value function 
NumIter=round((iter-1)/10);
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
  print(gcf,'-dpng',[plotpath 'ValueFunction.png'])
%  
% 

% 
% SOLVE THE T-0 PROBLEM given btild(-1)
btild_1=Para.btild_1;

  disp('Computed V...Now solving V0(btild_1) where btild_1 is')
disp(btild_1)
% c1 and c2 solve 
 options=optimset('Display','off');
[x,fval,exitflagv0,~,grad] = fminunc(@(x)  getValue0(x, btild_1,1,Para,c,V),[ .5 .5*mean(Para.RGrid)^(-1)],options);
if ~(exitflagv0==1)
    disp('Optimization failed for V0 once ..trying with fmincon')
opts = optimset('Algorithm', 'interior-point', 'Display','off', ...
    'GradObj','off','GradConstr','off',...
    'MaxIter',1000, ...
    'TolX', Para.ctol/10, 'TolFun', Para.ctol, 'TolCon', Para.ctol,'MaxTime',200);
lb=[0.001 0.001];
ub=[10 10];
%[x,fval,exitflagv0,output,lambda]  =fmincon(@(x) getValue0(x, btild_1,1,Para,c,V),[ x ],[],[],[],[],lb,ub,[],opts);
[x,fval,exitflagv0,output,lambda]  =fmincon(@(x) getValue0(x, btild_1,1,Para,c,V),[ 1 mean(Para.RGrid)^(-1)],[],[],[],[],lb,ub,[],opts);
end
c10 = x(1);
c20 = x(2);
R0=c10/c20;
TotalResources=(c10*n1+c20*n2+g(1));
FF=R0*theta_2/theta_1;
DenL2=n1*theta_1*FF+theta_2*n2;
l20=(TotalResources-n1*theta_1+n1*theta_1*FF)/(DenL2);
l10= 1-FF*(1-l20);
BracketTerm=l20/(1-l20)-(l10/(1-l10))*R0;
u2btildprime0=(((1-psi)/(psi))*BracketTerm+btild_1/(beta*psi)+R0-1)*psi;
btildprime0=u2btildprime0/(c20^-1*psi) ;
Rprime0=c20^(-1)/c10^(-1);


% RUN SIMULATION
NumSim=110;
u2btildHist=zeros(NumSim,1);
btildHist=zeros(NumSim,1);
RHist=zeros(NumSim,1);
btildHist(1)=btild_1;
TauHist=zeros(NumSim,1);
YHist=zeros(NumSim,1);
TransHist=zeros(NumSim,1);
GMul=zeros(NumSim,1);
GBCCheck=zeros(NumSim,1);

% INITIALIZE - t=0
u2btildHist(1)=u2btildprime0;
ul0=(1-psi)/(1-l20);
uc0=psi/c20;
TauHist(1)=1-(ul0/(theta_2*uc0));
TransHist(1)=c20-l20*ul0/uc0;
RHist(1)=Rprime0;  
YHist(1)=n1*c10+n2*c20+g(1);
sHist(1)=1;
n=1;

for i=1:NumSim
    % STATE (t) - x,R,s_
    u2btild=u2btildHist(i);
    R=RHist(i);
   s_=sHist(i);
  % SOLVE THE BELLMAN EQUATION  
 [PolicyRulesInit]=GetInitialApproxPolicy([u2btild R s_] ,x_state,PolicyRulesStore);  
 [PolicyRules, V_new,exitflag]=CheckGradNAG(u2btild,R,s_,c,V,PolicyRulesInit,Para);
  % GET THE POLICY RULES - 
 %PolicyRules=[c1_1 c1_2 c2_1 c2_2 l1(1) l1(2) l2(1) l2(2) btildprime c2_1^(-1)/c1_1^(-1) c2_2^(-1)/c1_2^(-1) u2btildprime(1) u2btildprime(2)]
 c1=PolicyRules(1:2);
 c2=PolicyRules(3:4);
 l1=PolicyRules(5:6);
 l2=PolicyRules(7:8);
 ul2=(1-psi)./(1-l2);
uc2=psi./c2;

% TAU - From the WAGE optimality of Agent 2
Tau=1-(ul2./(theta_2.*uc2));
% RPRIME - By definition u_c_2/u_c_1;
 Rprime=PolicyRules(end-3:end-2);
  % OUTPUT
 y=c1*n1+c2*n2+g;
  
 % TRANSFERS
 % These are transfers computed on the assumption that Agent 2 cannot
 % borrow and lend. The transfers are the difference between his
 % consumption and after tax earning (l . U_l/U_c)
 Trans=c2-l2.*ul2./uc2;
 
% G MULTIPLIER - Computed using (yh-yl)/(gh-gl)
 GMul(i)=(y(2)-y(1))/(g(2)-g(1));
  % x' - u_c_2* btildprime
 u2btildprime=PolicyRules(end-1:end);
  % btildprime - x'/u_c2
 btildprime=PolicyRules(9:10);
 ExitFlag(i)=exitflag;
  % DRAW THE s' ~ P(s,:)

if rand < Para.P(sHist(i),1)
    sHist(i+1)=1;
else
    
    sHist(i+1)=2;
end
% UPDATE THE SIMULATION HISTORY

   RHist(i+1)=Rprime(sHist(i+1));
     u2btildHist(i+1)=u2btildprime(sHist(i+1)) ;
     btildHist(i+1)=btildprime(sHist(i+1)) ;
      TauHist(i+1)=Tau(sHist(i+1)); 
      YHist(i+1)=y(sHist(i+1));
     TransHist(i+1)=Trans(sHist(i+1));
     
%  if exitflag==1
%      RHist(n)=Rprime(sHist(i+1));
%      u2btildHist(n)=u2btildprime(sHist(i+1)) ;
%  n=n+1;
%  end

end
figCombSimul=figure('Name','Simulations');
figCombSimul_b=figure('Name','Simulations');
  % SIMULATION 1  - btild_2
    X.Data=btildHist;
    X.sHist=sHist;
    X.name={'$\tilde{b}_2$'};
    PlotSimul(X);
    print(gcf,'-dpng',[plotpath 'Simulation_btild.png'])

    % SIMULATION 2  - btild_2 Trunc
    T=100;
    X.Data=btildHist(1:T);
    X.sHist=sHist;
    X.name={'$\tilde{b}_2$'};
    PlotSimul(X);
    print(gcf,'-dpng',[plotpath 'Simulation_btildTrunc.png'])
      % redraw the plot for the subplot (1/4)
    figure(figCombSimul)
    subplot(2,2,1)
    PlotSimul(X,1);
    title('$b_2-b_1$','Interpreter','Latex')

     % SIMULATION 2  - btild_2 Trunc (last 100 periods)
    T=100;
    X.Data=btildHist(end-T+1:end);
    X.sHist=sHist(end-T+1:end);
    X.name={'$\tilde{b}_2$'};
    PlotSimul(X);
    print(gcf,'-dpng',[plotpath 'Simulation_btildTrunc_b.png'])

    % redraw the plot for the subplot (1/4)
    figure(figCombSimul_b)
    subplot(2,2,1)
    PlotSimul(X,1);
    title('$b_2-b_1$','Interpreter','Latex')


    % SIMULATION 3 - Tau
    X.Data=TauHist;
    X.sHist=sHist;
    X.name={'$\tau$'};
    PlotSimul(X);
    print(gcf,'-dpng',[plotpath 'Simulation_Tau.png']);

    % SIMULATION 4 - Tau Trunc
    T=100;
    X.Data=TauHist(1:T);
    X.sHist=sHist;
    X.name={'$\tau$'};  
    PlotSimul(X);
    print(gcf,'-dpng',[plotpath 'Simulation_TauTrunc.png'])
      % redraw the plot for the subplot (3/4)
    figure(figCombSimul)
    subplot(2,2,3)
    PlotSimul(X,1);
    title('Labor Tax')
      % SIMULATION 4b - Tau Trunc
    T=100;
    X.Data=TauHist(end-T+1:end);
    X.sHist=sHist(end-T+1:end);
    X.name={'$\tau$'};  
    PlotSimul(X);
    print(gcf,'-dpng',[plotpath 'Simulation_TauTrunc_b.png'])

    % redraw the plot for the subplot (3/4)
    figure(figCombSimul_b)
    subplot(2,2,3)
    PlotSimul(X,1);
    title('Labor Tax')
    
   
    % SIMULATION 5 - x
    X.Data=u2btildHist;
    X.sHist=sHist;
    X.name={'$x$'};
    PlotSimul(X);
    print(gcf,'-dpng',[plotpath 'Simulation_u2btild.png'])

    % SIMULATION 6 - R
    X.Data=RHist;
    X.sHist=sHist;
    X.name={'$R$'};
    PlotSimul(X);
    print(gcf,'-dpng',[plotpath 'Simulation_R.png'])

    % SIMULATION 7 - Y
    T=100;
    X.Data=YHist(1:T);
    X.sHist=sHist;
    X.name={'$y$'};
    PlotSimul(X);
    print(gcf,'-dpng',[plotpath 'Simulation_YTrunc.png'])
      % redraw the plot for the subplot (4/4)
    figure(figCombSimul)
    subplot(2,2,4)
    PlotSimul(X,1);
    title('Output')
    
     % SIMULATION 7 - Y
    T=100;
    X.Data=YHist(end-T+1:end);
    X.sHist=sHist(end-T+1:end);
    X.name={'$y$'};
    PlotSimul(X);
    print(gcf,'-dpng',[plotpath 'Simulation_YTrunc_b.png'])
   
    % redraw the plot for the subplot (4/4)
    figure(figCombSimul_b)
    subplot(2,2,4)
    PlotSimul(X,1);
    title('Output')
    % SIMULATION 8 - T
    T=100;
    X.Data=TransHist(1:T);
    X.sHist=sHist;
    X.name={'$T$'};
    PlotSimul(X);
    print(gcf,'-dpng',[plotpath 'Simulation_TransTrunc.png'])
     % redraw the plot for the subplot (2/4)
    figure(figCombSimul)
    subplot(2,2,2)
    PlotSimul(X,1);
    title('Transfers')
    print(gcf,'-dpng',[plotpath 'Simulations.png'])
     % SIMULATION 8 - T
    T=100;
    X.Data=TransHist(end-T+1:end);
    X.sHist=sHist(end-T+1:end);
    X.name={'$T$'};
    PlotSimul(X);
    print(gcf,'-dpng',[plotpath 'Simulation_TransTrunc_b.png'])
    
    % redraw the plot for the subplot (2/4)
    figure(figCombSimul_b)
    subplot(2,2,2)
    PlotSimul(X,1);
    title('Transfers')
    print(gcf,'-deps',[plotpath 'Simulations_b1.ps'])
    print(gcf,'-depsc',[plotpath 'Simulations_b2.ps'])
    print(gcf,'-dpng',[plotpath 'Simulations_b.png'])


    % SIMULATION 9 - GMul Hist
    figure()
     hist(GMul,50)
    print(gcf,'-dpng',[plotpath 'Simulation_HistGMul.png'])


 % Policy Rules
 % Policy Rules
% Caption : fig:PolicyRules - This plot depicts the $\tilde{b}'_2$ as a function of $\tilde{b}_2$ 
 figu2BtildePrime =figure('Name','x');
 figBtildePrime =figure('Name','btild');
 figRprime=figure('Name','R');
 
 u2bdiffFineGrid=linspace(min(u2btildHist),max(u2btildHist),35);
 RList=quantile(RHist,[.5 .25 .5 .75]);
 s_=1;
for Rctr=1:4 
 for u2btildctr=1:length(u2bdiffFineGrid)
    R=RList(Rctr);
     u2btild=u2bdiffFineGrid(u2btildctr);
   [PolicyRulesInit]=GetInitialApproxPolicy([u2btild R s_] ,x_state,PolicyRulesStore);
    [PolicyRules, V_new,exitflag]=CheckGradNAG(u2btild,R,s_,c,V,PolicyRulesInit,Para);
    if exitflag==1
        IndxPrint(u2btildctr)=1;
    else
        IndxPrint(u2btildctr)=0;
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
 print(figu2BtildePrime,'-dpng',[plotpath 'u2BtildePrime.png'])
  print(figBtildePrime,'-dpng',[plotpath 'BtildePrime.png'])
%  
  print(figRprime,'-dpng',[plotpath 'RPrime.png'])
% 
% 

%% Policy Rules entire state space
% Caption : fig:PolicyRules - This plot depicts the $\tilde{b}'_2$ as a function of $\tilde{b}_2$ 
 figu2BtildePrime =figure('Name','x');
 figBtildePrime =figure('Name','btild');
 figRprime=figure('Name','R');
 
 u2bdiffFineGrid=linspace(min(Para.u2bdiffGrid),max(Para.u2bdiffGrid),35);
 RList=linspace(min(Para.RGrid),max(Para.RGrid),4);
 s_=1;
for Rctr=1:4 
 for u2btildctr=1:length(u2bdiffFineGrid)
    R=RList(Rctr);
     u2btild=u2bdiffFineGrid(u2btildctr);
   [PolicyRulesInit]=GetInitialApproxPolicy([u2btild R s_] ,x_state,PolicyRulesStore);
    [PolicyRules, V_new,exitflag]=CheckGradNAG(u2btild,R,s_,c,V,PolicyRulesInit,Para);
    if exitflag==1
        IndxPrint(u2btildctr)=1;
    else
        IndxPrint(u2btildctr)=0;
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
