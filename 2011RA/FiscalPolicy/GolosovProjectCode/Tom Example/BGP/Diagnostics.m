close all
clear all
SetParaStruc

clear numsolved
minIter = 1;
maxIter = 119;

load(['Data/c' num2str(maxIter) '.mat'])

theta_1=Para.theta_1;
theta_2=Para.theta_2;
Para.Sigma = 
sigma=Para.sigma;
gamma=Para.gamma;
g=Para.g;
alpha_1=Para.alpha_1;
alpha_2=Para.alpha_2;
P=Para.P;
n1=Para.n1;
n2=Para.n2;
beta=Para.beta;
u2btildLL=Para.u2btildLL;
u2btildUL=Para.u2btildUL;

btild_1=Para.btild_1;
  disp('Computed V...Now solving V0(btild_1) where btild_1 is')
disp(btild_1)
% c1 and c2 solve 
 options=optimset('Display','off');
[x,fval,exitflagv0,~,grad] = fminunc(@(x)  getValue0(x, btild_1,1,Para,c,V),[ 1 mean(Para.RGrid)^(-1/sigma)],options);
if ~(exitflagv0==1)
    disp('Optimization failed for V0 once ..trying with fmincon')
opts = optimset('Algorithm', 'interior-point', 'Display','off', ...
    'GradObj','off','GradConstr','off',...
    'MaxIter',1000, ...
    'TolX', Para.ctol/10, 'TolFun', Para.ctol, 'TolCon', Para.ctol,'MaxTime',200);
lb=[0.001 0.001];
ub=[10 10];
%[x,fval,exitflagv0,output,lambda]  =fmincon(@(x) getValue0(x, btild_1,1,Para,c,V),[ x ],[],[],[],[],lb,ub,[],opts);
[x,fval,exitflagv0,output,lambda]  =fmincon(@(x) getValue0(x, btild_1,1,Para,c,V),[ 1 mean(Para.RGrid)^(-1/sigma)],[],[],[],[],lb,ub,[],opts);
end
c10 = x(1);
c20 = x(2);
TotalResoucrces=(c10*n1+c20*n2+g(1));
DenL1=theta_1*n1+(c20/c10)^(-sigma/gamma)*(theta_2/theta_1)^(1/gamma)*n2*theta_2;
l10=TotalResoucrces/DenL1;
l20=(c20/c10)^(-sigma/gamma)*(theta_2/theta_1)^(1/gamma)*l10;
btildprime0=c20^(sigma)*btild_1/beta-(c20-c10-l20^(1+gamma)*c20^(sigma)+l10^(1+gamma)*c10^(sigma));
u2btildprime0=c20^(-sigma)*btildprime0;
Rprime0=c20^(-sigma)/c10^(-sigma);
NumSim=500;
u2btildHist=zeros(NumSim,1);
btildHist=zeros(NumSim,1);
btildHist(1)=btildprime0;
RHist=zeros(NumSim,1);
TauHist=zeros(NumSim,1);
u2btildHist(1)=u2btildprime0;
     RHist(1)=Rprime0;
     TauHist(1)=1-l10^(gamma)*c10^sigma/theta_1;
     TauHist(1)=1-l20^(gamma)*c20^sigma/theta_2;
sHist(1)=1;
n=1;
for i=1:NumSim
    
    u2btild=u2btildHist(i);
    R=RHist(i);
   s_=sHist(i);
 [PolicyRulesInit]=GetInitialApproxPolicy([u2btild R s_] ,x_state,PolicyRulesStore);  
 [PolicyRules, V_new,exitflag]=CheckGradNAG(u2btild,R,s_,c,V,PolicyRulesInit,Para);
 c1=PolicyRules(1:2);
 c2=PolicyRules(3:4);
 l1=PolicyRules(5:6);
 l2=PolicyRules(7:8);
 Tau=1-l1.^gamma.*c1.^(sigma)./theta_1;
 Rprime=PolicyRules(end-3:end-2);
 u2btildprime=PolicyRules(end-1:end);
 btildprime=PolicyRules(9:10);
 ExitFlag(i)=exitflag;
if rand < Para.P(sHist(i),1)
    sHist(i+1)=1;
else
    sHist(i+1)=2;
end
  RHist(i+1)=Rprime(sHist(i+1));
     u2btildHist(i+1)=u2btildprime(sHist(i+1)) ;
     btildHist(i+1)=btildprime(sHist(i+1)) ;
      TauHist(i+1)=Tau(sHist(i+1)); 
     
%  if exitflag==1
%      RHist(n)=Rprime(sHist(i+1));
%      u2btildHist(n)=u2btildprime(sHist(i+1)) ;
%  n=n+1;
%  end

end
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


% SIMULATION 4b - Tau Trunc from 1:T
T=100;
X.Data=TauHist(2:T);
X.sHist=sHist;
X.name={'$\tau$'};
PlotSimul(X);
print(gcf,'-dpng',[plotpath 'Simulation_TauTrunc.png'])


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
 




% CAPTION : fig:flagPoints - This plots the sucess of the optimizer to
% solve the FOC at the points selected in the state space for the final set of coeffecients. The red points
% denote failure. 

for iter=minIter:maxIter
    
    
    
    load(['Data/c' num2str(iter) '.mat'])
    numsolved(iter)=length(IndxSolved);

end
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


u2btildLL=Para.u2btildMin_1;
u2btildUL=Para.u2btildMax_1;
ucbtild_bounds = [u2btildLL,u2btildUL];
Rbounds=[min(Para.RGrid),max(Para.RGrid)];

%Caption : fig:FunctionalConvergence - This figure plots the value function
% with respect to $\tilde{b}_2$ across iterations. The four panels refer to
% vaules of R. The red line is the first iteration

NumIter=round((iter-1)/20);
MaxIter=iter;
ListIterations=(minIter:NumIter:maxIter);
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




 numtest=5;
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
 Vopt = CheckOpt(u2btild,R,s_,c,V,PolicyRulesInit,Para);
 Check2(n) = (Vopt-V_new)/V_new;
 [Vopt1,Vopt2] = CheckOpt(u2btild,R,s_,c,V,PolicyRules(1:3),Para);
 Check3(n) = (Vopt1-V_new)/V_new;
 Check4(n) = (Vopt2-V_new)/V_new;
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
 
 figure()
 plot(Check2);
 figure()
 plot(Check3);
 figure()
 plot(Check4);
 
% 
% % Caption : fig:ValueFunction - This plot depicts the value function 
NumIter=round((iter-1)/2);
MaxIter=iter;
ListIterations=(minIter:NumIter:iter);
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
 % Policy Rules
% Caption : fig:PolicyRules - This plot depicts the $\tilde{b}'_2$ as a function of $\tilde{b}_2$ 
 figu2BtildePrime =figure('Name','x');
 figBtildePrime =figure('Name','btild');
 figRprime=figure('Name','R');
 
 u2bdiffFineGrid=linspace(ucbtild_bounds(1),ucbtild_bounds(2),35);
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
% 
% simulation
     