% CAPTION : fig:flagPoints - This plots the sucess of the optimizer to
% solve the FOC at the points selected in the state space for the final set of coeffecients. The red points
% denote failure.
function  [sHist,u2btildHist,RHist,TauHist,YHist,TransHist,btildHist]=RunSimulations(startIter,endIter,btild0,NumSim,Para)
close all;
olddatapath=Para.datapath;
oldtexpath=Para.texpath;
oldplotpath=Para.plotpath;
plotpath=oldplotpath;
datapath=olddatapath;
for iter=startIter:endIter
    load([datapath 'c' num2str(iter) '.mat'])
    numsolved(iter)=length(IndxSolved);
%[Tau0,Rprime0,u2btildprime0]=SolveTime0(c,V,1,Para)
end
disp('Govt Exp')
g=Para.g

n1=Para.n1;
n2=Para.n2;
alpha_1=Para.alpha_1;
alpha_2=Para.alpha_2;
%disp('Govt Exp')
%g=Para.g
theta_1=Para.theta_1;
theta_2=Para.theta_2;
psi=Para.psi;
beta=Para.beta;

xSolved=x_state(IndxSolved,:);
xUnSolved=x_state(IndxUnSolved,:);

    

u2btildLL=Para.u2btildMin;
u2btildUL=Para.u2btildMax;
ucbtild_bounds = [u2btildLL,u2btildUL];
Rbounds=[min(Para.RGrid),max(Para.RGrid)];


% SOLVE THE T-0 PROBLEM given btild(-1)
btild_1=btild0;
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
[x,~,exitflagv0,output,lambda]  =fmincon(@(x) getValue0(x, btild_1,1,Para,c,V),[ 1 mean(Para.RGrid)^(-1)],[],[],[],[],lb,ub,[],opts);
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
btildHist(1)=u2btildprime0/uc0;
TauHist(1)=1-(ul0/(theta_2*uc0));
TransHist(1)=c20-l20*ul0/uc0;
RHist(1)=Rprime0;  
YHist(1)=n1*c10+n2*c20+g(1);
sHist(1)=1;
n=1;
tic
for i=1:NumSim
    if mod(i,100)==0
        disp('Running Simulation, t=')
        disp(i)
        toc
        tic
    end
    % STATE (t) - x,R,s_
    u2btild=u2btildHist(i);
    R=RHist(i);
   s_=sHist(i);
  % SOLVE THE BELLMAN EQUATION  
 [PolicyRulesInit]=GetInitialApproxPolicy([u2btild R s_] ,x_state,PolicyRulesStore);  
 [PolicyRules, V_new,exitflag,~]=CheckGradNAG(u2btild,R,s_,c,V,PolicyRulesInit,Para,0);
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
 y(1)=c1(1)*n1+c2(1)*n2+g(1);
 y(2)=c1(2)*n1+c2(2)*n2+g(2);
  
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


end