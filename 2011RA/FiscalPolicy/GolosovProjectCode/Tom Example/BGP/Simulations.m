% This program computes the 4 simulation. We need the following -
% 1. For each Case- Para,Coeff, functional space, domain, PolicyRulesStore
clear all
close all
SetParaStruc
btild_1=3;
NumSim=10;
ExName={'GammaHigh_n1High','GammaLow_n1Low','GammaHigh_n1Low','GammaLow_n1High'};
ExPath={'Data/GammaHigh_n1High_2/c298.mat','Data/GammaLow_nHigh/c298.mat'};
Sty={'k',':k'};
plotpath='Graphs/Simulations/';
  pwdd=pwd;

if strcmp(computer,'PCWIN')
    sl='\';
    coresize=2;
else
    sl='/';
    coresize=8;
end

compeconpath=[pwdd(1:end-length('\Tom Example\completelyStochasticCase')) sl 'compecon2011' sl];
addpath(genpath(compeconpath))

mkdir(plotpath);

% DRAW sHIST
sHist(1)=1;
for i=1:NumSim
   s_=sHist(i);
 if rand < Para.P(sHist(i),1)
    sHist(i+1)=1;
else
    sHist(i+1)=2;
end
end
 
figBtildSimul=figure('Name','btild');
figTauSimul=figure('Name','Tau');

for ExInd=1:2
% LOAD THE CASE
disp('Loading stuff for case number ')

disp(ExName{ExInd})
load(ExPath{ExInd})

theta_1=Para.theta_1;
theta_2=Para.theta_2;
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

% Initialize the simulation by solving the time0 problem
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
u2btildHist=zeros(NumSim,1);
btildHist=zeros(NumSim,1);
btildHist(1)=btildprime0;
RHist=zeros(NumSim,1);
u2btildHist(1)=u2btildprime0;
     RHist(1)=Rprime0;
     TauHist(1)=1-l10^(gamma)*c10^sigma/theta_1;

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

  RHist(i+1)=Rprime(sHist(i+1));
     u2btildHist(i+1)=u2btildprime(sHist(i+1)) ;
     btildHist(i+1)=btildprime(sHist(i+1)) ;
      TauHist(i+1)=Tau(sHist(i+1)); 

end
SimulStore(ExInd).btildHist=btildHist;
SimulStore(ExInd).TauHist=TauHist;
bbT{1}=0:1;
Inx=find(sHist(1:NumSim)<2);
    for i=2:length(Inx)-1
        bbT{i}=Inx(i)-1:Inx(i);
    end


figure(figBtildSimul)


plot((1:NumSim),SimulStore(ExInd).btildHist(1:NumSim),Sty{ExInd},'LineWidth',2)
axis tight
hold on
 
figure(figTauSimul)
subplot(2,1,ExInd)
  plot((1:NumSim),SimulStore(ExInd).TauHist(1:NumSim),Sty{ExInd},'LineWidth',2)
  axis tight
  ShadePlotForEmpahsis( bbT,'r',.05);
  ylabel('$\tau$','Interpreter','Latex')
hold on
 


end
figure(figBtildSimul)
ShadePlotForEmpahsis( bbT,'r',.05);  
legend('Gamma High','Gamma Low')
ylabel('$\tilde{b}_2$','Interpreter','Latex')
%print(gcf,'-dpng',[plotpath 'Simulation_btild.png'])
figure(figTauSimul)
%legend('Gamma High','Gamma Low')
%print(gcf,'-dpng',[plotpath 'Simulation_Tau.png'])

