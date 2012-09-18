clc
clear all
close all
SetParaStruc
theta_1=3.3;
theta_2=1;
n1=1;
n2=1;
tau=.2;
g_Y=.12;
AvfFETarget=.5;
x=fsolve(@(x) GetCalibrationFrischElasticity (x,AvfFETarget,theta_1,theta_2,tau,g_Y,n1,n2), [1 1 ]);
gamma=x(1)
Y=x(2)
g=g_Y*Y;
psi=1/(1+gamma);

beta=.9;
Para.beta=.9;
Para.alpha_1=alpha_1;
Para.alpha_2=alpha_2;
Para.psi=psi;
Para.g=ones(1,2).*g;
Para.theta_1=theta_1;
Para.theta_2=theta_2;
Para.btild_1=0;

%% Build Grid for the state variables
% This setups up the functional space and the grid.
u2btildMin=-10;
u2btildMax=-u2btildMin;
u2btildGrid=linspace(u2btildMin,u2btildMax,Para.u2btildGridSize);

Para.u2btildGrid=u2btildGrid;
Para.u2btildLL=u2btildMin;
Para.u2btildUL=u2btildMax;
Rbar0=5;
s_=1;
for u2btild_ind=1:Para.u2btildGridSize
  findR=@(R) getValueC1(u2btildGrid(u2btild_ind),R,s_,Para);
  [Rbar(u2btild_ind),~,exitval(u2btild_ind)]=fzero(findR,Rbar0);
  Rbar0=Rbar(u2btild_ind);
end


% R=u_2/u_1 = c1/c2
RMin=max(Rbar)*1.05;
RMax=max(Rbar)*2;

RGrid=linspace(RMin,RMax,Para.RGridSize);

Para.RGrid=RGrid;
GridSize=Para.u2btildGridSize*Para.RGridSize*Para.sSize;
Para.GridSize=GridSize;
Para.u2btildMin=u2btildMin;
Para.u2btildMax=u2btildMax;
RMax=max(RGrid);
Para.RMax=RMax;
RMin=min(RGrid);
Para.RMin=RMin;

alpha1GridSize=20;
alpha1grid=linspace(0.1, 0.85,alpha1GridSize);
tau1Target=.2;
c10guess=1;
c20guess=1/(RMax);
RGrid=[];
plotpath='Graphs/';

for i=1:alpha1GridSize
  alpha_1=alpha1grid(i);
alpha_2=1-alpha_1;
alpha_1=alpha_1*Para.n1;
alpha_2=alpha_2*Para.n2;
Para.alpha_1=alpha_1;
Para.alpha_2=alpha_2;
    
   [res]=GetSteadyStateMoments(Para,c10guess,c20guess);
     tau1(i)=res.tau1;
     AvgFrishElasticity(i)=res.AvgFrishElasticity;
     AvgHrs(i)=res.AvgHrs;
     exitflagv0(i)=res.exitflagv0;
     c10guess=res.c10;
     c20guess=res.c20;
     l1(i)=res.l1;
     l2(i)=res.l2;
  end
figure()
subplot(1,2,1)
plot(alpha1grid,l1,'k','LineWidth',2);
xlabel('alpha_1')
ylabel('l1')


subplot(1,2,2)
plot(alpha1grid,l2,'k','LineWidth',2);
xlabel('alpha_1')
ylabel('l2')

figure()
subplot(2,1,1)
plot(alpha1grid,tau1,'k','LineWidth',2);
hold on
plot(alpha1grid,tau1Target*ones(alpha1GridSize,1),':k','LineWidth',2);
xlabel('alpha_1')
ylabel('Labor Tax')

subplot(2,1,2)
plot(alpha1grid,AvgFrishElasticity,'k','LineWidth',2);
hold on
plot(alpha1grid,AvfFETarget*ones(alpha1GridSize,1),':k','LineWidth',2);
xlabel('alpha_1')
ylabel('Avg Frisch Elasticity')

print(gcf,'-depsc2 ',[plotpath 'CalibrationFE.eps'])
print(gcf,'-dpng ',[plotpath 'CalibrationFE.png'])



