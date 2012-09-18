clc
clear all
close all
SetParaStruc
x=fsolve(@(x) GetCalibration(x), [1 .5 1]);
theta_2=min(x(1),1);
gamma=x(2);
Y=x(3);
g=.12*Y;
psi=1/(1+gamma);
theta_1=max(x(1),1);
n1=1;
n2=1;
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
RMax=max(Rbar)*1.8;

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



alpha1GridSize=30;
alpha1grid=linspace(0.4, 0.9,alpha1GridSize);
tau1Target=.2;
GiniCoeffTarget=.38;
AvgHrsTarget=.23;
c10guess=.5;
c20guess=.5/(7);
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
     GiniCoeff(i)=res.GiniCoeff;
     AvgHrs(i)=res.AvgHrs;
     exitflagv0(i)=res.exitflagv0;
     VarHrs(i)=res.VarHrs;
  end


figure()
subplot(3,1,1)
plot(alpha1grid,tau1,'k','LineWidth',2);
hold on
plot(alpha1grid,tau1Target*ones(alpha1GridSize,1),':k','LineWidth',2);
xlabel('alpha_1')
ylabel('Labor Tax')
subplot(3,1,2)
plot(alpha1grid,GiniCoeff,'k','LineWidth',2);
hold on
plot(alpha1grid,GiniCoeffTarget*ones(alpha1GridSize,1),':k','LineWidth',2);
%%
xlabel('alpha_1')
ylabel('GiniCoeff')
subplot(3,1,3)
plot(alpha1grid,AvgHrs,'k','LineWidth',2);
hold on
plot(alpha1grid,AvgHrsTarget*ones(alpha1GridSize,1),':k','LineWidth',2);
xlabel('alpha_1')
ylabel('AvgHrs')
print(gcf,'-depsc2 ',[plotpath 'Calibration.eps'])
print(gcf,'-dpng ',[plotpath 'Calibration.png'])



clc
clear all
close all
SetParaStruc
x=fsolve(@(x) GetCalibration2(x), [1 1]);
theta_2=1;
gamma=x(1);
Y=x(2);
g=.12*Y;
psi=1/(1+gamma);
theta_1=3.3;
n1=1;
n2=1;
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
RMax=max(Rbar)*1.8;

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

alpha1GridSize=25;
alpha1grid=linspace(0.6, .95,alpha1GridSize);
tau1Target=.2;
GiniCoeffTarget=.38;
AvgHrsTarget=.23;
VarLogWagesTarget=.35;
c10guess=.5;
c20guess=.5/(7);
RGrid=[];
plotpath='Graphs/';

for i=1:alpha1GridSize
  alpha_1=alpha1grid(i);
alpha_2=1-alpha_1;
alpha_1=alpha_1*Para.n1;
alpha_2=alpha_2*Para.n2;
Para.alpha_1=alpha_1;
Para.alpha_2=alpha_2;
    
     [res]=GetSteadyStateMoments(Para,c10guess,c20guess)
     tau1(i)=res.tau1;
     GiniCoeff(i)=res.GiniCoeff;
     AvgHrs(i)=res.AvgHrs;
     exitflagv0(i)=res.exitflagv0;
     VarLogWages(i)=res.VarLogWages;
  end


figure()
subplot(3,1,1)
plot(alpha1grid,tau1,'k','LineWidth',2);
hold on
plot(alpha1grid,tau1Target*ones(alpha1GridSize,1),':k','LineWidth',2);
xlabel('alpha_1')
ylabel('Labor Tax')
subplot(3,1,2)
% plot(alpha1grid,GiniCoeff,'k','LineWidth',2);
% hold on
% plot(alpha1grid,GiniCoeffTarget*ones(alpha1GridSize,1),':k','LineWidth',2);
 plot(alpha1grid,VarLogWages,'k','LineWidth',2);
 hold on
 plot(alpha1grid,VarLogWagesTarget*ones(alpha1GridSize,1),':k','LineWidth',2);
%%
xlabel('alpha_1')
ylabel('VarLogWages')
subplot(3,1,3)
plot(alpha1grid,AvgHrs,'k','LineWidth',2);
hold on
plot(alpha1grid,AvgHrsTarget*ones(alpha1GridSize,1),':k','LineWidth',2);
xlabel('alpha_1')
ylabel('AvgHrs')
print(gcf,'-depsc2 ',[plotpath 'Calibration3.eps'])
print(gcf,'-dpng ',[plotpath 'Calibration3.png'])
