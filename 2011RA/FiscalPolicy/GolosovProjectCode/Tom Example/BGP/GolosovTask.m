SetParaStruc
alpha_1=.6
x=fsolve(@(x) GetCalibration(x), [1 .5 1]);
theta_2=min(x(1),1);
gamma=x(2);
Y=x(3);
g=.12*Y;
gamma = .42;
psi=1/(1+gamma);
theta_1=max(x(1),1);
n1=1;
n2=1;
alpha_2=1-alpha_1;
alpha_1=alpha_1*n1;
alpha_2=alpha_2*n2;
Para.alpha_1=alpha_1;
Para.alpha_2=alpha_2;
Para.psi=psi;
Para.g=ones(1,2).*.5;
Para.theta_1 =3.3;
Para.theta_2=1;
Para.btild_1=0;


%% Compute the undistorted FB
s_=1;

%% Build Grid for the state variables
% This setups up the functional space and the grid.
u2btildMin=-10;
u2btildMax=-u2btildMin;
u2btildGrid=linspace(u2btildMin,u2btildMax,Para.u2btildGridSize);

Para.u2bdiffGrid=u2btildGrid;
Para.u2btildLL=u2btildMin;
Para.u2btildUL=u2btildMax;
Rbar0=5;
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
Para.RMax=RMax;
Para.RMin=RMin;
%% Define the funtional space
V(1) = fundefn(Para.ApproxMethod,[Para.OrderOfAppx_u2btild Para.OrderOfApprx_R ] ,[u2btildMin RMin],[u2btildMax RMax]);
V(2) = V(1);
GridPoints=[Para.u2btildLL Para.u2btildUL;RMin RMax];
rowLabels = {'$x$','$R$'};
columnLabels = {'Lower Bound','Upper Bounds'};

    %% INITIALIZE THE COEFF
    %  This function computes c1,c2,l1,l2 and the value for an arbitrary x, R.
    % This section solves for V i.e the value function at the end of period 1
    % with g_t=g for all t >1. since the value function is static we need to
    % solve a equation in c_1 for each x,R. Th function getValueC1 does the job
    % by solving for the two roots of this equation and using the one that
    % supports the highest utility
    tic
    gTrue=Para.g;
    Para.g=mean(gTrue)*ones(2,1);
    for s_=1:Para.sSize
        n=1;
        if s_==1
            
            
            for u2btildctr=1:Para.u2btildGridSize
                for Rctr=1:Para.RGridSize
                    
                    u2btild_=u2btildGrid(u2btildctr);
                    R_=RGrid(Rctr);
                    %if R_>Rbar(u2btildctr)
                    x_state_(s_,n,:)=[u2btild_ R_ ];
                    % Solve for  c1
                    
                    c1_=max(getValueC1(u2btild_,R_,s_,Para ),.0001);
                    
                  
                    if c1_<.001
                        ExitFlagT(n)=0;
                    else
                        ExitFlagT(n)=1;
                    end
                    % compute c2
                    c2_=R_^(-1)*c1_;
                    TotalResources=(c1_*Para.n1+c2_*Para.n2+Para.g(s_));
                    FF=R_*Para.theta_2/Para.theta_1;
                    DenL2=Para.n1*Para.theta_1*FF+Para.theta_2*Para.n2;
                    l2_=(TotalResources-Para.n1*Para.theta_1+Para.n1*Para.theta_1*FF)/(DenL2);
                    l1_= 1-FF*(1-l2_);
                    u2btildPrime_=u2btild_;
                    V0(s_,n)=(Para.alpha_1*uBGP(c1_,l1_,Para.psi)+Para.alpha_2*uBGP(c2_,l2_,Para.psi))/(1-Para.beta);
                
                    xInit_0(s_,n,:)=[c1_ c2_ l1_ l2_ u2btildPrime_/(Para.psi*c2_^(-1)) R_ u2btildPrime_];
                    n=n+1;
                    %end
                end
            end
            c0(s_,:)=funfitxy(V(s_),squeeze(x_state_(s_,logical(ExitFlagT==1),:)),squeeze(V0(s_,logical(ExitFlagT==1)))' );
        else
            c0(s_,:)=c0(s_-1,:);
            V0(s_,:)=V0(s_-1,:);
            xInit_0(s_,:)=xInit_0(s_-1,:);
        end
    end
    disp('Number of points solved in initialization')
    sum(ExitFlagT)
    disp('Number of points solved out of a total of ')
    length(ExitFlagT)
    
    Para.g=gTrue;
    x_state=vertcat([squeeze(x_state_(1,:,:)) ones(length(x_state_),1)] ,[squeeze(x_state_(1,:,:)) 2*ones(length(x_state_),1)]);
    scatter(x_state(:,1),x_state(:,2))
    c=c0
 
    % slicing the state space for parfor loop later
    u2btild_slice=x_state(:,1) ;
    R_slice=x_state(:,2) ;
    s_slice=x_state(:,3) ;
    % This stores the values of the policy functions and multipliers that last
    % worked
    
    PolicyRulesWorked=[xInit_0(1,1,1) xInit_0(2,1,1) xInit_0(1,1,2)];
    
    % This stores the policy rules for each point in the state
    % space.
    PolicyRulesStore1=[squeeze(xInit_0(1,:,1))' squeeze(xInit_0(1,:,1))' ...
        squeeze(xInit_0(1,:,2))' squeeze(xInit_0(1,:,2))'...
        squeeze(xInit_0(1,:,3))' squeeze(xInit_0(1,:,3))' ...
        squeeze(xInit_0(1,:,4))' squeeze(xInit_0(1,:,4))' ....
        squeeze(xInit_0(1,:,5))' squeeze(xInit_0(1,:,5))' ....
        squeeze(xInit_0(1,:,6))' squeeze(xInit_0(1,:,6))' ....
        squeeze(xInit_0(1,:,7))' squeeze(xInit_0(1,:,7))' ....
        ];
    PolicyRulesStore2=[squeeze(xInit_0(2,:,1))' squeeze(xInit_0(2,:,1))' ...
        squeeze(xInit_0(2,:,2))' squeeze(xInit_0(2,:,2))'...
        squeeze(xInit_0(2,:,3))' squeeze(xInit_0(2,:,3))' ...
        squeeze(xInit_0(2,:,4))' squeeze(xInit_0(2,:,4))' ....
        squeeze(xInit_0(2,:,5))' squeeze(xInit_0(2,:,5))' ....
        squeeze(xInit_0(2,:,6))' squeeze(xInit_0(2,:,6))' ....
        squeeze(xInit_0(2,:,7))' squeeze(xInit_0(2,:,7))' ....
        ];
    PolicyRulesStore=vertcat(PolicyRulesStore1,PolicyRulesStore2);












% SOLVE THE T-0 PROBLEM given btild(-1)
btild_1=Para.btild_1;
disp('Computed V...Now solving V0(btild_1) where btild_1 is')
disp(btild_1)
% c1 and c2 solve
options=optimset('Display','off');
[x,~,exitflagv0,~,~] = fminunc(@(x)  getValue0(x, btild_1,1,Para,c,V),[ .5 .5*mean(Para.RGrid)^(-1)],options);
if ~(exitflagv0==1)
    disp('Optimization failed for V0 once ..trying with fmincon')
    opts = optimset('Algorithm', 'interior-point', 'Display','off', ...
        'GradObj','off','GradConstr','off',...
        'MaxIter',1000, ...
        'TolX', Para.ctol/10, 'TolFun', Para.ctol, 'TolCon', Para.ctol,'MaxTime',200);
    lb=[0.001 0.001];
    ub=[10 10];
    %[x,fval,exitflagv0,output,lambda]  =fmincon(@(x) getValue0(x, btild_1,1,Para,c,V),[ x ],[],[],[],[],lb,ub,[],opts);
    [x,~,exitflagv0,output,lambda]  =fmincon(@(x) getValue0(x, btild_1,1,Para,c,V),[ 1.8380 0.6526],[],[],[],[],lb,ub,[],opts);
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
ul20=(1-psi)/(1-l20);
ul10=(1-psi)/(1-l10);
uc20=psi/c20;
uc10=psi/c10;
tau0=1-(ul20/(theta_2*uc20))

 u2btild_=u2btildprime0;
                    R_=Rprime0;
                    c1=max(getValueC1(u2btild_,R_,s_,Para ),.0001);
                   
                    % compute c2
                    c2=R_^(-1)*c1;
                    TotalResources=(c1*Para.n1+c2_*Para.n2+Para.g(s_));
                    FF=R_*Para.theta_2/Para.theta_1;
                    DenL2=Para.n1*Para.theta_1*FF+Para.theta_2*Para.n2;
                    l2=(TotalResources-Para.n1*Para.theta_1+Para.n1*Para.theta_1*FF)/(DenL2);
                    l1= 1-FF*(1-l2);
                   ul2=(1-psi)/(1-l2);
ul1=(1-psi)/(1-l1);

uc1=psi/c1;
uc2=psi/c2;
tau=1-(ul2/(theta_2*uc2))
 Trans=c2-l2.*ul2./uc2;
    
  
    
     % Income
    AfterTaxWageIncome_Agent2=l2.*ul2./uc2+Trans;
    AfterTaxWageIncome_Agent1=l1.*ul1./uc1+Trans;
    % Gini Coeff
    GiniCoeff=(AfterTaxWageIncome_Agent2 +2*AfterTaxWageIncome_Agent1)./(AfterTaxWageIncome_Agent2+AfterTaxWageIncome_Agent1)-3/2;
   

AvgHrs=(l1+l2)/2;

Para.g(1)/(theta_1*l1+theta_2*l2)
