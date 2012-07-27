% Inputs - xInit, state variables - u2btild,,R,s_  coeff, value
% function, para
function [PolicyRules, V_new,exitflag]=OptimUnc(u2bdiff,RR,s,c,VV,xInit,Para)
global V Vcoef R u2btild Par s_ 

%Get the initial guess for the uconstraint problem. With the simplification
%we need only c1_1,c1_2and c2_1

xInit=xInit(1:3);
Para.theta=[Para.theta_1 Para.theta_2];
Para.alpha=[Para.alpha_1 Para.alpha_2];
Par=Para;
u2btild=u2bdiff;
R=RR;
Vcoef{1}=c(1,:)';
Vcoef{2}=c(2,:)';
V=VV;
s_=s;
u2btildLL=Para.u2btildLL;
u2btildUL=Para.u2btildUL;
n1=Para.n1;
n2=Para.n2;
ctol=Para.ctol;

%% Now solve the unconstraint problem FOC using NAG
% use the last solution
         options = optimset('GradObj','off'); % indicate gradient is provided 

   [x ,fval,exitflag] = fminunc(@(x) ValueFunction(x),xInit,options);

psi= Par.psi;
beta =  Par.beta;
P = Par.P;
theta_1 = Par.theta(1);
theta_2 = Par.theta(2);
g = Par.g;
alpha = Par.alpha;

sigma = 1;
    frac = (R*P(s_,1)*x(1)^(-sigma)+R*P(s_,2)*x(2)^(-sigma)-P(s_,1)*x(3)^(-sigma))...
        /( P(s_,2) );
     c1_1=x(1);
     c1_2=x(2);
     c2_1=x(3);
     
     
     
    %compute components from unconstrained guess
frac = (R*P(s_,1)*c1_1^(-sigma)+R*P(s_,2)*c1_2^(-sigma)-P(s_,1)*c2_1^(-sigma))...
        /( P(s_,2) );
    c2_2 = frac^(-1/sigma);
    
     %Compute l1 form formula
    l1_1den = n1*theta_1+n2*c2_1*theta_1/c1_1;
    l1_1num = (n1*c1_1+n2*c2_1+g(1) + n2*(c2_1*theta_1-c1_1*theta_2)/c1_1);
    l1(1) = l1_1num/l1_1den;
    l1_2den = n1*theta_1+n2*c2_2*theta_1/c1_2;
    l1_2num = (n1*c1_2+n2*c2_2+g(2) + n2*(c2_2*theta_1-c1_2*theta_2)/c1_2);
    l1(2) = l1_2num/l1_2den;
      %compute l2 from formula
    l2_1den = n2*theta_2+n1*c1_1*theta_2/c2_1;
    l2_1num = n1*c1_1+n2*c2_1+g(1)+n1*(c1_1*theta_2-c2_1*theta_1)/c2_1;
    l2(1) = l2_1num/l2_1den;
    l2_2den = n2*theta_2+n1*c1_2*theta_2/c2_2;
    l2_2num = n1*c1_2+n2*c2_2+g(2)+n1*(c1_2*theta_2-c2_2*theta_1)/c2_2;
    l2(2) = l2_2num/l2_2den;
    %get expected value of marginal utility of agent 2
    Eu2 = P(s_,1)*c2_1^(-1)+P(s_,2)*c2_2^(-1);
    
    %compute btildeprime from formula
    btildprime(1) = u2btild/(beta*psi*Eu2)...
        +c1_1-c2_1-(1-psi)*c1_1*l1(1)/(psi*(1-l1(1)))+(1-psi)*c2_1*l2(1)/(psi*(1-l2(1)));

    %Compute btildprime(2) from formula
    btildprime(2) = u2btild/(beta*psi*Eu2)...
        +c1_2-c2_2-(1-psi)*c1_2*l1(2)/(psi*(1-l1(2)))+(1-psi)*c2_2*l2(2)/(psi*(1-l2(2)));

u2btildprime=psi*[c2_1^(-1) c2_2^(-1)].*btildprime;


X(1,:) = [psi*c2_1^(-1)*btildprime(1),c2_1^(-1)/c1_1^(-1)];%state next period
X(2,:) = [psi*c2_2^(-1)*btildprime(2),c2_2^(-1)/c1_2^(-1)];%state next period
%compute objective
if ~isreal(X)
    X=real(X);
end
Vobj = P(s_,1)*(alpha(1)*uBGP(c1_1,l1(1),psi)+alpha(2)*uBGP(c2_1,l2(1),psi)...
    +beta*funeval(Vcoef{1},V(1),X(1,:)));

Vobj = Vobj + P(s_,2)*(alpha(1)*uBGP(c1_2,l1(2),psi)+alpha(2)*uBGP(c2_2,l2(2),psi)...
    +beta*funeval(Vcoef{2},V(2),X(2,:)));

V_new=Vobj;
PolicyRules=[c1_1 c1_2 c2_1 c2_2 l1(1) l1(2) l2(1) l2(2) btildprime c2_1^(-1)/c1_1^(-1) c2_2^(-1)/c1_2^(-1) u2btildprime(1) u2btildprime(2)];



end

% Knitro verification
function [errValue,errMultiplier]=CompareFOCwithOptimizer(u2btild,R,s_,Vcoef,V,Para,PolicyRulesFOC,~,MuL,MuU)
opts = optimset('Algorithm', 'interior-point', 'Display','off', ...
    'GradObj','off','GradConstr','off',...
    'MaxIter',1000, ...
    'TolX', Para.ctol/10, 'TolFun', Para.ctol, 'TolCon', Para.ctol,'MaxTime',200);
lb=[0.001 0.001 0.001 .001 Para.g(1)/Para.theta_1 Para.g(2)/Para.theta_1 Para.u2btildLL Para.u2btildLL];
ub=[Inf Inf Inf Inf Inf Inf Para.u2btildUL Para.u2btildUL];

[PolicyRules,fval,exitflag,output,lambda]  =fmincon(@(x) ValueFunction(x,u2btild,R,Vcoef,s_,V,Para),PolicyRulesFOC,[],[],[],[],lb,ub,@(x) NonLinearConstraints(x,u2btild,R,s_,Para),opts);
if exitflag==1
    disp('Percentage error in value')
errValue=(ValueFunction(PolicyRulesFOC,u2btild,R,Vcoef,s_,V,Para)-fval)/ValueFunction(PolicyRulesFOC,u2btild,R,Vcoef,s_,V,Para);
disp(errValue)
errMultiplier=[-sum(MuL)-sum(lambda.lower) sum(MuU)-sum(lambda.upper)];
else
    errValue=NaN;
    errMultiplier=NaN;
end
end
