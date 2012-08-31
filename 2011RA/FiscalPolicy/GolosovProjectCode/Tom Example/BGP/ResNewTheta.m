function [ res] = ResNewTheta( ThetaGuess,ThetaSpread, AverageOutput,Para)
% This function computes the residual 
theta_1=ThetaGuess(1);
theta_2=ThetaGuess(2);
res(1)=(theta_1-theta_2)-ThetaSpread;
Para.theta_1=theta_1;
Para.theta_2=theta_2;
[c1FB c2FB l1FB l2FB yFB g_yFB_h Agent1WageShareFB_h]=getFB(Para,2);
Output(1)=c1FB*Para.n1+c2FB*Para.n2+Para.g(1);
[c1FB c2FB l1FB l2FB yFB g_yFB_l Agent1WageShareFB_l]=getFB(Para,1);
Output(2)=c1FB*Para.n1+c2FB*Para.n2+Para.g(2);
res(2)=sum(Output)/2-AverageOutput;
end

