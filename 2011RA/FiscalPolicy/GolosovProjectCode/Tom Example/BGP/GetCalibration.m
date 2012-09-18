function [ res] = GetCalibration (x)
theta_2=x(1);
gamma=x(2);
Y=x(3);

res(1)= (1+1.2*gamma)*Y-1-theta_2;
res(2)=7.3*(0.8-0.08*gamma*Y)+0.08*gamma*Y-0.8*theta_2;
res(3)=0.23*(1+gamma)*theta_2+0.05*gamma*Y*(1+theta_2)-theta_2;


end

