function [ res] = GetCalibration2 (x)
theta_2=3.3;
gamma=x(1);
Y=x(2);
res(1)= (1+1.2*gamma)*Y-1-theta_2;
res(2)=0.23*(1+gamma)*theta_2+0.05*gamma*Y*(1+theta_2)-theta_2;
end

