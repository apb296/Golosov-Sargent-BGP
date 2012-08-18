function [ taul,tauh,Tl,Th,btildl,btildh,Rl,Rh]...
    = getTaxSchedule( cl,ch,ll,lh,btild_1,Para )
%GETTAXSCHEDULE Computes tax schedule from consumption labor choices
%   Detailed explanation goes here
theta_1 = Para.theta_1;
theta_2 = Para.theta_2;
gamma = Para.gamma;
sigma = Para.sigma;
beta=Para.beta;
taul = zeros(1,3);
tauh = zeros(1,3);
Rl = zeros(1,3); Rl(1) = cl(2,1)^(sigma)/beta;
Rh = zeros(1,3); Rh(1) = ch(2,1)^(sigma)/beta;
btildl = zeros(1,4); btildl(1) = btild_1;
btildh = zeros(1,4); btildh(1) = btild_1;
Tl = zeros(1,3);
Th = zeros(1,3);
for t= 1:3
    %get taxes
    taul(t) = 1-(ll(1,t)^gamma)/(theta_1*cl(1,t)^(-sigma));
    tauh(t) = 1-(lh(1,t)^gamma)/(theta_1*ch(1,t)^(-sigma));
    %get prices
    if(t > 1)
        Rl(t) = cl(1,t-1)^(-sigma)/(beta*cl(1,t)^(-sigma));
        Rh(t) = ch(1,t-1)^(-sigma)/(beta*ch(1,t)^(-sigma));
    end
    %get asset difference
    btildl(t+1) = Rl(t)*btildl(t)- (cl(2,t)-cl(1,t)-(1-taul(t))*theta_2*ll(2,t)...
        +(1-taul(t))*theta_1*ll(1,t));
    btildh(t+1) = Rh(t)*btildh(t)- (ch(2,t)-ch(1,t)-(1-tauh(t))*theta_2*lh(2,t)...
        +(1-tauh(t))*theta_1*lh(1,t));
    %get Transfers assuming b_{2,t} =0, all governmnet debt is held by
    %agent1
    Tl(t) = cl(2,t)-(1-taul(t))*theta_2*ll(2,t);
    %Th(t) = cl(2,t)-(1-taul(t))*theta_2*lh(2,t);
    Th(t) = ch(2,t)-(1-tauh(t))*theta_2*lh(2,t);
end


end

