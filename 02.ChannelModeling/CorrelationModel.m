function [y] = CorrelationModel(tau,mu,sigma,Ts)
t       = (0:1:1000)*Ts ;
X       = mu + sigma.*randn(1,1e5) ; %random(makedist('normal',mu,sigma),1,1e5) ;
h       = 2^0.75.*sqrt(1/tau).*besseli(0,t/tau)/pi^0.75 ;
Y       = filter(h,1,X) ;
y       = mu + (Y - mean(Y))*sigma/std(Y) ;
y       = y(1,1e5) ;
end