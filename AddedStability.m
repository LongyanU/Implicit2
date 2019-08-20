clear;
clc
close all

global v ratio M h tau



h=10;
tau=0.001;
M=3;

ratio=0.815;

x0=0.001*ones(1,M+2);%系数的初值是0

options = optimset('Algorithm','levenberg-marquardt','TolFun',10^-10,'TolX',10^-10,'MaxFunEvals',2000,'MaxIter',200);
ii=1;

aa=-2*ones(1,M+1);
bb=2*ones(1,M+1);




for v=100:100:10000
    r=v*tau/h
    [x,resnorm] = lsqnonlin(@myfun7,x0,aa,bb,options)    % Invoke optimizer
    
    temp=0;
    for i=1:M
        temp=temp+2*x(i)*(-1)^(i-1);
    end
    
    temp=temp-4*x(M+1); %%%%%%%%%%%%%%%%
    temp=2*(1-4*x(M+2))^2/temp^2;
    
    r1(ii)=r;
    s1(ii)=sqrt(temp);
    ii=ii+1;
    
end

figure;plot(r1(1:90),real(s1(1:90)),'k');

x0=0.001*ones(1,M+2);%系数的初值是0

options = optimset('Algorithm','levenberg-marquardt','TolFun',10^-10,'TolX',10^-10,'MaxFunEvals',2000,'MaxIter',200);
aa=-2*ones(1,M+1);
bb=2*ones(1,M+1);

ii=1;
for v=100:100:10000
    r=v*tau/h
    [x,resnorm] = lsqnonlin(@myfun2,x0,aa,bb,options);   % Invoke optimizer
    
    temp=0;
    for i=1:M
        temp=temp+2*x(i)*(-1)^(i-1);
    end
    
    temp=temp-4*x(M+1); %%%%%%%%%%%%%%%%
    temp=(1-4*x(M+2))/temp;
    
    r1(ii)=r;
    s2(ii)=sqrt(temp);
    ii=ii+1;
    
end

hold on;plot(r1(1:90),real(s2(1:90)),'k');
hold on;plot(r1(1:90),r1(1:90),'r');


legend('r','previous implicit FD scheme', 'implicit-explicit FD scheme')
xlabel('r')
ylabel('Stability')
grid on

save Stability.mat