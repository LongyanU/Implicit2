clear;
clc
% close all

global v ratio M h tau

v=1500;

h=20;
tau=0.0025;
M=3;
r=v*tau/h;
ratio=0.815;

x0=0.001*ones(1,M+2);%系数的初值是0

options = optimset('Algorithm','levenberg-marquardt','TolFun',10^-10,'TolX',10^-10,'MaxFunEvals',2000,'MaxIter',200);
aa=-2*ones(1,M+1);
bb=2*ones(1,M+1);

tic
[x,resnorm] = lsqnonlin(@myfun2,x0,aa,bb,options)   % Invoke optimizer
toc

x=real(x);
k=linspace(0,pi/h,1000); %

F=zeros(5,1000);

for i=1:5  % F的意思是5个角度,30个波数点，在这5个角度,30个波数点 满足时间-空间域频散关系。实际就是在这些特定的角度和波数点，公式的误差最小。
    xita=(i-1)*pi/16;
    
    for m=1:M
        F(i,:)=2*x(m)*( cos(m*k*h *cos(xita) )-cos((m-1)*k*h*cos(xita)) )+ F(i,:);
    end
    F(i,:)=F(i,:)+x(M+1)* 4*cos(k*sin(xita)*h).*(cos(k*h*cos(xita))-1);
    F(i,:)=F(i,:).*(1+2*x(M+2)*(cos(k*h*sin(xita))-1));
    
    
    temp=0;
    for m=1:M
        temp=2*x(m)*( cos(m*k*h *sin(xita) )-cos((m-1)*k*h*sin(xita)) )+ temp;
    end
    temp=temp+x(M+1)*4*cos(k*cos(xita)*h).*(cos(k*h*sin(xita))-1); 
    temp=temp.*(1+2*x(M+2)*(cos(k*h*cos(xita))-1));
    
    F(i,:)=F(i,:)+temp;
    F(i,:)=F(i,:)*r^2;
    
    F(i,:)=1/2*F(i,:)./  ( (1+2*x(M+2)*(cos(k*h*cos(xita))-1)) .*(1+2*x(M+2)*(cos(k*h*sin(xita))-1)) )  +1;
    F(i,:)=(acos(F(i,:))./ (tau*k*v));
    a1=(h/v*(1./F(i,:)-1));
    if (i==1)
        figure;plot(v*k*h/(2*pi*h),a1,'m','linewidth',2)
        hold on
    elseif i==2
        plot(v*k*h/(2*pi*h),a1,'r--','linewidth',2)
    elseif i==3
        plot(v*k*h/(2*pi*h),a1,'c:','linewidth',2)
    elseif i==4
        plot(v*k*h/(2*pi*h),a1,'k-.','linewidth',2)
    else
        plot(v*k*h/(2*pi*h),a1,'b','linewidth',2)
    end
end
grid on

axis([0 v/(2*h) -3*10^-5 7*10^-5])
xlabel('f(hz)')
legend('\theta=0', '\theta=π/16','\theta=2π/16','\theta=3π/16','\theta=4π/16')

ylabel('\delta (\theta)')
