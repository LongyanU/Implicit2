function F = myfun7(x)

global v ratio M h tau;

r=v*tau/h;

k=linspace(1*pi/(1000),ratio*pi/h,30); %
F=zeros(5,30);

for i=1:5  % F的意思是5个角度,30个波数点，在这5个角度,30个波数点 满足时间-空间域频散关系。实际就是在这些特定的角度和波数点，公式的误差最小。
    xita=(i)*pi/16;
    
    for m=1:M
        F(i,:)=2*x(m)*( sin((m-0.5)*k*h *cos(xita) ))+ F(i,:);
    end
    F(i,:)=F(i,:)+4*x(M+1)*cos(k*sin(xita)*h).*(sin( k*h/2*cos(xita) ) );
    F(i,:)=F(i,:).*(1+2*x(M+2)*(cos(k*h*sin(xita))-1));
    F(i,:)= F(i,:).^2;
    
    temp=0;
    for m=1:M
        temp=2*x(m)*( sin( (m-0.5)*k*h *sin(xita) ) )+ temp;
    end
    temp=temp+4*x(M+1)*cos(k*cos(xita)*h).*(sin( k*h/2*sin(xita) ) );
    temp=temp.*(1+2*x(M+2)*(cos(k*h*cos(xita))-1));
    
    F(i,:)=F(i,:)+temp.^2;
    F(i,:)=F(i,:)*r^2;
    
    F(i,:)=-1/2*F(i,:)./  ( (1+2*x(M+2)*(cos(k*h*cos(xita))-1)).^2.*(1+2*x(M+2)*(cos(k*h*sin(xita))-1)).^2 )  +1;
    F(i,:)=(acos(F(i,:))./ (tau*k*v));
    
%     F(i,:)=F(i,:)- (1+2*x(M+2)*(cos(k*h*cos(xita))-1)).^2.*(1+2*x(M+2)*(cos(k*h*sin(xita))-1)).^2 .* (2*cos(k*v*tau)-2);
end
F=F-1;
% F=F(1,:)*F(1,:)'+F(2,:)*F(2,:)'+F(3,:)*F(3,:)'+F(4,:)*F(4,:)'+F(5,:)*F(5,:)';