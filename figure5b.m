%May 7,2016
% authors Liang wenquan, Wang Yanfei, Yang Changchun

%May 7,2016
% Elapsed time is 159.406186 seconds.

%  180.344233 April 23, 2018
% 时间已过 179.028720 秒。

% 时间已过 182.199654 秒。 13 September,2018
clear
clc %%%%%%%
close all
nt=1602;    % number of time steps
eps=.6;     % stability
isnap=5;    % snapshot sampling
load('vv')

c1=flipud(c);

v=c1;
nx=800;
nx=nx+45*2;
nz=475;
nz=nz+45*2;

vv=zeros(nz,nx);
for ii=1:nz-90
    for jj=1:nx-90
        vv(ii+45,jj+45)=v(ii,jj);
    end
end

for ii=1:nz-90  %%left
    for jj=1:45
        vv(ii+45,jj)=v(ii,1);
    end
end

for ii=1:nz-90  %%right
    for jj=nx-45:nx
        vv(ii+45,jj)=v(ii,800);
    end
end


for ii=1:45  %%top
    for jj=1:nx
        vv(ii,jj)=vv(46,jj);
    end
end

for ii=nz-44:nz  %%bottom
    for jj=1:nx
        vv(ii,jj)=vv(nz-45,jj);
    end
end


clear v
v=vv;
% v=round(v);


dx=20;  %calculate space increment
h=dx;
x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
dt=0.0025; % calculate time step from stability criterion
tau=dt;


f0=40;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^8*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=(diff(src))/dx^2;				% time derivative to obtain gaussian


zs=46;
xs=600-150+25;

seis_record=zeros(nt,nx);
p=zeros([nz nx]); pnew=p; pold=p;

VxzAddedpoint=p;

r=v*dt/h;

% load ('coeffJune16.mat')
load ('coeffApril19.mat')
coeff=zeros(nz,nx,5);
for ii=1:nz
    for jj=1:nx
        
        coeff(ii,jj,:)=April19(floor((v(ii,jj)-1486))+1,:);
        
    end
end
tic


b=1-2*coeff(:,:,end);
a = coeff(:,:,end);
cc = a;

for j=1:nx,
    A1{j} = (gallery('tridiag',a(1:nz-1,j),b(:,j),cc(1:nz-1,j)));
end


% a = coeff(end)* ones(nx-1,1);
% cc = a;
A2=cell(nz,1);
for k=1:nz,
    A2 {k}= (gallery('tridiag',a(k,1:nx-1),b(k,:),cc(k,1:nx-1) ));      %%解三对角阵，直接matlab解了
end


Vx=zeros(nz,nx);
Vz=zeros(nz,nx);
d2pzz=p;d2pxx=p;

taper=ones(nz,nx);
for i=1:43
    for j=1:nx
        taper(i,j)=0.5-0.5*cos(pi*(i-1)/(43-1));
        taper(nz-i+1,j)=taper(i,j);
    end
end
for i=1:nz
    for j=1:43
        taper(i,j)=taper(i,j)*(0.5-0.5*cos(pi*(j-1)/(43-1)));
        taper(i,nx-j+1)=taper(i,j);
    end
end

for it=1:nt-2,
%     it;
    
    d2px11=Vx-circshift(Vx,[0 1]);
    d2px12=(circshift(Vx,[0 -1])-circshift(Vx,[0 2]));
    d2px13=(circshift(Vx,[0 -2])-circshift(Vx,[0 3]));
    
%                     d2px14=(circshift(Vx,[0 -3])-circshift(Vx,[0 4]));
%                     d2px15=(circshift(Vx,[0 -4])-circshift(Vx,[0 5]));
    %                 d2px16=(circshift(Vx,[0 -5])-circshift(Vx,[0 6]));
    %                 d2px17=(circshift(Vx,[0 -6])-circshift(Vx,[0 7]));
    
    addPoint=(circshift(Vx,[-1 0])-circshift(Vx,[-1 1]));
    addPoint=addPoint+(circshift(Vx,[1 0])-circshift(Vx,[1 1]));
    d2px=coeff(:,:,1).*d2px11+coeff(:,:,2).*d2px12+coeff(:,:,3).*d2px13 ;
    d2px=d2px+addPoint.*coeff(:,:,4);
    
    d2pz11=Vz-circshift(Vz,[1 0]);
    d2pz12=(circshift(Vz,[-1 0])-circshift(Vz,[2 0]));
    d2pz13=(circshift(Vz,[-2 0])-circshift(Vz,[3 0]));
%                     d2pz14=(circshift(Vz,[-3 0])-circshift(Vz,[4 0]));
%                     d2pz15=(circshift(Vz,[-4 0])-circshift(Vz,[5 0]));
%                     d2pz16=(circshift(Vz,[-5 0])-circshift(Vz,[6 0]));
    %                 d2pz17=(circshift(Vz,[-6 0])-circshift(Vz,[7 0]));
    addPoint=(circshift(Vz,[0 -1])-circshift(Vz,[1 -1]));
    addPoint=addPoint+(circshift(Vz,[0 1]) -circshift(Vz,[1 1]));
    
    
    d2pz=coeff(:,:,1).*d2pz11+coeff(:,:,2).*d2pz12+coeff(:,:,3).*d2pz13...
                             ;
    d2pz=d2pz+addPoint.*coeff(:,:,4);
    
    for j=1:nx,
        d2pzz(:,j)=A1{j}\d2pz(:,j);   %%解三对角阵，直接matlab解了
    end
    for k=1:nz,
        d2pxx(k,:)=A2{k}\d2px(k,:)';       %%解三对角阵，直接matlab解了
    end
    p=p-dt*v.^2.*(d2pxx+d2pzz)/h;
    
    p(zs,xs)= p(zs,xs)+src(it)*dt^2;
    
    d2px=(circshift(p,[0 -1])-circshift(p,[0 0]));
    %     d2px2=(circshift(p,[0 -2])-circshift(p,[0 1]));
    %     d2px3=(circshift(p,[0 -3])-circshift(p,[0 2]));
    %
    %     addPoint=(circshift(p,[-1 -1])-circshift(p,[-1 0]));
    %     addPoint=addPoint+(circshift(p,[1 -1])-circshift(p,[1 0]));
    
    %     d2px=coeff(:,:,1).*d2px1+coeff(:,:,2).*d2px2+coeff(:,:,3).*d2px3;
    %     d2px=d2px+addPoint.*coeff(:,:,4);
    
    d2pz=(circshift(p,[-1])-circshift(p,[0]));
    %     d2pz2=(circshift(p,[-2])-circshift(p,[1]));
    %     d2pz3=(circshift(p,[-3])-circshift(p,[2]));
    %
    %     addPoint=(circshift(p,[-1 -1])-circshift(p,[0 -1]));
    %     addPoint=addPoint+(circshift(p,[-1 1]) -circshift(p,[0 1]));
    %
    %     d2pz=coeff(:,:,1).*d2pz1+coeff(:,:,2).*d2pz2+coeff(:,:,3).*d2pz3;
    %     d2pz=d2pz+addPoint.*coeff(:,:,4);
    
    %     for j=1:nx,
    %         d2pz(:,j)=A1{j}\d2pz1(:,j);        %%解三对角阵，直接matlab解了
    %     end
    %     for k=1:nz,
    %         d2px(k,:)=A2{k}\d2px1(k,:)';       %%解三对角阵，直接matlab解了
    %     end
    
    Vx=Vx-dt*d2px/h;
    Vz=Vz-dt*d2pz/h;
    
    [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,45,45,0.009);
    
%     Vx=Vx.*taper;
%     Vz=Vz.*taper;
    
    if rem(it,isnap)== 0,
        imagesc(x,z,p,[-0.0007 0.0007]), axis equal
        colormap gray
        xlabel('x'),ylabel('z')
        title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(p))))
        drawnow
    end
    
    
    
    seis_record(it,:)=p(zs,:);
    seis_record2(it,:)=Vx(zs,:);
%     if it==500
%         pp1=p;
%     elseif it==1000
%         pp2=p;
%     elseif it==1500
%         pp3=p;
%     elseif it==2000
%         pp4=p;
%     elseif it==2500
%         pp5=p;
%     end
    
end


toc
save('Figure5bSaltModel.mat')
figure;imagesc(v(45:end-45,45:end-45))
xlabel('x/dx')
ylabel('z/dz')
h = colorbar;
set(get(h,'title'),'string','m/s');
hold on ;plot(xs,zs-45,'*r')