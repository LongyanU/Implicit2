clear
clc %%%%%%%
close all
% Elapsed time is 33.130840 seconds..
% Elapsed time is 29.353771 seconds.
% new implicit method
nt=403;    % number of time steps
isnap=20;    % snapshot sampling

nx=300;
nz=250;

v=ones(nz,nx)*2800;
% v(1:nz/2,:)=1500;




dx=20;  %calculate space increment
h=dx;
x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
dt=0.0025; % calculate time step from stability criterion
tau=dt;


f0=65;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^8*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=(diff(src))/dx^2;				% time derivative to obtain gaussian


zs=46;
xs=nx/2;

seis_record=zeros(nt,nx);
p=zeros([nz nx]); pnew=p; pold=p;

VxzAddedpoint=p;

r=v*dt/h;



coeff=zeros(nz,nx,5);
for ii=1:nz
    for jj=1:nx
        
        coeff(ii,jj,:)=[ 0.535641, 0.161028, -0.00732291, 0.00868726, 0.199816];
        
    end
end


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
tic
for it=1:nt-2,
    
    d2px11=Vx-circshift(Vx,[0 1]);
    d2px12=(circshift(Vx,[0 -1])-circshift(Vx,[0 2]));
    d2px13=(circshift(Vx,[0 -2])-circshift(Vx,[0 3]));
    %         d2px14=(circshift(Vx,[0 -3])-circshift(Vx,[0 4]));
    %                 d2px15=(circshift(Vx,[0 -4])-circshift(Vx,[0 5]));
    %                 d2px16=(circshift(Vx,[0 -5])-circshift(Vx,[0 6]));
    %                 d2px17=(circshift(Vx,[0 -6])-circshift(Vx,[0 7]));
    addPoint=(circshift(Vx,[-1 0])-circshift(Vx,[-1 1]));
    addPoint=addPoint+(circshift(Vx,[1 0])-circshift(Vx,[1 1]));
    d2px=coeff(:,:,1).*d2px11+coeff(:,:,2).*d2px12+coeff(:,:,3).*d2px13;
    d2px=d2px+addPoint.*coeff(:,:,4);
    
    d2pz11=Vz-circshift(Vz,[1 0]);
    d2pz12=(circshift(Vz,[-1 0])-circshift(Vz,[2 0]));
    d2pz13=(circshift(Vz,[-2 0])-circshift(Vz,[3 0]));
    %         d2pz14=(circshift(Vz,[-3 0])-circshift(Vz,[4 0]));
    %                 d2pz15=(circshift(Vz,[-4 0])-circshift(Vz,[5 0]));
    %                 d2pz16=(circshift(Vz,[-5 0])-circshift(Vz,[6 0]));
    %                 d2pz17=(circshift(Vz,[-6 0])-circshift(Vz,[7 0]));
    addPoint=(circshift(Vz,[0 -1])-circshift(Vz,[1 -1]));
    addPoint=addPoint+(circshift(Vz,[0 1]) -circshift(Vz,[1 1]));
    
    d2pz=coeff(:,:,1).*d2pz11+coeff(:,:,2).*d2pz12+coeff(:,:,3).*d2pz13;
    d2pz=d2pz+addPoint.*coeff(:,:,4);
    
    for j=1:nx,
        d2pzz(:,j)=A1{j}\d2pz(:,j);
    end
    for k=1:nz,
        d2pxx(k,:)=A2{k}\d2px(k,:)';
    end
    p=p-dt*v.^2.*(d2pxx+d2pzz)/h;
    
    p(zs,xs)= p(zs,xs)+src(it)*dt^2;
    
    %% time saved
    d2px1=circshift(p,[0 -1])-p;
    d2pz1=circshift(p,[-1])-p;
    
    Vx=Vx-dt*d2px1/h;
    Vz=Vz-dt*d2pz1/h;
    %% time saved
    [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,45,45,0.009);
    %     if these lines are commented, then the new
    %     method is two times faster
    %     if rem(it,isnap)== 0,
    %         imagesc(x,z,p), axis equal
    %         colormap gray
    %         xlabel('x'),ylabel('z')
    %         title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(p))))
    %         drawnow
    %     end
    
    if rem(it,isnap)== 0,
        imagesc(x,z,p), axis equal
        colormap gray
        xlabel('x'),ylabel('z')
        title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(p))))
        drawnow
    end
    
    
    
    seis_record(it,:)=p(zs,:);
    %     if these lines are commented, then the new
    %     method is two times faster
    %     if it==400
    %         pp1=p;
    %     end
    if it==400
        pp1=p;
    end
end


toc
save('Figure3bSimulation.mat')
figure;imagesc(pp1(45:end-45,45:end-45))
xlabel('x/dx')
ylabel('z/dz')
h = colorbar;