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
src=diff(diff(src))/dx^2;				% time derivative to obtain gaussian


zs=46;
xs=nx/2;

seis_record=zeros(nt,nx);
p=zeros([nz nx]); pnew=p; pold=p;



coeff=zeros(nz,nx,5);
for ii=1:nz
    for jj=1:nx
        
        coeff(ii,jj,:)=[ 0.535641, 0.161028, -0.00732291, 0.00868726, 0.199816];
        
    end
end

coeff2=cat(3, -coeff(:,:,1)*2  ,coeff(:,:,1)-coeff(:,:,2) ,coeff(:,:,2)-coeff(:,:,3), coeff(:,:,3),coeff(:,:,4) );



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

tic
for it=1:nt-2,
    
    
    
    
    d2pz11=(circshift(p,[ -1])+circshift(p,[ 1]));
    d2pz12=(circshift(p,[ -2])+circshift(p,[ 2]));
    d2pz13=(circshift(p,[ -3])+circshift(p,[ 3]));
    %    d2pzAddedpoints= circshift(p,[-1 -1])+circshift(p,[1  -1])-2*circshift(p,[0  -1]) ...
    %        + circshift(p,[-1 1])+circshift(p,[1  1])-2*circshift(p,[0  1]);
    
    d2px11=(circshift(p,[0  -1])+circshift(p,[0  1]));
    d2px12=(circshift(p,[0  -2])+circshift(p,[0  2]));
    d2px13=(circshift(p,[0  -3])+circshift(p,[0  3]));
    %     d2pxAddedpoints=circshift(p,[-1 -1])+circshift(p,[ -1 1])-2*circshift(p,[-1 0 ]) ...
    %         + circshift(p,[1 -1])+circshift(p,[1  1])-2*circshift(p,[1 0 ]);
    
    commonPoints=circshift(p,[-1 -1])+circshift(p,[1  1]) +circshift(p,[1 -1])+circshift(p,[-1 1]);
    d2pzAddedpoints= -2*d2px11 +commonPoints;
    d2pxAddedpoints=-2*d2pz11+commonPoints;
    
    d2px=coeff2(:,:,1).*p+coeff2(:,:,2).*d2px11+coeff2(:,:,3).*d2px12+coeff2(:,:,4).*d2px13+coeff2(:,:,5).*d2pxAddedpoints;
    d2pz=coeff2(:,:,1).*p+coeff2(:,:,2).*d2pz11+coeff2(:,:,3).*d2pz12+coeff2(:,:,4).*d2pz13 +coeff2(:,:,5).*d2pzAddedpoints;
    
    for j=1:nx,
        d2pzz(:,j)=A1{j}\d2pz(:,j);   %%解三对角阵，直接matlab解了
    end
    for k=1:nz,
        d2pxx(k,:)=A2{k}\d2px(k,:)';       %%解三对角阵，直接matlab解了
    end
    
    pnew=2*p-pold+v.*v.*(d2pzz +d2pxx)*dt^2/dx^2;
    pnew(zs,xs)=pnew(zs,xs)+src(it)*dt^2;
    pold=p;											% time lev(k,j)els
    p=pnew;
    
    [p,pold]=spongeABC(p,pold,nx,nz,45,45,0.009);
    
    if rem(it,isnap)== 0,
        it
        imagesc(x,z,p,[-0.0007 0.0007]), axis equal
        colormap gray
        xlabel('x'),ylabel('z')
        title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(p))))
        drawnow
    end
    
      if it==400
        pp1=p;
    end
    
end


toc
save('Figure3cSimulation.mat')
figure;imagesc(pp1(45:end-45,45:end-45))
xlabel('x/dx')
ylabel('z/dz')
h = colorbar;