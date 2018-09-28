 clear;
 clc
 close all

load('Figure3aSimulation.mat')
figure;imagesc(pp1(45:end-45,45:end-45),[-12*10^-3 6.6*10^-3])
pp11=pp1(45:end-45,45:end-45);
xlabel('x/dx')
ylabel('z/dz')
h = colorbar;
colormap gray


load('Figure3bSimulation.mat')
figure;imagesc(pp1(45:end-45,45:end-45),[-12*10^-3 6.6*10^-3])
pp22=pp1(45:end-45,45:end-45);
xlabel('x/dx')
ylabel('z/dz')
h = colorbar;
colormap gray

load('Figure3cSimulation.mat')
figure;imagesc(pp1(45:end-45,45:end-45),[-12*10^-3 6.6*10^-3])
pp22=pp1(45:end-45,45:end-45);
xlabel('x/dx')
ylabel('z/dz')
h = colorbar;
colormap gray

% figure;imagesc(pp11-pp22,[-12*10^-3 6.6*10^-3]);
% h = colorbar;
% xlabel('x/dx')
% ylabel('z/dz')
% colormap gray
% axis([])