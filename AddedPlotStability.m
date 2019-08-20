clear;
clc
close all;

load('Stability.mat')

figure;plot(r1(1:90),real(s1(1:90)),'b','linewidth',1.5);
hold on;plot(r1(1:90),real(s2(1:90)),'k','linewidth',1.5);
hold on;plot(r1(1:90),r1(1:90),'r','linewidth',1.5);

legend('previous implicit FD scheme', 'implicit-explicit FD scheme','r')
xlabel('r')
ylabel('s')
grid on
