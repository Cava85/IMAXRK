% test the new idea to solve the nonlinear set of equations

%  min f(x)
%  C1(x) <= eps

clc, clear all, close all;

fun =@(x) norm(x);
con =@(x) x(2,:) - x(1,:).^2 - 0.5;
t =0:0.01:1;
tt = [t;t];
tt = [t;zeros(size(t))];
cc = -con(tt);
figure(1); hold on;
plot([0,1], [1,1], '-')
plot(t,cc, 'k-')