function [ FUN,CON ] = IMAXRK_Solver3
% clear all
% close all
% clc


%% Define variables

syms bI1 bI2 bI3 bI4 bI5
syms bE1 bE2 bE3 bE4
syms     c2  c3  c4

bi = [bI1, bI2, bI3, bI4, bI5];
be = [bE1, bE2, bE3, bE4, 0];

Ai = [0, 0, 0, 0, 0; bI1, c2-bI1, 0, 0, 0;
      bI1, bI2, c3-bI1-bI2, 0, 0; bI1, bI2, bI3, c4-bI1-bI2-bI3, 0;
      bI1, bI2, bI3, bI4, bI5];

Ae = [0, 0, 0, 0, 0; c2, 0, 0, 0, 0; bE1, c3-bE1, 0, 0, 0; 
      bE1, bE2, c4-bE1-bE2, 0, 0; bE1, bE2, bE3, bE4, 0];

c = diag([0 c2 c3 c4 1]);

e = ones(5,1);


%% Define constraints

t11i = bi*e - 1;
t11e = be*e - 1; % Eq1

t21i = bi*c*e - 1/2;
t21e = be*c*e - 1/2; % Eq2

% t31i = bi*c*c*e/2 - 1/6;
t31e = be*c*c*e/2 - 1/6; % Eq3

t32ee = be*Ae*c*e - 1/6; % Eq4 -> 2nd order
t32ei = be*Ai*c*e - 1/6; % Eq
t32ie = bi*Ae*c*e - 1/6;
t32ii = bi*Ai*c*e - 1/6; % Eq9 --> bI5


%% Define objective functions

t44eee = be*Ae*Ae*c*e - 1/24;
Lstab = - (bI1 * (bI1 + bI2 - c2) * (bI1 + bI2 + bI3 - c3) * (bI1 + bI2 + bI3 + bI4 - c4)) / ...
        (bI5 * (bI1 - c2) * (bI1 + bI2 - c3) * (bI1 + bI2 + bI3 - c4));

t42ei = be*c*Ai*c*e - 3/24;
t42ee = be*c*Ae*c*e - 3/24;
t43ie = bi*Ae*c*c*e/2 - 1/24;
t43ee = be*Ae*c*c*e/2 - 1/24;
t44iii = bi*Ai*Ai*c*e - 1/24;
t44iie = bi*Ai*Ae*c*e - 1/24;
t44iei = bi*Ae*Ai*c*e - 1/24;
t44iee = bi*Ae*Ae*c*e - 1/24;
t44eii = be*Ai*Ai*c*e - 1/24;
t44eie = be*Ai*Ae*c*e - 1/24;
t44eei = be*Ae*Ai*c*e - 1/24;

tau4 = sqrt(t42ei^2 + t42ee^2 + t43ie^2 + t43ee^2 + t44iii^2 + t44iie^2 + ...
       t44iei^2 + t44iee^2 + t44eii^2 + t44eie^2 + t44eei^2 + t44eee^2);
% keyboard
   
%% Solve nonlinear system

SolC = [1/2, 9/10, 7/10];

t11eC = subs(t11e, [c2, c3, c4], SolC);
t21eC = subs(t21e, [c2, c3, c4], SolC);
t31eC = subs(t31e, [c2, c3, c4], SolC);
SolE1 = solve(t11eC, t21eC, t31eC, bE1, bE2, bE3);

t32eeEC = subs(t32ee, [c2, c3, c4], SolC);
t32eeEC = subs(t32eeEC, [bE1, bE2, bE3], [SolE1.bE1, SolE1.bE2, SolE1.bE3]);
SolE2 = solve(t32eeEC, bE4);

t11iC = subs(t11i, [c2, c3, c4], SolC);
t21iC = subs(t21i, [c2, c3, c4], SolC);
% t31iC = subs(t31i, [c2, c3, c4], SolC);
t32ieEC = subs(t32ie, [c2, c3, c4], SolC);
t32ieEC = subs(t32ieEC, [bE1, bE2, bE3], [SolE1.bE1, SolE1.bE2, SolE1.bE3]);
t32eiEC = subs(t32ei, [c2, c3, c4], SolC);
t32eiEC = subs(t32eiEC, [bE1, bE2, bE3], [SolE1.bE1, SolE1.bE2, SolE1.bE3]);
SolI1 = solve(t11iC, t21iC, t32ieEC, t32eiEC, bI1, bI2, bI3, bI4);

t32iiIEC = subs(t32ii, [c2, c3, c4], SolC);
t32iiIEC = subs(t32iiIEC, [bI1, bI2, bI3, bI4], [SolI1.bI1, SolI1.bI2, SolI1.bI3, SolI1.bI4]);
SolI2 = solve(t32iiIEC, bI5);


%% Evaluate objective functions

t44eeeSOL = subs(t44eee, [c2, c3, c4], SolC);
t44eeeSOL = subs(t44eeeSOL, [bE1, bE2, bE3], [SolE1.bE1, SolE1.bE2, SolE1.bE3]);
t44eeeSOL = subs(t44eeeSOL, bE4, SolE2(1));
t44eeeSOL_val = vpa(t44eeeSOL); % Must satisfy: -6/1000 < t44eeeSOL < 0

LstabSOL = subs(Lstab, [c2, c3, c4], SolC);
LstabSOL = subs(LstabSOL, [bI1, bI2, bI3, bI4], [SolI1.bI1, SolI1.bI2, SolI1.bI3, SolI1.bI4]);
LstabSOL = subs(LstabSOL, bI5, SolI2(2));
LstabSOL = subs(LstabSOL, bE4, SolE2(1));
LstabSOL_val = vpa(LstabSOL); % Must satisfy: -1/10 < LstabSOL < 1/10

tau4SOL = subs(tau4, [c2, c3, c4], SolC);
tau4SOL = subs(tau4SOL, [bE1, bE2, bE3], [SolE1.bE1, SolE1.bE2, SolE1.bE3]);
tau4SOL = subs(tau4SOL, [bI1, bI2, bI3, bI4], [SolI1.bI1, SolI1.bI2, SolI1.bI3, SolI1.bI4]);
tau4SOL = subs(tau4SOL, bI5, SolI2(2));
tau4SOL = subs(tau4SOL, bE4, SolE2(1));
tau4SOL_fun = vpa(tau4SOL); % Must satisfy: tau4SOL < 1/10


%% generate constraints functions
CON{1}.field = strcat('LstabSOL');
CON{1}.val = double(vpa(LstabSOL));
CON{1}.lowb = -1/10;
CON{1}.upb = 1/10;
%
CON{2}.field = strcat('t44eeeSOL');
CON{2}.val = double(vpa(t44eeeSOL));
CON{2}.lowb = -6/1000;
% extra boundary to inforce the 0 thereshold
CON{2}.upb = -0.0001;
%% generate obj. functions

FUN.tau4SOL = vpa(tau4SOL);