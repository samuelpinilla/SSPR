clc;
clear all;
close all;
clear classes;

%% variables
x = randn(1000,1);
n = length(x);

if exist('Params')                == 0,  Params.n           = n;                end
if isfield(Params, 'L')           == 0,  Params.L           = 2;                end         
if isfield(Params, 'T')           == 0,  Params.T           = 100;              end         
if isfield(Params, 'p')           == 0,  Params.p           = 6;                end                           
if isfield(Params, 'y1')          == 0,  Params.y1          = 0.9;              end         
if isfield(Params, 'u0')          == 0,  Params.u0          = 60;               end         
if isfield(Params, 'y')           == 0,  Params.y           = 0.08;             end
if isfield(Params, 'npower_iter') == 0,  Params.npower_iter = 60;               end
if isfield(Params, 'mu')          == 0,  Params.mu          = 1.5/1000;         end
if isfield(Params, 'eta')         == 0,  Params.eta         = 10;               end
if isfield(Params, 'alpha')       == 0,  Params.alpha       = 0.5;              end

        
L           = Params.L; 
m           = floor(n*L)-1;  
Params.m    = m;

display(Params)

Amatrix = randn(m,n);

% Make linear operators; 
A = @(I)  Amatrix*I;
At = @(I) Amatrix'*I;

%% Make signal and data (noiseless)
y = abs(A(x)).^2;

tic
[z0,z,Relerrs] = SSPR(x,y,Params, A, At,Amatrix);
toc
T = Params.T;

%% results
fprintf('Relative error after initialization: %f\n', Relerrs(1))
fprintf('Relative error after %d iterations: %f\n', T, Relerrs(end))

figure;
semilogy(0:length(Relerrs)-1,Relerrs);

xlabel('Iteration'), ylabel('Relative error (log10)'), ...
title('Relative error vs. iteration count')