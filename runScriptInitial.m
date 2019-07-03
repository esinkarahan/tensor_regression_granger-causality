clear;
close all;
clc
% LOAD THE DATA
% Replace this with your own directory
maindir = pwd;
datadir = fullfile(maindir,'data');
load(fullfile(datadir,'TestCase.mat'));
load(fullfile(datadir,'Lap.mat'));


% model order (number of atoms per mode)
R = 3;

% smoothness/sparsity/orthogonality/nonnegativity on spatial sigantures
% smoothness/sparsity/ on temporal siganture
opts.alphaorth{1}=10;
opts.alphaorth{2}=10;
opts.alphaorth{3}=10;

opts.alphaL1{1} =1;
opts.alphaL1{2} = 1;
opts.alphaL1{3} = 10;

opts.alphaL2{1} = 0.3162;
opts.alphaL2{2} = 0.3162;
opts.alphaL2{3} =10;
opts.maxiters    = 10;
opts.const       = [2 1 1];
opts.nonnegative = [0 1 1];
opts.verbose = 1;
opts.accel = 0;
opts.tol = 1e-8; 
opts.cmoda        = [2 3];
opts.cmodx        = [1 2];
szy = size(Y); Ny = ndims(Y);
opts = calculateEigenValueA(A,szy, Ny,opts);
opts.algorithmOpts{1}.L = sparse(Lap);
opts.algorithmOpts{2}.L = sparse(Lap);
opts.algorithmOpts{3}.L = speye(Nlag);


[Xf,Uf,output] = tensorRegressionCPpenSL(Y,A,R,opts);







