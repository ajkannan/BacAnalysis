% this code implements minimum probability flow learning for a fully connected Ising model.  see http://arxiv.org/abs/0906.4779

% Author: Jascha Sohl-Dickstein (2009)
% Web: http://redwood.berkeley.edu/wiki/Jascha_Sohl-Dickstein
% This software is made available under the Creative Commons
% Attribution-Noncommercial License.
% (http://creativecommons.org/licenses/by-nc/3.0/)
clear;clc;commandwindow;close all;
addpath(genpath('./minFunc_2012'));

% initialize
Xall = csvread('~/BacAnalysis/Data/bacteria_data_pared.csv');
Xall = Xall';
d = size(Xall,1); % number of units
maxlinesearch = 10000; % this number is excessive just to be safe!!!!!! learning works fine if this is just a few hundred
independent_steps = 10*d; % the number of gibbs sampling steps to take between samples

minf_options = [];
%options.display = 'none';
minf_options.maxFunEvals = maxlinesearch;
minf_options.maxIter = maxlinesearch;
minf_options.useMex = 0;

% make the weight matrix repeatable
rand('twister',355672);
randn('state',355672);

% randomly initialize the parameter matrix we're going to try to learn
% note that the bias units lie on the diagonal of J
Jest = randn( d, d ) / sqrt(d) / 100;
Jest = Jest + Jest';
Jest = Jest/2;

% perform parameter estimation
fprintf( '\nRunning minFunc for up to %d learning steps...\n', maxlinesearch );
t_min = tic();

%%%%%%%%%%% choose one of these two Ising model objective functions %%%%%%%%%%%%%
% K_dK_ising is slightly faster, and includes connectivity only to states
% which differ by a single bit flip.
%Jnew = minFunc( @K_dK_ising, Jnew(:), minf_options, Xall );
% K_dK_ising_allbitflipextension corresponds to a slightly modified choice
% of connectivity for MPF. This modified connectivity includes an
% additional connection between each state and the state which is reached
% by flipping all bits.  This connectivity pattern performs better in cases
% (like neural spike trains) where activity is extremely sparse.
Jest = minFunc( @K_dK_ising_allbitflipextension, Jest(:), minf_options, Xall );

Jest = reshape(Jest, [d,d]);
t_min = toc(t_min);
fprintf( 'parameter estimation in %f seconds \n', t_min );

% Calculate partition function to normalize probabilities
%Jest_tri = triu(Jest);Jest_tri(1:d+1:end) = 0;
bits = 0:(2^5 - 1);
bits = dec2bin(bits);
bits = bits - '0';
Z = 0;
for x = bits'
    Z = Z + exp(-1 * x' * Jest * x);
end

% Calculate log likelihood of seeing data given the model parameters found
nsamples = size(Xall);
nsamples = nsamples(2);
log_likli = 0;
sample_probabilities = zeros(size(Xall,2),1);
for i = 1:(nsamples)
    curr_sample = Xall(:,i);
    sample_probabilities(i) = exp(-1 * curr_sample' * Jest * curr_sample - log(Z));
    log_likli = log_likli + -1 * curr_sample' * Jest * curr_sample - log(Z);
end
sample_probabilities
csvwrite('sample_probabilities.csv', sample_probabilities);
log_likli
