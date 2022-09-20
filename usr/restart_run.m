% script to restart run
% Note: to restart run, the folder has to be in the ../out/ directory
% 

% clear workspace
clear; close all; clc;

% load parameters file
RunIDin = 'conv100_1lamb';
load(['../out/' RunIDin '/' RunIDin '_par.mat']);


% reassign some variables
RunID   = RunIDin;      % RunID (sometimes I change it after)

restart = 300;          % restart file index
rsstep  = restart*nop;  % current step

NtMax   = 3000;         % max time steps

% run model
run('../src/pantarhei');