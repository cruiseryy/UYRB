clear;clc;close all

addpath('utils_matlab');
addpath('data_matlab');

% read parameters including
% demand: monthly water demands
% kkcsa: csa ratios for distributing lateral inflow to small reservoir between LYX and LJX
% kkhead: approx waterhead difference used in estimating hydropower production from small reservoirs between LYX and LJX
% smin1, smax1: dead storage capacity and maximum storage capacity for LYX
% smin2, smax2: dead storage capacity and maximum storage capacity for LJX
% see Figure S1: Lat = Inf * b2 + b1 + e_l where e_l ~ N(0, db1 + db2 * Inf)
% inf_mm: climatologies of logarithmic values of Inf
% std_mm: standard deviations of logarithmic values of Inf
% rmin1, rmax1: lower and upper bounds of LYX release decisions
% rmin2, rmax2: lower and upper bounds of LJX release decisions 
% lower bounds = 80% of historical lowest; upper bounds = 120% of historical highest
load('paras.mat')

% read historic 10-day flow data 2001-2009 
% XXX_in: inflow pre-XXX
% XXX_out: (actual) outflow post-XXX
% XXX_s: XXX storage
load('flow_data_new.mat');

