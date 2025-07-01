close all, clear all,clc
%%
global vars
temp = dir('*.dat');
vars.datfile = temp(2).name;
vars.reps = 180;


vPIVOT__Read_vPIVOT_datXA30A(vars.datfile,vars.reps);
% vPIVOT__Read_vPIVOT_datVE11C(vars.datfile,vars.reps);

%%
vPIVOT__Inf_SvO2_reorganization
vPIVOT__Inf_T2star_analysis
vPIVOT__Inf_SvO2_analysis
%%
% % %%
vPIVOT__Sup_PC_reorganization
vPIVOT__Sup_SvO2_analysis
vPIVOT__Sup_PC_analysis

%%
vPIVOT__ASL_ksp2image
vPIVOT__ASL_perfusion_quantification
%%
init.datfile = vars.datfile;
init.reps = vars.reps;
varsASL = vars.ASL;
varsPC = vars.PC;
varsInfGRE = vars.InfGRE;


%%
Perf_TC = varsASL.Perf.AverageTimecourse_IMG_mean;
NormT2star_TC = varsInfGRE.T2star.AverageTimecourse.NormTEs2_5;
T2star_TC = varsInfGRE.T2star.AverageTimecourse.TEs2_5;
SvO2_TC = varsInfGRE.SvO2.AverageTimecourse.SvO2;

%%

save(['20250630_Perf_Analysis_',vars.datfile(end-15:end-4),'.mat'],'varsASL','init','varsPC','varsInfGRE')
