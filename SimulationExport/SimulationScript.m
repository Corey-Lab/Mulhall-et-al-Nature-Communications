freq=[[1:2:10],[20:20:100],[200:200:1000],10000]; % Different frequencies to simulation
amp=[0,1,2,4,8,16,30,60,120]; % Different amplitude this is the amplitude at the tiplink

%% Define Bond parameters
drag=0;
lifetime=[];
Rebind=[];
dF=[];

CN=0.4799/1000; % Effective Concentration in M
fcN=10.4;  % Decay force for force dependent Concentration in pN
model=3; % Uses an exp(-(F/fcN)^2) model for Effective concentration, 0 is constant concentration
fbN=14.3; % Force sensitivity of a single bond in pN
koffN=0.58; % Off rate of single bond 1/s
konN=70000; % On Rate of a single bond 1/M/s

Nreps=1000; % Number of Replicates for each simulation


adaptation=0; % set to 1 to enable motor adaptation and 0 to disable.

%%
[lifetime_avg,Forceout,RB] = TipLinkExperimentFunction(CN,fcN,model,fbN,koffN,konN,Nreps,adaptation,freq,amp,drag);
%%
%save('NoAdaptation_Rebinding.mat');