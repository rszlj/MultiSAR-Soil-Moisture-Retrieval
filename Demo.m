%Demo
% load the forward LUTs
load("C:\Users\zljnj\Dropbox\Matlab_Codes\EnsembleLearning_soil moisture\Cubes4BareLC\forwardCube_IIEM");
% load the synthetic data
load("C:\Users\zljnj\Dropbox\Matlab_Codes\EnsembleLearning_soil moisture\SyntheticInput\CaseA.mat");
bso = inputCase{1}; % taking the first of the 200 simulations of set A as an example
mv_groundTruth = mvGT{1};% the simulated soil moisture for the first simulation
sig_groundTruth = sig(1);% the simulated rms height for the first simulation
X_groundTruth = X(1);% the simulated correlation length for the first simulation
% set random seed
rng(88,'twister');

% Snapshot ensemble (the SSM ensemble)
ensembleMode=1; % 1: snapshot ensemble
Ne=10;
Nc=3;
Nt=1;
feMode=2;
dryDown=0;
[SSM_mvR,SSM_sig,SSM_X,SSM_cov_ave] = ensembleRetrievalFramework(bso,ensembleMode,forwardCube_IIEM,Ne,Nc,Nt,feMode,dryDown);


% Multi-temporal without dry down (the MT ensemble)
ensembleMode=2; % 1: Multi-temporal without dry down
Ne=10;
Nc=3;
Nt=8; % the full time series used, being 8 in the set A
feMode=2;
dryDown=0;
[MT_mvR,MT_sig,MT_X,MT_cov_ave] = ensembleRetrievalFramework(bso,ensembleMode,forwardCube_IIEM,Ne,Nc,Nt,feMode,dryDown);


% Multi-temporal with dry down (the MTD ensemble)
ensembleMode=3; % 1: Multi-temporal with dry down
Ne=10;
Nc=3;
Nt=8; % the full time series used, being 8 in the set A
feMode=2;
dryDown=1;
[MTD_mvR,MTD_sig,MTD_X,MTD_cov_ave] = ensembleRetrievalFramework(bso,ensembleMode,forwardCube_IIEM,Ne,Nc,Nt,feMode,dryDown);

% display the retrieved soil moisture
len_timeSeries = length(MTD_mvR);
plot(1:len_timeSeries, [SSM_mvR', MT_mvR', MTD_mvR',mv_groundTruth'],'lineWidth',2)
ylabel('Soil moisture [m^3/m^3]','fontsize',11)
xlabel('Time','fontsize',11)
legend('SSM ensemble','MT ensemble','MTD ensemble','GT')



