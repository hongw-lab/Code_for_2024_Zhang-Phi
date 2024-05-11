% dependencies -- deriveBehavs

addpath(genpath(fullfile('Z:/NR/NR_Toolbox/'))) % add simulateROC and dependencies to path

clearvars -except E

%%%% set most important parameters

params.pairSet = [1:5]; % add pair idx to run ROC on
params.sessionSet = {'exp'}; % 'sep', 'exp', 'exp#' (e.g. 'exp2'), or any combination (will skip if session does not exist)

%%%% set ROC file information

info.dateTime = datetime; % write only

info.structPath = '/PATH/TO/E/Structure/'; % E structure path
info.structName = 'Struct Name'; % E structure filename

info.rocPath = '/PATH/TO/ROC/OUT/'; % ROC file save path
info.loadROC = true;

%%%% set ROC run parameters

% ROC parameters
params.rocTypes = {'self','partner','derived_self','derived_partner'}; % 'self', 'partner', 'derived_self', 'derived_partner', or any combination
params.runAgainst = 'all'; % 'all' or 'other' as baseline
params.nSims = 2000; % number of simulations for null distribution
params.setSeed = true; % set seed (will use previous seed set for animal if updating)

% other parameters
params.nWorkers = 4; % number of workers for parallel pool

params.fullBvList = {'attack','chasing','sniff','dig'}; % list of behaviors to run ROC against              
params.deriveBvIdx = {[1:3]}; % list of derived behavior idx
                          
if ~exist('E','var')
    disp(['Loading E structure ',info.structName,' from ',info.structPath,'......'])
    load([info.structPath,info.structName], 'E');
end

%%%% run ROC

simulateROC(E,info,params);