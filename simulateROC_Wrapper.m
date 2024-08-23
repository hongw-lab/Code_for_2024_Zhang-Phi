%%%% set most important parameters

params.pairSet = [1:5]; % add pair idx to run ROC on
params.sessionSet = {'exp'}; % 'sep', 'exp', or any combination (will skip if session does not exist)

%%%% set ROC file information

info.dateTime = datetime; % write only

info.structPath = '/PATH/TO/E/Structure/'; % E structure path
info.structName = 'Struct Name'; % E structure filename

info.rocPath = '/PATH/TO/ROC/OUT/'; % ROC file save path
info.loadROC = false; % true: update existing; false: create new

%%%% set ROC run parameters

% ROC parameters
params.rocTypes = {'self','partner','derived_self','derived_partner'}; % 'self', 'partner', 'derived_self', 'derived_partner', or any combination
params.runAgainst = 'all'; % 'all' or 'other' as baseline
params.nSims = 2000; % number of simulations for null distribution
params.setSeed = true; % set seed (will use previous seed set for animal if updating)

% other parameters
params.nWorkers = 4; % number of workers for parallel pool

params.fullBvList = {'attack','chasing','tussling','threaten','escape','defend','flinch','general-sniffing',...
    'sniff_face','sniff_genital','approach','follow','interaction',...
    'socialgrooming','mount','dig','selfgrooming','climb','exploreobj',...
    'biteobj','stand','attention','human_interfere','other'}; % list of full behavior, must match E struct in sorted order            
params.deriveBvIdx = {[1:3]}; % list of derived behavior idx
                          
if ~exist('E','var')
    disp(['Loading E structure ',info.structName,' from ',info.structPath,'......'])
    load([info.structPath,info.structName], 'dlx_E'); % change variable name to match the file
end

%%%% run ROC

simulateROC(dlx_E,info,params); % change the variable name from dlx_E desired value