% add paths

restoredefaultpath; 

home_path = [pwd filesep];
%addpath([home_path '/DE_solver/']);
addpath([home_path '/infer_kernel/']);
addpath([home_path '/Generate_data/']);
%addpath([home_path '/performance/']);
addpath(genpath([home_path '/Examples/']));
addpath(genpath([home_path '/Plotting/']));
addpath(genpath([home_path '/Realdata/']));


%addpath([home_path '/outputs/']);

if ispc
    SAVE_DIR = [getenv('USERPROFILE'), '\outputs\LearningDynamics'];
else
    SAVE_DIR = [home_path,'outputs/'];
    % SAVE_DIR = [getenv('HOME'), 'LearningDynamics/stochasticCase/Outputs/'];
end

if ~exist('Params','var'), Params = [];     end
if ~exist(SAVE_DIR,'dir'), mkdir(SAVE_DIR); end
addpath(SAVE_DIR); 

fprintf('Setting up paths for Learning Dynamics... \n') 
fprintf(' Data to be saved to folder  \n %s...\n', SAVE_DIR);

