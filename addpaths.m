% add paths
function SAVE_DIR = addpaths()

    restoredefaultpath; 
    
    home_path = [pwd filesep];
    addpath([home_path '/infer_kernel/']);
    addpath([home_path '/Accel/']);
    addpath([home_path '/Generate_data/']);
    addpath(genpath([home_path '/Examples/']));
    addpath(genpath([home_path '/Plotting/']));
    addpath([home_path '/outputs/']);
    
    SAVE_DIR = [home_path,'outputs\'];
    
    if ~exist('Params','var'), Params = [];     end
    if ~exist(SAVE_DIR,'dir'), mkdir(SAVE_DIR); end
    addpath(SAVE_DIR); 
    
    fprintf('Setting up paths for Learning Dynamics... \n') 
    fprintf(' Data to be saved to folder  \n %s...\n', SAVE_DIR);

