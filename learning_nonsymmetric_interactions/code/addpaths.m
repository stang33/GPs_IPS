% Authored by Jinchao Feng and Sui Tang
% addpaths.m
% Reproducibility-package path setup for the JSC paper
%   "A Sparse Bayesian Learning Algorithm for Estimation of Interaction
%    Kernels in Motsch--Tadmor Model"  (Feng & Tang).
% Run from this 'code/' directory (the working directory must be code/).
% After this, all helper folders (Examples, Generate_data, infer_kernel,
% FSBL_tools, SR_Tool, random_feature_util, Plotting) are on the path.

restoredefaultpath;

home_path = [pwd filesep];

addpath([home_path 'basis_func/']);      % piecewise / B-spline / random-feature bases + phi_piecewise_const
addpath([home_path 'Generate_data/']);
addpath(genpath([home_path 'Examples/']));
addpath(genpath([home_path 'Plotting/']));
addpath(genpath([home_path 'SR_Tool/']));
addpath(genpath([home_path 'FSBL_tools/']));
addpath([home_path 'iSINDy/']);          % implicit-SINDy comparison (Mangan et al. 2016)
addpath(home_path);                      % ensures runners visible

% Output directories live one level UP, alongside main_SBL.tex
SAVE_DIR  = '../';                       % .tex outputs (table_OD.tex, table_CS.tex)
DATA_DIR  = '../trajectory_data/';       % .mat result files (input to make_*)
IMAGE_DIR = '../images/';                % PDF figures (output of make_figure)

if ~exist(SAVE_DIR,'dir'),  mkdir(SAVE_DIR);  end
if ~exist(DATA_DIR,'dir'),  mkdir(DATA_DIR);  end
if ~exist(IMAGE_DIR,'dir'), mkdir(IMAGE_DIR); end

if ~exist('Params','var'), Params = []; end

fprintf('Reproducibility paths configured (cwd = %s).\n', pwd);
fprintf('  data dir : %s\n', DATA_DIR);
fprintf('  image dir: %s\n', IMAGE_DIR);
fprintf('  tex dir  : %s\n', SAVE_DIR);
