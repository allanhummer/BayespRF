function prf_file=BAYESPRF_Run_pRF_analysis(Subject,U,glm_dir,TR,TE,sess,hemi,UseParpool)
% Example script to run a simple pRF analysis on a voxel-by-voxel basis 
% and view the results.
%
% Neuronal model:      Single Gaussian function
% Receptive field:     Circular (isotropic)
% Input specification: Polar coordinates
%
% Data is from SamSrf

% Settings

% Directory of the downloaded example dataset
%data_root_dir = '/net/mri.meduniwien.ac.at/departments/physics/fmrilab/lab/BayespRF/ret15-02/data/';
%data_dir      = fullfile(data_root_dir,'Example','pRF');
%data_dir      = '/net/mri.meduniwien.ac.at/departments/physics/fmrilab/lab/BayespRF/ret15-02/data/';


% Directory of GLM
%glm_dir  = '/net/mri.meduniwien.ac.at/departments/physics/fmrilab/lab/BayespRF/ret15-02/GLM';

%TR = 1.5;         % Repetition time
%TE = 0.0362;     % Echo time

% Which sessions to include
%sess = 1;
num_sess = length(sess);

%% Prepare inputs

% Build a structure containing which stimulus pixels were illuminated at
% each time step.
%load(fullfile(BayesPRF_Script_dir,[StimulusName,'.mat']));
%U = BAYESPRF_prepare_inputs_polar(eightbars,TR,Max_X);

% The timeseries from each session are stored in a VOI_xx.mat file. Build a
% cell array of the VOI files for each session.
% if strcmp(hemi,'both')
%    
%     xY = cell(1,2*num_sess);
%     for i = 1:2*num_sess
%         if i<=num_sess
%             hemi='lh';
%             session=sess(i);
%         else
%             hemi='rh';
%             session=sess(i-num_sess);
%         end
%         filename{i} = sprintf('VOI_%s_prf_mask_%d.mat',hemi,session);
%         xY{i}    = fullfile(glm_dir,filename{i});
%     end
%
% else

xY = cell(1,num_sess);
for i = 1:num_sess
    filename = sprintf('VOI_%s_prf_mask_%d.mat',hemi,sess(i));
    xY{i}    = fullfile(glm_dir,filename);
end
%end



%% Specify pRF model

% Load SPM for timing information / image dimensions
SPM = load(fullfile(glm_dir,'SPM.mat'));
SPM = SPM.SPM;

% Update SPM path as we don't know where this example will be saved
SPM.swd = glm_dir;

%PRF Model Name
%prf_filename = [Subject,'_',hemi,'_eightbars_tr',strrep(num2str(TR),'.',''),'_r',strrep(num2str(sess),' ','')];
prf_filename = [Subject,'_eightbars_tr',strrep(num2str(TR),'.',''),'_r',strrep(num2str(sess),' ','')];

% Set pRF specification options
options = struct('TE', TE,...
                 'voxel_wise', true,...
                 'name', prf_filename,...
                 'model', 'spm_prf_fcn_gaussian_polar',...
                 'B0',3);
             
% Specify pRF model (.mat file will be stored in the GLM directory)
PRF = spm_prf_analyse('specify',SPM,xY,U,options);

%% Estimate one voxel as an example (voxel 1)
%voxel = 1;

% Model to estimate
%prf_file = fullfile(glm_dir,['PRF_',hemi,'_',prf_filename,'.mat']);

% Estimation options
%options  = struct('voxels',voxel);
%options  = struct('voxels',3439);
%options  = struct('voxels',1:6669);

% Estimate
%PRF_est = spm_prf_analyse('estimate',prf_file,options);

% Review
%spm_prf_review(prf_file, voxel);

%% Estimate all voxels (slow)

% % Estimation options
options = ...
    struct('use_parfor', UseParpool, ... % Parallelization
    'nograph', true, ...    % If true, disables plots
    'usegpu', false);

% options = ...
%     struct('use_parfor', true, ... % Parallelization
%     'nograph', true, ...    % If true, disables plots
%     'voxels', 1, ... 
%     'usegpu', false);

% options = ...
%     struct('init', 'GLM_P', ...     % Initialization. See below
%     'use_parfor', true, ... % Parallelization
%     'nograph', true, ...    % If true, disables plots
%     'voxels', 3439, ...       % Voxels indices (optional)
%     'usegpu', true);

GPUScriptPath='/net/mri.meduniwien.ac.at/departments/physics/fmrilab/lab/BayespRF/MyBayesPRF/MySPMINT';

if options.usegpu
    addpath(GPUScriptPath)
else
    if contains(path,GPUScriptPath)
        rmpath(GPUScriptPath)
    end
end

% Estimate
% Model to estimate

prf_file = fullfile(glm_dir,['PRF_',prf_filename,'.mat']);

PRF_est = spm_prf_analyse('estimate',prf_file,options);

%% Review (single voxel)
%prf_file = fullfile(glm_dir,'PRF_ret15-02_eightbars_tr15_r1.mat');

%spm_prf_review(prf_file, 50);