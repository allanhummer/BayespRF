function BAYESPRF_SINGLESUBJECT(Study,Subject,NumberOfRuns,varargin)
%BAYES_SINGLESUBJECT Summary of this function goes here
%   Detailed explanation goes here

clear GLOBAL
spm fmri 

%% Inputs
p = inputParser;

p.addRequired('Study',@isstr);
p.addRequired('Subject',@isstr);
p.addRequired('NumberOfRuns',@isnumeric);

p.addParameter('StudiesDir','/sacher/melange/fmri/projects/',@isstr);
p.addParameter('Session','trio',@isstr);
p.addParameter('CalcFirstLevelGLM',1,@isnumeric);
p.addParameter('Hemisphere','both',@isstr); %or lh, rh
p.addParameter('TR',1.5,@isnumeric);
p.addParameter('TE',0.0362,@isnumeric);
p.addParameter('AnatomyName','lpivols.nii',@isstr);
p.addParameter('StimFolderName','eightbars_tr15',@isstr);
p.addParameter('StimulusName','eightbars_tr15',@isstr);
p.addParameter('RunSplit','folders',@isstr);
p.addParameter('NiftiName','ravols',@isstr);
p.addParameter('RPName','rp_avols',@isstr);
p.addParameter('FolderNameComment','',@isstr);
p.addParameter('MaxStimSize',9.4,@isnumeric);
p.addParameter('UseParpool',true,@islogical);

p.parse(Study,Subject,NumberOfRuns,varargin{:})

%Show inputs
inputs = p.Results

%% Settings

%Test Dataset Settings:
% NumberOfRuns=10
% TR = 1;
% TE = 0.055;

%Eightbars_tr15
% Repetition time
%TR    = 1.5;
% Echo Time
%TE = 0.0362; 

%Max visual angle degrees - Stimulus Diameter is calculated based on Max_X
Max_X=MaxStimSize;
Max_Y=MaxStimSize;

% Number of sessions
%nsess_firstlevel = 2;
%WhichRuns_PRFAnalysis = 1:2;
WhichRuns_PRFAnalysis=1:inputs.NumberOfRuns;

%Directories
BayesPRF_dir='/z/fmrilab/lab/BayespRF/MyBayesPRF';
BayesPRF_Script_dir=fullfile(BayesPRF_dir,'Scripts');
BayesPRF_Stimuli_dir=fullfile(BayesPRF_dir,'Stimuli');

study_dir=fullfile(inputs.StudiesDir,inputs.Study);
data_dir=fullfile(study_dir,'subjects',inputs.Subject,inputs.Session);
freesurfer_dir=fullfile(data_dir,'segmentation','freesurfer');
surf_dir=fullfile(freesurfer_dir,'surf');

% If it doesn't exist, create Coregistration Matrix for Anatomy

CoregisterMatrixPath=fullfile(surf_dir,'Coregistration.txt');

if ~exist(CoregisterMatrixPath,'file')
    
    disp(['[',mfilename,'] No Coregistration.txt file for registration of anatomy to freesurfer found. Creating one...']);
    
    [~,CoregisterMatrix{1}]=system(['tkregister2 --targ ',fullfile(freesurfer_dir,'mri','orig.mgz'),' --mov ',fullfile(data_dir,'anatomy','nifti','lpivols.nii'),' --s freesurfer --regheader --noedit --reg register.dat | grep Tmov -A 4 | grep -v Tmov']);
    [~,CoregisterMatrix{2}]=system(['tkregister2 --targ ',fullfile(freesurfer_dir,'mri','orig.mgz'),' --mov ',fullfile(data_dir,'anatomy','nifti','lpivols.nii'),' --s freesurfer --regheader --noedit --reg register.dat | grep RegMat -A 4 | grep -v RegMat']);
    
    CoregisterMatrix=cellfun(@(x) strrep(x,';',''),CoregisterMatrix,'UniformOutput',false);
    dlmwrite(CoregisterMatrixPath,CoregisterMatrix,'');
    
else
    
    disp(['[',mfilename,'] Coregistration.txt file for registration of Anatomy to freesurfer found.']);
    
end

% Structural image
Anatomy=fullfile(data_dir,'anatomy','nifti',inputs.AnatomyName);
%Anatomy=fullfile(freesurfer_dir,'mri','brain.nii');

% Directory for creating first level GLM & pRF Analysis
BayesPRF_Subject_dir=fullfile(data_dir,'BayesPRF_',inputs.StimFolderName);
if ~isempty(inputs.FolderNameComment); BayesPRF_Subject_dir=[BayesPRF_Subject_dir,'_',inputs.FolderNameComment]; end
if ~exist(BayesPRF_Subject_dir,'dir'); mkdir(BayesPRF_Subject_dir); end

%% Do It

% Prepare Stim Structure
StimulusStruct=load(fullfile(BayesPRF_Stimuli_dir,[inputs.StimulusName,'.mat']));
Stimulus=StimulusStruct.(cell2mat(fieldnames(StimulusStruct)));
StimStruct = BAYESPRF_prepare_inputs_polar(Stimulus,inputs.TR,Max_X);

if inputs.CalcFirstLevelGLM
BAYESPRF_Run_first_level(BayesPRF_Script_dir,data_dir,surf_dir,Anatomy,inputs.NiftiName,inputs.RPName,inputs.StimFolderName,StimStruct,BayesPRF_Subject_dir,inputs.NumberOfRuns,inputs.RunSplit,inputs.TR,Max_X,Max_Y,inputs.Hemisphere);
end

prf_file = BAYESPRF_Run_pRF_analysis(inputs.Subject,StimStruct,BayesPRF_Subject_dir,inputs.TR,inputs.TE,WhichRuns_PRFAnalysis,inputs.Hemisphere,inputs.UseParpool);

%% Review (single voxel)
%spm_prf_review(prf_file);

%spm_prf_review(prf_file, 50);
%prf_file = fullfile(glm_dir,'PRF_ret15-02_eightbars_tr15_r1.mat');
%spm_prf_review(prf_file, 50);

%% Review (all voxels)

%BAYESPRF_Review
%spm_prf_review(prf_file);


% glm_dir=BayesPRF_Subject_dir;
% % Convert the left V1 label to a nifti mask so we can do some ROI analyses
% label_file = fullfile(data_dir, 'lh_V1.label');
% spm_prf_import_label( label_file, glm_dir );
% 
% % Plot the summed pRF response in right V1
% % Load estimated pRF file
% %prf_file = fullfile(glm_dir,'PRF_ret15-02_eightbars_tr15_r1.mat');
% load(prf_file);
% 
% % Load VOI (imported using spm_prf_import_label)
% roi = fullfile(glm_dir, [hemi '_V1.nii']);
% 
% figure('Color','w');
% spm_prf_summarise(PRF,roi);
% title('Region of interest','FontSize',16);
% 
% % Compute a negative entropy map (certainty of the pRF location)
% % Load estimated pRF file
% %prf_file = fullfile(glm_dir,'PRF_ret15-02_eightbars_tr15_r1.mat');
% load(prf_file);
% 
% % Compute and plot
% spm_prf_plot_entropy(PRF,{'dist','angle'},'dist_angle',true);

end

