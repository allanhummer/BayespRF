function BAYESPRF_Run_first_level(BayesPRF_Script_dir,data_dir,surf_dir,Anatomy,NiftiName,RPName,StimFolderName,U,glm_dir,nsess,RunSplit,TR,Max_X,Max_Y,hemi)
%% Runs GLM analysis and extracts timeseries

% Settings

% Directory into which to download example dataset
%data_root_dir = '/z/fmrilab/lab/BayespRF/ret15-02/data';
%data_dir      = fullfile(data_root_dir,'Example','pRF');
%surf_dir      = fullfile(data_root_dir,'example','surf');
%data_dir = '/z/fmrilab/lab/BayespRF/ret15-02/data';

%NiftiName='ravols.nii';

% Directory for creating GLM
%glm_dir  = fullfile('/z/fmrilab/lab/BayespRF/ret15-02/GLM');

% Number of sessions
%nsess = 3;

% Repetition time
%TR    = 1.5;
%% Prepare onsets
%StimulusStruct=load(fullfile(BayesPRF_Script_dir,'Stimuli',[StimulusName,'.mat']));
%Stimulus=StimulusStruct.(cell2mat(fieldnames(StimulusStruct)));
%U = BAYESPRF_prepare_inputs_polar(Stimulus,TR,Max_X);

bins_x = [-Max_X 0 Max_X];
bins_y = [Max_Y 0 -Max_Y];

% Build time x screen bins matrix (3x3 screen bins)
onsets_matrix = zeros(length(U), length(bins_x) .^ 2);
for t = 1:length(U)
    % Loop over pixels activated at this time point
    for activated_pixel = 1:length(U(t).dist)
        % Get location
        dist  = U(t).dist(activated_pixel);
        angle = U(t).angle(activated_pixel);
        
        % Polar->cartesian
        x = dist * cos(angle);
        y = dist * sin(angle);

        % Identify closest bin
        [~,binx] = min(abs(bins_x-x));
        [~,biny] = min(abs(bins_y-y));

        % Binned coordintes -> index
        bin_idx = sub2ind([length(bins_x) length(bins_x)],biny,binx);

        onsets_matrix(t,bin_idx) = onsets_matrix(t,bin_idx) + 1;
    end
end

% Remove empty bins
onsets_matrix = onsets_matrix(:,any(onsets_matrix));
num_regressors = size(onsets_matrix,2);

% SPM inputs
names = cell(1,num_regressors); 
onsets = cell(1,num_regressors); 
durations = cell(1,num_regressors); 

for t = 1:num_regressors
    names{t} = ['Bin' num2str(t)];
    onsets{t} = (find( onsets_matrix(:,t) ) - 1) * TR;
    durations{t} = 0;
end

cd(glm_dir);
save('onsets.mat', 'names', 'onsets', 'durations');

%% Specify first level design

%start_dir = pwd;

% Make output directory
%if ~exist(glm_dir,'file')
%    mkdir(glm_dir);
%end

load(fullfile(BayesPRF_Script_dir,'first_level_batch.mat'));

% Session-specific options



for i = 1:nsess
    
    if strcmp(RunSplit,'folders')
        NiftiDir = fullfile(data_dir,[StimFolderName,'_r',num2str(i)],'nifti');
        movement = spm_select('FPList',NiftiDir,RPName,'.txt');
        epis     = spm_select('ExtFPList',NiftiDir,['^',NiftiName,'.nii'], 1:999);
    elseif strcmp(RunSplit,'files')
        NiftiDir = fullfile(data_dir,StimFolderName,'nifti');
        movement = spm_select('FPList',NiftiDir,[RPName,num2str(i),'.txt']);
        epis     = spm_select('ExtFPList',NiftiDir,['^',NiftiName,num2str(i),'.nii'], 1:999);
    end
    
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi_reg = cellstr(movement);
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).scans     = cellstr(epis);
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi     = cellstr('onsets.mat');
    
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).hpf = 128;    
    
    clear NiftiDir movement epis    
end

% Model spec options
matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(glm_dir);
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;

% Run
spm_jobman('run',matlabbatch);

%% Import cortical surface
%
% Creates images: GLM/lh_surface.nii and GLM/rh_surface.nii
% and .mat files: GLM/lh_Srf.mat and GLM/rh_Srf.mat

% Left hemisphere
spm_prf_import_surface(glm_dir, Anatomy, surf_dir, 'lh');

% Right hemisphere
spm_prf_import_surface(glm_dir, Anatomy, surf_dir, 'rh');

%% Build a mask of voxels which survive p < 0.001
clear matlabbatch;
matlabbatch{1}.spm.stats.results.spmmat = cellstr(fullfile(glm_dir,'SPM.mat'));
matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
matlabbatch{1}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.001;
matlabbatch{1}.spm.stats.results.conspec.extent = 0;
%matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
%matlabbatch{1}.spm.stats.results.export{1}.binary.basename = 'mask_uncorrected';
matlabbatch{1}.spm.stats.results.write.tspm.basename = 'mask_uncorrected';
%spm_jobman('interactive',matlabbatch);
spm_jobman('run',matlabbatch);
%% Remove voxels anterior to y = 0
cd(glm_dir);

% Read
V = spm_vol('spmF_0001_mask_uncorrected.nii');
[Y,XYZmm] = spm_read_vols(V);

% Threshold
i = XYZmm(2,:) > 0;

% Write
Y(i) = 0;
spm_write_vol(V,Y);

%cd(start_dir);

%% Extract timeseries from surface voxels which survive p < 0.001

% Identify masks
spm_F_mask   = fullfile(glm_dir,'spmF_0001_mask_uncorrected.nii');
if strcmp(hemi,'both')
    surface_mask{1} = fullfile(glm_dir,'lh_surface.nii');
    surface_mask{2} = fullfile(glm_dir,'rh_surface.nii');
else
    surface_mask{1} = fullfile(glm_dir,[hemi '_surface.nii']);
end

% Prepare batch

clear matlabbatch;
MatlabBatchTemplate=load(fullfile(BayesPRF_Script_dir,'extract_timeseries_batch.mat'));

matlabbatch{1}=MatlabBatchTemplate.matlabbatch{1};
matlabbatch{1}.spm.util.voi.name   = [hemi,'_prf_mask'];
matlabbatch{1}.spm.util.voi.spmmat = cellstr(fullfile(glm_dir,'SPM.mat'));
matlabbatch{1}.spm.util.voi.roi{1}.mask.image = cellstr(spm_F_mask);
matlabbatch{1}.spm.util.voi.roi{2}.mask.image = cellstr(surface_mask{1});
if strcmp(hemi,'both')
    matlabbatch{1}.spm.util.voi.roi{3}.mask.image = cellstr(surface_mask{2});
    matlabbatch{1}.spm.util.voi.expression = 'i1 & (i2 | i3)';
else
    matlabbatch{1}.spm.util.voi.expression = 'i1 & i2';
end

% Run batch
spm_jobman('run',matlabbatch);


end