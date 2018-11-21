function U = BAYESPRF_prepare_inputs_polar(ApFrm,TR,Max_X)
% Produces the input structure needed for spm_prf_analyse() given a 3D
% stimuli matrix with dimensions [x,y,t].
%
% Inputs:
%
% ApFrm       - [x,y,t] binary matrix indiciating which pixel locations
%               were illuminated at each time point t
% TR          - Scanner repetition time (TR)
%
% Returns:
%
% U           - Input structure to feed to spm_prf_analyse.m

% Settings
nmicrotime    = 16;     % Bins per TR
%stim_duration = 1;      % Duration of stimuli (secs) == TR in our case
stim_diameter = 2*Max_X;     % Diameter of stimuli in degrees

n = size(ApFrm,3);

for t = 1:n
    
    % Read
    im = ApFrm(:,:,t);
    
    % Skip if a baseline volume
    if ~any(im)
        continue;
    end
    
    % Rescale to 41 x 41 resolution
    res = [41 41];
    im  = imresize(im,res);
    
    % Binarize
    im = im > 0.01;
    
    % Extract stimulated coordinates
    [y,x] = ind2sub(res,find(im));
    
    % Rescale [1,41] to to units of visual angle (degrees)
    r = (stim_diameter/2);
    x = rescale(x, 1, res(1), -r, r);
    y = rescale(y, 1, res(1), -r, r);
    
    % Flip the y axis so positive is up
    y = -y;
            
    % Convert from x,y to distance >= 0 and angle (-pi,pi]
    dist  = sqrt( (x .^ 2) + (y .^ 2) );
    angle = atan2(y,x);

    U(t).dist  = dist;         % Distance
    U(t).angle = angle;        % Angle
    U(t).ons = TR * (t-1);     % Onset (secs)
    %U(t).dur = stim_duration;  % Duration (secs)
    U(t).dur = TR;  % Duration (secs)
    U(t).dt  = TR/nmicrotime;  % 1 bin=1/16 second, i.e. 16 bins per second
    U(t).pmax = stim_diameter; % Stimulus diameter
    U(t).pmin = 0.5;           % Minimum PRF size

    t = t + 1;
end

% -------------------------------------------------------------------------
function new_X = rescale(X, current_min, current_max, new_min, new_max)   
% Rescale X to the range [new_min new_max]
scale_factor = (current_max - current_min) / (new_max - new_min);
new_X = new_min + (X - current_min) / scale_factor;    