function M=BAYESPRFCreateStimFrommrVista(params)
%%Get params file from mrVista Analysis (fFit)

xRange = unique(params.analysis.X);
yRange = unique(params.analysis.Y);

M = [];

[X Y] = meshgrid(xRange, yRange);
nFrames = size(params.stim.images, 2);

% compute the indices I which correspond to the sampling locations
% in params.stim(:).images_org ...
modelCoords = [params.analysis.Y(:) params.analysis.X(:)];
imageCoords = [Y(:) X(:)];
[commonCoords Ia Ib] = intersectCols(imageCoords', modelCoords');

rng = 1:size(params.stim.images_org,2);
presentedImages  = params.stim.images_org(:,rng);

for f = 1:nFrames
    img = zeros(size(X));
    img(Ia) = presentedImages(Ib,f);
    M = cat(3, M, single(img));
end


end