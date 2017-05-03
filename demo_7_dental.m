clear
close all

addpath('Toolbox/');

dataset = 'Dental';
do_write_obj = 0;	% Set to 0 to de-active OBJ writing, to 1 to activate it

%%% Set dataset
% File containing the parameters of the LEDs 
light_file = 'Datasets/light.mat';
% File containing the camera's intrinsics
camera_file = 'Datasets/camera.mat';
% Folder containing the photometric stereo images
images_folder = sprintf('Datasets/%s',dataset);

%%% Read calibration files
disp('Reading calibration files');
load(light_file);
load(camera_file);
if(K(1,3)==0) K = K';end % Matlab uses a transposed definition of intrinsics matrix
% Store in a compact structure
calib.S = S; clear S; % Locations (nimgs x 3)
calib.Phi = Phi; clear Phi; % Intensities (nimgs x 3)
calib.Dir = Dir; clear Dir; % Orientations (nimgs x 3)
calib.mu = mu; clear mu; % Anisotropy factors (nimgs x 1)
calib.K = K; clear K; % Intrinsics (3 x 3)

%%% Read dataset
disp('Reading photometric stereo images');
% Get number of images from light calibration
nimgs = size(calib.S,1);
% Read ambient light image
Iamb = double(imread(sprintf('%s/photometric_sample_raw_ambient.png',images_folder)));
[nrows,ncols,nchannels] = size(Iamb);
% Read mask
mask = double(imread(sprintf('%s/photometric_sample_mask_raw.png',images_folder)));
mask = mask(:,:,1)>0;
% Create the cos^4 alpha map to compensate for darkening at the borders
[xx,yy] = meshgrid(1:ncols,1:nrows); % pixel coordinates
xx = xx-calib.K(1,3); % x coordinates w.r.t. principal point
yy = yy-calib.K(2,3); % y coordinates w.r.t. principal point
ff = mean([calib.K(1,1),calib.K(2,2)]); % focal length
cos4a = (ff./sqrt(xx.^2+yy.^2+ff^2)).^4; % cosine-forth law of attenuation
clear xx yy ff
% Read each image, substract ambient light, and compensate for cos^4 alpha
I = zeros(nrows,ncols,nchannels,nimgs);
for i = 1:nimgs
	Ii = double(imread(sprintf('%s/photometric_sample_raw_%04d.png',images_folder,i)));
	Ii = bsxfun(@rdivide,Ii-Iamb,cos4a);
	I(:,:,:,i) = max(0,Ii);
end
clear i Ii cos4a Iamb light_file camera_file images_folder nrows ncols nchannels
% Store in a compact structure
data.I = I; clear I; % Images (nrows x ncols x nchannels x nimgs)
data.mask = mask; clear mask; % Images (nrows x ncols x nchannels x nimgs)

%%% Show the input images
disp('Displaying data');
figure(1001)
subplot(3,1,1)
imshow(uint8(255*data.I(:,:,:,1)./max(data.I(:))))
title('$$I^1$$','Interpreter','Latex','Fontsize',18);
subplot(3,1,2)
imshow(uint8(255*data.I(:,:,:,2)./max(data.I(:))))
title('$$I^2$$','Interpreter','Latex','Fontsize',18);
subplot(3,1,3)
imshow(uint8(255*data.mask))
title('$$\Omega$$','Interpreter','Latex','Fontsize',18);
drawnow

%%% Set parameters
disp('Setting parameters');
params.precond = 'cmg'; % Use multigrid preconditioner
params.z0 = 700*double(data.mask); % Initial depth map: a plane at 700mm from camera
params.estimator = 'Lp'; % Lp norm optimization
params.lambda = 0.5; % p-norm that we used (<1 is nonconvex)
params.indices = 1:nimgs-2; % Use all data
params.ratio = 1; % Downsampling 
params.self_shadows = 1; % Explicitly take into account self-shadows
params.display = 1; % Display result at each iteration

%%% Solve photometric stereo
disp('Solving photometric stereo');
[XYZ,N,rho,Phi,mask] = near_ps(data,calib,params);
% Scale albedo for visualization
rho = rho./max(rho(:));
rho = uint8(255*rho);

%%% Show the result
disp('Displaying results');
figure(1002);
surfl(-XYZ(:,:,1),-XYZ(:,:,2),-XYZ(:,:,3),[0 90]);
shading flat;
colormap gray
view(-220,30);
axis ij
axis image
title('Shape')

figure(1003)
imagesc(rho)
if(size(rho,3)==1)
	colormap gray
end
axis image
axis off
title('Albedo')

figure(1004)
Ndisp = N;
Ndisp(:,:,3) = -Ndisp(:,:,3);
Ndisp = 0.5*(Ndisp+1);
imagesc(Ndisp)
axis image
axis off
title('Normals')
drawnow

%%% Save to an obj file
if(do_write_obj)
	disp('Saving results');
	export_obj(XYZ,N,rho,mask,dataset);
end
