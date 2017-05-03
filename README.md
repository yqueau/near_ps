# Codes_Near_PS
Code for solving photometric stereo under calibrated or semi-calibrated near point light source illumination (e.g., LEDs).

## Introduction

These Matlab codes implement the method for photometric stereo under point light source illumination advocated in [1,2]. Given a set of photometric stereo images of a still scene acquired from a still, calibrated, pinhole camera, and the light sources parameters, this algorithm estimates depth, normals, albedo and, optionally, lighting intensities. It can output a colored mesh in the .obj format. 

Features:
- Several calibrated datasets
- Graylevel or RGB-valued images
- Various robust estimation techniques (pixelwise image selection, least-squares or M-estimation, ...)
- Optional automatic estimation of lighting intensities (semi-calibrated setup)
- Isotropic or anisotropic (imperfect Lambertian) sources

[1] "Modeling, Calibrating and Solving LED-based Photometric Stereo", Yvain Quéau et al., 2017. 

[2] "Semi-calibrated Near-Light Photometric Stereo", Yvain Quéau et al., Proceedings of the international conference on Scale-Space and Variational Methods for computer vision (SSVM 2017). 

Please cite the above works if using the provided codes and/or datasets for your own research. 

Author: Yvain Quéau, Technical University Munich, yvain.queau@tum.de

## Datasets

- The `Datasets/` folder contains all the datasets used in [1,2]: a plaster statuette, two human faces, a dental plaster cast, a comic book and a box. Each dataset contains 21 files:
  * `photometric_sample_raw_0001.png -> photometric_sample_raw_0008.png`  : RAW images obtained under 8 different nearby LEDs.
  * `photometric_sample_mask_raw.png` : binary mask of the area to reconstruct.
  * `photometric_sample_raw_ambient.png` : image acquired without any LED on, to be substracted from the PS images in order to remove most of the additive bias.
  * `photometric_sample_all_raw.png` : image acquired with all 8 LEDs on. Not used in this work.
  * `photometric_sample_jpg_0001.png -> photometric_sample_jpg_0008.png`, `photometric_sample_jpg_ambient.png`, `photometric_sample_all_jpg.png` : JPG images corresponding to the above RAW ones. They are provided only to ease visualization of the data (because RAW images may be tedious to open with standard image viewers), and they should **NOT** be used for PS.

- Camera's intrinsics (calibrated using Matlab's computer vision toolbox) are provided in the `camera.mat` file, under the forme of a matrix `K = [fx,0,x0 ; 0,fy,y0; 0,0,1]`, with `(fx,fy)` the focal length scaled by the aspect ratio, and `(x0,y0)` the principal point.

- Calibrated LEDs' parameters are provided in the file `light.mat`, which contains the following variables:
  * S (8 x 3): XYZ location of the 8 sources, in mm and w.r.t. the camera center (REQUIRED). 
  * Phi (8 x 3): RGB intensities of the 8 sources, w.r.t. the color of the calibration target we used (assumed to have albedo equal to 1 in each wavelength). If these intensities were not calibrated, they can be estimated automatically using the semi-calibrated option.
  * mu (8 x 1): anisotropy parameter of the 8 sources. mu(i) = 0 means the it-th source is isotropic, mu =1 means it is a primary Lambertian source, and mu > 1 means it is anisotropic (imperfect Lambertian source model).  Set mu(:) = 0 if you have no idea about this parameter.
  * Dir (8 x 3): XYZ orientation of the 8 sources, w.r.t. the optical axis. Each row must have unit-length. Required if mu > 0, not required if mu=0.

- A laser-scan 3D-reconstruction of the statuette is provided for quantitative evaluation in Statuette_GT/. In [1-2], our 3D-reconstructions were first roughly aligned with this ground-truth using the manual point picking tool from the `CloudCompare` software. Then, ICP with farthest point removal was carried out to refine the alignment. The distances in [1,2] refer to the `cloud to cloud` distance tool from this software. 

## Usage

The main fuction is `Toolbox/near_ps.m` (see header for details).

Outputs:
- a gridded point cloud (the third dimension represents depth)
- normal map
- albedo
- lighting intensities
- binary mask
- evolution of energy w.r.t. iterations 

Inputs: 

- a data structure such that:
  * data.I contains the 3D or 4D image stack (**REQUIRED**)
  * data.mask contains a binary 2D mask 

- a calib structure such that:
  * calib.S contains the sources locations (**REQUIRED**)
  * calib.K contains the camera's intrinsics (**REQUIRED**)
  * calib.Phi contains the sources intensities
  * calib.mu contains the sources anisotropy factors
  * calib.Dir contains the sources orientations

- a param structure containing optional parameters. We strongly recommend to play a bit with the following parameters:
  * params.z0 is the initial depth. We advise to roughly (visually) estimate the distance from camera to object in mm, and set z0 to a constant matrix with this rough distance 
  * params.estimator sets the estimator. LS (least-squares) is a good initial choice, but robust M-estimators may be more accurate, though they require the parameter lambda to be set
  * params.lambda sets the M-estimator parameter. For Cauchy estimator, we use 0.1 in our datasets. L1 norm optimization is achieved by setting Lp as estimator, and 1 for lambda
  * params.self_shadows includes self-shadows in  the model or not
  * params.indices can be used to automatically remove brightest or darkest levels in each pixel. This can be useful in the presence of specularities or strong shadowing. Warning: this is not implemented for the semi-calibrated case yet
  * params.semi_calibrated enables automatic intensity refinement

For fast debugging or proof of concept, it may be useful to reduce the size of the data, to limit the number of iterations, or to display live surface, albedo, normals and energy:
  * params.ratio downsamples images by a factor of ratio
  * params.maxit sets the maximum number of iterations
  * params.tol sets the relative stopping criterion on the energy
  * params.display displays the result at each iteration

The inner conjugate gradient iterations can be controlled by:
  * params.maxit_pcg sets the max number of CG iterations within each global iteration
  * params.tol_pcg sets its relative stopping criterion
  * params.precond sets the preconditioner. We strongly recommend to use `cmg` (see Dependencies), but if you want to stick to Matlab's builtin function, use `ichol`


## Demo

The following demo files are provided: 
- `demo_1_calibrated_gray_LS_ps.m` : demo of fast calibrated PS using graylevel-converted images, least-squares estimator, no self-shadows modeling, and image downsampling. If downsampling is removed and CMG preconditioning is used, this recreates the result from Fig. 15 in [1]. 
- `demo_2_calibrated_color_robust_ps.m` : demo of calibrated PS using full-size RGB images, Cauchy's robust M-estimator and explicit self-shadows modeling. This recreates the result from Fig. 16 in [1].
- `demo_3_semicalibrated_color_robust_ps.m` : same, but automatically inferring the lighting intensities. Instead of the optimization over the rank-1 matrix manifold as avised in [2], this script performs simple alternating optimization. Convergence guarantees are thus lost, but in practice the same results are obtained. This script generates result akin to Fig. 3 in [2], though in color. 

To reproduce results from [2], select the appropriate dataset, convert images and light intensities (using the mean over the channels for instance) and use the folowing options:
- Cauchy's estimator with lambda = 0.1
- self_shadows = 1
- indices = 1:nimgs
- semi_calibrated = 1
- CMG preconditioner

To reproduce results from Section 5 in [1], use RGB images with the following options: 
- Cauchy's estimator with lambda = 0.1
- self_shadows = 1
- indices = 2:nimgs-2
- semi_calibrated = 0
- CMG preconditioner

Warning: for now it is not possible to use semi_calibrated PS if the "indices" variable is not equal to "1:nimgs". 
 

## Dependencies

We strongly recommend to use the CMG preconditioner from Koutis et al., which can be downloaded here: 
http://www.cs.cmu.edu/~jkoutis/cmg.html

If CMG it is not installed, set the "precond" parameter to "ichol". Standard incomplete Cholesky will be used, which should prevent any error message in Matlab, but may also be super slow or even non-convergent. 





