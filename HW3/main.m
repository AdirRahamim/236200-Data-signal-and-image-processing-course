clear; clc;

%% Part a
%  compute the DFT representation coeffitions of each row in
% for original and distorted image.

image_original = imread('mandril_original.png');
image_distorted = imread('mandril_distorted.png');

% normalize
image_original = double(image_original)/255;
image_distorted = double(image_distorted)/255;

DFT_original = fftshift(fft2(image_original));
DFT_distorted = fftshift(fft2(image_distorted));

figure(1);
hold on;
imagesc(log(abs(real(DFT_original).^2)));
title('DFT coefficients of original image');
colorbar;
hold off;

figure(2);
hold on;
title('DFT coefficients of distorted image');
imagesc(log(abs(real(DFT_distorted).^2)));
colorbar;
hold off;
%% part b
% Compute c_hat by solving least squares of:
% argmin_c s.t ||DFT-distorted-DFT_original*C||^2
r_orig= rank(DFT_original);
r_distorted = rank(DFT_distorted);

if(r_orig ~= size(DFT_original,1))
    fprintf('The matrix is not full rank!');
end
C_hat = inv(DFT_original)*DFT_distorted;

%% part c
% Distort the original image using the approx c_hat we got, and compute
% the PSNR value between the distored image and the approx distorted image.
approx_distored = DFT_original * C_hat;
approx_distored = ifft2(ifftshift(approx_distored));
imshow(double(approx_distored));
title('Distorted image using aprrox functional map');
psnr_val = compute_psnr(image_distorted, approx_distored);
fprintf('PSNR value = %d\n' ,psnr_val);
% PSNR value we got: 1.552137 * 10^2

%% part d
% Distort the butterfly image using obtained c_hat

 butterfy = imread('Butterfly_.png');
 % normalize:
 butterfy = double(butterfy)/255;
 DFT_butterfly = fftshift(fft2(butterfy));
 distorted_butterfly = DFT_butterfly * C_hat;
 approx_butterfly_distorted = ifft2(ifftshift(distorted_butterfly));
 imshow(double(approx_butterfly_distorted));
 title('Distorted butterfly image using C hat');
 

%% ~~~Q2~~~ 

% Read skycasle and skycastle distorted version.
 skycastle = audioread('skycastle.wav');
 skycastle_distorted = audioread('skycastle-distortion.wav');
 
% Reshape the signal according to the distortion period(512)
 skycastle_re=reshape(skycastle,[],2814);
 skycastle_re=skycastle_re';
 skycastle_dis_re = reshape(skycastle_distorted,[],2814);
 skycastle_dis_re=skycastle_dis_re';
 
 n=size(skycastle_re,2);
 
% Calculate the discrete fourier transform
 DFT_skycastle = skycastle_re*dftmtx(n);
 DFT_skycastle_distorted=skycastle_dis_re*dftmtx(n); 
  
%Same steps for totoro and totoro-distorted signals
 totoro = audioread('totoro.wav');
 totoro_distorted=audioread('totoro-distortion.wav');
 
 totoro_re = reshape(totoro,[],2814);
 totoro_re=totoro_re';
 totoro_dis_re = reshape(totoro_distorted,[],2814);
 totoro_dis_re=totoro_dis_re';
 
 DFT_totoro = totoro_re*dftmtx(n);
 DFT_totoro_distorted = totoro_dis_re*dftmtx(n);

% Compute the denoising matrix for skycaslte-distorted
 c = pinv(DFT_skycastle_distorted)*DFT_skycastle;
 skycastle_denoised = DFT_skycastle_distorted * c;
 fprintf('MSE skycatle = %d\n', norm(DFT_skycastle-skycastle_denoised));
% MSE error value = 6.342565*10^-2

% Compute denoised totoro by multiply with c
 totoro_temp = DFT_totoro_distorted*c;
 
% Compute inverse fourier transform of denoised signal 
 totoro_denoised = totoro_temp * inv(dftmtx(n));
 totoro_denoised = totoro_denoised';
 
% Reshape again to vector 
 totoro_denoised = reshape(totoro_denoised,[],1);

x = 1:size(totoro_denoised,1);
fprintf('MSE error totoro = %d' , norm(totoro-totoro_denoised));
figure(1);
plot(x , totoro_denoised);
title('Denoised totoro');
figure(2);
plot(x,totoro);
title('Original totoro');
% MSE value: 5.831843*10^-2