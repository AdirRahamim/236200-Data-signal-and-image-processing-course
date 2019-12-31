function [ psnr_val ] = compute_psnr(y_original, y_estimated)
% COMPUTE_PSNR Computes the PSNR between two images
%
% Input:
%  y_original  - The original image
%  y_estimated - The estimated image
%
% Output:
%  psnr_val - The Peak Signal to Noise Ratio (PSNR) score

y_original = y_original(:);
y_estimated = y_estimated(:);

dynamic_range = max(y_original)-min(y_original);

mse_val = mean((y_estimated-y_original).^2);

psnr_val = 10*log10(dynamic_range.^2/(mse_val+eps));

end

