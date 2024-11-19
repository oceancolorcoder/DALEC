
function [abs_QWIP_score,QCI,AVW] = AVW_QWIP_2D_fun(Rrs,wavelengths,AVW_coefs,sensor)
% function [abs_QWIP_score,QCI,AVW] = AVW_QWIP_2D_fun(Rrs,wavelengths,AVW_coefs,sensor)
%
%%% 1) Compute AVW for a 2D array of Rrs spectra (NOT AN IMAGE)
%%% 2) Determine the QWIP score for all the data points in an image
% based on original by Ryan Vandermuelen
%
%% Adapted by D.A.Aurin (Not for imagery)
%   2022-03-29: DAA
%       From 3D image-based approach (lat x lon x lambda) of RV to 2D matrix
%       Rrs: mxn array [sr^1] where m is the number of spectra to be
%           evaluated as a group and n is the length of wavelength.
%       wavelengths: 1xn [nm] (i.e. all spectra conform)
%
%   2022-09-23: DAA
%       Further altered to adjust multiband AVW using coefficients
%           derived by Ryan for each sensor.
%       AVW_coefs: structure with fields named for sensors, or 'NONE'
%       AVW_coef.(sensor).coefs: 6 polynomial coefficients
%       AVW_coef.(sensor).bands: actual (not nominal) bands of multispectral
%           sensor used in polynomial fit
%       Set sensor to 'NONE' to forego multispectral correction
%       Anything with < 20 wavebands will be considered multispectral and
%       should have AVW coefs in the input structure.

m = size(Rrs,1); % number of spectra

% Define the visible range of wavelenghts.
% Note: AVW doesn't use the UV/NIR in the calculation.
min_wave = 400;
max_wave = 700;
[~,min_index] = min(abs(wavelengths-min_wave));
[~,max_index] = min(abs(wavelengths-max_wave));
multi = 20; % Max number of wavebands defining a multispectral Rrs

RRS = Rrs(:,min_index:max_index);
wave_vis = wavelengths(:,min_index:max_index);
wave_1nm = (400:700);

for i = m:-1:1 % <- reverse loop for speed and to avoid preallocation        
    if length(wave_vis) > multi
        % For hyperspectral data only, spline to 1nm
        interim_arr = RRS(i,:);
        interim_arr(isnan(interim_arr))=0;
        RRS_1nm = interp1(wave_vis,interim_arr,wave_1nm,'spline');
        AVW_uncal(i) = sum(RRS_1nm)./sum(RRS_1nm./wave_1nm);
    else
        % Multispectral; no interpolation
        AVW_uncal(i) = sum(RRS(i,:))./sum(RRS(i,:)./wave_vis);
    end

end

% If Rrs is multispectral with coefficients defined apply polynomial fit correction
if length(wave_vis) > multi || ~isfield(AVW_coefs,sensor)
    % fprintf('Hyperspectral data or sensor %s not found in library. No calibration.\n',sensor)
    AVW = AVW_uncal;
else
    c0 = AVW_coefs.(sensor).coefs(1);
    c1 = AVW_coefs.(sensor).coefs(2);
    c2 = AVW_coefs.(sensor).coefs(3);
    c3 = AVW_coefs.(sensor).coefs(4);
    c4 = AVW_coefs.(sensor).coefs(5);
    c5 = AVW_coefs.(sensor).coefs(6);
    AVW = c0*AVW_uncal.^5 + c1*AVW_uncal.^4 + c2*AVW_uncal.^3 + c3*AVW_uncal.^2 + c4*AVW_uncal + c5;

%     % Test calibration
%     figure
%     plot(AVW_uncal,AVW,'.')
%     p11
%     xlabel('Uncalibrated AVW [nm]')
%     ylabel('Calibrated AVW [nm]')
%     title(sprintf('Instrument: %s', sensor))
%     grid on
%     print(sprintf('plt/calibration_%s.png',sensor),'-dpng')
end


% Calculate the Quality Control Index (QCI)
%   (Rrs_665-Rrs_490)/((Rrs_665+Rrs_490)
%   Find the closest wavelength match and compute the polynomial fit.
%   This is not ideal for multispectral data, but should be corrected
%   for with the polynomials above.
[~,index_490] = min(abs(wave_vis-490));
[~,index_665] = min(abs(wave_vis-665));
QCI = (RRS(:,index_665) - RRS(:,index_490))./(RRS(:,index_665) + RRS(:,index_490));
p = [-8.399885e-09,1.715532e-05,-1.301670e-02,4.357838,-5.449532e02];

% Generate an array of the "predicted" QCI based on the AVW
%   values and subtract this from the actual QCI calculated above to yield 
%   the QWIP score.
QCI_pred = (p(1)*AVW.^4 + p(2)*AVW.^3 + p(3)*AVW.^2 + p(4)*AVW.^1 + p(5));
QCI = QCI';
QWIP_score = QCI_pred-QCI;
abs_QWIP_score = abs(QWIP_score);





