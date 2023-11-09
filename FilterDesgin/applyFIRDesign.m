function [ h ] = applyFIRDesign( w, f, fs )
% Convert frequency-domain response to time-domain FIR filter using
% frequency sampling method
%
% INPUT
%   w - filter response in frequency-domain
%   f - linear frequency points
%   fs - sampling rate
%
% OUTPUT
%   h - FIR filter
%
% 2020-08-02 SHEN Yuchen

% conjugate-padding
Spectrum = [ w; conj( w( end-1 : -1 : 2 ) ) ];

% resample coeff
% [up, down] = rat( 0.5 * fs / ( f(end) + ( f(2) - f(1)) / 2 ) );
[up, down] = rat( 0.5 * fs / f(end) );
% convert to time-domain
h = real( ifftshift( ifft( Spectrum' ) ) );
% resample to sampling rate
% h = resample( tmp, up, down);

end
