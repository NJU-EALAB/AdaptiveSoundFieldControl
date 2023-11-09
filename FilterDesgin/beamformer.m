function [file_name, p_beam] = beamformer(D, f, beam_ang, beam_wid, h_eq)

if nargin == 4
    h_eq = zeros(1, 2048);
    h_eq(1) = 1;
end



% [ h ] = beamformingR( D, f);
[ h, p_rpc, H_rpc ] = beamformingHJ( D, f);
% load('h_24SPK_WFS(1).mat');
[ MPS, Theta, Phi, Frequency, theta_step, phi_step, p_beam  ] = DAP_beam3D(h, h_eq);
[file_name] = DAP_DAFF(MPS, Theta, Phi, Frequency, theta_step, phi_step, beam_ang, beam_wid);

end

