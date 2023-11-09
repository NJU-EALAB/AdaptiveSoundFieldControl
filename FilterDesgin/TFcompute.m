function [G0] = TFcompute(Data, f_BF)
    N_omega_BF = length(f_BF);
    fs = Data{1,1}.fs;
    nfft = 2^16;
    f_TD = 0 : fs / nfft : (nfft - 1) * fs / nfft;
    f_TD = f_TD( 1: end/2);
    N_theta = length(Data);
    theta_interval = 360 / N_theta;
    theta = 0 : theta_interval : 360 - theta_interval;
    N_spk = size(Data{1,1}.IR,2);
    G0 = zeros(N_theta, N_spk, length(f_BF));
    f_num = zeros(1, length(f_BF));
    for i1 = 1 : length(f_BF)
        [~, N_tmp] = min(abs(f_TD - f_BF(i1)));
        f_num(i1) = N_tmp;
    end
    
    for i1 = 1 : N_theta
        for i2 = 1 : N_spk
            IR_tmp = Data{i1, 1}.IR{1, i2};
            FR_tmp = fft(IR_tmp, nfft);
            G0(i1, i2, :) = FR_tmp(f_num);
        end
    end
    
%     spk_plot = [1 : 3 : 10];
%     for i1 = 1 : length(spk_plot)
%         G0_tmp = squeeze(G0(:, spk_plot(i1), :));
%         G0_tmp = G0_tmp ./ max(abs(G0_tmp));
%         SPL_plot = mag2db(abs(G0_tmp));
%         meshTheta = repmat(theta',1,N_omega_BF);
%         meshF = repmat(f_BF, N_theta,1);
%         figure
%         surf(meshF,meshTheta,SPL_plot);
%         view(0,90);
%         set(gca,'xscale','log');
%         set(gca, 'FontSize', 14, 'FontName', 'Times New Roma');
%         xlim([200,2e4]);
%         xlabel('Frequency / Hz'); 
%         ylabel('Angle / degree');
%         caxis([-20,0]);
%         colorbar;
%         colormap('jet');
%         shading interp;
%     end
%     
end

