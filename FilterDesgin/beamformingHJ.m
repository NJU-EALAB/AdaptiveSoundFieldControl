function [ h, p_rpc, H_rpc ] = beamformingHJ( D, f)
% output
% h: FIR filter
% N_spk * N_order
%
% input
% D: beamforming target
% f: frequencies corresonding to the preset target

% posi: position input of loudspeaker array, 1 means the 1st speaker is on
% the ceiling side, 0 means the 1st speaker is on the floor side.
% 
addpath('C:\Users\Administrator\Documents\RoomSimulation\SourceBeamCtrl\beamformer\Data');
addpath(genpath('C:\Users\Administrator\Documents\RoomSimulation\SourceBeamCtrl\beamformer\SpaRSA_2.0'));

posi = 1;

N_omega_BF = 65;
fs = 48e3;
fh = fs / 2;
fl = 0;
f_BF = linspace(fl,fh,N_omega_BF);

N_theta_D = size(D, 1);
theta_interval_D = 360 / N_theta_D;
theta_D = 0 : theta_interval_D : 360 - theta_interval_D;

%% CDPS data loading
load('Data.mat');
N_theta = length(Data);
theta_interval = 360 / N_theta;
theta = 0 : theta_interval : 360 - theta_interval;
N_spk = size(Data{1,1}.IR,2);
fs_data = Data{1,1}.fs;
if fs ~= fs_data
    error('sampling frequency mismatch');
end
[G0] = TFcompute(Data, f_BF);

%% angle calculation
direc_ang = zeros(1, N_spk);
ang_wid = 20;
num_wid = round(ang_wid / 2 / theta_interval);
for i1 = 1 : N_spk
    power_distr = zeros(1, N_theta);
    G0_tmp = squeeze(G0(:, i1, :));
    for i2 = 1 : N_theta
        theta_sel = i2 - num_wid : i2 + num_wid;
        for i3 = 1 : length(theta_sel)
            if theta_sel(i3) <= 0
                theta_sel(i3) = N_theta + theta_sel(i3);
            elseif theta_sel(i3)> N_theta
                theta_sel(i3) = theta_sel(i3) - N_theta;
            end
        end
        power_tmp = sum(sum((abs(squeeze(G0_tmp(theta_sel, :)))).^2));
        power_distr(i2) = power_tmp;
    end
    [~, N_tmp] = max(power_distr);
    direc_ang(i1) = N_tmp;
end

main_direc = mean(direc_ang); % revision could be made, what if N_theta~=360
if direc_ang(1) > main_direc
    posi_array = 1;
elseif direc_ang(1) < main_direc
    posi_array = 0;
end

%% target set
tar = zeros(N_theta, N_omega_BF); 
for i1 = 1 : N_omega_BF
    for i2 = 1 : N_theta
        [~, N_tmp1] = min(abs(f - f_BF(i1)));
        [~, N_tmp2] = min(abs(theta_D - theta(i2)));
        tar(i2, i1) = squeeze(D(N_tmp2, N_tmp1));
    end
end

%% target adjustment
if posi_array ~= posi
    G0 = flip(G0, 1);
    main_direc = N_theta - main_direc + 1;
end

angle_delta = main_direc - (N_theta / 2 + 1);
angle_delta = floor(angle_delta);

tar_tmp = tar;
for i1 = 1 : N_theta
    num_new = i1 + angle_delta;
    if num_new > N_theta
        num_new = num_new - N_theta;
    elseif num_new <= 0
        num_new = N_theta + num_new;
    end
    tar(num_new, :) = tar_tmp(i1, :);
end

tar = [tar(end-2:end,:); tar(1:end-3,:)];
% tar = zeros(N_theta, N_omega_BF); 
% tar(175:185, :) = 1;

%% beamforming in frequency domain
H = ones(N_spk, N_omega_BF);
for i1 = 2 : N_omega_BF
    tau = 0;
    G_tmp = squeeze(G0(:, :, i1));
    tar_tmp = squeeze(tar(:, i1));
    % Gradient Descent
    [w,~,~,...
          ~,~,~]= ...
        SpaRSA(tar_tmp,G_tmp,tau,...
        'Monotone',1,...
        'Debias',0,...
        'Initialization',0,...
        'StopCriterion',1,...
        'ToleranceA',0.001, ...
        'MaxiterA',2000);
    % G0 represents acoustic field, d0 represents target field
    H(:,i1) = w;
end


[~,max_idx] = max(rms(abs(H),2));
H = H./(H(max_idx,:));

fn = logspace(log10(8e3),log10(fs/2),N_spk/2);
Wn = fn/fs*2;
Wn(end) = 1; max_idx_win = N_spk/2;
Wn = [Wn, flip(Wn)];

% syn = max_idx_win - spk_idx;
syn = 0;
if syn<0
    Wn = [Wn(end-abs(syn)+1:end) Wn(1:end-abs(syn))];
elseif syn>0
    Wn = [Wn(abs(syn)+1:end) Wn(1:abs(syn))];
end
idx = find(Wn == 1);
Wn = [sort(Wn(1:idx(1))), flip(sort(Wn(idx(2):end)))];

win_freq = ones(N_spk, N_omega_BF);
idx = find((Wn-1));
% figure
for ii = idx
    win_time(ii,:) = fir1(N_omega_BF - 1,Wn(ii));
    temp = fft(win_time(ii,:), N_omega_BF*2);
    win_freq(ii,:) = abs(temp(1:end/2));
%     win_freq(array.N_speaker-ii+1,:) =  win_fre(ii,:);
%     semilogx(array.f,win_freq(ii,:));
%     hold on;
end
H = H.*win_freq;
%
for ii = 1:N_spk
%     [~,idx] = find(abs(H(ii,:))<0.5);
%     H(ii,idx) = 0;
    [~,idx] = find(abs(H(ii,:))>=1);
    H(ii,idx) = H(ii,idx)./abs(H(ii,idx));
end

%% get the filter in time domain
H_conj = [H,conj(H(:,end-1:-1:2))];
for i = 1:N_spk
    h(i,:) = real(ifftshift(ifft(H_conj(i,:))));
%     h(i,:) = ifftshift(real(ifft(H_conj(i,:))));
end


% h = zeros( N_spk, 2 * N_omega_BF-2 );
% for i1 = 1 : size( H, 1 )
%     h( i1, : ) = applyFIRDesign( H( i1, : )', f_BF, fs);
% end
% h = h / max( max( h ) );

%% visualization
nfft = 2^16;
f_TD = 0 : fs / nfft : (nfft - 1) * fs / nfft;
f_TD = f_TD( 1: end/2);
f_num = zeros(1, length(f_BF));
for i1 = 1 : length(f_BF)
    [~, N_tmp] = min(abs(f_TD - f_BF(i1)));
    f_num(i1) = N_tmp;
end

H_rpc = zeros(N_spk, N_omega_BF);
for i1 = 1 : N_spk
    H_tmp = fft(squeeze(h(i1, :)), nfft);
    H_rpc(i1, :) = H_tmp(f_num);
end 
    
p_rpc = zeros(N_theta, N_omega_BF);
for i1 = 1 : N_omega_BF
    G_tmp = squeeze(G0(:, :, i1));
    p_rpc(:, i1) = G_tmp * H_rpc(:, i1);
%     p_rpc(:, i1) = G_tmp * H(:, i1);
end

p_rpc_plot = p_rpc ./ max(abs(p_rpc));
SPL_plot = mag2db(abs(p_rpc_plot));
% meshTheta = repmat(theta',1,N_omega_BF);
% meshF = repmat(f_BF, N_theta,1);
% figure
% surf(meshF,meshTheta,SPL_plot);
% view(0,90);
% set(gca,'xscale','log');
% set(gca, 'FontSize', 14, 'FontName', 'Times New Roma');
% xlim([0,2e4]); ylim([0,360]);
% xlabel('Frequency / Hz'); 
% ylabel('Angle / degree');
% caxis([-20,0]);
% colorbar;
% colormap('jet');
% shading interp;


    
end

