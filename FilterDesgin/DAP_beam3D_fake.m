function [ MPS, Theta, Phi, Frequency, theta_step, phi_step, p_beam1  ] = DAP_beam3D(h, h_eq, D)

%% Matrix G
fs = 48e3;
% load('C:\Users\Administrator\Documents\RoomSimulation\SourceBeamCtrl\beamformer\Data\G0_cell.mat');
% G0: N_speaker * 1, N_angle * N_f
global G0;

N_spk = size(G0, 1);
load('C:\Users\Administrator\Documents\RoomSimulation\SourceBeamCtrl\beamformer\Data\f_1_24_oct.mat');
Frequency = f'; N_f = length(f);
theta_step = 1; phi_step = 1;
theta = 0:theta_step:180; phi = 0:phi_step:360-phi_step;
[~, idx1] = find(theta == 0); [~, idx2] = find(theta == 180);
N_theta = length(theta); N_phi = length(phi);
N_angle = (N_theta - length(idx1) - length(idx2)) * N_phi + length(idx1) + length(idx2);
Theta = []; Phi = [];
for ii = 1:N_theta
    if(theta(ii) == 0||theta(ii) == 180)
        Theta = [Theta, theta(ii)];
        Phi = [Phi, 0];
    else
        Theta = [Theta, theta(ii)*ones(1, length(phi))];
        Phi = [Phi, phi];
    end
end

%% filter transformation
fs = 48e3; Ts = 1/fs;
fl = 0; fh = fs/2;
fftsize = 2^15;
f_lin = linspace(fl, fh, fftsize/2+1);
H_tmp = zeros(N_spk, length(f_lin));
H = zeros(N_spk, N_f);
for ii = 1:N_spk
    h_tmp = h(ii,:);
    fft_tmp = fft(h_tmp,fftsize);
    H_tmp(ii,:) = fft_tmp(1:end/2+1);
end
fft_eq_tmp = fft(h_eq,fftsize);
H_eq_tmp = fft_eq_tmp(1:end/2+1);

for ii = 1:N_f
    [~, idx] = min(abs(f_lin - f(ii)));
    H(:,ii) = H_tmp(:,idx) * H_eq_tmp(idx);
end

%% beam formation
Beam = zeros(N_angle, N_f);
for ii = 1:N_spk
    SpkDir = G0{ii,1}.*H(ii,:);
    Beam = Beam + SpkDir;
end

% reverse theta
Beam1 = Beam;
if theta(1) == 0
    Beam1 = [zeros(N_phi-1, N_f); Beam1];
end
if theta(end) == 180
    Beam1 = [Beam1; zeros(N_phi-1, N_f)];
end

Beam2 = zeros(size(Beam1));
for ii = 1:N_theta
    Beam2((ii-1)*N_phi+1:ii*N_phi,:) = Beam1((N_theta-ii)*N_phi+1:(N_theta-ii+1)*N_phi,:);
end

if(theta(end) == 180)
    Beam2(end-N_phi+1:end-1,:) = [];
end
if(theta(1) == 0)
    Beam2(2:N_phi, :) = [];
end

Beam = Beam2./max(abs(Beam2));
MPS = Beam.';

%% octave filter
load('C:\Users\Administrator\Documents\RoomSimulation\SourceBeamCtrl\beamformer\Data\f_1_3_oct_right.mat');
OctfiltFreq = zeros(length(f), length(Frequency));
f_num2 = zeros(1, length(Frequency));
fftnum = 2^15;
f_lin2 = 0 : pi/fftnum : pi - pi/fftnum;
FrequencyRad = Frequency * 2 * pi / fs;
for i1 = 1 : length(Frequency)
    [~, N_tmp] = min(abs(f_lin2 - FrequencyRad(i1)));
    f_num2(i1) = N_tmp;
end

for i1 = 1 : length(f)
    octFilt_tmp = octaveFilter(f(i1),'1/3 octave','SampleRate',fs);
    [octFiltFreq_tmp, ~] = freqz(octFilt_tmp, fftnum);
    OctfiltFreq(i1, :) = octFiltFreq_tmp(f_num2) / max(abs(octFiltFreq_tmp(f_num2)));
end

% figure
% for i1 = 1 : length(f)
%     semilogx(Frequency, abs(OctfiltFreq(i1, :)));
%     hold on
% end
% xlabel('frequency / Hz');
% grid on

Beam3 = abs(Beam) * abs(OctfiltFreq.');
MPS = [];
MPS = Beam3.';
Frequency = f.';
N_f = length(f);

N_theta_D = size(D, 1);
theta_D = 0 : 360 / N_theta_D : 360 - 360 / N_theta_D;
theta_MPS = 90 : theta_step : 270;
N_f_D = size(D, 2);
f_D = linspace(0, fs/2, N_f_D);

for i1 = 1 : N_f
    [~, N_tmp] = min(abs(f_D - Frequency(i1)));
    D_tmp = squeeze(D(:, N_tmp));
    num_st = 1;

    for i2 = 1 : length(theta_MPS)
        [~, N_tmp2] = min(abs(theta_D - theta_MPS(i2)));
        if i2 == 1 || i2 == length(theta_MPS)
            num_end = num_st;
        else
            num_end = num_st + N_phi - 1;
        end
        if D_tmp(N_tmp2) == 0           
            MPS(i1, num_st : num_end) = 0;
        else
            MPS(i1, num_st + 91 / phi_step : phi_step : num_end - 90 / phi_step) = 0;
        end
        num_st = num_end + 1;
    end
end

    


%% visualization
theta1 = [theta(1:end-1), theta(1:end-1)+180];
N_theta1 = length(theta1);
beam = zeros(N_theta1, N_f);

for ii = 1:N_theta
    if(theta(ii) == 0)
        beam(ii,:) = MPS(:, 1).';
    elseif(theta(ii) == 180)
        beam(ii,:) = MPS(:, end).';
    else
        beam(ii,:) = MPS(:, N_phi*(ii-2)+2).';
        beam(N_theta1-ii+2,:) = MPS(:, N_phi*(ii-2)+2+N_phi/2).';
    end
end
tmp = beam(1:270,:);
beam = [beam(271:end,:); tmp];

theta0 = repmat(theta1',1,N_f); f0 = repmat(f, N_theta1,1);

p_beam1 = beam./max(abs(beam));

p_beam = mag2db(abs(p_beam1));
% figure
% surf(f0,theta0-180,p_beam);
% view(0,90);
% set(gca,'xscale','log');
% % xlim([min(f), max(f)]); ylim([-180,180]);
% xlim([100, max(f)]); ylim([-180,180]);
% caxis([-20,0]); %zlim([-10,0]);
% colorbar;
% shading interp;
% colormap('jet');
% 
% 
% 
% %%
% % beam = zeros(N_theta1, N_f);
% % p = G0{24,1};
% % for ii = 1:N_theta
% %     if(theta(ii) == 0)
% %         beam(ii,:) = p(1,:);
% %     elseif(theta(ii) == 180)
% %         beam(ii,:) = p(end,:);
% %     else
% %         beam(ii,:) = p(N_phi*(ii-2)+2,:);
% %         beam(N_theta1-ii+2,:) = p(N_phi*(ii-2)+2+N_phi/2,:);
% %     end
% % end
% % tmp = beam(1:270,:);
% % beam = [beam(271:end,:); tmp];
% % 
% % theta0 = repmat(theta1',1,N_f); f0 = repmat(f, N_theta1,1);
% % p_beam = mag2db(abs(beam./max(abs(beam))));
% % figure
% % surf(f0,theta0-180,p_beam);
% % view(0,90);
% % set(gca,'xscale','log');
% % xlim([min(f), max(f)]); ylim([-180,180]);
% % caxis([-6,0]); %zlim([-10,0]);
% % colorbar;
% % shading interp;
% % colormap('jet');
% 

    
    

end

