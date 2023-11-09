ccx;
tic;
%%
% addpath 'D:\Documents\RoomSimulation\EQ&D\Beam';
addpath 'D:\Documents\RoomSimulation\SourceBeamCtrl\beamformer'
%%
% G0: N_speaker * 1, N_angle * N_f
if ~exist('G0','var')
    load('D:\Documents\RoomSimulation\SourceBeamCtrl\beamformer\Data\G0_cell.mat');
    global G0;
end
size(G0)

fs = 48e3;
N_Omega = 1024;
%%
% Rcv_Pos = [2.000 1.200 -4.400;
% 6.000     1.5000  -1.800;     
% 10.0000   1.2000  -4.400;
% 2.000     1.2000  -5.400;
% 6.000     1.2000  -5.400;
% 10.0000   1.2000  -5.400;
% 2.000     1.2000  -6.400;
% 6.000     1.2000  -6.400;
% 10.000    1.2000  -6.400;
% 2.0000    1.3111   -7.3000;
% 6.0000    1.3111   -7.3000;
% 10.0000    1.3111   -7.3000;
% 2.0000    1.4222   -8.2000;
% 6.0000    1.4222   -8.2000;
% 10.0000    1.4222   -8.2000;
% 2.0000    1.5333   -9.1000;
% 6.0000    1.5333   -9.1000;
% 10.0000    1.5333   -9.1000;
% 2.0000    1.6444   -10.0000;
% 6.0000    1.6444   -10.0000;
% 10.0000    1.6444   -10.0000;
% 2.0000    1.7556   -10.9000;
% 6.0000    1.7556   -10.9000;
% 10.0000    1.7556   -10.9000;
% 2.0000    1.8667   -11.8000;
% 6.0000    1.8667   -11.8000;
% 10.0000    1.8667   -11.8000;
% 2.0000    1.9778   -12.7000;
% 6.0000    1.9778   -12.7000;
% 10.0000    1.9778   -12.7000;
% 2.0000    2.0889   -13.6000;
% 6.0000    2.0889   -13.6000;
% 10.0000    2.0889   -13.6000;
% 2.0000    2.2000   -14.5000;
% 6.0000    2.2000   -14.5000;
% 10.0000    2.2000   -14.5000;];


%浦口103设置
Rcv_Pos =[4 1.2 -3;
    4 1.2 -4.75;
    4 1.2 -6.5;
    4 1.2 -8.25;
    4 1.2 -10;
    2 1.2 -3;
    2 1.2 -4.75;
    2 1.2 -6.5;
    2 1.2 -8.25;
    2 1.2 -10;
    0 1.2 -3;
    0 1.2 -4.75;
    0 1.2 -6.5;
    0 1.2 -8.25;
    0 1.2 -10;
    -2 1.2 -3;
    -2 1.2 -4.75;
    -2 1.2 -6.5;
    -2 1.2 -8.25;
    -2 1.2 -10;
    -4 1.2 -3;
    -4 1.2 -4.75;
    -4 1.2 -6.5;
    -4 1.2 -8.25;
    -4 1.2 -10;];
Spk_Pos = [3.5, 2 -1.5; -3.5 2 -1.5];

%浦口109设置
% Rcv_Pos = [7.45 1.2 -7.3;
%     7.45 1.2 -10.165;
%     7.45 1.2 -13.03;
%     7.45 1.2 -15.895;
%     7.45 1.2 -18.76;
%     3.725 1.2 -7.3;
%     3.725 1.2 -10.165;
%     3.725 1.2 -13.03;
%     3.725 1.2 -15.895;
%     3.725 1.2 -18.76;
%     0 1.2 -7.3;
%     0 1.2 -10.165;
%     0 1.2 -13.03;
%     0 1.2 -15.895;
%     0 1.2 -18.76;
%     -3.725 1.2 -7.3;
%     -3.725 1.2 -10.165;
%     -3.725 1.2 -13.03;
%     -3.725 1.2 -15.895;
%     -3.725 1.2 -18.76;
%     -7.45 1.2 -7.3;
%     -7.45 1.2 -10.165;
%     -7.45 1.2 -13.03;
%     -7.45 1.2 -15.895;
%     -7.45 1.2 -18.76;
%     6.43 5.2 -16.24;
%     6.43 5.2 -17.89;
%     6.43 5.2 -19.51;
%     3.215 5.2 -16.24;
%     3.215 5.2 -17.89;
%     3.215 5.2 -19.51;
%     0 5.2 -16.24;
%     0 5.2 -17.89;
%     0 5.2 -19.51;
%     -3.215 5.2 -16.24;
%     -3.215 5.2 -17.89;
%     -3.215 5.2 -19.51;
%     -6.43 5.2 -16.24;
%     -6.43 5.2 -17.89;
%     -6.43 5.2 -19.51;];
% Spk_Pos = [5.8, 4.2 -3.8; -5.8 4.2 -3.8];


%  Spk_Pos = [11, 2.45 -0.1; 0.5 2.49 -0.1];  % Case1; 
%Spk_Pos = [2.9 2.4 -0.1; 8.7 2.4 -0.1];  % Case2; 
%Spk_Pos = [3.87 2.4 -0.1; 7.73 2.4 -0.1];  % Case3; 
% Spk_Pos = [0.35 2.4 -3.7; 11.21 2.4 -3.7];  % Case4; Case5;Case6;Case7;Case8;

% point1=[0.42 1.2 -4.4];point2=[0.624 1.2 -6.4];point3=[1.445 2.2 -14.5]; %Case5
%point1=[0.5146 1.2 -4.4];point2=[0.985 1.2 -6.4];point3=[2.89 2.2 -14.5]; %Case6
%point1=[0.6083 1.2 -4.4];point2=[1.346 1.2 -6.4];point3=[4.335 2.2 -14.5]; %Case7
% point1=[0.712 1.2 -4.4];point2=[1.707 1.2 -6.4];point3=[5.78 2.2 -14.5]; %Case8

aera_front = min(abs(Rcv_Pos(:,3)));
room_length = max(abs(Rcv_Pos(:,3)));
Loudspeaker_x = abs(Spk_Pos(1,3));
Loudspeaker_y =Spk_Pos(1,1);
Loudspeaker_z = Spk_Pos(1,2);
area_frontheight = min(abs(Rcv_Pos(:,2)));
area_backheight = max(abs(Rcv_Pos(:,2)));

%D = FilterCal_D(aera_front, room_length, Loudspeaker_x, Loudspeaker_z, area_frontheight, area_backheight, 2, N_Omega);
D = FilterCal_D36(aera_front, room_length, Loudspeaker_x, Loudspeaker_z, area_frontheight, area_backheight, 2, N_Omega); %Case1; Case2; Case3; Case4;
% D = FilterCal_D36_angle(point1, point2, point3,Loudspeaker_x,Loudspeaker_y, Loudspeaker_z, 2, N_Omega);%Case5;Case6;Case7;Case8;


f_D = linspace(0, fs/2, N_Omega);
y_D = linspace(6.4,13.6,9);
%%
D_flag = 1;
D_count = 0 + 180;
Ang_step = 1; Ang_cir = 360;
theta = 0:Ang_step:(Ang_cir-Ang_step);

dif_D_log = {};
D_log = {};

while(D_flag)
    D_flag = 0;
    D_count = D_count + 1;
%     [new_filename, p_beam] = FAKEbeamformer(D, fs, D_count,1);
%     new_filename = 'WFS_DAP_2022_1x1_STEER1_wid1.ms.daff';
%     new_filename = 'Pro_Fil_DAP_2022_1x1_STEER1_wid1.ms.daff';
%     new_filename = 'noProFil_DAP_2022_1x1_STEER1_wid1.ms.daff';
    [new_filename, p_beam] = beamformer(D, fs, D_count,1);
    run('ita_raven_demo.m')
    testpoint_num = size(mono_ir,2);                                               % the number of test points
    for iter_mic = 1:testpoint_num

        temp = mono_ir(1, iter_mic).timeData(2049:2048+round(end/4)) + mono_ir(2, iter_mic).timeData(2049:2048+round(end/4));
        f_lin = linspace(1,fs/2,length(temp)/2+1);
        fft_temp = fft(temp);
        mag_dB = smoothSpectrum(mag2db(abs(fft_temp(1 : end/2 + 1,1))), f_lin', 6);
%         mag_dB = mag2db(abs(fft_temp(1 : end/2 + 1,1)));
        for kk = 1:N_Omega
            [~,freq_id] = min(abs(f_lin - f_D(kk)));
            
%             per_freq(kk,ii) = db2mag(mag_dB(freq_id));
            per_freq(kk,iter_mic) = db2mag(mag_dB(freq_id));
        end
    end
% Every line
    per_num = 3;
    for p = 1:round(testpoint_num/per_num)
        for kk = 1:N_Omega
            mag_D_freq(p,kk) = 0;
                for pp =1:per_num
                    mag_D_freq(p,kk) = mag_D_freq(p,kk) + abs(per_freq(kk,pp + (p - 1)*per_num));
                end
            mag_D_freq(p,kk) = mag_D_freq(p,kk)/per_num;
        end
%         H_stu = linspace(area_frontheight, area_backheight, round(testpoint_num/per_num));
%         z_tmp = Loudspeaker_z(1)-H_stu;
%         x_tmp = Loudspeaker_x(1)-y_D(p);
        H_stu = Rcv_Pos(pp + (p - 1)*per_num, 2);
        z_tmp = Spk_Pos(1,2)-H_stu;
        x_tmp = Spk_Pos(1,3) - Rcv_Pos(pp + (p - 1)*per_num,3);
        theta_tmp(p) = 180-((pi - acos(z_tmp/sqrt(z_tmp^2+x_tmp^2)))*180/pi+90-180); % in degree
        [~, N_theta1_D(p)] = min(abs(theta - theta_tmp(p)));
    end
%250 - 8k
    [~,kk_id1] = min(abs(f_D - 300));
    [~,kk_id2] = min(abs(f_D - 8e3));
    for kk = kk_id1:kk_id2   
%         dif_D_freq(:,kk) = mag_D_freq(:,kk)/median(mag_D_freq(:,kk));
        dif_D_freq(:,kk) = mag_D_freq(:,kk)/mean(mag_D_freq(:,kk));
        angle_id = find(abs(mag2db(dif_D_freq(:,kk))) > 3);
        D_flag = D_flag || ~isempty(angle_id);
        if ~isempty(angle_id)
            disp(['kk = ',num2str(kk)]);
            for dd = 1:length(angle_id)
% %                 sign(dif_D_freq(angle_id(dd),kk))
                D(N_theta1_D(angle_id(dd)),kk) = D(N_theta1_D(angle_id(dd)),kk)/sqrt(dif_D_freq(angle_id(dd),kk));
            end
        end
    end
    
    dif_D_log{D_count-180} = dif_D_freq;
    D_log{D_count-180} = D;
    
%     if D_count>180+3
        disp('No Way..')
        break;
%     end
end
figure
semilogx(f_D(kk_id1:kk_id2), mag2db(abs(dif_D_freq(:,kk_id1:kk_id2))).');
% figure
% semilogx(f_D, db((per_freq.')./mean(per_freq.',2)));
% xlim([250 8e3]);

SPLplot(p_beam);
%%
testpoint_num = size(mono_ir,2);                                               % the number of test points
loudspeaker_num = size(mono_ir,1);
% IRs_minphase = cell(testpoint_num,1);                                       % save cutted minimum phase IRs
% IRs_minphase_windowed = cell(subband_num,1);                                % save IRs after been windowed.

for iter_mic = 1:testpoint_num
    for iter_loudspeaker = 1:loudspeaker_num
        temp = mono_ir(iter_loudspeaker, iter_mic).timeData;

        Data{iter_mic,1}.IR{1,iter_loudspeaker} = temp(2049:2048+round(end/4));
    end
end
IR_length = length(Data{1, 1}.IR{1, 1});                                               % the length of IRs
fs = 48e3;
run('C:\Users\Administrator\Documents\RoomSimulation\EQ&D\EQ\whatIRs_room')
for iter_loudspeaker = 1:loudspeaker_num
    h_eq(:,iter_loudspeaker) = IR_prototype_filter{1,iter_loudspeaker};
    eq_filter(iter_loudspeaker) = itaAudio();
    eq_filter(iter_loudspeaker).samplingRate = fs;
    eq_filter(iter_loudspeaker).signalType = 'energy';
    eq_filter(iter_loudspeaker).timeData = h_eq(:,iter_loudspeaker);
end
%% Sim by Convolution h_eq
for iter_mic = 1:testpoint_num
    for iter_loudspeaker = 1:loudspeaker_num  
         sim_mono_ir(iter_loudspeaker,iter_mic) = conv(mono_ir(iter_loudspeaker, iter_mic),eq_filter(iter_loudspeaker));    
    end
end

%% Assess STI + 声场不均匀度 + 传输频率特性 + 早后期声能比ita_roomacoustics
fs = mono_ir(1, 1).samplingRate;
f = [20;25;31.5;40;50;62.5;80;100;125;155;200;250;315;400;500;630;800;1000;1250;1600;2000;2500;3150;4000;5000;6350;8000;10000;12500;16000;20000];
band_range1 = 9:24;                                                         % Freq range: 125 ~ 4000 Hz
band_range2 = 6:27;                                                         % Freq range: 63 ~ 8000 Hz
freqRange = [f(6) f(27)];
for ii = 1:size(f,1)
    cpb_correction(ii) = 2^((ii-1)/3); % Constant percentage bandwidth
end

for iter_mic = 1:testpoint_num
    ir(iter_mic) = sim_mono_ir(1,iter_mic) + sim_mono_ir(2,iter_mic);
    sti_mic(iter_mic) = STI( ir(iter_mic).timeData, fs );
    conved_C80(iter_mic) = ita_roomacoustics(ir(iter_mic), 'freqRange', freqRange, 'bandsPerOctave', 3, 'C80','D50' );
    oct_ir = ita_filter_fractional_octavebands(ir(iter_mic));
    oct_mag = sum(abs(oct_ir.freqData),1)./cpb_correction;                                      % 1/3 oct magnitude
    oct_mag_points(:,iter_mic) = oct_mag(1,:).';
end
%for
oct_SPL_mag = mean(oct_mag_points,2);
oct_SPL_mean_db = mag2db(mean(oct_SPL_mag(band_range1)));

oct_TFR_db = mag2db(oct_SPL_mag(band_range2,1)) - oct_SPL_mean_db;                                           % Transmission frequency response
xTick = [20, 31.5, 63, 125, 250, 500, 1000, 2000, 4000, 8000, 16000, 20000];
yTick = [-12, -8, -4, 0, 4, 8, 12];
figure
semilogx(f(band_range2), oct_TFR_db, 'LineWidth', 1.5, 'color', 'c');
set(gca, 'xTick', xTick);
set(gca, 'XTickLabel', string(xTick));
set(gca, 'yTick', yTick);
axis([20, 20000, -16, 12]);
xlabel('频率（Hz）');
ylabel('相对声压级（dB）');
title({[ '仿真  ', '传输频率特性(平均数平移)'];['(标注普通教室、合班教室、阶梯教室一级和二级传输频率特性)']});
grid on;
%标注普通教室、合班教室、阶梯教室一级和二级传输频率特性
hold on;
semilogx(f(band_range2), [linspace(-12, -6, 4)'; (-6)*ones(14,1); linspace(-6, -12, 4)'], 'LineWidth', 1.5, 'color', 'r');
text(63, -8, '-6 dB/oct');
semilogx(f(band_range2), [linspace(-15, -6, 4)'; (-6)*ones(14,1); linspace(-6, -15, 4)'], 'LineWidth', 1.5, 'color', 'b');
text(100, -10, '-9 dB/oct');
semilogx(f(band_range2), 4*ones(22,1), 'LineWidth', 1.5, 'color', 'k');
semilogx(f(band_range1), zeros(16,1), 'LineWidth', 1.5, 'color', 'k');
hold off;
legend('实测结果','一级传输频率特性','二级传输频率特性');

oct_SD_1000 = max(db(oct_mag_points(18,:))) - min(db(oct_mag_points(18,:)));        % Sound distribution
oct_SD_4000 = max(db(oct_mag_points(24,:))) - min(db(oct_mag_points(24,:)));
disp(['稳态声场不均匀度（1000 Hz）： ', char(vpa(oct_SD_1000,2)),' dB'])
disp(['稳态声场不均匀度（4000 Hz）： ', char(vpa(oct_SD_4000,2)),' dB'])

the_STI = mean(sti_mic);                                               % STI
disp(['STI： ', char(vpa(the_STI,2))])


%%
toc;
