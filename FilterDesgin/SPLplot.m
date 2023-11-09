function [] = SPLplot(p_beam)

load('C:\Users\Administrator\Documents\RoomSimulation\SourceBeamCtrl\beamformer\Data\f_1_3_oct_right.mat');
if ~isempty(find(p_beam == 0))
    p_beam = p_beam + 1e-20 * ones(size(p_beam));
end


SPl_beam = mag2db(abs(p_beam));
theta = -180 : 1 : 179;
theta0 = repmat(theta',1,length(f)); f0 = repmat(f, length(theta),1);
figure
surf(f0,theta0,SPl_beam);
view(0,90);
set(gca,'xscale','log');
xlim([min(f), max(f)]); ylim([-180,180]);
% xlim([100, max(f)]); ylim([-180,180]);
caxis([-20,0]); %zlim([-10,0]);
colorbar;
shading interp;
colormap('jet');

end

