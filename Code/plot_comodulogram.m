function plot_comodulogram(PAC_mat,high,low,hs,ls)
% This function plot the comodulogram
%   INPUTS:
%   PAC_mat   : PAC values for required frequency
%   high      : High frequency range for plotting the comodulogram
%   low       : Low frequency range for plotting the comodulogram

%% plot comodulogram
% figure;
xticks = low; yticks = high;
pcolor(low,high,PAC_mat);
shading(gca,'interp'); 
colormap(jet);
set(gca,'FontSize',6);
% xlabel('Phase Frequency (Hz)','FontSize',10);ylabel('Amplitude Frequency (Hz)','FontSize',10);
set(gca,'FontName','Arial');
set(gca,'XTick',xticks);

end
