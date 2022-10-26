%% Dewey 2022 JASA EL Supplemental plot code
% fig_s7_8_swept_f2_fixed_dp_ave_plot.m

% Plots average BM or TM vibrations for f2 = ~10-32 kHz, with f1 varied 
% to fix either f2-f1 or 2f1-f2 at ~9kHz, L2 = 60 dB SPL, and 
% L1 = 30-80 dB SPL. Average data are plotted vs. f2 for 
% L1 = 50, 60, and 70 dB SPL.

% All individual data are loaded from 'bm_swept_f2_fixed_qdp.mat', 
% 'bm_swept_f2_fixed_cdp.mat', 'tm_swept_f2_fixed_qdp.mat', or
% 'tm_swept_f2_fixed_cdp.mat', which include
% vibration magnitudes (dB re 1 nm RMS) and phases (cycles).

% Set the location to 'bm' or 'tm' to plot the data from that measurement
% location.

close all;clear;clc;

%% Set measurement location here:
location = 'bm'; % either 'bm' or 'tm'

%% Load relevant data file
switch location 
    case 'bm'
        qdp=load('bm_swept_f2_fixed_qdp.mat');
        cdp=load('bm_swept_f2_fixed_cdp.mat');
    case 'tm'
        qdp=load('tm_swept_f2_fixed_qdp.mat');
        cdp=load('tm_swept_f2_fixed_cdp.mat');
end

f2s = qdp.all.f2s; % f2 (Hz)
L1s = qdp.all.L1s; % L1 (dB SPL)
plot_L1s = 50:10:70; % L1s to plot
min_n = 3; % min # mice with clean data required to plot average

%% Plot max mag
for L_i = 1:length(plot_L1s)
    L1 = plot_L1s(L_i);
    [~,L1_i] = ismember(L1,L1s);

    h=figure('units','normalized','position',[.3 .1 .14 .15]);
    ax1 = axes('position',[.15 .15 .75 .75],'box','off','LineWidth',1.2,'FontSize',13); hold on;
    xlim([9 34]); set(gca,'Xscale','log','Xtick',[10 20 30 40]);
    ylim([-35 20]);
    set(gca,'YTick',-40:20:20);
    ax1.YAxis.MinorTick = 'on';
    ax1.YAxis.MinorTickValues = -50:10:100;
    ax1.XAxis.TickDirection = 'out';
    ax1.XAxis.TickLength = [.015 .015];
    ax1.YAxis.TickLength = [.015 .015];
    
    for i=1:2
        if i==1
            all = cdp.all;
            cx = 'dp2f1_m1f2';
            clr = [0 126 224]/255;
            mrk = 'o';
            mrkClr = clr;

            if L_i==3
                plot([19.5 21.7],[-30 -30], '-', 'LineWidth', 1.4, 'Color', clr);
                plot(20.55,-30, 'Marker', mrk, 'LineWidth', 1.4, 'Color', clr, 'MarkerFaceColor', mrkClr, 'MarkerSize', 6);
            end
        else
            all = qdp.all;
            cx = 'dp1f2_m1f1';
            clr = 'r';
            mrk = 'd';
            mrkClr = 'w';
            if L_i==3
                plot([19.5 21.7],[-22 -22], '-', 'LineWidth', 1.4, 'Color', clr);
                plot(20.55, -22, 'Marker', mrk, 'LineWidth', 1.4, 'Color', clr, 'MarkerFaceColor', mrkClr, 'MarkerSize', 6);
            end
        end
    
        vib_magC = squeeze(all.vib.(genvarname(cx)).magC(:,L1_i,:));
        vib_magC_ave = nanmean(vib_magC,2);
        vib_magC_sd = nanstd(vib_magC,0,2);
        vib_magC_n = sum(~isnan(vib_magC),2);
        vib_magC_se = vib_magC_sd ./ sqrt(vib_magC_n);
        vib_magC_aveC = vib_magC_ave; % cleaned average (only show average when clean data are available from 3 mice
        vib_magC_aveC(vib_magC_n < min_n) = NaN;
    
        % Plot
        errorbar(f2s/1000, vib_magC_aveC, vib_magC_se, '-','Color', clr,'LineWidth', 1, 'CapSize', 4)
        plot(f2s/1000, vib_magC_aveC, 'Color', clr,'LineWidth', 1.4, 'Marker', mrk, 'MarkerFaceColor', mrkClr, 'MarkerSize', 6);
    end
end

