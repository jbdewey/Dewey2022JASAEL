%% Dewey 2022 JASA EL Supplemental plot code
% fig_s3_5_L2_equals_L1_ave_plot.m

% Plots average vibrations for f2 = 9 kHz, L2 = L1 = 20-85 dB SPL, with f1 varied.
% Average data are plotted vs. L2=L1 for each f2/f1 ratio (~1.07-1.67 in 0.1 steps).

% Set the location to 'ohc', 'bm', or 'tm' to plot the data from that measurement
% location.

% Individual data are loaded from 'ohc_fixed_L2_equals_L1.mat',
% 'bm_fixed_L2_equals_L1.mat', or 'tm_fixed_L2_equals_L1.mat', 
% which include vibration magnitudes (dB re 1 nm RMS) and phases (cycles).

close all;clear;clc;

%% Set measurement location here:
location = 'tm'; % either 'ohc', 'bm', or 'tm'

%% Load relevant data file
switch location 
    case 'ohc'
        load 'ohc_L2_equals_L1';
    case 'bm'
        load 'bm_L2_equals_L1';
    case 'tm'
        load 'tm_L2_equals_L1';
end

%% Stimulus parameters 
f1s = all.f1s; % f1 (Hz)
f2s = all.f2s; % f2 (Hz)
L1s = all.L1s; % L1 (dB SPL)
L2s = all.L2s; % L2 (dB SPL)
f1N = length(f1s); % # f1s
L1N = length(L1s); % # L1s


%% Plotting details
min_n = 3; % min # mice with clean data required for averaging
components = [{'f2'};{'f1'};{'dp2f1_m1f2'};{'dp1f2_m1f1'}]; % components to plot
componentN=length(components);           
component_clrs = [{'k'};{[.6 .6 .6]};{[0 126 224]/255};{'r'}]; % line color        
mrks = [{'s'};{'s'};{'o'};{'d'}]; % marker
mrkClrs = [component_clrs(1);component_clrs(2);{[0 126 224]/255};{'w'}]; % marker color
mrkSzs = [6 6 6 6 6]; % marker size

%% Plot component magnitudes vs. L1 for each f2/f1 ratio
for f_i = 1:f1N
    f1=f1s(f_i);
    f2=f2s(f_i);

    % Display component frequencies
    disp(['f2 = ' num2str(f2) ' Hz, f1 = ' num2str(f1) ' Hz, f2/f1 = ' num2str(f2./f1) ]);
    disp(['f2-f1 = ' num2str(f2-f1) ' Hz. 2f1-f2 = ' num2str(2*f1-f2) ' Hz. 2f2-f1 = ' num2str(2*f2-f1) ' Hz.']);
    
    %% Vibration data
    h = figure('units','normalized','position',[.4 .1 .1 .18]);
    ax1 = axes('position', [.15 .15 .75 .75]); set(gca,'FontSize',13,'LineWidth', 1.2, 'box','off'); hold on; % mic f2
    ax2 = axes('position', [.15 .12 .75 .75]); set(gca,'FontSize',13,'LineWidth', 1.2, 'box','off','Color','none'); hold on; % mic f2
 
    axes(ax1);
    xlim([17.5 87.5]); set(gca,'Xtick',-20:20:100);
    ax1.XAxis.MinorTick = 'on';
    ax1.XAxis.MinorTickValues = -20:10:100;
    ax1.XAxis.TickLength = [.015 .015];
    ax1.XAxis.Visible = 'off';
 
    ylim([-40 32.5]); set(gca,'Ytick', -60:20:100); 
    ax1.YAxis.MinorTick = 'on';
    ax1.YAxis.MinorTickValues = -60:10:100;
    ax1.YAxis.TickLength = [.02 .02];
    
    axes(ax2); % x-axis
    xlim([17.5 87.5]); set(gca,'Xtick',-20:20:100);
    ax2.XAxis.MinorTick = 'on';
    ax2.XAxis.MinorTickValues = -20:10:100;
    ax2.XAxis.TickLength = [.02 .02];
    ax2.YAxis.Visible = 'off';
    ax2.XAxis.TickDirection = 'out';

    for c_i = 1:componentN 
        c_x = components{c_i};   
        clr = component_clrs{c_i};
        mrkClr = mrkClrs{c_i};
        mrk = mrks{c_i};
        mrkSz = mrkSzs(c_i);

        %% Average vibrations
        vib_mag = squeeze(all.vib.(genvarname(c_x)).mag(f_i,:,:)); % all data (not used here)
        vib_magC = squeeze(all.vib.(genvarname(c_x)).magC(f_i,:,:)); % clean data

        vib_magC_ave = nanmean(vib_magC,2);
        vib_magC_sd = nanstd(vib_magC,0,2);
        vib_magC_n = sum(~isnan(vib_magC),2)
        vib_magC_se = vib_magC_sd ./ sqrt(vib_magC_n);

        vib_magC_aveC = vib_magC_ave; % cleaned average (only show average when clean data are available from 3 mice
        vib_magC_aveC(vib_magC_n < min_n) = NaN;
 
        %% Clear isolated points
        if f_i>1
            for L_i = 1:L1N
                if L_i == 1
                    if ~isnan(vib_magC_aveC(L_i)) && isnan(vib_magC_aveC(L_i+1))
                        vib_magC_aveC(L_i)=NaN;
                    end
                elseif L_i == L1N
                    if ~isnan(vib_magC_aveC(L_i)) && isnan(vib_magC_aveC(L_i-1))
                        vib_magC_aveC(L_i)=NaN;
                    end                  
                else
                    if ~isnan(vib_magC_aveC(L_i)) && isnan(vib_magC_aveC(L_i-1)) && isnan(vib_magC_aveC(L_i+1))
                        vib_magC_aveC(L_i)=NaN;
                    end       
                end        
            end
        end

        %% Plot
        axes(ax1);
        plot(L1s, vib_magC_aveC, '-', 'LineWidth', 1.6, 'Color', clr);
        errorbar(L1s, vib_magC_aveC, vib_magC_se, '-', 'LineWidth', 1.6, 'Color', clr);
        plot(L1s, vib_magC_aveC, '-', 'LineWidth', 1.4, 'Color', clr, 'Marker', mrk, 'MarkerFaceColor', mrkClr,'MarkerSize', mrkSz); 
    end
    pause; close all;
end




