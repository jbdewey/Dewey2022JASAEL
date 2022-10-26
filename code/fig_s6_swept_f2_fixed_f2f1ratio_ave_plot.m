%% Dewey 2022 JASA EL Supplemental plot code
% fig_s6_swept_f2_fixed_f2f1ratio_ave_plot.m

% Plots average TM vibrations for with f2 swept and f2/f1 fixed at 1.57, 
% L2 = L1 = 50-70 dB SPL.
% Average data are plotted vs. L2=L1 for each f2/f1 ratio (~1.07-1.67 in 0.1 steps).

% Individual data are loaded from 'tm_swept_f2_fixed_f2f1ratio.mat',
% which include vibration magnitudes (dB re 1 nm RMS) and phases (cycles).

close all;clear;clc;
load tm_swept_f2_fixed_f2f1ratio.mat

f2s = all.f2s;
f2N = length(f2s);
L2s = all.L2s;
L2N = length(L2s);


%% Plotting details
min_n = 3; % min # mice with clean data required for averaging
components = [{'f2'};{'f1'};{'dp2f1_m1f2'};{'dp1f2_m1f1'}]; % components to plot
componentN=length(components);           
component_clrs = [{'k'};{[.6 .6 .6]};{[0 126 224]/255};{'r'}]; % line color        
mrks = [{'s'};{'s'};{'o'};{'d'}]; % marker
mrkClrs = [component_clrs(1);component_clrs(2);{[0 126 224]/255};{'w'}]; % marker color
mrkSzs = [6 6 6 6 6]; % marker size

for L_i = 1:L2N
    %% Figure
    h1 = figure('units','normalized','position',[.1 .1 .14 .18]);
    ax1 = axes('position', [.15 .15 .75 .75]); set(gca,'FontSize',13,'LineWidth',1.2, 'box','off'); hold on;

    axes(ax1); % Magnitude vs. f2
    xlim([1.8 50]); set(gca,'Xscale','log','Xtick',[.1 .2 .5 1 2 5 10 20 50]);
    ylim([-45 30]); set(gca,'Ytick',-40:20:40);
    ax1.XAxis.TickDirection = 'out';
    ax1.XAxis.TickLength = [.012 .012];
    ax1.YAxis.TickLength = [.012 .012];
    ax1.YAxis.MinorTick = 'on';
    ax1.YAxis.MinorTickValues = -50:10:50;

    %% Plot each component
    for c_i = 1:componentN
        cx = components{c_i};
        clr = component_clrs{c_i};
        mrk = mrks{c_i};
        mrkClr = mrkClrs{c_i};

        vib_magC = squeeze(all.vib.(genvarname(cx)).magC(:,L_i,:));
        vib_phiC = squeeze(all.vib.(genvarname(cx)).phiC(:,L_i,:));

        vib_magC_ave = nanmean(vib_magC,2);
        vib_magC_sd = nanstd(vib_magC,0,2);
        vib_magC_n = sum(~isnan(vib_magC),2);
        vib_magC_se = vib_magC_sd ./ sqrt(vib_magC_n);

        vib_magC_aveC = vib_magC_ave; % cleaned average (only show average when clean data are available from 3 mice
        vib_magC_aveC(vib_magC_n < min_n) = NaN;

        %% Clear isolated points
        for f_i = 1:f2N
            if f_i == 1
                if ~isnan(vib_magC_aveC(f_i)) && isnan(vib_magC_aveC(f_i+1))
                    vib_magC_aveC(f_i)=NaN;
                end
            elseif f_i == f2N
                if ~isnan(vib_magC_aveC(f_i)) && isnan(vib_magC_aveC(f_i-1))
                    vib_magC_aveC(f_i)=NaN;
                end                  
            else
                if ~isnan(vib_magC_aveC(f_i)) && isnan(vib_magC_aveC(f_i-1)) && isnan(vib_magC_aveC(f_i+1))
                    vib_magC_aveC(f_i)=NaN;
                end       
            end        
        end

        axes(ax1);
        plot(f2s/1000, vib_magC_aveC, '-', 'LineWidth', 1.4, 'Color', clr);
        errorbar(f2s/1000, vib_magC_aveC, vib_magC_se, '-', 'LineWidth', .8, 'Color', clr, 'CapSize', 4);
        plot(f2s/1000, vib_magC_aveC, '-', 'LineWidth', 1.2, 'Color', clr, 'Marker', mrk, 'MarkerFaceColor', mrkClr,'MarkerSize', 5);
    end
end


