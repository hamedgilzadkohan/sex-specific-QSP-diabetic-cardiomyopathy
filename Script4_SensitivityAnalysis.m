%% ========================================================================
%  Script4_SensitivityAnalysis.m
%  Sex-Specific QSP Model of Diabetic Cardiomyopathy
%  SENSITIVITY ANALYSIS — All 3 Layers, All 6 States
%
%  Authors : Soheili M, Gilzad-Kohan H, Lotfi AS
%
%  Method  : One-at-a-time (OAT) sensitivity analysis
%            Each parameter perturbed ±20% from nominal value
%            Output metric: normalized sensitivity index (SI)
%            SI = (ΔOutput/Output) / (ΔParam/Param)
%
%  Focus states : female_post_T2DM (worst), male_T2DM (comparison)
%  Focus outputs: CFR, IMR, E/e', PGC-1α, Ceramide (primary endpoints)
%
%  Output files:
%    - SA_TornadoPlots.png         (tornado plots, publication figure)
%    - SA_HeatMap.png              (sensitivity heatmap all params)
%    - SA_Results.xlsx                  (full sensitivity table)
%    - SA_RunLog.txt
%
%  Version : 1.0  |  Date: 2026-03-09
%  MATLAB  : R2019b or newer
%% ========================================================================

clear; close all; clc;

script_version = '1.0';
run_timestamp  = datestr(now, 'yyyy-mm-dd HH:MM:SS');
script_path = fileparts(mfilename('fullpath'));
if isempty(script_path), script_path = pwd; end
outdir = script_path;

flog = fopen(fullfile(outdir,'SA_RunLog.txt'),'w');
fprintf(flog,'=======================================================\n');
fprintf(flog,' Sensitivity Analysis — Reproducibility Log\n');
fprintf(flog,'=======================================================\n');
fprintf(flog,' Version   : %s\n', script_version);
fprintf(flog,' Timestamp : %s\n', run_timestamp);
fprintf(flog,' MATLAB    : %s\n', version);
fprintf(flog,' Method    : One-at-a-time (OAT), perturbation = +/-20%%\n');
fprintf(flog,'=======================================================\n\n');

fprintf('=======================================================\n');
fprintf(' Script 4: Sensitivity Analysis\n');
fprintf(' %s\n', run_timestamp);
fprintf('=======================================================\n\n');

%% --- Settings ------------------------------------------------------------
perturb     = 0.20;      % 20% perturbation
focus_states = {'female_post_T2DM','male_T2DM','female_pre_healthy'};
focus_labels = {'Female Post T2DM','Male T2DM','Female Pre Healthy'};
tspan   = [0 168];
options = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.1);

%% --- Define parameters to test -------------------------------------------
% Format: {name, layer, nominal_value, description}
param_list = {
    % Layer 1 kinetic parameters
    'k_PGC_ERb',    1, 0.25,  'ERb->PGC coupling';
    'k_PGC_deg',    1, 0.15,  'PGC-1a degradation';
    'k_Mito_form',  1, 0.12,  'Mito biogenesis rate';
    'k_Mito_deg',   1, 0.05,  'Mito degradation rate';
    'k_ROS_prod',   1, 0.20,  'ROS basal production';
    'k_ROS_scav',   1, 0.30,  'ROS scavenging (MnSOD)';
    'k_Psi_form',   1, 3.00,  'DeltaPsi formation';
    'k_Psi_ROS',    1, 1.50,  'ROS-driven Psi loss';
    % State-specific (tested at female_post_T2DM values)
    'k_ERb_sens',   1, 0.50,  'ERb receptor sensitivity';
    'T2DM_factor',  1, 1.80,  'T2DM stress amplifier';
    'mito_stress',  1, 1.70,  'Mito fission stress';
    'ROS_drive',    1, 2.20,  'ROS production drive';
    'E2_ss',        1, 5.0,   'Estradiol steady state';
    % Layer 2 kinetic parameters
    'k_FAO_AMPK',   2, 0.20,  'AMPK->FAO coupling';
    'k_FAO_PGC',    2, 0.15,  'PGC->FAO coupling';
    'k_Cer_FA',     2, 0.08,  'FA overflow->Ceramide';
    'k_Cer_ROS',    2, 0.10,  'ROS->Ceramide synthesis';
    'k_Col_ROS',    2, 0.06,  'ROS->Collagen crosslink';
    'k_Col_cer',    2, 0.04,  'Ceramide->Fibrosis';
    'k_Ee_Col',     2, 0.10,  'Collagen->E/e stiffness';
    'k_Ee_Cer',     2, 0.06,  'Ceramide->E/e stiffness';
    'k_Ee_dPsi',    2, 0.04,  'dPsi->Diastolic dysfn';
    'T2DM_ins_res', 2, 2.60,  'Insulin resistance factor';
    'FA_overflow',  2, 2.20,  'FA overflow severity';
    'estrogen_prot',2, 0.35,  'Estrogen vascular prot';
    % Layer 3 kinetic parameters
    'k_NO_basal',   3, 0.50,  'eNOS basal NO prod';
    'k_NO_cer_inh', 3, 0.08,  'Ceramide->eNOS inhibition';
    'k_NO_ROS_scav',3, 0.20,  'ROS->NO scavenging';
    'k_ET1_ROS',    3, 0.12,  'ROS->ET1 transcription';
    'k_ET1_cer',    3, 0.08,  'Ceramide->ET1';
    'k_MRes_ET1',   3, 0.15,  'ET1->Vasoconstriction';
    'k_MRes_colx',  3, 0.10,  'Fibrosis->Struct resist';
    'k_CFR_MRes',   3, 0.20,  'Resistance->CFR loss';
    'k_IMR_MRes',   3, 0.12,  'Resistance->IMR rise';
    'k_IMR_Ee',     3, 0.04,  'Filling pressure->IMR';
    'eNOS_activity',3, 0.35,  'eNOS activity (sex/E2)';
    'vasc_estrogen',3, 0.30,  'Vascular estrogen prot';
};

n_params = size(param_list,1);
n_focus  = length(focus_states);

% Output variables to track
out_names = {'CFR','IMR','Ee','PGC1a','Ceramide'};
n_out = length(out_names);

%% --- Run Baseline for each focus state -----------------------------------
fprintf('Running baseline simulations...\n');
baseline = struct();

for fs = 1:n_focus
    sex_state = focus_states{fs};
    [p1,p2,p3] = get_all_params(sex_state);
    baseline(fs).CFR     = run_3layer(p1,p2,p3,tspan,options,'CFR');
    baseline(fs).IMR     = run_3layer(p1,p2,p3,tspan,options,'IMR');
    baseline(fs).Ee      = run_3layer(p1,p2,p3,tspan,options,'Ee');
    baseline(fs).PGC1a   = run_3layer(p1,p2,p3,tspan,options,'PGC1a');
    baseline(fs).Ceramide= run_3layer(p1,p2,p3,tspan,options,'Ceramide');
    fprintf('  Baseline %s: CFR=%.3f IMR=%.1f E/e=%.2f\n',...
        focus_labels{fs},baseline(fs).CFR,baseline(fs).IMR,baseline(fs).Ee);
end

%% --- One-At-A-Time Sensitivity -------------------------------------------
fprintf('\nRunning OAT sensitivity (%d parameters x %d states x 2 perturbations)...\n',...
    n_params, n_focus);

% SI(param, output, state, direction): +1=up, -1=down
SI_up   = zeros(n_params, n_out, n_focus);
SI_down = zeros(n_params, n_out, n_focus);

for i = 1:n_params
    pname  = param_list{i,1};
    player = param_list{i,2};
    pnom   = param_list{i,3};

    for fs = 1:n_focus
        sex_state = focus_states{fs};

        for dir = [1 -1]   % +20% then -20%
            p_new = pnom * (1 + dir*perturb);

            [p1,p2,p3] = get_all_params(sex_state);

            % Apply perturbation to correct layer
            if player == 1
                if isfield(p1, pname), p1.(pname) = p_new; end
            elseif player == 2
                if isfield(p2, pname), p2.(pname) = p_new; end
            else
                if isfield(p3, pname), p3.(pname) = p_new; end
            end

            % Run model
            for oi = 1:n_out
                out_pert = run_3layer(p1,p2,p3,tspan,options,out_names{oi});
                out_base = baseline(fs).(out_names{oi});
                if out_base ~= 0
                    si = ((out_pert - out_base)/out_base) / (dir*perturb);
                else
                    si = 0;
                end
                if dir == 1
                    SI_up(i,oi,fs)   = si;
                else
                    SI_down(i,oi,fs) = si;
                end
            end
        end
    end

    if mod(i,5)==0
        fprintf('  Parameter %d/%d complete\n', i, n_params);
    end
end

fprintf('Sensitivity analysis complete.\n\n');

%% --- Generate Tornado Plots ----------------------------------------------
fprintf('Generating tornado plots...\n');

state_colors_focus = [0.55 0.75 0.95; 0.85 0.20 0.20; 0.20 0.45 0.80];

fig_tornado = figure('Name','Sensitivity Tornado Plots',...
    'Units','inches','Position',[0.5 0.5 16 14],'Color','white');

output_titles = {'CFR (Coronary Flow Reserve)',...
                 'IMR (Microcirculatory Resistance)',...
                 'E/e'' (Diastolic Stiffness)',...
                 'PGC-1\alpha (Mitochondrial Biogenesis)',...
                 'Ceramide (Lipotoxicity)'};

focus_plot = 1;   % female_post_T2DM = index 1
top_n = 10;       % show top 10 parameters

for oi = 1:n_out
    subplot(3,2,oi);

    % Get SI values for this output and focus state
    si_vals = SI_up(:,oi,focus_plot);
    [~, idx] = sort(abs(si_vals),'descend');
    idx = idx(1:min(top_n,n_params));

    si_plot  = si_vals(idx);
    labels_plot = param_list(idx,4);

    % Horizontal bar (tornado)
    colors_bar = arrayfun(@(x) select_color(x), si_plot,'UniformOutput',false);
    colors_bar = vertcat(colors_bar{:});

    h = barh(si_plot, 0.7);
    h.FaceColor = 'flat';
    for b = 1:length(si_plot)
        if si_plot(b) > 0
            h.CData(b,:) = [0.85 0.20 0.20];  % red = increases output
        else
            h.CData(b,:) = [0.20 0.45 0.80];  % blue = decreases output
        end
    end
    set(gca,'YTickLabel', labels_plot,'FontSize',8,'TickDir','out');
    xline(0,'k-','LineWidth',1.0);
    xlabel('Sensitivity Index (SI)','FontSize',9);
    title(output_titles{oi},'FontSize',10,'FontWeight','bold');
    box on;
end

% Legend subplot
subplot(3,2,6); axis off;
patch(NaN,NaN,[0.85 0.20 0.20]); hold on;
patch(NaN,NaN,[0.20 0.45 0.80]);
legend({'Positive SI (parameter increase raises output)',...
        'Negative SI (parameter increase lowers output)'},...
    'Location','best','FontSize',9);
title('State: Female Post T2DM','FontSize',11,'FontWeight','bold');
text(0.5,0.3,sprintf('Perturbation: \\pm%.0f%%', perturb*100),...
    'Units','normalized','HorizontalAlignment','center','FontSize',10);
text(0.5,0.15,'Method: One-At-A-Time (OAT)',...
    'Units','normalized','HorizontalAlignment','center','FontSize',9,...
    'Color',[0.4 0.4 0.4]);

sgtitle({'Sensitivity Analysis — Top 10 Parameters per Output',...
    'Female Post T2DM | QSP Model of Diabetic Cardiomyopathy'},...
    'FontSize',12,'FontWeight','bold');
annotation('textbox',[0 0.96 1 0.04],'String',...
    sprintf('QSP-DCM | %s | v%s', run_timestamp, script_version),...
    'EdgeColor','none','HorizontalAlignment','center',...
    'FontSize',7,'Color',[0.5 0.5 0.5]);

save_figure(fig_tornado, fullfile(outdir,'SA_TornadoPlots'));
fprintf('  Saved: SA_TornadoPlots\n');

%% --- Sensitivity Heatmap -------------------------------------------------
fig_heat = figure('Name','Sensitivity Heatmap',...
    'Units','inches','Position',[0.5 0.5 14 12],'Color','white');

% Use female_post_T2DM SI for heatmap
SI_heat = SI_up(:,:,1);   % params x outputs

% Normalize each column to [-1,1] for display
SI_norm = SI_heat;
for oi = 1:n_out
    maxval = max(abs(SI_heat(:,oi)));
    if maxval > 0
        SI_norm(:,oi) = SI_heat(:,oi) / maxval;
    end
end

imagesc(SI_norm');
colormap(redblue_colormap());
colorbar;
caxis([-1 1]);
set(gca,'XTick',1:n_params,'XTickLabel',param_list(:,4),...
    'XTickLabelRotation',45,'YTick',1:n_out,'YTickLabel',out_names,...
    'FontSize',7,'TickDir','out');
xlabel('Parameters','FontSize',10);
ylabel('Model Outputs','FontSize',10);
title({'Sensitivity Heatmap — Female Post T2DM',...
    'Red = positive sensitivity | Blue = negative sensitivity'},...
    'FontSize',11,'FontWeight','bold');
annotation('textbox',[0 0.96 1 0.04],'String',...
    sprintf('QSP-DCM | %s | v%s', run_timestamp, script_version),...
    'EdgeColor','none','HorizontalAlignment','center',...
    'FontSize',7,'Color',[0.5 0.5 0.5]);

save_figure(fig_heat, fullfile(outdir,'SA_HeatMap'));
fprintf('  Saved: SA_HeatMap\n');

%% --- Export Results Table ------------------------------------------------
fprintf('Exporting sensitivity results...\n');

xlsx_file = fullfile(outdir,'SA_Results.xlsx');

% Sheet 1: Full SI table for Female Post T2DM
param_names_col = param_list(:,1);
param_desc_col  = param_list(:,4);
layer_col       = cell2mat(param_list(:,2));

T_sa = table(param_names_col, param_desc_col, layer_col,...
    SI_up(:,1,1), SI_up(:,2,1), SI_up(:,3,1), SI_up(:,4,1), SI_up(:,5,1),...
    'VariableNames',{'Parameter','Description','Layer',...
    'SI_CFR','SI_IMR','SI_Ee','SI_PGC1a','SI_Ceramide'});

writetable(T_sa, xlsx_file, 'Sheet','Female Post T2DM');

% Sheet 2: Male T2DM
T_sa2 = table(param_names_col, param_desc_col, layer_col,...
    SI_up(:,1,2), SI_up(:,2,2), SI_up(:,3,2), SI_up(:,4,2), SI_up(:,5,2),...
    'VariableNames',{'Parameter','Description','Layer',...
    'SI_CFR','SI_IMR','SI_Ee','SI_PGC1a','SI_Ceramide'});
writetable(T_sa2, xlsx_file, 'Sheet','Male T2DM');

% Metadata
meta = {'Field','Value';
        'Method','One-At-A-Time (OAT)';
        'Perturbation','±20%';
        'SI formula','(dY/Y)/(dP/P)';
        'Focus states','female_post_T2DM, male_T2DM, female_pre_healthy';
        'Outputs','CFR, IMR, E/e, PGC1a, Ceramide';
        'Parameters tested', num2str(n_params);
        'Script Version', script_version;
        'Run Timestamp', run_timestamp};
writecell(meta, xlsx_file, 'Sheet','Metadata');
fprintf('  Saved: %s\n', xlsx_file);

%% --- Print Top 5 Most Sensitive Parameters per Output -------------------
fprintf('\n=======================================================\n');
fprintf(' TOP 5 MOST SENSITIVE PARAMETERS — Female Post T2DM\n');
fprintf('=======================================================\n');
for oi = 1:n_out
    si_vals = SI_up(:,oi,1);
    [sorted_si, idx] = sort(abs(si_vals),'descend');
    fprintf('\n  %s:\n', out_names{oi});
    for k = 1:5
        fprintf('    %d. %-25s  SI = %+.3f\n',...
            k, param_list{idx(k),4}, si_vals(idx(k)));
    end
end
fprintf('\n=======================================================\n');

fprintf(flog,'\nSensitivity analysis complete: %s\n', datestr(now));
fclose(flog);
fprintf('\nScript 4 complete. All sensitivity files saved.\n\n');

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function [p1,p2,p3] = get_all_params(sex_state)
% Returns parameter structs for all 3 layers
    p1 = get_L1_params(sex_state);
    p2 = get_L2_params(sex_state, p1);
    p3 = get_L3_params(sex_state, p2);
end

function out = run_3layer(p1,p2,p3,tspan,options,output_name)
% Runs all 3 layers sequentially and returns the requested output SS value
    options_fast = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',0.5);

    % Layer 1
    y0_1 = [p1.E2_0; p1.ERb_0; p1.PGC1a_0; p1.Mito_0;
            p1.ROS_0; p1.MnSOD_0; p1.DeltaPsi_0];
    try
        [~,y1] = ode15s(@(t,y) L1_odes(t,y,p1), tspan, y0_1, options_fast);
    catch
        out = NaN; return;
    end
    p2.PGC1a  = y1(end,3);
    p2.Mito   = y1(end,4);
    p2.ROS    = y1(end,5);
    p2.dPsi   = y1(end,7);
    p2.MnSOD  = y1(end,6);

    if strcmp(output_name,'PGC1a'), out = y1(end,3); return; end

    % Layer 2
    y0_2 = [p2.AMPK_0; p2.FAO_0; p2.GlucOx_0;
            p2.Cer_0; p2.ColX_0; p2.Ee_0];
    try
        [~,y2] = ode15s(@(t,y) L2_odes(t,y,p2), tspan, y0_2, options_fast);
    catch
        out = NaN; return;
    end
    p3.Ceramide = y2(end,4);
    p3.ColX     = y2(end,5);
    p3.Ee       = y2(end,6);
    p3.ROS      = p2.ROS;
    p3.dPsi     = p2.dPsi;
    p3.PGC1a    = p2.PGC1a;
    p3.Mito     = p2.Mito;

    if strcmp(output_name,'Ee'),      out = y2(end,6); return; end
    if strcmp(output_name,'Ceramide'),out = y2(end,4); return; end

    % Layer 3
    y0_3 = [p3.NO_0; p3.ET1_0; p3.MRes_0;
            p3.CFR_0; p3.IMR_0; p3.WS_0];
    try
        [~,y3] = ode15s(@(t,y) L3_odes(t,y,p3), tspan, y0_3, options_fast);
    catch
        out = NaN; return;
    end

    if strcmp(output_name,'CFR'), out = y3(end,4); return; end
    if strcmp(output_name,'IMR'), out = y3(end,5); return; end
    out = NaN;
end

function c = select_color(val)
    if val >= 0, c = [0.85 0.20 0.20];
    else,        c = [0.20 0.45 0.80];
    end
end

function cmap = redblue_colormap()
    n = 64;
    r = [linspace(0.2,1,n/2), ones(1,n/2)];
    g = [linspace(0.45,1,n/2), linspace(1,0.2,n/2)];
    b = [ones(1,n/2), linspace(1,0.2,n/2)];
    cmap = [r(:), g(:), b(:)];
end

function save_figure(fig, basepath)
    print(fig, basepath, '-dpng', '-r300');
end

%% --- Embedded Layer 1 ODEs (self-contained) ------------------------------
function p = get_L1_params(sex_state)
    p.k_ERb_on=0.05; p.k_ERb_off=0.10; p.k_PGC_basal=0.08;
    p.k_PGC_ERb=0.25; p.k_PGC_deg=0.15; p.k_Mito_form=0.12;
    p.k_Mito_deg=0.05; p.k_ROS_prod=0.20; p.k_ROS_mito=0.15;
    p.k_ROS_scav=0.30; p.k_MnSOD_basal=2.00; p.k_MnSOD_ROS=8.00;
    p.k_MnSOD_PGC=5.00; p.k_MnSOD_deg=0.05; p.k_Psi_form=3.00;
    p.k_Psi_ROS=1.50; p.k_Psi_deg=0.08;
    switch sex_state
        case 'female_pre_healthy'
            p.E2_ss=150;p.k_ERb_sens=1.00;p.T2DM_factor=1.00;
            p.mito_stress=1.00;p.ROS_drive=1.00;
            p.E2_0=150;p.ERb_0=0.20;p.PGC1a_0=0.50;p.Mito_0=30;
            p.ROS_0=1.00;p.MnSOD_0=60;p.DeltaPsi_0=1.00;
        case 'female_pre_T2DM'
            p.E2_ss=150;p.k_ERb_sens=0.65;p.T2DM_factor=1.60;
            p.mito_stress=1.40;p.ROS_drive=1.80;
            p.E2_0=150;p.ERb_0=0.15;p.PGC1a_0=0.40;p.Mito_0=25;
            p.ROS_0=1.50;p.MnSOD_0=70;p.DeltaPsi_0=0.80;
        case 'female_post_healthy'
            p.E2_ss=15;p.k_ERb_sens=0.80;p.T2DM_factor=1.00;
            p.mito_stress=1.20;p.ROS_drive=1.20;
            p.E2_0=15;p.ERb_0=0.15;p.PGC1a_0=0.40;p.Mito_0=25;
            p.ROS_0=1.20;p.MnSOD_0=65;p.DeltaPsi_0=0.90;
        case 'female_post_T2DM'
            p.E2_ss=5;p.k_ERb_sens=0.50;p.T2DM_factor=1.80;
            p.mito_stress=1.70;p.ROS_drive=2.20;
            p.E2_0=5;p.ERb_0=0.10;p.PGC1a_0=0.35;p.Mito_0=20;
            p.ROS_0=1.80;p.MnSOD_0=75;p.DeltaPsi_0=0.70;
        case 'male_healthy'
            p.E2_ss=30;p.k_ERb_sens=0.60;p.T2DM_factor=1.00;
            p.mito_stress=1.00;p.ROS_drive=1.10;
            p.E2_0=30;p.ERb_0=0.12;p.PGC1a_0=0.45;p.Mito_0=28;
            p.ROS_0=1.10;p.MnSOD_0=55;p.DeltaPsi_0=0.95;
        case 'male_T2DM'
            p.E2_ss=30;p.k_ERb_sens=0.45;p.T2DM_factor=1.70;
            p.mito_stress=1.60;p.ROS_drive=2.00;
            p.E2_0=30;p.ERb_0=0.10;p.PGC1a_0=0.35;p.Mito_0=22;
            p.ROS_0=1.70;p.MnSOD_0=72;p.DeltaPsi_0=0.75;
        otherwise
            error('Unknown state: %s',sex_state);
    end
end

function dydt = L1_odes(~,y,p)
    E2=max(y(1),0);ERb=max(y(2),0);PGC1a=max(y(3),0);
    Mito=max(y(4),0);ROS=max(y(5),0);MnSOD=max(y(6),0);DeltaPsi=max(y(7),0);
    dE2=p.k_ERb_on*(p.E2_ss-E2);
    dERb=p.k_ERb_on*p.k_ERb_sens*E2*(1-ERb)-p.k_ERb_off*ERb;
    dPGC1a=p.k_PGC_basal+p.k_PGC_ERb*ERb-p.k_PGC_deg*p.T2DM_factor*PGC1a;
    dMito=p.k_Mito_form*PGC1a*30-p.k_Mito_deg*p.mito_stress*Mito;
    ROS_prod=p.k_ROS_prod*p.ROS_drive+p.k_ROS_mito*p.T2DM_factor*(1/(PGC1a+0.1));
    ROS_scav=p.k_ROS_scav*MnSOD/100*ROS;
    dROS=ROS_prod-ROS_scav;
    dMnSOD=p.k_MnSOD_basal+p.k_MnSOD_ROS*ROS+p.k_MnSOD_PGC*PGC1a-p.k_MnSOD_deg*MnSOD;
    dDeltaPsi=p.k_Psi_form*Mito/30*PGC1a-p.k_Psi_ROS*ROS*DeltaPsi-p.k_Psi_deg*DeltaPsi;
    dydt=[dE2;dERb;dPGC1a;dMito;dROS;dMnSOD;dDeltaPsi];
end

function p = get_L2_params(sex_state, L1)
    p.PGC1a=L1.PGC1a_0; p.Mito=L1.Mito_0; p.ROS=L1.ROS_0;
    p.dPsi=L1.DeltaPsi_0; p.MnSOD=L1.MnSOD_0;
    p.k_AMPK_basal=0.10;p.k_AMPK_dPsi=0.08;p.k_AMPK_ROS=0.06;p.k_AMPK_deg=0.12;
    p.dPsi_ref=50.0;p.k_FAO_AMPK=0.20;p.k_FAO_PGC=0.15;p.k_FAO_deg=0.10;
    p.k_FAO_cer_inh=0.05;p.k_Gluc_basal=0.12;p.k_Gluc_ins_res=0.08;
    p.k_Gluc_deg=0.10;p.k_Cer_FA=0.08;p.k_Cer_ROS=0.10;p.k_Cer_deg=0.04;
    p.Cer_ref=6.0;p.k_Col_ROS=0.06;p.k_Col_cer=0.04;p.k_Col_deg=0.02;
    p.Col_ref=1.0;p.k_Ee_Col=0.10;p.k_Ee_Cer=0.06;p.k_Ee_dPsi=0.04;
    p.k_Ee_deg=0.03;
    switch sex_state
        case 'female_pre_healthy'
            p.T2DM_ins_res=1.00;p.FA_overflow=1.00;p.estrogen_prot=1.00;
            p.AMPK_0=1.0;p.FAO_0=1.0;p.GlucOx_0=1.0;p.Cer_0=4.0;p.ColX_0=1.0;p.Ee_0=7.0;
        case 'female_pre_T2DM'
            p.T2DM_ins_res=2.20;p.FA_overflow=1.80;p.estrogen_prot=0.65;
            p.AMPK_0=0.9;p.FAO_0=1.2;p.GlucOx_0=0.7;p.Cer_0=7.0;p.ColX_0=1.2;p.Ee_0=10.0;
        case 'female_post_healthy'
            p.T2DM_ins_res=1.30;p.FA_overflow=1.30;p.estrogen_prot=0.70;
            p.AMPK_0=0.95;p.FAO_0=1.0;p.GlucOx_0=0.9;p.Cer_0=5.5;p.ColX_0=1.1;p.Ee_0=9.0;
        case 'female_post_T2DM'
            p.T2DM_ins_res=2.60;p.FA_overflow=2.20;p.estrogen_prot=0.35;
            p.AMPK_0=0.7;p.FAO_0=1.5;p.GlucOx_0=0.5;p.Cer_0=10.0;p.ColX_0=1.5;p.Ee_0=13.0;
        case 'male_healthy'
            p.T2DM_ins_res=1.10;p.FA_overflow=1.10;p.estrogen_prot=0.60;
            p.AMPK_0=1.0;p.FAO_0=1.0;p.GlucOx_0=1.0;p.Cer_0=4.5;p.ColX_0=1.0;p.Ee_0=7.5;
        case 'male_T2DM'
            p.T2DM_ins_res=2.10;p.FA_overflow=1.70;p.estrogen_prot=0.45;
            p.AMPK_0=0.8;p.FAO_0=1.3;p.GlucOx_0=0.65;p.Cer_0=8.5;p.ColX_0=1.3;p.Ee_0=11.5;
        otherwise, error('Unknown state: %s',sex_state);
    end
end

function dydt = L2_odes(~,y,p)
    AMPK=max(y(1),0);FAO=max(y(2),0);GlucOx=max(y(3),0);
    Cer=max(y(4),0);ColX=max(y(5),0);Ee=max(y(6),0);
    dPsi_norm=max(1-p.dPsi/p.dPsi_ref,0);
    dAMPK=p.k_AMPK_basal+p.k_AMPK_dPsi*dPsi_norm+p.k_AMPK_ROS*p.ROS-p.k_AMPK_deg*AMPK;
    dFAO=p.k_FAO_AMPK*AMPK*p.FA_overflow+p.k_FAO_PGC*p.PGC1a-p.k_FAO_deg*FAO-p.k_FAO_cer_inh*Cer*FAO;
    dGlucOx=p.k_Gluc_basal-p.k_Gluc_ins_res*p.T2DM_ins_res*GlucOx-p.k_Gluc_deg*GlucOx*FAO/(FAO+0.5);
    dCer=p.k_Cer_FA*p.FA_overflow*(1-p.estrogen_prot*0.4)+p.k_Cer_ROS*p.ROS-p.k_Cer_deg*Cer;
    dColX=p.k_Col_ROS*p.ROS*(1-p.estrogen_prot*0.3)+p.k_Col_cer*Cer/10-p.k_Col_deg*ColX;
    dEe=p.k_Ee_Col*ColX+p.k_Ee_Cer*Cer/10+p.k_Ee_dPsi*dPsi_norm*5-p.k_Ee_deg*Ee;
    dydt=[dAMPK;dFAO;dGlucOx;dCer;dColX;dEe];
end

function p = get_L3_params(sex_state, L2)
    p.Ceramide=L2.Cer_0; p.ROS=L2.ROS; p.ColX=L2.ColX_0;
    p.Ee=L2.Ee_0; p.dPsi=L2.dPsi; p.PGC1a=L2.PGC1a; p.Mito=L2.Mito;
    p.k_NO_basal=0.50;p.k_NO_cer_inh=0.08;p.k_NO_ROS_scav=0.20;p.k_NO_deg=0.15;
    p.k_ET1_basal=0.05;p.k_ET1_ROS=0.12;p.k_ET1_cer=0.08;p.k_ET1_NO_inh=0.10;p.k_ET1_deg=0.08;
    p.k_MRes_ET1=0.15;p.k_MRes_colx=0.10;p.k_MRes_NO_dil=0.20;p.k_MRes_deg=0.05;
    p.k_CFR_MRes=0.20;p.k_CFR_dPsi=0.10;p.k_CFR_deg=0.30;
    p.k_IMR_MRes=0.12;p.k_IMR_Ee=0.04;p.k_IMR_deg=0.03;
    p.k_WS_Ee=0.04;p.k_WS_IMR=0.002;p.k_WS_colx=0.06;p.k_WS_deg=0.03;
    switch sex_state
        case 'female_pre_healthy'
            p.eNOS_activity=1.00;p.vasc_estrogen=1.00;p.ET1_basal_mult=1.00;
            p.NO_0=0.80;p.ET1_0=1.00;p.MRes_0=1.00;p.CFR_0=3.20;p.IMR_0=15.0;p.WS_0=1.00;
        case 'female_pre_T2DM'
            p.eNOS_activity=0.65;p.vasc_estrogen=0.65;p.ET1_basal_mult=1.50;
            p.NO_0=0.55;p.ET1_0=1.40;p.MRes_0=1.30;p.CFR_0=2.40;p.IMR_0=28.0;p.WS_0=1.30;
        case 'female_post_healthy'
            p.eNOS_activity=0.75;p.vasc_estrogen=0.70;p.ET1_basal_mult=1.20;
            p.NO_0=0.65;p.ET1_0=1.15;p.MRes_0=1.15;p.CFR_0=2.80;p.IMR_0=20.0;p.WS_0=1.10;
        case 'female_post_T2DM'
            p.eNOS_activity=0.35;p.vasc_estrogen=0.30;p.ET1_basal_mult=2.00;
            p.NO_0=0.30;p.ET1_0=1.80;p.MRes_0=1.70;p.CFR_0=1.80;p.IMR_0=42.0;p.WS_0=1.70;
        case 'male_healthy'
            p.eNOS_activity=0.80;p.vasc_estrogen=0.60;p.ET1_basal_mult=1.05;
            p.NO_0=0.72;p.ET1_0=1.05;p.MRes_0=1.05;p.CFR_0=3.00;p.IMR_0=16.0;p.WS_0=1.05;
        case 'male_T2DM'
            p.eNOS_activity=0.55;p.vasc_estrogen=0.45;p.ET1_basal_mult=1.60;
            p.NO_0=0.48;p.ET1_0=1.55;p.MRes_0=1.45;p.CFR_0=2.10;p.IMR_0=32.0;p.WS_0=1.45;
        otherwise, error('Unknown state: %s',sex_state);
    end
end

function dydt = L3_odes(~,y,p)
    NO=max(y(1),0.001);ET1=max(y(2),0);MRes=max(y(3),0);
    CFR=max(y(4),0);IMR=max(y(5),0);WS=max(y(6),0);
    dPsi_def=max(1-p.dPsi/50,0);
    dNO=p.k_NO_basal*p.eNOS_activity*p.vasc_estrogen-p.k_NO_cer_inh*p.Ceramide*NO-p.k_NO_ROS_scav*p.ROS*NO-p.k_NO_deg*NO;
    dET1=p.k_ET1_basal*p.ET1_basal_mult+p.k_ET1_ROS*p.ROS+p.k_ET1_cer*p.Ceramide/10-p.k_ET1_NO_inh*NO*ET1-p.k_ET1_deg*ET1;
    dMRes=p.k_MRes_ET1*ET1+p.k_MRes_colx*p.ColX-p.k_MRes_NO_dil*NO*MRes-p.k_MRes_deg*MRes;
    dCFR=-p.k_CFR_MRes*max(MRes-1.0,0)*CFR-p.k_CFR_dPsi*dPsi_def*CFR+p.k_CFR_deg*(3.5-CFR);
    dIMR=p.k_IMR_MRes*MRes+p.k_IMR_Ee*p.Ee-p.k_IMR_deg*IMR;
    dWS=p.k_WS_Ee*p.Ee/10+p.k_WS_IMR*IMR+p.k_WS_colx*p.ColX-p.k_WS_deg*WS;
    dydt=[dNO;dET1;dMRes;dCFR;dIMR;dWS];
end
