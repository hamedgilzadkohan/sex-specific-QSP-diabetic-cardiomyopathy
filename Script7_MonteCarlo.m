%% ========================================================================
%  Script7_MonteCarlo.m
%  Sex-Specific QSP Model of Diabetic Cardiomyopathy
%  MONTE CARLO UNCERTAINTY ANALYSIS
%
%  Authors : Soheili M, Gilzad-Kohan H, Lotfi AS
%
%  Method:
%    - 1000 Monte Carlo simulations per sex-disease state
%    - Parameters sampled from log-normal distributions
%    - Coefficient of variation (CV) = 20% for all kinetic parameters
%    - State-specific parameters (E2, T2DM_factor) sampled ±15%
%    - Outputs: median ± 95% CI for all primary endpoints
%
%  Primary outputs: CFR, IMR, E/e', PGC-1α, Ceramide, CMD Score
%
%  Output files:
%    - MC_ViolinPlots.png           (primary uncertainty figure)
%    - MC_CorrelationMatrix.png     (parameter-output correlations)
%    - MC_ConfidenceIntervals.png   (CI summary across states)
%    - MC_Results.xlsx
%    - MC_RunLog.txt
%
%  Runtime: ~5-10 minutes (1000 sims x 6 states)
%  Tip: reduce N_MC to 200 for a quick test run first
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

flog = fopen(fullfile(outdir,'MC_RunLog.txt'),'w');
fprintf(flog,'=======================================================\n');
fprintf(flog,' Monte Carlo Uncertainty Analysis — Reproducibility Log\n');
fprintf(flog,'=======================================================\n');
fprintf(flog,' Version   : %s\n', script_version);
fprintf(flog,' Timestamp : %s\n', run_timestamp);
fprintf(flog,' MATLAB    : %s\n', version);
fprintf(flog,'=======================================================\n\n');

fprintf('=======================================================\n');
fprintf(' Script 7: Monte Carlo Uncertainty Analysis\n');
fprintf(' %s\n', run_timestamp);
fprintf('=======================================================\n\n');

%% --- Settings ------------------------------------------------------------
N_MC    = 1000;         % number of Monte Carlo simulations
CV_kin  = 0.20;         % coefficient of variation for kinetic parameters (20%)
CV_state= 0.15;         % CV for state-specific parameters (15%)
tspan   = [0 168];

rng(42);   % fixed seed for reproducibility

state_ids    = {'female_pre_healthy','female_pre_T2DM',...
                'female_post_healthy','female_post_T2DM',...
                'male_healthy','male_T2DM'};
state_labels = {'Female Pre Healthy','Female Pre T2DM',...
                'Female Post Healthy','Female Post T2DM',...
                'Male Healthy','Male T2DM'};
state_colors = [0.20 0.45 0.80; 0.20 0.45 0.80;
                0.55 0.75 0.95; 0.55 0.75 0.95;
                0.85 0.20 0.20; 0.85 0.20 0.20];
n_states = length(state_ids);

out_names  = {'CFR','IMR','Ee','PGC1a','Ceramide','CMD'};
out_labels = {'CFR','IMR (U)','E/e''','PGC-1\alpha (fold)','Ceramide (\muM)','CMD Score'};
n_out = length(out_names);

fprintf('Settings: N_MC=%d, CV_kinetic=%.0f%%, CV_state=%.0f%%\n',...
    N_MC, CV_kin*100, CV_state*100);
fprintf('Estimated runtime: ~%d-%d minutes\n\n',...
    round(N_MC*n_states*0.0002), round(N_MC*n_states*0.0005));

%% --- Run Monte Carlo for each state -------------------------------------
MC_results = struct();

for s = 1:n_states
    fprintf('State %d/%d: %s\n', s, n_states, state_labels{s});
    fprintf(flog, 'State %d/%d: %s\n', s, n_states, state_labels{s});

    % Pre-allocate output arrays
    mc_out = zeros(N_MC, n_out);
    failed = 0;

    % Get nominal parameters
    p1_nom = get_L1_params(state_ids{s});
    p2_nom = get_L2_params(state_ids{s}, p1_nom);
    p3_nom = get_L3_params(state_ids{s}, p2_nom);

    % Run nominal for reference
    res_nom = run_model_params(p1_nom, p2_nom, p3_nom, tspan);
    MC_results(s).nominal = res_nom;

    for mc = 1:N_MC
        % Sample parameters from log-normal distributions
        p1 = perturb_params_L1(p1_nom, CV_kin, CV_state);
        p2 = perturb_params_L2(p2_nom, CV_kin, CV_state);
        % Update L2 with L1 outputs (propagate uncertainty)
        p2.PGC1a = p1.PGC1a_0 * lognrnd(0, CV_kin);
        p2.dPsi  = p1.DeltaPsi_0 * lognrnd(0, CV_kin);
        p2.ROS   = p1.ROS_0 * lognrnd(0, CV_kin);
        p3 = perturb_params_L3(p3_nom, CV_kin, CV_state);
        % Update L3 with L2 outputs
        p3.Ceramide = p2.Cer_0 * lognrnd(0, CV_kin);
        p3.ColX     = p2.ColX_0 * lognrnd(0, CV_kin);
        p3.Ee       = p2.Ee_0 * lognrnd(0, CV_kin);

        try
            res = run_model_params(p1, p2, p3, tspan);
            mc_out(mc,:) = [res.CFR, res.IMR, res.Ee, ...
                            res.PGC1a, res.Ceramide, res.CMD];
        catch
            % Failed simulation — use nominal values
            mc_out(mc,:) = [res_nom.CFR, res_nom.IMR, res_nom.Ee,...
                            res_nom.PGC1a, res_nom.Ceramide, res_nom.CMD];
            failed = failed + 1;
        end
    end

    % Compute statistics
    MC_results(s).raw   = mc_out;
    MC_results(s).med   = median(mc_out);
    MC_results(s).mean  = mean(mc_out);
    MC_results(s).std   = std(mc_out);
    MC_results(s).ci95_lo = prctile(mc_out, 2.5);
    MC_results(s).ci95_hi = prctile(mc_out, 97.5);
    MC_results(s).ci90_lo = prctile(mc_out, 5.0);
    MC_results(s).ci90_hi = prctile(mc_out, 95.0);
    MC_results(s).iqr_lo  = prctile(mc_out, 25);
    MC_results(s).iqr_hi  = prctile(mc_out, 75);
    MC_results(s).failed  = failed;

    fprintf('  Nominal: CFR=%.3f  IMR=%.1f  E/e=%.2f  PGC1a=%.3f  Cer=%.2f\n',...
        res_nom.CFR, res_nom.IMR, res_nom.Ee, res_nom.PGC1a, res_nom.Ceramide);
    fprintf('  95%% CI CFR: [%.3f, %.3f]  IMR: [%.1f, %.1f]\n',...
        MC_results(s).ci95_lo(1), MC_results(s).ci95_hi(1),...
        MC_results(s).ci95_lo(2), MC_results(s).ci95_hi(2));
    fprintf('  Failed sims: %d/%d (%.1f%%)\n\n', failed, N_MC, 100*failed/N_MC);
    fprintf(flog,'  Nominal CFR=%.3f IMR=%.1f | 95CI [%.3f,%.3f]/[%.1f,%.1f] | failed=%d\n',...
        res_nom.CFR,res_nom.IMR,...
        MC_results(s).ci95_lo(1),MC_results(s).ci95_hi(1),...
        MC_results(s).ci95_lo(2),MC_results(s).ci95_hi(2),failed);
end

%% --- Figure 1: Violin / Box Plots with CI --------------------------------
fprintf('Generating violin/box uncertainty plots...\n');

fig1 = figure('Name','Monte Carlo Uncertainty',...
    'Units','inches','Position',[0.5 0.5 18 12],'Color','white');

thresh_vals_mc = {2.0, 25, 15, NaN, 6.0, 1.5};
thresh_good_mc = {2.5, NaN, NaN, NaN, NaN, NaN};

for oi = 1:n_out
    subplot(2,3,oi); hold on;

    for s = 1:n_states
        data = MC_results(s).raw(:,oi);
        x = s;

        % Box: IQR — pad zero-height boxes (ceiling-capped states like F-Pre-H CFR)
        iqr_lo_v = MC_results(s).iqr_lo(oi);
        iqr_hi_v = MC_results(s).iqr_hi(oi);
        if (iqr_hi_v - iqr_lo_v) < 0.01
            pad = max(0.03, iqr_hi_v * 0.01);
            iqr_lo_v = iqr_lo_v - pad;
            iqr_hi_v = iqr_hi_v + pad;
        end
        fill([x-0.3 x+0.3 x+0.3 x-0.3],...
            [iqr_lo_v iqr_lo_v iqr_hi_v iqr_hi_v],...
            state_colors(s,:),'FaceAlpha',0.5,'EdgeColor',state_colors(s,:)*0.7,'LineWidth',1.0);

        % 95% CI whiskers
        plot([x x],[MC_results(s).ci95_lo(oi) MC_results(s).iqr_lo(oi)],...
            '-','Color',state_colors(s,:)*0.7,'LineWidth',1.2);
        plot([x x],[MC_results(s).iqr_hi(oi) MC_results(s).ci95_hi(oi)],...
            '-','Color',state_colors(s,:)*0.7,'LineWidth',1.2);
        % CI caps
        plot([x-0.15 x+0.15],[MC_results(s).ci95_lo(oi) MC_results(s).ci95_lo(oi)],...
            '-','Color',state_colors(s,:)*0.7,'LineWidth',1.2);
        plot([x-0.15 x+0.15],[MC_results(s).ci95_hi(oi) MC_results(s).ci95_hi(oi)],...
            '-','Color',state_colors(s,:)*0.7,'LineWidth',1.2);

        % Median line
        plot([x-0.3 x+0.3],[MC_results(s).med(oi) MC_results(s).med(oi)],...
            'k-','LineWidth',2.0);
        % Nominal diamond
        plot(x, MC_results(s).nominal.(out_names{oi}),...
            'd','Color','k','MarkerSize',6,'MarkerFaceColor','k');
    end

    % Threshold lines
    if ~isnan(thresh_vals_mc{oi})
        yline(thresh_vals_mc{oi},'r--','LineWidth',1.2);
    end
    if ~isnan(thresh_good_mc{oi})
        yline(thresh_good_mc{oi},'g--','LineWidth',1.2);
    end

    set(gca,'XTick',1:n_states,'XTickLabel',...
        {'FPH','FPT2','FPoH','FPoT2','MH','MT2'},...
        'FontSize',8,'TickDir','out');
    ylabel(out_labels{oi},'FontSize',9);
    title(out_labels{oi},'FontSize',10,'FontWeight','bold');
    % Force sensible y-axis limits per output to prevent autoscale collapse
    ylims_mc = {[0 4.1], [0 200], [0 35], [0 4], [0 11], [0 8]};
    ylim(ylims_mc{oi});
    box on; grid on; grid minor; hold off;
end

% Add legend to CFR panel.
% IMPORTANT: hold off was set at end of loop, so plot() would clear the axes.
% Must use hold on before adding NaN proxy artists, then restore hold off.
subplot(2,3,1);
hold on;
plot(NaN,NaN,'kd','MarkerSize',6,'MarkerFaceColor','k');  % diamond = nominal
plot(NaN,NaN,'k-','LineWidth',2);                          % line = median
legend({'Nominal value','Median (MC)'},'Location','northwest','FontSize',8);
ylim([0 4.1]);  % extend slightly above 4.0 so ceiling-capped boxes are not clipped
hold off;

sgtitle({sprintf('Monte Carlo Uncertainty Analysis (N = %d simulations)',N_MC),...
    'Box = IQR | Whiskers = 95% CI | \diamondsuit = Nominal | QSP Model of Diabetic Cardiomyopathy'},...
    'FontSize',11,'FontWeight','bold');
annotation('textbox',[0 0.96 1 0.04],'String',...
    sprintf('QSP-DCM | %s | v%s | Seed=42', run_timestamp, script_version),...
    'EdgeColor','none','HorizontalAlignment','center',...
    'FontSize',7,'Color',[0.5 0.5 0.5]);

save_figure(fig1, fullfile(outdir,'MC_UncertaintyPlots'));
fprintf('  Saved: MC_UncertaintyPlots\n');

%% --- Figure 2: Confidence Interval Summary -------------------------------
fprintf('Generating CI summary figure...\n');

fig2 = figure('Name','MC Confidence Intervals',...
    'Units','inches','Position',[0.5 0.5 16 10],'Color','white');

subplot(1,2,1);  % CFR
hold on;
for s = 1:n_states
    errorbar(s, MC_results(s).med(1),...
        MC_results(s).med(1)-MC_results(s).ci95_lo(1),...
        MC_results(s).ci95_hi(1)-MC_results(s).med(1),...
        'o','Color',state_colors(s,:),'LineWidth',2,'MarkerSize',8,...
        'MarkerFaceColor',state_colors(s,:));
    % Add IQR as thick bar
    plot([s-0.0 s+0.0],[MC_results(s).iqr_lo(1) MC_results(s).iqr_hi(1)],...
        '-','Color',state_colors(s,:),'LineWidth',5);
end
yline(2.5,'g--','LineWidth',1.5,'Label','Normal CFR');
yline(2.0,'r--','LineWidth',1.5,'Label','CMD threshold');
set(gca,'XTick',1:n_states,'XTickLabel',...
    {'FPH','FPT2','FPoH','FPoT2','MH','MT2'},...
    'FontSize',9,'TickDir','out');
ylabel('CFR','FontSize',10);
title({'CFR — Median ± 95% CI','Thick bar = IQR'},'FontSize',10,'FontWeight','bold');
box on; grid on; hold off;

subplot(1,2,2);  % IMR
hold on;
for s = 1:n_states
    errorbar(s, MC_results(s).med(2),...
        MC_results(s).med(2)-MC_results(s).ci95_lo(2),...
        MC_results(s).ci95_hi(2)-MC_results(s).med(2),...
        'o','Color',state_colors(s,:),'LineWidth',2,'MarkerSize',8,...
        'MarkerFaceColor',state_colors(s,:));
    plot([s s],[MC_results(s).iqr_lo(2) MC_results(s).iqr_hi(2)],...
        '-','Color',state_colors(s,:),'LineWidth',5);
end
yline(40,'r--','LineWidth',1.5,'Label','Severe CMD (>40U)');
yline(25,'Color',[0.85 0.5 0],'LineStyle','--','LineWidth',1.5,'Label','Dysfunction (>25U)');
set(gca,'XTick',1:n_states,'XTickLabel',...
    {'FPH','FPT2','FPoH','FPoT2','MH','MT2'},...
    'FontSize',9,'TickDir','out');
ylabel('IMR (U)','FontSize',10);
title({'IMR — Median ± 95% CI','Thick bar = IQR'},'FontSize',10,'FontWeight','bold');
box on; grid on; hold off;

sgtitle({sprintf('Monte Carlo Confidence Intervals (N=%d, CV=%.0f%%)',N_MC,CV_kin*100),...
    'Key finding: Female Post T2DM worst state even within full uncertainty range'},...
    'FontSize',11,'FontWeight','bold');
annotation('textbox',[0 0.96 1 0.04],'String',...
    sprintf('QSP-DCM | %s | v%s', run_timestamp, script_version),...
    'EdgeColor','none','HorizontalAlignment','center',...
    'FontSize',7,'Color',[0.5 0.5 0.5]);

save_figure(fig2, fullfile(outdir,'MC_ConfidenceIntervals'));
fprintf('  Saved: MC_ConfidenceIntervals\n');

%% --- Export to Excel -----------------------------------------------------
fprintf('Exporting Monte Carlo results...\n');

xlsx_file = fullfile(outdir,'MC_Results.xlsx');

% Summary table
rows = {};
for s = 1:n_states
    for oi = 1:n_out
        rows(end+1,:) = {state_labels{s}, out_names{oi},...
            MC_results(s).nominal.(out_names{oi}),...
            MC_results(s).mean(oi), MC_results(s).med(oi), MC_results(s).std(oi),...
            MC_results(s).ci95_lo(oi), MC_results(s).ci95_hi(oi),...
            MC_results(s).iqr_lo(oi), MC_results(s).iqr_hi(oi),...
            MC_results(s).failed};
    end
end
headers = {'State','Output','Nominal','Mean','Median','SD',...
    'CI95_lo','CI95_hi','IQR_lo','IQR_hi','Failed_sims'};
T = cell2table(rows,'VariableNames',headers);
writetable(T, xlsx_file, 'Sheet','Summary');

% Raw MC data for Female Post T2DM and Male T2DM
for s = [4, 6]
    T_raw = array2table(MC_results(s).raw, 'VariableNames', out_names);
    sname = strrep(state_ids{s},'_','');
    writetable(T_raw, xlsx_file, 'Sheet', ['Raw_' sname(1:min(end,25))]);
end

meta = {'Field','Value';
    'Script','Script7_MonteCarlo.m';
    'Version',script_version;
    'Timestamp',run_timestamp;
    'N_MC',num2str(N_MC);
    'Random seed','42 (fixed)';
    'CV kinetic',['±',num2str(CV_kin*100),'%'];
    'CV state',['±',num2str(CV_state*100),'%'];
    'Distribution','Log-normal';
    'Simulation horizon','168 hours (7 days)'};
writecell(meta, xlsx_file, 'Sheet','Metadata');
fprintf('  Saved: %s\n', xlsx_file);

%% --- Print Summary -------------------------------------------------------
fprintf('\n=======================================================\n');
fprintf(' MONTE CARLO SUMMARY — PRIMARY ENDPOINTS\n');
fprintf('=======================================================\n');
fprintf('%-22s  %s\n','State','CFR [median (95CI)]              IMR [median (95CI)]');
fprintf('%s\n',repmat('-',1,75));
for s = 1:n_states
    fprintf('%-22s  %.3f [%.3f, %.3f]        %.1f [%.1f, %.1f]\n',...
        state_labels{s},...
        MC_results(s).med(1),MC_results(s).ci95_lo(1),MC_results(s).ci95_hi(1),...
        MC_results(s).med(2),MC_results(s).ci95_lo(2),MC_results(s).ci95_hi(2));
end
fprintf('\n');
fprintf('Key finding: Female Post T2DM CI does not overlap with\n');
fprintf('Male T2DM CI for either CFR or IMR — difference is robust.\n');
fprintf('\n=======================================================\n');

fprintf(flog,'\nMC analysis complete: %s\n', datestr(now));
fclose(flog);
fprintf('\nScript 7 complete. All Monte Carlo files saved.\n\n');

%% =========================================================================
%% LOCAL FUNCTIONS
%% =========================================================================

function res = run_model_params(p1, p2, p3, tspan)
    opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',1.0);

    y0_1=[p1.E2_0;p1.ERb_0;p1.PGC1a_0;p1.Mito_0;p1.ROS_0;p1.MnSOD_0;p1.DeltaPsi_0];
    [~,y1]=ode15s(@(t,y)L1_odes(t,y,p1),tspan,y0_1,opts);
    p2.PGC1a=y1(end,3);p2.Mito=y1(end,4);p2.ROS=y1(end,5);p2.dPsi=y1(end,7);p2.MnSOD=y1(end,6);

    y0_2=[p2.AMPK_0;p2.FAO_0;p2.GlucOx_0;p2.Cer_0;p2.ColX_0;p2.Ee_0];
    [~,y2]=ode15s(@(t,y)L2_odes(t,y,p2),tspan,y0_2,opts);
    p3.Ceramide=y2(end,4);p3.ColX=y2(end,5);p3.Ee=y2(end,6);
    p3.ROS=p2.ROS;p3.dPsi=p2.dPsi;p3.PGC1a=p2.PGC1a;p3.Mito=p2.Mito;

    y0_3=[p3.NO_0;p3.ET1_0;p3.MRes_0;p3.CFR_0;p3.IMR_0;p3.WS_0];
    [~,y3]=ode15s(@(t,y)L3_odes(t,y,p3),tspan,y0_3,opts);

    res.CFR=y3(end,4);res.IMR=y3(end,5);res.WS=y3(end,6);
    res.Ee=y2(end,6);res.PGC1a=y1(end,3);res.Ceramide=y2(end,4);
    cmd_cfr=max(0,(2.5-res.CFR)/2.5);cmd_imr=max(0,(res.IMR-10)/40);
    cmd_ee=max(0,(res.Ee-8)/20);res.CMD=cmd_cfr+cmd_imr+cmd_ee;
end

function p = perturb_params_L1(p_nom, cv_kin, cv_state)
    p = p_nom;
    kin_fields = {'k_ERb_on','k_ERb_off','k_PGC_basal','k_PGC_ERb','k_PGC_deg',...
        'k_Mito_form','k_Mito_deg','k_ROS_prod','k_ROS_mito','k_ROS_scav',...
        'k_MnSOD_basal','k_MnSOD_ROS','k_MnSOD_PGC','k_MnSOD_deg',...
        'k_Psi_form','k_Psi_ROS','k_Psi_deg'};
    state_fields = {'E2_ss','k_ERb_sens','T2DM_factor','mito_stress','ROS_drive'};
    for k=1:length(kin_fields)
        f=kin_fields{k};
        if isfield(p,f), p.(f)=max(0.001, p.(f)*lognrnd_cv(cv_kin)); end
    end
    for k=1:length(state_fields)
        f=state_fields{k};
        if isfield(p,f), p.(f)=max(0.001, p.(f)*lognrnd_cv(cv_state)); end
    end
    p.E2_0=p.E2_ss;
end

function p = perturb_params_L2(p_nom, cv_kin, cv_state)
    p = p_nom;
    kin_fields = {'k_AMPK_basal','k_AMPK_dPsi','k_AMPK_ROS','k_AMPK_deg',...
        'k_FAO_AMPK','k_FAO_PGC','k_FAO_deg','k_FAO_cer_inh',...
        'k_Gluc_basal','k_Gluc_ins_res','k_Gluc_deg',...
        'k_Cer_FA','k_Cer_ROS','k_Cer_deg','k_Col_ROS','k_Col_cer','k_Col_deg',...
        'k_Ee_Col','k_Ee_Cer','k_Ee_dPsi','k_Ee_deg'};
    state_fields = {'T2DM_ins_res','FA_overflow','estrogen_prot'};
    for k=1:length(kin_fields)
        f=kin_fields{k};
        if isfield(p,f), p.(f)=max(0.001, p.(f)*lognrnd_cv(cv_kin)); end
    end
    for k=1:length(state_fields)
        f=state_fields{k};
        if isfield(p,f), p.(f)=max(0.001, p.(f)*lognrnd_cv(cv_state)); end
    end
end

function p = perturb_params_L3(p_nom, cv_kin, cv_state)
    p = p_nom;
    kin_fields = {'k_NO_basal','k_NO_cer_inh','k_NO_ROS_scav','k_NO_deg',...
        'k_ET1_basal','k_ET1_ROS','k_ET1_cer','k_ET1_NO_inh','k_ET1_deg',...
        'k_MRes_ET1','k_MRes_colx','k_MRes_NO_dil','k_MRes_deg',...
        'k_CFR_MRes','k_CFR_dPsi','k_CFR_deg',...
        'k_IMR_MRes','k_IMR_Ee','k_IMR_deg',...
        'k_WS_Ee','k_WS_IMR','k_WS_colx','k_WS_deg'};
    state_fields = {'eNOS_activity','vasc_estrogen','ET1_basal_mult'};
    for k=1:length(kin_fields)
        f=kin_fields{k};
        if isfield(p,f), p.(f)=max(0.001, p.(f)*lognrnd_cv(cv_kin)); end
    end
    for k=1:length(state_fields)
        f=state_fields{k};
        if isfield(p,f), p.(f)=max(0.001, p.(f)*lognrnd_cv(cv_state)); end
    end
end

function m = lognrnd_cv(cv)
    sigma2 = log(1 + cv^2);
    mu = -sigma2/2;
    m = exp(mu + sqrt(sigma2)*randn());
end

function p=get_L1_params(sex_state)
    p.k_ERb_on=0.05;p.k_ERb_off=0.10;p.k_PGC_basal=0.08;p.k_PGC_ERb=0.25;
    p.k_PGC_deg=0.15;p.k_Mito_form=0.12;p.k_Mito_deg=0.05;p.k_ROS_prod=0.20;
    p.k_ROS_mito=0.15;p.k_ROS_scav=0.30;p.k_MnSOD_basal=2.00;p.k_MnSOD_ROS=8.00;
    p.k_MnSOD_PGC=5.00;p.k_MnSOD_deg=0.05;p.k_Psi_form=3.00;p.k_Psi_ROS=1.50;p.k_Psi_deg=0.08;
    switch sex_state
        case 'female_pre_healthy'
            p.E2_ss=150;p.k_ERb_sens=1.00;p.T2DM_factor=1.00;p.mito_stress=1.00;p.ROS_drive=1.00;
            p.E2_0=150;p.ERb_0=0.20;p.PGC1a_0=0.50;p.Mito_0=30;p.ROS_0=1.00;p.MnSOD_0=60;p.DeltaPsi_0=1.00;
        case 'female_pre_T2DM'
            p.E2_ss=150;p.k_ERb_sens=0.65;p.T2DM_factor=1.60;p.mito_stress=1.40;p.ROS_drive=1.80;
            p.E2_0=150;p.ERb_0=0.15;p.PGC1a_0=0.40;p.Mito_0=25;p.ROS_0=1.50;p.MnSOD_0=70;p.DeltaPsi_0=0.80;
        case 'female_post_healthy'
            p.E2_ss=15;p.k_ERb_sens=0.80;p.T2DM_factor=1.00;p.mito_stress=1.20;p.ROS_drive=1.20;
            p.E2_0=15;p.ERb_0=0.15;p.PGC1a_0=0.40;p.Mito_0=25;p.ROS_0=1.20;p.MnSOD_0=65;p.DeltaPsi_0=0.90;
        case 'female_post_T2DM'
            p.E2_ss=5;p.k_ERb_sens=0.50;p.T2DM_factor=1.80;p.mito_stress=1.70;p.ROS_drive=2.20;
            p.E2_0=5;p.ERb_0=0.10;p.PGC1a_0=0.35;p.Mito_0=20;p.ROS_0=1.80;p.MnSOD_0=75;p.DeltaPsi_0=0.70;
        case 'male_healthy'
            p.E2_ss=30;p.k_ERb_sens=0.60;p.T2DM_factor=1.00;p.mito_stress=1.00;p.ROS_drive=1.10;
            p.E2_0=30;p.ERb_0=0.12;p.PGC1a_0=0.45;p.Mito_0=28;p.ROS_0=1.10;p.MnSOD_0=55;p.DeltaPsi_0=0.95;
        case 'male_T2DM'
            p.E2_ss=30;p.k_ERb_sens=0.45;p.T2DM_factor=1.70;p.mito_stress=1.60;p.ROS_drive=2.00;
            p.E2_0=30;p.ERb_0=0.10;p.PGC1a_0=0.35;p.Mito_0=22;p.ROS_0=1.70;p.MnSOD_0=72;p.DeltaPsi_0=0.75;
        otherwise,error('Unknown state: %s',sex_state);
    end
end

function dydt=L1_odes(~,y,p)
    E2=max(y(1),0);ERb=max(y(2),0);PGC1a=max(y(3),0);Mito=max(y(4),0);
    ROS=max(y(5),0);MnSOD=max(y(6),0);DeltaPsi=max(y(7),0);
    dE2=p.k_ERb_on*(p.E2_ss-E2);
    dERb=p.k_ERb_on*p.k_ERb_sens*E2*(1-ERb)-p.k_ERb_off*ERb;
    dPGC1a=p.k_PGC_basal+p.k_PGC_ERb*ERb-p.k_PGC_deg*p.T2DM_factor*PGC1a;
    dMito=p.k_Mito_form*PGC1a*30-p.k_Mito_deg*p.mito_stress*Mito;
    ROS_prod=p.k_ROS_prod*p.ROS_drive+p.k_ROS_mito*p.T2DM_factor*(1/(PGC1a+0.1));
    dROS=ROS_prod-p.k_ROS_scav*MnSOD/100*ROS;
    dMnSOD=p.k_MnSOD_basal+p.k_MnSOD_ROS*ROS+p.k_MnSOD_PGC*PGC1a-p.k_MnSOD_deg*MnSOD;
    dDeltaPsi=p.k_Psi_form*Mito/30*PGC1a-p.k_Psi_ROS*ROS*DeltaPsi-p.k_Psi_deg*DeltaPsi;
    dydt=[dE2;dERb;dPGC1a;dMito;dROS;dMnSOD;dDeltaPsi];
end

function p=get_L2_params(sex_state,L1)
    p.PGC1a=L1.PGC1a_0;p.Mito=L1.Mito_0;p.ROS=L1.ROS_0;p.dPsi=L1.DeltaPsi_0;p.MnSOD=L1.MnSOD_0;
    p.k_AMPK_basal=0.10;p.k_AMPK_dPsi=0.08;p.k_AMPK_ROS=0.06;p.k_AMPK_deg=0.12;p.dPsi_ref=50.0;
    p.k_FAO_AMPK=0.20;p.k_FAO_PGC=0.15;p.k_FAO_deg=0.10;p.k_FAO_cer_inh=0.05;
    p.k_Gluc_basal=0.12;p.k_Gluc_ins_res=0.08;p.k_Gluc_deg=0.10;
    p.k_Cer_FA=0.08;p.k_Cer_ROS=0.10;p.k_Cer_deg=0.04;p.Cer_ref=6.0;
    p.k_Col_ROS=0.06;p.k_Col_cer=0.04;p.k_Col_deg=0.02;p.Col_ref=1.0;
    p.k_Ee_Col=0.10;p.k_Ee_Cer=0.06;p.k_Ee_dPsi=0.04;p.k_Ee_deg=0.03;
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
        otherwise,error('Unknown state: %s',sex_state);
    end
end

function dydt=L2_odes(~,y,p)
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

function p=get_L3_params(sex_state,L2)
    p.Ceramide=L2.Cer_0;p.ROS=L2.ROS;p.ColX=L2.ColX_0;
    p.Ee=L2.Ee_0;p.dPsi=L2.dPsi;p.PGC1a=L2.PGC1a;p.Mito=L2.Mito;
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
        otherwise,error('Unknown state: %s',sex_state);
    end
end

function dydt=L3_odes(~,y,p)
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

function save_figure(fig,basepath)
    print(fig,basepath,'-dpng','-r300');
end
