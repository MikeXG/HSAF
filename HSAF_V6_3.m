%% ========================================================================
%% HSAF_V6.m - 基于诊断结果优化的完整Hankel滤波系统
%% ========================================================================
%
% 目标：将残差RMS从62.9%降低到<30%
% 
% 优化策略：
%   1. 增加k_max (24→28) - 捕获更多模态
%   2. 增加win_min (60→90) - 提高窗口稳定性
%   3. 放宽band_scale ([0.5,3.5]→[0.4,4.0]) - 更宽容的带宽
%   4. 增加win_overlap (0.85→0.90) - 更平滑的OLA
%   5. 添加二次迭代选项
%
% 作者：Optimized based on 20251223 diagnosis
% 日期：2025-01-06
% 版本：V6.0 (Iterative + Enhanced Parameters)
%
%% ========================================================================

%% ========================================================================
%% 第一部分：环境初始化
%% ========================================================================

fprintf('\n╔════════════════════════════════════════╗\n');
fprintf('║  HSAF V6.0 - 二次优化迭代版本         ║\n');
fprintf('║  Target: Residual RMS < 30%%           ║\n');
fprintf('╚════════════════════════════════════════╝\n');
fprintf('Start time: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

%% 1.1) 环境检测
is_cluster = ~usejava('desktop');
if is_cluster
    fprintf('[ENV] Running on HPC cluster (nodisplay)\n');
    set(0, 'DefaultFigureVisible', 'off');
else
    fprintf('[ENV] Running on local workstation\n');
end

%% 1.2) 路径配置
base_path = '/home/um202370130/HSAF';
data_path = fullfile(base_path, 'Data');
script_path = fullfile(base_path, 'Software');
output_base = fullfile(data_path, 'EWH_Output', 'CSR_EWH', 'HSAF_Adaptive');

% 添加必要路径
addpath(fullfile(data_path, 'EWH_Output'));
addpath(fullfile(data_path, 'Auxiliary_Data'));
addpath(fullfile(script_path, 'HSA_Filter'));
addpath(fullfile(script_path, 'Tool_Functions', 'gmt'));
addpath(fullfile(script_path, 'Tool_Functions'));
addpath(fullfile(script_path, 'Tool_Functions', 'm_map'));

% 创建输出目录
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
output_dir = fullfile(output_base, timestamp);
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
fprintf('[DIR] Output: %s\n', output_dir);

%% 1.3) 初始化日志
log_file = fullfile(output_dir, 'hsaf_v6.log');
diary(log_file); diary on;

%% ========================================================================
%% 第二部分：数据加载
%% ========================================================================

fprintf('\n[LOAD] Loading GRACE data...\n');
try
    load(fullfile(data_path, 'EWH_Output', 'CSR_EWH_data.mat'), ...
         'CSR_EWH', 'csr_lon', 'csr_lat', 'dateTime');
    fprintf('[LOAD] ✓ Data loaded: %s\n', mat2str(size(CSR_EWH.None)));
catch ME
    error('[ERROR] Failed to load data: %s', ME.message);
end

EWH_all = CSR_EWH.None;
lon = csr_lon(:);
lat = csr_lat(:);

%% 2.2) 时间范围选择
t0 = 19; t1 = 150;  % 2004.01-2016.01 (132 months)
EWH = EWH_all(:,:,t0:t1);
dates = dateTime(t0:t1);
[Nlon, Nlat, Nt] = size(EWH);

fprintf('[DATA] Period: %s to %s\n', ...
        datestr(dates(1), 'yyyy-mm'), datestr(dates(end), 'yyyy-mm'));
fprintf('       Grid: %d lon × %d lat × %d months\n', Nlon, Nlat, Nt);

dlon = mean(diff(lon));
fprintf('       Spatial resolution: %.3f°\n', dlon);

%% ========================================================================
%% 第三部分：参数配置（基于诊断优化）
%% ========================================================================

fprintf('\n[CONFIG] ═══════════════════════════════════════\n');
fprintf('         Parameter Configuration (Optimized)\n');
fprintf('         ═══════════════════════════════════════\n');

opt = struct();

%% ─────────────────────────────────────────────────────────────
%% 【核心参数组1】PSD分析参数
%% ─────────────────────────────────────────────────────────────
opt.psd_win   = 64;          % Welch窗口长度
opt.psd_ov    = 32;          % 重叠点数 (50%)
opt.psd_nfft  = 512;         % FFT点数

fprintf('\n┌─ PSD Parameters ───────────────────────────┐\n');
fprintf('│ Window: %d | Overlap: %d | NFFT: %d       │\n', ...
        opt.psd_win, opt.psd_ov, opt.psd_nfft);
fprintf('└────────────────────────────────────────────┘\n');

%% ─────────────────────────────────────────────────────────────
%% 【核心参数组2】条带主峰识别
%% ─────────────────────────────────────────────────────────────
opt.peak_wl_range_deg = [2, 20];     % 主峰搜索范围
opt.default_wl_peak_deg = 7.5;       % 默认峰值波长

fprintf('\n┌─ Stripe Peak Detection ────────────────────┐\n');
fprintf('│ Search range: [%.1f, %.1f]°              │\n', ...
        opt.peak_wl_range_deg(1), opt.peak_wl_range_deg(2));
fprintf('│ Default peak: %.1f°                        │\n', ...
        opt.default_wl_peak_deg);
fprintf('└────────────────────────────────────────────┘\n');

%% ─────────────────────────────────────────────────────────────
%% 【核心参数组3】嵌入维数P规则 ⚙️ 关键参数
%% ─────────────────────────────────────────────────────────────
% 物理意义：P决定Hankel矩阵的"记忆长度"
% 规则：P = p_factor × (λ_peak / dlon)
%       即：让P覆盖主波长的p_factor个采样周期
% 影响：P↑ → 可捕获更长周期，但计算量↑、易过拟合
%       P↓ → 计算快，但无法分解长周期信号

opt.p_factor = 10;                   % P倍数因子
opt.p_min = 30;                      % P的下界
opt.p_max = 120;                     % P的上界

fprintf('\n┌─ Embedding Dimension P (Critical) ─────────┐\n');
fprintf('│ Formula: P = %.0f × (λ_peak / %.3f°)       │\n', ...
        opt.p_factor, dlon);
fprintf('│ Range: [%d, %d]                            │\n', ...
        opt.p_min, opt.p_max);
fprintf('│ Physical meaning:                           │\n');
fprintf('│   P = "memory" of Hankel matrix             │\n');
fprintf('│   Larger P → capture longer periods         │\n');
fprintf('│              but ↑ computation & overfitting│\n');
fprintf('└────────────────────────────────────────────┘\n');

%% ─────────────────────────────────────────────────────────────
%% 【核心参数组4】模态数K规则 🎯 最关键参数
%% ─────────────────────────────────────────────────────────────
% 物理意义：K是保留的主要频率分量数（成对共轭极点）
% 诊断结果：当前k_max=24不足，残差62.9%
% 优化方案：k_max 24→28，k_min 8→10

opt.k_determination = 'eigenvalue_gap';  % 'eigenvalue_gap' 或 'energy'
opt.energy = 0.85;                   % 能量阈值法的累积能量
opt.k_min  = 10;                     % ✅ 8→10 (防止±30°区域的突降)
opt.k_max  = 38;                     % ✅ 24→28 (诊断建议，捕获更多模态)
opt.force_even_k = true;             % 强制偶数K (保证共轭配对)
opt.sv_gap_threshold = 0.18;         % ✅ 0.20→0.18 (更容易选到k_max)

fprintf('\n┌─ Mode Number K (MOST CRITICAL) ────────────┐\n');
fprintf('│ Method: %s                  │\n', opt.k_determination);
fprintf('│ Range: [%d, %d] (was [8, 24])              │\n', ...
        opt.k_min, opt.k_max);
fprintf('│ SV gap threshold: %.2f (was 0.20)          │\n', ...
        opt.sv_gap_threshold);
fprintf('│                                             │\n');
fprintf('│ 🎯 Physical meaning:                        │\n');
fprintf('│   K = # of frequency components to keep     │\n');
fprintf('│   K↓ → underfitting, stripes remain         │\n');
fprintf('│   K↑ → overfitting, lose real signals       │\n');
fprintf('│                                             │\n');
fprintf('│ 📊 Optimization rationale:                  │\n');
fprintf('│   - Current residual: >50%% (target <30%%)  │\n');
fprintf('│   - Diagnosis suggests: k_max→28            │\n');
fprintf('│   - Strategy: Allow more modes to capture   │\n');
fprintf('│     all stripe wavelengths (4.8-24.4°)      │\n');
fprintf('└────────────────────────────────────────────┘\n');

%% ─────────────────────────────────────────────────────────────
%% 【核心参数组5】波长带宽 📏 直接控制去除范围
%% ─────────────────────────────────────────────────────────────
% 物理意义：[λ_min, λ_max] 是要去除的条带噪声波长范围
% 诊断结果：残差主导波长在4.8-24.4°
% 当前带宽：[3.5, 30]° 已覆盖残差，但效果不佳
% 优化方案：放宽band_scale以更宽容地捕获边缘条带

opt.band_scale = [0.4, 4.0];         % ✅ [0.5,3.5]→[0.4,4.0]
opt.band_abs_deg = [3.0, 35];        % ✅ [3.5,30]→[3.0,35]

fprintf('\n┌─ Wavelength Bandwidth (Direct Control) ────┐\n');
fprintf('│ Adaptive scale: [%.1f, %.1f] × λ_peak      │\n', ...
        opt.band_scale(1), opt.band_scale(2));
fprintf('│ Absolute limit: [%.1f, %.1f]°              │\n', ...
        opt.band_abs_deg(1), opt.band_abs_deg(2));
fprintf('│                                             │\n');
fprintf('│ 📏 Physical meaning:                        │\n');
fprintf('│   This range defines which wavelengths      │\n');
fprintf('│   will be REMOVED as stripe noise           │\n');
fprintf('│                                             │\n');
fprintf('│ 📊 Optimization rationale:                  │\n');
fprintf('│   - Residual wavelengths: 4.8-24.4°         │\n');
fprintf('│   - Old band: [3.5, 30]° (covered but %%62.9)│\n');
fprintf('│   - New band: [3.0, 35]° (wider safety margin)│\n');
fprintf('│   - Strategy: More tolerant capture of      │\n');
fprintf('│     edge stripes (especially at 24.4°)      │\n');
fprintf('└────────────────────────────────────────────┘\n');

%% ─────────────────────────────────────────────────────────────
%% 【核心参数组6】预处理参数
%% ─────────────────────────────────────────────────────────────
opt.taper_alpha  = 0.03;             % Tukey窗锥化系数
opt.detrend_mode = 'constant';       % 'constant'或'linear'
opt.smooth_lat   = 10;               % 纬向平滑窗口(度)

fprintf('\n┌─ Preprocessing ────────────────────────────┐\n');
fprintf('│ Taper: %.2f | Detrend: %s             │\n', ...
        opt.taper_alpha, opt.detrend_mode);
fprintf('│ Latitude smoothing: %d°                    │\n', ...
        opt.smooth_lat);
fprintf('└────────────────────────────────────────────┘\n');

%% ─────────────────────────────────────────────────────────────
%% 【核心参数组7】极区处理
%% ─────────────────────────────────────────────────────────────
opt.polar_lat_threshold = 60;        % |lat|>80°视为极区
opt.polar_fixed_params  = true;      % 极区采用固定参数

fprintf('\n┌─ Polar Region Handling ────────────────────┐\n');
fprintf('│ Threshold: |lat| > %d°                     │\n', ...
        opt.polar_lat_threshold);
fprintf('│ Fixed params: %s                           │\n', ...
        mat2str(opt.polar_fixed_params));
fprintf('└────────────────────────────────────────────┘\n');

%% ─────────────────────────────────────────────────────────────
%% 【核心参数组8】共轭配对
%% ─────────────────────────────────────────────────────────────
opt.force_conjugate_pairs = true;    % 强制扣除共轭模态对
opt.freq_pair_tol = 0.02;            % 配对容差(cycles/deg)
opt.pair_tol = opt.freq_pair_tol;      % 兼容字段：确保local_htls_destripe生效

fprintf('\n┌─ Conjugate Pair Enforcement ───────────────┐\n');
fprintf('│ Force pairs: %s | Tolerance: %.2f       │\n', ...
        mat2str(opt.force_conjugate_pairs), opt.freq_pair_tol);
fprintf('└────────────────────────────────────────────┘\n');

%% ─────────────────────────────────────────────────────────────
%% 【核心参数组9】滑动窗口 🪟 时变信号关键
%% ─────────────────────────────────────────────────────────────
% 物理意义：将信号分段处理，每段视为指数衰减正弦
% 诊断结果：win_min=60可能偏短，无法充分平滑
% 优化方案：win_min 60→90，overlap 0.85→0.90

opt.use_sliding = true;
opt.win_len = [];                    % 自动确定
opt.win_min = 60;                    % ✅ 60→90 (诊断建议)
opt.win_overlap = 0.90;              % ✅ 0.85→0.90 (更平滑OLA)
opt.step = [];                       % 自动计算
opt.circular = true;                 % 经度环形连接
opt.p_cap_ratio = 1/3;             % P在窗口内的上限比例
opt.p_min_win = 16;
opt.k_per_window = true;            % 不逐窗自适应K
opt.k_energy = 0.95;
opt.min_mode_energy_ratio = 0.0;
opt.protect_wl_gt_deg = 25;          % 保护>35°的长波
opt.ola_window = 'hann';             % OLA窗函数类型
opt.ola_tukey_alpha = 0.20;

fprintf('\n┌─ Sliding Window (Quasi-periodicity) ───────┐\n');
fprintf('│ Min window: %d (was 60)                    │\n', opt.win_min);
fprintf('│ Overlap: %.0f%%%% (was 85%%)                 │\n', ...
        opt.win_overlap*100);
fprintf('│ Circular: %s | OLA window: %s          │\n', ...
        mat2str(opt.circular), opt.ola_window);
fprintf('│ Protect wavelength: >%.0f°                 │\n', ...
        opt.protect_wl_gt_deg);
fprintf('│                                             │\n');
fprintf('│ 🪟 Physical meaning:                        │\n');
fprintf('│   Window = segment for local HTLS analysis  │\n');
fprintf('│   Overlap = smooth transition (OLA)         │\n');
fprintf('│   win_min↑ → more stable, less adaptive     │\n');
fprintf('│                                             │\n');
fprintf('│ 📊 Optimization rationale:                  │\n');
fprintf('│   - Increase win_min for better stability   │\n');
fprintf('│   - Increase overlap for smoother OLA       │\n');
fprintf('│   - Protect long waves (>35°) from removal  │\n');
fprintf('└────────────────────────────────────────────┘\n');

%% ─────────────────────────────────────────────────────────────
%% 【新增】二次迭代参数 🔄
%% ─────────────────────────────────────────────────────────────
opt.enable_iteration = true;         % ✅ 启用二次迭代
opt.max_iterations = 5;              % 最大迭代次数
opt.iter_residual_threshold = 0.08;  % 残差<30%时停止

fprintf('\n┌─ Iterative Filtering (NEW) ────────────────┐\n');
fprintf('│ Enabled: %s                                │\n', ...
        mat2str(opt.enable_iteration));
fprintf('│ Max iterations: %d                          │\n', ...
        opt.max_iterations);
fprintf('│ Stop threshold: %.0f%% residual               │\n', ...
        opt.iter_residual_threshold*100);
fprintf('│                                             │\n');
fprintf('│ 🔄 Strategy:                                │\n');
fprintf('│   1st pass: Remove bulk of stripes          │\n');
fprintf('│   2nd pass: Target remaining fine stripes   │\n');
fprintf('│   Advantage: Handles nonlinear residuals    │\n');
fprintf('└────────────────────────────────────────────┘\n');

fprintf('\n[CONFIG] ═══════════════════════════════════════\n');
fprintf('         Parameter Summary\n');
fprintf('         ═══════════════════════════════════════\n');
fprintf('Key optimizations vs. previous version:\n');
fprintf('  • k_max:      24 → %d (more modes)\n', opt.k_max);
fprintf('  • band:       [3.5,30]° → [%.1f,%.1f]°\n', ...
        opt.band_abs_deg(1), opt.band_abs_deg(2));
fprintf('  • win_min:    60 → %d (more stable)\n', opt.win_min);
fprintf('  • overlap:    85%% → %.0f%% (smoother)\n', opt.win_overlap*100);
fprintf('  • iteration:  OFF → ON (2-pass filtering)\n');
fprintf('\nExpected outcome: Residual RMS 62.9%% → <30%%\n');
fprintf('═══════════════════════════════════════════════\n\n');

%% ========================================================================
%% 第四部分：参数标定
%% ========================================================================

fprintf('\n[CALIB] Calibrating latitude-adaptive parameters...\n');
nCal = min(10, Nt);
t_list = unique(round(linspace(1, Nt, nCal)));
fprintf('[CALIB] Using %d sample months: %s\n', numel(t_list), mat2str(t_list));

tic;
param = calibrate_lat_params_v5(EWH, lon, lat, t_list, opt);
t_calib = toc;
fprintf('[CALIB] ✓ Calibration completed in %.1f sec\n', t_calib);

%% 可视化参数
fprintf('\n[PLOT] Generating parameter plots...\n');
fig = figure('Position', [100 100 1000 700], 'Visible', 'off');

subplot(4,1,1);
plot(lat, param.wl_peak_deg, 'b-', 'LineWidth', 1.5); grid on;
ylabel('λ_{peak} (°)', 'FontSize', 11, 'FontWeight', 'bold');
title('Stripe Peak Wavelength', 'FontSize', 12);
xlim([-90 90]); ylim([0 30]);

subplot(4,1,2);
plot(lat, param.p, 'r-', 'LineWidth', 1.5); grid on;
ylabel('p (embed dim)', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('P range: [%d, %d]', min(param.p), max(param.p)), 'FontSize', 10);
xlim([-90 90]); ylim([0 130]);

subplot(4,1,3);
plot(lat, param.k, 'g-', 'LineWidth', 1.5); grid on;
ylabel('k (#modes)', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('K range: [%d, %d] (optimized)', min(param.k), max(param.k)), ...
      'FontSize', 10);
xlim([-90 90]); ylim([0, max(param.k)+2]);

subplot(4,1,4);
plot(lat, param.wl_band_deg(:,1), 'c--', 'LineWidth', 1.2, 'DisplayName', 'Lower bound');
hold on;
plot(lat, param.wl_band_deg(:,2), 'm--', 'LineWidth', 1.2, 'DisplayName', 'Upper bound');
grid on; ylabel('λ_{band} (°)', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Latitude (°)', 'FontSize', 11, 'FontWeight', 'bold');
legend('Location', 'best');
xlim([-90 90]); ylim([0 40]);

param_fig = fullfile(output_dir, sprintf('params_%s.png', timestamp));
print(fig, param_fig, '-dpng', '-r300');
fprintf('[PLOT] ✓ Parameter plot: %s\n', param_fig);
close(fig);

%% ========================================================================
%% 第五部分：并行计算配置
%% ========================================================================

fprintf('\n[PARALLEL] Configuring parallel pool...\n');
max_workers = 56;

poolobj = gcp('nocreate');
if isempty(poolobj)
    try
        parpool('local', min(max_workers, feature('numcores')));
        fprintf('[PARALLEL] ✓ Started %d workers\n', gcp().NumWorkers);
    catch ME
        warning('[WARN] Parallel pool failed: %s', ME.message);
        fprintf('[PARALLEL] Running in serial mode\n');
    end
else
    fprintf('[PARALLEL] ✓ Using existing pool (%d workers)\n', poolobj.NumWorkers);
end

%% ========================================================================
%% 第六部分：迭代滤波主循环
%% ========================================================================

fprintf('\n[FILTER] ═══════════════════════════════════════\n');
fprintf('         Starting Iterative HSAF Filtering\n');
fprintf('         ═══════════════════════════════════════\n');

EWH_current = EWH;
iteration_stats = struct();

for iter = 1:opt.max_iterations
    fprintf('\n┌────────────────────────────────────────────┐\n');
    fprintf('│ ITERATION %d / %d                            │\n', ...
            iter, opt.max_iterations);
    fprintf('└────────────────────────────────────────────┘\n');
    
    t_iter_start = tic;
    EWH_filt = nan(size(EWH_current), 'like', EWH_current);
    
    %% 逐月逐纬度滤波
    for t = 1:Nt
        t_month_start = tic;
        
        Xt = EWH_current(:,:,t);
        Xo = nan(size(Xt), 'like', Xt);
        
        parfor j = 1:Nlat
            x = Xt(:,j);
            if all(isnan(x))
                Xo(:,j) = x;
                continue;
            end
            
            p = param.p(j);
            k = param.k(j);
            wl_band_deg = param.wl_band_deg(j,:);
            if iter >= 2
                shrink = 0.85;  % 二次迭代：收窄带宽，减少斑点并专治残余细条带
                mid = mean(wl_band_deg);
                half = (wl_band_deg(2)-wl_band_deg(1))/2 * shrink;
                wl_band_deg = [mid-half, mid+half];
            end
            Ts = dlon;
            
            Xo(:,j) = hankel_destripe_oneprofile_sw(x, Ts, p, k, wl_band_deg, opt);
        end
        
        EWH_filt(:,:,t) = Xo;
        
        % 进度报告
        if mod(t, 20) == 0 || t == 1 || t == Nt
            elapsed = toc(t_iter_start);
            eta = elapsed / t * (Nt - t);
            fprintf('[Iter%d] Month %3d/%d (%.1f%%) | ETA: %.1fmin\n', ...
                    iter, t, Nt, 100*t/Nt, eta/60);
        end
    end
    
    t_iter_total = toc(t_iter_start);
    fprintf('[Iter%d] ✓ Completed in %.1f sec (%.2f sec/month)\n', ...
            iter, t_iter_total, t_iter_total/Nt);
    
    %% 计算残差统计
    residual = EWH_current - EWH_filt;
    t_check = round(Nt/2);
    
    residual_rms = sqrt(mean(residual(:,:,t_check).^2, 'all', 'omitnan'));
    original_rms = sqrt(mean(EWH(:,:,t_check).^2, 'all', 'omitnan'));
    residual_pct = 100 * residual_rms / original_rms;
    
    fprintf('[Iter%d] Residual RMS: %.2f cm (%.1f%% of original)\n', ...
            iter, residual_rms, residual_pct);
    
    % 保存迭代统计
    iteration_stats(iter).residual_rms = residual_rms;
    iteration_stats(iter).residual_pct = residual_pct;
    iteration_stats(iter).time_sec = t_iter_total;
    
    % 检查是否达标
    if residual_pct < opt.iter_residual_threshold * 100
        fprintf('[Iter%d] 🎉 SUCCESS: Residual < %.0f%%, stopping iteration\n', ...
                iter, opt.iter_residual_threshold*100);
        break;
    end
    
    % 准备下一次迭代
    if iter < opt.max_iterations
        fprintf('[Iter%d] ⚠️  Residual still %.1f%%, proceeding to Iteration %d...\n', ...
                iter, residual_pct, iter+1);
        EWH_current = EWH_filt;  % 使用滤波结果作为下次输入
    else
        fprintf('[Iter%d] ⚠️  Max iterations reached, residual = %.1f%%\n', ...
                iter, residual_pct);
    end
end

fprintf('\n[FILTER] ═══════════════════════════════════════\n');
fprintf('Iterative filtering completed.\n');
fprintf('Final residual: %.1f%% (target: <30%%)\n', residual_pct);
if residual_pct < 30
    fprintf('✅ TARGET ACHIEVED!\n');
else
    fprintf('⚠️  Target not fully met, consider:\n');
    fprintf('   - Further increase k_max\n');
    fprintf('   - Manual inspection of residual patterns\n');
end
fprintf('═══════════════════════════════════════════════\n\n');

%% 回填到原始维度
out = nan(size(EWH_all), 'like', EWH_all);
out(:,:,t0:t1) = EWH_filt;
CSR_EWH.HankelLatAdaptive = out;

%% ========================================================================
%% 第七部分：质量评估与诊断
%% ========================================================================

fprintf('\n[DIAGNOSE] Analyzing final residual characteristics...\n');

% 创建诊断日志
diag_log = fullfile(output_dir, sprintf('final_diagnosis_%s.log', timestamp));
diary(diag_log); diary on;

% 选择代表性纬度
lat_test = [-60, -30, 0, 30, 60];
lat_idx_valid = arrayfun(@(x) find(abs(lat - x) < 1, 1), lat_test);

% 创建诊断图
fig_diag = figure('Position', [100 100 1600 1000], 'Visible', 'off');

for i = 1:length(lat_test)
    lat_val = lat_test(i);
    lat_idx = lat_idx_valid(i);
    
    x_resid = residual(:, lat_idx, t_check);
    x_resid = fillmissing(x_resid, 'linear', 'EndValues', 'nearest');
    x_resid = x_resid(:) - mean(x_resid);
    
    % 功率谱
    [pxx, f] = pwelch(x_resid, hamming(64), 32, 512, 1/dlon);
    wl = 1 ./ f; wl(1) = inf;
    mask = (wl >= 2) & (wl <= 60) & isfinite(pxx);
    wl_plot = wl(mask);
    pxx_plot = pxx(mask);
    
    subplot(2, 5, i);
    semilogx(wl_plot, 10*log10(pxx_plot + eps), 'LineWidth', 1.5);
    xlabel('Wavelength (°)'); ylabel('Power (dB)');
    title(sprintf('Lat=%d° Residual Spectrum', lat_val), 'FontSize', 10);
    grid on; xlim([2 60]);
    
    subplot(2, 5, i+5);
    plot(lon, EWH(:, lat_idx, t_check), 'k-', 'DisplayName', 'Before');
    hold on;
    plot(lon, EWH_filt(:, lat_idx, t_check), 'b-', 'LineWidth', 1.5, ...
         'DisplayName', 'After');
    plot(lon, x_resid, 'r--', 'DisplayName', 'Residual');
    xlabel('Longitude (°)'); ylabel('EWH (cm)');
    title(sprintf('Lat=%d° Profiles', lat_val), 'FontSize', 10);
    legend('Location', 'best', 'FontSize', 8); grid on;
end

diag_fig = fullfile(output_dir, sprintf('final_diagnosis_%s.png', timestamp));
print(fig_diag, diag_fig, '-dpng', '-r300');
fprintf('[DIAGNOSE] ✓ Diagnosis plot: %s\n', diag_fig);
close(fig_diag);

diary off;

%% ========================================================================
%% 第八部分：数据保存
%% ========================================================================

fprintf('\n[SAVE] Saving results...\n');

% 主数据文件
save_file = fullfile(output_dir, 'CSR_EWH_HSAF_filtered_optimized.mat');
save(save_file, 'CSR_EWH', 'param', 'opt', 'iteration_stats', ...
     'dateTime', 'dates', 'csr_lon', 'csr_lat', '-v7.3');
fprintf('[SAVE] ✓ Main data: %s (%.1f MB)\n', save_file, ...
        dir(save_file).bytes / 1024^2);

% 统计表格
stats_tbl = struct2table(iteration_stats);
stats_csv = fullfile(output_dir, sprintf('iteration_stats_%s.csv', timestamp));
writetable(stats_tbl, stats_csv);
fprintf('[SAVE] ✓ Statistics: %s\n', stats_csv);

% 脚本快照
src = mfilename('fullpath');
if ~isempty(src)
    dst = fullfile(output_dir, 'script_snapshot.m');
    copyfile([src '.m'], dst);
    fprintf('[SAVE] ✓ Script snapshot: %s\n', dst);
end


fprintf('[DATA] Filtered data saved to CSR_EWH.HankelLatAdaptive\n');

%% ============================================================
%% 8) 质量检查图（保存3个代表性月份）
%% ============================================================
fprintf('\n[PLOT] Generating result check plots...\n');
check_months = [1, round(Nt/2), Nt];  % 首、中、末月份

for idx = 1:length(check_months)
    t_check = check_months(idx);

    fig = figure('Position', [100 100 1400 500], 'Visible', 'off');
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot_map(EWH(:,:,t_check), lon, lat, 1);
    title(sprintf('Before (Month %d)', t0+t_check-1), 'FontSize', 12, 'FontWeight', 'bold');
    clim([-15 15]); colorbar;

    nexttile;
    plot_map(EWH_filt(:,:,t_check), lon, lat, 1);
    title('After HSAF (V6)', 'FontSize', 12, 'FontWeight', 'bold');
    clim([-15 15]); colorbar;

    nexttile;
    plot_map(EWH(:,:,t_check) - EWH_filt(:,:,t_check), lon, lat, 1);
    title('Removed Stripe Noise', 'FontSize', 12, 'FontWeight', 'bold');
    clim([-15 15]); colorbar;

    check_fig = fullfile(output_dir, sprintf('result_check_month%03d_%s.png', ...
                                             t0+t_check-1, timestamp));
    print(fig, check_fig, '-dpng', '-r300');
    fprintf('[PLOT] Check plot saved: %s\n', check_fig);
    close(fig);
end

%% ============================================================
%% 9) 球谐系数转换（可选）
%% ============================================================
fprintf('\n[SHC] Converting to spherical harmonic coefficients...\n');
try
    SC_Hankel = zeros(61, 121, Nt);

    for t = 1:Nt
        % EWH (cm) -> m -> SHC
        cs = gmt_grid2cs(EWH_filt(:,:,t)' / 100, 60);
        gc = gmt_mc2gc(cs);
        sc = gmt_cs2sc(gc);
        SC_Hankel(:,:,t) = sc;

        if mod(t, 20) == 0
            fprintf('[SHC] Processed %d/%d months\n', t, Nt);
        end
    end

    SC_Hankel(SC_Hankel == 0) = NaN;
    fprintf('[SHC] SHC conversion completed.\n');
catch ME
    warning('[WARN] SHC conversion failed: %s', ME.message);
    SC_Hankel = [];
end

%% ============================================================
%% 10) 保存结果
%% ============================================================
fprintf('\n[SAVE] Saving filtered data...\n');
save_file = fullfile(output_dir, 'CSR_EWH_HSAF_filtered.mat');

save(save_file, 'CSR_EWH', 'SC_Hankel', 'param', 'opt', ...
     'dateTime', 'dates', 'csr_lon', 'csr_lat', '-v7.3');

fprintf('[SAVE] Data saved: %s (%.1f MB)\n', save_file, ...
        dir(save_file).bytes / 1024^2);

% ========================================================================
% 残余条带频谱诊断（健壮版 - 添加数据有效性检查）
% ========================================================================
fprintf('\n[DIAGNOSE] Analyzing residual stripe characteristics...\n');

if ~exist(output_dir, 'dir'); mkdir(output_dir); end

log_file = fullfile(output_dir, sprintf('residual_diagnosis_%s.log', timestamp));
if exist(log_file, 'file'); delete(log_file); end
diary(log_file); diary on;
fprintf('[LOG] Console output -> %s\n', log_file);

% 计算残余场
residual = EWH - EWH_filt;

% 检查残余场是否有效
if all(isnan(residual(:)))
    error('[ERROR] Residual field is all NaN. Check if filtering was successful.');
end

% 选择中间月份
t_check = round(Nt/2);
fprintf('[DIAGNOSE] Using Month %d (index %d) for analysis\n', t0+t_check-1, t_check);

% 自动选择3个有数据的典型纬度
lat_candidates = [-60, -30, 0, 30, 60];
lat_test = [];
lat_idx_valid = [];

for lat_val = lat_candidates
    idx = find(abs(lat - lat_val) < 1, 1);  % 容差1度
    if ~isempty(idx)
        x_test = residual(:, idx, t_check);
        % 检查是否有足够的有效数据（至少50%非NaN）
        if sum(~isnan(x_test)) > length(x_test)/2
            lat_test(end+1) = lat_val; %#ok<AGROW>
            lat_idx_valid(end+1) = idx; 
        end
    end
end

if isempty(lat_test)
    warning('[WARN] No valid latitudes found for diagnosis. Using fallback method.');
    % 回退方案：使用全球RMS统计
    rms_all = squeeze(std(residual(:,:,t_check), 0, 1, 'omitnan'));
    [~, best_idx] = max(rms_all);
    lat_test = lat(best_idx);
    lat_idx_valid = best_idx;
end

fprintf('[DIAGNOSE] Selected latitudes for analysis: %s\n', mat2str(lat_test));

% 创建诊断图
fig_diag = figure('Position', [100 100 1400 900], 'Visible', 'off');
colors = lines(length(lat_test));

% 初始化存储主导波长
dominant_wl = zeros(length(lat_test), 3);

for i = 1:length(lat_test)
    lat_val = lat_test(i);
    lat_idx = lat_idx_valid(i);

    % 提取残余剖面
    x_resid = residual(:, lat_idx, t_check);

    % 数据清洗：填补NaN（如果有少量缺失）
    if any(isnan(x_resid))
        x_resid = fillmissing(x_resid, 'linear', 'EndValues', 'nearest');
    end

    % 去均值（PWelch要求）
    x_resid = x_resid(:) - mean(x_resid);

    % 再次检查（防御性编程）
    if all(isnan(x_resid)) || numel(x_resid) < 10
        warning('[WARN] Lat=%d° has invalid data, skipping...', lat_val);
        continue;
    end

    % 计算功率谱
    try
        [pxx, f] = pwelch(x_resid, hamming(64), 32, 512, 1/dlon);
    catch ME
        warning('[WARN] PWelch failed for Lat=%d°: %s', lat_val, ME.message);
        continue;
    end

    wl = 1 ./ f;  % 波长（度）
    wl(1) = inf;  % 避免除以零

    % 聚焦于条带相关波长范围（2-60°）
    mask = (wl >= 2) & (wl <= 60) & isfinite(pxx);
    wl_plot = wl(mask);
    pxx_plot = pxx(mask);

    if isempty(wl_plot)
        warning('[WARN] No valid frequency data for Lat=%d°', lat_val);
        continue;
    end

    % 子图1：功率谱（对数坐标）
    subplot(2, 3, 1); hold on;
    h = semilogx(wl_plot, 10*log10(pxx_plot + eps), 'Color', colors(i,:), ...
                 'LineWidth', 1.5, 'DisplayName', sprintf('Lat=%d°', lat_val));
    xlabel('Wavelength (°)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Power (dB)', 'FontSize', 11, 'FontWeight', 'bold');
    title('Residual Stripe Power Spectrum', 'FontSize', 12);
    grid on; xlim([2 60]); set(gca, 'XTick', [2 5 10 20 40 60]);

    % 子图2：累积能量分布
    subplot(2, 3, 2); hold on;
    energy_cumsum = cumsum(pxx_plot.^2) / sum(pxx_plot.^2);
    semilogx(wl_plot, energy_cumsum, 'Color', colors(i,:), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('Lat=%d°', lat_val));
    xlabel('Wavelength (°)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Cumulative Energy Ratio', 'FontSize', 11, 'FontWeight', 'bold');
    title('Residual Energy Distribution', 'FontSize', 12);
    grid on; xlim([2 60]); ylim([0 1]); set(gca, 'XTick', [2 5 10 20 40 60]);

    % 找出主导波长（能量最大的3个峰值）
    [pks, locs] = findpeaks(pxx_plot, 'SortStr', 'descend', 'NPeaks', 3);
    if length(locs) >= 3
        dominant_wl(i,:) = [wl_plot(locs(1)), wl_plot(locs(2)), wl_plot(locs(3))];
        fprintf('[DIAGNOSE] Lat=%+4d° - Top 3 residual wavelengths: %.1f°, %.1f°, %.1f°\n', ...
                lat_val, dominant_wl(i,1), dominant_wl(i,2), dominant_wl(i,3));
    else
        fprintf('[DIAGNOSE] Lat=%+4d° - Insufficient peaks detected\n', lat_val);
    end

    % 子图3：经向剖面（原始 vs 滤波 vs 残余）
    subplot(2, 3, 3); hold on;
    plot(lon, EWH(:, lat_idx, t_check), 'k-', 'LineWidth', 1, ...
         'DisplayName', sprintf('Before (Lat=%d°)', lat_val));
    plot(lon, EWH_filt(:, lat_idx, t_check), 'Color', colors(i,:), ...
         'LineWidth', 1.5, 'DisplayName', sprintf('After (Lat=%d°)', lat_val));
    plot(lon, x_resid, '--', 'Color', colors(i,:), 'LineWidth', 1, ...
         'DisplayName', sprintf('Residual (Lat=%d°)', lat_val));
    xlabel('Longitude (°)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('EWH (cm)', 'FontSize', 11, 'FontWeight', 'bold');
    title('Zonal Profiles', 'FontSize', 12);
    grid on; xlim([-180 180]);
end

% 添加图例
subplot(2, 3, 1); legend('Location', 'best', 'FontSize', 9);
subplot(2, 3, 2); legend('Location', 'best', 'FontSize', 9);
subplot(2, 3, 3); legend('Location', 'best', 'FontSize', 8);

% 子图4：空间分布（滤波前）
subplot(2, 3, 4);
plot_map(EWH(:,:,t_check), lon, lat, 1);
title('Before Filtering', 'FontSize', 12, 'FontWeight', 'bold');
clim([-15 15]); colorbar;

% 子图5：空间分布（滤波后）
subplot(2, 3, 5);
plot_map(EWH_filt(:,:,t_check), lon, lat, 1);
title('After HSAF', 'FontSize', 12, 'FontWeight', 'bold');
clim([-15 15]); colorbar;

% 子图6：残余场分布
subplot(2, 3, 6);
plot_map(residual(:,:,t_check), lon, lat, 1);
title('Residual Stripes', 'FontSize', 12, 'FontWeight', 'bold');
clim([-15 15]); colorbar;

% 保存诊断图
diag_fig = fullfile(output_dir, sprintf('residual_diagnosis_%s.png', timestamp));

diag_png = fullfile(output_dir, sprintf('residual_diagnosis_%s.png', timestamp));
print(fig_diag, diag_png, '-dpng', '-r300');

diag_figfile = fullfile(output_dir, sprintf('residual_diagnosis_%s.fig', timestamp));
savefig(fig_diag, diag_figfile);

diag_pdf = fullfile(output_dir, sprintf('residual_diagnosis_%s.pdf', timestamp));
exportgraphics(fig_diag, diag_pdf, 'ContentType', 'vector');

fprintf('[SAVE] Figures:\n  %s\n  %s\n  %s\n', diag_png, diag_figfile, diag_pdf);


fprintf('[DIAGNOSE] Diagnosis plot saved: %s\n', diag_fig);
close(fig_diag);

% ========================================================================
% 生成参数调整建议
% ========================================================================
fprintf('\n========================================\n');
fprintf('   Parameter Tuning Recommendation\n');
fprintf('========================================\n');

% 分析主导波长范围
all_wl = dominant_wl(dominant_wl > 0);
if ~isempty(all_wl)
    wl_min = min(all_wl);
    wl_max = max(all_wl);
    wl_mean = mean(all_wl);

    fprintf('Dominant wavelength range: %.1f - %.1f° (mean: %.1f°)\n', ...
            wl_min, wl_max, wl_mean);
    fprintf('Current filter bandwidth: %.1f - %.1f°\n', ...
            opt.band_abs_deg(1), opt.band_abs_deg(2));

    % 判断是否需要调整
    if wl_min < opt.band_abs_deg(1) - 1
        fprintf('\n⚠️  WARNING: Residual wavelengths BELOW current lower bound!\n');
        fprintf('   Recommended action: Decrease band_abs_deg(1) to %.1f°\n', ...
                floor(wl_min));
    end

    if wl_max > opt.band_abs_deg(2) + 1
        fprintf('\n⚠️  WARNING: Residual wavelengths ABOVE current upper bound!\n');
        fprintf('   Recommended action: Increase band_abs_deg(2) to %.1f°\n', ...
                ceil(wl_max));
    end

    if (wl_min >= opt.band_abs_deg(1)) && (wl_max <= opt.band_abs_deg(2))
        fprintf('\n✅ Residual wavelengths are WITHIN current bandwidth.\n');
        fprintf('   Possible issues:\n');
        fprintf('   1. k_max too small (current: %d) - try increasing to %d\n', ...
                opt.k_max, opt.k_max + 4);
        fprintf('   2. Window too short (current min: %d) - try increasing to %d\n', ...
                opt.win_min, opt.win_min + 30);
        fprintf('   3. Need second-pass iterative filtering\n');
    end
else
    fprintf('⚠️  No dominant wavelengths detected. Check data quality.\n');
end

% 计算残余能量比例
residual_rms = sqrt(mean(residual(:,:,t_check).^2, 'all', 'omitnan'));
original_rms = sqrt(mean(EWH(:,:,t_check).^2, 'all', 'omitnan'));
filtered_rms = sqrt(mean(EWH_filt(:,:,t_check).^2, 'all', 'omitnan'));

fprintf('\nRMS Statistics (cm):\n');
fprintf('  Original:  %.2f\n', original_rms);
fprintf('  Filtered:  %.2f\n', filtered_rms);
fprintf('  Residual:  %.2f (%.1f%% of original)\n', ...
        residual_rms, 100*residual_rms/original_rms);

wl_tbl = table(lat_test(:), dominant_wl(:,1), dominant_wl(:,2), dominant_wl(:,3), ...
    'VariableNames', {'Lat_deg','Top1_wl_deg','Top2_wl_deg','Top3_wl_deg'});
wl_csv = fullfile(output_dir, sprintf('dominant_wavelengths_%s.csv', timestamp));
writetable(wl_tbl, wl_csv);

rms_tbl = table(original_rms, filtered_rms, residual_rms, 100*residual_rms/original_rms, ...
    'VariableNames', {'Original_RMS_cm','Filtered_RMS_cm','Residual_RMS_cm','Residual_pct_of_original'});
rms_csv = fullfile(output_dir, sprintf('rms_stats_%s.csv', timestamp));
writetable(rms_tbl, rms_csv);

mat_file = fullfile(output_dir, sprintf('residual_diagnosis_%s.mat', timestamp));
save(mat_file, 'lat_test','lat_idx_valid','dominant_wl', ...
    'original_rms','filtered_rms','residual_rms', 't_check','timestamp','opt', ...
    'wl_tbl','rms_tbl', '-v7.3');

fprintf('[SAVE] Tables:\n  %s\n  %s\n  %s\n', wl_csv, rms_csv, mat_file);
fprintf('========================================\n\n');
fprintf('[LOG] Done. Log saved: %s\n', log_file);
diary off;

%% ============================================================
%% 11) 保存脚本快照
%% ============================================================

% script_snapshot = fullfile(output_dir, 'HSAF_V6_snapshot.m');
% 
% src = [mfilename('fullpath') 'HSAF_V6.m'];   % 关键：补后缀
% try
%     copyfile(src, script_snapshot);
%     fprintf('[SNAP] Script snapshot saved: %s\n', script_snapshot);
% catch ME
%     warning('[WARN] Failed to save script snapshot: %s (%s)', ME.message, ME.identifier);
% end
% 
src = mfilename('fullpath');
if ~endsWith(src, '.m'); src = [src '.m']; end
dst = fullfile(output_dir, 'HSAF_V6_snapshot.m');

if ~exist(output_dir,'dir'); mkdir(output_dir); end

txt = fileread(src);              % 读整个脚本
fid = fopen(dst,'w');             % 覆盖写入
assert(fid>0, 'Cannot open dst for writing: %s', dst);

fwrite(fid, txt, 'char');
fclose(fid);

fprintf('[SNAP] Script snapshot saved: %s\n', dst);

% 统计总滤波耗时（迭代求和）
if exist('iteration_stats','var') && ~isempty(iteration_stats)
    t_filter_total = sum([iteration_stats.time_sec]);
else
    t_filter_total = 0;
end

%% ============================================================
%% 12) 生成运行报告
%% ============================================================
fprintf('\n========================================\n');
fprintf('   HSAF V6 Run Summary\n');
fprintf('========================================\n');
fprintf('End time:        %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('Total runtime:   %.1f minutes\n', (t_calib + t_filter_total)/60);
fprintf('Calibration:     %.1f seconds\n', t_calib);
fprintf('Filtering:       %.1f seconds (%.2f sec/month)\n', ...
        t_filter_total, t_filter_total/Nt);
fprintf('Output dir:      %s\n', output_dir);
fprintf('Data file:       CSR_EWH_HSAF_filtered.mat\n');
fprintf('Log file:        hsaf_v6.log\n');
fprintf('========================================\n');

diary off;

%% ========================================================================
%% 第九部分：最终报告
%% ========================================================================

fprintf('\n╔════════════════════════════════════════╗\n');
fprintf('║  HSAF V6.0 Run Summary                 ║\n');
fprintf('╚════════════════════════════════════════╝\n');
fprintf('End time:        %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('Total runtime:   %.1f minutes\n', ...
        (t_calib + sum([iteration_stats.time_sec]))/60);
fprintf('\nIteration breakdown:\n');
for iter = 1:length(iteration_stats)
    fprintf('  Iter%d: %.1fs | Residual: %.1f%%\n', ...
            iter, iteration_stats(iter).time_sec, ...
            iteration_stats(iter).residual_pct);
end
fprintf('\nOutput directory: %s\n', output_dir);
fprintf('Main log:         hsaf_v6.log\n');
fprintf('Diagnosis log:    final_diagnosis_%s.log\n', timestamp);
fprintf('════════════════════════════════════════════\n');

diary off;
%% QUAILTY AND FORMATTING OUTPUT

%% ========================================================================
%% HSAF_Quality_Control.m - HSAF V6配套质量评估模块
%% ========================================================================
%
% 功能：对HSAF_V6滤波结果进行全面质量评估
% 
% 评估指标：
%   1. CC (相关系数)
%   2. SNR (信噪比)
%   3. RMSE (均方根误差)
%   4. PSNR (峰值信噪比)
%   5. MAE (平均绝对误差)
%   6. NSC (Nash-Sutcliffe系数)
%   7. SRMSE (空间均方根误差)
%   8. Basin RMSE (流域尺度误差)
%   9. 球谐系数误差谱
%   10. 大地水准面阶误差
%
% 使用方法：
%   1. 先运行 HSAF_V6.m 生成滤波结果
%   2. 修改下方 output_dir 为HSAF_V6的输出路径
%   3. 运行本脚本
%
% 作者：HSAF Team
% 日期：2025-01-06
% 版本：V6.0
%
%% ========================================================================

clc; clear; warning('off', 'verbose');

fprintf('\n╔════════════════════════════════════════╗\n');
fprintf('║  HSAF Quality Control Module          ║\n');
fprintf('║  For HSAF V6 Filtered Data            ║\n');
fprintf('╚════════════════════════════════════════╝\n');
fprintf('Start time: %s\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

%% ========================================================================
%% 第一部分：路径配置与数据加载
%% ========================================================================

%% 1.1) 路径配置
base_path = '/home/um202370130/HSAF';
data_path = fullfile(base_path, 'Data');
script_path = fullfile(base_path, 'Software');

% 添加必要路径
addpath(fullfile(script_path, 'Tool_Functions', 'gmt'));
addpath(fullfile(script_path, 'Tool_Functions'));
addpath(fullfile(script_path, 'Tool_Functions', 'm_map'));
addpath(fullfile(script_path, 'Tool_Functions', 'GRACE-filter-master', 'src', 'matlab'));

%% 1.2) 用户输入：指定HSAF_V6输出路径
fprintf('[INPUT] Please specify HSAF V6 output directory:\n');
fprintf('        Example: /home/.../HSAF_Adaptive/20251223_153608\n');
fprintf('        Or press Enter to use latest output\n');

% user_input = input('Output directory: ', 's');
% 
% if isempty(user_input)
%     % 自动查找最新输出
%     output_base = fullfile(data_path, 'EWH_Output', 'CSR_EWH', 'HSAF_Adaptive');
%     dirs = dir(output_base);
%     dirs = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
%     [~, idx] = max([dirs.datenum]);
%     hsaf_output_dir = fullfile(output_base, dirs(idx).name);
%     fprintf('[AUTO] Using latest output: %s\n', hsaf_output_dir);
% else
%     hsaf_output_dir = user_input;
% end

output_base = fullfile(data_path, 'EWH_Output', 'CSR_EWH', 'HSAF_Adaptive');
dirs = dir(output_base);
dirs = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
[~, idx] = max([dirs.datenum]);
hsaf_output_dir = fullfile(output_base, dirs(idx).name);
fprintf('[AUTO] Using latest output: %s\n', hsaf_output_dir);


% 验证路径存在
if ~exist(hsaf_output_dir, 'dir')
    error('[ERROR] Output directory not found: %s', hsaf_output_dir);
end

% 创建统计输出子文件夹
stats_output_dir = fullfile(hsaf_output_dir, 'Statistics');
if ~exist(stats_output_dir, 'dir'), mkdir(stats_output_dir); end
fprintf('[DIR] Statistics output: %s\n\n', stats_output_dir);

%% 1.3) 加载HSAF滤波结果
fprintf('[LOAD] Loading HSAF filtered data...\n');
try
    load(fullfile(hsaf_output_dir, 'CSR_EWH_HSAF_filtered.mat'), ...
         'CSR_EWH', 'param', 'opt', 'dateTime', 'dates', 'csr_lon', 'csr_lat');
    fprintf('[LOAD] ✓ HSAF filtered data loaded\n');
catch ME
    error('[ERROR] Failed to load HSAF output: %s', ME.message);
end

%% 1.4) 加载辅助数据
fprintf('[LOAD] Loading auxiliary data...\n');
try
    % 辅助数据路径
    aux_path = fullfile(data_path, 'Auxiliary_Data');
    
    % 土地掩膜
    load(fullfile(aux_path, 'land_mask.mat'), 'land_mask_r');
    
    % 300km掩膜（用于SNR计算）
    [lon_msk, lat_msk, bool_msk] = textread(fullfile(aux_path, 'msk_300.xyz'), ...
                                             '%f %f %u', 'headerlines', 0);
    msk = reshape(bool_msk, [360, 180]);
    land_mask_300_r(1:180, :) = msk(181:360, :);
    land_mask_300_r(181:360, :) = msk(1:180, :);
    land_mask_300_r = fliplr(land_mask_300_r);
    land_mask_r_300_025 = ones(360, 180);
    land_mask_r_300_025(land_mask_300_r == 0) = NaN;
    
    % 流域数据
    load(fullfile(aux_path, 'rivers1.mat'), 'rivers_new');
    
    % 时间数据
    load(fullfile(aux_path, 'time.mat'), 'years');
    
    fprintf('[LOAD] ✓ Auxiliary data loaded\n');
catch ME
    error('[ERROR] Failed to load auxiliary data: %s', ME.message);
end

%% 1.5) 数据准备
% 提取关键变量
lon = csr_lon(:);
lat = csr_lat(:);
[Nlon, Nlat, ~] = size(CSR_EWH.None);

% 确定时间范围（与HSAF_V6一致）
t0 = 19; t1 = 150;  % 2004.01-2016.01
start_idx = t0;
end_idx = t1;
time_size = end_idx - start_idx + 1;
dates = dateTime(start_idx:end_idx);

fprintf('[DATA] Processing period: %s to %s\n', ...
        datestr(dates(1), 'yyyy-mm'), datestr(dates(end), 'yyyy-mm'));
fprintf('       Grid: %d lon × %d lat × %d months\n\n', Nlon, Nlat, time_size);

%% 1.6) 构建EWH结构体（统一接口）
fprintf('[PREP] Constructing unified EWH structure...\n');

EWH = struct();
EWH.None = CSR_EWH.None(:, :, start_idx:end_idx);           % 原始数据
EWH.Hankel = CSR_EWH.HankelLatAdaptive(:, :, start_idx:end_idx);  % HSAF滤波结果
EWH.Mascon = CSR_EWH.Mascon(:, :, start_idx:end_idx);      % 参考数据

% 如果有其他滤波方法，也加载进来
if isfield(CSR_EWH, 'Gaussian')
    EWH.Gaussian = CSR_EWH.Gaussian(:, :, start_idx:end_idx);
end
if isfield(CSR_EWH, 'DDK4')
    EWH.DDK4 = CSR_EWH.DDK4(:, :, start_idx:end_idx);
else
    % 7. 计算 DDK4 滤波数据
    EWH.DDK4 = DDKs_Filter(CSR_EWH.None(:,:,start_idx:end_idx), 'DDK4', 1);

end
if isfield(CSR_EWH, 'Fan_Decorrelation')
    EWH.Fan_Decorrelation = CSR_EWH.Fan_Decorrelation(:, :, start_idx:end_idx);
end

filter_methods = fieldnames(EWH);
num_methods = length(filter_methods);
fprintf('[PREP] ✓ %d filtering methods prepared:\n', num_methods);
for i = 1:num_methods
    fprintf('       - %s\n', filter_methods{i});
end
fprintf('\n');

%% ========================================================================
%% 第二部分：统计指标计算
%% ========================================================================

fprintf('[STATS] ═══════════════════════════════════════\n');
fprintf('        Computing Statistical Metrics\n');
fprintf('        ═══════════════════════════════════════\n\n');

%% 2.1) 初始化存储结构
CC_results = struct();
SNR_results = struct();
RMSE_results = struct();
PSNR_results = struct();
MAE_results = struct();
NSC_results = struct();
Basin_RMSE_results = struct();
SRMSE_accumulator = struct();

for i = 1:num_methods
    method = filter_methods{i};
    CC_results.(method) = zeros(time_size, 1);
    SNR_results.(method) = zeros(time_size, 1);
    RMSE_results.(method) = zeros(time_size, 1);
    PSNR_results.(method) = zeros(time_size, 1);
    MAE_results.(method) = zeros(time_size, 1);
    NSC_results.(method) = zeros(time_size, 1);
    Basin_RMSE_results.(method) = zeros(112, 1);  % 假设112个流域
    SRMSE_accumulator.(method) = zeros(Nlon, Nlat);
end

%% 2.2) 逐月计算全局统计指标
fprintf('[STATS] Computing global metrics (逐月计算)...\n');
tic;

for t = 1:time_size
    Ft = EWH.Mascon(:, :, t);  % 参考数据
    mean_Ft = mean(Ft(:), 'omitnan');
    
    for i = 1:num_methods
        method = filter_methods{i};
        Fo = EWH.(method)(:, :, t);
        mean_Fo = mean(Fo(:), 'omitnan');
        
        % 1) 相关系数 (CC)
        numerator = sum((Fo(:) - mean_Fo) .* (Ft(:) - mean_Ft), 'omitnan');
        denominator = sqrt(sum((Fo(:) - mean_Fo).^2, 'omitnan')) * ...
                      sqrt(sum((Ft(:) - mean_Ft).^2, 'omitnan'));
        CC_results.(method)(t) = numerator / denominator;
        
        % 2) Nash-Sutcliffe系数 (NSC)
        numerator_NSC = sum((Fo(:) - Ft(:)).^2, 'omitnan');
        denominator_NSC = sum((Ft(:) - mean_Ft).^2, 'omitnan');
        NSC_results.(method)(t) = 1 - numerator_NSC / denominator_NSC;
        
        % 3) 均方根误差 (RMSE)
        RMSE_results.(method)(t) = sqrt(mean((Fo(:) - Ft(:)).^2, 'omitnan'));
        
        % 4) 平均绝对误差 (MAE)
        MAE_results.(method)(t) = mean(abs(Fo(:) - Ft(:)), 'omitnan');
        
        % 5) 峰值信噪比 (PSNR)
        max_Ft = max(Ft(:), [], 'omitnan')^2;
        PSNR_results.(method)(t) = 10 * log10(max_Ft / (RMSE_results.(method)(t)^2 + eps));
        
        % 6) 累积SRMSE误差平方
        SRMSE_accumulator.(method) = SRMSE_accumulator.(method) + (Fo - Ft).^2;
        
        % 7) 信噪比 (SNR) - 陆地RMS vs 海洋RMS
        land_EWH = Fo(land_mask_r_300_025 == 1);
        ocean_EWH = Fo(isnan(land_mask_r_300_025));
        RMS_land = sqrt(mean(land_EWH.^2, 'omitnan'));
        RMS_ocean = sqrt(mean(ocean_EWH.^2, 'omitnan'));
        SNR_results.(method)(t) = 10 * log10(RMS_land / (RMS_ocean + eps));
    end
    
    % 进度报告
    if mod(t, 20) == 0 || t == 1 || t == time_size
        fprintf('[STATS] Processed %3d/%d months (%.1f%%)\n', ...
                t, time_size, 100*t/time_size);
    end
end

% 计算SRMSE（时间平均）
for i = 1:num_methods
    method = filter_methods{i};
    SRMSE_results.(method) = sqrt(SRMSE_accumulator.(method) / time_size);
end

t_global = toc;
fprintf('[STATS] ✓ Global metrics completed in %.1f sec\n\n', t_global);

%% 2.3) 球谐系数转换
fprintf('[STATS] Converting to spherical harmonic coefficients...\n');
tic;

SC_results = struct();
for method = filter_methods'
    SC_results.(method{1}) = zeros(61, 121, time_size);
end

for t = 1:time_size
    for i = 1:num_methods
        method = filter_methods{i};
        % EWH (cm) -> m -> SHC
        cs = gmt_grid2cs(EWH.(method)(:, :, t)' / 100, 60);
        gc = gmt_mc2gc(cs);
        sc = gmt_cs2sc(gc);
        SC_results.(method)(:, :, t) = sc;
    end
    
    if mod(t, 20) == 0
        fprintf('[STATS] SHC conversion: %d/%d months\n', t, time_size);
    end
end

% 处理NaN值
for i = 1:num_methods
    method = filter_methods{i};
    SC_results.(method)(SC_results.(method) == 0) = NaN;
end

t_shc = toc;
fprintf('[STATS] ✓ SHC conversion completed in %.1f sec\n\n', t_shc);

%% 2.4) 大地水准面阶误差
fprintf('[STATS] Computing geoid degree errors...\n');
tic;

lmax = 60;
GeoidError = struct();
degree = (0:lmax)';

% 参考球谐系数（Mascon）
C1_ref = SC_results.Mascon(:, lmax+1:2*lmax+1, 1);
S1_ref = [zeros(lmax+1, 1), SC_results.Mascon(:, 1:lmax, 1)];
C1_ref(isnan(C1_ref)) = 0;
S1_ref(isnan(S1_ref)) = 0;

for t = 1:time_size
    C1 = SC_results.Mascon(:, lmax+1:2*lmax+1, t);
    S1 = [zeros(lmax+1, 1), SC_results.Mascon(:, 1:lmax, t)];
    C1(isnan(C1)) = 0;
    S1(isnan(S1)) = 0;
    
    for i = 1:num_methods
        method = filter_methods{i};
        C_current = SC_results.(method)(:, lmax+1:2*lmax+1, t);
        S_current = [zeros(lmax+1, 1), SC_results.(method)(:, 1:lmax, t)];
        C_current(isnan(C_current)) = 0;
        S_current(isnan(S_current)) = 0;
        
        [GeoidError.(method).Degree(:, :, t), GeoidError.(method).Cumulative(:, :, t)] = ...
            Function_GeoidErrorRealEGM(C_current, S_current, C1, S1, lmax, lmax);
    end
end

% 计算平均误差
mean_errors = struct();
for i = 1:num_methods
    method = filter_methods{i};
    mean_errors.(method) = mean(GeoidError.(method).Degree(:, 2, :), 3);
end

t_geoid = toc;
fprintf('[STATS] ✓ Geoid errors completed in %.1f sec\n\n', t_geoid);

%% 2.5) 全球相关系数分布
fprintf('[STATS] Computing global correlation coefficient map...\n');
tic;

global_cc = struct();

for i = 1:num_methods
    method = filter_methods{i};
    global_cc.(method) = nan(Nlon, Nlat);
    
    for x = 1:Nlon
        for y = 1:Nlat
            ts_method = squeeze(EWH.(method)(x, y, :));
            ts_mascon = squeeze(EWH.Mascon(x, y, :));
            
            if sum(~isnan(ts_method)) > 10 && sum(~isnan(ts_mascon)) > 10
                [global_cc.(method)(x, y), ~] = corr(ts_method, ts_mascon, ...
                                                      'rows', 'complete');
            end
        end
    end
end

t_cc_map = toc;
fprintf('[STATS] ✓ Global CC map completed in %.1f sec\n\n', t_cc_map);

%% 2.6) 流域尺度分析
fprintf('[STATS] Computing basin-scale metrics...\n');
tic;

basin = struct('Lon', [], 'Lat', [], 'mask', []);
num_basins = min(112, length(rivers_new));

for i = 1:num_basins
    [basin(i).Lon, basin(i).Lat, basin(i).mask] = mkmask(rivers_new(i), 1);
    basin(i).name = rivers_new(i).DRAINAGE;
    basin(i).AREA = rivers_new(i).AREA;
    
    for j = 1:num_methods
        method = filter_methods{j};
        
        % 提取流域数据
        basin(i).(['ewh_', method]) = EWH.(method)( ...
            find(csr_lon == basin(i).Lon(1)):find(csr_lon == basin(i).Lon(end)), ...
            find(csr_lat == basin(i).Lat(1)):find(csr_lat == basin(i).Lat(end)), ...
            :) .* basin(i).mask';
        
        % 调和分析
        [basin(i).(['Annual_Amplitude_', method]), ~, ~, ~, ...
         basin(i).(['Semi_Annual_Amplitude_', method]), ~, ~, ~, ...
         basin(i).(['Trend_', method]), ~, ~, ...
         basin(i).(['residual_', method]), ~] = ...
            gmt_harmonic(years(start_idx:end_idx), [], basin(i).(['ewh_', method]));
        
        % 流域时间序列
        [~, ~, ~, ~, basin(i).(['time_series_', method]), basin(i).(['b_', method])] = ...
            Basin_Analysis(basin(i).(['ewh_', method]), years(start_idx:end_idx)', ...
                          basin(i).Lon, basin(i).Lat);
    end
    
    % 计算流域RMSE
    for j = 1:num_methods
        method = filter_methods{j};
        if strcmp(method, 'Mascon'), continue; end
        
        Basin_RMSE_results.(method)(i) = rmse( ...
            basin(i).(['time_series_', method]), ...
            basin(i).time_series_Mascon);
    end
    
    if mod(i, 20) == 0
        fprintf('[STATS] Processed %d/%d basins\n', i, num_basins);
    end
end

t_basin = toc;
fprintf('[STATS] ✓ Basin metrics completed in %.1f sec\n\n', t_basin);

%% 2.7) 全球时间序列调和分析
fprintf('[STATS] Computing global time-series harmonic analysis...\n');
tic;

global_ewh = struct();

for i = 1:num_methods
    method = filter_methods{i};
    
    [global_ewh.(['Amplitude1_', method]), ~, global_ewh.(['Phase1_', method]), ~, ...
     global_ewh.(['Amplitude2_', method]), ~, global_ewh.(['Phase2_', method]), ~, ...
     global_ewh.(['Trend_', method]), ~, global_ewh.(['Trend_line_', method]), ...
     global_ewh.(['Resid_', method]), ~] = ...
        gmt_harmonic(years(start_idx:end_idx), [], EWH.(method));
end

t_harmonic = toc;
fprintf('[STATS] ✓ Harmonic analysis completed in %.1f sec\n\n', t_harmonic);

%% ========================================================================
%% 第三部分：结果可视化
%% ========================================================================

fprintf('[PLOT] ═══════════════════════════════════════\n');
fprintf('       Generating Visualization Plots\n');
fprintf('       ═══════════════════════════════════════\n\n');

%% 3.1) 统计指标时间序列图
fprintf('[PLOT] Plotting statistical time series...\n');

fig = figure('Position', [100 100 1600 900], 'Visible', 'off');
tiledlayout(2, 3, "TileSpacing", "compact", "Padding", "compact");

colors = lines(num_methods);
markers = {'o', 's', '^', 'd', 'v', '>', '<', 'p', 'h', '*'};

% 子图1: 相关系数
nexttile;
hold on;
for i = 1:num_methods
    method = filter_methods{i};
    if strcmp(method, 'None'), continue; end
    plot(dates, CC_results.(method), ...
         'Color', colors(i, :), 'Marker', markers{mod(i-1, 10)+1}, ...
         'MarkerIndices', 1:5:time_size, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('%s (μ=%.4f)', method, mean(CC_results.(method))));
end
ylabel('Correlation Coefficient', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9, 'Interpreter', 'none');
grid on; xtickformat('yyyy-MM');

% 子图2: SNR
nexttile;
hold on;
for i = 1:num_methods
    method = filter_methods{i};
    if strcmp(method, 'None'), continue; end
    plot(dates, SNR_results.(method), ...
         'Color', colors(i, :), 'Marker', markers{mod(i-1, 10)+1}, ...
         'MarkerIndices', 1:5:time_size, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('%s (μ=%.2f dB)', method, mean(SNR_results.(method))));
end
ylabel('SNR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9, 'Interpreter', 'none');
grid on; xtickformat('yyyy-MM');

% 子图3: RMSE
nexttile;
hold on;
for i = 1:num_methods
    method = filter_methods{i};
    if strcmp(method, 'None') || strcmp(method, 'Mascon'), continue; end
    plot(dates, RMSE_results.(method), ...
         'Color', colors(i, :), 'Marker', markers{mod(i-1, 10)+1}, ...
         'MarkerIndices', 1:5:time_size, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('%s (μ=%.2f cm)', method, mean(RMSE_results.(method))));
end
ylabel('RMSE (cm)', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9, 'Interpreter', 'none');
grid on; xtickformat('yyyy-MM');

% 子图4: PSNR
nexttile;
hold on;
for i = 1:num_methods
    method = filter_methods{i};
    if strcmp(method, 'None') || strcmp(method, 'Mascon'), continue; end
    plot(dates, PSNR_results.(method), ...
         'Color', colors(i, :), 'Marker', markers{mod(i-1, 10)+1}, ...
         'MarkerIndices', 1:5:time_size, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('%s (μ=%.2f dB)', method, mean(PSNR_results.(method))));
end
ylabel('PSNR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9, 'Interpreter', 'none');
grid on; xtickformat('yyyy-MM');

% 子图5: MAE
nexttile;
hold on;
for i = 1:num_methods
    method = filter_methods{i};
    if strcmp(method, 'None') || strcmp(method, 'Mascon'), continue; end
    plot(dates, MAE_results.(method), ...
         'Color', colors(i, :), 'Marker', markers{mod(i-1, 10)+1}, ...
         'MarkerIndices', 1:5:time_size, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('%s (μ=%.2f cm)', method, mean(MAE_results.(method))));
end
ylabel('MAE (cm)', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9, 'Interpreter', 'none');
grid on; xtickformat('yyyy-MM');

% 子图6: NSC
nexttile;
hold on;
for i = 1:num_methods
    method = filter_methods{i};
    if strcmp(method, 'None') || strcmp(method, 'Mascon'), continue; end
    plot(dates, NSC_results.(method), ...
         'Color', colors(i, :), 'Marker', markers{mod(i-1, 10)+1}, ...
         'MarkerIndices', 1:5:time_size, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('%s (μ=%.4f)', method, mean(NSC_results.(method))));
end
ylabel('NSC', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9, 'Interpreter', 'none');
grid on; xtickformat('yyyy-MM');

sgtitle('Global Statistical Metrics Time Series', 'FontSize', 16, 'FontWeight', 'bold');
print(fig, fullfile(stats_output_dir, 'Global_Statistics_TimeSeries.png'), '-dpng', '-r300');
fprintf('[PLOT] ✓ Statistical time series saved\n');
close(fig);

%% 3.2) 大地水准面阶误差图
fprintf('[PLOT] Plotting geoid degree errors...\n');

fig = figure('Position', [100 100 1200 700], 'Visible', 'off');
hold on;

for i = 1:num_methods
    method = filter_methods{i};
    if strcmp(method, 'None'), continue; end
    semilogy(degree, mean_errors.(method), ...
             'Color', colors(i, :), 'LineWidth', 2, ...
             'DisplayName', method);
end

xlabel('Degree', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Geoid Degree Error (m)', 'FontSize', 14, 'FontWeight', 'bold');
title('Geoid Degree Error (Averaged Over All Months)', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 12, 'Interpreter', 'none');
grid on; grid minor;
set(gca, 'YScale', 'log', 'FontSize', 12);

print(fig, fullfile(stats_output_dir, 'Geoid_Degree_Error.png'), '-dpng', '-r300');
fprintf('[PLOT] ✓ Geoid degree error plot saved\n');
close(fig);

%% 3.3) 流域时间序列图（前12个主要流域）
fprintf('[PLOT] Plotting basin time series...\n');

fig = figure('Position', [100 100 1600 1000], 'Visible', 'off');
tiledlayout(3, 4, "TileSpacing", "compact", "Padding", "compact");

for i = 1:min(12, num_basins)
    nexttile;
    hold on;
    for j = 1:num_methods
        method = filter_methods{j};
        if strcmp(method, 'None'), continue; end
        plot(dates, basin(i).(['time_series_', method]), ...
             'Color', colors(j, :), 'LineWidth', 1.5, ...
             'DisplayName', method);
    end
    title(basin(i).name, 'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'none');
    grid on; xtickformat('yyyy');
    if i == 1
        legend('Location', 'best', 'FontSize', 8, 'Interpreter', 'none');
    end
end

sgtitle('Basin-Scale EWH Time Series (Top 12 Basins)', 'FontSize', 16, 'FontWeight', 'bold');
print(fig, fullfile(stats_output_dir, 'Basin_TimeSeries_Top12.png'), '-dpng', '-r300');
fprintf('[PLOT] ✓ Basin time series plot saved\n');
close(fig);

%% ========================================================================
%% 第四部分：保存统计结果
%% ========================================================================

fprintf('\n[SAVE] ═══════════════════════════════════════\n');
fprintf('       Saving Statistical Results\n');
fprintf('       ═══════════════════════════════════════\n\n');

%% 4.1) 保存MAT文件
fprintf('[SAVE] Saving .mat file...\n');
save_file = fullfile(stats_output_dir, 'HSAF_Quality_Statistics.mat');
save(save_file, 'CC_results', 'SNR_results', 'RMSE_results', 'PSNR_results', ...
     'MAE_results', 'NSC_results', 'SRMSE_results', 'Basin_RMSE_results', ...
     'SC_results', 'GeoidError', 'mean_errors', 'global_cc', 'global_ewh', ...
     'basin', 'filter_methods', 'dates', 'param', 'opt', '-v7.3');
fprintf('[SAVE] ✓ MAT file saved (%.1f MB): %s\n', ...
        dir(save_file).bytes / 1024^2, save_file);

%% 4.2) 保存CSV汇总表
fprintf('[SAVE] Generating CSV summary tables...\n');

% 全局统计汇总
summary_tbl = table();
summary_tbl.Method = filter_methods;

for i = 1:num_methods
    method = filter_methods{i};
    if strcmp(method, 'Mascon') || strcmp(method, 'None')
        summary_tbl.CC_mean(i) = NaN;
        summary_tbl.SNR_mean_dB(i) = NaN;
        summary_tbl.RMSE_mean_cm(i) = NaN;
        summary_tbl.PSNR_mean_dB(i) = NaN;
        summary_tbl.MAE_mean_cm(i) = NaN;
        summary_tbl.NSC_mean(i) = NaN;
    else
        summary_tbl.CC_mean(i) = mean(CC_results.(method));
        summary_tbl.SNR_mean_dB(i) = mean(SNR_results.(method));
        summary_tbl.RMSE_mean_cm(i) = mean(RMSE_results.(method));
        summary_tbl.PSNR_mean_dB(i) = mean(PSNR_results.(method));
        summary_tbl.MAE_mean_cm(i) = mean(MAE_results.(method));
        summary_tbl.NSC_mean(i) = mean(NSC_results.(method));
    end
end

csv_summary = fullfile(stats_output_dir, 'Global_Statistics_Summary.csv');
writetable(summary_tbl, csv_summary);
fprintf('[SAVE] ✓ Global summary: %s\n', csv_summary);

% HSAF性能汇总（重点关注Hankel方法）
if ismember('Hankel', filter_methods)
    hsaf_perf = table();
    hsaf_perf.Metric = {'CC'; 'SNR_dB'; 'RMSE_cm'; 'PSNR_dB'; 'MAE_cm'; 'NSC'};
    hsaf_perf.Mean = [
        mean(CC_results.Hankel);
        mean(SNR_results.Hankel);
        mean(RMSE_results.Hankel);
        mean(PSNR_results.Hankel);
        mean(MAE_results.Hankel);
        mean(NSC_results.Hankel)
    ];
    hsaf_perf.Std = [
        std(CC_results.Hankel);
        std(SNR_results.Hankel);
        std(RMSE_results.Hankel);
        std(PSNR_results.Hankel);
        std(MAE_results.Hankel);
        std(NSC_results.Hankel)
    ];
    hsaf_perf.Min = [
        min(CC_results.Hankel);
        min(SNR_results.Hankel);
        min(RMSE_results.Hankel);
        min(PSNR_results.Hankel);
        min(MAE_results.Hankel);
        min(NSC_results.Hankel)
    ];
    hsaf_perf.Max = [
        max(CC_results.Hankel);
        max(SNR_results.Hankel);
        max(RMSE_results.Hankel);
        max(PSNR_results.Hankel);
        max(MAE_results.Hankel);
        max(NSC_results.Hankel)
    ];
    
    csv_hsaf = fullfile(stats_output_dir, 'HSAF_Performance_Summary.csv');
    writetable(hsaf_perf, csv_hsaf);
    fprintf('[SAVE] ✓ HSAF performance: %s\n', csv_hsaf);
end

%% 4.3) 生成HTML报告
fprintf('[SAVE] Generating HTML quality report...\n');

html_file = fullfile(stats_output_dir, 'Quality_Report.html');
fid = fopen(html_file, 'w');

fprintf(fid, '<html><head><title>HSAF Quality Control Report</title>\n');
fprintf(fid, '<style>body{font-family:Arial;margin:40px;} table{border-collapse:collapse;width:100%%;} th,td{border:1px solid #ddd;padding:8px;text-align:center;} th{background-color:#4CAF50;color:white;}</style>\n');
fprintf(fid, '</head><body>\n');
fprintf(fid, '<h1>HSAF V6 Quality Control Report</h1>\n');
fprintf(fid, '<p><b>Generated:</b> %s</p>\n', datestr(now));
fprintf(fid, '<p><b>HSAF Output Directory:</b> %s</p>\n', hsaf_output_dir);
fprintf(fid, '<p><b>Time Period:</b> %s to %s (%d months)</p>\n', ...
        datestr(dates(1), 'yyyy-mm'), datestr(dates(end), 'yyyy-mm'), time_size);

fprintf(fid, '<h2>1. Global Statistics Summary</h2>\n');
fprintf(fid, '<table><tr><th>Method</th><th>CC</th><th>SNR (dB)</th><th>RMSE (cm)</th><th>PSNR (dB)</th><th>MAE (cm)</th><th>NSC</th></tr>\n');
for i = 1:height(summary_tbl)
    fprintf(fid, '<tr><td>%s</td><td>%.4f</td><td>%.2f</td><td>%.2f</td><td>%.2f</td><td>%.2f</td><td>%.4f</td></tr>\n', ...
            summary_tbl.Method{i}, summary_tbl.CC_mean(i), summary_tbl.SNR_mean_dB(i), ...
            summary_tbl.RMSE_mean_cm(i), summary_tbl.PSNR_mean_dB(i), ...
            summary_tbl.MAE_mean_cm(i), summary_tbl.NSC_mean(i));
end
fprintf(fid, '</table>\n');

if ismember('Hankel', filter_methods)
    fprintf(fid, '<h2>2. HSAF (Hankel) Performance Details</h2>\n');
    fprintf(fid, '<table><tr><th>Metric</th><th>Mean</th><th>Std</th><th>Min</th><th>Max</th></tr>\n');
    for i = 1:height(hsaf_perf)
        fprintf(fid, '<tr><td>%s</td><td>%.4f</td><td>%.4f</td><td>%.4f</td><td>%.4f</td></tr>\n', ...
                hsaf_perf.Metric{i}, hsaf_perf.Mean(i), hsaf_perf.Std(i), ...
                hsaf_perf.Min(i), hsaf_perf.Max(i));
    end
    fprintf(fid, '</table>\n');
end

fprintf(fid, '<h2>3. Visualization</h2>\n');
fprintf(fid, '<img src="Global_Statistics_TimeSeries.png" width="100%%"><br>\n');
fprintf(fid, '<img src="Geoid_Degree_Error.png" width="100%%"><br>\n');
fprintf(fid, '<img src="Basin_TimeSeries_Top12.png" width="100%%"><br>\n');

fprintf(fid, '</body></html>\n');
fclose(fid);
fprintf('[SAVE] ✓ HTML report: %s\n', html_file);

%% ========================================================================
%% 第五部分：最终报告
%% ========================================================================

fprintf('\n╔════════════════════════════════════════╗\n');
fprintf('║  Quality Control Summary              ║\n');
fprintf('╚════════════════════════════════════════╝\n\n');

fprintf('Processing Time:\n');
fprintf('  - Global metrics:  %.1f sec\n', t_global);
fprintf('  - SHC conversion:  %.1f sec\n', t_shc);
fprintf('  - Geoid errors:    %.1f sec\n', t_geoid);
fprintf('  - CC map:          %.1f sec\n', t_cc_map);
fprintf('  - Basin analysis:  %.1f sec\n', t_basin);
fprintf('  - Harmonic:        %.1f sec\n', t_harmonic);
fprintf('  Total:             %.1f sec (%.2f min)\n\n', ...
        t_global+t_shc+t_geoid+t_cc_map+t_basin+t_harmonic, ...
        (t_global+t_shc+t_geoid+t_cc_map+t_basin+t_harmonic)/60);

if ismember('Hankel', filter_methods)
    fprintf('HSAF Performance Highlights:\n');
    fprintf('  - Correlation Coefficient:  %.4f ± %.4f\n', ...
            mean(CC_results.Hankel), std(CC_results.Hankel));
    fprintf('  - Signal-to-Noise Ratio:    %.2f ± %.2f dB\n', ...
            mean(SNR_results.Hankel), std(SNR_results.Hankel));
    fprintf('  - Root Mean Square Error:   %.2f ± %.2f cm\n', ...
            mean(RMSE_results.Hankel), std(RMSE_results.Hankel));
    fprintf('  - Nash-Sutcliffe Coeff:     %.4f ± %.4f\n\n', ...
            mean(NSC_results.Hankel), std(NSC_results.Hankel));
end

fprintf('Output Files:\n');
fprintf('  - Statistics:     %s\n', stats_output_dir);
fprintf('  - MAT file:       HSAF_Quality_Statistics.mat\n');
fprintf('  - Summary CSV:    Global_Statistics_Summary.csv\n');
fprintf('  - HTML report:    Quality_Report.html\n');
fprintf('  - Plots:          3 PNG figures\n\n');

fprintf('════════════════════════════════════════════\n');
fprintf('Quality control completed successfully!\n');
fprintf('End time: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('════════════════════════════════════════════\n\n');


%% ========================================================================
%% HSAF_Visualization.m - HSAF V6配套可视化模块
%% ========================================================================
%
% 功能：对质量评估结果生成高质量可视化图表
% 
% 图表类型：
%   1. 滤波前后对比图（多月多方法）
%   2. 球谐误差谱对比
%   3. 全球趋势分布图
%   4. 流域EWH分布图
%   5. 流域趋势对比图
%   6. 空间相关系数分布图
%   7. 全局统计指标对比柱状图
%
% 使用方法：
%   先运行 HSAF_Quality_Control.m 生成统计结果，再运行本脚本
%
% 作者：HSAF Team
% 日期：2025-01-06
% 版本：V6.0
%
%% ========================================================================

clc; clear; warning('off', 'verbose');

fprintf('\n╔════════════════════════════════════════╗\n');
fprintf('║  HSAF Visualization Module            ║\n');
fprintf('║  Generate Publication-Ready Figures   ║\n');
fprintf('╚════════════════════════════════════════╝\n');
fprintf('Start time: %s\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

%% ========================================================================
%% 第一部分：路径配置与数据加载
%% ========================================================================

%% 1.1) 路径配置
base_path = '/home/um202370130/HSAF';
data_path = fullfile(base_path, 'Data');
script_path = fullfile(base_path, 'Software');

% 添加必要路径
addpath(fullfile(script_path, 'Tool_Functions', 'gmt'));
addpath(fullfile(script_path, 'Tool_Functions'));
addpath(fullfile(script_path, 'Tool_Functions', 'm_map'));

%% 1.2) 定位统计结果文件夹
fprintf('[INPUT] Locating statistics directory...\n');

output_base = fullfile(data_path, 'EWH_Output', 'CSR_EWH', 'HSAF_Adaptive');
dirs = dir(output_base);
dirs = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
[~, idx] = max([dirs.datenum]);
hsaf_output_dir = fullfile(output_base, dirs(idx).name);
stats_output_dir = fullfile(hsaf_output_dir, 'Statistics');

if ~exist(stats_output_dir, 'dir')
    error('[ERROR] Statistics directory not found. Run HSAF_Quality_Control.m first.');
end

% 创建可视化输出子文件夹
vis_output_dir = fullfile(stats_output_dir, 'Figures');
if ~exist(vis_output_dir, 'dir'), mkdir(vis_output_dir); end

fprintf('[DIR] Statistics: %s\n', stats_output_dir);
fprintf('[DIR] Figures output: %s\n\n', vis_output_dir);

%% 1.3) 加载统计结果
fprintf('[LOAD] Loading statistical results...\n');
try
    load(fullfile(stats_output_dir, 'HSAF_Quality_Statistics.mat'));
    fprintf('[LOAD] ✓ Statistics loaded\n');
catch ME
    error('[ERROR] Failed to load statistics: %s', ME.message);
end

% 加载原始HSAF输出
fprintf('[LOAD] Loading HSAF filtered data...\n');
try
    load(fullfile(hsaf_output_dir, 'CSR_EWH_HSAF_filtered.mat'), ...
         'CSR_EWH', 'csr_lon', 'csr_lat', 'dateTime');
    fprintf('[LOAD] ✓ HSAF data loaded\n\n');
catch ME
    error('[ERROR] Failed to load HSAF output: %s', ME.message);
end

% 加载辅助数据
aux_path = fullfile(data_path, 'Auxiliary_Data');
load(fullfile(aux_path, 'rivers1.mat'), 'rivers_new');

%% 1.4) 提取关键变量
lon = csr_lon(:);
lat = csr_lat(:);
num_methods = length(filter_methods);
[Nlon, Nlat, ~] = size(CSR_EWH.None);

% 月份名称
month_names = arrayfun(@(d) datestr(d, 'yyyy-mm'), dateTime, 'UniformOutput', false);

%% ========================================================================
%% 第二部分：滤波前后对比图
%% ========================================================================

fprintf('[PLOT] ═══════════════════════════════════════\n');
fprintf('       Generating Comparison Figures\n');
fprintf('       ═══════════════════════════════════════\n\n');

%% 2.1) 多月多方法对比图
fprintf('[PLOT] Creating multi-month multi-method comparison...\n');

% 选择代表性月份
check_months = [2, 12, 22, 32];  % 第2, 12, 22, 32个月

% 选择要显示的方法（排除None和Mascon）
display_methods = filter_methods(~ismember(filter_methods, {'None'}));

fig = figure('Position', [100 100 1800 1000], 'Visible', 'off');
tiledlayout(length(check_months), length(display_methods), ...
            "TileSpacing", "compact", "Padding", "compact");

for i = 1:length(check_months)
    month_idx = check_months(i);
    
    for j = 1:length(display_methods)
        method = display_methods{j};
        
        nexttile;
        
        if strcmp(method, 'Mascon')
            data = CSR_EWH.Mascon(:, :, month_idx);
        elseif strcmp(method, 'Hankel')
            data = CSR_EWH.HankelLatAdaptive(:, :, month_idx);
        else
            if isfield(CSR_EWH, method)
                data = CSR_EWH.(method)(:, :, month_idx);
            else
                continue;
            end
        end
        
        plot_map(data, lon, lat, 1);
        
        % 第一行添加方法名
        if i == 1
            title(method, 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');
        end
        
        % 第一列添加月份标签
        if j == 1
            ylabel(month_names{month_idx}, 'FontSize', 12, 'FontWeight', 'bold', ...
                   'Units', 'Normalized', 'Position', [-0.15 0.5]);
        end
    end
end

sgtitle('EWH Comparison: Different Filtering Methods', ...
        'FontSize', 18, 'FontWeight', 'bold');
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'EWH (cm)';
cb.FontSize = 12;

print(fig, fullfile(vis_output_dir, 'MultiMonth_MultiMethod_Comparison.png'), '-dpng', '-r300');
fprintf('[PLOT] ✓ Multi-month comparison saved\n');
close(fig);

%% 2.2) 单月详细对比（Before/After/Residual）
fprintf('[PLOT] Creating before-after-residual comparison...\n');

target_month = round(length(dates) / 2);  % 中间月份

fig = figure('Position', [100 100 1800 600], 'Visible', 'off');
tiledlayout(1, 3, "TileSpacing", "compact", "Padding", "compact");

% Before
nexttile;
data_before = CSR_EWH.None(:, :, target_month);
plot_map(data_before, lon, lat, 1);
title('Before Filtering (Original)', 'FontSize', 14, 'FontWeight', 'bold');
clim([-15 15]); colorbar;

% After
nexttile;
data_after = CSR_EWH.HankelLatAdaptive(:, :, target_month);
plot_map(data_after, lon, lat, 1);
title('After HSAF Filtering', 'FontSize', 14, 'FontWeight', 'bold');
clim([-15 15]); colorbar;

% Residual (Removed Noise)
nexttile;
data_residual = data_before - data_after;
plot_map(data_residual, lon, lat, 1);
title('Removed Stripe Noise', 'FontSize', 14, 'FontWeight', 'bold');
clim([-15 15]); colorbar;

sgtitle(sprintf('HSAF Filtering Result (%s)', month_names{target_month}), ...
        'FontSize', 16, 'FontWeight', 'bold');

print(fig, fullfile(vis_output_dir, 'HSAF_BeforeAfter_Residual.png'), '-dpng', '-r300');
fprintf('[PLOT] ✓ Before-after-residual comparison saved\n');
close(fig);

%% ========================================================================
%% 第三部分：球谐误差谱图
%% ========================================================================

fprintf('[PLOT] Creating spherical harmonic error spectra...\n');

time_idx = round(length(dates) / 2);  % 选择中间月份
m1 = -60:60;
l1 = 0:60;
[m, l] = meshgrid(m1, l1);

fig = figure('Position', [100 100 1800 800], 'Visible', 'off');

% 动态布局
n_display = sum(~strcmp(filter_methods, 'None'));
cols = ceil(sqrt(n_display));
rows = ceil(n_display / cols);

tiledlayout(rows, cols, "TileSpacing", "compact", "Padding", "compact");

plot_idx = 1;
for i = 1:num_methods
    method = filter_methods{i};
    if strcmp(method, 'None'), continue; end
    
    nexttile;
    
    sc_data = flipud(SC_results.(method)(:, :, time_idx));
    sc_data(sc_data == 0) = NaN;
    
    pcolor(m, l, abs(sc_data) ./ 1e-11);
    shading flat;
    clim([0, 1]);
    colormap('jet');
    
    title(method, 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'none');
    xlabel('Order (m)', 'FontSize', 10);
    ylabel('Degree (l)', 'FontSize', 10);
    set(gca, 'YDir', 'reverse', 'FontSize', 10);
    xticks([-60, -30, 0, 30, 60]);
    yticks([0, 30, 60]);
    grid on;
    
    plot_idx = plot_idx + 1;
end

% 统一色标
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = '10^{-11}';
cb.FontSize = 12;

sgtitle(sprintf('Spherical Harmonic Error Spectra (%s)', month_names{time_idx}), ...
        'FontSize', 16, 'FontWeight', 'bold');

print(fig, fullfile(vis_output_dir, 'SC_Error_Spectra_Comparison.png'), '-dpng', '-r300');
fprintf('[PLOT] ✓ SC error spectra saved\n');
close(fig);

%% ========================================================================
%% 第四部分：趋势分布图
%% ========================================================================

fprintf('[PLOT] Creating global trend maps...\n');

fig = figure('Position', [100 100 1800 1000], 'Visible', 'off');
tiledlayout(2, 3, "TileSpacing", "compact", "Padding", "compact");

plot_idx = 1;
for i = 1:num_methods
    method = filter_methods{i};
    if strcmp(method, 'None'), continue; end
    
    nexttile;
    
    trend_data = global_ewh.(['Trend_', method]);
    plot_map(trend_data, lon, lat, 1);
    
    title(method, 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');
    clim([-2 2]);
    colorbar;
    
    plot_idx = plot_idx + 1;
    if plot_idx > 6, break; end
end

sgtitle('Global EWH Trend Distribution (cm/year)', ...
        'FontSize', 16, 'FontWeight', 'bold');

print(fig, fullfile(vis_output_dir, 'Global_Trend_Comparison.png'), '-dpng', '-r300');
fprintf('[PLOT] ✓ Global trend maps saved\n');
close(fig);

%% ========================================================================
%% 第五部分：流域分析图
%% ========================================================================

fprintf('[PLOT] Creating basin analysis figures...\n');

%% 5.1) 流域EWH分布（选择主要方法）
main_method = 'Hankel';
if isfield(basin(1), ['ewh_', main_method])
    
    fig = figure('Position', [100 100 1800 1200], 'Visible', 'off');
    tiledlayout(3, 4, "TileSpacing", "compact", "Padding", "compact");
    
    for i = 1:min(12, length(basin))
        nexttile;
        
        data = basin(i).(['ewh_', main_method])(:, :, 1);
        plot_map(data, basin(i).Lon, basin(i).Lat, 1);
        hold on;
        m_plot(rivers_new(i).X, rivers_new(i).Y, 'LineWidth', 2.0, 'Color', 'r');
        
        title(basin(i).name, 'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'none');
    end
    
    sgtitle(sprintf('Basin EWH Distribution (%s Filter)', main_method), ...
            'FontSize', 16, 'FontWeight', 'bold');
    
    print(fig, fullfile(vis_output_dir, 'Basin_EWH_Distribution_Hankel.png'), '-dpng', '-r300');
    fprintf('[PLOT] ✓ Basin EWH distribution saved\n');
    close(fig);
end

%% 5.2) 流域趋势对比
fig = figure('Position', [100 100 1800 1200], 'Visible', 'off');
tiledlayout(3, 4, "TileSpacing", "compact", "Padding", "compact");

for i = 1:min(12, length(basin))
    nexttile;
    
    if isfield(basin(i), 'Trend_Hankel')
        trend_data = basin(i).Trend_Hankel;
        plot_map(trend_data, basin(i).Lon, basin(i).Lat, 1);
        hold on;
        m_plot(rivers_new(i).X, rivers_new(i).Y, 'LineWidth', 2.0, 'Color', 'r');
        
        title(basin(i).name, 'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'none');
        clim([-1 1]);
        colorbar;
    end
end

sgtitle('Basin EWH Trend (cm/year, HSAF Filtered)', ...
        'FontSize', 16, 'FontWeight', 'bold');

print(fig, fullfile(vis_output_dir, 'Basin_Trend_Hankel.png'), '-dpng', '-r300');
fprintf('[PLOT] ✓ Basin trend maps saved\n');
close(fig);

%% ========================================================================
%% 第六部分：统计指标对比柱状图
%% ========================================================================

fprintf('[PLOT] Creating statistical metrics bar charts...\n');

fig = figure('Position', [100 100 1400 800], 'Visible', 'off');
tiledlayout(2, 3, "TileSpacing", "compact", "Padding", "compact");

% 筛选非None和Mascon方法
plot_methods = filter_methods(~ismember(filter_methods, {'None', 'Mascon'}));
colors_bar = lines(length(plot_methods));

% CC
nexttile;
cc_means = arrayfun(@(i) mean(CC_results.(plot_methods{i})), 1:length(plot_methods));
bar(cc_means, 'FaceColor', 'flat', 'CData', colors_bar);
set(gca, 'XTickLabel', plot_methods, 'XTickLabelRotation', 45);
ylabel('Mean CC', 'FontSize', 12, 'FontWeight', 'bold');
title('Correlation Coefficient', 'FontSize', 14);
grid on;

% SNR
nexttile;
snr_means = arrayfun(@(i) mean(SNR_results.(plot_methods{i})), 1:length(plot_methods));
bar(snr_means, 'FaceColor', 'flat', 'CData', colors_bar);
set(gca, 'XTickLabel', plot_methods, 'XTickLabelRotation', 45);
ylabel('Mean SNR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('Signal-to-Noise Ratio', 'FontSize', 14);
grid on;

% RMSE
nexttile;
rmse_means = arrayfun(@(i) mean(RMSE_results.(plot_methods{i})), 1:length(plot_methods));
bar(rmse_means, 'FaceColor', 'flat', 'CData', colors_bar);
set(gca, 'XTickLabel', plot_methods, 'XTickLabelRotation', 45);
ylabel('Mean RMSE (cm)', 'FontSize', 12, 'FontWeight', 'bold');
title('Root Mean Square Error', 'FontSize', 14);
grid on;

% PSNR
nexttile;
psnr_means = arrayfun(@(i) mean(PSNR_results.(plot_methods{i})), 1:length(plot_methods));
bar(psnr_means, 'FaceColor', 'flat', 'CData', colors_bar);
set(gca, 'XTickLabel', plot_methods, 'XTickLabelRotation', 45);
ylabel('Mean PSNR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
title('Peak Signal-to-Noise Ratio', 'FontSize', 14);
grid on;

% MAE
nexttile;
mae_means = arrayfun(@(i) mean(MAE_results.(plot_methods{i})), 1:length(plot_methods));
bar(mae_means, 'FaceColor', 'flat', 'CData', colors_bar);
set(gca, 'XTickLabel', plot_methods, 'XTickLabelRotation', 45);
ylabel('Mean MAE (cm)', 'FontSize', 12, 'FontWeight', 'bold');
title('Mean Absolute Error', 'FontSize', 14);
grid on;

% NSC
nexttile;
nsc_means = arrayfun(@(i) mean(NSC_results.(plot_methods{i})), 1:length(plot_methods));
bar(nsc_means, 'FaceColor', 'flat', 'CData', colors_bar);
set(gca, 'XTickLabel', plot_methods, 'XTickLabelRotation', 45);
ylabel('Mean NSC', 'FontSize', 12, 'FontWeight', 'bold');
title('Nash-Sutcliffe Coefficient', 'FontSize', 14);
grid on;

sgtitle('Statistical Metrics Comparison (Mean Values)', ...
        'FontSize', 16, 'FontWeight', 'bold');

print(fig, fullfile(vis_output_dir, 'Statistics_BarChart_Comparison.png'), '-dpng', '-r300');
fprintf('[PLOT] ✓ Statistics bar chart saved\n');
close(fig);

%% ========================================================================
%% 第七部分：空间相关系数分布图
%% ========================================================================

fprintf('[PLOT] Creating spatial correlation coefficient maps...\n');

fig = figure('Position', [100 100 1800 600], 'Visible', 'off');
tiledlayout(1, min(3, num_methods), "TileSpacing", "compact", "Padding", "compact");

plot_idx = 1;
for i = 1:num_methods
    method = filter_methods{i};
    if strcmp(method, 'None') || strcmp(method, 'Mascon'), continue; end
    
    nexttile;
    
    if isfield(global_cc, method)
        plot_map(global_cc.(method), lon, lat, 1);
        title(method, 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');
        clim([0 1]);
        colorbar;
    end
    
    plot_idx = plot_idx + 1;
    if plot_idx > 3, break; end
end

sgtitle('Global Spatial Correlation Coefficient (vs Mascon)', ...
        'FontSize', 16, 'FontWeight', 'bold');

print(fig, fullfile(vis_output_dir, 'Global_Correlation_Maps.png'), '-dpng', '-r300');
fprintf('[PLOT] ✓ Correlation maps saved\n');
close(fig);

%% ========================================================================
%% 第八部分：SRMSE分布图
%% ========================================================================

fprintf('[PLOT] Creating SRMSE distribution maps...\n');

fig = figure('Position', [100 100 1800 600], 'Visible', 'off');
tiledlayout(1, min(3, num_methods), "TileSpacing", "compact", "Padding", "compact");

plot_idx = 1;
for i = 1:num_methods
    method = filter_methods{i};
    if strcmp(method, 'None') || strcmp(method, 'Mascon'), continue; end
    
    nexttile;
    
    if isfield(SRMSE_results, method)
        plot_map(SRMSE_results.(method), lon, lat, 1);
        title(method, 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');
        clim([0 10]);
        colorbar;
    end
    
    plot_idx = plot_idx + 1;
    if plot_idx > 3, break; end
end

sgtitle('Spatial Root Mean Square Error (SRMSE, cm)', ...
        'FontSize', 16, 'FontWeight', 'bold');

print(fig, fullfile(vis_output_dir, 'SRMSE_Distribution.png'), '-dpng', '-r300');
fprintf('[PLOT] ✓ SRMSE distribution saved\n');
close(fig);

%% ========================================================================
%% 第九部分：生成可视化索引
%% ========================================================================

fprintf('\n[INDEX] Generating figure index...\n');

index_file = fullfile(vis_output_dir, 'Figure_Index.txt');
fid = fopen(index_file, 'w');

fprintf(fid, '═══════════════════════════════════════════════════\n');
fprintf(fid, '  HSAF V6 Visualization Figure Index\n');
fprintf(fid, '═══════════════════════════════════════════════════\n\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now));

fprintf(fid, 'Figure List:\n');
fprintf(fid, '  1. MultiMonth_MultiMethod_Comparison.png\n');
fprintf(fid, '     - Multi-month comparison of different filtering methods\n\n');
fprintf(fid, '  2. HSAF_BeforeAfter_Residual.png\n');
fprintf(fid, '     - Before/After/Residual comparison for HSAF\n\n');
fprintf(fid, '  3. SC_Error_Spectra_Comparison.png\n');
fprintf(fid, '     - Spherical harmonic error spectra\n\n');
fprintf(fid, '  4. Global_Trend_Comparison.png\n');
fprintf(fid, '     - Global EWH trend distribution\n\n');
fprintf(fid, '  5. Basin_EWH_Distribution_Hankel.png\n');
fprintf(fid, '     - Basin-scale EWH distribution\n\n');
fprintf(fid, '  6. Basin_Trend_Hankel.png\n');
fprintf(fid, '     - Basin-scale trend maps\n\n');
fprintf(fid, '  7. Statistics_BarChart_Comparison.png\n');
fprintf(fid, '     - Statistical metrics bar chart comparison\n\n');
fprintf(fid, '  8. Global_Correlation_Maps.png\n');
fprintf(fid, '     - Global spatial correlation coefficient\n\n');
fprintf(fid, '  9. SRMSE_Distribution.png\n');
fprintf(fid, '     - Spatial RMSE distribution\n\n');

fprintf(fid, 'Total: 9 figures\n');
fprintf(fid, 'Location: %s\n', vis_output_dir);

fclose(fid);
fprintf('[INDEX] ✓ Figure index saved: %s\n\n', index_file);

%% ========================================================================
%% 第十部分：最终报告
%% ========================================================================

fprintf('╔════════════════════════════════════════╗\n');
fprintf('║  Visualization Summary                ║\n');
fprintf('╚════════════════════════════════════════╝\n\n');

fprintf('Generated Figures:\n');
fprintf('  1. Multi-month comparison\n');
fprintf('  2. Before-after-residual\n');
fprintf('  3. SC error spectra\n');
fprintf('  4. Global trend maps\n');
fprintf('  5. Basin EWH distribution\n');
fprintf('  6. Basin trend maps\n');
fprintf('  7. Statistics bar charts\n');
fprintf('  8. Correlation maps\n');
fprintf('  9. SRMSE distribution\n\n');

fprintf('Output Location:\n');
fprintf('  %s\n\n', vis_output_dir);

fprintf('════════════════════════════════════════════\n');
fprintf('Visualization completed successfully!\n');
fprintf('End time: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('════════════════════════════════════════════\n\n');


%% ========================================================================
%% 局部函数定义
%% ========================================================================

function param = calibrate_lat_params_v5(EWH, lon, lat, t_list, opt)
[Nlon, Nlat, Nt] = size(EWH);
dlon = mean(diff(lon));
Fs = 1/dlon;

wl_peak = nan(Nlat,1);
p_arr   = nan(Nlat,1);
k_arr   = nan(Nlat,1);
wl_band = nan(Nlat,2);

for j = 1:Nlat
    % 多月平均PSD
    pxx_sum = 0; cnt = 0;
    for tt = t_list(:)'
        if tt<1 || tt>Nt, continue; end
        x = squeeze(EWH(:,j,tt));
        if all(isnan(x)), continue; end
        x = preprocess_1d_v5(x, opt);
        [pxx, f] = pwelch(x, hamming(opt.psd_win), opt.psd_ov, opt.psd_nfft, Fs);
        pxx_sum = pxx_sum + pxx;
        cnt = cnt + 1;
    end

    if cnt==0
        wl_peak(j) = opt.default_wl_peak_deg;
    else
        pxx_m = pxx_sum / cnt;
        wl = 1./f; wl(1) = inf;
        mask = (wl >= opt.peak_wl_range_deg(1)) & (wl <= opt.peak_wl_range_deg(2));
        if ~any(mask)
            wl_peak(j) = opt.default_wl_peak_deg;
        else
            wl_seg = wl(mask); p_seg = pxx_m(mask);
            [~, im] = max(p_seg);
            wl_peak(j) = wl_seg(im);
        end
    end

    % 计算P
    period_samples = max(2, wl_peak(j)/dlon);
    p0 = round(opt.p_factor * period_samples);
    p0 = max(opt.p_min, min(p0, opt.p_max));
    p0 = min(p0, Nlon-2);
    p_arr(j) = p0;

    % 计算K
    xbar = squeeze(mean(EWH(:,j,t_list), 3, 'omitnan'));
    if all(isnan(xbar))
        k0 = opt.k_min;
    else
        xbar = preprocess_1d_v5(xbar, opt);
        L = Nlon + 1 - p0;
        H = hankel(xbar(1:L), xbar(L:Nlon));
        s = svd(H, 'econ');

        if strcmpi(opt.k_determination, 'eigenvalue_gap')
            s_norm = s / (s(1)+eps);
            s_log = log10(s_norm + 1e-12);
            gaps = -diff(s_log);
            k0 = find(gaps > opt.sv_gap_threshold, 1, 'first');
            if isempty(k0), k0 = opt.k_min; end
        else
            e = cumsum(s.^2) / sum(s.^2);
            k0 = find(e >= opt.energy, 1, 'first');
            if isempty(k0), k0 = opt.k_min; end
        end

        k0 = max(opt.k_min, min(k0, opt.k_max));
        k0 = min(k0, p0-1);
        if opt.force_even_k && mod(k0,2)==1
            k0 = min(k0+1, p0-1);
        end
    end
    k_arr(j) = k0;

    % 计算波长带宽
    b1 = wl_peak(j) * opt.band_scale(1);
    b2 = wl_peak(j) * opt.band_scale(2);
    b1 = max(b1, opt.band_abs_deg(1));
    b2 = min(b2, opt.band_abs_deg(2));
    if b2 <= b1
        b1 = opt.band_abs_deg(1);
        b2 = opt.band_abs_deg(2);
    end
    wl_band(j,:) = [b1, b2];
end

% 纬向平滑
wl_peak = smoothdata(wl_peak, 'movmedian', opt.smooth_lat);
p_arr   = round(smoothdata(double(p_arr), 'movmedian', opt.smooth_lat));
% K：两步平滑 + 限坡（避免纬向跳变导致分段伪影）
k_s = smoothdata(double(k_arr), 'movmedian', opt.smooth_lat);
k_s = smoothdata(k_s, 'movmean', opt.smooth_lat);

dk_max = 2;  % 相邻纬度K变化不超过2（可调：2~3）
for jj = 2:numel(k_s)
    k_s(jj) = min(max(k_s(jj), k_s(jj-1)-dk_max), k_s(jj-1)+dk_max);
end
k_arr = round(k_s);

% 极区固定
if opt.polar_fixed_params
    lat_abs = abs(lat);
    polar = lat_abs > opt.polar_lat_threshold;
    ref   = (lat_abs >= 60) & (lat_abs <= opt.polar_lat_threshold);
    if any(polar) && any(ref)
        wl_peak(polar) = median(wl_peak(ref));
        p_arr(polar)   = round(median(p_arr(ref)));
        k_arr(polar)   = round(median(k_arr(ref)));
        wl_band(polar,1) = median(wl_band(ref,1));
        wl_band(polar,2) = median(wl_band(ref,2));
    end
end

% 约束
p_arr = max(opt.p_min, min(p_arr, opt.p_max));
p_arr = min(p_arr, Nlon-2);
k_arr = max(opt.k_min, min(k_arr, opt.k_max));
k_arr = min(k_arr, p_arr-1);
if opt.force_even_k
    odd = mod(k_arr,2)==1;
    k_arr(odd) = min(k_arr(odd)+1, p_arr(odd)-1);
end

param.wl_peak_deg = wl_peak;
param.wl_band_deg = wl_band;
param.p = p_arr;
param.k = k_arr;
end

function x = preprocess_1d_v5(x, opt)
x = x(:);
if any(isnan(x))
    x = fillmissing(x, 'linear', 'EndValues', 'nearest');
end
x = x - mean(x);
if isfield(opt,'taper_alpha') && opt.taper_alpha > 0
    x = x .* tukeywin(numel(x), opt.taper_alpha);
end
end

% [其他辅助函数：hankel_destripe_oneprofile_sw, local_htls_destripe, 等]
% [完整实现请参考原HSAF_V6.m中的相应函数]
%--------------------------------------------------------------------------
% function plot_map(data, lon, lat, interp_factor)
% % 简易地图绘制（集群兼容）
% if nargin < 4, interp_factor = 1; end
% 
% [LON, LAT] = meshgrid(lon, lat);
% 
% if interp_factor > 1
%     lon2 = linspace(min(lon), max(lon), length(lon)*interp_factor);
%     lat2 = linspace(min(lat), max(lat), length(lat)*interp_factor);
%     [LON2, LAT2] = meshgrid(lon2, lat2);
%     data2 = interp2(LON, LAT, data', LON2, LAT2, 'linear');
%     pcolor(LON2, LAT2, data2); shading interp;
% else
%     pcolor(LON, LAT, data'); shading flat;
% end
% 
% axis equal tight;
% xlabel('Longitude (°)', 'FontSize', 10, 'FontWeight', 'bold');
% ylabel('Latitude (°)', 'FontSize', 10, 'FontWeight', 'bold');
% set(gca, 'FontSize', 10);
% end

%--------------------------------------------------------------------------
% 以下函数与HSAF_V4保持一致（滑动窗口去条带）
%--------------------------------------------------------------------------
function x_clean = hankel_destripe_oneprofile_sw(x, Ts, p_lat, k_lat, wl_band_deg, opt)
% 滑动窗口HTLS去条带（与V4版本一致）
% （完整代码见HSAF_V4.txt，此处省略以节省篇幅）

    % 默认参数
    if nargin < 6, opt = struct(); end
    opt = set_default(opt, 'detrend_mode', 'constant');
    opt = set_default(opt, 'taper_alpha', 0.02);
    opt = set_default(opt, 'force_conjugate_pairs', true);
    opt = set_default(opt, 'pair_tol', 0.01);
    opt = set_default(opt, 'use_sliding', true);
    opt = set_default(opt, 'win_len', []);
    opt = set_default(opt, 'win_min', 30);
    opt = set_default(opt, 'win_overlap', 0.75);
    opt = set_default(opt, 'step', []);
    opt = set_default(opt, 'circular', true);
    opt = set_default(opt, 'p_cap_ratio', 1/3);
    opt = set_default(opt, 'p_min_win', 24);
    opt = set_default(opt, 'k_per_window', false);
    opt = set_default(opt, 'k_energy', 0.95);
    opt = set_default(opt, 'k_min', 6);
    opt = set_default(opt, 'k_max', 20);
    opt = set_default(opt, 'min_mode_energy_ratio', 0.0);
    opt = set_default(opt, 'protect_wl_gt_deg', inf);
    opt = set_default(opt, 'ola_window', 'hann');
    opt = set_default(opt, 'ola_tukey_alpha', 0.25);

    x_in = x;
    x = x(:);
    N = numel(x);

    if any(isnan(x))
        x = fillmissing(x, 'linear', 'EndValues', 'nearest');
    end

    if ~opt.use_sliding
        x_clean = local_htls_destripe(x, Ts, p_lat, k_lat, wl_band_deg, opt);
        x_clean = reshape_like(x_clean, x_in);
        return;
    end

    % 窗口长度
    wl2 = wl_band_deg(2);
    T1 = round(3 * p_lat);
    T2 = ceil(1.5 * (wl2 / Ts));
    if isempty(opt.win_len)
        T = max([opt.win_min, T1, T2]);
        T = min(T, N);
    else
        T = min(max(opt.win_len, opt.win_min), N);
    end

    if isempty(opt.step)
        step = max(1, round(T * (1 - opt.win_overlap)));
    else
        step = max(1, round(opt.step));
    end

    w = make_ola_window(T, opt);
    w = w(:);

    acc  = zeros(N,1);
    wsum = zeros(N,1);

    for s = 1:step:N
        idx = mod((s-1):(s+T-2), N) + 1;
        seg = x(idx);

        p_win = min(p_lat, floor(T * opt.p_cap_ratio));
        p_win = max(p_win, opt.p_min_win);
        p_win = min(p_win, T-2);
        k_win = min(k_lat, p_win-1);

        if opt.k_per_window
            k_win = estimate_k_svd(seg, p_win, opt);
        else
            k_win = min(max(k_win, opt.k_min), opt.k_max);
        end

        seg_clean = local_htls_destripe(seg, Ts, p_win, k_win, wl_band_deg, opt);

        acc(idx)  = acc(idx)  + seg_clean(:) .* w;
        wsum(idx) = wsum(idx) + w;
    end

    x_clean = acc ./ max(wsum, eps);
    x_clean = reshape_like(x_clean, x_in);
end

function y = local_htls_destripe(x, Ts, p, k, wl_band_deg, opt)
    x = x(:);

    if strcmpi(opt.detrend_mode,'linear')
        n = (0:numel(x)-1)';
        pp = polyfit(n, x, 1);
        trend = polyval(pp, n);
    else
        trend = mean(x) * ones(size(x));
    end
    xd = x - trend;

    if opt.taper_alpha > 0
        xd = xd .* tukeywin(numel(xd), opt.taper_alpha);
    end

    try
        [~, ~, freq, ~, Y, ~] = H_RCs(xd, Ts, p, k);
    catch
        y = x;
        return;
    end

    freq = freq(:);
    wl = 1 ./ abs(freq);
    wl(~isfinite(wl)) = inf;

    idx = (wl >= wl_band_deg(1)) & (wl <= wl_band_deg(2));
    idx = idx & (wl <= opt.protect_wl_gt_deg);

    if opt.force_conjugate_pairs && any(idx)
        idx2 = false(size(idx));
        for i = find(idx)'
            idx2(i) = true;
            [dmin, j] = min(abs(freq - (-freq(i))));
            if j ~= i && dmin < opt.pair_tol
                idx2(j) = true;
            end
        end
        idx = idx2;
    end

    if ~any(idx)
        y = x;
        return;
    end

    if opt.min_mode_energy_ratio > 0
        e = mean(abs(Y).^2, 2);
        thr = opt.min_mode_energy_ratio * median(e);
        idx = idx & (e >= thr);
        if ~any(idx)
            y = x;
            return;
        end
    end

    noise = sum(Y(idx, :), 1).';
    y = (xd - noise) + trend;
end

function k = estimate_k_svd(x, p, opt)
    x = x(:);
    N = numel(x);
    L = N - p + 1;
    if L < 2
        k = opt.k_min; return;
    end
    H = hankel(x(1:L), x(L:N));
    s = svd(H, 'econ');
    e = cumsum(s.^2) / max(sum(s.^2), eps);
    k = find(e >= opt.k_energy, 1, 'first');
    if isempty(k), k = opt.k_min; end
    k = min(max(k, opt.k_min), min(opt.k_max, p-1));
    if mod(k,2)==1, k = k + 1; end
    k = min(k, p-1);
end

function w = make_ola_window(T, opt)
    switch lower(opt.ola_window)
        case 'tukey'
            w = tukeywin(T, opt.ola_tukey_alpha);
        otherwise
            w = hann(T, 'periodic');
    end
    w = max(w, eps);
end

function opt = set_default(opt, name, val)
    if ~isfield(opt, name) || isempty(opt.(name))
        opt.(name) = val;
    end
end

function y = reshape_like(ycol, xlike)
    if isrow(xlike)
        y = ycol(:).';
    else
        y = ycol(:);
    end
end 