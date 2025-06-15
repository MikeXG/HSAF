%clc; clear; warning('off', 'verbose');
fprintf('\n📌 [Hankel_Filter] 进行 Hankel 滤波...\n');

% % 1. 加载必要数据
EWH_data_address = 'G:\HSAF\Data\EWH_Output\';
addpath 'G:\HSAF\Data\Auxiliary_Data'
addpath 'G:\HSAF\Software\HSA_Filter\'
data_path = 'G:\HSAF\Data\';
script_path = 'G:\HSAF\Software\';

addpath(fullfile(script_path, 'Tool_Functions', 'gmt'));
addpath(fullfile(script_path, 'Tool_Functions'));

% opengl("save", 'software');
% 
% % 2. 询问用户选择数据类型
% disp('请选择要处理的数据类型:');
% disp('1. CSR 数据 (GRACE 观测)');
% disp('2. HIS 数据 (模拟数据)');
% data_choice = input('请输入数字 1 或 2: ');

% % **根据用户选择加载相应数据**
% if data_choice == 1
%     load(fullfile(EWH_data_address, 'CSR_EWH_data.mat'));  % HIS 数据
%     EWH = CSR_EWH;
%     data_type = 'CSR';
% 
% elseif data_choice == 2
%     load(fullfile(EWH_data_address, 'HIS_EWH_data.mat'));  % HIS 数据
%     EWH = HIS_EWH;
%     data_type = 'HIS';
% else
%     error('输入无效！请重新运行程序，并输入 1 (CSR) 或 2 (HIS)。');
% end


% 3. 用户输入初始参数
Ts = 30;  % 采样周期
window_size = input('请输入滑动窗口大小 window_size: ');
p = input('请输入嵌套维数 p: ');
order = input('请输入阶数 order: ');
buffer = input('请输入缓冲参数 buffer: ');

% % **用户选择起始时间**
start_year = input('请输入起始年份 (如 2004年): ');
start_month = input('请输入起始月份 (如 1月): ');
end_year = input('请输入结束年份 (如 2016年): ');
end_month = input('请输入结束月份 (如 1月): ');

% start_year = 2004;
% start_month = 1;
% end_year = 2016;
% end_month = 1;


% **确定起始索引和结束索引**
start_idx = find(year(dateTime) == start_year & month(dateTime) == start_month, 1);
end_idx = find(year(dateTime) == end_year & month(dateTime) == end_month, 1);
time_size = end_idx - start_idx + 1;  % 计算时间步长

dates = dateTime(start_idx:end_idx); % 选定时间范围

% 生成输出目录
time_now = string(strrep(datestr(now,15), ':', '-'));
folder_name = sprintf('window_size=%dp=%dorder=%d,window_size_s=%dp_s=%dorder_s=%d', window_size, p, order, window_size, p, order);
Dir_images_output = fullfile(data_path, 'EWH_Output', [data_type, '_EWH'], 'EWH_HSA_Filter', time_now, folder_name);
mkdir(Dir_images_output);

% 备份当前脚本
script_name = 'G:\HSAF\Software\HSAF_HUST.m'; 
snapshot_name = sprintf('Snapshot_window_size=%d_p=%d_order=%d.txt', window_size, p, order);
copyfile(script_name, fullfile(Dir_images_output, snapshot_name));
disp(['脚本快照已保存为: ', snapshot_name]);

% 5. 预分配存储变量
filter_methods = fieldnames(EWH); % 获取所有滤波方法的字段名
filter_results = struct();

% **球谐系数存储 (60×60×time_size)**
SC_results = struct();
for method = filter_methods
    SC_results.(method{1}) = zeros(61,121,time_size);
end




% 6. 进行第一次 Hankel 滤波
tic;
parfor t = 1:time_size
    month_idx = start_idx + (t - 1); % **自动计算当前月份索引**
    % Hankel 滤波计算
    Hankel_Mode= HSA(EWH.Decorrelation(:,:,month_idx), Ts, window_size, p,order, buffer);
    filter_results_Hankel(:,:,t) = EWH.Decorrelation(:,:,month_idx) - (sum(Hankel_Mode(:,:,1:6),3) - sum(Hankel_Mode(:,:,3:4),3));
end
EWH.Hankel = filter_results_Hankel;
plot_map(EWH.Hankel(:,:,1),csr_lon,csr_lat,1);colorbar;
toc;

% 6. 进行第一次 Hankel 滤波

% 3. 询问用户是否进行二次迭代
% iter_flag = input('是否进行二次迭代？(1: 是, 0: 否) ');

% if iter_flag == 1
%     % 允许用户修改参数
%     window_size = input('请输入新的滑动窗口大小 window_size: ');
%     p = input('请输入新的嵌套维数 p: ');
%     order = input('请输入新的阶数 order: ');
%     buffer = input('请输入新的缓冲参数 buffer: ');
% 
%     % 4. 进行第二次 Hankel 滤波
%     tic;
%     parfor t = 1:time_size
%         % 以第一次滤波后的结果作为输入
%         Hankel_Mode = HSA(EWH.Hankel(:,:,t), Ts, window_size, p, order, buffer);
%         filter_results_Hankel(:,:,t) = EWH.Hankel(:,:,t) - ...
%             (sum(Hankel_Mode(:,:,1:6),3) - sum(Hankel_Mode(:,:,3:4),3));
%     end
%     EWH.Hankel = filter_results_Hankel;
%     toc;
% end

disp('Hankel 滤波完成。');
%plot_map_layout(Hankel_Mode,csr_lon,csr_lat,1,2,3,month_names)

% 7. 变换到球谐系数 (Hankel -> Spherical Harmonics)
for t = 1:time_size
    for i = 1:length(filter_methods)
        method = filter_methods{i}; % 当前滤波方法的名称
        if method == "Hankel"
            start_id = t;
        else
            start_id = start_idx + (t - 1);
        end
        SC_results.(method)(:,:,t) = gmt_cs2sc(gmt_mc2gc(gmt_grid2cs(EWH.(method)(:,:,start_id) .\ 100', 60)));
    end
end

% 8. 处理 NaN 值
for i = 1:length(filter_methods)
    method = filter_methods{i};
    SC_results.(method)(SC_results.(method) == 0) = NaN;
end

% 9. 保存数据
% save(fullfile(EWH_data_address, [data_type, '_EWH_Hankel_Filtered.mat']), 'EWH', 'SC_results', 'dateTime', 'years');
% disp([data_type, ' Hankel 滤波完成，并已保存！']);
