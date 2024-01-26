%自动识别峰波区间，但需要手动选择峰间距
%使用findpeaks函数寻找峰值点
%去除峰波区间进行插值获取基线
%基线校正后的电流数据在y_new，电压数据在x_raw
clear all;
clc;

%% 选择Excel文件,注意文件后缀是.csv
[file, path] = uigetfile('*.csv', '选择Excel文件');

if isequal(file, 0)
    disp('未选择任何文件');% 检查用户是否取消选择
else
    % 构建文件路径
    filepath = fullfile(path, file);
    % 读取Excel文件中的数据到data
    data = xlsread(filepath);
end

%% 电压x,电流y初值，初始参数选择
x_raw= data(:, 1);
y_raw= data(:, 2);
windowSize=input('滤波窗口长度选择，不处理选1：');%滑动平均参数选择

% 使用smoothdata函数进行滑动平均处理
y_raw1= smoothdata(y_raw, 'movmean', windowSize);
x_raw1=x_raw;%不处理

y=y_raw1;
x=x_raw1;

%% 查找峰值
% 绘制曲线和峰值
figure(1);
plot(x, y);
title('最小峰距参数选择');
xlabel('电压');
ylabel('电流');
%legend('Original Data');
%legend('Location', 'best');
hold on;

minPeakDistance=input('相邻峰最小间距（小一点0.1,大一点0.3）：');%最小峰距参数选择
minPeakHeight = 0;%峰值最低高度
[peaks, locs, widths, proms] = findpeaks(y, x, 'MinPeakHeight', minPeakHeight, 'MinPeakDistance', minPeakDistance);

% 绘制曲线和峰值
plot(locs, peaks, 'ro');
%legend('Original Data','peak');
%legend('Location', 'best');
hold on;

%% 寻找峰区间
% 遍历每个峰值
for i = 1:length(peaks)
    % 当前峰值的位置和高度
    peak_pos = find(x >= locs(i), 1, 'first');%找到峰值对应的位置
    
    % 寻找左侧峰谷
          left_idx(i) = peak_pos - 1;
    while left_idx(i) > 3 && y(left_idx(i)) >= y(left_idx(i)-3)  %判断条件，-3原因是提高容错性
          left_idx(i) = left_idx(i) - 1;
    end
          leftBound(i) = x(left_idx(i)-1);  %弥补-3，提高准确率
    
    % 寻找右侧峰谷
          right_idx(i) = peak_pos + 1;
    while right_idx(i) <= length(y)-3 && y(right_idx(i)) >= y(right_idx(i)+3)  %判断条件,+3原因是提高容错性
          right_idx(i) = right_idx(i) + 1;
    end
          rightBound(i) = x(right_idx(i)+1);  %弥补+3，提高准确率

    line([x(left_idx(i)), x(right_idx(i))], [y(left_idx(i)), y(right_idx(i))], 'Color', 'g', 'LineStyle', '--');%显示波形区间
 
end
legend('Original Data','peak','data');
legend('Location', 'best');

title('识别峰值与峰区间');
hold off;

%% 基线插值

%去除峰波区间
for i = 1:length(peaks)
    % 在原始数据中将区间内的数据值设置为 NaN
    x(left_idx(i):right_idx(i)) = NaN;
    y(left_idx(i):right_idx(i)) = NaN;
end
x = x(~isnan(x));% 去除NAN值
y = y(~isnan(y));

Intertype=input('线性插值选1，样条插值选2，最近邻插值选3：');%插值方法选择
% 对数据进行插值
if Intertype == 1
    % 线性插值
    yInterpolated = interp1(x, y, x_raw, 'linear');

elseif Intertype == 2
    % 样条插值
    yInterpolated = interp1(x, y, x_raw, 'spline');

elseif Intertype == 3
    % 最相邻插值
    yInterpolated = interp1(x, y, x_raw, 'nearest');

else
    disp('未知的命令数值');
    return;
end

% 绘制原始数据和插值后的曲线
figure(2);
plot(x_raw, y_raw1, 'b-', 'DisplayName', 'Original Data');
hold on;
plot(x_raw, yInterpolated, 'r-', 'DisplayName', 'baseline');
xlabel('电压');
ylabel('电流');
title('原始曲线与基线');
legend('Location', 'best');
hold on;

%% 基线校正
y_new=y_raw1-yInterpolated;% 基线校正

%寻找峰值
num_peaks = length(peaks);
disp('峰值信息（从左到右）：');

for i = 1:length(peaks);
    peakdata_y=y_new(left_idx(i):right_idx(i));
    peakdata_x=x_raw(left_idx(i):right_idx(i));
    [max_peak, max_index] = max(peakdata_y);%区间内最大的峰值点

    peak_info(i, 1) = max_peak; % 峰值大小
    peak_info(i, 2) = peakdata_x(max_index); % 峰值对应的坐标值

    % 显示多个峰值的大小和对应的坐标值
    disp(['峰波', num2str(i), '：峰值大小=', num2str(peak_info(i, 1)), ', 峰值电压=', num2str(peak_info(i, 2))]);
end

%输出校正后的峰值以及位置图像
figure(3);
plot(x_raw, y_new, 'r-', 'DisplayName', 'baseline correction');
hold on;
%plot(peak_info(:, 2), peak_info(:, 1), 'ro', 'MarkerSize', 10);
xlabel('电压');
ylabel('电流');
title('基线校正后的曲线');
legend('Location', 'best');
