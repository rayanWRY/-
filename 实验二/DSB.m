
clear,clc
% 生成DSB信号
Fs = 100000; % 采样率
T = 1/Fs; % 采样周期
L = 1000; % 信号长度
t = (0:L-1)*T; % 时间向量

f_c = 10000; % 载波频率
f_m = 1000; % 调制信号频率

A_c = 0.5; % 载波幅度
A_m = 1; % 调制信号幅度

carrier = A_c * cos(2*pi*f_c*t); % 载波信号
modulating = A_m * cos(2*pi*f_m*t); % 调制信号

dsb_signal = carrier .* modulating; % DSB信号

% 绘制时域波形图
subplot(2,1,1);
plot(t, dsb_signal);
title('DSB信号时域波形图');
hold on 
plot(t, modulating,'r');
xlabel('时间（秒）');
ylabel('幅度');

% 计算并绘制频谱图
Y = fft(dsb_signal);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

subplot(2,1,2);
plot(f, P1);
title('DSB信号频谱图');
xlabel('频率（Hz）');
ylabel('幅度');

% 计算统计特性
a=[0.5,0.8,1,1.5];
for i=1:4
    Fs = 100000; % 采样率
T = 1/Fs; % 采样周期
L = 1000; % 信号长度
t = (0:L-1)*T; % 时间向量

f_c = 10000; % 载波频率
f_m = 1000; % 调制信号频率

A_c = a(i); % 载波幅度
A_m = 1; % 调制信号幅度

carrier = A_c * cos(2*pi*f_c*t); % 载波信号
modulating = A_m * cos(2*pi*f_m*t); % 调制信号

dsb_signal = carrier .* modulating; % DSB信号
    
rms_value = rms(dsb_signal);                   % 均方根值
variance = var(dsb_signal);                    % 方差
autocorr = xcorr(dsb_signal, 'unbiased');        % 自相关函数
power_density = pwelch(dsb_signal); % 功率谱密度
max_value = max(dsb_signal);     %最大值
min_value = min(dsb_signal);       %最小值
std_deviation = std(dsb_signal);   %标准差

disp(['调制深度：', num2str(a(i))]);
disp(['均方根值：', num2str(rms_value)]);
disp(['方差：', num2str(variance)]);
disp(['最大值：', num2str(max_value)]);
disp(['最小值：', num2str(min_value)]);
disp(['标准差：', num2str(std_deviation)]);
disp(['-----------------']);
figure(3)
hold on
subplot(2,2,i);
plot(t, dsb_signal);
title('DSB信号时域波形图');
hold on 
plot(t, modulating,'r');
xlabel('时间（秒）');
ylabel('幅度');
end

% 画出无偏估计的自相关函数图形
figure(2);
plot(-length(dsb_signal)+1:length(dsb_signal)-1, autocorr);
title('DSB信号无偏估计的自相关函数');
xlabel('延迟');
ylabel('自相关值');