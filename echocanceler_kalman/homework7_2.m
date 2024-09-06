clc;
clear;
close all;
%% %%%%%%%%%%%%%%作业六 时域回声抵消器%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 载入系统响应
%fs = 8e3;%采样频率为8kHz
[h ,fs] = audioread('echo_path.wav');

Time = 40;%语音信号长度为40s
dotnumber = Time/(1/fs);%采样的点数为320000
t = [0:1:dotnumber-1].*(1/fs);
%% 载入远端信号

%x为远端语音信号x(n)
[x,fs]=audioread('farend_signal.wav');
timeDelay_x = length(x)*(1/fs);
%% 载入近端语音信号v2
[v2,fs] = audioread('nearend_signal.wav');
v2 = [zeros(dotnumber/2,1);v2;zeros(dotnumber/2-length(v2),1)];

%% 远端信号与空间响应卷积产生回声信号

%x为远端语音信号，此时与系统响应h(n)进行卷积
y = conv(x,h) ;
audiowrite('audio_withEcho.wav', y, fs);
%取回声信号40s 远端信号40s（题目要求）
y = y(1:dotnumber);
x = x(1:dotnumber);
%计算回声信号功率
y_power= sum(y.^2) / length(y);

% 计算噪声功率，以达到 30 dB 的 SNR
noise_power = y_power / (10^(30/10));

%生成背景噪声v1
v1 = sqrt(noise_power) * randn(length(y),1);

%d(n) = y(n)+v1(n) +v2(n) 麦克风接收的信号
d = y+v1+v2;           
%% 绘制各个信号波形
figure(8);
subplot(3,1,1);
plot(t,x);
title('远端语音信号x(n)');xlabel('t/s');ylabel('amplitude');grid on;axis([0 40 -1 1]);
subplot(3,1,2);
plot(t,v1 );
title('背景噪声v_1(n)');xlabel('t/s');ylabel('amplitude');grid on;axis([0 40 -1 1]);
subplot(3,1,3);
plot(t,v2);
title('近端语音信号v_2(n)');xlabel('t/s');ylabel('amplitude');grid on;axis([0 40 -1 1]);
figure(9)
subplot(2,1,1);
plot(t,y);
title('回声信号y(n)');xlabel('t/s');ylabel('amplitude');grid on;axis([0 40 -1 1]);
subplot(2,1,2);
plot(t,d);
title('麦克风接收的信号d(n)');xlabel('t/s');ylabel('amplitude');grid on;axis([0 40 -1 1]);



%% Kalman滤波
%考虑仅存在远端的状况
len = 1200;%len为自适应滤波器的阶数

w_record = zeros(len,dotnumber);%w_SD 记录最速下降法每次滤波器的系数
%调试1
%vector_u = zeros(len,1);%vector_input_filter为每次滤波器的输入向量
vector_u = x(len:-1:1);
filter_output = [];
e = zeros(dotnumber,1);
%{
a = 0.9999;

lambda = 0.9999;
%======初始化========
%w自适应滤波器的实时系数
w = ones(len,1);

%e为误差信号
e(1) = d(1);

%过程噪声方差估计sigma^2_delta
E_w2 = 0;

%测量噪声方差估计
sigma2_v = mean(e(1:len).^2);

%p为系统距离============》初始值存在问题
p = (1-a^2)*(1-lambda)*norm(w).^2/len;

%k为卡尔曼增益===========>初始值存在问题
k = p*vector_u/(p*vector_u.'*vector_u+   (1-lambda)*e(1).^2);

x = [x;zeros(len-1,1)];

for n = 1:dotnumber  %

    vector_u = [x(n);vector_u(1:end-1)];        %每次滤波器抽头的输入向量,*100是为了防止浮点数溢出
    filter_output(n) = w.'*vector_u;       %与滤波器系数相乘

      %得到n时刻的误差信号
    e(n) = d(n) - filter_output(n);

    %估计测量噪声sigam^2_v
    sigma2_v = lambda*sigma2_v + (1-lambda)*e(n).^2;

    %卡尔曼增益
    k = p*vector_u/(p* (vector_u.') *vector_u + sigma2_v );
    
    %滤波器系数进行迭代
    w = a*(w + k*e(n));   

    %估计E{||w(n)||^2}和 测量噪声sigma^2_delta
    E_w2 = lambda*E_w2 + (1-lambda)*norm(w).^2;
    sigma2_delta = (1-a^2) * E_w2 /len;

    %系统距离更新
    p = a^2*(1 - vector_u.'*k/len)*p+sigma2_delta;

    %记录此时滤波器系数
    w_record(:,n) = w;  
end

%% 绘制输出后的信号

figure(9)
subplot(2,1,1);
plot(t,d);
title('麦克风输入信号');xlabel('t/s');ylabel('amplitude');grid on;axis([0 40 -1 1]);hold on;
subplot(2,1,2);
plot(t,e);
title('消除回声信号后的信号');xlabel('t/s');ylabel('amplitude');grid on;axis([0 40 -1 1]);hold on;

%保存音频
audiowrite('canceled_audio_file.wav', e, fs);
%}
%% 计算MSD曲线
E_d2 = zeros(1,dotnumber);
E_e2 = zeros(1,dotnumber);
w0 = h;
i =0;
lambda = 0.9999;
% 循环不同的 mu 值
for a = [0.9999 0.99995 0.999975]
      
    disp(['a=',num2str(a)]);
    w = ones(len,1); % 重置滤波器系数为初始值
    vector_u = x(len:-1:1); % 初始化滤波器输入向量
    MSD = zeros(1, dotnumber); % 初始化MSD向量
   % w = ones(len,1);

    %e为误差信号
    e(1) = d(1);

    %过程噪声方差估计sigma^2_delta
    E_w2 = 0;

    %测量噪声方差估计
    sigma2_v = mean(e(1:len).^2);

    %  p为系统距离============》初始值存在问题
    p = (1-a^2)*(1-lambda)*norm(w).^2/len;

    %k为卡尔曼增益===========>初始值存在问题
    k = p*vector_u/(p*vector_u.'*vector_u+   (1-lambda)*e(1).^2);

    for n = 1:dotnumber
        disp(['这是第',num2str(n),'次循环']);
       
    vector_u = [x(n);vector_u(1:end-1)];        %每次滤波器抽头的输入向量,*100是为了防止浮点数溢出
    filter_output(n) = w.'*vector_u;       %与滤波器系数相乘

      %得到n时刻的误差信号
    e(n) = d(n) - filter_output(n);

    %估计测量噪声sigam^2_v
    sigma2_v = lambda*sigma2_v + (1-lambda)*e(n).^2;

    %卡尔曼增益
    k = p*vector_u/(p* (vector_u.') *vector_u + sigma2_v );
    
    %滤波器系数进行迭代
    w = a*(w + k*e(n));   

    %估计E{||w(n)||^2}和 测量噪声sigma^2_delta
    E_w2 = lambda*E_w2 + (1-lambda)*norm(w).^2;
    sigma2_delta = (1-a^2) * E_w2 /len;

    %系统距离更新
    p = a^2*(1 - vector_u.'*k/len)*p+sigma2_delta;

    %记录此时滤波器系数
    w_record(:,n) = w;  
        MSD(n) = 10 * log10( norm(w - w0).^2 / norm(w0).^2); % 计算MSD
    end
    figure(10);
    subplot(3,1,i+1);
    plot(t,e);
    titleStr = sprintf('消除回声后的信号(a=%.6f)', a);
    title(titleStr);
    xlabel('t/s');ylabel('amplitude');axis([0 40 -1 1]);
    grid on;

    figure(11); % 创建绘图窗口
    hold on; % 保持图形，以便在同一图上绘制多条曲线
    grid on; % 开启网格
    plot(t, MSD); % 绘制MSD曲线

    figure(12);
    title('ERLE曲线'); % 设置图表标题
    hold on;
    grid on;
    lambda = 0.9999;
    E_d2(1) = (1-lambda)*d(1).^2;
    E_e2(1) = (1-lambda)*e(1).^2;
  
    ERLE(1) = 10*log10(E_d2(1)/E_e2(1));
    for n = 1:dotnumber-1
    E_d2(n+1)  = lambda*E_d2(n)+(1-lambda)*d(n+1).^2;
    E_e2(n+1) = lambda*E_e2(n)+(1-lambda)*e(n+1).^2;
    ERLE(n+1) = 10*log10(E_d2(n+1)/E_e2(n+1));
    end
    plot(t,ERLE);axis([0 40 -35 40]);
    xlabel('t/s');ylabel('ERLE/dB');
     i = i+1;
end
figure(11)
title('MSD曲线'); % 设置图表标题
xlabel('t/s'); % 设置x轴标签
ylabel('MSD/dB'); % 设置y轴标签
legend('a=0.9999', 'a=0.99995', 'a=0.999975'); % 添加图例
figure(12)
legend('a=0.9999', 'a=0.99995', 'a=0.999975'); % 添加图例
hold off; % 取消持续绘图状态
