clc;
clear;
close all;
%% %%%%%%%%%%%%%%作业六 时域回声抵消器%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 载入系统响应
%fs = 8e3;%采样频率为8kHz
[h ,fs] = audioread('echo_path.wav');
%计算时长
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
%{
figure(6);
subplot(3,1,1);
plot(t,x);
title('远端语音信号x(n)');xlabel('t/s');ylabel('amplitude');grid on;axis([0 40 -1 1]);
subplot(3,1,2);
plot(t,v1 );
title('背景噪声v_1(n)');xlabel('t/s');ylabel('amplitude');grid on;axis([0 40 -1 1]);
subplot(3,1,3);
plot(t,v2);
title('近端语音信号v_2(n)');xlabel('t/s');ylabel('amplitude');grid on;axis([0 40 -1 1]);
figure(7)
subplot(2,1,1);
plot(t,y);
title('回声信号y(n)');xlabel('t/s');ylabel('amplitude');grid on;axis([0 40 -1 1]);
subplot(2,1,2);
plot(t,d);
title('麦克风接收的信号d(n)');xlabel('t/s');ylabel('amplitude');grid on;axis([0 40 -1 1]);
%}

%% 滤波器参数设置

len = 1200;%len为自适应滤波器的阶数

w_record = zeros(len,dotnumber);%w_SD 记录最速下降法每次滤波器的系数
e = zeros(dotnumber,1);
vector_u = zeros(len,1);%vector_input_filter为每次滤波器的输入向量



%% 计算算法的MSD和ERLE曲线
mu = [0.01;0.1];
kesi = [];
wo = h;
for i = 1:2
    MSD = [];

    w = zeros(len,1);%w自适应滤波器的实时系数
    vector_u = [zeros(len,1)]; % 初始化滤波器输入向量
%for mu = 0.01:0,05:0.05*5  %mu为LMS算法步长 0<u<0.6667为了计算步长的影响
    for n = 1:dotnumber  %
         disp(['这是第',num2str(n),'次循环']);
        vector_u = [x(n);vector_u(1:end-1)]; % 更新滤波器输入向量
         filter_output(n) = w' * vector_u; % 滤波器输出
         e(n) = d(n) - filter_output(n); % 误差
     
    
    %%检测双端对讲，当存在双端对讲的时候Mu为0
        beta = 0.9992;
        if n==1
            E_d2(1) = (1-beta)*d(1).^2;
            E_e2(1) = (1-beta)*e(1).^2;
            E_de(1) = (1-beta)*d(1)*e(1);
        else 
            E_d2(n) = beta*E_d2(n-1)+(1-beta)*d(n).^2;
            E_e2(n) = beta*E_e2(n-1)+(1-beta)*e(n).^2;
            E_de(n) = beta*E_de(n-1)+(1-beta)*e(n)*d(n);
        end
        c_de(n) = E_de(n)/sqrt(E_d2(n)*E_e2(n) );
        if (n>1000)
             if c_de(n)>0.5
                w = w;   %不对滤波器系数进行迭代
             else
    
                w = w + mu(i)*conj(e(n))*vector_u/(norm(vector_u).^2+eps);   %滤波器系数进行迭代
             end
        else 
            w = w + mu(i)*conj(e(n))*vector_u/(norm(vector_u).^2+eps);   %滤波器系数进行迭代
        end
        MSD(n) = 10*log10( norm(w - wo).^2 / norm(wo).^2 );%norm()计算向量的2范数
    end


    for n = 1:dotnumber
        ERLE(n) = 10*log10(E_d2(n)/E_e2(n));
    end
%此时ERLE不用再平滑处理
figure(12);
    subplot(2,1,i);
    if i==1
        plot(t,e);
    else 
        plot(t,e,'r');
    end
   % titleStr = sprintf('消除回声后的信号(\mu=%.6f)', mu);
    title(['消除回声后的信号（\mu=',num2str(mu(i)),')']);
    xlabel('t/s');ylabel('amplitude');axis([0 40 -1 1]);
    grid on;
figure(13)
grid on;
hold on;
 plot(t,MSD);
axis([0  40 -10 30]);
    
figure(14)
grid on;
hold on;
plot(t,ERLE);
axis([0  40 -35 40]);
   
end
figure(13)
legend('\mu=0.01', '\mu=0.1');
title('MSD曲线');xlabel('t/s');ylabel('MSD');
figure(14)
legend('\mu=0.01', '\mu=0.1');
title('ERLE曲线');xlabel('t/s');ylabel('ERLE');



