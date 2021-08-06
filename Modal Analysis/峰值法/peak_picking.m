%峰值法

clear
close all
clc

%提取测试数据
load modal_test1;
data=signal(1,:); 
%所选的时域样本内数据个数
[m,n]=size(data);
ntdata = n; 



%【1】正常信号处理:
Cvtdata = data;
%【1】正常信号处理end;


% %【2】Narada信号转换:
% for j=1:length(data(:,1))
%     data(j,1) = data(j,1) * 5.0 / 65535.0;  %Convert the data 转化为电压信号
% end       
% 
% %增益系数
% gain = 10;
% for j=1:ntdata      
%     Cvtdata(j,1)= data(j,1)/(gain*200/1000);
% end
% Avgdata = zeros(1,1);
% Avgdata = mean( Cvtdata(:,1)  ); %求数据平均值
% %【2】Narada信号转换end;

% 峰值法识别
f = 500;  %采样频率
fs = f;
fc = 0.02;  %高通截断频率
% dummyV = IdealHighPass (DSet, fs, fc) ;     %滤波
% nDSet= dummyV(1:size(DSet,1));   
t = 0:1/f:(ntdata-1)/f;


%%%  FFT 快速傅里叶变换  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Ns,dummy]= size(nDSet(:,1));    %Sample size 样本大小
mag=abs(fft(Cvtdata,ntdata));   %Magnitude 
mag(1:1,1)= mag(1:1,1)/20;
freqdomain= fs*(0:ntdata-1)'/ntdata;   % Frequency series

figure;
plot(t,Cvtdata);
xlabel(['Time (s)']);
ylabel(['Accel (g)']); 
grid on
% annotation('textbox',[0.2,0.8,0.9,0.2],'string',strcat( Narada_name,'  Time:  ',dataFolder),'EdgeColor','none','FontSize',16);

figure;
hold on
title('fft-data');
plot(freqdomain(1:ntdata/2),mag(1:ntdata/2));
xlabel('Frequency (Hz)','FontSize',14);
% ylabel('Magnitude','FontSize',14);
ylabel('Magnitude','FontSize',14);
hold off

%%%Obtain psd estimate using welch's method
if 1==1
figure;
% hold on
title('pwelch-data');
[psd_c,f_c]=pwelch(Cvtdata,ntdata,0,[],fs,'onesided');  %pls set the stablized time of zeropadding as 10s.
pwelch(Cvtdata,ntdata,0,[],fs,'onesided');
figure;
loglog(f_c,psd_c);
xlabel('Frequency (Hz)');
ylabel('Power Spectrum Density ');
% hold off
end



















