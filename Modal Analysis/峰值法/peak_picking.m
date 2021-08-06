%��ֵ��

clear
close all
clc

%��ȡ��������
load modal_test1;
data=signal(1,:); 
%��ѡ��ʱ�����������ݸ���
[m,n]=size(data);
ntdata = n; 



%��1�������źŴ���:
Cvtdata = data;
%��1�������źŴ���end;


% %��2��Narada�ź�ת��:
% for j=1:length(data(:,1))
%     data(j,1) = data(j,1) * 5.0 / 65535.0;  %Convert the data ת��Ϊ��ѹ�ź�
% end       
% 
% %����ϵ��
% gain = 10;
% for j=1:ntdata      
%     Cvtdata(j,1)= data(j,1)/(gain*200/1000);
% end
% Avgdata = zeros(1,1);
% Avgdata = mean( Cvtdata(:,1)  ); %������ƽ��ֵ
% %��2��Narada�ź�ת��end;

% ��ֵ��ʶ��
f = 500;  %����Ƶ��
fs = f;
fc = 0.02;  %��ͨ�ض�Ƶ��
% dummyV = IdealHighPass (DSet, fs, fc) ;     %�˲�
% nDSet= dummyV(1:size(DSet,1));   
t = 0:1/f:(ntdata-1)/f;


%%%  FFT ���ٸ���Ҷ�任  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Ns,dummy]= size(nDSet(:,1));    %Sample size ������С
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



















