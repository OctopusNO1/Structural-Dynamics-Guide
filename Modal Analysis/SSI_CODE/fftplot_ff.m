function [yfft,f] =fftplot_ff(y,Fs,NFFT,ratio,show)
%y,measured response 
%Fs,sampling frequency
%NFFT,the NO. of data for FFT
%ratio,
if (nargin==4), show='off'; end 

L0=length(y);
if nargin<=2
 L=L0;NFFT=2^nextpow2(L);
 Y=fft(y,NFFT)/L;
% Fs=1/ts;
 f=Fs/2*linspace(0,1,NFFT/2+1);
 yfft=2*abs(Y(1:NFFT/2+1));
else
    
 NN=floor(L0/NFFT);
 hd=NFFT*(1-ratio);
 f=Fs/2*linspace(0,1,NFFT/2+1);
%  f=[1:L/2]*Fs/L;

 ytempfft=zeros(NN,NFFT/2+1);
  for i=1:NN
    tempY=y(hd*(i-1)+1:NFFT+hd*(i-1));
    tempYfft=fft(tempY,NFFT)/NFFT;
    ytempfft(i,:)=2*abs(tempYfft(1:NFFT/2+1));
  end 
  yfft=sum(ytempfft)/NN;
  
end 

TF = strcmp(show,'on');
if TF
figure;
% subplot(212);
plot(f,yfft);
set(gca,'FontName','Times New Roman','FontSize',11);
xlabel('Frequency (Hz)','Fontsize',11,'Fontname','Times');
% ylabel('|Y(f)|')
ylabel('Spectrum')
% title('Single-Sided Amplitude Spectrum')
set(gcf,'Position',[100 100 0.6*560 0.6*420]);
end

end