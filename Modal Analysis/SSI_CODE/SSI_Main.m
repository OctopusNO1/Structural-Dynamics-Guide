%
%  The modal parameters are obtained from the system matrices A and C,
%  which are idntified by data-driven stochastic system identification
%  method
%************************************************************************%

clear 
close all
clc;
% % input data
%     fname = input('Data filename to be read ? ','s') ; 
%     signal = load(fname) ;
%     fs = input('sampling frequency ? ') ; 
%     nch = 11; % nch is the number of channels
path='C:\Users\shidada\Desktop\东宝河大桥 SSI\mt\'; %修改到数据所在文件见
% filename='Data_20181215.dat'; 
% filename='Data_20181216.dat'; 
% filename='Data_20181217.dat'; 
filename='Data_20181215.dat'; 
% filename='Data_20181219.dat';
signal =load([path,filename]); 
Vsignal=signal(:,[3,1,5,7,19,17,13,12]);
Hsignal=signal(:,[4,2,6,10,20,18,15,11]);

fs = 50;% sampling frequency 
nch = 8;% the number of channels
% system matrices identification using SSI
Hour=1;
Datasct=Hour*60*60*fs+1:1:(Hour+1)*60*60*fs;

yy=Vsignal(Datasct,1:end);

NFFT=512*8;
for kk=1:size(yy,2)
  fftplot_ff(yy(:,kk),fs,NFFT,0.5,'on');
end



 nbr=20;          % The number of block rows used in block hankel matrices. 
                 % The maximal order that can be estimated is nbr times the number of outputs.
% W = input('weighting matrices? (CVA/PC/UPC)');
 W = 'CVA';
%  W = 'PC';
%  W = 'UPC';
fa=zeros(1,2);
dampa=zeros(1,2);
phia=zeros(nch,2); 
f=zeros(1,2);
damp=zeros(1,2);
phi=zeros(nch,2); 
mm=2;
r=1;
for n=4:30
    [A,C,G,L0,ss] = sto_ac(yy,nbr,n,W);
    [V,D] = eig(A);   % eigenvalue decomposition
    mu = diag(D);
    lambda = log(mu)*fs;
    mode = length(lambda);
  for i=1:mode
    omega(i) = sqrt((real(lambda(i)))^2+(imag(lambda(i)))^2);
    fn(i) = omega(i)/(2*pi);
    dampn(i) = -real(lambda(i))/omega(i);
  end
  phin = real(C*V);

  % sort the modal parameters in ascending order according the frequencies
  [fx,ix] = sort(fn,'ascend');
  fn = fx ;
  for i=1:mode
  dampx(i) = dampn(ix(i));
  phix(:,i) = phin(:,ix(i));  
  end
  dampn = dampx;
  phin = phix;
  
  i=1;
  m=1;
 while i<=mode-1
      if dampn(i)>0 && dampn(i)<1
        if fn(i)==fn(i+1)
         f(m)=fn(i);
         damp(m)=dampn(i);
         phi(:,m)=phin(:,i);
         i=i+2;
         m=m+1;
        else 
         f(m)=fn(i);
         damp(m)=dampn(i);
         phi(:,m)=phin(:,i);
         i=i+1;
         m=m+1;
        end
      else
          i=i+1;
      end
 end
  
  for i=1:mm
  if fa(i)~=0 && dampa(i)~=0 && (phi(:,i)'*phi(:,i))*(phia(:,i)'*phia(:,i))~=0
      nn=length(f);
      for j=1:nn
       if (f(j)-fa(i))/fa(i)<0.01 && (damp(j)-dampa(i))/dampa(i)<0.05 && 1-((phi(:,j)'*phia(:,i))^2/((phi(:,j)'*phi(:,j))*(phia(:,i)'*phia(:,i))))<0.02 
       xf(r)=fa(i);
       xdamp(r)=dampa(i);
       xphi(:,r)=phia(:,i);
       xn(r)=n-1;
       r=r+1;
       figure(1);
       plot(fa(i),n-1,'*b');
       title(' stabilization diagram ');
       ylabel('order ');
       xlabel('frequency（Hz）');
       hold on;
       end
      end
  end
  end
  fa = f;
  mm=length(fa);
  dampa = damp;
  phia = phi;
end

Nfft=1024*1;
% L=length(y);Nfft=2^nextpow2(L);
x = yy(:,5);  % select one channel signal to plot the PSD
[Pxx,ff]=pwelch(x,hanning(Nfft),Nfft/2,Nfft,fs,'onesided');
figure(1);
plot(ff,abs(Pxx)*5E0,'r');  % 5E7 is the scale to amplify the PSD (50E-2)

nsa=input('No. of the stabilization axises to pick ? ');
freq=zeros(nsa);
nsp=length(xf);
Tmp=zeros(nsa,2+nch);
for i=1:nsa
    [xx,yy]=ginput(2) ; 
    nsum=0;
    sumf=0;
    sumdamp=0;
    sumphi=zeros(nch,1);
    for j=1:nsp
      if  xf(j)>xx(1) && xf(j)<xx(2)
          sumf = sumf+xf(j);
          sumdamp = sumdamp+xdamp(j);
          sumphi = sumphi+xphi(:,j);
          nsum=nsum+1;
      end
    end
    freq = sumf/nsum;
    damping = sumdamp/nsum;
    modeshape = sumphi/nsum;
    No_mode = input('Mode_number : ') ; 
    fname = ['mode' num2str(No_mode)]; 
    tmp = [freq  damping  modeshape']' ; 
    Tmp(i,:)= tmp ;
    eval(['save ',fname,'.txt tmp -ascii']) ;  
end
save ModePara.txt Tmp -ascii

% eval(['save ','ModePara','.txt Tmp -ascii']) 

