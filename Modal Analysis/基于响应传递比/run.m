A1=load('C:\Users\shidada\Desktop\基于响应传递比\data.txt');   %%% 读取实测数据存放在A1矩阵中   
Fs=1000;%%%采样频率为1000Hz
nfft=1024;%%%fft的分段长度
[m,n]=size(A1);%%%得到矩阵A1的行数m和列数n
for i=1:n
for j=1:n
    [Pff2,f_disc]=cpsd(A1(:,i),A1(:,j),hamming(nfft),nfft/2,nfft,Fs);
    G0(i,j,:)=Pff2;
end 
end
G0_real=G0;  %%%计算各个通道的互谱值并存储在三维数组G0_real中
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%找到残差函数，残差函数的列式见Engineering Structures 102(2015) 108-119
%%%的文章公式（15），通过残差函数极点对应的频率值即为频率值。
 k=0;
 for i=1:n                                                           
 for j=1:n
 if j~=i
     G0_trans1=(G0_real(j,j,:)./G0_real(i,j,:));
     G0_trans2=(G0_real(j,i,:)./G0_real(i,i,:));
     CSD_change=G0_trans1-G0_trans2;
     k=k+1;
     PSDT_inverse(:,k)=(1./abs(CSD_change));%%%把前面的表达转换成k列了
 end
 end
 end
    
PSDT_inverse0=sum((PSDT_inverse),2); %%%对每一行进行求和,生成一个数组PSDT_inverse0  
%%%%绘制图PSDT_inverse0%%%%绘图有问题需要重新绘制
figure(1)
plot(f_disc,PSDT_inverse0);%%%拾取该图中峰值对应的频率即可