fc=20000; %载波频率
fs=40000; %采样速率
k=2;
code_size=15*round(k*fs/fc); %信息码元长度
t0=5.5; %信号长度
Ns=256; %采样点个数
fd=125; %符号速率
ts=1/fs; %采样周期
M=64; %码元个数
ti=1/fd; %码元间隔
N=ti/ts;
t=0:ts:t0;
select=menu('调制方式','2ASK','2FSK','2PSK','4ASK','4FSK','4PSK');
switch select
case 1, %                             2ASK signal
    x=randi(1,M);
    m=sin(2*pi*fc*t);
    y=ones(1,M*N);
    for i=1:M
        for j=1:N
            y((i-1)*N+j)=x(i)*m(j);
        end
    end
    T=zeros(6,50);
    T(1,1:50)=1;    
case 2, %                           2FSK signal
    x=randi(1,M);
    m1=sin(2*pi*fc*t);
    m2=sin(2*pi*2*fc*t);
    y=zeros(1,M*N); 
    for i=1:M
        if x(i)==1;
        for j=1:N;
        y((i-1)*N+j)=x(i)*m1(j);
        end
        elseif x(i)==0;
            for j=1:N;
                y((i-1)*N+j)=(1-x(i))*m2(j);
            end
        end
    end
    T=zeros(6,50);
    T(2,1:50)=1;
case 3, %                       2PSK signal,
    x=randi(1,M);
    m1=sin(2*pi*fc*t);
    m2=sin(2*pi*fc*t+pi);
    y=zeros(1,M*N);
    for i=1:M
        if x(i)==1;
        for j=1:N;
            y((i-1)*N+j)=x(i)*m1(j);
        end
        elseif x(i)==0;
            for j=1:N;
                y((i-1)*N+j)=(1-x(i))*m2(j);
            end
        end
    end
    T=zeros(6,50);
    T(3,1:50)=1;
case 4, % 4ASK signal 
    x1=randi(1,M);
    x2=randi(1,M);
    x3=randi(1,M);
    x=x1+x2+x3;
    m=sin(2*pi*fc*t);
    y=ones(1,M*N);
    for i=1:M
        if x(i)==0;
        for j=1:N 
            y((i-1)*N+j)=x(i)*m(j);
        end
        elseif x(i)==1;
            for j=1:N
                y((i-1)*N+j)=x(i)*m(j);
            end
        elseif x(i)==2;
            for j=1:N
                y((i-1)*N+j)=x(i)*m(j);
            end
        elseif x(i)==3;
            for j=1:N
                y((i-1)*N+j)=x(i)*m(j);
            end
        end
    end
    T=zeros(6,50);
    T(4,1:50)=1;
case 5, %           4FSK signal
    x1=randi(1,M);
    x2=randi(1,M);
    x3=randi(1,M);
    x=x1+x2+x3;
    m1=sin(2*pi*fc*t);
    m2=sin(2*pi*2*fc*t);
    m3=sin(2*pi*3*fc*t);
    m4=sin(2*pi*4*fc*t);
    y=zeros(1,M*N);
    for i=1:M
        if x(i)==0;
        for j=1:N
            y((i-1)*N+j)=(1-x(i))*m1(j);
        end 
        elseif x(i)==1;
        for j=1:N
            y((i-1)*N+j)=x(i)*m2(j);
        end
        elseif x(i)==2;
        for j=1:N
            y((i-1)*N+j)=(x(i)-1)*m3(j);
        end 
        elseif x(i)==3;
        for j=1:N
            y((i-1)*N+j)=(x(i)-2)*m4(j);
        end
        end
    end
    T=zeros(6,50);
    T(5,1:50)=1;
case 6, %4PSK signal
    x1=randi(1,M);
    x2=randi(1,M);
    x3=randi(1,M);
    x=x1+x2+x3;
    m1=sin(2*pi*fc*t);
    m2=sin(2*pi*fc*t+pi/2);
    m3=sin(2*pi*fc*t+pi);
    m4=sin(2*pi*fc*t+3*pi/2);
    y=zeros(1,M*N);
    for i=1:M
        if x(i)==0;
        for j=1:N
            y((i-1)*N+j)=(1-x(i))*m1(j);
        end
        elseif x(i)==1;
            for j=1:N
                y((i-1)*N+j)=x(i)*m2(j);
            end
        elseif x(i)==2;
            for j=1:N
                y((i-1)*N+j)=(x(i)-1)*m3(j);
            end
        elseif x(i)==3;
            for j=1:N
                y((i-1)*N+j)=(x(i)-2)*m4(j);
            end
        end
    end
    T=zeros(6,50);
    T(6,1:50)=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNR=20; %定义信噪比，单位 DB 
sigpow=mean(abs(y).^2); %power of input signal
noisefac=10^(-SNR/10);
noise=randn(1,size(y,2));
noise=noise*(sqrt(sigpow*noisefac)/sqrt(mean(noise.^2))); %产生所需的高斯噪声
ynoise=noise+y; %加噪后的信号
for n=1:1:50
     m=n*Ns;
     x=(n-1)*Ns;
for i=x+1:m; %提取信号段
    y0(i)=ynoise(i);
end
Y=fft(y0); %调制信号的傅立叶变换
y1=hilbert(y0); %实信号的解析式
z=abs(y0); %实信号的瞬时幅度
phase=angle(y1); %实信号的瞬时相位
add=0; %求 Rmax
for i=x+1:m;
    add=add+z(i);
end
ma=add/Ns; %瞬时幅度的平均值
y2=z./ma ; %幅度比,即为文献中的 an(i)
y3=y2-1; %归一化瞬时幅度
y4=max(abs(y3));
y2n=y3./y4; % 即为文献中的 acn(i)
s=fft(y2n);
R=abs(s);
Rmax=max((R)/Ns).^2; %零中心归一化瞬时幅度的谱密度的最大值
Xcn=0;
Ycn=0;
for i=x+1:m;
     Xcn=Xcn+y2n(i).*y2n(i);
     Ycn=Ycn+abs(y2n(i));
end
Xcnav=Xcn/Ns;
Ycnav=(Ycn/Ns).*(Ycn/Ns);
deltaaa=sqrt(Xcnav-Ycnav); %零中心归一化瞬时幅度绝对值得标准偏差
if phase(2+x)-phase(1+x)>pi; %修正相位序列
     Ck(1+x)=-2*pi;
elseif phase(1+x)-phase(2+x)>pi;
    Ck(1+x)=2*pi;
else Ck(1+x)=0;
end
for i=x+2:m-1;
     if phase(i+1)-phase(i)>pi;
         Ck(i)=Ck(i-1)-2*pi;
     elseif phase(i)-phase(i+1)>pi
        Ck(i)=Ck(i-1)+2*pi;
     else
        Ck(i)=Ck(i-1);
     end
end
if -phase(m)>pi;
    Ck(m)=Ck(m-1)-2*pi;
elseif phase(m)>pi;
     Ck(m)=Ck(m-1)+2*pi;
else Ck(m)=Ck(m-1);
end
phase1=phase+Ck; %去相位卷叠后的相位序列
phasenl=phase1-2*pi*fc*i/fs; %非线性相位
 at=1; %判决门限电平
 a=0; %求取零中心非弱信号段瞬时相位非线性分量绝对值的标准偏差和零中心非弱信号段瞬时相位非线性分量的标准偏差
 b=0;
 d=0;
 c=0;
for i=x+1:m;
     if y2(i)>at
         c=c+1;
         phasesquare(i)=phasenl(i).*phasenl(i);
         a=a+phasesquare(i);
         phaseabs(i)=abs(phasenl(i));
         b=b+phaseabs(i);
         d=d+phasenl(i)
     end
end
a1=a/c;
b1=(b/c).*(b/c);
d1=(d/c).*(d/c);
deltaap=sqrt(a1-b1); %零中心非弱信号段瞬时相位非线性分量绝对值的标准偏差
deltadp=sqrt(a1-d1); %零中心非弱信号段瞬时相位非线性分量的标准偏差
freqN(i)=phase1(i)-phase1(i-1);
for i=x+1:m; 
     if i>at;
         c=c+1;
         freqNsquare(i)=freqN(i).*freqN(i);
         a=a+freqNsquare(i);
         b=b+freqN(i);
     end
end
a1=a/c;
b1=(b/c)^2;
deltaaf=sqrt(a1-b1); %零中心归一化非弱信号段瞬时频率绝对值得标准偏差
end
if (Rmax>80)
    if(deltaaa>0.2)
        type=menu('输入信号是','2PSK 信号');
    else
        type=menu('输入信号是','2ASK 信号');
    end
else
    if (deltaaa<0.23)
        type=menu('输入信号是','4ASK 信号');
    else
        if(deltaap>50)
            type=menu('输入信号是','4FSK 信号');
        else 
            if(Rmax>80)
                type=menu('输入信号是','2FSK 信号');
            else
                type=menu('输入信号是','4PSK 信号');
            end
        end
    end
end