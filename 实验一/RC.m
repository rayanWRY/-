R = 5.1 * 10^3;
C = 1983.3333333 * 10^(-12);
f0 = 1 / (2 * pi * R * C);
f = [0:10:1000*10^3];
w = 2 * pi * f;
H = @(w) 1 ./ (1 + 1i * R * C * w);




figure(1);
sys=tf([1], [R*C 1]);
bode(sys);

figure(2);
subplot(2, 1, 1);
hold on;
grid on;
%xlabel("f/10^xHz");
xlabel("f/Hz");
ylabel("|H(jw)|");
title("RC幅频特性曲线");
%plot(log10(f), abs(H(w)), '-b');
%plot(log10(f0), abs(H(2*pi*f0)), 'r*');
%plot(log10(156556), 0.1, 'r*');
plot(f, abs(H(w)), '-b');
plot(f0, abs(H(2*pi*f0)), 'r*');
plot(156556, 0.1, 'r*');


subplot(2, 1, 2);
hold on;
grid on;
xlabel("f/10^xHz");
ylabel("|H(jw)|/dB");
title("RC幅频特性曲线(dB)");
plot(log10(f), 20 * log10(abs(H(w))), '-b');
plot(log10(f0), 20 * log10(abs(H(2*pi*f0))), 'r*');
plot(log10(156556), 20 * log10(0.1), 'r*');
hold off;

N=10000;
fHz=[5 10 20 50 200 300]*10^3;
fs=100*fHz;
w=2*pi*fHz;
for i=1:length(w)
    t=1/fs(i)*[0:N-1];
    u=1.08*sin(w(i)*t);
    %u=2*square(w(i)*t,50);
    figure;    
    y=lsim(sys,u,t);
    p=plot(t,u,t,y);
    legend('original signal','pass filter signal');
    title(strcat('频率',num2str(fHz(i)/10^3),'kHz','时域'));
    xlabel("t/s");
    ylabel("Amplitude");
    axis([0 3/fHz(i) -3 3]);
    
    f=fs(i)/N*[0:N/2];
    U=fft(u,N);
    Y=fft(y,N);
    UA=abs(U);
    YA=abs(Y); 
    UA=UA(1:N/2+1);
    YA=YA(1:N/2+1);
    figure;
    plot(f,UA,f,YA);
    legend('original signal','pass filter signal');
    title(strcat('频率',num2str(fHz(i)/10^3),'kHz','频域'));
    xlabel("f/Hz");
    ylabel("Amplitude");
    axis([0 20*fHz(i) 0 inf]);
    
    PU=1/N*UA.^2;
    PY=1/N*YA.^2;
    figure;
    plot(f,PU,f,PY);
    legend('original signal','pass filter signal');
    title(strcat('频率',num2str(fHz(i)/10^3),'kHz','功率谱'));
    xlabel("f/Hz");
    ylabel("W/Hz");
    axis([0 20*fHz(i) 0 inf]);
    
    [Ru,maxlagsu]=xcorr(u,'unbiased');
    [Ry,maxlagsy]=xcorr(y,'unbiased');
    figure;
    plot(maxlagsu/fs(i),Ru/max(Ru),maxlagsy/fs(i),Ry/max(Ry));
    legend('original signal','pass filter signal');
    title(strcat('频率',num2str(fHz(i)/10^3),'kHz','自相关函数'));
    xlabel("t/s");
    ylabel("R(t)");
    axis([ -10/fHz(i) 10/fHz(i) -1.2 1.2]);
end