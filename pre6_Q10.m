clc,clear all,close all
[x,fs]=audioread('sound.wav');
x=x'
h1=fir1(50,0.2,hamming(51));
l_x=length(x);

dft_x=dftadd(x,h1,l_x);
t=[1:l_x]/fs;
[fr_x,w_x]=freqz(x);
[fr_dftx,w_o]=freqz(dft_x);
fr_x=[flipud(fr_x);fr_x];
fr_dftx=[flipud(fr_dftx);fr_dftx];
freq=linspace(-1,1,length(fr_x))
figure,subplot(3,2,1),plot(t,x);
title('Input Signal'),xlabel('Time(s)'),ylabel('Magnitude');
subplot(3,2,2),plot(freq,abs(fr_x)/l_x);
title('Magnitude Spectrum of Input Signal'),xlabel('Frequency'),ylabel('Magnitude');
subplot(3,2,3),plot(t,dft_x);
title('Output Signal=Filtered Input Signal Using Overlap Add Method with DFT')
xlabel('Time(s)'),ylabel('Magnitude');
subplot(3,2,4),plot(freq,abs(fr_dftx)/l_x)
title('Magnitude Spectrum of Output Signal(DFT)'),xlabel('Frequency'),ylabel('Magnitude'); 
conv_x=convadd(x,h1,l_x);
t=[1:l_x]/fs;
[fr_x,w_x]=freqz(x);
[fr_convx,w_o]=freqz(conv_x);
fr_x=[flipud(fr_x);fr_x];
fr_convx=[flipud(fr_convx);fr_convx];
freq=linspace(-1,1,length(fr_x));
subplot(3,2,5),plot(t,x);
title('Output Signal=Filtered Input Signal Using Overlap Add Method with Circular Convolution.');
xlabel('Time(s)'),ylabel('Magnitude');
subplot(3,2,6),plot(freq,abs(fr_convx)/l_x)
title('Magnitude Spectrum of Output Signal(Cir)'),xlabel('Frequency'),ylabel('Magnitude');
gtext('Ziya Tuna Bölükbaşı')
 function y=convadd(x,h,L)
  P=length(h);
  y=convfilt(x(1:L),h,L);
    for i=1:length(x)/L-1
        k=convfilt(x(i*L+1:(i+1)*L),h,L);
        y(end-P+1:end)=y(end-P+1:end)+k(1:P);
        y=cat(2,y,k(P+1:end));
    end
  end
 function y=dftadd(x,h,L)
  P=length(h);
  y=dftfilt(x(1:L),h,L);
    for i=1:length(x)/L-1
        k=dftfilt(x(i*L+1:(i+1)*L),h,L+P-1);
        y(end-P+1:end)=y(end-P+1:end)+k(1:P);
        y=cat(2,y,k(P+1:end));
    end
  end


function y=dftfilt(x,h,N)
X=fft(x,N);
H=fft(h,N);
Y=X.*H;
y=ifft(Y,N);
end

function  y  = convfilt( x, h, N )
x = [x , zeros(1,N-length(x))];
r = [h , zeros(1,N-length(h))];
c = [h(1) , fliplr(r(2:end))];
H = toeplitz(c, r);
y = x*H;
end