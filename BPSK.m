clc;
clear;

%% 1) BPSK 
%% Define transmitted signal (BPSK)
N=8; %Number of bits
X_input= randi([0, 1],1,N);  %Binary signal                    
Tb=1; %Bit duration = 1 msec
X_digit=[]; 
nb=10000;  %Number of points between two symbols (it's used to convert the symbols into continuous digital signal)
for i=1:N
    if X_input(i)==1
       x_temp=ones(1,nb);
    else
        x_temp=zeros(1,nb);
    end
    X_digit=[X_digit x_temp];
end
t_sig = Tb/nb : Tb/nb : N*Tb; %Time vector of continuous digital signal
%Plotting the input message signal
figure();
plot(t_sig,X_digit, 'LineWidth',2,'Color','black'); grid on; xlim([0 Tb*N]); ylim([-0.5 1.5]);
xlabel('Time(sec)'); ylabel('Amplitude(volt)'); title('Input Message Signal');

%% BPSK Modulation (BPSK)
Ac=sqrt(2/Tb);   %Carrier Signal Amplitude
fc = 4*(1/Tb);  %Carrier Signal Frequency
phi_1 = 0;       %Phase Shift for bit '1'
phi_0 = pi;      %Phase Shift for bit '0'
t_cycle = Tb/nb : Tb/nb : Tb;  %Time of one symbol to be used to calculate the carrier signal to be multiplied by the input signal
X_BPSK=[];
x_con_mod = [];
%Multiplying the input message signal by the carrier based on the message symbol
for i=1:N
    if X_input(i)==1
        x_temp_mod = Ac*cos(2*pi*fc*t_cycle + phi_1);
        x_con_mod = [x_con_mod, 1];
    else
        x_temp_mod = Ac*cos(2*pi*fc*t_cycle + phi_0);
        x_con_mod = [x_con_mod, -1];
    end
    X_BPSK=[X_BPSK x_temp_mod];
end
t_mod = Tb/nb : Tb/nb : N*Tb; %Time of the modulated signal
%Plotting the BPSK Signal
figure();
plot(t_mod,X_BPSK,'Color','blue'); grid on;
xlabel('Time(sec)'); ylabel('Amplitude(volt)'); title('BPSK Modulated Signal');

%% Noise in the Communication Channel (BPSK)
Y = awgn(X_BPSK,0.000001,'measured');  %Adds white Gaussian noise to the Modulated signal
%Plotting the Noisy Signal
figure();
plot(t_mod,Y,'Color','red'); grid on;
xlabel('Time(sec)'); ylabel('Amplitude(volt)'); title('Noisy Signal');

%% BPSK Demodulation (BPSK)
X_demod=[];
x_con_dem = [];
for i=nb:nb:length(Y)
  t1_dem = Tb/nb:Tb/nb:Tb;                %Time of one symbol
  Ac_dem = sqrt(2/Tb);                    %Carrier Amplitude
  C_dem = Ac_dem*cos(2*pi*fc*t1_dem);     %Carrier that will be used in Coherent Detection 
  %Correlator 
  y_temp_dem = C_dem.*Y((i-(nb-1)):i);    %Multiply with the carrier
  y_corr=trapz(t1_dem,y_temp_dem);             %Integrate over the time period of the symbol
  x_con_dem = [x_con_dem y_corr];
  
  %Decision Making Device
  if(y_corr>0) %If the value > threshold (0) 
    S=1;       %Then symbol = 1
  else         %If the value < threshold (0) 
    S=0;       %Then symbol = 0
  end
  X_demod=[X_demod S];
end

%Convert the symbols into continuous digital signal
X_dem_sig = [];
for i=1:length(X_demod)
    if X_demod(i)==1
       x_temp_dem=ones(1,nb);
    else
        x_temp_dem=zeros(1,nb);
    end
    X_dem_sig=[X_dem_sig x_temp_dem];
end
t_sig_dem = Tb/nb : Tb/nb : length(X_demod)*Tb; %Time vector of continuous digital signal
figure();
plot(t_sig_dem,X_dem_sig, 'LineWidth',2,'Color','blue'); grid on; xlim([0 Tb*length(X_demod)]); ylim([-0.5 1.5]);
xlabel('Time(sec)'); ylabel('Amplitude(volt)'); title('BPSK Demodulated Signal');


%% Graphs in Subplot (BPSK)

%Plotting the input message signal
figure();
subplot(4,1,1);
plot(t_sig,X_digit, 'LineWidth',2,'Color','black'); grid on; xlim([0 Tb*N]); ylim([-0.5 1.5]);
xlabel('Time(sec)'); ylabel('Amplitude(volt)'); title('Input Message Signal');

%Plotting the BPSK Signal
subplot(4,1,2);
plot(t_mod,X_BPSK,'Color','blue'); grid on;
xlabel('Time(sec)'); ylabel('Amplitude(volt)'); title('BPSK Modulated Signal');

%Plotting the Noisy Signal
subplot(4,1,3);
plot(t_mod,Y,'Color','red'); grid on;
xlabel('Time(sec)'); ylabel('Amplitude(volt)'); title('Noisy Signal');

%Plotting the Demodulated BPSK Signal
subplot(4,1,4);
plot(t_sig_dem,X_dem_sig, 'LineWidth',2,'Color','blue'); grid on; xlim([0 Tb*length(X_demod)]); ylim([-0.5 1.5]);
xlabel('Time(sec)'); ylabel('Amplitude(volt)'); title('BPSK Demodulated Signal');


%% Constellation Diagram (BPSK)
figure();
scatter(x_con_dem, zeros(1,length(x_con_dem)),40, 'MarkerEdgeColor', [0 .5 .5],'MarkerFaceColor',[0 .7 .7], 'LineWidth',1.5); grid on;
hold on
scatter(x_con_mod, zeros(1,length(x_con_mod)), 70, 'red', '+', 'LineWidth',2);
hold off
ylim([-1 1]); xlim([-2 2]);
xlabel('Phi1(t)'); title('Constellation Diagram'); legend('Noisy BPSK Signal','Transmitted BPSK Signal')


%% PSD (BPSK)
psd_BPSK = fft(X_BPSK);  %FFT for the Modulated Signal
psds_BPSK = fftshift(psd_BPSK);  %FFT Shift
psd_neg = flip(psds_BPSK);   %The negative side frequencies
psd = [psd_neg psds_BPSK];   %Combining the frequencies of both sides
f = linspace(-2*fc, 2*fc, length(psd));  %The frequency Vector
%Plotting PSD
figure();
plot(f, abs(psd), 'LineWidth',2,'Color','black'); grid on;
xlabel('Frequency(Hz)'); ylabel('PSD(Magnitude)'); title('PSD of the Transmitted BPSK Signal');


%% BER (BPSK) 
Bit_Error_Rate = biterr(X_demod, X_input)


