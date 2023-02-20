%% BER (PBSK)
clc;
clear;
Eb_N0_dB = 0:0.1:10;
Eb_N0 = 10.^(Eb_N0_dB/10);
x= sqrt(Eb_N0);
BER = 1/2.*erfc(x);
figure();
semilogy(Eb_N0_dB,BER); grid on; ylabel('BER'); xlabel('Eb/N0 (dB)'); title('Theoretical BER for BPSK');

%% BER (QPSK)
clc;
clear;
Eb_N0_dB = 0:0.1:10;
Eb_N0 = 10.^(Eb_N0_dB/10);
x= sqrt(Eb_N0);
BER = erfc(x);
figure();
semilogy(Eb_N0_dB,BER); grid on; ylabel('BER'); xlabel('Eb/N0 (dB)'); title('Theoretical BER for QPSK');