close;
clear;

load ber.log;

i=1:1:11;
semilogy(ber(:,1), ber(:,3), '-ro');

grid on;
title('BER Performance of Convolutional Code (2,1,2) in AWGN Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');

legend('Viterbi Algorithm (SOVA), 1-path AWGN', 'Generator Polynomial: \{5,7\} in Octal', 'd_f_r_e_e = 5', 3);
%print -djpeg100 convolutional_SOVA.jpg;
