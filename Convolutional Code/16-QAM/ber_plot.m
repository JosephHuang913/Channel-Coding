close;
clear;

load ber_16QAM.log;
i=1:1:16;
semilogy(ber_16QAM(:,1), ber_16QAM(:,2), '-ro');
hold on;

load ber_HDVA_16QAM.log
i=1:1:11;
semilogy(ber_HDVA_16QAM(:,1), ber_HDVA_16QAM(:,2), '-b*');

grid on;
title('BER Performance of Convolutional Code (2,1,2) with 16-QAM in AWGN Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');

legend('16-QAM', 'Viterbi Algorithm (HDVA)', 'Generator Polynomial: \{5,7\} in Octal', 'd_f_r_e_e = 5', 3);
%print -djpeg100 convolutional_SDVA.jpg;
