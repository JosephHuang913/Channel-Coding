close;
clear;

ber1 = fopen('ber.log', 'r');
ber2 = fscanf(ber1, '%e');
i=1:1:11;
Eb_No(i) = ber2(2*i-1);
err_rate(i) = ber2(2*i);
semilogy(Eb_No, err_rate, '-ro');
fclose(ber1);

grid on;
title('BER Performance of Convolutional Code (2,1,2) in AWGN Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');

legend('Viterbi Algorithm (Hard Decision), 1-path AWGN', 'Generator Polynomial: \{5,7\} in Octal', 'd_f_r_e_e = 5', 3);
%print -djpeg100 convolutional_HDVA.jpg;
