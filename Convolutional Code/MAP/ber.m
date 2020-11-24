close;
clear;

load ber_ite_CC.log;
i=1:1:11;
semilogy(ber_ite_CC(i,1), ber_ite_CC(i,2), '-rs');
hold on;
semilogy(ber_ite_CC(i,1), ber_ite_CC(i,3), '-bo');
semilogy(ber_ite_CC(i,1), ber_ite_CC(i,4), '-g^');
grid on;
title('BER Performance of Convolutional Code (2,1,2) in AWGN Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');

legend('BCJR Algorithm, 1_s_t Iteration', '2_n_d Iteration', '3_r_d Iteration', 'Generator Polynomial: \{5,7\} in Octal', 'd_f_r_e_e = 5', 3);
%print -djpeg100 convolutional_1_awgn.jpg;
