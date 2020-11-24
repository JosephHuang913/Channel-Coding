close;
clear;

set(gca, 'fontsize', 12)
load ber.log;
i=1:1:8;
semilogy(ber(:,1), ber(:,3), '-ro', 'LineWidth', 1.6, 'MarkerSIze', 8);
grid on;
hold on;

r=0:1:10;
Pb=0.5.*erfc(sqrt(10.^(r./10)));
semilogy(r, Pb, '-.k', 'LineWidth', 1.6, 'MarkerSIze', 8);

title('BER Performance of Convolutional Code (2, 1, 4) in AWGN Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');
axis([0 10 1e-6 1]);

legend('Viterbi Algorithm (SOVA)', 'BPSK, AWGN', 'Generator Polynomial: \{23,35\} in Octal', 'd_f_r_e_e = 7', 3);
%print -djpeg100 convolutional_SOVA.jpg;
 