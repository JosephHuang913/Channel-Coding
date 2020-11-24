close;
clear;

set(gca, 'fontsize', 14)
load ber_tail_biting.log;
semilogy(ber_tail_biting(:,1), ber_tail_biting(:,3), '-ro', 'LineWidth', 2.5, 'MarkerSIze', 10);
grid on;
hold on;

load ber_non_tail_biting.log;
semilogy(ber_non_tail_biting(:,1), ber_non_tail_biting(:,3), '-b*', 'LineWidth', 2.5, 'MarkerSIze', 12);

r=0:1:10;
Pb=0.5.*erfc(sqrt(10.^(r./10)));
semilogy(r, Pb, '-.k', 'LineWidth', 2.5, 'MarkerSIze', 12);

title('BER Performance of Convolutional Code (2, 1, 6) in AWGN Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');
axis([0 10 1e-6 1]);

legend('Viterbi Algorithm (SOVA), Tail-Biting', 'Viterbi Algorithm (SOVA)', 'BPSK, AWGN', 'Generator Polynomial: \{171,133\} in Octal', 'd_f_r_e_e = 10', 1);
%print -djpeg100 convolutional_171_133_SOVA.jpg;
