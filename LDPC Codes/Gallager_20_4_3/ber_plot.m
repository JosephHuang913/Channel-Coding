close; 
clear;

set(gca, 'fontsize', 12)
load ber_LDPC_20_5.log;
i=1:1:11;
semilogy(ber_LDPC_20_5(i,1), ber_LDPC_20_5(i,2), '-ro', 'LineWidth', 1.6, 'MarkerSIze', 8);
hold on;
%grid on;
semilogy(ber_LDPC_20_5(i,1), ber_LDPC_20_5(i,3), '-b*',  'LineWidth', 1.6, 'MarkerSIze', 8);
% semilogy(ber_LDPC_20_5(i,1), ber_LDPC_20_5(i,4), '-.ro');
% semilogy(ber_LDPC_20_5(i,1), ber_LDPC_20_5(i,5), '-.ro');
% semilogy(ber_LDPC_20_5(i,1), ber_LDPC_20_5(i,6), '-.ro');
% semilogy(ber_LDPC_20_5(i,1), ber_LDPC_20_5(i,7), '-.ro');
% semilogy(ber_LDPC_20_5(i,1), ber_LDPC_20_5(i,8), '-.ro');
% semilogy(ber_LDPC_20_5(i,1), ber_LDPC_20_5(i,9), '-.ro');
% semilogy(ber_LDPC_20_5(i,1), ber_LDPC_20_5(i,10), '-.gs');
semilogy(ber_LDPC_20_5(i,1), ber_LDPC_20_5(i,11), '-gs',  'LineWidth', 1.6, 'MarkerSIze', 8);

r=0:1:10;
Pb=0.5.*erfc(sqrt(10.^(r./10)));
semilogy(r, Pb, '-.k',  'LineWidth', 1.6, 'MarkerSIze', 8);

title('BER Performance of Galager (20,5) Regular LDPC Code in AWGN Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');
%axis([0 10 1e-6 1]);

legend('1_s_t Iteration', '2_n_d Iteration', '10_t_h Iteration', 'BPSK, AWGN', 'Parity Check Matrix: (20,15)', 'Row Weight: 4', 'Column Weight: 3', 3);
%print -djpeg100 LDPC_20_5.jpg;
