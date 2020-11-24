close;
clear;

set(gca, 'fontsize', 12);
load ber_Galager_20_4_3.log;
i=1:1:11;
semilogy(0,1,'w');
hold on;
grid on;
semilogy(ber_Galager_20_4_3(i,1), ber_Galager_20_4_3(i,2), '-rv', 'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(ber_Galager_20_4_3(i,1), ber_Galager_20_4_3(i,3), '-y*', 'LineWidth', 2.2, 'MarkerSIze', 12);
semilogy(ber_Galager_20_4_3(i,1), ber_Galager_20_4_3(i,6), '-g+', 'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(ber_Galager_20_4_3(i,1), ber_Galager_20_4_3(i,11), '-cs', 'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(ber_Galager_20_4_3(i,1), ber_Galager_20_4_3(i,21), '-m^', 'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(ber_Galager_20_4_3(i,1), ber_Galager_20_4_3(i,51), '-bo', 'LineWidth', 2.0, 'MarkerSIze', 10);
r=0:1:10;
Pb=0.5.*erfc(sqrt(10.^(r./10)));
semilogy(r, Pb, '-.k',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(0,1,'w');
semilogy(0,1,'w');
semilogy(0,1,'w');
semilogy(0,1,'w');
%axis([0 5 1e-4 1]);

title('BER Performance of Gallager (20, 7) Regular LDPC Code in AWGN Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');

legend('Sum Product Algorithm', '1_s_t Iteration', '2_n_d Iteration', '5_t_h Iteration', '10_t_h Iteration', '20_t_h Iteration', '50_t_h Iteration', 'BPSK, AWGN', 'Parity Check Matrix: 15x20', 'Row Weight: 4', 'Column Weight: 3', 3);
%print -djpeg100 Galager_awgn.jpg;
