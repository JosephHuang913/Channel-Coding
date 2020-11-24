close;
clear;

ber1 = fopen('ber_BCJR.log', 'r');
ber2 = fscanf(ber1, '%e');
i=1:1:11;
Eb_No(i) = ber2(2*i-1);
err_rate(i) = ber2(2*i);
semilogy(Eb_No, err_rate, '-ro');
fclose(ber1);
hold on;
grid on;

clear;
ber1 = fopen('ber_VA_S.log', 'r');
ber2 = fscanf(ber1, '%e');
i=1:1:11;
Eb_No(i) = ber2(2*i-1);
err_rate(i) = ber2(2*i);
semilogy(Eb_No, err_rate, '-b^');
fclose(ber1);

clear;
ber1 = fopen('ber_VA_H.log', 'r');
ber2 = fscanf(ber1, '%e');
i=1:1:11;
Eb_No(i) = ber2(2*i-1);
err_rate(i) = ber2(2*i);
semilogy(Eb_No, err_rate, '-k.');
fclose(ber1);

clear;
load ber_SOVA.log;
i=1:1:11;
semilogy(ber_SOVA(i,1), ber_SOVA(i,2), '-.gs');
%semilogy(ber_SOVA(i,1), ber_SOVA(i,3), '-.gs');

title('BER Performance of Convolutional Code (2,1,2) in AWGN Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');

legend('BCJR Algorithm', 'Viterbi Algorithm (Soft Decision)', 'Viterbi Algorithm (Hard Decision)', 'SOVA', 'Generator Polynomial: \{5,7\} in Octal', 'd_f_r_e_e = 5', 3);
%print -djpeg100 convolutional_1_awgn.jpg;
