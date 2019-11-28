% 2. Rice fading channel K = 0, K = 10, K = 30
clear
clc

N = 100000; % number of symbols
Es_N0 = 0: 1: 30; %SNR

% generate QPSK symbols
% 1 + i, 1 - i, -1 + i, -1 - i
Tx = (2 *(rand(1, N) > 0.5) - 1) + j * (2 * (rand(1, N) > 0.5) - 1);
Tx_normal = Tx/sqrt(2); % normalize transmitted signal

n = (randn(1, N) + j * randn(1, N))/sqrt(2); % noise
K = [0 10 30];

% fading channel
h0 = (randn(1, N) + j * randn(1, N))/sqrt(2) + (1 + j)*sqrt(K(1)/2); % K = 0
h10 = (randn(1, N) + j * randn(1, N))/sqrt(2) + (1 + j)*sqrt(K(2)/2); % K = 10
h30 = (randn(1, N) + j * randn(1, N))/sqrt(2) + (1 + j)*sqrt(K(3)/2); % K = 30

Rx_output = zeros(3, length(Es_N0));
P_s = zeros(3, length(Es_N0)); % probability of symbol error

for i = 1: length(Es_N0)
    Rx_0 = abs(h0) .* Tx_normal * 10^(Es_N0(i)/20) + n; % received signal with K = 0
    Rx_10 = (1/sqrt(K(2))) * abs(h10) .* Tx_normal * 10^(Es_N0(i)/20) + n; % with
    Rx_30 = (1/sqrt(K(3))) * abs(h30) .* Tx_normal * 10^(Es_N0(i)/20) + n; % with
    
    Rx_re_0 = real(Rx_0);
    Rx_im_0 = imag(Rx_0);
    
    Rx_re_10 = real(Rx_10);
    Rx_im_10 = imag(Rx_10);
    
    Rx_re_30 = real(Rx_30);
    Rx_im_30 = imag(Rx_30);
    
    Rx_output(1, find(Rx_re_0 > 0 & Rx_im_0 > 0)) = 1 + j;
    Rx_output(1, find(Rx_re_0 > 0 & Rx_im_0 < 0)) = 1 - j;
    Rx_output(1, find(Rx_re_0 < 0 & Rx_im_0 > 0)) = -1 + j;
    Rx_output(1, find(Rx_re_0 < 0 & Rx_im_0 < 0)) = -1 - j;
    
    Rx_output(2, find(Rx_re_10 > 0 & Rx_im_10 > 0)) = 1 + j;
    Rx_output(2, find(Rx_re_10 > 0 & Rx_im_10 < 0)) = 1 - j;
    Rx_output(2, find(Rx_re_10 < 0 & Rx_im_10 > 0)) = -1 + j;
    Rx_output(2, find(Rx_re_10 < 0 & Rx_im_10 < 0)) = -1 - j;
    
    Rx_output(3, find(Rx_re_30 > 0 & Rx_im_30 > 0)) = 1 + j;
    Rx_output(3, find(Rx_re_30 > 0 & Rx_im_30 < 0)) = 1 - j;
    Rx_output(3, find(Rx_re_30 < 0 & Rx_im_30 > 0)) = -1 + j;
    Rx_output(3, find(Rx_re_30 < 0 & Rx_im_30 < 0)) = -1 - j;
    
    count_error_0 = size(find(Tx - Rx_output(1, :)), 2);
    count_error_10 = size(find(Tx - Rx_output(2, :)), 2);
    count_error_30 = size(find(Tx - Rx_output(3, :)), 2);
    
    P_s(1, i) = count_error_0/N;
    P_s(2, i) = count_error_10/N;
    P_s(3, i) = count_error_30/N;
end
figure(3)

semilogy(Es_N0, P_s(1, :), '-b+');
hold on
semilogy(Es_N0, P_s(2, :), '-bo');
semilogy(Es_N0, P_s(3, :), '-b<');
grid on
xlabel('E_s/N_0(dB)')
ylabel('Probability of symbol error (%)')
hold off
legend('Rice K = 0', 'Rice K = 10', 'Rice K = 30')
axis([1 30 10^-3 1])
