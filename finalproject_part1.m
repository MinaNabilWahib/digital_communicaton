close all
clc
clear all
M = 16;                     % Size of signal constellation
k = log2(M);                % Number of bits per symbol
NumSNR= 10 ^(6);
SNR=0:2:30; %range of SNR
RandomSignal=randi([0,1],NumSNR,1); % generate stream of random bits
numSamplesPerSymbol = 1;  


%% modulations

OOK_mod = RandomSignal;% OOK modulation
OOK_mod_fn = (pammod(RandomSignal,2)+1)./2; % OOK modulation using fn 

PRK_mod = 2.*RandomSignal - 1;%PRK modulation
PRK_mod_fn = pammod(RandomSignal,2);%PRK modulation using fn 

FSk_H_i = find(RandomSignal==1);
FSk_L_i = find(RandomSignal==0);
FSK_modv(FSk_H_i)=1j;%FSK modulation
FSK_modv(FSk_L_i)=1;%FSK modulation
FSK_mod=FSK_modv.';%transpose 
FSK_mod_fn = sqrt(pskmod(RandomSignal,2));

% qam modulation 
b_data_block = reshape(RandomSignal,length(RandomSignal)/k,k);   % Reshape data into binary k-tuples, k = log2(M)
block_data_symbol = bi2de(b_data_block);                 % Convert to integers

QAM_mod_fn = qammod(block_data_symbol,M,'bin');         % Binary coding, phase offset = 0



%% calculating power 
PTX_OOK= mean(OOK_mod.^(2)); %ook power 
PTX_OOK_fn= mean(OOK_mod_fn.^(2)); %ook power 

PTX_PRK= mean(PRK_mod.^(2)); %prk power
PTX_PRK_fn= mean(PRK_mod_fn.^(2)); %prk power

PTX_FSK= mean(abs(FSK_mod.^(2))); %fsk power
PTX_FSK_fn= mean(abs(FSK_mod_fn.^(2))); %fsk power



%% adding noise , demodulation and calculating bit error 
for n=1:length(SNR)
    
     snr_i=10^(SNR(n)/10);% converting from db to linear 
     %calculating noise 
     noise_OOK=sqrt(PTX_OOK/(2 * snr_i) *( randn(size(OOK_mod)) + 1j * randn(size(OOK_mod)) ) ) ; %OOK noise
     noise_OOK_fn=sqrt(PTX_OOK_fn/(2 * snr_i) *( randn(size(OOK_mod_fn)) + 1j * randn(size(OOK_mod_fn)) ) ) ; %OOK noise
     
     
     noise_PRK=sqrt(PTX_PRK/(2 * snr_i) *( randn(size(PRK_mod)) + 1j * randn(size(PRK_mod)) ) ) ; %PRK noise
     noise_PRK_fn=sqrt(PTX_PRK_fn/(2 * snr_i) *( randn(size(PRK_mod_fn)) + 1j * randn(size(PRK_mod_fn)) ) ) ; %PRK noise
     
     
     noise_FSK=sqrt(PTX_FSK/(2 * snr_i) *( randn(size(FSK_mod)) + 1j * randn(size(FSK_mod)) ) ) ; %FSK noise
     noise_FSK_fn=sqrt(PTX_FSK_fn/(2 * snr_i) *( randn(size(FSK_mod_fn)) + 1j * randn(size(FSK_mod_fn)) ) ) ; %FSK noise
     
     %recieved signals 
     RX_OOK = OOK_mod + noise_OOK ; % OOK recieved signal
     RX_OOK_fn = OOK_mod_fn + noise_OOK_fn ; % OOK recieved signal
     

     RX_PRK = PRK_mod +noise_PRK ;  % PRK recieved signal
     RX_PRK_fn = PRK_mod_fn + noise_PRK_fn;

     RX_FSK = FSK_mod + noise_FSK ; % FSK recieved signal
     RX_FSK_fn = FSK_mod_fn + noise_FSK_fn;
     
     snr_i_qam = SNR(n) + 10*log10(k) - 10*log10(numSamplesPerSymbol);
     RX_QAM_fn = awgn(QAM_mod_fn,snr_i_qam,'measured');% qam  recieved signal
     
     %demodulation 
     %OOK_demodulation 
     RX_OOK_abs=abs(RX_OOK);
     RX_OOK_H_i = find(RX_OOK_abs>0.5);
     RX_OOK(RX_OOK_H_i) = 1;
     RX_OOK_L_i = find(RX_OOK_abs<0.5);
     RX_OOK(RX_OOK_L_i) = 0;
     
     RX_OOK_fn=pamdemod((2.*RX_OOK_fn)-1,2);
     
     %PRK demodulation
     RX_PRK_H_i = find(real(RX_PRK)>0);
     RX_PRK(RX_PRK_H_i) = 1;
     RX_PRK_L_i = find(real(RX_PRK)<0);
     RX_PRK(RX_PRK_L_i) = 0;
     
     RX_PRK_fn=pamdemod(RX_PRK_fn,2);
     
     %FSK demodulation
     RX_FSK_R_abs=abs(real(RX_FSK));
     RX_FSK_I_abs=abs(imag(RX_FSK));
     RX_FSK_L_i = find (RX_FSK_R_abs>=RX_FSK_I_abs);
     RX_FSK(RX_FSK_L_i)=0;
     RX_FSK_H_i = find (RX_FSK_R_abs<RX_FSK_I_abs);
     RX_FSK(RX_FSK_H_i)=1;
     RX_FSK_fn = pskdemod(RX_FSK_fn.^2,2);
     
      if SNR(n)==0
         RX_QAM_fn_0=RX_QAM_fn;
     elseif SNR(n)==16
         RX_QAM_fn_15=RX_QAM_fn;
     elseif SNR(n)==30
         RX_QAM_fn_30=RX_QAM_fn;
     end
     %QAM demodulation 
     symbol_data_block = qamdemod(RX_QAM_fn,M,'bin');  
     block_data_binary = de2bi(symbol_data_block,k);
     RX_QAM_fn = block_data_binary(:);% Return data in column vector 
    
     %error detection 
    [number_OOK,ratio_OOK]=biterr(OOK_mod,RX_OOK); %OOK_error detection
    [number_OOK_fn,ratio_OOK_fn]=biterr(OOK_mod_fn,RX_OOK_fn); %OOK_error detection
    
    [number_PRK,ratio_PRK]=biterr(OOK_mod,RX_PRK); %PRK_error detection
    [number_PRK_fn,ratio_PRK_fn]=biterr(OOK_mod,RX_PRK_fn); %PRK_error detection

    [number_FSK,ratio_FSK]=biterr(OOK_mod,RX_FSK); %FSK_error detection
    [number_FSK_fn,ratio_FSK_fn]=biterr(OOK_mod,RX_FSK_fn); %FSK_error detection
    
    [number_QAM_fn,ratio_QAM_fn]=biterr(OOK_mod,RX_QAM_fn); %QAM_error detection
    

     Error_OOK(n)=ratio_OOK;
     Error_OOK_fn(n)=ratio_OOK_fn;
     
     Error_PRK(n)=ratio_PRK;
     Error_PRK_fn(n)=ratio_PRK_fn;     
     
     Error_FSK(n)=ratio_FSK;
     Error_FSK_fn(n)=ratio_FSK_fn;
     
     Error_QAM_fn(n)=ratio_QAM_fn;


end

%% plotting 
figure 
semilogy(SNR,Error_OOK_fn,'r','linewidth',3,'DisplayName','OOK') %plotting BER VS SNR
hold on
semilogy(SNR,Error_FSK_fn,'b','linewidth',3,'DisplayName','FSK') %plotting BER vs SNR
hold on
semilogy(SNR,Error_PRK_fn,'g','linewidth',3,'DisplayName','PRK') %plotting BER vs SNR
hold on
semilogy(SNR,Error_QAM_fn,'k','linewidth',3,'DisplayName','QAM') %plotting BER vs SNR
ylabel('BER'), xlabel('SNR') ,grid on;
hold off
legend; 

figure 
semilogy(SNR,Error_OOK,'r','linewidth',3,'DisplayName','OOK') %plotting BER VS SNR
hold on
semilogy(SNR,Error_FSK,'b','linewidth',3,'DisplayName','FSK') %plotting BER vs SNR
hold on
semilogy(SNR,Error_PRK,'g','linewidth',3,'DisplayName','PRK') %plotting BER vs SNR
hold on
semilogy(SNR,Error_QAM_fn,'k','linewidth',3,'DisplayName','QAM') %plotting BER vs SNR
ylabel('BER'), xlabel('SNR') ,grid on;
hold off
legend; 


figure
subplot(4,2,1);
semilogy(SNR,Error_OOK);
xlabel('Range(dB)'), ylabel('Error'),grid on;
title('Error against SNR OOK','FontSize',12);
subplot(4,2,2);
semilogy(SNR,Error_OOK_fn);
xlabel('Range(dB)'), ylabel('Error'),grid on;
title('Error against SNR OOK-fn','FontSize',12);
subplot(4,2,3);
semilogy(SNR,Error_PRK);
xlabel('Range(dB)'), ylabel('Error'),grid on;
title('Error against SNR PRK','FontSize',12);
subplot(4,2,4);
semilogy(SNR,Error_PRK_fn);
xlabel('Range(dB)'), ylabel('Error'),grid on;
title('Error against SNR PRK-fn','FontSize',12);
subplot(4,2,5)
semilogy(SNR,Error_FSK);
xlabel('Range(dB)'), ylabel('Error'),grid on;
title('Error against SNR FSK','FontSize',12);
subplot(4,2,6)
semilogy(SNR,Error_FSK_fn);
subplot(4,2,7:8)
semilogy(SNR,Error_QAM_fn);
xlabel('Range(dB)'), ylabel('Error'),grid on;
title('Error against SNR FSK-fn','FontSize',12);

% constillation diagram 
sPlotFig = scatterplot(RX_QAM_fn_15,1,0,'g.');
hold on
scatterplot(QAM_mod_fn,1,0,'k*',sPlotFig)

sPlotFig2 = scatterplot(RX_QAM_fn_30,1,0,'g.');
hold on
scatterplot(QAM_mod_fn,1,0,'k*',sPlotFig2)

sPlotFig3 = scatterplot(RX_QAM_fn_0,1,0,'g.');
hold on
scatterplot(QAM_mod_fn,1,0,'k*',sPlotFig3)
hold off






