clear all
clc
%close all
NumSNR= 10 ^(1);
SNR=0:2:30; %range of SNR
RandomSignal=randi([0,1],1,NumSNR); % generate stream of random bits

%% Encoding :
    %first repetiton
    %second linear block code
    %third convolution block code

%% first repetition
% using repmat to repeat bit then reshaping it using reshape to make it in
% one row
%{
encoded_data_3=reshape(repmat(RandomSignal,3,1),[],1);% take data repeat every bit 3 times
encoded_data_5=reshape(repmat(RandomSignal,5,1),[],1);% take data repeat every bit 5 times
encoded_data_11=reshape(repmat(RandomSignal,11,1),[],1);% take data repeat every bit 11 times
%}
%% second Linear block code
%{
k = 4;        % Data length
m = 7;        % Code length
% first create polynomial using cyclpoly then generate parity matrix
% then generate G matrix
% then get the decoding table syndrome

Polynomial = cyclpoly(m,k);
P_M = cyclgen(m,Polynomial);% parity matrix
G_M = gen2par(P_M);%generator matrix
d_t = syndtable(P_M);%decoding table using syndrome to detect error

% Encode the message sequence by using the generator matrix.
data_mat = vec2mat (RandomSignal,k); %convert data from vector form to matrix form of column length k
[G,U] = size(data_mat);
encoded_data_LBC = [];

for i  = 1 : G
    output = encode(data_mat(i,:),m,k,'linear/binary',G_M);

    encoded_data_LBC = [encoded_data_LBC , output] ;
end;
%}
encoded_data_LBC = encode(RandomSignal,7,4,'hamming/binary');
%% third Convolutional Codes
constLength = 9;
traceBack = 5 * constLength;
polynomial = [657 435];
trellis = poly2trellis(constLength , polynomial);

encoded_data_Conv = convenc(RandomSignal , trellis);

%% Modulation :
%{

%% first repetition
%{
Rep3_Mod_Data = pskmod(encoded_data_3,2);
PTX_R_3=mean(Rep3_Mod_Data.^2);
Rep3_Mod_Data = Rep3_Mod_Data*sqrt(1/3);

Rep5_Mod_Data = pskmod(encoded_data_5,2);
PTX_R_5=mean(Rep5_Mod_Data.^2);
Rep5_Mod_Data = Rep5_Mod_Data*sqrt(1/5);

Rep11_Mod_Data = pskmod(encoded_data_11,2);
PTX_R_11=mean(Rep11_Mod_Data.^2);
Rep11_Mod_Data = Rep11_Mod_Data*sqrt(1/11);
%}
%% second Linear block code

LBC_Mod_Data= pskmod(encoded_data_LBC,2);
PTX_LBC= mean(LBC_Mod_Data.^(2));
LBC_Mod_Data = LBC_Mod_Data*sqrt(1/1.75);

%% third Convolutional Codes

Conv_Mod_Data = pskmod(encoded_data_Conv,2);
PTX_conv= mean(Conv_Mod_Data.^(2));
Conv_Mod_Data = Conv_Mod_Data*sqrt(1/2);

%% fourth uncoded

uncoded_Mod_Data=pskmod(RandomSignal,2);
PTX_uncoded = mean(uncoded_Mod_Data.^2);
%}
%% adding noise + demodulation + decoding + error detection


%% adding noise


for n=1:length(SNR)

     snr_i=10^(SNR(n)/10);
     % noise repetition
     %{
     noise_rep3=sqrt(1/(2*snr_i))*( randn(size(Rep3_Mod_Data))+1j*randn(size(Rep3_Mod_Data)));
     noise_rep5=sqrt(1/(2*snr_i))*( randn(size(Rep5_Mod_Data))+1j*randn(size(Rep5_Mod_Data)));
     noise_rep11=sqrt(1/(2*snr_i))*( randn(size(Rep11_Mod_Data))+1j*randn(size(Rep11_Mod_Data)));
     %}
     %noise lbc
     noise_LBC=sqrt(1/(2*snr_i))*( randn(size(LBC_Mod_Data))+1j*randn(size(LBC_Mod_Data)));
     %noise conv
     noise_conv=sqrt(1/(2*snr_i))*( randn(size(Conv_Mod_Data))+1j*randn(size(Conv_Mod_Data)));
     %uncoded
     %noise_BPSK=sqrt(1/(2 * snr_i) *( randn(size(RandomSignal)) + 1j * randn(size(RandomSignal)) ) ) ; %PRK noise


     %recieved repetition
     %{
     Rx_Rep3 = Rep3_Mod_Data + noise_rep3;
     Rx_Rep5 = Rep5_Mod_Data + noise_rep5;
     Rx_Rep11 = Rep11_Mod_Data + noise_rep11;
     %}
     %recieved repetition
     %Rx_LBC = awgn(LBC_Mod_Data ,SNR (n), 'measured');
     Rx_LBC = LBC_Mod_Data +  noise_LBC;
     %recieved repetition
     Rx_conv = Conv_Mod_Data + noise_conv;
     %recieved uncoded
     Rx_uncoded= awgn(uncoded_Mod_Data ,SNR (n), 'measured');



     %% demodulation


         %% first repetition
    %{
     Rep3_Demod_Data = pskdemod(Rx_Rep3,2);
     Rep5_Demod_Data = pskdemod(Rx_Rep5,2);
     Rep11_Demod_Data = pskdemod(Rx_Rep11,2);
     %}    
     %% second Linear block code
     LBC_Demod_Data = pskdemod(Rx_LBC,2);

         %% third Convolutional Codes

     Conv_Demod_Data = pskdemod(Rx_conv,2);

         %% fourth Convolutional Codes

     uncoded_Demod = pskdemod(Rx_uncoded,2);



     %% Decoding

          %% first repetition
%{
     for i=1:NumSNR
         Rep3_Dec_Data(i)= sum(Rep3_Demod_Data(3*i-2:3*i))>=2;
         Rep5_Dec_Data(i)= sum(Rep5_Demod_Data(5*i-4:5*i))>=3;
         Rep11_Dec_Data(i)= sum(Rep11_Demod_Data(11*i-10:11*i))>=6;


     end
 %}         
          %% second Linear block code
%{
      mat2 = vec2mat (LBC_Demod_Data,m);
     [y,r] = size(mat2);
     LBC_dec_Data =[];

     for j  = 1 : y
        p = decode(mat2(j,:),m,k,'linear/binary',G_M,d_t);
        LBC_dec_Data = [LBC_dec_Data , p] ;
     end
%}
 LBC_dec_Data= decode(LBC_Demod_Data,7,4,'hamming/binary');
          %% third Convolutional Codes

     Conv_Dec_Data = vitdec(Conv_Demod_Data,trellis,traceBack,'trunc','hard');

%% error and ber

          %% first repetition
%{
     [numoferrors_rep3(n),ratio_rep_3(n)]=biterr(RandomSignal,Rep3_Dec_Data);
     [numoferrors_rep5(n),ratio_rep_5(n)]=biterr(RandomSignal,Rep5_Dec_Data);
     [numoferrors_rep11(n),ratio_rep_11(n)]=biterr(RandomSignal,Rep11_Dec_Data);
%}
          %% second Linear block code

     [numoferrors_LBC(n),ratio_LBC(n)]=biterr(RandomSignal,LBC_dec_Data);



          %% third Convolutional Codes

     [numoferrors_conv(n),ratio_conv(n)]= biterr(RandomSignal,Conv_Dec_Data);

          %% fourth Convolutional Codes
          
     [numoferrors_uncoded(n),ratio_uncoded(n)]= biterr(RandomSignal,uncoded_Demod);




end

figure
%{
semilogy(SNR,ratio_rep_3,'r','linewidth',1.5,'DisplayName','rep3') %plotting BER VS SNR
hold on
semilogy(SNR,ratio_rep_5,'b','linewidth',1.5,'DisplayName','rep5') %plotting BER vs SNR
hold on
semilogy(SNR,ratio_rep_11,'g','linewidth',1.5,'DisplayName','rep11') %plotting BER vs SNR
hold on
%}
semilogy(SNR,ratio_conv,'k','linewidth',1.5,'DisplayName','conv') %plotting BER vs SNR
hold on
semilogy(SNR,ratio_uncoded,'m','linewidth',1.5,'DisplayName','uncoded') %plotting BER vs SNR
hold on
semilogy(SNR,ratio_LBC,'y','linewidth',1.5,'DisplayName','linear') %plotting BER vs SNR

ylabel('BER'), xlabel('SNR') ,grid on;
hold off
legend;
