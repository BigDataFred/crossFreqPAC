function [SNR] = compute_SNR(s1,s2)
%%
SNR= 10*log10(mean(s1.^2,2)./(mean(s2.^2,2)));