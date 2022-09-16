function Sf = myFFTSmooth(S,n)
% very simple low pass filter
    S = S(:);
    fft_data_v = fft(S);
    s_fft_data_v = zeros(1,length(S));
    s_fft_data_v(1:n) = fft_data_v(1:n);
    s_fft_data_v(end-n:end) = fft_data_v(end-n:end);
    Sf = real(ifft(s_fft_data_v));
    Sf = Sf(:);
end