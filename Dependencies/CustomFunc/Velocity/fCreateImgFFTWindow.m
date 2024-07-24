function [FFT_Window_Output] = fCreateImgFFTWindow(size_Img)

M = size_Img(1);
N = size_Img(2);
w1 = cos(linspace(-pi/2, pi/2, M));
w2 = cos(linspace(-pi/2, pi/2, N));
FFT_Window_Output = w1' * w2; 
