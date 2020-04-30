freqs =[500:50:3000]; %[[500:50:1000]];%1500], [1600:100:3000]];
for i = 1:length(freqs)
    freq = freqs(i) ;
    t = 0:1/48000:3;
    x = cos(2*pi*t*freq);
    x = x./(max(abs(x))+0.01);
    wavwrite(x, 48000, 16, ['D:\ImpedanceTube\Processing\2018.4.6\sins\sin', num2str(freq), '.wav']);
end    