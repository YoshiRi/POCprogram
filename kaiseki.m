Data = Pp1 ;
for i = 1:256
    for j = 1:256
        if  Data(j,i) <0.05
            Data(j,i) = 0;
        end
    end
end

fData = fft2(Data);
afData = imag(fData);

figure;
plot(afData(128,:))

%%
fPp1 = fft2(Pp1);
afPp1 = imag(fPp1);

figure;
plot(afPp1(128,:))

%% gfilter

gfilt = [1 2 1; 2 4 2 ; 1 2 1];
figure;
mesh(conv2(Pp1,gfilt));
figure;
mesh(conv2(Pp2,gfilt));
