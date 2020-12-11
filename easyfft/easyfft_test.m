clear


% generate data
dt = 1/600;
N = 2500;
t = 0:dt:(N-1)*dt;                     
x = 5.9*sin(2*pi*11*t) + 35.4*sin(2*pi*40*t);      

y = fft(x);                  
f = (0:length(y)-1)*(1/dt)/length(y);        
                       

for i = 1:length(y)
    
    % things i've tried
%     y(i) = y(i)*(dt);
    y(i) = y(i)*2/N;

    % filter out nyquist
    if f(i)>(1/dt)/2
        y(i) = 0;
    end

end

% get magnitude
m = abs(y);  

% plot it
figure
hold on 
plot(f,m)
title('Magnitude')
hold off


power_t = sum(dot(x,x))*dt
power_f = sum(abs(y).^2)


%%
n = 0:319;
x = cos(pi/4*n)+randn(size(n));
[pxx,w] = periodogram(x);
figure
plot(w,(pxx))

y = fft(x);

p = (abs(y).^2)/length(x);
figure
plot(p)