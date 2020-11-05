clear


load test.mat

%%
newf = 0:1e-6:.5;

for i = 1:80
    y_avg(i) = 0 ;
end

newf = newf(1:10000);
y_avg = y_avg(1:10000);

figure
hold on 
plot(newf,y_avg)
set(gca, 'XScale', 'log')
hold off 

sqry_avg = zeros(size(y_avg));
for i = 1:length(y_avg)
    sqry_avg(i) = sqrt(y_avg(i));
end
figure
hold on 
plot(newf,sqry_avg)
set(gca, 'XScale', 'log')
hold off 

%% FL STUDIO 
y_avg = sqry_avg;
col1 = [85, 65, 113]/255;
col2 = [255, 140, 65]/255;
delta_color = col2-col1;

scaled_y_avg = y_avg / max(abs(y_avg));

figure
hold on 

for i = 1:1:length(y_avg)
    
    col = col1 + delta_color*scaled_y_avg(i);
    
    plot([newf(i),newf(i)],[0 1],'color',col)
end
xlim([1e-4,100])
set(gca, 'XScale', 'log')
set(gca,'Color',col1)
hold off 

%% 
y_avg = sqry_avg;
col1 = [1, 1, 1];
col2 = [0, 0, 0];
delta_color = col2-col1;

scaled_y_avg = y_avg / max(abs(y_avg));

figure
hold on 

for i = 1:1:length(y_avg)
    
    col = col1 + delta_color*scaled_y_avg(i);
    
    plot([newf(i),newf(i)],[0.001 .999],'color',col)
    
    if i == length(y_avg)
        back_color = col;
    end
end
xlim([1e-4,100])
xlabel('Frequency (hz)')
ylabel('Normalized Magnitude')
set(gca, 'XScale', 'log')
set(gca,'Color',back_color)
hold off 

    
    
    
    
% figure
% hold on 
% plot(newf)
% hold off 
