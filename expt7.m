fprintf('----------------------------------------------');
fprintf('Sampling and Reconstruction');
fprintf('----------------------------------------------\n');
F=input('Input1: Frequenncy of analog sinusoid in Hz:');
fn=2*F;
fprintf('\n Output1: for sinusoid of freq. %d Hz, Nyquist rate is %d Hz \n \n',F, fn);
T=1/F;
t=0:0.001*T:4*T;
x=2*cos(2*pi*F*t);

%% UnderSampling
fprintf('\n\n----------------------------------------------');
fprintf('Case 1: UnderSampling');
fprintf('----------------------------------------------\n');
fsun=input('Input2: Choose appropriate Sampling frequency in Hz: ');
if fsun >= fn
    error('For Undersampling, sampling frequency must be samller than Nyquist rate')
endif
t1=0:(1/fsun):4*T;
xun=2*cos(2*pi*F*t1);

%% NyquistSampling
fprintf('\n\n----------------------------------------------');
fprintf('Case 2: NyquistSampling');
fprintf('----------------------------------------------\n');
fsn=input('Choose appropriate Sampling frequency in Hz:');
if fsn ~= fn
    error('For Nyquistsampling, sampling frequency must be equal to Nyquist rate')
endif
t2=0:(1/fsn):4*T;
xn=2*cos(2*pi*F*t2);

%% OverSampling
fprintf('\n\n----------------------------------------------');
fprintf('Case 3: OverSampling');
fprintf('----------------------------------------------\n');
fsov=input('Choose appropriate Sampling frequency in Hz:');
if fsov <= fn
    error('For Oversampling, sampling frequency must be greater than Nyquist rate')
endif
t3=0:(1/fsov):4*T;
xov=2*cos(2*pi*F*t3);

%% ReconstructedAnlog Signal
% Case 1: From UnderSampled signal
xr1=zeros(1,length(t));
for k=1:length(t)
    for p= 1:length(t1)
        xr1(k)=xr1(k)+xun(p)*sinc((k-1)*fsun*0.001*T-(p-1));
    end
end

% Case 2: From NyquistSampled signal
xr2=zeros(1,length(t));
for k=1:length(t)
    for p= 1:length(t2)
        xr2(k)=xr2(k)+xn(p)*sinc((k-1)*fsn*0.001*T-(p-1));
    end
end

% Case 3: From OverSampled signal
xr3=zeros(1,length(t));
for k=1:length(t)
    for p= 1:length(t3)
        xr3(k)=xr3(k)+xov(p)*sinc((k-1)*fsov*0.001*T-(p-1));
    end
end

%% Plotting
figure
subplot(4,2,1);
plot(t,x,'-','LineWidth',2);
title(['Sinusoid analog signal (' num2str(F) ' Hz) to be sampled']);
xlabel('time');
ylabel('amplitude');

subplot(4,2,2);
plot(t,x,'-','LineWidth',2);
title(['Sinusoid analog signal (' num2str(F) ' Hz) to be reconstructed']);
xlabel('time');
ylabel('amplitude');

subplot(4,2,3);
plot(t,x,'--','LineWidth',1);
hold on;
stem(t1,xun,'k','LineWidth',2,'MarkerFaceColor','red', 'MarkerEdgeColor','yellow');
title(['UnderSampled signal with ' num2str(fsun) ' Hz Sampling rate']);
xlabel('time');
ylabel('amplitude');
%legend('Original analog', 'UnderSampled');

subplot(4,2,5);
plot(t,x,'--','LineWidth',1);
hold on;
stem(t2,xn,'k','LineWidth',2,'MarkerFaceColor','red', 'MarkerEdgeColor','yellow');
title(['NyquistSampled signal with ' num2str(fn) ' Hz Sampling rate']);
xlabel('time');
ylabel('amplitude');
%legend('Original analog', 'NyquistSampled');

subplot(4,2,7);
plot(t,x,'--','LineWidth',1);
hold on;
stem(t3,xov,'k','LineWidth',2,'MarkerFaceColor','red', 'MarkerEdgeColor','yellow');
title(['OverSampled signal with '  num2str(fsov) ' Hz Sampling rate']);
xlabel('time');
ylabel('amplitude');
%legend('Original analog', 'OverSampled');

subplot(4,2,4);
plot(t,x,'--','LineWidth',1);
hold on;
plot(t,xr1,'-','LineWidth',2);
stem(t1,xun,'k','LineWidth',2,'MarkerFaceColor','red', 'MarkerEdgeColor','yellow');
title('Reconstructed signal (red) from UnderSampled signal');
ylim([-2.5,2.5]);
xlabel('time');
ylabel('amplitude');
%legend('Original analog','Reconstructed signal','UnderSampled signal',"location", "southoutside");

subplot(4,2,6);
plot(t,x,'--','LineWidth',1);
hold on;
plot(t,xr2,'-','LineWidth',2);
stem(t2,xn,'k','LineWidth',2,'MarkerFaceColor','red', 'MarkerEdgeColor','yellow');
title('Reconstructed signal (red) from NyquistSampled signal');
ylim([-2.5,2.5]);
xlabel('time');
ylabel('amplitude');
%legend('Original analog','Reconstructed signal','NyquistSampled signal');

subplot(4,2,8);
plot(t,x,'--','LineWidth',1);
hold on;
plot(t,xr3,'-','LineWidth',2);
stem(t3,xov,'k','LineWidth',2,'MarkerFaceColor','red', 'MarkerEdgeColor','yellow');
title('Reconstructed signal(red) from OverSampled signal');
ylim([-2.5,2.5]);
xlabel('time');
ylabel('amplitude');
%legend('Original analog','Reconstructed signal','OverSampled signal');