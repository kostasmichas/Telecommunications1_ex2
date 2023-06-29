%A1

T=10^(-2);
over=10;
Ts=T/over;
A=4;
a=0.5;
figure('Name','Srrc pulse');
grid on;
hold on;
[phi,t] = srrc_pulse(T, over, A, a);
plot(t,phi);
legend('a=0.5')
hold off;

Nf=2048;
Fs=1/Ts;
DT=Fs/Nf;
F=-Fs/2:DT:(Fs/2)-1/Nf;
f_axis=[-0.5:1/Nf:0.5-1/Nf];
F_axis=Fs*f_axis;
figure('Name',' |Î¦(F)|');
grid on;
hold on;

X_f=fftshift(fft(phi,Nf))*Ts;
mf=abs(X_f);
plot(F_axis,mf);
xlim([-length(phi)-100 length(phi)+100]);
    
legend('a=0.5');
hold off;

figure('Name','Power Spectral Density with semilogy ');
  
  
psd=abs(X_f).^2;
semilogy(F_axis,psd);
hold on;
grid on;
legend('a=0.5');

%A2
for i=1:5
N=100;
b=(sign(randn(N,1))+1)/2;



X_n=bits_to_2PAM(b);
delays = [0:Ts:T];

phi_shifted=zeros(length(X_n), length(t) + (length(X_n) -1)*length(delays));

for n =1:length(X_n)
    if n==1
        phi_shifted(n,1:length(t))=X_n(1).*phi;
    else 
        phi_shifted(n, ((n-1)*length(delays):(n-1)*length(delays)+length(t)-1)+1)=X_n(n).*phi;
    end
end

X=sum(phi_shifted,1);
t_x=linspace(min(X), max(X), length(X));


figure('Name',sprintf('Signal X(t) with 2 PAM for i=%d',i));
plot(t_x,X);
xlabel('Time');
ylabel('Amplitude');

grid on;

Sx=(var(X_n)^2/T)*psd;


%A3
%a

P=fftshift((fft(X, Nf)))*Ts;
Ttotal=t(end);
P_x=(abs(P).^2)/Ttotal;
figure('Name',sprintf('Px(F) for i=%d',i));
plot(F_axis,P_x);
semilogy(F_axis,P_x);
grid on;
end

%b
K=500;
for i=1:K
   b=(sign(randn(N,1))+1)/2;
   X_n=bits_to_2PAM(b);

   P_x=(abs(P).^2)/Ttotal;
   P_x_array(:, i) = P_x;
   
end

P_x_mean = mean(P_x_array(:))
P_x_mean2 = mean(P_x_array, 2);
Sx=(var(X_n)^2/T)*psd

figure('Name',' Theoritical and Estimated Power Spectral Density with semilogy');
semilogy(F_axis, P_x_mean2, 'b-', 'LineWidth', 2);
hold on;
semilogy(F_axis, Sx, 'r-', 'LineWidth', 2);
grid on;
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density');
legend('Estimated PSD', 'Theoretical PSD');




%A4
for i=1:5


X_n=bits_to_4PAM(b);

phi_shifted=zeros(length(X_n), length(t) + (length(X_n) -1)*length(delays));
for n =1:(length(X_n)/2)
    if n==1
        phi_shifted(n,1:length(t))=X_n(1).*phi;
    else 
        phi_shifted(n, ((n-1)*length(delays):(n-1)*length(delays)+length(t)-1)+1)=X_n(n).*phi;
    end
end

X=sum(phi_shifted,1);
t_x=linspace(min(X), max(X), length(X));
figure('Name',sprintf('Signal X(t) with 4 PAM for i=%d',i));
plot(t_x,X);
xlabel('Time');
ylabel('Amplitude');

grid on;

Sx=(var(X_n)^2/T)*psd;

P=fftshift((fft(X, Nf)))*Ts;
Ttotal=t(end);
P_x=(abs(P).^2)/Ttotal;
figure('Name',sprintf('Px(F) with 4 PAM for i=%d',i));
plot(F_axis,P_x);
semilogy(F_axis,P_x);
grid on;
end

K=500;
for i=1:K
   b=(sign(randn(N,1))+1)/2;
   X_n=bits_to_4PAM(b);
   P_x=(abs(P).^2)/Ttotal;
   P_x_array(:, i) = P_x;
   
end

P_x_mean = mean(P_x_array(:))
P_x_mean2 = mean(P_x_array, 2);
Sx=(var(X_n)^2/T)*psd

figure('Name','Estimated and Theoritical Power Spectral Density with semilogy with 4 PAM');
semilogy(F_axis, P_x_mean2, 'b-', 'LineWidth', 2);
hold on;
semilogy(F_axis, Sx, 'r-', 'LineWidth', 2);
grid on;
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density');
legend('Estimated PSD', 'Theoretical PSD');


%A5

T_new=2*T;
over_new=2*over;

[phi,t] = srrc_pulse(T_new, over_new, A, a);

X_f=fftshift(fft(phi,Nf))*Ts;
mf=abs(X_f);
plot(F_axis,mf);
xlim([-length(phi)-100 length(phi)+100]);

psd=abs(X_f).^2;
for i=1:5
X_n=bits_to_2PAM(b);
delays = [0:Ts:T];

phi_shifted=zeros(length(X_n), length(t) + (length(X_n) -1)*length(delays));

for n =1:length(X_n)
    if n==1
        phi_shifted(n,1:length(t))=X_n(1).*phi;
    else 
        phi_shifted(n, ((n-1)*length(delays):(n-1)*length(delays)+length(t)-1)+1)=X_n(n).*phi;
    end
end

X=sum(phi_shifted,1);
t_x=linspace(min(X), max(X), length(X));

Sx=(var(X_n)^2/T)*psd;



P=fftshift((fft(X, Nf)))*Ts;
Ttotal=t(end);
P_x=(abs(P).^2)/Ttotal;
figure('Name',sprintf('Px(F) with 2T for i=%d',i));
plot(F_axis,P_x);
semilogy(F_axis,P_x);
grid on;
end


K=500;
for i=1:K
   b=(sign(randn(N,1))+1)/2;
   X_n=bits_to_2PAM(b);
   P_x=(abs(P).^2)/Ttotal;
   P_x_array(:, i) = P_x;
   
end

P_x_mean = mean(P_x_array(:))
P_x_mean2 = mean(P_x_array, 2);
Sx=(var(X_n)^2/T)*psd

figure('Name',sprintf('Estimated and Theoritical Power Spectral Density with semilogy with 2T for i=%d',i));
semilogy(F_axis, P_x_mean2, 'b-', 'LineWidth', 2);
hold on;
semilogy(F_axis, Sx, 'r-', 'LineWidth', 2);
grid on;
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density');
legend('Estimated PSD A5', 'Theoretical PSD A5');