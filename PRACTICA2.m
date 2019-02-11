%% PSAVC Pràctica 2 - Cancelació d'interfèrencia
%Constants
L = 10^3*2;   %Nombre experiments
mu = 0.01; %Constant de canvi

%SENYAL D'INTERES
%Cas soroll blanc
Pw=1;
%s_interes = randn(L,1)*sqrt(Pw);
%Cas constant (OPCIONAL)
s_interes = zeros(L,1)+10;

%SENYAL INTERFERENT
%Cas soroll autocorrelat (OPCIONAL)
Pw = 1;
v_autocorr = [1,0.5,0.3];
s_interf = randn(L,length(v_autocorr))*v_autocorr'*sqrt(Pw);

%Cas senyal sinusoidal (P3)
P = 0;
n = (1:L);
fo = 1/4;
valor_continua = 0;
s_interf = (sqrt(2*P)*cos(2*pi*fo*n))'+valor_continua;

%Cas primer zeros i després senyal sinusodial (P4)
P = 10;
%s_interf = ([zeros(L/2,1);repmat(sqrt(2*P),L/2,1)].*cos(2*pi*fo*n'));

%Senyal x
x = s_interes + s_interf;

%Inicialització de vectors
h0 = zeros(L,1)+4;  
h1 = zeros(L,1)+4;
y  = zeros(L,1);
error = zeros(L,1);

%ALGORITME ADAPTATIU
for n = (3:L) 
y(n) = x(n)-h0(n)*x(n-1)-h1(n)*x(n-2);

error(n) = (y(n)-x(n))^2;

h0(n+1)=h0(n)+mu*x(n-1)*y(n);
h1(n+1)=h1(n)+mu*x(n-2)*y(n);
end

%VISUALITZACIONS
close all
%Plot de les dos constants
figure(4)
plot(h0,h1)
title('Plot 2D de les constants h0 i h1')
axis([-2 4 -2 4])
ylabel('h1')
xlabel('h0')

%Plot fft
n_samples = L/2;
h = [[1;-h0(n);-h1(n)];zeros(n_samples,1)];
len = length(h);
figure(3)
plot((1:len)/len,abs(fft(h)))
title('Transformada de fourier del filtre amb Zero Padding')

%Plot error quadràtic
figure(1)
plot(error)
title('Error quadràtic')
ylabel('Error quadràtic')
xlabel('N')

%Plot H0 i H1 
figure(2)
plot(h0,'DisplayName','h0')
hold on
plot(h1,'DisplayName','h1')
title('Valors h0 i h1')
legend
ylabel('Valor')
xlabel('N')

%Plot senyal original vs filtrada
figure(5)

hold on
plot(y,'DisplayName','senyal filtrada')
hold on
plot(x,'DisplayName','senyal original')
hold on
plot(s_interes,'DisplayName','senyal interes')
legend
title('Senyal original i senyal filtrada')
