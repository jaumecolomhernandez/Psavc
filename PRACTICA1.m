%% PRÀCTICA 1 PSAVC DETECCIÓ - JAUME COLOM I JR HERNANDEZ
%% Exercici 1 - Generar senyal i observacions
%definició de constants i variables
a=0.0707;
b=0.0707;
fo=1/400;
var=1;
N=2000;

%Generació de seyal
v=0:N-1;    %vector enters ascendent
H=[cos(2*pi*fo*v.'),sin(2*pi*fo*v.')]; %matriu 2x1 amb senyals sinusoidals
O=[a,b];    %Vector constants
s=H*O.';    %senyal calculada

%Hipotesis
h1=s+randn(N,1);
h0=randn(N,1);

plot(h1)

%% Exercici 2 - Càlcul Pfa
L=10^5;             %nombre d'experiments
umbral=7.35;      %valor de l'umbral
N=2000;

ho=randn(N,L);  %hipotesis sense senyal
t=s.'*ho;
num = t>umbral;

Pfa1=sum(num)/L           %falses alarmes entre nombre d'experiments(L)

%% Exercici 3 - Càlcul Pfa amb diferent umbral
L=10^5;             %nombre d'experiments
umbral=9.22;      %valor de l'umbral
N=2000;

ho=randn(N,L);  %hipotesis sense senyal
t=s.'*ho;
num = t>umbral;

Pfa2=sum(num)/L           %falses alarmes entre nombre d'experiments(L)

%% Exercici 4 - Escombrat umbrals i gràfica
L=10^4;                 %nombre d'experiments
n_umbrals=10;
N=2000;

ho=randn(N,L);  %hipotesis sense senyal
t=s.'*ho;

count=t'>linspace(0,12,n_umbrals) ;
v=sum(count)/L;       

semilogy(linspace(0,12,10),v);  %plot logaritmic de pfa

%% Exercici 5 - Càlcul Pd
L=10^4;             %nombre d'experiments
umbral=7.35;      %valor de l'umbral
N=2000;

h1=repmat(s,1,L)+randn(N,L);  %hipotesis sense senyal
t=s.'*h1;
num = t>umbral;

Pd=sum(num)/L

%% Exercici 6 - Càlcul ROC
L=10^4;                 %nombre experiments

%Càlcul senyals
soroll=randn(N,L);
h1=repmat(s,1,L)+soroll;%hipotesis senyal 
h0=soroll;  %hipotesis sense senyal
t=s.'*h1;
t2=s.'*h0;

%Càlcul probabilitats
umbrals = linspace(0,12,100);
c_pd=t'>umbrals;
c_pfa=t2'>umbrals;
p_d=sum(c_pd)/L;     
p_fa=sum(c_pfa)/L;

plot(p_fa,p_d);         %plot de la ROC
axis([0 0.5 0.5 1]);

%% Exercici 7 - ROC amb amplituds aleatòries
%Declaració de variables
L=10^4;
N=2000;

%Càlcul senyal i estimadors
O=randn(2,1)*sqrt(0.01);
s=H*O;
soroll=randn(N,L);
h1=repmat(s,1,L)+soroll;%hipotesis senyal 
h0=soroll;              %hipotesis sense senyal
t=diag(h1.'*H*H.'*h1);
t2=diag(h0.'*H*H.'*h0);

%Càlcul probabilitats
umbrals = linspace(1,200000,1000);
c_pd=t>umbrals;
c_pfa=t2>umbrals;
p_d=sum(c_pd)/L;     
p_fa=sum(c_pfa)/L;

plot(p_fa,p_d);
%axis([0 0.5 0.5 1]);
%%  Exercici 8 - Càlcul teòric
%Considerant el vector O conegut
pot=1;
umbral=4;
%tret de la solució del examen
var = sqrt(pot*(a^2+b^2)/(2*N));
pfa1=qfunc(umbral/var);
pd=qfunc(qfuncinv(pfa)-sqrt(N*(a^2+b^2)/(2*pot)));

%Considerant vector O desconegut, vector aleatori
potr=0.01;
var=sqrt(potr*N/2+pot);
pfa2=qfunc(umbral/var);

%NO ENS COINCIDEIXEN AMB EL VALOR EXPERIMENTAL!!

%% PRÀCTICA 1 PSAVC DETECCIÓ - JAUME COLOM I JR HERNANDEZ
%% Exercici 1 - Generar senyal
%definició de constants
A=1.5;  %Amplitud del pols
B=1;    %Amplitud constant
Pw=1;   %Potència del soroll
N=9000; %Duració total
M=3000; %Duració del pols

%Creació del vector S
s = [cat(1,ones(M,1),zeros(N-M,1)),ones(N,1)];
%Vector de constants a estimar
O=[A;B];
%Generació de senyal
x=s*O+randn(N,1)*sqrt(Pw);
plot(x);
%% Exercici 2 - Estimador Òptim
%Estimador òptim multivariable
est= (s.'*s)^-1*s.'*x;

%% Exercici 3 - Estimació mitjana i variança
%constants i variables
L=10^4;                 %nombre d'experiments
v_est_a=[];             %vector d'estimacions de a
v_est_b=[];             %vector d'estimacions de b
%loop d'operacions
x=repmat(s*O,1,L)+randn(N,L)*sqrt(Pw);
est = (s.'*s)^-1*s.'*x;

%Estimació de mitja i variància de A
est_mean_a=1/L*est(1,:)*ones(L,1);
est_var_a=1/L*(est(1,:)-est_mean_a).^2*ones(L,1);

%Estimacio de mitja i variança de B
est_mean_b=1/L*est(2,:)*ones(L,1);
est_var_b=1/L*(est(2,:)-est_mean_b).^2*ones(L,1);

%% Exercici 4 - Gràfica variància
N = 9000;       %Duracio total
est_mean_a=[];
est_var_a=[];
soroll=randn(N,L);

for M = [100, 500, 1000, 2000, 3000, 4000, 7000, 8000,8500, 8900]; 
    v_est_a=[]; %Vector estimacions de a
    %Generem s per cada valor de M
    s = [cat(1,ones(M,1),zeros(N-M,1)),ones(N,1)];
    vector=s*O;
    operador=(s.'*s)^-1*s.';
    %loop d'operacions
    x=repmat(vector,1,L)+soroll*sqrt(Pw);
    est = (s.'*s)^-1*s.'*x;

    %Estimació de mitja i variància de A
    est_mean_a(end+1)=1/L*est(1,:)*ones(L,1);
    est_var_a(end+1)=1/L*(est(1,:)-est_mean_a(end)).^2*ones(L,1);
end
figure(1);
plot(est_mean_a)
figure(2);
plot(est_var_a);

%% Exercici 5 - Valor de salt
N = 9000;           %Duració total
soroll=randn(N,L);

%vector resultats
est_var_c=[];
for M = [100, 500, 1000, 2000, 3000, 4000, 7000, 8000, 8500, 8900];
    %Càlcul senyal i estimador
    s = [cat(1,ones(M,1),zeros(N-M,1)),ones(N,1)];
    vector=s*O;
    operador=(s.'*s)^-1*s.';
    x=repmat(vector,1,L)+soroll*sqrt(Pw);
    est = (s.'*s)^-1*s.'*x;
    v_est_c=est(1,:)-est(2,:);
    
    %Càlcul estadístiques mostrals
    est_mean_c(end+1)=1/L*v_est_c*ones(L,1);   
    est_var_c(end+1)=1/L*(v_est_c-est_mean_c(end)).^2*ones(L,1); 
end
%plots
plot(est_var_c);

%% Exercici 6 - Càlcul teòric
var_est_a=[];
N=9000;
Pw=1;
M = [100, 500, 1000, 2000, 3000, 4000, 7000, 8000,8500, 8900];
var_est_a=N./(M.*(N-M)).*Pw;
plot(var_est_a);
hold on;
var_est_c=[];
M = [100, 500, 1000, 2000, 3000, 4000, 7000, 8000,8500, 8900];
var_est_c=Pw.*(N+3*M)./(M.*(N-M));
plot(var_est_c);

%OBSERVEM QUE ELS DOS COINCIDEIXEN AMB ELS VALORS EXPERIMENTALS

%% Exercici 7 - Estimació de soroll
%Càlcul projector ortogonal
%Proj=eye(3000)-s*(s.'*s)^-1*s.';
%Estimació de potència de la senyal
%est_pow=1/3000*(Proj*x)'*(Proj*x);
L=10000;

%Vectores resultados
est_mean_p=[];
est_var_p=[];

%WARNING SI FEM SERVIR M>2000 SALTA ERRROR OUT OF MEMORY
%SI M<2000 FUNCIONA BASTANT RÀPID
for M_usr = linspace(100,1000,8)
    %definició de M i N
    M=round(M_usr);
    N=3*M;
    
    %Generació senyal x
    s = [cat(1,ones(M,1),zeros(N-M,1)),ones(N,1)];
    O=[1.5;1];  %Haurien de ser els valors estimats
    senyal=s*O;
    %codi vectoritzat, aquí és calcula totes les mostres de cop 
    mostra_soroll = randn(N,L);
    x=repmat(senyal,1,L)+mostra_soroll; 
    
    %Càlcul projector ortogonal
    Proj=eye(N)-s*(s.'*s)^-1*s.';
    
    %Càlcul estimador
    m_est_p= 1/N*(Proj*x)'*(Proj*x);
    v_est_p=diag(m_est_p);  %els valors diagonals són les estimacions

    %Càlcul mitjana i variància estimador de p
    est_mean_p(end+1)=1/L*v_est_p'*ones(L,1);
    est_var_p(end+1)=1/L*(v_est_p-repmat(est_mean_p(end),L,1)).^2'*ones(L,1);
end
%plots
figure(1)
plot(est_mean_p);
figure(2)
plot(est_var_p);

%% Exercici 8 - Estimació SNR
%Estimació SNR
%est_snr=est(1)^2/est_pow;

%Gràfica de SNR en funció de M en el cas N=3M
est_mean_snr=[];
est_var_snr=[];
var_crb=[];
for M = linspace(500,1500,8)
    %Constants
    M=round(M);
    N=3*M;   %nou valor de N per cada valor de M
    
    %Generació senyal 
    s = [cat(1,ones(M,1),zeros(N-M,1)),ones(N,1)];
    senyal=s*O;
    mostra_soroll = randn(N,L);
    x=repmat(senyal,1,L)+mostra_soroll;
    
    %Càlcul projector
    Proj=eye(N)-s*(s.'*s)^-1*s.';
    
    %Estimació de potència del soroll
    m_est_p= 1/N*(Proj*x)'*(Proj*x);
    v_est_p=diag(m_est_p);  %els valors diagonals són les estimacions
    
    %Càlcul estimació constants A i B
    calcul_previ=(s.'*s)^-1*s.';
    est = calcul_previ*x;

    %Càlcul estimació SNR
    v_est_snr=est(1)^2./v_est_p;

    %Càlcul de mitjana i variància mostrals
    est_mean_snr(end+1)=1/L*v_est_snr'*ones(L,1);
    est_var_snr(end+1)=1/L*(v_est_snr-est_mean_snr(end)).^2'*ones(L,1);
    var_crb(end+1) = 9/N;
end

%Plots
figure(1)
plot(est_mean_snr);
figure(2)
plot(est_var_snr);
hold on;
plot(var_crb);


