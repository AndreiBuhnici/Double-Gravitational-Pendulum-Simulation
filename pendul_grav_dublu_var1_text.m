% Pendulul gravitational dublu

clear; close all; clc;

RT=1; % selector real time(1) / slow motion(0)

g=9.80665; % m/s^2; acceleratia gravitationala terestra

% Parametrii fizici ai sistemului mecanic:

L1=2.2; L2=1.3; % m; lungimile tijelor

m1=1.2; m2=1.7; % kg; masele corpurilor
count=0;

% Conditiile initiale - orice unghiuri intre -180 si +180!

theta10=100; theta20=-10; % grade; unghiurile initiale

theta10=theta10*pi/180; theta20=theta20*pi/180; % conversia in rad

OM10=-30; OM20=20; % grade/s; vitezele unghiulare initiale

OM10=OM10*pi/180; OM20=OM20*pi/180; % conversia in rad/s

% Definirea duratelor caracteristice:

omega1=sqrt(g/L1); omega2=sqrt(g/L2); % pulsatii proprii ale componentelor

T1=2*pi/omega1; T2=2*pi/omega2; % perioade proprii ale componentelor

T=max(T1,T2); % "timp caracteristic" al miscarii pendulului dublu

ti=0; tf=5*T; N=200000; t=linspace(ti,tf,N); dt=t(2)-t(1); % timpul discret
tic;
while count<3
    
% Prealocare si valori de start:

theta1=zeros(1,N); theta2=theta1; % prealocare unghiuri

OM1=zeros(1,N); OM2=OM1; % prealocare viteze unghiulare

theta1(1)=theta10; theta2(1)=theta20; % unghiuri de start pas 1

theta1(2)=theta10+OM10*dt; theta2(2)=theta20+OM20*dt; % unghiuri de start pas 2

OM1(1)=OM10; OM2(1)=OM20; % valori de start ale vitezelor unghiulare


% Notatii ajutatoare:


miu=1+m1/m2; % coeficient adimensional

r=L2/L1; % coeficient adimensional

a11=miu; a22=r; % coeficienti diagonala principala (constanti)



for i=2:N-1 % ciclul recurentelor

    aux=theta2(i)-theta1(i);

    a21=cos(aux); a12=a21*r; % coeficienti diagonala secundara (variabili)

    % Pentru vitezele unghiulare curente folosim derivatele la stanga:

    OM1(i)=(theta1(i)-theta1(i-1))/dt; % viteza corpului 1 la pasul i

    OM2(i)=(theta2(i)-theta2(i-1))/dt; % viteza corpului 2 la pasul i

    b1=r*OM2(i)^2*sin(aux)-g/L1*miu*sin(theta1(i)); % termen "liber" 1

    b2=-OM1(i)^2*sin(aux)-g/L1*sin(theta2(i)); % termen "liber" 2

    A=[a11,a12;a21,a22]; B=[b1;b2]; % matrice sistem si coloana termeni liberi

    E=A\B; % rezolvarea sistemului liniar in forma matriceala

    eps1=E(1); eps2=E(2); % acceleratiile unghiulare curente

    % Recurentele de ordinul II:

    theta1(i+1)=2*theta1(i)-theta1(i-1)+dt^2*eps1; % corp 1

    theta2(i+1)=2*theta2(i)-theta2(i-1)+dt^2*eps2; % corp 2

end;

OM1(N)=(theta1(N)-theta1(N-1))/dt; % viteza unghiulara corp 1 la pasul N

OM2(N)=(theta2(N)-theta2(N-1))/dt; % viteza unghiulara corp 2 la pasul N



% Coordonate carteziene ale corpurilor:

x1=L1*sin(theta1); x2=x1+L2*sin(theta2); % coordonate orizontale

y1=-L1*cos(theta1); y2=y1-L2*cos(theta2); % coordonate verticale

if count==0
a1=x1;
a2=x2;
c1=y1;
c2=y2;
elseif count==1
a3=x1;
a4=x2;
c3=y1;
c4=y2;
elseif count==2
a5=x1;
a6=x2;
c5=y1;
c6=y2;
end
count=count+1;
m2=m2+0.3;
end

toc; % afiseaza timpul de calcul al solutiei numerice
% Energiile cinetica, potentiala, si respectiv totala:

%T=1/2*(m1*L1^2*OM1.^2+m2*(L1^2*OM1.^2+L2^2*OM2.^2+2*L1*L2*OM1.*OM2.*cos(theta2-theta1)));

%U=-g*((m1+m2)*L1*cos(theta1)+m2*L2*cos(theta2));

%H=T+U; % hamiltonianul sistemului - energia totala

figure(1);

Lmax=L1+L2; % semilatura cadrului grafic

coef=30; % controleaza dimensiunile grafice ale corpurilor

rg1=coef*m1^(1/3); rg2=coef*m2^(1/3); % raze "grafice"

tic; simt=0; % porneste cronometrul si initializeaza timpul simularii

while simt<=tf % ciclul grafic

  j=abs(t-simt)==min(abs(t-simt)); % cauta cel mai apropiat t din discretizare

  plot([0 a1(j) a2(j)],[0 c1(j) c2(j)],'-g','LineWidth',3); hold on; % tija 1 oscilator 1
  plot([0 a3(j) a4(j)],[0 c3(j) c4(j)],'-g','LineWidth',3); hold on; % tija 1 oscilator 2
  plot([0 a5(j) a6(j)],[0 c5(j) c6(j)],'-g','LineWidth',3); hold on; % tija 1 oscilator 3
  xlabel('x/m'); ylabel('y/m');

  plot(0,0,'.k','MarkerSize',10); % articulatie de suspensie

  plot(a1(j),c1(j),'.r','MarkerSize',rg1); % corpul 1 oscilator 1
  plot(a3(j),c3(j),'.r','MarkerSize',rg1); % corpul 1 oscilator 2
  plot(a5(j),c5(j),'.r','MarkerSize',rg1); % corpul 1 oscilator 3
  plot(a1(j),c1(j),'.k','MarkerSize',10); % articulatie corp 1 oscilator 1
  plot(a3(j),c3(j),'.k','MarkerSize',10); % articulatie corp 1 oscilator 2
  plot(a5(j),c5(j),'.k','MarkerSize',10); % articulatie corp 1 oscilator 3
  plot(a2(j),c2(j),'.b','MarkerSize',rg2); % corpul 2 oscilator 1
  plot(a4(j),c4(j),'.y','MarkerSize',rg2); % corpul 2 oscilator 2
  plot(a6(j),c6(j),'.m','MarkerSize',rg2); % corpul 2 oscilator 3
  
  axis([-Lmax Lmax -Lmax Lmax]); axis square; % cadrul grafic

 

  if RT==1 % real time(1) / slow motion(0)

    simt=toc; % actualizeaza timpul simularii cu ceasul sistemului

    text(3/5*Lmax,4/5*Lmax,['t = ',num2str(round(t(j))),' s']);

  else

    simt=simt+1e-2; % incrementeaza cu o centisecunda

    text(3/5*Lmax,4/5*Lmax,['t=',num2str(round(t(j)*100)),' cs']);

  end

  pause(1e-6); hold off

end