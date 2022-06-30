
%jonathan obredor
%proyecto automatica 2017
%constantes para calcular en simulink
clear all
close all
clc

%momento de inercia referido al eje de salida del tren de trasmision
Jl=0.2520%nominal
%Jl=0.1260%minimo
%Jl=0.3780 %maximo

bl=0%0.0630%; % amortiguamiento referido al eje de salida del tren de trasmision
Tl=1.57; %torque de carga
r=314.3008; %relacion de reduccion total

Jm=3.1*10^-6;%momento de inercia (motor y caja)
bm=1.5*10^-5; %coef de friccion viscosa (motor y caja)
Pp=3; %pares de polos magneticos
lambda=0.01546; %flujo concatenado en el bobinado del estator
Lq=5.8*10^-3; %inductancia del estator en cuadratura
Ld=6.6*10^-3; %inductancia del estator en el eje directo
%Lts
Rs=1.02; %resistencia del estator por fase
J_m=Jm+(Jl/r^2); %inercia referido al eje del motor
b_m=bm+(bl/r^2); %amortiguamiento referido al eje del motor
T_l=Tl/r; %carga referida al eje del motor



%%

%calcular ceros y polos del sistema

syms s
disp('G1')
[polos, ordenes, residuos]=poles(       (3*Pp*lambda)/(     (2*J_m*Lq*s^2 + 2*(J_m*Rs+b_m*Lq)*s + 2*b_m*Rs + 3*(Pp*lambda)^2    )     *s)              );
[polos]=vpa(polos, 6)
%[residuos]=vpa(residuos,6);

%ver funcion de transferencia
lazo_abierto=tf(3*Pp*lambda,    [2*J_m*Lq       2*(J_m*Rs+b_m*Lq)   2*b_m*Rs + 3*(Pp*lambda)^2    0 ]            )



disp('G2')
[polos2, ordenes2, residuos2]=poles(       ( -2*(Lq*s+Rs) )/(        (2*J_m*Lq*s^2 + 2*(J_m*Rs+b_m*Lq)*s + 2*b_m*Rs + 3*(Pp*lambda)^2    )     *s)          );

[polos2]=vpa(polos2, 6)
[residuos]=vpa(residuos,6);
num2=[0 0 -2*Lq -2*(Rs)];
sys2=tf(   num2   ,      [2*J_m*Lq       2*(J_m*Rs+b_m*Lq)   2*b_m*Rs + 3*(Pp*lambda)^2    0 ]      )


%graficar polos y ceros de las G
% figure(1)
% iopzmap(lazo_abierto,'r',sys2,'b')
% hold on
% grid on

discriminante=((J_m* Rs+b_m*Lq)/(2*J_m*Lq ))^2  -  (  ( 2*b_m*Rs+3*(Pp*lambda)^2)  /  (2*J_m*Lq )   )
%frecuencia natural de lazo abierto
frec_natural= sqrt(             (      2*b_m* Rs+3*(Pp*lambda )^2     )/(2*J_m*Lq )           )
sita_natural=(J_m*Rs+b_m*Lq)/(2*J_m*Lq * frec_natural)


%control solo proporcional para el lazo de realimentacion de corriente
p1i=-5000;
R1=-p1i*Lq;
R2=-p1i*Ld;
control_Corriente=tf(1 , [  Lq/R1 1 ]  )
[polosCorriente, ordenesCorriente, residuosCorriente]=poles(       1/(s *Lq/R1+1)         );



%%
%constantes del controlador PID
%metodo de sintonia serie
n=2.5;
wpos=800; %este es el ancho de banda
Kd=J_m*n*wpos;
Kp=J_m*n*wpos^2;
Ki=J_m*wpos^3;
sitaPID=(n-1)/2;
num3=[Kd  Kp  Ki ];
den3=[ J_m   Kd   Kp   Ki   ];
[polosPID, ordenesPID, residuosPID]=poles(     (Kd *s^2+Kp *s+Ki)/(J_m*s^3+Kd *s^2+Kp* s+Ki)     )
controlador_PID=tf(  num3 ,den3   )
figure
iopzmap(controlador_PID,'r',control_Corriente,'g',lazo_abierto,'b')


%ganancias para el observador de posicion y velocidad
Ketita=2*3200-b_m/J_m;
Kew=3200^2-Ketita*b_m/J_m;

Ketita_corr=8000;
Kew_corr=25600e3;
Ke_int=3.2768e10

consigna_pos=[0:100];

%sistema termico
Rs_cero=1.02 %resistencia del estator
Rts=55 %resistencia termica
Cts=1.091 %capacitancia termica
Ta=40 %temperatura ambiente
alfa_cu=3.9*10^-3
Ts0=40


%%
%calcular puntos de operacion
figure
iopzmap(lazo_abierto,'g')
xlim([-7000 1000])
ylim([-4000 4000])
title('traslación de los polos segun variaciones de Id')

hold on
%id positiva
idd=0
Vq=19.596
iteraciones=5;
for i=1:iteraciones
    idd=idd+i/10;



    coef1=3/2*Pp/(J_m) *(lambda+(Ld-Lq )*idd);
    coef2=-(bm)/(J_m);
    coef3=-1/(J_m)*T_l;

    coef4=-Rs/Lq;
    coef5=  -(lambda+Ld *idd )*Pp/Lq ;
    coef6=1/Lq  *Vq;

    coef7=-Rs/Ld*idd;
    coef8=(Lq*Pp)/Ld;
    coef9=1/Ld;

    A=[coef1 coef2;
        coef4 coef5];
    b=[coef3;
        coef6];
    Ab=[A b];
    R=rref(Ab);

    iqq=R(1,3)
    omegaa=R(2,3)

    vdd=(-coef7-coef8*iqq*omegaa)*Ld

    %la matriz del modelo lpv
    A=zeros(5,5);
    A(1,2)=1;
    A(2,2)=-((b_m)/(J_m));
    A(2,3)=3/2 * (Pp *lambda)/(J_m )+3/2 * (Pp *(Ld-Lq ))/(J_m ) * idd;
    A(2,4)=3/2 *Pp* (Ld-Lq ) * iqq;

    A(3,2)=-Pp/Lq *(lambda+Ld * idd );
    A(3,3)=-(Rs/Lq );
    A(3,4)=-(Ld*  Pp)/Lq *omegaa;

    A(4,2)=(Lq * Pp)/Ld *iqq;
    A(4,3)=(Lq *Pp)/Ld *omegaa;
    A(4,4)=-(Rs/Ld );
    Lls=.8e-3;
    A(5,5)=-(Rs/Lls );
    %encontrar los polos del sistema y comparar con el sistema lazo abierto
    EIG=eig(A);
    p2=plot(EIG,'xr')
    pause
end

%calculo de id negativa
iteraciones=5;
idd=0;

for i=1:iteraciones
    idd=idd-i/10;
    Vq=19.596;

    coef1=3/2*Pp/(J_m) *(lambda+(Ld-Lq )*idd);
    coef2=-(bm)/(J_m);
    coef3=-1/(J_m)*T_l;

    coef4=-Rs/Lq;
    coef5=  -(lambda+Ld *idd )*Pp/Lq ;
    coef6=1/Lq  *Vq;

    coef7=-Rs/Ld*idd;
    coef8=(Lq*Pp)/Ld;
    coef9=1/Ld;

    A=[coef1 coef2;
        coef4 coef5];
    b=[coef3;
        coef6];
    Ab=[A b];
    RR=rref(Ab);

    iqq=RR(1,3)
    omegaa=RR(2,3)

    vdd=(-coef7-coef8*iqq*omegaa)*Ld

    %la matriz del modelo lpv
    A=zeros(5,5);
    A(1,2)=1;
    A(2,2)=-((b_m)/(J_m));
    A(2,3)=3/2 * (Pp *lambda)/(J_m )+3/2 * (Pp *(Ld-Lq ))/(J_m ) * idd;
    A(2,4)=3/2 *Pp* (Ld-Lq ) * iqq;

    A(3,2)=-Pp/Lq *(lambda+Ld * idd );
    A(3,3)=-(Rs/Lq );
    A(3,4)=-(Ld*  Pp)/Lq *omegaa;

    A(4,2)=(Lq * Pp)/Ld *iqq;
    A(4,3)=(Lq *Pp)/Ld *omegaa;
    A(4,4)=-(Rs/Ld );
    Lls=.8e-3;
    A(5,5)=-(Rs/Lls );
    %encontrar los polos del sistema y comparar con el sistema lazo abierto
    EIG1=eig(A);
    p3=plot(EIG1,'xb')
    pause
end

legend([p2 p3],'LPV con Id>0','LPV con Id<0')