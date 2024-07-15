clc
clear
close all

myVideo = VideoWriter('myVideoFile'); %Crea el archivo del video
myVideo.FrameRate = 60;  %Selecciona los fps
open(myVideo) %Abre el archivo para editarlo

posicion=[0,0]; %Posicion inicial de mi carga
posicion_ret=[0,0]; %Posicion retardada
x=linspace(-30,30,61); %Hago mi vector x para manejar con mayor facilidad y rapidez el calculo del tiempo y radio retardado
y=linspace(-30,30,61); %Hago mi vector y para manejar con mayor facilidad y rapidez el calculo del tiempo y radio retardado
[X,Y]=meshgrid(x,y); %Creo mi espacio vectorial
q=5; %Valor de la carga
w=4; %Frecuencia angular
c=15; %Velocidad de la luz
A=5; %Amplitud
t=0:0.01:30; %Vector tiempo
k=1; %Esto es 4*pi*epsilon0
V=zeros(length(x),length(y)); %Espacio de mi campo escalar potencial electrico

pot_a_x=zeros(length(x),length(y)); %Espacio de mi campo vectorial potencial magnetico en x
pot_a_y=zeros(length(x),length(y)); %Espacio de mi campo vectorial potencial magnetico en y
a_pot_x=cell(length(t),1); %Tensor que almacenará mi campo vectorial potencial magnetico en x
a_pot_y=cell(length(t),1); %Tensor que almacenará mi campo vectorial potencial magnetico en y
a_posicion_ret=cell(length(y),1); %Tensor que almacenará mis posiciones retardadas
tr=zeros(length(t),1); %Lista que almacenará mis tiempos retardados

for i=1:length(t) %Ciclo que calculará la posicion de la carga, potencial magnetico, gradiente del potencial electrico, y el campo electrico

    posicion(2)=A*cos(w*t(i)); %Aqui declaro el movimiento que hará mi carga

    for j=1:length(x) %Estos dos ciclos serán los encargados de calcular mi campo escalar potencial electrico, el tiempo retardado y el radio retardado
        for k=1:length(y)
            if i==1 %Se declara que en la primera iteración, la posicion retardada será cero, y se hace el calculo del tiempo retardado

                posicion_ret(2)=0;
                a_posicion_ret{i}=posicion_ret(2);
                tr(i)=t(i)-(sqrt(((x(j))^2+(y(k)-0)^2)))/c;

            else %Se calcula tiempo retardado y posicion retardada

                tr(i)=t(i)-(sqrt(((x(j))^2+(y(k)-a_posicion_ret{i-1})^2)))/c;
                posicion_ret(2)=A*cos(w*tr(i));
                a_posicion_ret{i}=posicion_ret(2);

            end

            rad_ret=[x(j),(y(k)-posicion_ret(2))]; %Se calcula el radio retardado
            velocidad_ret=[0,-w*A*sin(w*tr(i))]; %Se calcula la velocidad retardada (se deriva la funcion de movimiento y se evalúa en el tiempo retardado)
            V(k,j)=k*(q*c)/(sqrt(rad_ret(1)^2+rad_ret(2)^2)*c-dot(rad_ret,velocidad_ret)); %Se calcula el potencial electrico

        end
    end

    pot_a_x=(velocidad_ret(1)/c^2).*V; %Calculo de los potenciales magneticos
    pot_a_y=(velocidad_ret(2)/c^2).*V;
    a_pot_x{i}=pot_a_x; %Se almacena en el tensor el potencial magnetico
    a_pot_y{i}=pot_a_y;

    %En este if se deriva el potencial magnetico
    if i == 1
        t_pot_x=a_pot_x{i};
        t_pot_y=a_pot_y{i};
    else
        t_pot_x=(a_pot_x{i}-a_pot_x{i-1});
        t_pot_y=(a_pot_y{i}-a_pot_y{i-1});
    end


    [gradVx,gradVy]=gradient(V);%Se calcula el gradiente del potencial electrico
    Ex=-gradVx-t_pot_x; %Se calcula el campo electrico
    Ey=-gradVy-t_pot_y;
    mag=sqrt((Ex).^2+(Ey).^2); %Se busca la magnitud del campo electrico

    quiver(x,y,Ex./mag,Ey./mag,'ShowArrowHead','off');    axis([-30,30,-30,30]); drawnow;   %Se grafica
    hold on
    plot(posicion(1),posicion(2),'.');
    hold off

    pause(0.00001)
    frame = getframe(gcf);
    writeVideo(myVideo, frame);
 

end
close(myVideo)