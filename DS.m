clear all
N = 10; % число частиц
c = 200; % жесткость
m = 1; % масса
L0 = 1; % начальная дистанция
dt = 0.01; % шаг времени
betta = 0.1;
g=9.8;
L1 = L0/(N-1);

AA = zeros(1,1);
A = [];
deltaY = zeros(1,1);
Y_N = zeros(1,1);
YF = zeros(1,1);

% Создание массива
a = 2;
for i = 1:N
    X(i) = L1*(i-1)-L0/2;
    Y(i) = a*(exp(X(i)/a)+exp(-X(i)/a))/2;
end
Y=Y-Y(1);

Vx = zeros(1,N);
Vy = zeros(1,N);

% начальные условия
Vx(1) = 0;
Vy(1) = 0;
Vy(N) = 0;
Vx(N) = 0;

X_free = 2*L0;
Y_free = 0;

figure
hold on
grid on
for st = 1:2500
    
    for i = 1:N-1
        
        L = sqrt((X(i+1)-X(i))^2+(Y(i+1)-Y(i))^2);
        e_x = (X(i+1)-X(i))/L;
        e_y = (Y(i+1)-Y(i))/L;
        
        delta_L = L-L1;
        %         if delta_L >= 0
        F_x = -c*(L-L1)*e_x;
        F_y = -c*(L-L1)*e_y-m*g;
        %         else
        %             F_x = 0;
        %             F_y = 0;
        %         end
        
        Vx(i) = Vx(i)-F_x*dt/m;
        Vx(i+1) = Vx(i+1)+F_x*dt/m;
        Vy(i) = Vy(i)-F_y*dt/m;
        Vy(i+1) = Vy(i+1)+F_y*dt/m;
    end
    
    Vx=Vx-(betta*dt/m).*Vx;
    Vy=Vy-(betta*dt/m).*Vy;
    
    Vx(1) = 0;
    Vy(1) = 0;
    
    X = X+Vx .* dt;
    Y = Y+Vy .* dt;
    X(1) = 0;
    Y(1) = 0;
    
    X_free = X_free;
    Y_free = Y_free - g*dt^2/2;
    
    A = sqrt(F_x^2+F_y^2)/m;
    AA = [AA A];
    DY = abs(Y(N)-Y_free);
    deltaY = [deltaY DY];
    
    Y_N = [Y_N (Y(N))];
    
    YF = [YF (Y_free)];
    
    if mod(st,20) == 0
        
        clf
        subplot(4,1,1)
        hold on
        grid on
        xlim([-L0 3*L0])
        ylim([-3*L0 L0/2])
        plot(X,Y,'r.-','MarkerSize',20)
        plot(X_free,Y_free,'r.-','MarkerSize',20)
        title ('Двумерная цепочка')
        xlabel('x')
        ylabel('y')
        %         pause(0.1)
        
        subplot(4,1,2)
        hold on
        grid on
        plot(AA,'r','MarkerSize',3)
        xlabel('t')
        ylabel('A')
        title ('Ускорение')
        %         pause(0.1)
        
        subplot(4,1,3)
        hold on
        grid on
        plot(deltaY,'r','MarkerSize',3)
        xlabel('t')
        ylabel('|delta Y|')
        title ('Разность координат')
        %         pause(1)
        
        subplot(4,1,4)
        hold on
        grid on
        plot(Y_N,'r','MarkerSize',3)
        plot(YF,'g','MarkerSize',3)
        xlabel('t')
        ylabel('Y')
        title ('Координата')
        legend('Последняя частица цепочки','Свободная частица');
        pause(1)
    end
    
end
