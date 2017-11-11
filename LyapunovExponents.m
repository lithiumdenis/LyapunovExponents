close all; clear; clc;

y = zeros(12,1);            % Инициализируем y0
y(1:3) = [0.01 0.01 0.01];  % Начальные условия
y(4:12) = eye(3);           % Единичная матрица

tstart = 0;                 % Начальное значение времени
iteratetimes = 500;         % Число итераций
tincre = 1.0;               % Шаг в момент времени
sum = zeros(3,1);

% Инициализируем три показателя Ляпунова и будем хранить их для каждой
% итерации
expo = zeros(iteratetimes,3);

for i=1:iteratetimes
    %Переопределим конец
    tend = tstart + tincre;
    [~,Y] = ode45(@lorenz_ode,[tstart,tend], y);

    % Возьмем значение последнего момента, полученного интегрированием
    y = Y(size(Y,1),:);
    % Переопределим начало
    tstart = tend;
    y0 = [y(4) y(7) y(10);
          y(5) y(8) y(11);
          y(6) y(9) y(12)];
    % Ортогонализация по Граму-Шмидту
    [y0,znorm] = GS(y0);
    sum = sum + log(znorm);
    y(4:12) = y0;

    % Получим показатели Ляпунова для времени итерации
    for j=1:3
        expo(i,j) = sum(j) / tstart;
    end
end

% Выведем окончательные показатели Ляпунова
lyap = expo(length(expo),:);
for i=1:length(lyap)
    fprintf(1, 'Lambda %d = %10.5f\n', i, lyap(i));
end

% Визуализируем показатели Ляпунова по времени
i = 1:iteratetimes;
plot(i,expo(:,1),'r-',i,expo(:,2),'g-',i,expo(:,3),'b-','LineWidth',1.5)
xlabel('Время');
ylabel('Показатели Ляпунова')
title('Зависимость показателей Ляпунова от времени');
legend('\Lambda_1','\Lambda_2','\Lambda_3','Location','Best')