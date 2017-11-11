close all; clear; clc;

y = zeros(12,1);            % �������������� y0
y(1:3) = [0.01 0.01 0.01];  % ��������� �������
y(4:12) = eye(3);           % ��������� �������

tstart = 0;                 % ��������� �������� �������
iteratetimes = 500;         % ����� ��������
tincre = 1.0;               % ��� � ������ �������
sum = zeros(3,1);

% �������������� ��� ���������� �������� � ����� ������� �� ��� ������
% ��������
expo = zeros(iteratetimes,3);

for i=1:iteratetimes
    %������������� �����
    tend = tstart + tincre;
    [~,Y] = ode45(@lorenz_ode,[tstart,tend], y);

    % ������� �������� ���������� �������, ����������� ���������������
    y = Y(size(Y,1),:);
    % ������������� ������
    tstart = tend;
    y0 = [y(4) y(7) y(10);
          y(5) y(8) y(11);
          y(6) y(9) y(12)];
    % ��������������� �� �����-������
    [y0,znorm] = GS(y0);
    sum = sum + log(znorm);
    y(4:12) = y0;

    % ������� ���������� �������� ��� ������� ��������
    for j=1:3
        expo(i,j) = sum(j) / tstart;
    end
end

% ������� ������������� ���������� ��������
lyap = expo(length(expo),:);
for i=1:length(lyap)
    fprintf(1, 'Lambda %d = %10.5f\n', i, lyap(i));
end

% ������������� ���������� �������� �� �������
i = 1:iteratetimes;
plot(i,expo(:,1),'r-',i,expo(:,2),'g-',i,expo(:,3),'b-','LineWidth',1.5)
xlabel('�����');
ylabel('���������� ��������')
title('����������� ����������� �������� �� �������');
legend('\Lambda_1','\Lambda_2','\Lambda_3','Location','Best')