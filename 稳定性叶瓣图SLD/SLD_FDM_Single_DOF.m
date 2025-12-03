close all; clear all; clc;
tic

% 铣削工艺参数
N = 2;                   % 刀具齿数
Kt = 6e8;           % 切向切削力系数 (N/m^2)
Kn = 2e8;            % 法向切削力系数 (N/m^2)
aD = 0.05;               % 径向切深比

up_or_down = -1;        % 1: 逆铣; -1: 顺铣
if up_or_down == 1      % 逆铣
    fist = 0;           % 切入角 (rad)
    fiex = acos(1-2*aD); % 切出角 (rad)
elseif up_or_down == -1 % 顺铣 
    fist = acos(2*aD-1); % 切入角 (rad)
    fiex = pi;          % 切出角 (rad)
end

% 单自由度模态参数
modal_freq = 922 * 2 * pi;     % 固有频率 (rad/s)
modal_damping = 0.011;           % 阻尼比
modal_mass = 0.03993;     % 模态质量 (kg)


% ========== 仿真参数 ==========
stx = 400;              % 主轴转速步数
sty = 200;              % 切削深度步数
w_st = 0;               % 起始切削深度 (m)
w_fi = 0.01;            % 最终切削深度 (m)  
o_st = 0.5e4;           % 起始主轴转速 (rpm)
o_fi = 2.5e4;           % 最终主轴转速 (rpm)
m = 40;                 % 离散区间数

% ========== 初始化离散映射矩阵 ==========
D = zeros(m+2, m+2);
d = ones(m+1, 1);
d(1:2) = 0;
D = D + diag(d, -1);
D(3, 1) = 1;

% ========== 切削力系数离散化 ==========
h = zeros(m+1, 1);
dtr = 2*pi/N/m;         % 角度步长

for i = 1:m+1
    for j = 1:N         % 遍历每个刀齿
        fi = i*dtr + (j-1)*2*pi/N;
        
        % 判断刀齿是否参与切削
        if (fi >= fist) && (fi <= fiex)
            g = 1;
        else
            g = 0;
        end
        
        % 计算时变切削力系数
        h(i) = h(i) + g*(Kt*cos(fi) + Kn*sin(fi))*sin(fi);
    end
end

% ========== 核心算法部分 ==========
A0 = [-modal_damping*modal_freq, 1/modal_mass;
      modal_mass*((modal_damping*modal_freq)^2 - modal_freq^2), -modal_damping*modal_freq];
I = eye(size(A0));
invA0 = inv(A0);

% 初始化结果矩阵
ss = zeros(stx+1, sty+1);
dc = zeros(stx+1, sty+1);
ei = zeros(stx+1, sty+1);

% 主循环
for x = 1:stx+1
    o = o_st + (x-1)*(o_fi-o_st)/stx;  % 当前主轴转速
    tau = 60/o/N;                      % 时滞周期
    dt = tau/m;                        % 时间步长
    
    % 计算状态转移矩阵
    Fi0 = expm(A0*dt);
    Fi1 = invA0*(Fi0 - I);
    Fi2 = invA0*(Fi0*dt - Fi1);
    Fi3 = invA0*(Fi0*dt^2 - 2*Fi2);
    Fi4 = invA0*(Fi0*dt^3 - 3*Fi3);
    Fi5 = invA0*(Fi0*dt^4 - 4*Fi4);
    
    % 计算加权系数矩阵
    Q_minus_2_1k = (Fi2/(3*dt) - 5*Fi3/(6*dt^2) + 2*Fi4/(3*dt^3) - Fi5/(6*dt^4));
    Q_minus_2_0k = (Fi3/(3*dt^2) - Fi4/(2*dt^3) + Fi5/(6*dt^4));
    Q_minus_1_1k = (-3*Fi2/(2*dt) + 7*Fi3/(2*dt^2) - 5*Fi4/(2*dt^3) + Fi5/(2*dt^4));
    Q_minus_1_0k = (-3*Fi3/(2*dt^2) + 2*Fi4/(dt^3) - Fi5/(2*dt^4));
    Q_0_1k = (3*Fi2/(dt) - 11*Fi3/(2*dt^2) + 3*Fi4/(dt^3) - Fi5/(2*dt^4));
    Q_0_0k = (3*Fi3/(dt^2) - 5*Fi4/(2*dt^3) + Fi5/(2*dt^4));
    Q_1_1k = (Fi1 - 17*Fi2/(6*dt) + 17*Fi3/(6*dt^2) - 7*Fi4/(6*dt^3) + Fi5/(6*dt^4));
    Q_1_0k = (Fi2/dt - 11*Fi3/(6*dt^2) + Fi4/(dt^3) - Fi5/(6*dt^4));
    
    H_0_1j = (Fi2/(3*dt) + Fi3/(6*dt^2) - Fi4/(3*dt^3) - Fi5/(6*dt^4));
    H_0_0j = (Fi3/(3*dt^2) + Fi4/(2*dt^3) + Fi5/(6*dt^4));
    H_minus_1_1j = (Fi1 - Fi2/(2*dt) - 3*Fi3/(2*dt^2) + Fi4/(2*dt^3) + Fi5/(2*dt^4));
    H_minus_1_0j = (Fi2/dt + Fi3/(2*dt^2) - Fi4/(dt^3) - Fi5/(2*dt^4));
    H_minus_2_1j = (-Fi2/(dt) + 3*Fi3/(2*dt^2) - Fi5/(2*dt^4));
    H_minus_2_0j = (-Fi3/(dt^2) + Fi4/(2*dt^3) + Fi5/(2*dt^4));
    H_minus_3_1j = (Fi2/(6*dt) - Fi3/(6*dt^2) - Fi4/(6*dt^3) + Fi5/(6*dt^4));
    H_minus_3_0j = (Fi3/(6*dt^2) - Fi5/(6*dt^4));
    
    for y = 1:sty+1
        w = w_st + (y-1)*(w_fi-w_st)/sty;  % 当前切削深度
        Fi = eye(m+2);
        
        for i = 1:m
            % 构建时变切削力矩阵
            hi0 = w * h(i);
            hi1 = w * h(i+1);
            
            % 单自由度系统矩阵
            Ak0 = [0, 0; -hi0, 0];
            Ak1 = [0, 0; -hi1, 0];
            Bk0 = [0, 0; hi0, 0];
            Bk1 = [0, 0; hi1, 0];
            
            % 计算加权矩阵
            Q_minus_2 = Q_minus_2_1k*Ak1 + Q_minus_2_0k*Ak0;
            Q_minus_1 = Q_minus_1_1k*Ak1 + Q_minus_1_0k*Ak0;
            Q_0 = Q_0_1k*Ak1 + Q_0_0k*Ak0;
            Q_1 = Q_1_1k*Ak1 + Q_1_0k*Ak0;
            
            H_0 = H_0_1j*(-Bk1) + H_0_0j*(-Bk0);
            H_minus_1 = H_minus_1_1j*(-Bk1) + H_minus_1_0j*(-Bk0);
            H_minus_2 = H_minus_2_1j*(-Bk1) + H_minus_2_0j*(-Bk0);
            H_minus_3 = H_minus_3_1j*(-Bk1) + H_minus_3_0j*(-Bk0);
            
            inv0fimk1 = inv(Q_1 - I);
            
            % 更新离散映射矩阵D
            D(1:2, 1:2) = -inv0fimk1*(Q_0 + Fi0);
            D(1:2, 3) = -inv0fimk1*Q_minus_1(:, 1);
            D(1:2, 4) = -inv0fimk1*Q_minus_2(:, 1);
            
            D(1:2, m-1) = inv0fimk1*H_minus_3(:, 1);
            D(1:2, m) = inv0fimk1*H_minus_2(:, 1);
            D(1:2, m+1) = inv0fimk1*H_minus_1(:, 1);
            D(1:2, m+2) = inv0fimk1*H_0(:, 1);
            
            Fi = D * Fi;
        end
        
        ss(x,y) = o;
        dc(x,y) = w;
        ei(x,y) = max(abs(eig(Fi)));
    end
    
    fprintf('进度: %d/%d, 转速: %.1f rpm\n', x, stx+1, o);
end

toc

% ========== 结果可视化 ==========
figure;
contour(ss, dc, ei, [1, 1], 'k');
xlabel('主轴转速 (rpm)');
ylabel('轴向切深 (m)');
title('铣削稳定性叶瓣图');
grid on;