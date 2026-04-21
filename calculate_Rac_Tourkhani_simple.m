% ---------------------------------------------------------
% This script implements the analytical AC resistance model 
% proposed in the following paper:
% [1] F. Tourkhani and P. Viarouge, “Accurate analytical model of winding losses in round Litz wire windings,” IEEE Transactions on Magnetics, vol. 37, no. 1, pp. 538–543, Jan. 2001, doi: 10.1109/20.914375.

% ---------------------------------------------------------
% 计算变压器利兹线交流电阻（考虑集肤效应和邻近效应）

clear; clc; close all;

% --- 输入参数 (参考施耐德样机参数) ---
f_s = 70e3;             % 开关频率 70 kHz
rho = 1.724e-8;          % 铜的电阻率 (室温下，单位: Ω·m)
mu = 4 * pi * 1e-7;     % 真空磁导率 (H/m)

N = 20;                 % 原边绕组匝数
N0 = 250;               % 原边单根导线包含的股线数
d_o = 0.2e-3;           % 单根股线直径 0.2mm
m = 2;                  % 绕组层数 (假定示例为2层)
beta = 0.5776;             % 填充系数 (假定示例)第一次原边0.5115 副边0.5103  第二次原边0.5779 副边0.5776
l_single = 160e-3;       % 每匝绕组长度 (m)

% --- 1. 计算趋肤深度和基础电阻 ---
delta = sqrt(rho / (pi * mu * f_s)); 
zeta = d_o / delta;     

% 计算基础电阻 Rb
Rb = l_single* (N * rho) / (N0 * pi * delta^2); 
Rdc = l_single * (4 * N * rho) / (N0 * pi * d_o^2); % 直流电阻
Awdg = 1/(4/(pi * d_o^2 *N0));

% 计算最优 zeta (zeta_op)
zeta_op = 4 * (3 / (1 + (pi^2 * N0 * beta / 4) * (16 * m^2 - 1 + 24/pi^2)))^(1/4);
fprintf('最优参数 zeta_op = %g\n', zeta_op);

% --- 2. 泰勒展开法近似计算 ---
% 根据文中给出的公式：
psi1_taylor = 2 * sqrt(2) * (1/zeta + (1/(3*2^8))*zeta^3 - (1/(3*2^14))*zeta^5);
psi2_taylor = (1/sqrt(2)) * (-(1/2^5)*zeta^3 + (1/2^12)*zeta^7);

% 计算交流电阻系数 Kr Kd
term2_taylor = (pi^2 * N0 * beta / 24) * (16*m^2 - 1 + 24/pi^2) * psi2_taylor;
Kr_taylor = (sqrt(2) / zeta) * (psi1_taylor - term2_taylor);
Kd_taylor = Kr_taylor * zeta^2 / 4; % 交流电阻与直流电阻的比值

Rac_taylor = Rb * Kr_taylor;

% --- 3. 原始开尔文函数(Kelvin functions) 精确计算 ---
% Kelvin 函数与 Bessel 函数的关系:
% ber_n(x) + j*bei_n(x) = j^(-n) * J_n(x * exp(j*3*pi/4))
x = zeta / sqrt(2);
z = x * exp(3i*pi/4);

% 0阶
J0 = besselj(0, z);
ber = real(J0);
bei = imag(J0);

% J0的导数 (用于计算 ber' 和 bei')
% J0'(z) = -J1(z)
% d/dx [ber(x) + j*bei(x)] = d/dx [J0(x*exp(3j*pi/4))] = -J1(x*exp(3j*pi/4)) * exp(3j*pi/4)
J1_term = -besselj(1, z) * exp(3i*pi/4);
ber_p = real(J1_term); % ber'
bei_p = imag(J1_term); % bei'

% 2阶
%J2 = (1 / 1i^2) * besselj(2, z); % j^(-2) * J_2(...)
J2 = besselj(2, z); % 注意：这里直接使用 J2 的值，因为在计算 psi2_exact 时需要 ber^2 + bei^2 的分母，而 ber 和 bei 已经是 J0 的实部和虚部了。
ber2 = real(J2);
bei2 = imag(J2);

% 组装精确的 psi1 和 psi2
psi1_exact = (ber * bei_p - bei * ber_p) / (ber_p^2 + bei_p^2);
psi2_exact = (ber2 * ber_p + bei2 * bei_p) / (ber^2 + bei^2);

term2_exact = (pi^2 * N0 * beta / 24) * (16*m^2 - 1 + 24/pi^2) * psi2_exact;
Kr_exact = (sqrt(2) / zeta) * (psi1_exact - term2_exact);
kd_exact = Kr_exact * zeta^2 / 4; % 交流电阻与直流电阻的比值

Rac_exact = Rb * Kr_exact;

% --- 输出对比结果 ---
fprintf('=== 交流电阻计算结果 ===\n');
fprintf('集肤深度 delta = %g mm\n', delta * 1000);
fprintf('参数 zeta = %g\n', zeta);
fprintf('基础交流阻抗 Rb = %g Ω\n', Rb);
fprintf('直流电阻 Rdc = %g Ω\n\n', Rdc);

fprintf('【泰勒展开近似法】:\n');
fprintf('psi1_taylor = %g, psi2_taylor = %g\n', psi1_taylor, psi2_taylor);
fprintf('Kr_taylor = %g\n', Kr_taylor);
fprintf('Rac_taylor = %g Ω\n\n', Rac_taylor);

fprintf('【原始 psi 函数法】:\n');
fprintf('psi1_exact = %g, psi2_exact = %g\n', psi1_exact, psi2_exact);
fprintf('Kr_exact = %g\n', Kr_exact);
fprintf('Rac_exact = %g Ω\n\n', Rac_exact);

fprintf('相对误差 = %g %%%\n', abs(Rac_exact - Rac_taylor)/Rac_exact * 100);
fprintf('导线有效截面积： %g m^2 \n',Awdg);

% --- 附加：绘制 K_r 和 K_d 随 zeta 变化的曲线 ---
zeta_vec = linspace(0.01, 3, 500);
Kr_taylor_vec = zeros(size(zeta_vec));
Kd_taylor_vec = zeros(size(zeta_vec));
Kr_exact_vec = zeros(size(zeta_vec));
Kd_exact_vec = zeros(size(zeta_vec));

for k = 1:length(zeta_vec)
    z_val = zeta_vec(k);
    
    % --- 泰勒展开法 ---
    psi1_t = 2 * sqrt(2) * (1/z_val + (1/(3*2^8))*z_val^3 - (1/(3*2^14))*z_val^5);
    psi2_t = (1/sqrt(2)) * (-(1/2^5)*z_val^3 + (1/2^12)*z_val^7);
    
    term2_t = (pi^2 * N0 * beta / 24) * (16*m^2 - 1 + 24/pi^2) * psi2_t;
    Kr_t = (sqrt(2) / z_val) * (psi1_t - term2_t);
    Kd_t = Kr_t * z_val^2 / 4;
    
    Kr_taylor_vec(k) = Kr_t;
    Kd_taylor_vec(k) = Kd_t;
    
    % --- 原始精确法 ---
    x_val = z_val / sqrt(2);
    z_cplx = x_val * exp(3i*pi/4);
    
    J0_val = besselj(0, z_cplx);
    ber_v = real(J0_val);
    bei_v = imag(J0_val);
    
    J1_term_v = -besselj(1, z_cplx) * exp(3i*pi/4);
    ber_p_v = real(J1_term_v);
    bei_p_v = imag(J1_term_v);
    
    J2_val = besselj(2, z_cplx);
    ber2_v = real(J2_val);
    bei2_v = imag(J2_val);
    
    psi1_e = (ber_v * bei_p_v - bei_v * ber_p_v) / (ber_p_v^2 + bei_p_v^2);
    psi2_e = (ber2_v * ber_p_v + bei2_v * bei_p_v) / (ber_v^2 + bei_v^2);
    
    term2_e = (pi^2 * N0 * beta / 24) * (16*m^2 - 1 + 24/pi^2) * psi2_e;
    Kr_e = (sqrt(2) / z_val) * (psi1_e - term2_e);
    Kd_e = Kr_e * z_val^2 / 4;
    
    Kr_exact_vec(k) = Kr_e;
    Kd_exact_vec(k) = Kd_e;
end

% 绘图1：Kr 随 zeta 变化
figure('Position', [100, 400, 400, 300]);
plot(zeta_vec, Kr_exact_vec, 'b-', 'LineWidth', 1.5); hold on;
plot(zeta_vec, Kr_taylor_vec, 'r--', 'LineWidth', 1.5);
plot(zeta, Kr_taylor, 'o', 'Color', [0.85 0.325 0.098], 'MarkerSize', 6, 'MarkerFaceColor', [0.85 0.325 0.098]); % 标注当前zeta
text(zeta + 0.15, Kr_taylor, sprintf('当前 \\zeta=%.3f, K_r=%.2f', zeta, Kr_taylor), 'FontSize', 11, 'Color', [0.85 0.325 0.098]);

% 标注最优 zeta_op (基于原方程精确曲线插值)
Kr_op_exact = interp1(zeta_vec, Kr_exact_vec, zeta_op, 'linear', 'extrap');
plot(zeta_op, Kr_op_exact, 'p', 'Color', [0.4660 0.6740 0.1880], 'MarkerSize', 8, 'MarkerFaceColor', [0.4660 0.6740 0.1880]); 
text(zeta_op + 0.15, Kr_op_exact-20, sprintf('最优 \\zeta_{op}=%.3f', zeta_op), 'FontSize', 11, 'Color', [0.4660 0.6740 0.1880]);

xlabel('\zeta', 'FontSize', 12);
ylabel('K_r', 'FontSize', 12);
%title('K_r 随 \zeta 变化曲线', 'FontSize', 14);
%legend('原方程 (Kelvin 函数)', '泰勒展开', '当前 \zeta 取值', '最优 \zeta_{op}', 'Location', 'Best', 'FontSize', 11);
axis([0 3 0 800]);
grid on;

% 绘图2：Kd 随 zeta 变化
figure('Position', [750, 400, 400, 300]);
plot(zeta_vec, Kd_exact_vec, 'b-', 'LineWidth', 1.5); hold on;
plot(zeta_vec, Kd_taylor_vec, 'r--', 'LineWidth', 1.5);
plot(zeta, Kd_taylor, 'o', 'Color', [0.85 0.325 0.098], 'MarkerSize', 6, 'MarkerFaceColor', [0.85 0.325 0.098]); % 标注当前zeta
text(zeta + 0.2, Kd_taylor + 3, sprintf('当前 \\zeta=%.3f, K_d=%.2f', zeta, Kd_taylor), 'FontSize', 11, 'Color', [0.85 0.325 0.098]);

% 标注最优 zeta_op (基于原方程精确曲线插值)
Kd_op_exact = interp1(zeta_vec, Kd_exact_vec, zeta_op, 'linear', 'extrap');
plot(zeta_op, Kd_op_exact, 'p', 'Color', [0.4660 0.6740 0.1880], 'MarkerSize', 8, 'MarkerFaceColor', [0.4660 0.6740 0.1880]); 
text(zeta_op + 0.4, Kd_op_exact + 3, sprintf('最优 \\zeta_{op}=%.3f', zeta_op), 'FontSize', 11, 'Color', [0.4660 0.6740 0.1880]);

xlabel('\zeta', 'FontSize', 12);
ylabel('K_d', 'FontSize', 12);
%title('K_d 随 \zeta 变化曲线', 'FontSize', 14);
legend('原方程 (Kelvin 函数)', '泰勒展开', '当前 \zeta 取值', '最优 \zeta_{op}', 'Location', 'Best', 'FontSize', 11);
axis([0 2 0 100]);
grid on;
