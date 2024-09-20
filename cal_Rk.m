function Rk = cal_Rk(N, d, Rcap, fill_efficiency)
 % calculateOuterRadius: 计算同心圆柱体的外圆半径
    % 输入参数:
    % N - 小圆球数量
    % d - 小圆球直径 (单位: 米)
    % r_inner - 内圆柱半径 (单位: 米)
    % fill_efficiency - 填充效率 (介于 0 到 1 之间)
    % 输出:
    % r_outer - 同心圆柱体外圆半径 (单位: 米)

    % 计算小球的半径
    r_small = d / 2;
    
    % 计算小球的体积 (单位: 立方米)
    V_small = (4/3) * pi * r_small^3;
    
    % 计算总有效体积 (单位: 立方米)
    V_effective = N * V_small;
    
    % 考虑填充效率，计算实际需要的总体积 (单位: 立方米)
    V_total = V_effective / fill_efficiency;
    
    % 计算内圆柱体的体积 (单位: 立方米)
    V_inner = pi * Rcap^2 * Rcap;
    
    % 设定方程并求解外圆半径 (单位: 米)
    % V_cylinder = V_total
    % pi * r_outer^3 = V_total + V_inner
    V_cylinder = V_total + V_inner;
    Rk = (V_cylinder / pi)^(1/3);
end


