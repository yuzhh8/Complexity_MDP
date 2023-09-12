function ComplexityInTheoryT_Single =ST_Theory(gamma,epsilon,beta,Q_0_infty)
r_infty=2;
epsilon_iter=(1-beta)*epsilon; % 迭代精度

% 迭代时间复杂度
ComplexityInTheoryT_Single = ceil(log( (1-gamma)*epsilon_iter/(r_infty+(1+gamma)*Q_0_infty))/log(gamma));