p = 0.75; choice = 2; year = 0;

%% 定义
j = sym('j');
t = sym('t');
mu1 = 0.5259; mu2 = 0.5465; mu3 = 0.5086; % 树干在生物量中的占比
eta1 = 0.1764; eta2 = 0.1632; eta3 = 0.1801;% 树枝在生物量中的占比
rho1 = 0.596; rho2 = 0.358; rho3 = 0.475; % 树干密度
average_D1 = 16.86-1.25; average_D2 = 23.01-1.42; average_D3 = 17.21-1.15; % 平均胸径
N1 = 12600; N2 = 11400; N3 = 8700; % 数目总数
T1 = 35; T2 = 40; T3 = 40; % 0时刻树木最大年龄
p1 = p; p2 = p; p3 = p; % 砍伐率
P1 = 40+year; P2 = 55+year; P3 = 60+year; % 轮伐期
q = 0.334; % 填埋率
r = 0.0912; % 分解速率
Kc = 0.5; % 生物量-固碳量转换系数
L1 = 15; L2 = 2; L3 = 1; % 产品使用寿命
if choice == 1
    alpha1 = 0.8*mu1; beta1 = 0; gamma1 = 0.2*mu1+eta1; % 产品策略
    alpha2 = 0.8*mu2; beta2 = 0; gamma2 = 0.2*mu2+eta2;
    alpha3 = 0.8*mu3; beta3 = 0; gamma3 = 0.2*mu3+eta3;
else
    alpha1 = 0; beta1 = 0.8*mu1; gamma1 = 0.2*mu1+eta1; % 产品策略
    alpha2 = 0; beta2 = 0.8*mu2; gamma2 = 0.2*mu2+eta2;
    alpha3 = 0; beta3 = 0.8*mu3; gamma3 = 0.2*mu3+eta3;
end
Pr1 = 5300; Pr2 = 3300; Pr3 = 900;
%% 胸径Dij数据拟合
fit_func = fittype('c/(1+a*exp(-b*t))', 'independent', 't', 'coefficients', {'a','b','c'});  % 拟合函数类型a*t/(a+b*t)
raw_D1 = [0 4.650302778 9.912501403 13.17246799 14.9256427 15.89324511 16.69899162 17.081377]';
j1 = [0:5:35]';
cfun_D1 = fit(j1, raw_D1, fit_func);
%     aaa=1:100;
%     bbb=cfun_D1(aaa);
%     plot(aaa,bbb,j1,raw_D1);

raw_D2 = [0 1.569559733 9.211390006 14.32563346 15.50488341 18.24382062 20.54041114 22.08200435 23.316527]';
j2 = [0:5:40]';
cfun_D2 = fit(j2, raw_D2, fit_func);

raw_D3 = [0 2.488493109 5.155219209 8.356016932 12.1941016 14.74688142 15.89991481 16.88967252 17.45319651]';
j3 = [0:5:40]';
cfun_D3 = fit(j3, raw_D3, fit_func);
%% 高度Hij数据拟合
raw_H1 = [0 6.21600482 12.83876406 14.10135112 14.79148078 14.95318602 15.40840566 15.64705269]';
j1 = [0:5:35]';
cfun_H1 = fit(j1, raw_H1, fit_func);

raw_H2 = [0 2.506633912 8.04676184 13.63045922 15.01685999 17.42665037 18.45308173 20.34652777 21.42269899]';
j2 = [0:5:40]';
cfun_H2 = fit(j2, raw_H2, fit_func);

raw_H3 = [0 4.797796035 7.280907079 10.66437618 15.37407734 19.12513471 21.86314442 23.28158495 23.70719823]';
j3 = [0:5:40]';
cfun_H3 = fit(j3, raw_H3, fit_func);
%% 计算生物量BMij,t
% 单位换算：高度m*胸径cm2*密度g/cm3=m*g/cm=100g=0.1kg=0.00001t
BM1 = @(j) 0.0001*(pi*(cfun_H1.c/(1+cfun_H1.a*exp(-cfun_H1.b*j)))*(cfun_D1.c/(1+cfun_D1.a*exp(-cfun_D1.b*j)))^2*rho1)/(4*mu1);
BM2 = @(j) 0.0001*(pi*(cfun_H2.c/(1+cfun_H2.a*exp(-cfun_H2.b*j)))*(cfun_D2.c/(1+cfun_D2.a*exp(-cfun_D2.b*j)))^2*rho2)/(4*mu2);
BM3 = @(j) 0.0001*(pi*(cfun_H3.c/(1+cfun_H3.a*exp(-cfun_H3.b*j)))*(cfun_D3.c/(1+cfun_D3.a*exp(-cfun_D3.b*j)))^2*rho3)/(4*mu3);
%% 计算死亡率mi
%     syms m1 m2 m3;
%     eqn1 = average_D1 * N1==symsum((cfun_D1.c/(1+cfun_D1.a*exp(-cfun_D1.b*j)))*N1*m1*(1-m1)^(j-1),j,1,T1-1)+(cfun_D1.c/(1+cfun_D1.a*exp(-cfun_D1.b*T1)))*N1*(1-m1)^(T1-1);
%     m1 = vpasolve(eqn1,m1,0.00558);
%     m1 = m1(imag(m1)==0 & real(m1)<1);
% 
%     eqn2 = average_D2 * N2==symsum((cfun_D2.c/(1+cfun_D2.a*exp(-cfun_D2.b*j)))*N2*m2*(1-m2)^(j-1),j,1,T2-1)+(cfun_D2.c/(1+cfun_D2.a*exp(-cfun_D2.b*T2)))*N2*(1-m2)^(T2-1);
%     m2 = vpasolve(eqn2,m2,0.00083);
%     m2 = m2(imag(m2)==0 & real(m2)<1);
% 
%     eqn3 = average_D3 * N3==symsum((cfun_D3.c/(1+cfun_D3.a*exp(-cfun_D3.b*j)))*N3*m3*(1-m3)^(j-1),j,1,T3-1)+(cfun_D3.c/(1+cfun_D3.a*exp(-cfun_D3.b*T3)))*N3*(1-m3)^(T3-1);
%     m3 = vpasolve(eqn3,m3,0.00447);
%     m3 = m3(imag(m3)==0 & real(m3)<1);
m1 = 0.00558; m2 = 0.00083; m3 = 0.00447;
%% CT
n1 = zeros(200,200);
CT1 = zeros(1,200);
for j = 1:T1
    if j == T1
        n1(j,0+1) = N1*(1-m1)^(j-1);
    else
        n1(j,0+1) = N1*m1*(1-m1)^(j-1);
    end
end
for t = 1:100
    for j = 1:T1+t
        if j == 1
            n1(j,t+1) = N1*m1+sum(n1(P1-1:T1+t,t))*(1-m1)*p1;
        elseif j >= P1
            n1(j,t+1) = n1(j-1,t)*(1-m1)*(1-p1);
        else
            n1(j,t+1) = n1(j-1,t)*(1-m1);
        end
    end
    for j = 1:T1+t
        CT1(1,t) = CT1(1,t) + Kc*BM1(j)*n1(j,t+1);
    end
%     CT1(1,t) = sum(Kc*BM1(1:j)*n1(1:j,t+1));
end

n2 = zeros(200,200);
CT2 = zeros(1,200);
for j = 1:T2
    if j == T2
        n2(j,0+1) = N2*(1-m2)^(j-1);
    else
        n2(j,0+1) = N2*m2*(1-m2)^(j-1);
    end
end
for t = 1:100
    for j = 1:T2+t
        if j == 1
            n2(j,t+1) = N2*m2+sum(n2(P2-1:T2+t,t))*(1-m2)*p2;
        elseif j >= P2
            n2(j,t+1) = n2(j-1,t)*(1-m2)*(1-p2);
        else
            n2(j,t+1) = n2(j-1,t)*(1-m2);
        end
    end
    for j = 1:T2+t
        CT2(1,t) = CT2(1,t) + Kc*BM2(j)*n2(j,t+1);
    end
%     CT2(1,t) = sum(Kc*BM2(1:j)*n2(1:j,t+1));
end

n3 = zeros(200,200);
CT3 = zeros(1,200);
for j = 1:T3
    if j == T3
        n3(j,0+1) = N3*(1-m3)^(j-1);
    else
        n3(j,0+1) = N3*m3*(1-m3)^(j-1);
    end
end
for t = 1:100
    for j = 1:T3+t
        if j == 1
            n3(j,t+1) = N3*m3+sum(n3(P3-1:T3+t,t))*(1-m3)*p3;
        elseif j >= P3
            n3(j,t+1) = n3(j-1,t)*(1-m3)*(1-p3);
        else
            n3(j,t+1) = n3(j-1,t)*(1-m3);
        end
    end
    for j = 1:T3+t
        CT3(1,t) = CT3(1,t) + Kc*BM3(j)*n3(j,t+1);
    end
%     CT3(1,t) = sum(Kc*BM3(1:j)*n3(1:j,t+1));
end

CT = CT1 + CT3 + CT3;
%% CP
CP1 = zeros(1,200);
H1 = zeros(1,200);
for t = 1:100
    for j = P1:t+T1
        H1(1,t+L1) = H1(1,t+L1)+n1(j-1,t)*(1-m1)*p1*BM1(j);
    end
    CP1(1,t) = Kc*(sum(alpha1*H1(1,t+L1-L1:t+L1))+sum(beta1*H1(1,t+L1-L2:t+L1))+sum(gamma1*H1(1,t+L1-L3:t+L1)));
%         sum(alpha1*q*H1(1,1:t+L1-L1)*(1-exp(-r*(t+L1-L1:1))))+sum(beta1*q*H1(1,1:t+L1-L2)*(1-exp(-r*(t+L1-L2:1))))+sum(gamma1*q*H1(1,1:t+L1-L3)*(1-exp(-r*(t+L1-L3:1)))));
    for j = 1:t+L1-L1
        CP1(1,t) = CP1(1,t) + Kc*alpha1*q*H1(1,j)*(exp(-r*(t-j-L1)));
    end
    for j = 1:t+L1-L2
        CP1(1,t) = CP1(1,t) + Kc*beta1*q*H1(1,j)*(exp(-r*(t-j-L2)));
    end
end

CP2 = zeros(1,200);
H2 = zeros(1,200);
for t = 1:100
    for j = P2:t+T2
        H2(1,t+L1) = H2(1,t+L1)+n2(j-1,t)*(1-m2)*p2*BM2(j);
    end
%     H2(1,t+L1) = n2(P2-1,t)*(1-m2)*p2*BM2(P2);
    CP2(1,t) = Kc*(sum(alpha2*H2(1,t+L1-L1:t+L1))+sum(beta2*H2(1,t+L1-L2:t+L1))+sum(gamma2*H2(1,t+L1-L3:t+L1)));
%         sum(alpha1*q*H1(1,1:t+L1-L1)*(1-exp(-r*(t+L1-L1:1))))+sum(beta1*q*H1(1,1:t+L1-L2)*(1-exp(-r*(t+L1-L2:1))))+sum(gamma1*q*H1(1,1:t+L1-L3)*(1-exp(-r*(t+L1-L3:1)))));
    for j = 1:t+L1-L1
        CP2(1,t) = CP2(1,t) + Kc*alpha2*q*H2(1,j)*(exp(-r*(t-j-L1)));
    end
    for j = 1:t+L1-L2
        CP2(1,t) = CP2(1,t) + Kc*beta2*q*H2(1,j)*(exp(-r*(t-j-L2)));
    end
end

CP3 = zeros(1,200);
H3 = zeros(1,200);
for t = 1:100
    for j = P3:t+T3
        H3(1,t+L1) = H3(1,t+L1)+n3(j-1,t)*(1-m3)*p3*BM3(j);
    end
%     H3(1,t+L1) = n3(P3-1,t)*(1-m3)*p3*BM3(P3);
    CP3(1,t) = Kc*(sum(alpha3*H3(1,t+L1-L1:t+L1))+sum(beta3*H3(1,t+L1-L2:t+L1))+sum(gamma3*H3(1,t+L1-L3:t+L1)));
%         sum(alpha1*q*H1(1,1:t+L1-L1)*(1-exp(-r*(t+L1-L1:1))))+sum(beta1*q*H1(1,1:t+L1-L2)*(1-exp(-r*(t+L1-L2:1))))+sum(gamma1*q*H1(1,1:t+L1-L3)*(1-exp(-r*(t+L1-L3:1)))));
    for j = 1:t+L1-L1
        CP3(1,t) = CP3(1,t) + Kc*alpha3*q*H3(1,j)*(exp(-r*(t-j-L1)));
    end
    for j = 1:t+L1-L2
        CP3(1,t) = CP3(1,t) + Kc*beta3*q*H3(1,j)*(exp(-r*(t-j-L2)));
    end
end

CP = CP1 + CP2 + CP3;

C = CT + CP;

mC = sum(C)/100;
mCT = sum(CT)/100;
mCP = sum(CP)/100;

ddd = mC-C(1,1);

% %% 作图
% plot(1:100,C(1,1:100),1:100,CT(1,1:100),1:100,CP(1,1:100));
% 
% %% 综合评价模型
% x = zeros(1,3);
% %% x1,固碳
% x(1,1) = sum(CT)/100;
% %% x2,经济
% EV1 = (alpha1*Pr1/mu1+beta1*Pr2/mu1+gamma1*Pr3)*sum(H1);
% EV2 = (alpha2*Pr1/mu2+beta2*Pr2/mu2+gamma2*Pr3)*sum(H2);
% EV3 = (alpha3*Pr1/mu3+beta3*Pr2/mu3+gamma3*Pr3)*sum(H3);
% EV = EV1 + EV2 + EV3;
% x(1,2) = EV/100;
% %% x3,稳定程度
% x(1,3) = sum((CT(1,1:100)-x(1,1)).^2)/(100-1);

