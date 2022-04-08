f = 850e6;
omega = 2 * pi * f;
Z0 = 50;
Zin = 0.196 - 1i*0.376;
Zout = 9.169 + 1i*8.367;

Qin = sqrt(Z0/abs(Zin) - 1);
Qout = sqrt(Z0/abs(Zout) - 1);

[Cin, Lin] = HP_Lmatch(Z0, real(Zin), f);
Cin_ss = 1 / (omega *  imag(Zin));
[Cout, Lout] = HP_Lmatch(Z0, real(Zout), f);
Cout_ss = 1 / (omega *  imag(Zout));






function [Cs,Lp] = HP_Lmatch(Rhigh, Rlow, fc)
    % Upward HP L-match
    w = 2 * pi * fc;
    Rp = Rhigh;
    Rs = Rlow;
    Q = sqrt(Rp/Rs - 1);
    Lp = Rp / w / Q;
    Cp = 1 / (w*w) / Lp;
    Cs = Cp * (1 + Q*Q) / (Q*Q);  
end