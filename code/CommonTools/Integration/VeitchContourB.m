function [C,dC]=VeitchContourB(Dp,Wp,Dm,Wm)

% Using the contour given by Veitch and Peake 2008. Intersections the real
% axis at 0. Gives C=C(t)

% On RHS of Im axis uses Dp,Wp, on LHS Dm, Wm.

% D is the height, W the width.

Dp=-Dp;
Dm=Dm; % So both heights above real axis. This is because I define everything backwards, partially.

Bp=3*Wp^4; Ap=4*1i*Dp*Wp^3;
Bm=3*Wm^4; Am=4*1i*Dm*Wm^3;

Cp=@(t) t-Ap*t./(Bp+t.^4);
dCp=@(t) 1-Ap*(Bp-3*t.^4)./(Bp+t.^4).^2;

Cm=@(t) t-Am*t./(Bm+t.^4);
dCm=@(t) 1-Am*(Bm-3*t.^4)./(Bm+t.^4).^2;

C=@(t) Cp(t).*heaviside(t)+Cm(t).*heaviside(-t);
dC=@(t) dCp(t).*heaviside(t)+dCm(t).*heaviside(-t);

end