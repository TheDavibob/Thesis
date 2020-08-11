function [Lu,Ld] = ImpedanceUpDown(Zu,Zd,omega,U0,sigma0)
%ImpedanceUpDown - computes operator for upstream/downstream impdendance
%solutions

C0 = @(k1) 1i*(omega - U0*k1);
dC0 = @(k1) -1i*sigma0*k1;

if Zu == inf
    Lu.l1 = @(k1) -C0(k1).^2;
    Lu.l2 = @(k1) -3*C0(k1).*dC0(k1);
else
    Lu.l1 = @(k1) 1i*omega*Zu*(-C0(k1).^2) - 0;
    Lu.l2 = @(k1) 1i*omega*Zu*(-3*C0(k1).*dC0(k1)) - (-C0(k1).^4);
end
    
if Zd == inf
    Ld.l1 = @(k1) -C0(k1).^2;
    Ld.l2 = @(k1) -3*C0(k1).*dC0(k1);
else
    Ld.l1 = @(k1) 1i*omega*Zd*(-C0(k1).^2) - 0;
    Ld.l2 = @(k1) 1i*omega*Zd*(-3*C0(k1).*dC0(k1)) - (-C0(k1).^4);
end
end

