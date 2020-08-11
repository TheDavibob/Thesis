% Some potential contours, nicely setup

% Straight lines
% contour.dCp = @(t) 1;
% contour.Cp = @(t) t - 0.5*1i;
% contour.Cm = @(t) t + 0.5*1i;
% contour.dCm = @(t) 1;

% TanhContourB
% [contour.Cp,contour.dCp] = TanhContour(-0.5,1,0,1);
% [contour.Cm,contour.dCm] = TanhContour(0.5,1,-1,0);

% Pair of separated Tanh contours
% [contour.Cp,contour.dCp] = TanhContour(-0.5,1,-1,1);
% [contour.Cm,contour.dCm] = TanhContour(0.5,1,-1,1);

% Pair of identical Tanh contours
[contour.Cp,contour.dCp] = TanhContour(0,1,-1,1);
[contour.Cm,contour.dCm] = TanhContour(0,1,-1,1);