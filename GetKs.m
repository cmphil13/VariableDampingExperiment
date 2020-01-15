function [Kp, Kn] = GetKs(xdot, xdotdot)
    
    % va = vel*accel
    va = xdot.*xdotdot;  % the quantity xdot*xdotdot is a scaled version
                         % of the change in kinetic energy

    s = 0.95;           % sensitivity
    va_max = max(va);   % vel-accel max
    va_min = min(va);   % vel-accel min
    
    % Calculate Kp, Kn
    Kp = -log((1-s)/(1+s))/va_max;
    Kn = -log((1+s)/(1-s))/va_min;
end