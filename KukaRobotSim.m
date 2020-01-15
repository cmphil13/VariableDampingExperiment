% Calculate end effector x position as a function of joint angles
% Make KUKA robot with Corke Robotic Toolbox 
alpha = [-pi/2, pi/2, -pi/2, -pi/2, pi/2, pi/2, -pi/2, 0];
a = [0, 0, 0, 0, 0, 0, 0, 0];
d = [1, 0.36, 0, 0.42, 0, 0.4, 0, 0.126];
theta=[-pi/2, 0 0 0 0 0 0 -0.958709];
dh = [theta' d' a' alpha'];

% BUILD ROBOT--------------------------------------------------------------
for i = 1:length(dh(:,1))
    L{i} = Link('d', dh(i,2), 'a', dh(i,3), 'alpha', dh(i,4));
end
robot = SerialLink([L{1} L{2} L{3} L{4} L{5} L{6} L{7} L{8}]);


q0 = [pi/2, -pi/2, pi/2, 0, pi/2, 0, -pi/2, 0];

orange = [ 0.9100 0.4100 0.1700];

fig = figure;
robot.plot(q0, 'floorlevel', 0, 'workspace', [-1.5, 0, -1.5, 1.5, 0, 1.5])
hold on
cyl('z', 0.15, 'c', [0, 0.3], [-1.125, 0, 0])
cyl('z', 0.15, orange, [0, 0.3], [-0.76, 1.125-0.76, 0])


function cyl(ax, r, color, extent, offset, varargin)
    
    if isempty(offset)
        offset = [0 0 0];
    end
    
    %fprintf('   cyl: %s, r=%f, extent=[%g, %g]\n', ax, r, extent);
    
    n = 20;
    
    r = [r;r];
    
    theta = (0:n)/n*2*pi;
    sintheta = sin(theta); sintheta(n+1) = 0;
    
    switch ax
        case 'x'
            y = r * cos(theta) + offset(2);
            z = r * sintheta + offset(3);
            x = extent(:) * ones(1,n+1) + offset(1);
        case 'y'
            x = r * cos(theta) + offset(1);
            z = r * sintheta + offset(3);
            y = extent(:) * ones(1,n+1) + offset(2);
        case 'z'
            x = r * cos(theta) + offset(1);
            y = r * sintheta + offset(2);
            z = extent(:) * ones(1,n+1) + offset(3);
    end
    
    % walls of the cylinder
    surf(x,y,z, 'FaceColor', color, 'EdgeColor', 'none', varargin{:})
    
    % put the ends on
    patch(x', y', z', color, 'EdgeColor', 'none', varargin{:});
end