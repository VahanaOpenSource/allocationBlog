% allocationBlog.m
% by Geoffrey Bower
% Dec 2018
%
% Example implementation of a multirotor allocation using MPT3 toolbox.
% Includes example for octocopter with equally spaced rotors both in
% nominal conditions and with a single failed rotor.
%
% Copyright (c) 2018 A³ by Airbus
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

clear all
clc

%% Check for mpt3 installation
if exist('mptdoc','file') ~= 2
    error('Follow installation instructions for MPT3: https://www.mpt3.org/Main/Installation');
end

%% Vehicle Parameters
m = 5;    % Vehicle mass [kg]
g = 9.8;  % Gravitational acceleration [m/s^2]
r = 0.35; % Rotor radius [m]
rotDir = [+1,-1,+1,-1,+1,-1,+1,-1]; % +ve is right handed rotation about thrust direction
n = length(rotDir); % Number of rotors

% Rotor locations in x forward, y right, z down coordinate frame with origin at CG
% Assumed to be axisymmetric about CG
phi = linspace(0,360,n+1)'; phi = phi(1:end-1); % Azimuth
xRotor = sind(phi);  % Rotor x location [m]
yRotor = cosd(phi);  % Rotor y location [m]
zRotor = zeros(n,1); % Rotor z locaton [m]

% Rotor parameters
Ct = 0.014; % Rotor thrust coefficient [T = rho * pi * r^4 * Ct * omega^2]
FOM = 0.7;  % Rotor figure of merit
Cq = Ct^1.5 / FOM / sqrt(2); % Rotor torque coefficient [Q = rho * pi * r^5 * Cq * omega^2]

% Motor sizing
thrustRatio = 1.4; % Maximum individual rotor thrust over average thrust in hover

%% Generate hover control Jacobian: dFM/dTau
% Ignore Fx and Fy (assume no cant on rotors)
dFMdTau = zeros(4,n);

% Thrust per torque
dTdTau = Ct / (r * Cq);

% Fz = -dTdTau
dFMdTau(1,:) = -dTdTau;

% Mx = dTdTau * -y
dFMdTau(2,:) = dTdTau * -yRotor;

% My = dTdTau * x
dFMdTau(3,:) = dTdTau * xRotor;

% Mz = rotDir
dFMdTau(4,:) = rotDir;

%% Allocation parameters 
% Force and moment limtis
rollPitchMomentLimits = [-15,15]; % Roll and pitch moment limits [Nm]
yawMomentLimits = [-3,3]; % Yaw moment limits [Nm]
thrustLimits = [-1.2 * m * g, -0.8 * m * g]; % Thrust limits [Nm]

% Force and moment weighting
WFM = diag([20,100,100,5]); % Fz, Mx, My, Mz

%% Nominal case with no failures

% Force/Moment trim command, assuming CG is at origin 
FMTrim = [-m * g; 0; 0; 0]; % 1g thrust command, 0 moment commands

% Compute trim actuation using Moore-Penrose pseudo-inverse
xTrim = pinv(dFMdTau) * FMTrim;

% Create optimization variables and parameters
x = sdpvar(n,1); % design variables (actuator commands)
xMin = zeros(size(xTrim));
xMax = thrustRatio * mean(xTrim) * ones(size(xTrim));
FMCmd = sdpvar(4,1);
FMCmdMin = [thrustLimits(1); rollPitchMomentLimits(1); rollPitchMomentLimits(1); yawMomentLimits(1)];
FMCmdMax = [thrustLimits(2); rollPitchMomentLimits(2); rollPitchMomentLimits(2); yawMomentLimits(2)];

% Objective function:
% First term: weighted square of errors in achieving Force/Moment commands
% Second term: Minimize deviation from trim actuation in the null-space
J = (dFMdTau * x - FMCmd)' * WFM * (dFMdTau * x - FMCmd) + sum((null(dFMdTau)'*(x-xTrim)).^2);

% Constraints
% Actuators within bounds, force moment commands within limits
C = [ x >= xMin, x <= xMax, ...
    FMCmd >= FMCmdMin, FMCmd <= FMCmdMax];

% Convert the problems into the MPT format:
p = Opt(C, J, FMCmd, x);

% Solve the parameteric problem
s = p.solve();


%% Motor 2 failure example

% Set maximum motor 2 torque to a small number
xMax(2) = 1e-4; 
xTrim2 = pinv(dFMdTau(:,setdiff(1:n,2))) * FMTrim;

% Adjust actuator Jacobian for motor 2 failure
dFMdTau2 = dFMdTau;
dFMdTau2(:,2) = [];

% Objective
% First term: weighted square of errors in achieving Force/Moment commands
% Second term: Minimize deviation from trim actuation in the null-space of
% remaining actuators
J2 = (dFMdTau * x - FMCmd)' * WFM * (dFMdTau * x - FMCmd) + sum((null(dFMdTau2)'*(x(setdiff(1:n,2))-xTrim2)).^2);

% Actuator constraints
C2 = [ x >= xMin, x <= xMax, ...
    FMCmd >= FMCmdMin, FMCmd <= FMCmdMax];

% Convert the problems into the MPT format:
p2 = Opt(C2, J2, FMCmd, x);

% Solve the parameteric problems
s2 = p2.solve();


%% Evaluate solution space

% Mx and My commands to loop over
MxCmd = rollPitchMomentLimits(1):rollPitchMomentLimits(2);
MyCmd = rollPitchMomentLimits(1):rollPitchMomentLimits(2);
FMCmdXY = zeros(length(MxCmd),length(MyCmd),4);

% Storage for solutions
xXY = zeros(length(MxCmd),length(MyCmd),n);
xXY2 = zeros(length(MxCmd),length(MyCmd),n);

% Storage for objective value
fvalXY = zeros(length(MxCmd),length(MyCmd));
fvalXY2 = zeros(length(MxCmd),length(MyCmd));

% Storage for realized force/moments
FMRetXY = zeros(length(MxCmd),length(MyCmd),4);
FMRetXY2 = zeros(length(MxCmd),length(MyCmd),4);

for i = 1:length(MxCmd)
    for j = 1:length(MyCmd)
        
        % Commanded force/moment
        FMTest = FMTrim; FMTest(2) = MxCmd(i); FMTest(3) = MyCmd(j);
        FMCmdXY(i,j,:) = FMTest;
        
        % Find optimal design variables
        xTmp = s.xopt.feval(FMTest, 'primal');
        xTmp2 = s2.xopt.feval(FMTest, 'primal');
        xXY(i,j,:) = xTmp(:,1);
        xXY2(i,j,:) = xTmp2(:,1);
        
        % Find objective value
        tmp = s.xopt.feval(FMTest, 'obj');
        tmp2 = s2.xopt.feval(FMTest, 'obj');
        fvalXY(i,j) = tmp(1);
        fvalXY2(i,j) = tmp2(1);

        % Compute returned force/moments for optimal design variables
        FMRetXY(i,j,:) = dFMdTau * squeeze(xXY(i,j,:));
        FMRetXY2(i,j,:) = dFMdTau * squeeze(xXY2(i,j,:));
    end
end


%% Plots

% Plot rotor geometry
figure(1); clf; hold on;
plot(yRotor,xRotor,'k.','MarkerSize',20);
theta = linspace(0,2*pi);
for i = 1:n
    if rotDir(i) > 0
        plot(yRotor(i) + r * sin(theta), xRotor(i) + r * cos(theta),'b')
    else
        plot(yRotor(i) + r * sin(theta), xRotor(i) + r * cos(theta),'r')
    end
    text(yRotor(i),xRotor(i),['   ',num2str(i)])
end
grid on
axis equal
ylabel('x [m]')
xlabel('y [m]')
title('Rotor Geometry')

% Plot lines showing commanded and realized Mx and My with Fz = -m*g and Mz = 0
figure(2); clf; hold on;
for i = 1:length(MxCmd)
    for j = 1:length(MyCmd)
        plot([FMCmdXY(i,j,2),FMRetXY(i,j,2)],[FMCmdXY(i,j,3),FMRetXY(i,j,3)],'b')
        plot(FMRetXY(i,j,2),FMRetXY(i,j,3),'r.')
    end
end
grid on
xlabel('Mx [N-m]')
ylabel('My [N-m]')
axis equal
title('FM Commanded to FM Realized as Mx and My Cmd vary with Fz = -m*g, Mz = 0')

% Surface plots of errors in Fz, Mx, My, Mz as a function of Mx and My
% commands with Fz = -m*g and Mz = 0
FMName = {'Fz','Mx','My','Mz'};
figure(3); clf; hold on;
for i = 1:4
    subplot(2,2,i); hold on;
    surface(MxCmd, MyCmd, squeeze(FMCmdXY(:,:,i) - FMRetXY(:,:,i))','EdgeAlpha',0,'FaceColor','interp');
    grid on;
    xlabel('Mx')
    ylabel('My')
    title([FMName{i},' Errors'])
    colorbar
    if i == 1
        caxis([-20,20])
        title('Fz Errors')
    elseif i == 2 
        caxis([-10,10])
        title('Mx Errors')
    elseif i == 3
        caxis([-10,10])
        title('My Errors')
    elseif i == 4
        caxis([-1,1])
        title('Mz Errors')
    end
end

% Plot projection of polytopes for nominal case into each 3D FM set.
figure(4); clf;
FM = {'Fz','Mx','My','Mz'};
for i = 1:4
    subplot(2,2,i);
    ind = setdiff(1:4,i);
    tmp = projection(s.xopt.Set,ind,'vrep');
    tmp.plot('Alpha',1);
    xlabel(FM{ind(1)})
    ylabel(FM{ind(2)})
    zlabel(FM{ind(3)})
end
title('Polytope projections for nominal case')

% Plot lines showing commanded and realized Mx and My with rotor 2 failed, Fz = -m*g and Mz = 0
figure(5); clf; hold on;
for i = 1:length(MxCmd)
    for j = 1:length(MyCmd)
        plot([FMCmdXY(i,j,2),FMRetXY2(i,j,2)],[FMCmdXY(i,j,3),FMRetXY2(i,j,3)],'b')
        plot(FMRetXY2(i,j,2),FMRetXY2(i,j,3),'r.')
    end
end
grid on
xlabel('Mx [N-m]')
ylabel('My [N-m]')
axis equal
title('FM Commanded to FM Realized as Mx and My Cmd vary with motor 2 failure, Fz = -m*g, and Mz = 0')

% Surface plots of errors in Fz, Mx, My, Mz as a function of Mx and My
% commands with rotor 2 failed, Fz = -m*g and Mz = 0
FMName = {'Fz','Mx','My','Mz'};
figure(6); clf; hold on;
for i = 1:4
    subplot(2,2,i); hold on;
    surface(MxCmd, MyCmd, squeeze(FMCmdXY(:,:,i) - FMRetXY2(:,:,i))','EdgeAlpha',0,'FaceColor','interp');
    grid on;
    xlabel('Mx')
    ylabel('My')
    title([FMName{i},' Errors'])
    colorbar
    if i == 1
        caxis([-20,20])
        title('Fz Errors with Motor 2 failure')
    elseif i == 2
        caxis([-10,10])
        title('Mx Errors with Motor 2 failure')
    elseif i == 3
        caxis([-10,10])
        title('My Errors with Motor 2 failure')
    elseif i == 4
        caxis([-1,1])
        title('Mz Errors with Motor 2 failure')
    end
end


% Plot projection of polytopes for rotor 2 failure case into each 3D FM set.
figure(7); clf;
FM = {'Fz','Mx','My','Mz'};
for i = 1:4
    subplot(2,2,i);
    ind = setdiff(1:4,i);
    tmp = projection(s2.xopt.Set,ind,'vrep');
    tmp.plot('Alpha',1);
    xlabel(FM{ind(1)})
    ylabel(FM{ind(2)})
    zlabel(FM{ind(3)})
end
title('Polytope projection for rotor 2 out');




