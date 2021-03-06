%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  16833 Robot Localization and Mapping  % 
%  Assignment #1                         %
%  EFK-SLAM                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clf;

%==== TEST: Setup uncertianty parameters (try different values!) ===
sig_x = 0.25;%0.25;
sig_y = 0.1;%0.1;
sig_alpha = 0.1;%0.1;
sig_beta = 0.01;%0.01;
sig_r = 0.08;%0.08;

%==== Generate sigma^2 from sigma ===
sig_x2 = sig_x^2;
sig_y2 = sig_y^2;
sig_alpha2 = sig_alpha^2;
sig_beta2 = sig_beta^2;
sig_r2 = sig_r^2;

%==== Open data file ====
fid = fopen('../data/data.txt');

%==== Read first measurement data ====
tline = fgets(fid);
arr = str2num(tline);
measure = arr';
t = 1;
 
%==== Setup control and measurement covariances ===
control_cov = diag([sig_x2, sig_y2, sig_alpha2]);
measure_cov = diag([sig_beta2, sig_r2]);

%==== Setup initial pose vector and pose uncertainty ====
pose = [0 ; 0 ; 0];
pose_cov = diag([0.02^2, 0.02^2, 0.1^2]);

%==== TODO: Setup initial landmark vector landmark[] and covariance matrix landmark_cov[] ====
%==== (Hint: use initial pose with uncertainty and first measurement) ====

% Write your code here...
k = length(measure)/2;
landmark = zeros(k*2,1);
for i = 1 : 2 : k*2
        landmark(i)= pose(1)+cos(measure(i)+pose(3))*measure(i+1);
        landmark(i+1)= pose(2) + sin(measure(i)+pose(3))*measure(i+1);
end

%Measurement Jacobian w.r.t pose
dzdp = @(i,j) [1, 0, -measure(j)*sin(measure(i));...
               0, 1, measure(j)*cos(measure(i))];
%Measurement Jacobian w.r.t landmark
dzdl = @(i,j) [-measure(j)*sin(measure(i)), cos(measure(i)); ...
                measure(j)*cos(measure(i)), sin(measure(i))]; 
landmark_cov = zeros(k*2, k*2);
for i = 1 : 2 : k*2
   j = i+1;
   landmark_cov(i:j,i:j) = dzdp(i,j)*pose_cov*dzdp(i,j)'+dzdl(i,j)*measure_cov*dzdl(i,j)';
end
%==== Setup state vector x with pose and landmark vector ====
x = [pose ; landmark];

%==== Setup covariance matrix P with pose and landmark covariances ====
P = [pose_cov zeros(3, 2*k) ; zeros(2*k, 3) landmark_cov];

%==== Plot initial state and conariance ====
last_x = x;
drawTrajAndMap(x, last_x, P, 0);

%==== Read control data ====
tline = fgets(fid);
while ischar(tline)
    arr = str2num(tline);
    d = arr(1);
    alpha = arr(2);
    
    %==== TODO: Predict Step ====
    %==== (Notice: predict state x_pre[] and covariance P_pre[] using input control data and control_cov[]) ====
    
    % Write your code here...
    
    % Predict state space
    theta_t = x(3);
    x_pre = x + [cos(theta_t)*d; sin(theta_t)*d; alpha; zeros(2*k,1)];

    % Predict Covariance
    %    State Jacobian w.r.t pose
    G = eye(3+2*k, 3+2*k);
    G(1:3,1:3) = [1, 0, -d*sin(theta_t);...
                  0, 1,  d*cos(theta_t);...
                  0, 0,  1];
    %    State Jacobian w.r.t input
    B = zeros(3+2*k, 3+2*k);
    B(1:3,1:3) = [cos(theta_t), -sin(theta_t) ,0;...
                  sin(theta_t), cos(theta_t), 0;...
                             0, 0, 1];

    C_cov = zeros(3+2*k, 3+2*k);
    C_cov(1:3,1:3) = control_cov;
    P_pre = G*P*G' + B*C_cov*B';
    
    %==== Draw predicted state x_pre[] and covariance P_pre[] ====
    drawTrajPre(x_pre, P_pre);
    
    %==== Read measurement data ====
    tline = fgets(fid);
    arr = str2num(tline);
    measure = arr';
    
    %==== TODO: Update Step ====
    %==== (Notice: update state x[] and covariance P[] using input measurement data and measure_cov[]) ====
    
    % Write your code here...
  
    B = zeros(k*2, k*2);
    for i = 1 : 2 : k*2
       B(i:i+1,i:i+1) = measure_cov;
    end
    
    % beta, alpha to lx, ly
    z = zeros(k*2,1);
    for i = 1 : 2 : k*2
            z(i)= x_pre(1)+cos(measure(i)+x_pre(3))*measure(i+1);
            z(i+1)= x_pre(2)+sin(measure(i)+x_pre(3))*measure(i+1);
    end
    
    %Landmark Jacobian w.r.t pose
    delta_x = @(i) x_pre(i+3)-x_pre(1);
    delta_y = @(j) x_pre(j+3)-x_pre(2);
    dldp = @(i,j) [delta_y(j)/(delta_x(i)^2+delta_y(j)^2), -delta_x(i)/(delta_x(i)^2+delta_y(j)^2), -1;...
                   -delta_x(i)/sqrt((delta_x(i)^2+delta_y(j)^2)) , -delta_y(j)/sqrt((delta_x(i)^2+delta_y(j)^2)), 0];
    %Landmark Jacobian w.r.t landmark
    dldl = @(i,j) [-delta_y(j)/(delta_x(i)^2+delta_y(j)^2), delta_x(i)/(delta_x(i)^2+delta_y(j)^2); ...
                delta_x(i)/sqrt((delta_x(i)^2+delta_y(j)^2)), delta_y(j)/sqrt((delta_x(i)^2+delta_y(j)^2))]; 

    H = zeros(2*k, 2*k+3);
    for i = 1 : 2 : 2*k
        H(i:i+1,1:3) = dldp(i,i+1);
        H(i:i+1,i+3:i+4) = dldl(i,i+1);
    end
    
    z_hat = zeros(2*k,1);
    for i = 1:2:2*k
        z_hat(i) = wrapToPi(atan2(x_pre(i+4)-x_pre(2),x_pre(i+3)-x_pre(1))-x_pre(3));
        z_hat(i+1) = sqrt((x_pre(i+3)-x_pre(1))^2 + (x_pre(i+4)-x_pre(2))^2);
    end

    K = P_pre*H'*(H*P_pre*H'+B)^-1;
    x = x_pre + K*(measure - z_hat);
    P = (eye(2*k+3,2*k+3)-K*H)*P_pre;
    %==== Plot ====   
    drawTrajAndMap(x, last_x, P, t);
    last_x = x;
    drawnow limitrate;
    %==== Iteration & read next control data ===
    t = t + 1;
    tline = fgets(fid);
    
end
%==== EVAL: Plot ground truth landmarks ====
landmark_true_x = [6,6,10,10,14,14];
landmark_true_y = [6,12,6,12,6,12];
% Write your code here...
plot(landmark_true_x, landmark_true_y, '.k','MarkerSize',20);

EUdistance = zeros(6,1);
MANdistance = zeros(6,1);
for i = 1:2:2*k
    errDis = [x(i+3)-landmark_true_x(floor(i/2+1));
              x(i+4)-landmark_true_y(floor(i/2+1))];
    EUdistance(floor(i/2+1)) = sqrt(dot(errDis,errDis));
    MANdistance(floor(i/2+1)) = sqrt(errDis'*P(i+3:i+4,i+3:i+4)^-1*errDis);

end

%==== Close data file ====
fclose(fid);
