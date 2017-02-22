%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  16833 Robot Localization and Mapping  % 
%  Assignment #1                         %
%  EFK-SLAM                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

%==== TEST: Setup uncertianty parameters (try different values!) ===
sig_x = 0.25;
sig_y = 0.1;
sig_alpha = 0.1;
sig_beta = 0.01;
sig_r = 0.08;

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
        landmark(i)= cos(measure(i))*measure(i+1);
        landmark(i+1)= sin(measure(i))*measure(i+1);
end

%Measurement Jacobian w.r.t pose
dz_dp = @(i,j) [landmark(j)/(landmark(i)^2+landmark(j)^2), -landmark(i)/(landmark(i)^2+landmark(j)^2), -1;...
               -landmark(i)/sqrt((landmark(i)^2+landmark(j)^2)) , -landmark(j)/sqrt((landmark(i)^2+landmark(j)^2)), 0];
%Measurement Jacobian w.r.t landmark
dz_dl = @(i,j) [-landmark(j)/(landmark(i)^2+landmark(j)^2), landmark(i)/(landmark(i)^2+landmark(j)^2); ...
            landmark(i)/sqrt((landmark(i)^2+landmark(j)^2)), landmark(j)/sqrt((landmark(i)^2+landmark(j)^2))]; 
landmark_cov = zeros(k*2, k*2);
for i = 1 : 2 : k*2
   j = i+1;
   landmark_cov(i:j,i:j) = dz_dp(i,j)*pose_cov*dz_dp(i,j)'+dz_dl(i,j)*measure_cov*dz_dl(i,j)';
   %landmark_cov(i:j,i:j) = dz_dl(i,j)*[landmark(i);landmark(i+1)]*dz_dl(i,j)';
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
    theta_t = x(3);
    x_pre = x + [cos(theta_t)*d; sin(theta_t*d); alpha; zeros(k*2,1)];
    
    % State Jacobian w.r.t pose
    G_robot = [1, 0, -d*sin(theta_t);...
               0, 1,  d*cos(theta_t);...
               0, 0,  1              ];
    G = zeros(3+2*k, 3+2*k);
    G(1:3,1:3) = G_robot;
    % State Jacobian w.r.t input
    B_robot = [cos(theta_t), 0 ,0;...
               sin(theta_t), 0, 0;...
                          0, 0, 1];
    B = zeros(3+2*k, 3+2*k);
    B(1:3,1:3) = B_robot;   

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
    Q = zeros(3+2*k, 3+2*k);
    for i = 1 : 2 : 2*k
        Q(i+3:i+4,i+3:i+4) = measure_cov;
    end
    
    z = zeros(k*2,1);
    for i = 1 : 2 : k*2
            z(i)= cos(measure(i))*measure(i+1);
            z(i+1)= sin(measure(i))*measure(i+1);
    end
    %Measurement Jacobian w.r.t pose
    delta_x = @(i) z(i)-x_pre(1);
    delta_y = @(i) z(i)-x_pre(2);
    dz_dp = @(i,j) [delta_y(j)/(delta_x(i)^2+delta_y(j)^2), -delta_x(i)/(delta_x(i)^2+delta_y(j)^2), -1;...
                   -delta_x(i)/sqrt((delta_x(i)^2+delta_y(j)^2)) , -delta_y(j)/sqrt((delta_x(i)^2+delta_y(j)^2)), 0];
    %Measurement Jacobian w.r.t landmark
    dz_dl = @(i,j) [-delta_y(j)/(delta_x(i)^2+delta_y(j)^2), delta_x(i)/(delta_x(i)^2+delta_y(j)^2); ...
                delta_x(i)/sqrt((delta_x(i)^2+delta_y(j)^2)), delta_y(j)/sqrt((delta_x(i)^2+delta_y(j)^2))]; 

    H = zeros(3+2*k, 3+2*k);
    H(1:3,1:3) = eye(3,3);
    for i = 1 : 2 : 2*k
        H(i+3:i+4,1:3) = dz_dp(i,i+1);
        H(i+3:i+4,i+3:i+4) = dz_dl(i,i+1);
    end
    
    z_hat = zeros(2*k,1);
    for i = 1:2:2*k
        %z_hat(i) = wrapToPi(atan2(x_pre(i+4)-x_pre(2),x_pre(i+3)-x_pre(1))-x_pre(3));
        z_hat(i) = (atan2(x_pre(i+4)-x_pre(2),x_pre(i+3)-x_pre(1))-x_pre(3));
        z_hat(i+1) = sqrt((x_pre(i+3)-x_pre(1))^2 + (x_pre(i+4)-x_pre(2))^2);
    end
    
    K = P_pre*H'*(H*P_pre*H'+Q)^-1;
    
    z = [0;0;0;z];
    z_hat = [0;0;0;z_hat];
    x = x_pre + K*(z - z_hat);
    P = (eye(2*k+3,2*k+3)-K*H)*P_pre;
    %==== Plot ====   
    drawTrajAndMap(x, last_x, P, t);
    last_x = x;
    
    %==== Iteration & read next control data ===
    t = t + 1;
    tline = fgets(fid);
end

%==== EVAL: Plot ground truth landmarks ====

% Write your code here...
    

%==== Close data file ====
fclose(fid);
