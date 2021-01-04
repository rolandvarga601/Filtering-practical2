%% Filtering & Identification Practical Assignment 2
% Jane Doe #123456
% John Doe #654321

% St. id.: 5230543

%% Assignment 1
%%% Input signal
% We chose the following input signal because ...
% The following plot shows the input signal , you can see that ...
% u = [zeros(20,1);10000*ones(80,1)];
u = 10000*ones(10000,1);

figure(1)
plot (u,'.')
xlabel ('Time (s)')
ylabel ('Input ')
title ('Input Signal ')
%%% Sampling Frequency
% We found the sampling frequency by ...
fs1 = 100;
Ts1 = 1/fs1;
fs2 = 50;
Ts2 = 1/fs2;
fs3 = 25;
Ts3 = 1/fs3;
%%% Other important design choises
% We use the comments to provide clear motivations of all our design choices .
% Everything we comment directly under this header will be published as
% text .
% The text under the following line will be published as a green comment .
excitement_level = 100;
% Excitement level

STUDENTID = 5230543;
y1 = exciteSystem(STUDENTID,10000*ones(5/Ts1,1),fs1);
y2 = exciteSystem(STUDENTID,10000*ones(5/Ts2,1),fs2);
y3 = exciteSystem(STUDENTID,10000*ones(5/Ts3,1),fs3);
y_DCoffset =  exciteSystem(STUDENTID,zeros(5/Ts3,1),fs3);
mean(spike_filter(y_DCoffset))

% Linearity check
% u = 10000*ones(10/Ts1,1);
% y1 = exciteSystem(STUDENTID,u,fs1);
% y2 = exciteSystem(STUDENTID,2*u,fs1);

y_data1 = timeseries(y1,linspace(0,(size(y1,1)-1)*Ts1,size(y1,1)));
y_data2 = timeseries(y2,0:Ts2:(size(y2,1)-1)*Ts2);
y_data3 = timeseries(y3,0:Ts3:(size(y3,1)-1)*Ts3);
%y_data2 = timeseries(y2,0:1/(fs/2):99*(2/fs));
figure(2)
plot(y_data1,'.','LineWidth',5)
hold on
plot(y_data2,'.','LineWidth',5)
plot(y_data3,'.','LineWidth',5)
%ylim([-1500 1500])
hold off
xlabel ('Time (s)')
ylabel ('Measured value')
title ('Step response frequency dependence')
legend('fs=100','fs=50','fs=25')

%% Filtering the spikes
y1_filt = spike_filter(y1);

figure
plot(y1_filt,'.');

%% Linearity check
y_lin1 = exciteSystem(STUDENTID,10000*ones(5/Ts2,1),fs2);
y_lin2 = exciteSystem(STUDENTID,2*10000*ones(5/Ts2,1),fs2);
y_lin3 = exciteSystem(STUDENTID,4*10000*ones(5/Ts2,1),fs2);

figure
plot(spike_filter(y_lin1),'.');
hold on
plot(spike_filter(y_lin2),'.');
plot(spike_filter(y_lin3),'.');
hold off
disp(max(spike_filter(y_lin2))/max(spike_filter(y_lin1)))
disp(max(spike_filter(y_lin3))/max(spike_filter(y_lin1)))

%% Delay
% Checked from the plot
delay = 0.9;

%% DC offset
for i = 1:2000
    y_all(i,:) = (spike_filter(exciteSystem(STUDENTID,0*10000*ones(5/Ts2,1),fs2)));
end
figure()
plot(sum(y_all,1)/2000)
mean(sum(y_all,1)/2000)
%hist(y_all)
%mean(y_all)
% y_data1 = timeseries(y1_filt,0:Ts1:(size(y1_filt,1)-1)*Ts1);
% %plot(y_data1,'.')
% 
% 
% for i=1:size(y_data1.Time,1)
%    if (y_data1.Time(i)>=0.8*time_delay)
%        
%    end
% end
% omega = 2*pi;
% sim_time = 20;
% t = 0:Ts2:sim_time-Ts2;
% u = [0*1000000*sin(omega*t(1:size(t,2)*0.8)) zeros(1,size(t,2)*0.2)];
% figure(5)
% hold on
% %plot(0:Ts2:sim_time-Ts2,u,'.')
% y = spike_filter(exciteSystem(STUDENTID,u,fs2));
% plot(t,y,'.')
% mean(y)

%% Persistency of excitation
omega = 2*pi/4;
sim_time = 60;
t = 0:Ts1:sim_time-Ts1;
t_length = size(t,2);
ramp = linspace(1,4,size(t_length*0.6,2));
%ramp = [linspace(1,4,size(t,2)/2) linspace(4,1,size(t,2)/2)];
u = [zeros(1,t_length*0.2)...
    sin(2*pi/(2*40)*t(1:t_length*0.6)).*...
    (100000*sin(omega*(ramp.*t(1:t_length*0.6))))...
    zeros(1,t_length*0.2)];
figure()
hold on
plot(0:Ts1:sim_time-Ts1,u,'.')

y = spike_filter(exciteSystem(STUDENTID,u,fs1));
plot(t,y,'.')
hold off

s = 50;
U_0sN = hankel(u(1:s),u(s:end));
Y_0sN = hankel(y(1:s),y(s:end));
rank(U_0sN)
figure()
semilogy(svd(Y_0sN),'x')

%% Bandwidth check
max_values = zeros(size(10:10:270,2),1);
min_values = zeros(size(10:10:270,2),1);
for i = 10:10:270
    %omega = 2*pi/2*i/1000;
    omega = 2/2*i/1000;
    h = 0.4/omega;
    sim_time = 40;
    t = 0:Ts2:sim_time-Ts2;
    u = 10000*sin(omega*t);
    y = spike_filter(exciteSystem(STUDENTID,u,1/h));
    max_values(i) = max(y);
    min_values(i) = min(y);
end

%%
for i = 500
    omega = 2*pi/2*i/100;
    h = 0.4/omega;
    sim_time = 40;
    t = 0:Ts2:sim_time-Ts2;
    u = 10000*sin(omega*t);
    y = spike_filter(exciteSystem(STUDENTID,u,1/h));
    max_values(i) = max(y);
    min_values(i) = min(y);
end

%% Assignment 2: Identification
%%% Model Estimation
% We use the following identification method because ...
% As you can see in our figure , we can expect that ...

delay = 0.9;
delay_samples = 0.9/Ts1;
y_nodelay = y(delay_samples:end);
u_nodelay = u(1:end-delay_samples+1)';

r = 0.3;
lim = floor(size(u_nodelay,1)*r);

method = 'po-moesp';
n = 2;
s = 50;
[A,B,C,D,x0,sv] = subspaceID(u_nodelay(1:lim),y_nodelay(1:lim),s,n,method);
[yhat,xhat] = simsystem(A,B,C,D,x0,u_nodelay);

%% Results
figure
hold on
plot(y_nodelay(lim+1:end),'.')
plot(yhat(lim+1:end),'.')
vaf(y_nodelay(lim+1:end),yhat(lim+1:end))


%% Assignment 3: Validation
%%% VAF
morecode = 1+1;
%%% Conlusions
% after going through the id - cycle many times we found that by doing this
% and that our results improved ....
scaler = 10;
uv=zeros(400*scaler,1);
uv(20*scaler:30*scaler,1)=100000;
uv(50*scaler:70*scaler,1)=200000;
uv(70*scaler:90*scaler,1)=-200000;
uv(130*scaler:140*scaler,1)=300000;
uv(180*scaler:190*scaler,1)=-100000;

yv = spike_filter(exciteSystem(STUDENTID,uv,fs1));
y_uf = exciteSystem(STUDENTID,uv,fs1);
yv_uf_nodelay = y_uf(delay_samples:end);

yv_nodelay = yv(delay_samples:end);
uv_nodelay = uv(1:end-delay_samples+1);

method = 'po-moesp';
n = 2;
s = 50;
[Av,Bv,Cv,Dv,x0v,svv] = subspaceID(uv_nodelay,yv_nodelay,s,n,method);
[yhatv,xhatv] = simsystem(Av,Bv,Cv,Dv,x0v,uv_nodelay);

figure
hold on
plot(yv_nodelay,'.')
%plot(yv_uf_nodelay(),'.')
plot(yhatv,'.')

%%
vaf(yv_nodelay,yhatv)

%% Functions
function [A,B,C,D,x0,sv] = subspaceID(u,y,s,n,method)
% Function INPUT 
% u         system input (matrix of size N x m)
% y         system output (matrix of size N x l)
% s         block size (scalar)
% n         model order (scalar)
% method    method (string e.g. 'moesp')
%
% Function OUTPUT
% A         System matrix A (matrix of size n x n)
% B         System matrix B (matrix of size n x m)
% C         System matrix C (matrix of size l x n)
% D         System matrix D (matrix of size l x m)
% x0        Initial state (vector of size n x one)
% sv        Singular values (vector of size n x one)
    
    switch method
        case 'moesp'
            % Computation of the Hankel matrices
            UY = [hankel(u(1:s),u(s:end)); hankel(y(1:s),y(s:end))];
            
            % LQ factorization to get L22
            [~,R] = qr(UY',0);
            L = R(1:2*s,1:2*s)';
            L22 = L(size(L,1)/2+1:end,size(L,2)/2+1:end);
            
            % SVD of L22
            [U,S,~] = svd(L22,0);
            
        case 'pi-moesp'
            % Computation of the Hankel matrices
            U = hankel(u(1:2*s),u(2*s:end));
            Up = U(1:s,:);
            Uf = U(s+1:2*s,:);
            Y = hankel(y(1:2*s),y(2*s:end));
            Yp = Y(1:s,:);
            Yf = Y(s+1:2*s,:);
            
            % LQ factorization to get L32
            [~,R] = qr([Uf' Up' Yf'],0);
            L = R(1:3*s,1:3*s)';
            L32 = L(2*s+1:end,s+1:2*s);
            
            % SVD of L32
            [U,S,~] = svd(L32,0); 
            
        case 'po-moesp'
            % Computation of the Hankel matrices
            U = hankel(u(1:2*s),u(2*s:end));
            Up = U(1:s,:);
            Uf = U(s+1:2*s,:);
            Y = hankel(y(1:2*s),y(2*s:end));
            Yp = Y(1:s,:);
            Yf = Y(s+1:2*s,:);
            
            % Constructing the instrumental variable
            Z = [Up;Yp];
            
            % LQ factorization to get L32
            [~,R] = qr([Uf; Z; Yf]',0);
            L = R(1:4*s,1:4*s)';
            L32 = L(3*s+1:end,s+1:3*s);
            
            % SVD of L32
            [U,S,~] = svd(L32,0);
    end
    
    % Storing singular values in output 'sv'
    sv = diag(S);
    %semilogy(sv,'x')
            
    % Computing A by linear least squares
    A = U(1:s-1,1:n)\U(2:s,1:n);

    % Extracting C
    C = U(1,1:n);
    
    % Construction of the phi matrix transpose
    m = size(u,2);
    l = size(y,2);
    phi_t = zeros(size(y,1),n+n*m+l*m);
    for k=0:size(y,1)-1
        phi_t1 = C*(A^k);
        phi_t2 = zeros(1,n*m);
        for j=0:k-1
            phi_t2 = phi_t2 + kron(u(j+1),C*(A^(k-1-j)));
        end
        phi_t3 = u(k+1);
        phi_t(k+1,:) = [phi_t1 phi_t2 phi_t3];
    end

    % Linear least squares to get the solution vector
    sol = phi_t\y;
    
    % Extracting information from the vector
    x0 = sol(1:n);
    B = sol(n+1:end-1);
    D = sol(end);
end


function y_filt = spike_filter(y)
    std_dev = std(y);
    m = mean(y);

    y_filt = y;

    ws_start = 0;
    ws_end = 0;

    for i=1:size(y,1)
       if (y(i)>m+std_dev*4)
           if (ws_start == 0)
               ws_start = i-1;
           else
               % Nothing to do
           end
       else
           if (ws_start ~= 0)
               ws_end = i;

               % Doing the interpolation
               x = [ws_start; ws_end];
               v = [y(ws_start); y(ws_end)];
               xq = ws_start:ws_end;
               vq1 = interp1(x,v,xq);

               % Substituting with the correct values
               y_filt(xq) = vq1;

               ws_start = 0;
               ws_end = 0;
           end
       end 
    end
    
    if (ws_start ~= 0)
        y_filt(ws_start:end) = y_filt(ws_start);
    end
    
end

function v = vaf(y, y_est)
    % Variance Accounted For (VAF) | Percentage value (%)
    % y     : measured output
    % y_est : estimated output
    
    v = var(y - y_est) / var(y) ;
    v = 100 * ( 1 - v );
    
    if ( v < 0 )
        v = 0;
    end
end



