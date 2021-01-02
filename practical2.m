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



%% DC offset
for i = 1:2000
    y_all(i) = mean(spike_filter(exciteSystem(STUDENTID,0*10000*ones(5/Ts2,1),fs2)));
end
hist(y_all)
% y_data1 = timeseries(y1_filt,0:Ts1:(size(y1_filt,1)-1)*Ts1);
% %plot(y_data1,'.')
% 
% 
% for i=1:size(y_data1.Time,1)
%    if (y_data1.Time(i)>=0.8*time_delay)
%        
%    end
% end

%% Persistency of excitation
omega = 2*pi/4;
sim_time = 40;
t = 0:Ts1:sim_time-Ts1;
ramp = linspace(1,4,size(t,2));
%ramp = [linspace(1,4,size(t,2)/2) linspace(4,1,size(t,2)/2)];
u = sin(2*pi/(2*40)*t).*(10000*sin(omega*(ramp.*t)));
figure()
hold on
plot(0:Ts1:sim_time-Ts1,u,'.')

y = spike_filter(exciteSystem(STUDENTID,u,fs1));
plot(t,y,'.')
hold off

s = 8;
U_0sN = hankel(u(1:s),u(s:end));
rank(U_0sN)

%% Bandwidth check
for i = 1:10:1000
    omega = 2*pi/2;
    sim_time = 40;
    t = 0:Ts2:sim_time-Ts2;
end

%% Assignment 2: Identification
%%% Model Estimation
% We use the following identification method because ...
% As you can see in our figure , we can expect that ...
code = rand ; % important code

%% Assignment 3: Validation
%%% VAF
morecode = 1+1;
%%% Conlusions
% after going through the id - cycle many times we found that by doing this
% and that our results improved ....


%% Functions
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




