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
fs1 = 1000;
Ts1 = 1/fs1;
fs2 = 100;
Ts2 = 1/fs2;
fs3 = 10;
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

% Linearity check
% u = 10000*ones(10/Ts1,1);
% y1 = exciteSystem(STUDENTID,u,fs1);
% y2 = exciteSystem(STUDENTID,2*u,fs1);

y_data1 = timeseries(y1,linspace(0,(size(y1,1)-1)*Ts1,size(y1,1)));
y_data2 = timeseries(y2,0:Ts2:(size(y2,1)-1)*Ts2);
y_data3 = timeseries(y3,0:Ts3:(size(y3,1)-1)*Ts3);
%y_data2 = timeseries(y2,0:1/(fs/2):99*(2/fs));
figure(2)
plot(y_data1,'.')
hold on
plot(y_data2,'x')
plot(y_data3,'x')
%ylim([-1500 1500])
hold off
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