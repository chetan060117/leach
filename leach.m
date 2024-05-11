% LEACH Protocol

% Clearing workspace and command window
clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Field Dimensions - Maximum x and y dimensions of the network area (in meters)
xm = 100; % x dimension
ym = 100; % y dimension

% Coordinates of the Sink - Position of the base station
sink.x = 0.5 * xm;
sink.y = 1.75 * ym;

% Number of Nodes in the field
n = 100;

% Optimal Election Probability of a node to become cluster head
p = 0.05;

% Energy Model (all values in Joules)
Eo = 0.5; % Initial Energy of nodes
ETX = 50 * 0.000000001; % Energy consumed during transmission
ERX = 50 * 0.000000001; % Energy consumed during reception
Efs = 10 * 0.000000000001; % Energy for free space model
Emp = 0.0013 * 0.000000000001; % Energy for multi-path fading model
EDA = 5 * 0.000000001; % Energy consumed during data aggregation

% Values for Heterogeneity 
m = 0.05; % Percentage of nodes that are advanced

% Parameter alpha for energy model
a = 0.1;

% Maximum number of rounds for the simulation
rmax = 2000;
%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%

% Computation of do - Calculating the optimal distance for energy calculations
do = sqrt(Efs / Emp);

% Creation of the random Sensor Network
figure(1);
for i = 1:1:n
    S(i).xd = rand(1, 1) * xm;
    XR(i) = S(i).xd;
    S(i).yd = rand(1, 1) * ym;
    YR(i) = S(i).yd;
    S(i).G = 0;
    S(i).type = 'N'; % initially there are no cluster heads only nodes
    
    temp_rnd0 = i;
    
    % Random Election of Normal Nodes
    if (temp_rnd0 >= m * n + 1) 
        S(i).E = Eo;
        S(i).ENERGY = 0;
        plot(S(i).xd, S(i).yd, 'o');
        hold on;
    end
    
    % Random Election of Advanced Nodes
    if (temp_rnd0 < m * n + 1)  
        S(i).E = Eo * (1 + a);
        S(i).ENERGY = 1;
        plot(S(i).xd, S(i).yd, '+');
        hold on;
    end
end

S(n + 1).xd = sink.x;
S(n + 1).yd = sink.y;
plot(S(n + 1).xd, S(n + 1).yd, 'x');

% First Iteration - Initializing the first iteration
countCHs = 0;
rcountCHs = 0;
cluster = 1;
countCHs;
rcountCHs = rcountCHs + countCHs;

flag_first_dead = 0;  % Initialize flag for the first dead node
DEAD = zeros(rmax + 1, 1);  % Initialize array to store number of dead nodes per round
DEAD_A = zeros(rmax + 1, 1);  % Initialize array to store number of dead advanced nodes per round
DEAD_N = zeros(rmax + 1, 1);  % Initialize array to store number of dead normal nodes per round

for r = 0:1:rmax
    disp(['Current round: ', num2str(r)]);
    if (mod(r, round(1/p)) == 0)
        for i = 1:1:n
            S(i).G = 0;
            S(i).cl = 0;
        end
    end
    
    % Number of dead nodes
    dead = 0;
    dead_a = 0;
    dead_n = 0;
    
    % Counter for bits transmitted to Base Station and to Cluster Heads
    packets_TO_BS = 0;
    packets_TO_CH = 0;
    
    clf; % Clear the current figure
    hold on; % Hold the plot
    
    for i = 1:1:n
        if (S(i).E <= 0)
            plot(S(i).xd, S(i).yd, 'red .');
            dead = dead + 1;
            if (S(i).ENERGY == 1)
                dead_a = dead_a + 1;
            end
            if (S(i).ENERGY == 0)
                dead_n = dead_n + 1;
            end
            hold on;    
        end
        if S(i).E > 0
            S(i).type = 'N';
            if (S(i).ENERGY == 0)  
                plot(S(i).xd, S(i).yd, 'o');
            end
            if (S(i).ENERGY == 1)  
                plot(S(i).xd, S(i).yd, '+');
            end
            hold on;
        end
    end
    plot(S(n + 1).xd, S(n + 1).yd, 'x');
    
    if (dead == 1)
        if (flag_first_dead == 0)
            first_dead = r;
            flag_first_dead = 1;
        end
    end
    countCHs = 0;
    cluster = 1;
    for i = 1:1:n
        if (S(i).E > 0)
            temp_rand = rand;     
            if ((S(i).G) <= 0)
                if (temp_rand <= (p / (1 - p * mod(r, round(1 / p)))))
                    countCHs = countCHs + 1;
                    packets_TO_BS = packets_TO_BS + 1;
                    
                    S(i).type = 'C';
                    S(i).G = round(1 / p) - 1;
                    plot(S(i).xd, S(i).yd, 'k*');
                    
                    distance = sqrt((S(i).xd - (S(n + 1).xd))^2 + (S(i).yd - (S(n + 1).yd))^2);
                    cluster = cluster + 1;
                    
                    if (distance > do)
                        S(i).E = S(i).E - ((ETX + EDA) * (4000) + Emp * 4000 * (distance^4)); 
                    end
                    if (distance <= do)
                        S(i).E = S(i).E - ((ETX + EDA) * (4000) + Efs * 4000 * (distance^2)); 
                    end
                end     
            end
        end
    end
    
    for i = 1:1:n
        if (S(i).type == 'N' && S(i).E > 0)
            if (cluster - 1 >= 1)
                min_dis = sqrt((S(i).xd - S(n + 1).xd)^2 + (S(i).yd - S(n + 1).yd)^2);
                min_dis_cluster = 1;
                for c = 1:1:cluster - 1
                    temp = min(min_dis, sqrt((S(i).xd - S(c).xd)^2 + (S(i).yd - S(c).yd)^2));
                    if (temp < min_dis)
                        min_dis = temp;
                        min_dis_cluster = c;
                    end
                end
                if (min_dis > do)
                    S(i).E = S(i).E - (ETX * (4000) + Emp * 4000 * (min_dis^4)); 
                end
                if (min_dis <= do)
                    S(i).E = S(i).E - (ETX * (4000) + Efs * 4000 * (min_dis^2)); 
                end
                if (min_dis > 0)
                    S(min_dis_cluster).E = S(min_dis_cluster).E - ((ERX + EDA) * 4000); 
                end
            end
        end
    end
    
    hold off;
    countCHs;
    rcountCHs = rcountCHs + countCHs;
    
    DEAD(r + 1) = dead;
    DEAD_N(r + 1) = dead_n;
    DEAD_A(r + 1) = dead_a;
    
    % Plot the graph for every round
    xlabel('X');
    ylabel('Y');
    title(['Random Sensor Network - Round: ', num2str(r)]);
    grid on;
    pause(0.1); % Pause to display the figure
end
