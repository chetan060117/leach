clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Field Dimensions - x and y maximum (in meters)
xm = 100;
ym = 100;

% x and y Coordinates of the Sink
sink.x = 0.5 * xm;
sink.y = 0.5 * ym;

% Number of Nodes in the field
n = 100;

% Optimal Election Probability of a node to become a cluster head
p = 0.05;

% Energy Model (all values in Joules)
% Initial Energy 
Eo = 5;
% Eelec=Etx=Erx
ETX = 50 * 0.000000001;
ERX = 50 * 0.000000001;
% Transmit Amplifier types
Efs = 10 * 0.000000000001;
Emp = 0.0013 * 0.000000000001;
% Data Aggregation Energy
EDA = 5 * 0.000000001;

% Values for Heterogeneity
% Percentage of nodes that are advanced
m = 0.005;
% Alpha
a = 1;

% Maximum number of rounds
rmax = 1200;
%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%

% Computation of do
do = sqrt(Efs / Emp);

% Preallocation
S = struct('xd', zeros(1, n), 'yd', zeros(1, n), 'G', zeros(1, n), 'type', char(zeros(1, n)), 'E', zeros(1, n), 'ENERGY', zeros(1, n));
XR = zeros(1, n);
YR = zeros(1, n);
C = struct('xd', [], 'yd', [], 'distance', [], 'id', []);
X = [];
Y = [];
STATISTICS = struct('DEAD', zeros(1, rmax));
DEAD = zeros(1, rmax);
DEAD_N = zeros(1, rmax);
DEAD_A = zeros(1, rmax);
ALIVE_NODE = zeros(1, rmax);
PACKETS_TO_CH = zeros(1, rmax);
PACKETS_TO_BS = zeros(1, rmax);
energy_res = zeros(n, rmax);
energy_moy = zeros(n, rmax);

% Creation of the random Sensor Network
figure(1);
for i = 1:n
    S(i).xd = rand(1, 1) * xm;
    XR(i) = S(i).xd;
    S(i).yd = rand(1, 1) * ym;
    YR(i) = S(i).yd;
    S(i).G = 0;
    % initially there are no cluster heads, only nodes
    S(i).type = 'N';

    temp_rnd0 = i;
    % Random Election of Normal Nodes
    if temp_rnd0 >= m * n + 1
        S(i).E = Eo;
        S(i).ENERGY = 0;
        plot(S(i).xd, S(i).yd, 'o');
        hold on;
    end
    % Random Election of Advanced Nodes
    if temp_rnd0 < m * n + 1
        S(i).E = Eo * (1 + a);
        S(i).ENERGY = 1;
        plot(S(i).xd, S(i).yd, '+');
        hold on;
    end
end

S(n + 1).xd = sink.x;
S(n + 1).yd = sink.y;
plot(S(n + 1).xd, S(n + 1).yd, 'x');

% First Iteration
figure(1);
title('LEACH-SWDN Simulation: Node Distribution and Dynamics');

for r = 0:rmax
    fprintf('Iteration: %d\n', r);
    % Operation for epoch
    if mod(r, round(1 / p)) == 0
        for i = 1:n
            S(i).G = 0;
            S(i).cl = 0;
        end
    end

    % Number of dead nodes
    dead = 0;
    % Number of dead Advanced Nodes
    dead_a = 0;
    % Number of dead Normal Nodes
    dead_n = 0;
    % Alive node
    alive = 0;

    % counter for bit transmitted to Base Station and to Cluster Heads
    packets_TO_BS = 0;
    packets_TO_CH = 0;
    % counter for bit transmitted to Base Station and to Cluster Heads
    % per round
    PACKETS_TO_CH(r + 1) = 0;
    PACKETS_TO_BS(r + 1) = 0;

    for i = 1:n
        % checking if there is a dead node
        if S(i).E <= 0
            plot(S(i).xd, S(i).yd, 'red .');
            dead = dead + 1;
            if S(i).ENERGY == 1
                dead_a = dead_a + 1;
            end
            if S(i).ENERGY == 0
                dead_n = dead_n + 1;
            end
            hold on;
        end
        if S(i).E > 0
            S(i).type = 'N';
            if S(i).ENERGY == 0
                plot(S(i).xd, S(i).yd, 'o');
            end
            if S(i).ENERGY == 1
                plot(S(i).xd, S(i).yd, '+');
            end
            hold on;
        end
    end
    plot(S(n + 1).xd, S(n + 1).yd, 'blue x');

    STATISTICS(r + 1).DEAD = dead;
    DEAD(r + 1) = dead;
    DEAD_N(r + 1) = dead_n;
    DEAD_A(r + 1) = dead_a;
    alive = n - dead;
    ALIVE_NODE(r + 1) = alive;

    % When the first node dies
    if dead == 1
        if flag_first_dead == 0
            first_dead = r;
            flag_first_dead = 1;
        end
    end

    countCHs = 0;
    cluster = 1;

    for i = 1:n
        if S(i).E > 0
            % Average energy
            energy_moy(i, r + 1) = (1 / r) * sum(energy_res(i, 1:(r + 1)));
            % Generate a random number based on energy average
            temp_rand = ((energy_moy(i) / Eo) * rand(1, 1));
            if S(i).G <= 0
                % k = (number of live nodes in the network * percentage of cluster heads)
                k = ALIVE_NODE(r + 1) * p;

                % Election of Cluster Heads
                if temp_rand <= (n * (S(i).E / Eo) / (n - k * mod(r, round(n / k))))
                    countCHs = countCHs + 1;
                    packets_TO_BS = packets_TO_BS + 1;
                    PACKETS_TO_BS(r + 1) = packets_TO_BS;

                    S(i).type = 'C';
                    S(i).G = round(1 / k) - 1;
                    C(cluster).xd = S(i).xd;
                    C(cluster).yd = S(i).yd;
                    plot(S(i).xd, S(i).yd, 'k*');

                    distance = sqrt((S(i).xd - S(n + 1).xd)^2 + (S(i).yd - S(n + 1).yd)^2);
                    C(cluster).distance = distance;
                    C(cluster).id = i;
                    X(cluster) = S(i).xd;
                    Y(cluster) = S(i).yd;
                    cluster = cluster + 1;

                    % Calculation of Energy dissipated
                    if distance > do
                        S(i).E = S(i).E - ((ETX + EDA) * (4000) + Emp * 4000 * (distance^4));
                        energy_res(i, r + 1) = S(i).E;
                    end
                    if distance <= do
                        S(i).E = S(i).E - ((ETX + EDA) * (4000) + Efs * 4000 * (distance^2));
                        energy_res(i, r + 1) = S(i).E;
                    end
                end
            end
        end
    end

    % Update the plot and pause
    drawnow;
    pause(0.1);
end
