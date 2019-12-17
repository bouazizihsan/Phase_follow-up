function model_based_deterministic_policy_iteration
% initialization
state = [0,1,2,3,4,5];                      % set of states
action = [-1,1];                            % set of actions
% make all (state,action) pairs
[A,B] = meshgrid(state,action);
C = cat(2,A',B');
D = reshape(C,[],2);                        % D includes all possible combination of (state,action) pairs
% make a new D to change the iteration from [(1,1),(2,1),...] tp [(1,1),(1,-1),(2,1),...]
% jjj = 1;
% for iii = 1:length(state)
%     [X,Y]=meshgrid(state(1,iii),action);
%     D(jjj:jjj+1,:) = [X Y];
%     jjj = jjj + 2;
% end
% initialize the policy
h = ones(1,length(state));
h_old = h;                                  % save a backup to compare later
% initialize the Q-function
Q = zeros(length(state),length(action));    % initial Q can be chosen arbitrariliy
Q_old = Q;                                  % save a backup to compare later
L = 100;                                    % maximum number of iterations
gamma = 0.5;                                % discounting factor
epsilon = 0.001;                            % final error to stop the algorithm (here we have the same epsilon for policy and q-function)
Lpe = 0;                                    % counting the number of policy evaluations
%% deterministic policy-iteration algorithm
for l = 1:L
    for ii = 1:size(D,1)
        [x,u] = getStateActionPair(ii,D);  % get a pair of state,action
        if x == 0 || x == 5
            %disp('terminal states reached!')
        else
        xn = model(x,u);                % find the next state xn=x+u
        r = reward(x,u);                % find the reward
        un = getNextAction(x,h);        % find the next action according to policy
        qn = r + gamma * getQValue(xn,un,Q);
        Q = updateQValue(x,u,Q,qn);
        end
    end
    if abs(sum(sum(Q - Q_old))) < epsilon    % evaluate the Q-function
        h = updatePolicy(Q);        % policy improvement
        disp('Policy updated');
        disp(h);
        Lpe = Lpe + 1;              % count the number of policy improvements
        if abs(sum(h - h_old)) < epsilon
            disp('Epsilon criteria satisfied for the policy!');
            break;
        else
            h_old = h;
            Q = zeros(length(state),length(action));
        end
    else
        Q_old = Q;
    end
end
disp('final results :');
% disp('final Q-function = ');
% disp(Q);
disp('final policy = ');
disp(h);
disp(['Number of updating Q-function = ' num2str(l)]);
disp(['Number of updating policy = ' num2str(Lpe)]);
end
%% this function updates the policy and we use it for policy improvement
function h = updatePolicy(Q)
[~,I] = max(Q,[],2);
h(find(I'==1)) = -1;
h(find(I'==2)) = 1;
end
%% this function gets a pair of state,action from a whole set and a iterator
function [x,u] = getStateActionPair(i,D)
% remember that Q = zeros(length(state),length(action));
x = D(i,1);
u = D(i,2);
end
%% this function calculates the next action using the policy
function un = getNextAction(x,h)
un = h(1,x+1);
end
%% this function gets the Q-function value
function q = getQValue(x,u,Q)
if u == -1
    iu = 1;
else
    iu = 2;
end
ix = x + 1;
q = Q(ix,iu);
end
%% This function updates the Q-function elements
function Q = updateQValue(x,u,Q,q)
if u == -1
    iu = 1;
else
    iu = 2;
end
ix = x + 1;
Q(ix,iu) = q;
end

%% this function is the transition model of the robot
% the inputs are: the current state, and the chosen action
% the output is the next state
% in the deterministic case the model function returns the next state
function f = model(x,u)
if (x <= 4 && x>=1)
    f = x + u;
else
    f = x;
end
end

%% this function is the reward function for the task
% the inputs are: the current state, and the chosen action
% the output is the expected reward
function r = reward(x,u)
if (x == 4 && u == 1)
    r = 5;
elseif (x == 1 && u == -1)
    r = 1;
else
    r = 0;
end
end