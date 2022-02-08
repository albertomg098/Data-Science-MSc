%% First of all we initialize all the variables

N=6;
pi=1/4;
pr=1/5;
pi_emergency=1/20;

%% CREATE TRANSITION AND REWARD MATRICES
%ACTION 1.................................

%Transition Matrix Related to Action 1
prob=0;
transition_1=zeros(N);

for i=1:(N+1)
    for j=1:(N+1)
        prob=0;
        for k=0:N
            prob=prob+binopdf(k,(i-1),pi)*binopdf((j-1)-(i-1)+k,N-(i-1),pr);
        end
        transition_1(i,j)=prob;
    end
end

%Reward Matrix Related to Action 1
reward_1=zeros(N);

for i=1:(N+1)
    for j=1:(N+1)
        reward_1(i,j)=50*(exp((i-1)/7)-1);
    end
end

%ACTION 2.............................

%Transition Matrix for Related to Action 2
transition_2=zeros(N);

for i=1:(N+1)
    for j=1:(N+1)
        prob=0;
        for k=0:N
            prob=prob+binopdf(k,(i-1),pi_emergency)*binopdf((j-1)-(i-1)+k,N-(i-1),pr);
        end
        transition_2(i,j)=prob;
    end
end

%Reward Matrix related to Action 2
reward_2=zeros(N);

for i=1:(N+1)
    for j=1:(N+1)
        if i==1
            reward_2(i,j)=0;
        else
            reward_2(i,j)=35*(exp(((i-1)-1)/7)-1);
        end
    end
end

%% SPECIFY A GENERIC MULTIDIMENSIONAL REWARD MATRIX and TRANSITION MATRIX

%When computing the different dynamic optimization algorithms it is helpful
%to make a tri-dimensional matrix unifying the two reward matrices so, the access is easier

total_reward=reward_1;
total_reward(:,:,2)=reward_2;

total_transition=transition_1;
total_transition(:,:,2)=transition_2;

%% OPTIMIZATION WITH Q-VALUE VALUE ITERATION

%Initialization
k=0;
Q=zeros([N+1,2]);
Q_new_dyn=Q;
epsilon= 1e-10;
lamda=0.995;
J=zeros([N+1,1]);
J_new_dyn=J;
diff=1000;

tol=epsilon*(1-lamda)/(2*lamda);
%%
tic
%Once we have initialized the algoritm we can get into the loop

while diff>=tol
    k=k+1;
    Q_old_dyn=Q_new_dyn;
    J_old_dyn=J_new_dyn;
    %We compute at each iteration the Matrix of Q-factors
    for i=1:(N+1)
        for a=1:2
            val=0;
            for j=1:(N+1)
                val=val+total_transition(i,j,a)*(total_reward(i,j,a)+(lamda*(max(Q_old_dyn(j,:)))));
            end
            Q_new_dyn(i,a)=val;
        end
    end
    
    %We compute the vector composed by the solutions of the different
    %Bellman optimmality equations
    for i=1:N+1
        J_new_dyn(i)=max(Q_new_dyn(i,:));
    end
    
    %We compute the infinite norm of the difference
    diff=max(J_new_dyn-J_old_dyn);
end

%We specify the optimal policy
optimal_policy_dyn=zeros([N+1,1]);
for i=1:(N+1)
    j=1;
    max_Q=max(Q_new_dyn(i,:));
    while Q_new_dyn(i,j)~=max_Q
        j=j+1;
    end
    optimal_policy_dyn(i)=j;
end
toc
%% %% OPTIMIZATION WITH REINFORCEMENT LEARNING: Q-LEARNING ALGORITHM

%Initialization
k=0;
Q_reinf=Q;
lamda=0.995;
k_max=1e+06;

%%
tic
%We select the initial state arbitrarily
in_state=2;
curr_state=in_state;
next_state=curr_state;

%Define alpha as a function handle
A=3000;
B=6000;
alpha=@(k) A/(B+k);

%Once we have initialized the algoritm we can get into the loop
while k<k_max
    
    curr_state=next_state;
    
    %We create a uniform random variable in order to pick a random action
    %among the two possible ones
    U=rand;
    if U<0.5  
        a=1;
    else
        a=2;
    end
    
    %We simulate the next state, given the current one
    Z=rand;
    count=0;
    for j=1:(N+1)
        count=count+total_transition(curr_state,j,a);
        if Z<count
            next_state=j;
            break
        else
        end
    end
    
    %We update the matrix of Q-Factors according the Reinforcement Learning
    %Rule
    
    Q_reinf(curr_state,a)=((1-alpha(k+1))*Q_reinf(curr_state,a))+alpha(k+1)*...
        (total_reward(curr_state,next_state,a)+(lamda*(max(Q_reinf(next_state,:)))));
    
    %We update the number of iterations
    k=k+1;
end

%We have built by simulations the optimal matrix of Q-Factors and now we have to translate that
%information into an optimal policy the optimal policy
 J_RL=zeros(N+1,1);
for i=1:N+1
        J_RL(i)=max(Q_reinf(i,:));
 end
    
optimal_policy_reinf=zeros([N+1,1]);

for i=1:(N+1)
    j=1;
    max_Q=max(Q_reinf(i,:));
    while Q_reinf(i,j)~=max_Q
        j=j+1;
    end
    optimal_policy_reinf(i)=j;
end
toc