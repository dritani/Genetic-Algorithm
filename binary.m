% % % %
% %
% Binary Genetic Algorithm
% minimizes the objective function designated in ff Before beginning, set all the parameters in parts
% I, II, and III Haupt & Haupt 2003

%%%%%%%%
%%%%%%%%

% I
%
% 1. Average value 2. Best value convergence
% % all ponts -5 to 5 heatmap:
%
% usually goes to 1,1
% occasionally goes to -1,-1
% occasionally goes to 3,3 (rarely)


%
% cluster straight line
% contour plot with path



% II
% a. popsize:
% convergence again
% heatmap
% b. chrolen
% convergence 
% heatmap

% III
% 0 mutation conv
% 0 mutation heatmap
% 0.9 mutation conv
% 0.9 mutation heatmap







%%%%%%%%
%%%%%%%%
clf
clc
clear
clearvars
size=5;
a = linspace(-size,size);
b = linspace(-size,size);
[A,B] = meshgrid(a,b);
C = -10*exp(-((A-1).^2 + 0.25*(B-1).^2)) -5*exp(-((A+1).^2 + 0.25*(B+1).^2)) - 5*exp(-((A-3).^2 +0.25*(B-3).^2));

%____________% I. Setup the GA
ff='testfunction'; % objective function
npar=2; % number of optimization variables
%____________
%           II. Stopping criteria
maxit=100; % max number of iterations
mincost=-9999999; % minimum cost

%_______________
%           III. GA parameters
popsize=12;% set population size
mutrate=.01;% set mutation rate
selection=0.5;% fraction of population kept
nbits=8;% number of bits in each % parameter
Nt=nbits*npar;% total number of bits in a chormosome
keep=floor(selection*popsize);% #population members that survive
%_________________ % Create the initial population
iga=0; % generation counter initialized
pop=round(rand(popsize,Nt)); % random population of % 1s and 0s
par=gadecode(pop,-5,5,nbits); % convert binary to continuous values
cost=feval(ff,par); % calculates population cost using ff
[cost,ind]=sort(cost); % min cost in element 1
par=par(ind,:);pop=pop(ind,:); % sorts population with % lowest cost first
minc(1)=min(cost);% minc contains min of % population
meanc(1)=mean(cost);% meanc contains mean of population

hold on
plot (-1,-1,'co')
plot(1,1,'co')
plot(3,3,'co')
kkk=zeros(maxit,2);
%__________ Iterate through generations
while iga<maxit
iga=iga+1;% increments generation counter
%___________ Pair and mate

M=ceil((popsize-keep)/2);% number of matings
prob=flipud([1:keep]'/sum([1:keep]));% weights chromosomes based upon position in list
odds=[0 cumsum(prob(1:keep))']; % probability distribution function
pick1=rand(1,M);% mate #1
pick2=rand(1,M);% mate #2
% ma and pa contain the indicies of the chromosomes % that will mate
ic=1;
while ic<=M
  for id=2:keep+1
    if pick1(ic)<=odds(id) & pick1(ic)>odds(id-1)
        ma(ic)=id-1; 
    end % if
    if pick2(ic)<=odds(id) & pick2(ic)>odds(id-1) 
        pa(ic)=id-1;
    end % if 
  end % id
    ic=ic+1; 
end % while
%_______ Performs mating using single point crossover
ix=1:2:keep; % index of mate #1
xp=ceil(rand(1,M)*(Nt-1)); % crossover point
pop(keep+ix,:)=[pop(ma,1:xp) pop(pa,xp+1:Nt)];% first offspring
pop(keep+ix+1,:)=[pop(pa,1:xp) pop(ma,xp+1:Nt)];% second offspring
%_____ Mutate the population 

nmut=ceil((popsize-1)*Nt*mutrate);% total number of mutations
mrow=ceil(rand(1,nmut)*(popsize-1))+1;% row to mutate
%nmut=ceil((popsize)*Nt*mutrate);% total number of mutations
%mrow=ceil(rand(1,nmut)*(popsize));% row to mutate
mcol=ceil(rand(1,nmut)*Nt);% column to mutate

for ii=1:nmut
    pop(mrow(ii),mcol(ii))=abs(pop(mrow(ii),mcol(ii))-1); % toggles bits
end % ii

%____________ The population is re-evaluated for cost
par(2:popsize,:)=gadecode(pop(2:popsize,:),-5,5,nbits);% decode
cost(2:popsize)=feval(ff,par(2:popsize,:));
%________ Sort the costs and associated parameters
[cost,ind]=sort(cost);
par=par(ind,:); pop=pop(ind,:);
%_________ Do statistics for a single nonaveraging run
minc(iga+1)=min(cost);
meanc(iga+1)=mean(cost);
kkk(iga,1)=par(1,1); 
kkk(iga,2)=par(1,2);
%___________ Stopping criteria

    thecolor=[1-iga/maxit 0 0];
k=0;

for k=1:popsize
plot(par(k,1),par(k,2),'Color',thecolor,'Marker','x')
end
if iga==10

end


if iga>maxit | cost(1)<mincost
    break
end
[iga cost(1)];







end %iga
hold off
axis([-5 5 -5 5])
xlabel('X')
ylabel('Y')
%xlabel('function value')
%ylabel('function value')
day=clock;
disp(datestr(datenum(day(1),day(2),day(3),day(4),day(5), day(6)),0))
disp(['optimized function is ' ff])
format short g
disp(['popsize = ' num2str(popsize) ' mutrate = ' num2str(mutrate) ' # par = ' num2str(npar)])
disp(['#generations=' num2str(iga) ' best cost=' num2str(cost(1))])
disp(['best solution'])
disp([num2str(par(1,:))])
disp('continuous genetic algorithm')
figure(2)
iters=0:length(minc)-1;
plot(iters,minc,iters,meanc,'r-');
xlabel('Generation number');ylabel('Function value');
legend('lowest', 'average')





