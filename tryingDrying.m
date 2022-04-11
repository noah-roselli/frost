%% Constants
clear; 
tic;
global nx ny

nx = 200; ny = 200;
radius_fixed = 50; radius = radius_fixed; max_r = 0;
stopExecution = false; % tells us to stop if the radius = nx or ny
dryLimitIterations = 200; dryIterations = 0;

global orderedPair state neighbors neighborsOn neighborsDry
orderedPair  = 1;
state        = 2;
neighbors    = 3; 
neighborsOn  = 4;
neighborsDry = 5;

limit_iterations = 3000;   % limit the amount of iterations per random walker    
particles = 7500;          % Maximum possible amount of particles
%% Make Playing Field
C = makePlayingField();
%% Make seed particle(s)
stuck = [];
seed = xyToEntry(0,0);

C = iJustGotStuck(C,seed, 1);
stuck = [stuck, seed];

prob_freezing = 0.5;
prob_drying = 0.75;
figure(); hold on; grid on; xlim([-nx, nx]); ylim([-ny, ny]);
dry = [];       % the particles that are part of the dry region
radii = [radius];

%% Looping
while (length(stuck) + length(dry) < particles) && (~stopExecution)
    position = randomWalkerStarting(radius); % Starting Position of random walker (RW)
    num = 0;
    
    while num <= limit_iterations % Limit steps of RW 
        if C{position,neighborsOn} > 0 % Check if any RW neighbors are ice
            r = rand(1);
            if r <= prob_freezing 
                C = iJustGotStuck(C,position,1);
                stuck = [stuck, position];

                max_r = max(max_r, ceil(rad(C{position,orderedPair}(1),C{position,orderedPair}(2))));
                radius = min([max_r + radius_fixed, nx]);

                radii = [radii radius];
                if ((length(radii)) > 100 && (radii(end) == radii(end-100)))
                    stopExecution = true;
                end

                if radius == nx
                    stopExecution = true;
                end
                
                break; % go onto next loop
            elseif (r <= prob_drying)
                C = iJustGotStuck(C,position,2);
                dry = unique([dry, position]);
                break;
            else
                temp_position = walk(C,position); % check where you are going next
                if (C{temp_position,state} == 1) || (C{temp_position,state} == 2)
                    continue; % continue the loop
                else
                    num = num + 1;
                    position = temp_position; % if it is liquid, move to it.
                end
            end
        else
            num = num + 1;
            position = walk(C,position); % check where you are going next
        end
    end
end

%% Drawing
% Draw dry zones
for i = 1:length(dry)
    plot(C{dry(i),orderedPair}(1),C{dry(i),orderedPair}(2), '.r', 'MarkerFaceColor','r')
end

% Draw stuck ice particles
for i = 1:length(stuck)
    plot(C{stuck(i),orderedPair}(1),C{stuck(i),orderedPair}(2), '.b','MarkerFaceColor','b')
end

set(gca,'Color','k');  % Make background black
N = length(dry)+length(stuck);
title(sprintf('p_{freeze} = %2.3f, p_{dry} = %2.3f, N = %4.0f', prob_freezing, prob_drying,N))
t = toc; fprintf('It took %2.3f minutes to complete\n', t/60)
%% functions
function status = checkExecutionStatus(C,stuck,dry)
    global state
    ld = length(dry); ls = length(stuck);

    x = NaN(ld+ls,1);
    y = NaN(ld+ls,1);
    
    for i = 1:ld
        x(i) = C{dry(i),1}(1);
        y(i) = C{dry(i),1}(2);
    end
    
    for i = 1:ls
        x(ld+i) = C{stuck(i),1}(1);
        y(ld+i) = C{stuck(i),1}(2);
    end
    
    k = boundary(x,y);
    
    status = true; %border is dry
    for i = 1:length(k)
        if k(i) > ld
            status = false; % border is not dry
            return;
        end

        if i > 1
            if rad(x(k(i))-x(k(i-1)),y(k(i))-y(k(i-1))) > sqrt(2)
                status = false;
                return;
            end
        end
    end
end

function C = makePlayingField()
    global orderedPair state neighbors neighborsOn neighborsDry nx ny

    [X,Y] = meshgrid(-nx:1:nx,-ny:1:ny);
    lenX = 2*nx + 1; lenY = 2*ny + 1;
    
    num = 1;
    for i = 1:lenY
        for j = 1:lenX
            C{num,orderedPair} = [X(i,j),Y(i,j)];
            C{num,state} = 0;
            C{num,neighbors} = findNeighbors(X(i,j),Y(i,j),false);
            C{num,neighborsOn} = 0;
            C{num,neighborsDry} = 0;
            
            num = num+1;
        end
    end
end
function list = findNeighbors(x,y,diagonal)
    global nx ny
    list = [];
    
    i = 1;
    if x-1 >= -nx %left
        list(i,:) = [x-1,y]; i = i+1;
        if diagonal
            if y-1 >= -ny % left down
                list(i,:) = [x-1,y-1]; i = i+1;
            end
    
            if y+1 <= ny % left up
                list(i,:) = [x-1,y+1]; i = i+1;
            end
        end
    end
    
    if x+1 <= nx % right
        list(i,:) = [x+1,y]; i = i+1;
        if diagonal
            if y-1 >= -ny % left down
                list(i,:) = [x+1,y-1]; i = i+1;
            end
    
            if y+1 <= ny % left up
                list(i,:) = [x+1,y+1]; i = i+1;
            end
        end
    end
    
    if y-1 >= -ny % down
        list(i,:) = [x,y-1]; i = i+1;
    end
    
    if y+1 <= ny % up
        list(i,:) = [x,y+1];
    end
end
function position = randomWalkerStarting(r)
    angle = rand(1)*2*pi;
    
    if ((angle >= 0) && (angle < pi/2))
        position = xyToEntry(ceil(r*cos(angle)),ceil(r*sin(angle)));
    elseif ((angle >= pi/2) && (angle < pi))
        position = xyToEntry(floor(r*cos(angle)),ceil(r*sin(angle)));
    elseif ((angle >= pi) && (angle < 3*pi/2))
        position = xyToEntry(ceil(r*cos(angle)),floor(r*sin(angle)));
    elseif ((angle >= 3*pi/2) && (angle < 2*pi))
        position = xyToEntry(floor(r*cos(angle)),floor(r*sin(angle)));
    end
end
function position = xyToEntry(x,y)
    global nx ny
    position = x + (nx + 1) + (2*nx + 1)*(y + ny);
end
function position = walk(C,pos)
    global orderedPair state nx ny

    x = C{pos,orderedPair}(1); y = C{pos,orderedPair}(2);
    r = rand(1); chances = 0:1/4:1;
    
    if r >= chances(1) && r < chances(2)     % walk right
        xNew = x + 1; yNew = y;
    elseif r >= chances(2) && r < chances(3) % walk up
        xNew = x;     yNew = y + 1;
    elseif r >= chances(3) && r < chances(4) % walk left
        xNew = x - 1; yNew = y;
    else                                     % walk down 
        xNew = x; yNew = y - 1;
    end
    
    % Check if the new position is out of the bounds
    if xNew == nx + 1
        xNew = -nx; % send it in the other side
    end
    
    if xNew == -nx - 1
        xNew = nx;
    end
        
    if yNew == ny + 1
        yNew = -ny; % send it in the other side
    end
    
    if yNew == -ny - 1
        yNew = ny;
    end
    
    newPos = xyToEntry(xNew,yNew); % calculate new position
    if C{newPos,state} == 0
        position = newPos;
    else
        position = walk(C,pos);
    end
end
function C = iJustGotStuck(C,pos, NEW_STATE)
    global state neighbors neighborsOn neighborsDry

    C{pos,state} = NEW_STATE; % Frozen or dry
    buddies = C{pos,neighbors};
    for i = 1:length(buddies)
        position = xyToEntry(buddies(i,1),buddies(i,2));
        if NEW_STATE == 1
            C{position,neighborsOn} = C{position,neighborsOn} + 1; % Add one to the neighborsOn for each neighbor
        else
            C{position,neighborsDry} = C{position,neighborsDry} + 1;
        end
    end
end
function d = rad(x,y)
    d = sqrt(x^2 + y^2);
end