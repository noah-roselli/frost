%% PREAMBLE

% Noah Roselli, Ricky Palaguachi, 
% Divjyot Singh, Joseph Torsiello
% MATH 451 Capstone Project
% Professor Kondic
% 11 May 2022

%% ALGORITHM
tic; % To start 
global nx ny

nx = 200; ny = 200; % size of the domain

radius_fixed = 30; % the buffer region between the furthest particle and where the walkers are released
radius = radius_fixed; % where random walkers are released
max_r = 0; % to track the maximum radius each iteration

stopExecution = false; % tells us to stop if the radius = nx or ny

global orderedPair state neighbors neighborsOn neighborsDry
orderedPair  = 1; % These are global variables that are used to 
state        = 2; %     keep track of specific values in the cell
neighbors    = 3; %     array that holds all of the data that we
neighborsOn  = 4; %     have for each of the grid points.
neighborsDry = 5;


limit_iterations = 3000;   % limit the amount of iterations per random walker    
particles = 5000;          % Maximum possible amount of particles
%% Make Playing Field
C = makePlayingField(); % Sets up the domain and fills it with liquid particles
%% Make seed particle(s)
stuck = []; % initalize an empty array of stuck particles
seed = xyToEntry(0,0); % place an ice particle at the origin
                       % the xyToEntry takes an x,y ordered pair and gives
                       % the location of this entry in the cell array.

C = iJustGotStuck(C, seed, 1); % change the cell array to update the fact that there is an ice particle now
stuck = [stuck, seed]; % append it to the stuck aray

freezing = 0.25; % probability of freezing
drying = 0.125; % probability of drying if it does not freeze

prob_freezing = freezing; % probability of freezing
prob_drying = freezing+(1-freezing)*drying; % probability of drying if it does not freeze

% initialize a figure
figure(); hold on; xlim([-nx, nx]); ylim([-ny, ny]); set(gca,'XColor','none','YColor','none');
dry = [];       % the particles that are part of the dry region
%% Looping
while (length(stuck) + length(dry) < particles) && (~stopExecution)
    % BIG while look that tracks whether or not to continue releasing
    % random walkers
    position = randomWalkerStarting(radius); % Starting Position of random walker (RW)
    num = 0; % traicks the iterations of the random walker
    
    while num <= limit_iterations % Limit steps of RW 
        if C{position,neighborsOn} > 0 % Check if any RW neighbors are ice
            r = rand(1); % random number
            if r <= prob_freezing % check if it freezes
                % IF IT FREEZES:
                C = iJustGotStuck(C,position,1); % update cell array
                stuck = [stuck, position]; % append to stuck

                % update the radius at which radom walkers are released
                max_r = max(max_r, ceil(rad(C{position,orderedPair}(1),C{position,orderedPair}(2))));
                radius = min([max_r + radius_fixed, nx]);

                % check if you should stop sending out random walkers
                if radius == nx
                    stopExecution = true;
                end
                
                break; % go onto next loop
            elseif (r <= prob_drying) % check if it dries instead
                % IF IT DRIES:
                C = iJustGotStuck(C,position,2); % update the cell array
                dry = unique([dry, position]); % parse for only the unique dry values
                break;
            else
                % IF IT DOES NOT FREEZE OR DRY:
                temp_position = walk(C,position); % check where you are going next
                if (C{temp_position,state} == 1) || (C{temp_position,state} == 2)
                    % If the next location is frozen or dry
                    continue; % continue the loop
                else
                    % otherwise walk
                    num = num + 1;
                    position = temp_position; % if it is liquid, move to it.
                end
            end
        else
            % THIS EXECUTES IF: none of the neighbors of the RW are ice
            num = num + 1;
            position = walk(C,position); % check where you are going next
        end
    end
end
t = toc;
%% Drawing
xy = zeros(length(dry),2); % XY VALUES OF DRY
for i = 1:length(dry)
    xy(i,:) = [C{dry(i),orderedPair}(1),C{dry(i),orderedPair}(2)];
end
% Draw dry zones
hex = "#83b1ab";
scatter(xy(:,1), xy(:,2), 5, 'o', ...
    'MarkerFaceColor', hex, ...
    'MarkerEdgeColor', hex, 'Color', hex)

xy = zeros(length(stuck),2); % XY VALUES OF ICE
for i = 1:length(stuck)
    xy(i,:) = [C{stuck(i),orderedPair}(1),C{stuck(i),orderedPair}(2)];
end
% Draw ice particles
hex = "#004840";
scatter(xy(:,1),xy(:,2),5,'o', ...
    'MarkerFaceColor', hex, ...
    'MarkerEdgeColor', hex, 'Color', hex)

N = length(dry)+length(stuck); % calculate total length
title(sprintf('p_{freeze} = %2.3f, p_{dry} = %2.3f, N = %4.0f', freezing, drying, N))
fprintf('It took %2.3f minutes to complete\n', t/60)
%% functions
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

    [X,Y] = meshgrid(-nx:1:nx,-ny:1:ny); % create a grid that is (2 nx + 1)(2 ny + 1) units
    lenX = 2*nx + 1; lenY = 2*ny + 1;
    
    num = 1;
    for i = 1:lenY
        for j = 1:lenX
            C{num,orderedPair} = [X(i,j),Y(i,j)]; % the ordered pair of the lattice point
            C{num,state} = 0; % the current state of lattice (initially liquid)
            C{num,neighbors} = findNeighbors(X(i,j),Y(i,j),false); % find neighbors in the cardinal directions
            C{num,neighborsOn} = 0; % the number of neighbors that are ice
            C{num,neighborsDry} = 0; % the number of neighbors that are dry
            
            num = num+1;
        end
    end
end
function list = findNeighbors(x,y, diagonal)
    % x, y: ordered pair you are finding neighbors of
    % diagonal: boolean that is whether or not you want to include
    % diagonals

    global nx ny
    list = []; % the list of neighbors
    
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
    % r: radius at which the walkers are released
    angle = rand(1)*2*pi; % find a random angle
    
    if ((angle >= 0) && (angle < pi/2)) % in first quadrant
        position = xyToEntry(ceil(r*cos(angle)),ceil(r*sin(angle))); % push out to grid point
    elseif ((angle >= pi/2) && (angle < pi)) % in second quadrant
        position = xyToEntry(floor(r*cos(angle)),ceil(r*sin(angle))); % push out to grid point
    elseif ((angle >= pi) && (angle < 3*pi/2)) % in third quadrant
        position = xyToEntry(floor(r*cos(angle)),floor(r*sin(angle))); % push out to grid point
    elseif ((angle >= 3*pi/2) && (angle < 2*pi)) % in fourth quadrant
        position = xyToEntry(ceil(r*cos(angle)),floor(r*sin(angle))); % push out to grid point
    end
end
function position = xyToEntry(x,y)
    % takes in an x,y ordered pair
    global nx ny
    position = x + (nx + 1) + (2*nx + 1)*(y + ny); 
    % calculates the position of it in the cell array
end
function position = walk(C,pos)
    % takes in the cell array and the current position and chooses a new
    % position for it to go tod
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
    % Takes in cell array, current position and the new state
    global state neighbors neighborsOn neighborsDry

    C{pos,state} = NEW_STATE; % Frozen or dry
    buddies = C{pos,neighbors}; % the list of neighbors
    for i = 1:length(buddies) % iterate over neighbors
        position = xyToEntry(buddies(i,1),buddies(i,2)); % find cell array position of neighbor
        if NEW_STATE == 1
            C{position,neighborsOn} = C{position,neighborsOn} + 1; % Add one to the neighborsOn for each neighbor
        else
            C{position,neighborsDry} = C{position,neighborsDry} + 1; % Add one to the neighborsDry for each neighbor
        end
    end
end
function d = rad(x,y) % just a simple distance function
    d = sqrt(x^2 + y^2);
end