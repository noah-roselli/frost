load('./dataSets/data-1-boundary.mat');
shrinks = [1];
for i = 1:length(shrinks)
    [x,y,k,status] = findBoundary(C,dry,stuck,shrinks(i));
    drawImage(C,x,y,k,stuck,dry,shrinks(i))
end
%% Functions
function d = rad(x,y)
    d = sqrt(x^2 + y^2);
end

function drawImage(C,x,y,k,stuck,dry,shrink)
    %% Drawing
    figure(); hold on;
    orderedPair = 1;
    % Draw dry zones
    for i = 1:length(dry)
        plot(C{dry(i),orderedPair}(1),C{dry(i),orderedPair}(2), '.r', 'MarkerFaceColor','r')
    end
    
    % Draw stuck ice particles
    for i = 1:length(stuck)
        plot(C{stuck(i),orderedPair}(1),C{stuck(i),orderedPair}(2), '.b','MarkerFaceColor','b')
    end
    
    plot(x(k),y(k),'w');
    
    set(gca,'Color','k');  % Make background black
    title(sprintf('Shrink Factor = %1.2f',shrink))
end

function [x,y,k,status] = findBoundary(C,dry,stuck,shrinkFactor)
    %% Creating Boundary
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
    
    k = alphaShape(x,y,shrinkFactor);
    k = boundaryFacets(k);
    status = true;
%     status = true; %border is dry
%     for i = 1:length(k)
%         if k(i) > ld
%             status = false; % border is not dry
%             return;
%         end
%     
%         if i > 1
%             if rad(x(k(i))-x(k(i-1)),y(k(i))-y(k(i-1))) > sqrt(2)
%                 status = false;
%                 return;
%             end
%         end
%     end
end