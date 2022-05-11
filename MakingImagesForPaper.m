r = 5;
figure(); hold on; grid on;
for i = 1:25
    [og, position] = randomWalkerStarting(r);
    plot(og(1), og(2),'ok','MarkerFaceColor','k');
    plot(position(1), position(2),'or','MarkerFaceColor','r');
end

a = linspace(0,2*pi);
plot(r.*cos(a), r.*sin(a),'k','LineWidth',2);
set(gca,'FontSize',18);
xticks([-r:1:r]); yticks([-r:1:r]); axis square;

function [og, position] = randomWalkerStarting(r)
    angle = rand(1)*2*pi;
    og = [r*cos(angle), r*sin(angle)];

    if ((angle >= 0) && (angle < pi/2))
        position = [ceil(r*cos(angle)),ceil(r*sin(angle))];
    elseif ((angle >= pi/2) && (angle < pi))
        position = [floor(r*cos(angle)),ceil(r*sin(angle))];
    elseif ((angle >= pi) && (angle < 3*pi/2))
        position = [floor(r*cos(angle)),floor(r*sin(angle))];
    elseif ((angle >= 3*pi/2) && (angle < 2*pi))
        position = [ceil(r*cos(angle)),floor(r*sin(angle))];
    end
end