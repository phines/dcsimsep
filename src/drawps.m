function [x,y] = drawps(ps,opt)
% draw the case defined in ps
% usage: [x,y] = drawps(ps,options)
% options are defined in psoptions

%cmap = colormap('jet');
cmap = colormap('Copper');
cmap = cmap(size(cmap,1):-1:1,:);
cmap(64,:) = [1 1 1];
n_color = size(cmap,1);
cmap = cmap((n_color:-1:1),:);
colormap(cmap);

%% process inputs
if nargin<2,
    opt=psoptions;
end
width_base=opt.draw.width;
show_bus_nos=opt.draw.bus_nos;
simple=opt.draw.simple;
fs = opt.draw.fontsize;

%% prep work and constants
GREY = [1 1 1]*0.3;
%ORANGE = [255 127 0]/255;
EPS = 1e-6;
C = psconstants;
ps = updateps(ps);
locs = ps.bus(:,C.bu.locX:C.bu.locY);
x = normalize(locs(:,1));
y = normalize(locs(:,2));
width_min     = width_base/100;
circle_size   = width_base*0.1;
triangle_size = width_base*0.05;
flow_max_factor = 1.8;

if all(all(locs==0))
    error('no bus location data provided');
end
    
n = size(ps.bus,1);
bus_sizes = zeros(n,1);

%% prepare the figure
clf; hold on;
set(gcf,'Color','w'); % set bg to white
axis off;

%% draw the branches
for i = 1:size(ps.branch,1)
    f = ps.bus_i(ps.branch(i,C.br.from));
    t = ps.bus_i(ps.branch(i,C.br.to));
    status = ps.branch(i,C.br.status)==1;
    flow = max(ps.branch(i,C.br.Imag_f:C.br.Imag_t)) * status;
    flow_max = ps.branch(i,C.br.rateA:C.br.rateC)/ps.baseMVA;
    flow_ratio = abs(flow)/flow_max(1);
    X = [x(f) x(t)];
    Y = [y(f) y(t)];
%<<<<<<< HEAD
    width = max(sqrt(flow)*width_base,width_min);
    color_ix = max(min( ceil(flow_ratio*n_color/flow_max_factor), n_color ),1);
%=======
%    width = max(flow*width_base,width_min);
%    color_ix = max(1,min( ceil(flow_ratio/2*n_color), n_color ));
%>>>>>>> 66be317674e4dfb07e029e4b3685e03a3fd0325a
    color = cmap(color_ix,:);
    if flow>flow_max
        drawline(X,Y,color,width,'k');
    elseif flow<EPS
        %drawline(X,Y,GREY,width_min);
    else
        drawline(X,Y,color,width);
    end
    %{
    if simple
        plot(X,Y,'k-');
    elseif ~status || flow==0
        %drawline(X,Y,GREY,width);
    elseif flow<(flow_max(1)+EPS)
        drawline(X,Y,GREY,width);
    elseif flow<(flow_max(2)+EPS)
        drawline(X,Y,'y',width);
    elseif flow<(flow_max(3)+EPS)
        drawline(X,Y,ORANGE,width);
    else
        drawline(X,Y,'r',width);
    end
    %}
    % note the bus size
    bus_sizes(f) = max(bus_sizes(f),width);
    bus_sizes(t) = max(bus_sizes(t),width);
end

%% draw the buses
for i = 1:n
    if simple
        drawcircle(x(i),y(i),bus_sizes(i),'w','k',1);
    else
        drawcircle(x(i),y(i),bus_sizes(i),'w','w',1);
    end
end
drawnow

%% draw the generators
green = [133 251 163]/256;
if ~simple
    for i = 1:size(ps.gen,1)
        P = ps.gen(i,C.ge.P)/ps.baseMVA;
        if P>0
            gen_bus_i = ps.bus_i(ps.gen(i,1));
            gx = x(gen_bus_i);
            gy = y(gen_bus_i);
            diameter = sqrt(abs(P))*circle_size;
            drawcircle(gx,gy,diameter,green,green,1);
            %if ps.bus(gen_bus_i,C.bu.type)==C.REF || ps.gen(i,C.ge.type)==C.REF
            %    text(gx,gy,'R','HorizontalAlignment','center');
            %end
        elseif P<0
            text(gx,gy,'?','HorizontalAlignment','center');
        end
    end
    drawnow
end

%% draw the loads
blue = [96 177 222]/256;
%blue = [64 64 256]/256;
if ~simple
    for i = 1:size(ps.shunt,1)
        P = ps.shunt(i,C.sh.P)/ps.baseMVA*ps.shunt(i,C.sh.factor);
        sh_bus_i = ps.bus_i(ps.shunt(i,1));
        sx = x(sh_bus_i);
        sy = y(sh_bus_i);
        drawtriangle(sx,sy,blue,abs(P)*triangle_size,1);
    end
    drawnow
end

%% show bus nos
if show_bus_nos
    for i=1:n
        text(x(i),y(i),num2str(ps.bus(i,1)),'HorizontalAlignment','center','FontSize',fs);
    end
end

%% Keep the x and y proportional
axis equal;
caxis([0 200]);
colorbar

return

%% functions
function xp = normalize(x)
% simple normalization to the range [0,1]
xp = x - min(x);
xp = xp/max(xp);

% draw a line with a given thickness
function drawline(X,Y,color,width,edgecolor)

dx = X(2) - X(1);
dy = Y(2) - Y(1);
ortho = [-dy dx]/sqrt(dx^2 + dy^2);
w = width/2;

p = [[X(1) Y(1)] + ortho*w;
     [X(1) Y(1)] - ortho*w;
     [X(2) Y(2)] - ortho*w;
     [X(2) Y(2)] + ortho*w;];

if nargin<5
    edgecolor = color;
end
patch(p(:,1),p(:,2),color,'EdgeColor',edgecolor);

% draw a circle
function drawcircle(x,y,diameter,face_color,border_color,alpha)

if nargin<5
    border_color = face_color;
end

angles = (0:pi/12:2*pi)';

r = diameter/2;
dy = sin(angles)*r;
dx = cos(angles)*r;

patch( x+dx, y+dy, face_color ,'EdgeColor',border_color,'FaceAlpha',alpha);

% draw a triangle
function drawtriangle(x,y,color,area,alpha)

cos30 = cos(pi/6);
w = sqrt(2*area/cos30);
h = w*cos30;

dx = [-w/2 0 w/2];
dy = [h/2 -h/2 h/2];

patch( x+dx, y+dy, color,'EdgeColor',color,'FaceAlpha',alpha);

