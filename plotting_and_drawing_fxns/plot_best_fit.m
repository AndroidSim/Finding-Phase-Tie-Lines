function plot_best_fit(varargin)
if (nargin == 2)
    best_tielines=varargin{1};
    figure;
    plot_tielines(best_tielines);
    best_surface=varargin{2};
    figure;
    plot_surface(best_surface);    
elseif (nargin == 1)
    best_trajs=varargin{1};
    [nbpts,ntpts,xy]=size(best_trajs);
    hold on;
    for i=1:nbpts
	    plot(best_trajs(i,:,1),best_trajs(i,:,2));
    end
    hold off;
else
    error('number of arguments not 1 or 2');
    break
end