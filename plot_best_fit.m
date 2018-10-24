[nbpts,ntpts,xy]=size(best_fit);
hold on;
for i=1:nbpts
	plot(best_fit(i,:,1),best_fit(i,:,2));
end
hold off;
