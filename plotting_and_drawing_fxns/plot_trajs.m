[ntrajs,ntpts,xy]=size(trajs);
hold on;
for i=1:ntrajs
	plot(trajs(i,:,1),trajs(i,:,2));
end
hold off;