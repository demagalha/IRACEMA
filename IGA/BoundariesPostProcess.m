function autovector = BoundariesPostProcess(autovector,boundaries)
boundaries = sort(boundaries,'ascend');
[sz1,sz2] = size(autovector);
ndof = sz1+numel(boundaries);


for i=1:numel(boundaries)
	if boundaries(i) == 1
		autovector = [zeros(1,sz2); autovector(1:end,:)];
    elseif boundaries(i) == ndof
		autovector = [autovector(1:end,:); zeros(1,sz2)];
	else
		autovector = [autovector(1:boundaries(i)-1,:); zeros(1,sz2); autovector(boundaries(i):end,:)];
	end
end

end
