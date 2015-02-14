function snaps = snapshot_data(data, time_s, is_y) 

[xy, t_end ,num_samp] = size(data);
if(is_y == true)
    snaps = data(2,time_s,:);
else 
    snaps = data(1,time_s,:);
end

snaps = reshape(snaps, length(time_s), num_samp);


end 