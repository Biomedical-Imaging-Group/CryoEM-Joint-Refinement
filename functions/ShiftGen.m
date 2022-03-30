function shifts = ShiftGen(sig_shiftX, sig_shiftY, num_samples)
%% Generates num_samples random shifts in x and y directions
%input sig_shiftX:  the maximum value of shift along x direction
%input sig_shiftY:  the maximum value of shift along y direction
%input num_samples: the number of samples/projections
%output shifts: the corresponding shift for each particle image (num_samples x 2)
% Mona Zehni, July 2018

shiftX = (rand([num_samples,1])-0.5)*sig_shiftX;
shiftY = (rand([num_samples,1])-0.5)*sig_shiftY;
shifts = [shiftX,shiftY];

end

