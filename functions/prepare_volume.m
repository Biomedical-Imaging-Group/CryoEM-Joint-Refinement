function vol = prepare_volume(vol_name, DL)
%% Reads the volume in .mrc format 

%input vol_name: mrc file path of the volume
%input DL: downsampling factor
%output vol: the volume in the shape of a 3D tensor
% Mona Zehni, Aug 2018

if strcmp(vol_name, 'emd_7795')
    [vol_ld, sz] = volGen('emd_7795_54.mrc', DL);
    vol = zeros(84,84,84);
    vol(:,14:71,12:73) = vol_ld;
    
elseif strcmp(vol_name, 'emd_3400') || strcmp(vol_name, 'test')
    [vol_ld,sz] = volGen('emd_3400_86.mrc', DL);
    vol = zeros(90, 90, 90);
    vol(:,12:79,4:87) = vol_ld;
    vol = vol/max(vol(:));
    
elseif strcmp(vol_name, 'test_3400') || strcmp(vol_name, 'test_7795')
    [vol, ~] = volGen([vol_name '.mrc'], DL);
 
else
    error('The volume does not exist.')
    
end


end