function interpolate_data = interpcontour(Lat,Lon,Value,Resolution,Ask_Plot)

arguments
    Lat double
    Lon double
    Value double
    Resolution double
    Ask_Plot logical = 0;
end

lat_v = min(Lat):Resolution:max(Lat);
lon_v = min(Lon):Resolution:max(Lon);

[grid_lon,grid_lat] = meshgrid(lon_v,lat_v);

Vq = griddata(Lon,Lat,Value, lon_v, lat_v.',"linear");
interpolate_data = table(grid_lat(:),grid_lon(:),Vq(:),VariableNames=["Latitude","Longitude","Value"]);
idx_nan = isnan(table2array(interpolate_data(:,end)));
interpolate_data(idx_nan,:) = [];

if Ask_Plot
    % sch = scatter(grid_lon(:), grid_lat(:), 35, Vq(:), 'fill');
    geoscatter(interpolate_data,"Latitude","Longitude","filled","ColorVariable","Value");
end

end