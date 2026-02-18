%% load data
clearvars;clc;
cd("C:\Users\Admin\Desktop\nuclear-wind");
NPPStation = readtable("NPPstation.xlsx");
NPPStation = NPPStation(NPPStation.status == "Operational",:);

NPPStation.Distance = deg2km(distance('gc',3.1319,101.6841,NPPStation.latitude,NPPStation.longitude));
NPPStation.Bearing = azimuth(3.1319,101.6841,NPPStation.latitude,NPPStation.longitude);

NPPStation = NPPStation(strcmp(NPPStation.status,'Operational'),:);
NPPStation = NPPStation(NPPStation.Distance<4500,:);
[~,idxNPP] = unique(NPPStation.powerplant);
NPPStation(idxNPP,:) = [];
cd("C:\Users\Admin\Desktop\nuclear-wind");

NPPStation.powerplant = arrayfun(@(x) regexprep(x, '[^a-zA-Z]', ''),NPPStation.powerplant);

NPPStation.Elevation = repelem(50,height(NPPStation)).';

[~,idxNPPunique] = unique(NPPStation.powerplant);
NPPStation = NPPStation(idxNPPunique,:);

%%
Location = NPPStation.powerplant;
Location = arrayfun(@(x) regexprep(x, '[^a-zA-Z]', ''),Location);
Lat = NPPStation.latitude;
Lon = NPPStation.longitude;
% Time = datetime(2020,1,1):days(1):datetime(2024,12,31);
Range_time = datetime(2022,1,1):days(0.5):datetime(2024,12,31);
Period = 24*5;
numClusters = 6;
idx_NPP = 11;%3,5,11,13,15
%%

% folderPath = sprintf('ERA5Trajectories/%d',year(Range_time(1)));
folderPath = 'ERA5Trajectories/20222024';

matFiles = dir(fullfile(folderPath, '*.mat'));

DataTraj = struct();
i = 0;
store_idx = [];

for namefiles = 1:numel(matFiles)

    if contains(matFiles(namefiles).name,string(NPPStation.powerplant([idx_NPP]))) %,
        i = i+1;
        store_idx(i,:) = namefiles;
    end
end

% Loop through each .mat file and load the contents
for k = store_idx.'
    % Get the full path of the .mat file
    filePath = fullfile(folderPath, matFiles(k).name);

    % Load the .mat file into a temporary structure
    tempStruct = load(filePath);
    tempStruct = tempStruct.(erase(matFiles(k).name,".mat"));

    % Combine the loaded data into the main structure
    % The field name will be the file name without the '.mat' extension
    [~, fileName] = fileparts(matFiles(k).name);
    DataTraj.(fileName) = tempStruct;

    fprintf("Loaded file of : %s\n", matFiles(k).name);
end


%% Classify to cell form
clc;warning off;
FullSummary = table;
FieldNumber = fieldnames(DataTraj);
Location = string(fieldnames(DataTraj));

All_Trajectories = [];


for SelectedLocation = 1:numel(fieldnames(DataTraj))

    str_time = "";
    for t = 1:numel(Range_time)
        str_time(t,:) = "T"+string(year(Range_time(t)))+sprintf("%02d",month(Range_time(t)))+sprintf("%02d",day(Range_time(t)))+sprintf("%02d",hour(Range_time(t)));
    end

    LocationL = string(Location(SelectedLocation));

    FieldN  = string(fieldnames(DataTraj.(LocationL)));
    DataTraj.(LocationL) = rmfield(DataTraj.(LocationL),FieldN(~ismember(FieldN,str_time)));


    idx_takcukup = [];

    for p =1:numel(fieldnames(DataTraj.(LocationL)))
        FieldN = string(fieldnames(DataTraj.(LocationL)));
        FieldN = FieldN(p);

        if (size(DataTraj.(LocationL).(FieldN).Trajectory,1))<Period+1
            idx_takcukup = [idx_takcukup;p];

        end

    end

    FieldN = string(fieldnames(DataTraj.(LocationL)));
    FieldN = FieldN(idx_takcukup);
    del_file = string(replace(FieldN,"T20","tdump"));
    DataTraj.(LocationL) = rmfield(DataTraj.(LocationL),FieldN);




    SelectedData = DataTraj.(Location(SelectedLocation));
    Field_Date = fieldnames(SelectedData);
    Time = cellfun(@(x) datetime(erase(x,'T'),'InputFormat', 'yyyyMMddHH'),Field_Date);
    FieldN = string(fields(SelectedData)); % input Data
    % trajectories = {};
    clear Trajectories;
    clear OldTrajectories;
    Trajectories =deal(cell(numel(FieldN),1));
    OldTrajectories =deal(cell(numel(FieldN),1));

    for f = 1:numel(FieldN)
        Trajectories{f,:} = [SelectedData.(FieldN(f)).Trajectory.Longitude,SelectedData.(FieldN(f)).Trajectory.Latitude,SelectedData.(FieldN(f)).Trajectory.Height];
    end

    %% Make trajectory sizes same
    sizes = cellfun(@(x) size(x, 1), Trajectories);
    mostCommonSize = max(sizes);

    for t = 1:height(Trajectories)
        OldTrajectories{t} = Trajectories{t};
        if size(Trajectories{t}) ~= mostCommonSize
            NaNTemp = zeros(mostCommonSize - size(Trajectories{t},1),size(Trajectories{t},2));
            Trajectories{t} = [Trajectories{t};NaNTemp];
        end
    end


    %% Step 1: Flatten Trajectories for Pairwise Distance Calculation

    % Each trajectory is reshaped into a 1D array for distance computation.
    numTrajectories = length(Trajectories);
    flattenedTrajectories = zeros(numTrajectories, numel(Trajectories{1}));
    for i = 1:numTrajectories
        flattenedTrajectories(i, :) = reshape(Trajectories{i}, 1, []);
    end

    % Step 2: Compute Pairwise Distance Matrix
    % Calculate Euclidean distance between all trajectories.
    distanceMatrix = pdist(flattenedTrajectories, 'euclidean');

    % Step 3: Perform Hierarchical Clustering
    % Use the "linkage" function to perform agglomerative clustering.
    linkageMatrix = linkage(distanceMatrix, 'ward');

    % Step 4: Visualize the Dendrogram
    % figure;
    % dendrogram(linkageMatrix, 'Labels', arrayfun(@(x) sprintf('Traj %d', x), 1:numTrajectories, 'UniformOutput', false));
    % title('Trajectory Clustering Dendrogram');
    % xlabel('Trajectories');
    % ylabel('Distance');

    % Step 5: Determine Cluster Assignments
    % Choose the number of clusters (e.g., 2) based on the dendrogram or criteria.
    
    colors = lines(numClusters); % Generate distinct colors for clusters

    clusterAssignments = cluster(linkageMatrix, 'maxclust', numClusters);
    clusterID = lines(numClusters);

    % idx_Traj = find(idxActual);
    Traj_Data = table(Trajectories);
    % Display Results
    for i = 1:numTrajectories
        % fprintf('Trajectory %d is in Cluster %d\n', idx_Traj(i), clusterAssignments(i));
        Traj_Data.Cluster(i) = clusterAssignments(i);
    end

    Traj_Data = groupsummary(Traj_Data,"Cluster");
    Traj_Data.Percentage = (Traj_Data.GroupCount / sum(Traj_Data.GroupCount))*100;
    Trajectories_table = table(Trajectories);
    Trajectories_table.OldTrajectories = OldTrajectories;
    Trajectories_table.Cluster = clusterAssignments;
    Trajectories_table.Location = repelem(Location(SelectedLocation),height(Trajectories_table),1);


    for clusterID = 1:numClusters
        Traj = Trajectories_table(Trajectories_table.Cluster == clusterID,:);

        % Traj = Traj(~idxZero,:);
        MeanTraj = cat(3, Traj.Trajectories{:});
        mean_trajectory = mean(MeanTraj, 3);
        % figure(10)
        % Plot each trajectory in the cluster
        for t = 1:height(Traj)
            idxZero = sum(Traj.Trajectories{t}==0,2) == 3;

            trajLon = Traj.Trajectories{t}(~idxZero, 1); % Assuming Var1 contains [Lon, Lat]
            trajLat = Traj.Trajectories{t}(~idxZero, 2);
            % geoplot(trajLat,trajLon, 'Color', colors(clusterID, :), 'LineWidth', 1.5);
            % hold on
        end
        % figure(11)
        %     geoplot(mean_trajectory(:,2),mean_trajectory(:,1))
        %     hold on

    end

    % title("Trajectory based on HYSPLIT approach")

    %%
    FieldN = fields(SelectedData);

    % Initialize the ClassBearing column
    Trajectories_table.ClassBearing = strings(height(Trajectories), 1);

    % Trajectories_table.Bearing = cell2mat(cellfun(@(x) mean(azimuth(x(1,2),x(1,1),x(:,2),x(:,1))),Trajectories_table.Trajectories,'UniformOutput',false));
    Trajectories_table.Bearing = cell2mat(cellfun(@(x) azimuth(x(1,2),x(1,1),x(end,2),x(end,1)),Trajectories_table.Trajectories,'UniformOutput',false));
    Trajectories_table.DiffBearing = cell2mat(cellfun(@(x) abs(azimuth(x(1,2),x(1,1),x(end,2),x(end,1))-azimuth(x(1,2),x(1,1),3.1319,101.6841)),Trajectories_table.OldTrajectories,'UniformOutput',false));


    for i = 1:height(Trajectories)
        Bearing = Trajectories_table.Bearing(i);

        % Classify bearing into compass directions
        if (Bearing >= 0 && Bearing <= 10) || (Bearing > 350 && Bearing <= 360)
            Trajectories_table.ClassBearing(i) = "N";
        elseif Bearing > 10 && Bearing < 80
            Trajectories_table.ClassBearing(i) = "NE";
        elseif Bearing >= 80 && Bearing < 100
            Trajectories_table.ClassBearing(i) = "E";
        elseif Bearing >= 100 && Bearing < 170
            Trajectories_table.ClassBearing(i) = "SE";
        elseif Bearing >= 170 && Bearing < 190
            Trajectories_table.ClassBearing(i) = "S";
        elseif Bearing >= 190 && Bearing < 260
            Trajectories_table.ClassBearing(i) = "SW";
        elseif Bearing >= 260 && Bearing < 280
            Trajectories_table.ClassBearing(i) = "W";
        elseif Bearing >= 280 && Bearing <= 350
            Trajectories_table.ClassBearing(i) = "NW";
        end
    end

    Trajectories_table.ClassBearing = categorical(Trajectories_table.ClassBearing);
    Trajectories_table.StartUTC = Time;
    All_Trajectories = [All_Trajectories;Trajectories_table];

    %% K means
    Centroid = [];
    Trajectories_Ktable = table(Trajectories);
    Trajectories_Ktable.Cluster = zeros(height(Trajectories),1);

    for t = 1:height(Trajectories)
        temp_Traj = Trajectories_Ktable.Trajectories{t}(Trajectories_Ktable.Trajectories{t}(:,1)~=0 & Trajectories_Ktable.Trajectories{t}(:,2)~=0,:);
        % [CLon,CLat] = centroid(polyshape(trajectories.trajectories{t}(:,1),trajectories.trajectories{t}(:,2)));
        Clon = temp_Traj(end,1); %(temp_Traj(end,1) - temp_Traj(1,1))/2;
        Clat = temp_Traj(end,2); %(temp_Traj(end,2) - temp_Traj(1,2))/2;
        Centroid = [Centroid;Clon,Clat];
    end

    % numCLuster = 5;
    idxCluster = kmeans(Centroid,numClusters);

    for i = 1:height(Trajectories_Ktable)
        Trajectories_Ktable.Cluster(i,1) = idxCluster(i);
    end

    numClusters = max(Trajectories_Ktable.Cluster); % Total number of clusters
    colors = lines(numClusters);


    %% Plot trajctory and calculate the percentage
    % groupCounts.Percentage = (groupCounts.GroupCount / sum(groupCounts.GroupCount)) * 100;

    figure
    for clusterID = 1:numClusters
        Traj = Trajectories_Ktable(Trajectories_Ktable.Cluster == clusterID,:);
        % Traj(Traj.Latitude == 0 & Traj.Longitude==0,:) = [];
        MeanTraj = cat(3, Traj.Trajectories{:});
        mean_trajectory = mean(MeanTraj, 3);
        % Plot each trajectory in the cluster
        for t = 1:height(Traj)
            % figure(4)
            % idxZero = (Traj.Trajectories{t}(:, 1) == 0) & (Traj.Trajectories{t}(:, 2) == 0);
            % Traj.Trajectories{t}(idxZero, :) = [];
            % trajLon = Traj.Trajectories{t}(:, 1); % Assuming Var1 contains [Lon, Lat]
            % trajLat = Traj.Trajectories{t}(:, 2);
            %
            % geoplot(trajLat,trajLon, 'Color', colors(clusterID, :), 'LineWidth', 1.5);
            % hold on
            % legend(string(1:numClusters));
        end
        
        geoplot(mean_trajectory(:,2),mean_trajectory(:,1))
        hold on

    end

    % title(sprintf("Clustering from %d trajcetories ",height(Trajectories_table)));
    % legend(string(1:numClusters));


    Trajectories_Ktable.ClassBearing = Trajectories_table.ClassBearing;
    Trajectories_Ktable.Bearing = Trajectories_table.Bearing;

    % Summarize the count of each cluster
    SummaryTraj_Kdata = groupsummary(Trajectories_Ktable, ["ClassBearing"]); %"Cluster",
    disp("Summary Trajectory at "+Location(SelectedLocation));
    SummaryTraj_Kdata.PercentageClustres = (SummaryTraj_Kdata.GroupCount / sum(SummaryTraj_Kdata.GroupCount))*100;
    SummaryTraj_Kdata.Location = repelem(Location(SelectedLocation),height(SummaryTraj_Kdata),1);

    FullSummary = [FullSummary;SummaryTraj_Kdata];
    Date_EachCluster= deal(cell(numClusters,1));
    Date_Consistent = deal(cell(numClusters,1));

    figure(Theme="light")
    for loc = 1:numClusters
        subplot(2,numClusters,loc)
        idx_takcukup = [];
        date_persistent = Trajectories_table.StartUTC(Trajectories_Ktable.Cluster == loc);
        Date_EachCluster{loc,:} = date_persistent;
        xdate = date_persistent;
        xdate.Year = 2020;
        ConsistentSeason = groupsummary(array2table(xdate),"xdate");
        ConsistentSeason = ConsistentSeason.xdate(ConsistentSeason.GroupCount == 3,:);
        
        try
            Date_Consistent{loc,:} = ConsistentSeason;
        catch 
            Date_Consistent{loc,:} = [];
        end

        for p =1:numel(date_persistent)
            % FieldN = string(fieldnames(DataTraj.Fuqing));
            % FieldN = FieldN(p);
            str_time = "T"+string(year(date_persistent(p)))+sprintf("%02d",month(date_persistent(p)))+sprintf("%02d",day(date_persistent(p)))+sprintf("%02d",hour(date_persistent(p)));
            geoplot(DataTraj.(LocationL).(str_time).Trajectory,"Latitude","Longitude","Color","r","Marker","x");

            % if (size(DataTraj.Fuqing.(FieldN).Trajectory,1))<240
            %     idx_takcukup = [idx_takcukup;p];
            % end
            hold on

        end

         geoplot(NPPStation.latitude(idx_NPP),NPPStation.longitude(idx_NPP),Marker="diamond",Color="b",MarkerFaceColor="auto");
            
        % Stem plot is better for individual dates

        subplot(2,numClusters,loc+numClusters)
        stem(date_persistent, ones(size(date_persistent)), ...
            'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 6, ...
            'MarkerFaceColor', 'r', 'Color', 'b');

        xlabel('Tarikh');
        ylabel('Occurrence');
        title(sprintf('Date : %s - %s ',string(min(date_persistent)),string(max(date_persistent))));
        grid on;
        % datetick('x', 'dd-mmm', 'keeplimits'); % Show day-month only

        % Optional: Add month separators
        hold on;
        month_starts = dateshift(min(date_persistent):calmonths(1):max(date_persistent), 'start', 'month');
        for i = 1:length(month_starts)
            xline(month_starts(i), '--', 'Color', [0.5, 0.5, 0.5], 'Alpha', 0.5);
        end
        

        % legend(string(Trajectories_table.StartUTC));
        % title(sprintf("%s - %s",string(date_persistent(1)),string(date_persistent(end))))
    end
    sgtitle(sprintf("Location : %s",LocationL));
    % FullSummary{SelectedLocation} = rows2vars(SummaryTraj_Kdata,"VariableNamesSource","ClassBearing");
    % FullSummary.Properties.VariableNames =

end

data = FullSummary; % Assuming the data is loaded in a table

% Define the compass directions we are interested in
directions = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"];

% Initialize the final table for the output
outputTable = [];

% Get the unique locations
locations = unique(data.Location);

% Loop through each location and populate the table
for i = 1:numel(locations)
    locationData = data(data.Location == locations{i}, :);

    % Initialize a row for the current location
    row = {locations{i}};

    % Loop through each direction (N, NE, E, etc.)
    for j = 1:length(directions)
        % Filter data for the current direction
        bearingData = locationData(locationData.ClassBearing == directions(j), :);

        % If data exists for this bearing, extract the GroupCount and PercentageClusters
        if ~isempty(bearingData)
            groupCount = bearingData.GroupCount;
            percentageClusters =  bearingData.PercentageClustres;
            % Concatenate GroupCount and PercentageClusters in the desired format
            bearingInfo = sprintf('%d (%.4f)', groupCount, percentageClusters);
            % bearingInfo = percentageClusters;
        else
            bearingInfo = 0; % No data for this bearing
        end

        % Append the bearing info to the row
        row = [row, bearingInfo];
    end

    % Add the row to the output table
    outputTable = [outputTable; row];

    
end

% % Set the variable names for the columns
outputTable = array2table(outputTable);
outputTable.Properties.VariableNames = ['Location', directions];

% Display the final table
disp(outputTable);

%% Polarhistogram
% figure(Name="Polar Histogram")
% % clf;
% store_idx = [];
% for l = 1
%     % figure(Name="Polar Histogram"+l)
%     angle_loc = All_Trajectories(All_Trajectories.Location == Location(l),:).DiffBearing;
% 
%     if sum(angle_loc<10) > 366/2
%         Polarfig = polarhistogram(deg2rad(angle_loc));
%         % Polarfig = polarhistogram(deg2rad(All_Trajectories.DiffBearing));
%         Polarfig.Normalization = "pdf";
%         hold on
%         store_idx = [store_idx;l];
%     end
% 
%     % title("Difference angle distribution of trajectories among at all NPP sites");
% end
% 
% legend(Location(store_idx));
% title("Angle distribution of trajectories among at all NPP sites");

%%
for l = 1:numel(Location)
    MeanDiff(l,:) = mean((All_Trajectories(All_Trajectories.Location == Location(l),:).DiffBearing));
end

%% PLot traj

%%
% NorthNation = NPPStation.powerplant(ismember(NPPStation.country,["China" "South Korea" "Taiwan"]),:);
% NorthNation = erase(string(NorthNation)," ");
% % NorthNation = outputTable(ismember(outputTable.Location,NorthNation),:);
%
% WestNation = NPPStation.powerplant(~ismember(NPPStation.country,["China" "South Korea" "Taiwan"]),:);

%% Find similarities between each trajectories
i = 0;
All_day = datetime(2020,1,1):days(1):datetime(2024,12,31);
% num_traj = [1,367];
%
for a = 1:numel(All_day)
    for l = 1:numel(Location)
        Dy = day(All_day(a));
        Mnth = month(All_day(a));
        num_traj{a,l} = find((month(All_Trajectories.StartUTC)==Mnth)&(day(All_Trajectories.StartUTC)==Dy) & All_Trajectories.Location == Location(l));
    end
end


%%
% i = 0;
% Minimum_size = [];
% R_Minimum_size = [];
% 
% for p = 1:numel(num_traj)
%     i = 0;
%     for sm_traj = 1:numel(num_traj{p})
%         i = i+1;
%         idx_Traj = num_traj{p}(sm_traj);
%         Minimum_size(i,:) = height(All_Trajectories.OldTrajectories{idx_Traj});
%     end
%     R_Minimum_size(p,:) = min(Minimum_size);
% end


%%
% for d = 1:numel(num_traj)
%     % diffs = All_Trajectories.OldTrajectories{num_traj{d}(1)}(1:R_Minimum_size(d),1:2) - All_Trajectories.OldTrajectories{num_traj{d}(2)}(1:R_Minimum_size(d),1:2);
%     % euclid_dist = sqrt(sum(diffs.^2, 2));  % Per-point distance
%     % mean_dist(d,:) = mean(euclid_dist);        % Overall similarity score
%     try
%         Dtw_score(d,:) = dtw(All_Trajectories.OldTrajectories{num_traj{d}(1)}(1:R_Minimum_size(d),1:2), All_Trajectories.OldTrajectories{num_traj{d}(2)}(1:R_Minimum_size(d),1:2));
%         dtw_lon = dtw(All_Trajectories.OldTrajectories{num_traj{d}(1)}(1:R_Minimum_size(d),1), All_Trajectories.OldTrajectories{num_traj{d}(2)}(1:R_Minimum_size(d),1));
%         dtw_lat = dtw(All_Trajectories.OldTrajectories{num_traj{d}(1)}(1:R_Minimum_size(d),2), All_Trajectories.OldTrajectories{num_traj{d}(2)}(1:R_Minimum_size(d),2));
%         Dtw_combine(d,:) = (dtw_lon + dtw_lat)/2;
%     catch
%         Dtw_score(d,:) = NaN;
%     end
% end
%
% Dtw_score = reshape(Dtw_score,[],17);
% Dtw_combine = reshape(Dtw_combine,[],17);

%% figure(Name="Minimum Distance")
% figure(Name = "DTW score")
% min_value = min(Dtw_combine(:,15));
% idx_min = find(Dtw_combine == min_value(end));
% idx_min = idx_min(1);
% for t = 1:5
%     geoplot(All_Trajectories.OldTrajectories{num_traj{idx_min}(t)}(1:R_Minimum_size(idx_min),2),All_Trajectories.OldTrajectories{num_traj{idx_min}(t)}(1:R_Minimum_size(idx_min),1))
% hold on
% end
% % geoplot(All_Trajectories.OldTrajectories{num_traj{idx_min}(2)}(1:R_Minimum_size(idx_min),2),All_Trajectories.OldTrajectories{num_traj{idx_min}(2)}(1:R_Minimum_size(idx_min),1))
% title(sprintf("DTW value : %.2f",min_value(end)));
%%

% figure(Name = "DTW score")
% plot(datetime(2020,1,1):datetime(2020,12,31),Dtw_score(:,3));

%% Compute each 5 year
% Sum_dtwscore = [];
% Sum_dtwcombinescore = [];
%
% for d = 1:numel(num_traj)
%     % diffs = All_Trajectories.OldTrajectories{num_traj{d}(1)}(1:R_Minimum_size(d),1:2) - All_Trajectories.OldTrajectories{num_traj{d}(2)}(1:R_Minimum_size(d),1:2);
%     % euclid_dist = sqrt(sum(diffs.^2, 2));  % Per-point distance
%     % mean_dist(d,:) = mean(euclid_dist);        % Overall similarity score
%     Dtw_tempscore = 0;
%     Dtw_tempsumscore = 0;
%     for l = 1 %:numel(num_traj{d})-1
%         Dtw_tempscore = Dtw_tempscore+dtw([All_Trajectories.OldTrajectories{num_traj{d}(l)}(1:R_Minimum_size(d),2),All_Trajectories.OldTrajectories{num_traj{d}(l)}(1:R_Minimum_size(d),1)], [All_Trajectories.OldTrajectories{num_traj{d}(l+1)}(1:R_Minimum_size(d),2),All_Trajectories.OldTrajectories{num_traj{d}(l+1)}(1:R_Minimum_size(d),1)]);
%         Sum_dtwscore(d,:) = Dtw_tempscore/numel(num_traj{d});
%         Dtw_tempsumscore = Dtw_tempsumscore + dtw(All_Trajectories.OldTrajectories{num_traj{d}(1)}(1:R_Minimum_size(d),1), All_Trajectories.OldTrajectories{num_traj{d}(2)}(1:R_Minimum_size(d),1)) + dtw(All_Trajectories.OldTrajectories{num_traj{d}(1)}(1:R_Minimum_size(d),2), All_Trajectories.OldTrajectories{num_traj{d}(2)}(1:R_Minimum_size(d),2));
%         Sum_dtwcombinescore(d,:) = Dtw_tempsumscore/(2*numel(num_traj{d}));
%     end
% end
%
% Sum_dtwscore = reshape(Sum_dtwscore,[],2);
% Sum_dtwcombinescore = reshape(Sum_dtwcombinescore,[],2);

%%
% Sum_fdistance = [];
%
% for d = 1:numel(num_traj)
%     % diffs = All_Trajectories.OldTrajectories{num_traj{d}(1)}(1:R_Minimum_size(d),1:2) - All_Trajectories.OldTrajectories{num_traj{d}(2)}(1:R_Minimum_size(d),1:2);
%     % euclid_dist = sqrt(sum(diffs.^2, 2));  % Per-point distance
%     % mean_dist(d,:) = mean(euclid_dist);        % Overall similarity score
%     fdistance = 0;
%     count = 0;
%     % Sum_fdistance = 0;
%     for l_1 = 1:numel(num_traj{d})
%         for l_2 = l_1+1:numel(num_traj{d})
%         % fdistance = fdistance + frechet_distance([All_Trajectories.OldTrajectories{num_traj{d}(l)}(1:R_Minimum_size(d),2),All_Trajectories.OldTrajectories{num_traj{d}(l)}(1:R_Minimum_size(d),1)]-[All_Trajectories.OldTrajectories{num_traj{d}(l+1)}(1:R_Minimum_size(d),2),All_Trajectories.OldTrajectories{num_traj{d}(l+1)}(1:R_Minimum_size(d),1)]);
%         Difference_dist = frechet_distance([All_Trajectories.OldTrajectories{num_traj{d}(l_1)}(1:R_Minimum_size(d),2),All_Trajectories.OldTrajectories{num_traj{d}(l_1)}(1:R_Minimum_size(d),1)],[All_Trajectories.OldTrajectories{num_traj{d}(l_2)}(1:R_Minimum_size(d),2),All_Trajectories.OldTrajectories{num_traj{d}(l_2)}(1:R_Minimum_size(d),1)]);
%         fdistance = fdistance + Difference_dist; %mean(sqrt(sum(Difference_dist.^2, 2)));
%         % fdistance = fdistance + distance(All_Trajectories.OldTrajectories{num_traj{d}(l)}(1:R_Minimum_size(d),2),All_Trajectories.OldTrajectories{num_traj{d}(l)}(1:R_Minimum_size(d),1), All_Trajectories.OldTrajectories{num_traj{d}(l+1)}(1:R_Minimum_size(d),2),All_Trajectories.OldTrajectories{num_traj{d}(l+1)}(1:R_Minimum_size(d),1));
%         count = count +1;
%         end
%     end
%     % Sum_fdistance(d,:) = sqrt(sum(fdistance.^2))/(l+1);
%     Sum_fdistance(d,:) = fdistance/count;
% end
%
% Sum_fdistance = reshape(Sum_fdistance,[],17);
% Sum_fdistance = Sum_fdistance(1:366,:);

%%
% for m = 1:2
%     figure(m)
%     plot(datetime(2020,1,1):datetime(2020,12,31),Sum_fdistance(:,m),LineStyle="-",LineWidth=1);
%     legend(Location(m));
% end
% [1,2,3,8,13,17]
% %
% % num_loc = 1;
% for num_loc = 1
% 
%     figure(Name="Minimum Distance")
%     min_value = mink(Sum_fdistance(:,num_loc),1);
%     idx_min = find(Sum_fdistance == min_value(end));
%     idx_min = idx_min(1);
%     for t = 1:numel(num_traj{idx_min})
%         geoplot(All_Trajectories.OldTrajectories{num_traj{idx_min}(t)}(:,2),All_Trajectories.OldTrajectories{num_traj{idx_min}(t)}(:,1),Marker="*");
%         hold on
%         % text(All_Trajectories.OldTrajectories{num_traj{idx_min}(t)}(end,2), All_Trajectories.OldTrajectories{num_traj{idx_min}(t)}(end,1), datestr(start_datetime+years(t-1)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right','FontWeight','bold');
%     end
% 
%     title(sprintf("Frechet distance between trajectories at %s : %.2f",Location(num_loc),min_value(end)));
% 
% end

%%
% plot([NewData.Changjiang.T2024010206.Sums.DOSE00100])
%% Frechet distance plot
% clf
% hold off
% 
% for f_d1 = 1:numel(num_traj{idx_min})
%     for f_d2 = f_d1+1:numel(num_traj{d})
%         [fd,cs,fm,pd] = frechetDistance([All_Trajectories.OldTrajectories{num_traj{idx_min}(f_d1)}(:,2),All_Trajectories.OldTrajectories{num_traj{idx_min}(f_d1)}(:,1)],[All_Trajectories.OldTrajectories{num_traj{idx_min}(f_d2)}(:,2),All_Trajectories.OldTrajectories{num_traj{idx_min}(f_d2)}(:,1)]);
%         geoplot(All_Trajectories.OldTrajectories{num_traj{idx_min}(f_d1)}(:,2),All_Trajectories.OldTrajectories{num_traj{idx_min}(f_d1)}(:,1),Marker="*",Color="b");
%         hold on
%         [xs,ys] = plotCouplingSeq([All_Trajectories.OldTrajectories{num_traj{idx_min}(f_d1)}(:,2),All_Trajectories.OldTrajectories{num_traj{idx_min}(f_d1)}(:,1)],[All_Trajectories.OldTrajectories{num_traj{idx_min}(f_d2)}(:,2),All_Trajectories.OldTrajectories{num_traj{idx_min}(f_d2)}(:,1)],cs, pd, curveLabels);
%         geoplot([xs(1) xs(2)],[ys(1) ys(2)],"k");
%         hold on
%     end
% 
% end

% curveLabels = {'28 May 2020', '28 May 2021'};
% WestNation = erase(string(WestNation)," ");
% WestNation = outputTable(ismember(outputTable.Location,WestNation),:);
%
% % Sample data
% categories = WestNation.Properties.VariableNames(2:end); % X-axis categories
% data = WestNation{:,2:end}; % Stacked bar data (rows = categories, columns = groups)
% data = data.';
% % Create the stacked bar chart
% figure;
% b = bar(data, 'stacked');
%
% % Set different colors for each stack
% numStacks = size(data, 2); % Number of stacks (columns)
% % colors = lines(numStacks); % Generate distinct colors
% colors = rand(numStacks,3);
% for k = 1:numStacks
%     b(k).FaceColor = 'flat';
%     b(k).CData = repmat(colors(k, :), size(data, 1), 1); % Apply unique colors
% end
%
% % Add labels and title
% xticks(1:length(categories));
% xticklabels(categories); % Add category labels
% xlabel('Direction');
% ylabel('Frequency');
% title('Frequency of trajectories direction from West NPP stations');
% legend(WestNation.Location, 'Location', 'best'); % Legend for stacks
% grid on;
%
% %%
% FieldN = fieldnames(Data.Changjiang);
% FieldN = string(FieldN);
% for f = 1:numel(FieldN)
%     geoplot(Data.Changjiang.(FieldN(f)).Trajectory.Latitude,Data.Changjiang.(FieldN(f)).Trajectory.Longitude,Marker="x")
%     hold on
% end
%
% title(sprintf("%d Trajectories from %s",f,"Changjiang"));

%% HYSPLIT cluster
%
% FieldN = string(FieldN);
% DiffLat = [];
% DiffLon = [];
% DiffDistance = [];
%
% for i = 1:height(trajectories)
%     Limiter = numel(trajectories.Var1{i}(:,1));
%
%     try
%         DiffLat(i,:) = sum(abs((trajectories.Var1{1}(:,1) - trajectories.Var1{i}(:,1))));
%         DiffLon(i,:) = sum(abs((trajectories.Var1{1}(:,2) - trajectories.Var1{i}(:,2))));
%         DiffDistance(i,:) = mean(distance('gc',trajectories.Var1{1}(:,1),trajectories.Var1{1}(:,2),trajectories.Var1{i}(:,1),trajectories.Var1{i}(:,2)));
%
%     catch
%         DiffLat(i,:) = sum(abs((trajectories.Var1{1}(1:Limiter,1) - trajectories.Var1{i}(:,1))));
%         DiffLon(i,:) = sum(abs((trajectories.Var1{1}(1:Limiter,2) - trajectories.Var1{i}(:,2))));
%         DiffDistance(i,:) = mean(distance('gc',trajectories.Var1{1}(1:Limiter,1),trajectories.Var1{1}(1:Limiter,2),trajectories.Var1{i}(:,1),trajectories.Var1{i}(:,2)));
%     end
%
%
% end

%%
% Step 1: Flatten Trajectories for Pairwise Distance Calculation
% Each trajectory is reshaped into a 1D array for distance computation.
% numTrajectories = length(trajectories);
% % flattenedTrajectories = zeros(numTrajectories, numel(trajectories{1}));
% for i = 1:numTrajectories
%     flattenedTrajectories{i, :} = reshape(trajectories{i}, 1, []);
% end
%
% distanceMatrix = cellfun(@(x) pdist(x, 'euclidean'),trajectories,'UniformOutput',false);
%
% linkageMatrix = cellfun(@(x) linkage(x, 'ward'),distanceMatrix,'UniformOutput',false);
%
% numClusters = 3;
% clusterAssignments = cluster(linkageMatrix{1}, 'maxclust', numClusters);

% %%
%
% numTrajectories = height(trajectories);
%
%
% % Start with each trajectory as its own cluster
% clusters = cell(numTrajectories, 1);
% for i = 1:numTrajectories
%     clusters{i} = {trajectories.Var1{i}};
% end
%
% % Iterative clustering
% TSVs = zeros(numTrajectories - 1, 1); % To store TSV values
% for iteration = 1:numTrajectories - 1
%     minTSV = inf;
%     bestMerge = [];
%
%     % Find the best pair of clusters to merge
%     for i = 1:length(clusters) - 1
%         for j = i + 1:length(clusters)
%             tempClusters = clusters;
%             % Merge clusters i and j
%             tempClusters{i} = [tempClusters{i}, tempClusters{j}];
%             tempClusters(j) = []; % Remove the merged cluster
%
%             % Compute TSV for this merge
%             currentTSV = computeTSV(tempClusters);
%             if currentTSV < minTSV
%                 minTSV = currentTSV;
%                 bestMerge = [i, j];
%             end
%         end
%     end
%
%     % Perform the best merge
%     clusters{bestMerge(1)} = [clusters{bestMerge(1)}, clusters{bestMerge(2)}];
%     clusters(bestMerge(2)) = [];
%     TSVs(iteration) = minTSV; % Store TSV for this iteration
% end
%
% figure;
% plot(1:length(TSVs), TSVs, '-o');
% xlabel('Number of Clusters');
% ylabel('Total Spatial Variance (TSV)');
% title('TSV Curve');
% grid on;
%
%
% function SV = computeSV(trajectory, clusterMean)
%     % trajectory: Nx2 matrix of points
%     % clusterMean: Nx2 matrix of mean trajectory points
%     SV = sum(sum((trajectory - clusterMean).^2, 2)); % Sum of squared differences
% end
%
% function meanTrajectory = computeClusterMean(cluster)
%     % cluster: cell array of trajectories
%     numTrajectories = length(cluster);
%     trajectoryLengths = cellfun(@(x) size(x, 1), cluster);
%     maxLength = max(trajectoryLengths);
%
%     % Interpolate trajectories to the same length
%     interpolatedTrajectories = cellfun(@(x) interp1(1:size(x, 1), x, ...
%                               linspace(1, size(x, 1), maxLength)), cluster, 'UniformOutput', false);
%     interpolatedTrajectories = cat(3, interpolatedTrajectories{:});
%
%     % Compute mean trajectory
%     meanTrajectory = mean(interpolatedTrajectories, 3, 'omitnan');
% end
%
% function TSV = computeTSV(clusters)
%     TSV = 0;
%     for i = 1:length(clusters)
%         cluster = clusters{i};
%         clusterMean = computeClusterMean(cluster);
%         for j = 1:length(cluster)
%             TSV = TSV + computeSV(cluster{j}, clusterMean);
%         end
%     end
% end
%
%


% % Bersihkan data (Buang NaN jika ada)
% idx_clean = ~isnan(all_lats) & ~isnan(all_lons);
% all_lats = all_lats(idx_clean);
% all_lons = all_lons(idx_clean);
% 
% fprintf('Selesai! Total titik terkumpul: %d\n', length(all_lats));
% 
% %% 3. KERNEL DENSITY CALCULATION (Binning & Smoothing)
% fprintf('Menjana Density Map...\n');
% 
% % Setting Resolusi Grid (0.25 darjah ~= 25km)
% grid_res = 0.25;
% 
% % Tentukan had peta (Bounding Box)
% % Kita bagi buffer sikit (+5 darjah) supaya tak rapat sangat ke tepi
% lat_edges = floor(min(all_lats)-5):grid_res:ceil(max(all_lats)+5);
% lon_edges = floor(min(all_lons)-5):grid_res:ceil(max(all_lons)+5);
% 
% % Kira Histogram 2D (Frequency)
% [N, Xedges, Yedges] = histcounts2(all_lons, all_lats, lon_edges, lat_edges);
% 
% % Transpose N untuk plotting
% N = N';
% 
% %% 5. VISUALIZATION MENGGUNAKAN GEOSCATTER
% fprintf('Menjana plot Geoscatter...\n');
% 
% % 1. Cari titik tengah bagi setiap kotak grid (Pixel)
% % Xedges/Yedges adalah garisan tepi, kita nak tengah-tengah kotak
% lon_centers = (Xedges(1:end-1) + Xedges(2:end)) / 2;
% lat_centers = (Yedges(1:end-1) + Yedges(2:end)) / 2;
% 
% % 2. Wujudkan Grid Coordinate untuk setiap pixel dalam N
% % Ini akan hasilkan matrix Lat dan Lon yang sama saiz dengan N
% [LON_GRID, LAT_GRID] = meshgrid(lon_centers, lat_centers);
% 
% % 3. "Ratakan" (Flatten) matrix jadi satu senarai panjang (Vector)
% % Sebab geoscatter cuma terima input vector (Nx1)
% lat_flat = LAT_GRID(:);
% lon_flat = LON_GRID(:);
% N_flat = N(:); % Nilai frequency
% 
% % 4. BUANG NILAI KOSONG (Zero Filtering) - PENTING!
% % Kita tak nak plot titik yang nilainya 0 (kawasan angin tak lalu).
% % Kalau tak buang, peta awak penuh titik biru kosong yang semakkan mata.
% idx_valid = N_flat > 0; 
% 
% lat_plot = lat_flat(idx_valid);
% lon_plot = lon_flat(idx_valid);
% N_plot = N_flat(idx_valid);
% 
% % 5. PLOT GEOSCATTER
% figure('Color','w');
% 
% % geoscatter(Lat, Lon, SaizMarker, WarnaVariable, 'filled')
% % SaizMarker: Kita boleh buat saiz berubah ikut frequency (lagi kerap = lagi besar)
% marker_size = rescale(N_plot, 10, 50); % Skala saiz antara 10 hingga 50 pts
% 
% h = geoscatter(lat_plot, lon_plot, marker_size, N_plot*(100)/numel(fieldnames(DataTraj.Fuqing)), 'filled');
% 
% % 6. KOSMETIK
% % Set Transparency supaya nampak bertindan cantik
% h.MarkerFaceAlpha = 0.6; 
% 
% % Tukar tema peta (Paling cantik untuk publication: 'colorterrain' atau 'streets')
% geobasemap('colorterrain'); 
% % Pilihan lain: 'satellite', 'streets', 'grayland'
% 
% % Warna & Label
% colormap(jet); % Warna panas
% c = colorbar;
% c.Limits = [0,100];
% c.Label.String = 'Frequency (Trajectory Count)';
% title('Trajectory Frequency Map (Geoscatter)');
% 
% % Zoom auto ke kawasan data
% geolimits([min(lat_plot)-2 max(lat_plot)], [min(lon_plot)-2 max(lon_plot)]);
% 
% %%
% 
% % Ambil Log10 supaya jurang nilai tak jauh sangat
% % Tambah 1 supaya tak error log(0)
% N_log = log10(N_plot + 1); 
% 
% figure(Theme="Light");
% % Plot guna nilai LOG
% h = geoscatter(lat_plot, lon_plot, 20, N_log, 'filled');
% 
% colormap(jet); 
% h.MarkerFaceAlpha = 0.5; % Biar nampak lutsinar sikit
% geobasemap('colorterrain');
% c = colorbar;
% 
% % Tukar label colorbar supaya orang faham ni skala log
% c.Label.String = 'Frequency';
% title('Density Map of trajectories from 5 consecutive years');
