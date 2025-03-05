%% Load previous data
clc;clearvars;
folderPath = '';  

matFiles = dir(fullfile(folderPath, '*.mat')); 

if ~strcmp(folderPath,'')

    Data = struct();

    % Loop through each .mat file and load the contents
    for k = 1:27
        % Get the full path of the .mat file
        filePath = fullfile(folderPath, matFiles(k).name);

        % Load the .mat file into a temporary structure
        tempStruct = load(filePath);
        tempStruct = tempStruct.(erase(matFiles(k).name,".mat"));

        % Combine the loaded data into the main structure
        % The field name will be the file name without the '.mat' extension
        [~, fileName] = fileparts(matFiles(k).name);
        Data.(fileName) = tempStruct;
    end

end


%% Location and time setup

load("WorldMap.mat");

currentfile = pwd;

Location = "Tanjung Malim"; 
Location = arrayfun(@(x) regexprep(x, '[^a-zA-Z]', ''),Location);
Lat = 3.6795;      
Lon = 101.5203; 
Time = datetime(2025,2,5); 
Period = 24*15;

Filepath = pwd; 

if exist('Data', 'var')
    disp("Reload existed data...")
else 
    Data = HYSPLITRun(Location,Lat,Lon,Time,Period,all_geotable,currentfile,"",""); % insert folder name at the end of function to save it in a folder
end
%% Time arrival and highest measurement **
Location = string(fieldnames(Data)).';

% Dose Unit = µSv/hr;
Dose_Ratecoefficient = ((3.89E-16)*(3600))./(1e6);

SelectedCountry = "Malaysia";
Dose_Ratecoefficient = 1;

[FirstReach,HighestReach,RecordReach,IdxAffectedregion,IdxUnaffectedregion] = ReachRegion(Data,Location,Time,all_geotable,SelectedCountry,Dose_Ratecoefficient);
%% Calculate probability from each sources
clc;

for l = 1:numel(Location)
    r = RecordReach.(Location(l));
    idxEmpty = cellfun(@(x) isempty(x),r,'UniformOutput',false);
    idxEmpty=cell2mat(idxEmpty);
    TotalSimulation = height(r);
    LocColumn = find(strcmp(r{:}.Properties.VariableNames,'Location'));
    Location_probability = cellfun(@(x) unique(x(:,LocColumn)),r(~idxEmpty),'UniformOutput',false);
    Location_probability = cat(1,Location_probability{:});
    Location_probability = groupsummary(Location_probability,"Location");
    Location_probability.Probability = (Location_probability.GroupCount./TotalSimulation).*100; % bahagi setiap simulations
end
%% Calculate probability from all sources
clc;
All_SummaryTable = [];

% Loop through each location
for l = 1:numel(Location)
    % Extract records for the current location
    r = RecordReach.(Location{l}); % Ensure proper indexing with {}


% Initialize an empty structure for location probabilities
Location_probability = struct();
    % Identify and exclude empty records
    idxEmpty = cellfun(@isempty, r); % Logical array of empty cells
    validRecords = r(~idxEmpty);     % Retain non-empty records

    % Calculate total number of simulations
    TotalSimulation = height(r);
    LocColumn = find(strcmp(r{:}.Properties.VariableNames,'Location'));

    % Extract the 3rd column from valid records and concatenate them
    locationData = cellfun(@(x) unique(x(:, LocColumn)), validRecords, 'UniformOutput', false);
    locationData = vertcat(locationData{:}); % Combine into a single array

    % Store the result in the output structure
    try  
        % Compute the grouped summary
        summaryTable = groupsummary(locationData, "Location");

        % Add probability column
        summaryTable.Probability = (summaryTable.GroupCount ./ TotalSimulation) * 100;
        All_SummaryTable = [All_SummaryTable;summaryTable];
        Location_probability.(Location{l}) = summaryTable;
    
    catch
        Location_probability.(Location{l}) = [];
    end

end

AffectedLocation = unique(All_SummaryTable.Location);

g = table;

for aff = 1:numel(AffectedLocation)
    g.Location(aff,:) = AffectedLocation(aff);
    g.GroupCount(aff,:) = sum(All_SummaryTable(All_SummaryTable.Location == AffectedLocation(aff),:).GroupCount);
end

g.Probability = (g.GroupCount./(TotalSimulation.*l)).*100;

%% Figure of Merit in space
try
    Contour1 = Data2.Hanbit.T2020030100.Concentration;
    Contour2 = Data2.Hanbit.T2024030100.Concentration;
    RecordedFieldName = [];
    s = 0;

    [a,b,c,f] = FMSpace(Contour1,Contour2,1);
    disp("FMS :"+f);

catch
end

%% Contour plot

num_loc =numel(Location);
SelectedTime = Time(1);
SelectedTime = "T"+sprintf("%d%02d%02d%02d",year(SelectedTime),month(SelectedTime),day(SelectedTime),hour(SelectedTime));

videoFile = VideoWriter('ContourPlotVideo.mp4', 'MPEG-4');
open(videoFile); % Open the video writer to prepare for frame writing

for n = 1:num_loc

SelectedData = Data.(Location(n)).(SelectedTime); % set the date
Conc = SelectedData.Concentration;
ldata = SelectedData.ldata;
wgs84 = wgs84Ellipsoid("km");

Area = zeros(height(ldata),1);

figure(Name = "Contour plot")
geobasemap colorterrain
geoplot(all_geotable);
hold on
for i = 1:height(ldata)
    % k = convhull(ldata.Concentration{i}.Fun_LON,ldata.Concentration{i}.Fun_LAT);
    % Area(i) = areaint(ldata.Concentration{i}.Fun_LAT,ldata.Concentration{i}.Fun_LON,wgs84);
    % title(sprintf("UTC : %s - %s",string(SelectedData.Concentration.UTC(1)),string(SelectedData.Concentration.UTC(end))));

    geoplot(all_geotable.Shape(ldata.RegionCell{i}(2:end)),"FaceColor","r");
    hold on
    geoscatter(ldata.Concentration{i},"Fun_LAT","Fun_LON","filled",ColorVariable="Fun_C13700000");
    title(sprintf("Affected population :"+"\n%s"),string(sum(all_geotable.Population(unique(cat(1,ldata.RegionCell{i}(2:end)))))));

    frame = getframe(gcf);  % Capture the figure's current frame
    writeVideo(videoFile, frame); % Write frame to video

    pause(0.01)

end

end

% Close the video file after the loop
close(videoFile);  % Finalize and save the video

disp('Video saved as ContourPlotVideo.mp4');


%% Affected and unaffected region 
num_loc = 1;

for num_loc = 1
    SelectedData = Data.(Location(num_loc)).(SelectedTime);
    ldata = SelectedData.ldata;
    Dose_Ratecoefficient = 1; %7.85e-18;

    IdxAffectedregion = (unique(SelectedData.Concentration.Region,"rows","sorted"));
    IdxAffectedregion = IdxAffectedregion(2:end,:);
    IdxAffectedregion = ismember(1:height(all_geotable),IdxAffectedregion);
    IdxAffectedregion = IdxAffectedregion.';

    fig = figure(Name = "Affected Area");

    l = groupsummary(cat(1,SelectedData.ldata.Concentration{:}),"Region","sum",["Fun_C13700000" "Fun_C13700010" "Fun_C13700020" "Fun_C13700030"]);
    l(1,:) = [];
    all_geotable.dose = zeros(height(all_geotable),1);
    AffectedregionIdx = unique(cat(1,ldata.RegionCell{:}));
    all_geotable.dose(IdxAffectedregion,:) = l.sum_Fun_C13700000*(Dose_Ratecoefficient);
    geoplot(all_geotable(IdxAffectedregion,:),"ColorVariable","dose");
    colormap(1-winter);
    colorbar;

    hold on
    geoplot(all_geotable(~IdxAffectedregion,:),FaceColor="b");
    hold on
    geoplot(SelectedData.climatedata.Fields3D.LAT(1),SelectedData.climatedata.Fields3D.LON(1),"r",Marker="diamond");
    legend(["Affected region" "Unaffected region" "NPP Station" ]);
    geobasemap colorterrain
    title(sprintf("%s\n UTC : %s - %s",Location(num_loc),string(SelectedData.Concentration.UTC(1)),string(SelectedData.Concentration.UTC(end))));
    hold off
    
    FileName = sprintf("%s%s.fig",Location(num_loc),SelectedTime);
    savefig(fig,FileName)

end

%% Time-series dose & accumulated dose

% SelectedData = Data.TanjungMalim.(SelectedTime);
% List = all_geotable.District(unique(cat(1,ldata.RegionCell{:})),:);
% [idxLocation,loc] = listdlg('ListString',List);
% idxLocation = find(strcmp(all_geotable.District,List(idxLocation)));
% ldata = SelectedData.ldata;
% 
% 
% for i = 1:height(ldata)
% 
%     idxReadDose = ldata.RegionCell{i}==idxLocation;
%     Dose(i,:) = median(ldata.Concentration{i}.Fun_C13700000(idxReadDose,:));
%     [Centre_Lat,Centre_Lon] = centroid(polyshape(ldata.Concentration{i}.Fun_LAT,ldata.Concentration{i}.Fun_LON));
% 
%     % if isnan(Dose(i,:)) 
%         % Dose(i,:) = griddata(ldata.Concentration{i}.Fun_LON,ldata.Concentration{i}.Fun_LAT,ldata.Concentration{i}.Fun_C13700000,Centre_Lon,Centre_Lat);
%         if isnan(Dose(i,:))
%             Dose(i,:) = 0 ;
%         end
%     % end
% 
% end
% 
% Total_Dose = cumsum(Dose);
% 
% figure
% plot(0:height(Dose),[0;Dose])
% hold on 
% plot(0:height(Dose),[0; Total_Dose]);
% legend(["Dose","Accumulated dose"]);
% title(["Location :" + all_geotable.District(idxLocation)]);
% hold off
% 
% %% Effective doses 
% clc
% DCFNewborn = 9.23e-18;
% DCFTeenager = 7.95e-18;
% DCFAdult = 7.85e-18;
% 
% DCF = [DCFNewborn;DCFTeenager;DCFAdult];
% d= [1;height(ldata)];
% 
% for f = 1:numel(DCF)
%     f
%     disp("Mean effective day " +d+" : "+SelectedData.ldata.mean_Fun_C13700010(d)*(DCF(f))*1e6 + " +- "+ SelectedData.ldata.std_Fun_C13700010([1 7])*(DCF(f))*1e6)
% 
% end
% %% Mean Concentration (High - Low / Area)
% i = 1;
% (SelectedData.ldata.max_Fun_C13700010(i) - min(SelectedData.ldata.Concentration{i}.Fun_C13700010))./(areaint(SelectedData.ldata.Concentration{i}.Fun_LAT,SelectedData.ldata.Concentration{i}.Fun_LON,wgs84))
% 
% %% Mortality risk based on inhalation 
% Activity = SelectedData.ldata.max_Fun_C13700010(end);
% Risk_Mortality = 2.19e-10;
% Life_span = 1;
% Breathing_rate = 19.2/24;%16.5,19.2
% 
% (Life_span)*(Breathing_rate)*(Activity)*(Risk_Mortality)*(100000)

%% Local function
function [AreaA,AreaB,AreaIntersact,FigureofMerit] = compareContour(Contour1,Contour2,time)

u = unique(Contour1.UTC);

for d = time
    ConcCell{d,:} = Contour1(Contour1.UTC==u(d),:);
    [k,Area(d,:)] = convhull(ConcCell{d,:}.Fun_LON,ConcCell{d,:}.Fun_LAT);
    PointConvA = [ConcCell{d}.Fun_LAT(k,:),ConcCell{d}.Fun_LON(k,:)];
    PolyA = polyshape(PointConvA);
end


u = unique(Contour2.UTC);

for d = time
    ConcCell{d,:} = Contour2(Contour2.UTC==u(d),:);
    [k,Area(d,:)] = convhull(ConcCell{d,:}.Fun_LON,ConcCell{d,:}.Fun_LAT);
    % geoplot(ConcCell{d,:}.Fun_LAT(k,:),ConcCell{d,:}.Fun_LON(k,:));
    % hold on 
    PointConvB = [ConcCell{d}.Fun_LAT(k,:),ConcCell{d}.Fun_LON(k,:)];
    % [~,AreaB] = convhull(PointConvB);
    % PointConvA = array2table(PointConvA);
    PolyB = polyshape(PointConvB); 
    % [c1,c2] = centroid(polyshape(ConcCell{d}.Fun_LAT,ConcCell{d}.Fun_LON)); %,Trajectory = Traj
    % CentroidCell{d,:} = [c1,c2];
end

% geoplot(ConcCell{d,:}.Fun_LAT(k,:),ConcCell{d,:}.Fun_LON(k,:));
% [c1,c2] = centroid(polyshape(Contour1.Fun_LAT,Contour1.Fun_LON))

% DifferencesMetric = [Area cell2mat(CentroidCell)]; 

IntersactPoly= intersect(PolyA,PolyB);
% [~,AreaIntersact] = convhull(IntersactPoly.Vertices);

wgs84 = wgs84Ellipsoid("km");



figure

geoplot(PolyA.Vertices(:,1),PolyA.Vertices(:,2),"b");
AreaA = areaint(PolyA.Vertices(:,1),PolyA.Vertices(:,2),wgs84);
% AreaA = deg2km(AreaA);
text(mean(PolyA.Vertices(:,1)), mean(PolyA.Vertices(:,2)), ...
     ['Area A: ', num2str(AreaA)], 'Color', 'k', 'FontSize', 12)
hold on 

% plot(PolyB, 'FaceColor', 'r')
geoplot(PolyB.Vertices(:,1),PolyB.Vertices(:,2),"r");
AreaB = areaint(PolyB.Vertices(:,1),PolyB.Vertices(:,2),wgs84);
% AreaB = deg2km(AreaB);
text(mean(PolyB.Vertices(:,1)), mean(PolyB.Vertices(:,2)), ...
     ['Area B: ', num2str(AreaB)], 'Color', 'k', 'FontSize', 12)
hold on 


% plot(IntersactPoly, 'FaceColor', 'g')
geoplot(IntersactPoly.Vertices(:,1),IntersactPoly.Vertices(:,2),"k");
AreaIntersact = areaint(IntersactPoly.Vertices(:,1),IntersactPoly.Vertices(:,2),wgs84);
% AreaIntersact = deg2km(AreaIntersact);
text(mean(IntersactPoly.Vertices(:,1)), mean(IntersactPoly.Vertices(:,2)), ...
     ['Intersect Area: ', num2str(AreaIntersact)], 'Color', 'k', 'FontSize', 12)


%Normal axes

figure


plot(PolyA, 'FaceColor', 'g') 
% text(mean(PolyA.Vertices(:,1)), mean(PolyA.Vertices(:,2)), ...
     % ['Area A: ', num2str(AreaA)], 'Color', 'k', 'FontSize', 12)
hold on 
plot(PolyB, 'FaceColor', 'r')
% text(mean(PolyB.Vertices(:,1)), mean(PolyB.Vertices(:,2)), ...
     % ['Area B: ', num2str(AreaB)], 'Color', 'k', 'FontSize', 12)
hold on 
plot(IntersactPoly, 'FaceColor', 'k')
% text(mean(IntersactPoly.Vertices(:,1)), mean(IntersactPoly.Vertices(:,2)), ...
     % ['Intersect Area: ', num2str(AreaIntersact)], 'Color', 'k', 'FontSize', 12)
hold off

FigureofMerit = (AreaIntersact./(AreaA + AreaB))*100*2;

text(mean(IntersactPoly.Vertices(:,1)), mean(IntersactPoly.Vertices(:,2)), ...
     ['FMS: ', num2str(FigureofMerit)], 'Color', 'k', 'FontSize', 12)

end

