%% Initiate NPP
clc;clearvars;clear all;
savedfiles = pwd;

NPPStation = readtable("C:\Users\Admin\Desktop\nuclear-wind\NPPstation.xlsx");
NPPStation = NPPStation(NPPStation.status == "Operational",:);

NPPStation.Distance = deg2km(distance('gc',3.1319,101.6841,NPPStation.latitude,NPPStation.longitude));
NPPStation.Bearing = azimuth(3.1319,101.6841,NPPStation.latitude,NPPStation.longitude);

NPPStation = NPPStation(strcmp(NPPStation.status,'Operational'),:);
NPPStation = NPPStation(NPPStation.Distance<4500,:);
[~,idxNPP] = unique(NPPStation.powerplant);
NPPStation(idxNPP,:) = [];

NPPStation.powerplant = arrayfun(@(x) regexprep(x, '[^a-zA-Z]', ''),NPPStation.powerplant);

[~,idxNPPunique] = unique(NPPStation.powerplant);
NPPStation = NPPStation(idxNPPunique,:);

%% Choose NPP
% % 
NPPStation = NPPStation([11],:);%,3,5,11,13,15
load("C:\Users\Admin\Desktop\nuclear-wind\WorldMap.mat");
Location = arrayfun(@(x) regexprep(x, '[^a-zA-Z]', ''),NPPStation.powerplant);
Location = string(Location)+"constant";

Lat = NPPStation.latitude;
Lon = NPPStation.longitude;


Altitude(1:height(NPPStation)) = 20;
Time = datetime(2022,1,1):days(1):datetime(2023,12,31);
Time = Time.';

str_time = "";

for t = 1:numel(Time)
    str_time(t,:) ="T"+string(year(Time(t)))+sprintf("%02d",month(Time(t)))+sprintf("%02d",day(Time(t)))+sprintf("%02d",hour(Time(t)));
end

rs=[];
rt = [];
rd = [];

FullAccident = [];
TotalLocation_probabilities = [];
Dose_Ratecoefficient = ((3.89E-16)*(3600))./(1e6);
all_Country = unique(all_geotable.Country);
Dose_Ratecoefficient = 1;
% M_SeasonalProbAccidents =zeros(12,numel(Location));
M_SeasonalProbAccident = zeros(12,numel(Location));


Location_probability = [];
ProbPerLocation = [];

for l = 1:numel(Location)
    clear Data
    files = dir(sprintf("ERA5simulations\\%s\\*.mat",Location(l))); % list all .mat files
    for i = 1:length(files)
        filepath = fullfile(files(i).folder, files(i).name);
        if contains(filepath,str_time)
            temp_data = load(filepath);
            Fieldnames = erase(string(files(i).name),".mat");

            try
                % optional: store each dataset in a structure array
                Data.(Location(l)).(Fieldnames) = temp_data.(Fieldnames);

                % optional: display progress
                fprintf("Loaded file %d of %d: %s\n", i, length(files), files(i).name);
            catch ME
                disp(sprintf("Error : %s",ME.message))
            end
        end


    end

    Dose_Ratecoefficient = ((3.89E-16)*(3600))./(1e6);

    SelectedCountry = "Malaysia";
    Dose_Ratecoefficient = 1;

    [FirstReach,HighestReach,RecordReach,DoseSumTable,RateTable,IdxAffectedregion,IdxUnaffectedregion] = ReachRegion(Data,Location(l),Time.',all_geotable,SelectedCountry,Dose_Ratecoefficient);

    r = RecordReach.(Location(l));
    rat = RateTable.(Location(l));
    d = DoseSumTable.(Location(l));
    idxEmpty = cellfun(@(x) isempty(x),r,'UniformOutput',false);
    idxEmpty=cell2mat(idxEmpty);
    TotalSimulation = height(r);
    r = r(~idxEmpty,:);
    rat = rat(~idxEmpty,:);
    d = d(~idxEmpty,:);

    for num_r = 1:numel(r)
        r{num_r} = addvars(r{num_r},repelem(Location(l),height(r{num_r})).');
    end

    for num_t = 1:numel(rat)
        rat{num_t} = addvars(rat{num_t},repelem(Location(l),height(rat{num_t})).');
    end

    % for num_d = 1:numel(d)
    %     d{num_d} = addvars(d{num_d},repelem(Location(l),height(rd{num_d})).');
    % end

    Location_probability = cellfun(@(x) unique(x(:,9)),r,'UniformOutput',false);
    Location_probability = cat(1,Location_probability{:});
    try
        rs = [rs;r];
        rt = [rt;rat];
        rd = [rd;d];
        % ProbPerLocation.(Location(l)) = groupsummary(Location_probability,"Location");
        % ProbPerLocation.(Location(l)).Probability = (ProbPerLocation.(Location(l)).GroupCount./(TotalSimulation)).*100;
        % TotalLocation_probabilities = [TotalLocation_probabilities;Location_probability];
    catch
        disp("Error");
        TotalLocation_probabilities = TotalLocation_probabilities;

    end


    % Location = string(fieldnames(Data)).';


    for m = 1:12   % 7 & 9 out of memory
        temp_Season = RecordReach(month(RecordReach.TimeStamp)==m,:);
        check_Touch = cellfun(@isempty,temp_Season.(Location(l)));

        if sum(~check_Touch)>=1
            M_SeasonalProbAccident(m,l) = 1;
        end

    end
    % M_SeasonalProbAccidents = [M_SeasonalProbAccidents,M_SeasonalProbAccident]

end

%%
FullAccident = cat(1,rs{:});
FullAccident.StartUTC = FullAccident.UTC - FullAccident.Period;
TotalLocation_probabilities = groupsummary(FullAccident,["StartUTC","Region","Location"]);
TotalLocation_probabilities = groupsummary(TotalLocation_probabilities,["Region","Location"]);
TotalLocation_probabilities.Probability = (TotalLocation_probabilities.GroupCount./(365*numel(Location))).*100;

%%
AccumulatedDose_Regional = groupsummary(cat(1,rd{:}),"Region","max","DOSE00030");


NPP_source = groupsummary(FullAccident,["StartUTC","Var16"],"mean","Period");
% NPP_source = groupsummary(FullAccident(FullAccident.Region == c_dx,:),["StartUTC","Var16"],"mean","Period");
NPP_source = groupsummary(NPP_source,["Var16"],"mean","mean_Period");
NPP_source.Probability = (NPP_source.GroupCount/366)*100;

Mean_timearrivalsource = groupsummary(cat(1,rs{:}),["Var16"],"mean","Period");


EffectiveDose_hr = cat(1,rt{:});
% EffectiveDose_hr =EffectiveDose_hr(EffectiveDose_hr.Region == c_dx,:);
% EffectiveDose_hr = EffectiveDose_hr(EffectiveDose_hr.StartUTC ==EffectiveDose_hr.StartUTC(1) & EffectiveDose_hr.Var7 == "Changjiang",:);
% plot(EffectiveDose_hr.DOSE00030);
% EffectiveDose_hr = groupsummary(EffectiveDose_hr,["StartUTC","Region"],"max",M_SeasonalProbAccident)
%%
figure()

subplot(2,1,1)

all_geotable.Probability = zeros(height(all_geotable),1);
all_geotable.Probability(TotalLocation_probabilities.Region) = TotalLocation_probabilities.Probability;

g = geoplot(all_geotable(all_geotable.Country == "Malaysia",:),ColorVariable="Probability",HandleVisibility="off");
clim([0,max(all_geotable.Probability)+2]);
cm = colormap(1-winter);
c = colorbar;
c.Label.String = 'Probability (%)';


title("Frequency of nuclear accident occurence reached by radionuclide");


subplot(2,1,2)
all_geotable.MaximumDose = zeros(height(all_geotable),1);
all_geotable.MaximumDose(AccumulatedDose_Regional.Region) = AccumulatedDose_Regional.max_DOSE00030;

geoplot(all_geotable(all_geotable.Country=="Malaysia" &  ...
     ~isnan(all_geotable.MaximumDose),:),ColorVariable="MaximumDose",HandleVisibility="off");
colormap(1-winter)
clim([0,max(all_geotable.MaximumDose)])
% clim([0,1]);
c = colorbar;
c.Label.String = "Maximum cumulative effective dose (mSv)";

% hold on
%
% geoscatter(Detector.Latitude,Detector.Longitude,"filled",Marker="o",MarkerEdgeColor="k",MarkerFaceColor="b");
%
% hold on
%
% for d = 1:height(Detector)
%     [lat_temp,lon_temp] = getCoordinates(Detector.Latitude(d),Detector.Longitude(d),50,0:360);
%     g1 = geoplot(lat_temp,lon_temp,LineWidth=2,Color="r",HandleVisibility="off");
%     hold on
% end
% legend("ERMS station","FontSize",12.7);


title("Maximum cummulative effective dose");


%%
figure(Theme="light")
% NPP_source.Var16(6) = Location(end);
% NPP_source.GroupCount(6) = 0;
% NPP_source.mean_mean_Period(6) = 0;
% NPP_source.Probability(6) = 0;


b = bar(NPP_source(1:numel(Location),:).Var16, floor(NPP_source(1:numel(Location),:).Probability.*(366/100)));

% Get bar positions and heights
x = b.XData;
y = b.YData;

% Add label for the first bar
text(x, y, [string(duration(NPP_source(1:numel(Location),:).mean_mean_Period,'Format','dd:hh:mm:ss'))+" days"], ...
    'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center', ...
    'FontSize',12);

title("Probability with mean arrival time for each NPP reached Malaysia from daily simulations in 2024","FontSize",12);
xlabel("Locations","FontSize",13);
ylabel("Number of days","FontSize",13);
ylim([0 40]);


%%
figure()
heatmap(double(M_SeasonalProbAccident),XData= Location,YData=unique(string(month(Time,"name")),"stable"),CellLabelColor="none",ColorLimits=[0,1],XLabel="NPP Location",YLabel="Start Time of accident (Month)",FontSize=10);
% title(["Probability of radionuclide reach Malaysia from nearest NPP"]);

%%
Table_overall = cat(1,rd{:});
Table_overall = groupsummary(Table_overall,["LAT","LON","Region"],"mean","DOSE00030");
Table_overall.Population = all_geotable.Population(Table_overall.Region);
Table_overall.Properties.VariableNames(end-1) = "DOSE00030";

[~,i] = max(Table_overall.DOSE00030);
MaxDose = Table_overall(i,:);
EffectiveDose_hrMax = EffectiveDose_hr(EffectiveDose_hr.LAT == MaxDose.LAT & EffectiveDose_hr.LON == MaxDose.LON,:);
plot(EffectiveDose_hrMax.DOSE00030);
hold on 
yline([1e-6])

 %% ginput 

% [x,y] = ginput(1);
% 
%  for p = 1:height(all_geotable)
%      c = find(inpolygon(y,x,all_geotable.Lon{p},all_geotable.Lat{p}));
%      if c ==1
%          c_dx = p;
%      end
%  end
% 
% 
% NPP_source = groupsummary(FullAccident(FullAccident.Region == c_dx,:),["StartUTC","Var16"],"mean","Period");
% NPP_source = groupsummary(NPP_source,["Var16"],"mean","mean_Period");
% NPP_source.Probability = (NPP_source.GroupCount/366)*100;
% 
% figure(Theme="light")
% % NPP_source.Var16(6) = Location(end);
% % NPP_source.GroupCount(6) = 0;
% % NPP_source.mean_mean_Period(6) = 0;
% % NPP_source.Probability(6) = 0;
% 
% % subplot(2,1,1)
% 
% 
% b = bar(NPP_source.Var16, floor(NPP_source.Probability.*(366/100)));
% 
% % Get bar positions and heights
% x = b.XData;
% y = b.YData;
% 
% % Add label for the first bar
% text(x, y, [string(duration(NPP_source.mean_mean_Period,'Format','dd:hh:mm:ss'))+" days"], ...
%     'VerticalAlignment','bottom', ...
%     'HorizontalAlignment','center', ...
%     'FontSize',12);
% 
% title("Probability with mean arrival time for each NPP reached Malaysia from daily simulations in 2024","FontSize",12);
% xlabel("Locations","FontSize",13);
% ylabel("Number of days","FontSize",13);
% ylim([0 40]);


%% contour plot %{

figure("Theme","light","Units","pixels","Position",[0 0 1062 1540]);
i= 0 ;
str_time = "";
% Range_time = [datetime(2024,6,1),datetime(2024,6,15),datetime(2024,7,1),datetime(2024,7,15),datetime(2024,8,1),datetime(2024,8,15)]';
Range_time = datetime(2024,12,18):days(1):datetime(2024,12,23);

str_time ="";
for t = 1:numel(Range_time)
    str_time(t,:) ="T"+string(year(Range_time(t)))+sprintf("%02d",month(Range_time(t)))+sprintf("%02d",day(Range_time(t)))+sprintf("%02d",hour(Range_time(t)));
end

for l = 1 %[1,2,4,5,8,11]
    for tt =1:numel(str_time)

    i = i+1;
    subplot(3,2,i)

    % Extract table
    T = cat(1, Data.(Location(l)).(str_time(tt)).ldata.Concentration{:});
    RegionAffected = unique(T.Region);
    RegionAffected =RegionAffected(2:end);
    T = groupsummary(T,["LAT" "LON"],"max","Cs137000");
    T = T(T.max_Cs137000>0.1,:);

    geoplot(all_geotable(RegionAffected,:),"r",HandleVisibility='off');
    hold on
    % Extract concentration values
    C = T.max_Cs137000;

    % Compute quartiles
    q = quantile(C, [0.25 0.50 0.75]);

    % Assign category (1 to 4)
    CatLevel = ones(size(C));
    CatLevel(C > q(1)) = 2;
    CatLevel(C > q(2)) = 3;
    CatLevel(C > q(3)) = 4;

    % Plot geoscatter with quartile bins
    geoscatter(T.LAT, T.LON, 40, CatLevel, 'filled',HandleVisibility='off');

    % Custom colormap (4 levels: green → yellow → orange → red)
    cmap = [
        0.00 0.40 0.00   % dark green (lowest)   % light green
        1.00 0.85 0.20   % yellow
        1.00 0.50 0.00   % orange
        1.00 0.20 0.20   % red (highest)
        ];

    cmap_I131 = [
        0.85 0.93 1.00   % very light blue
        0.40 0.70 1.00   % medium blue
        0.35 0.20 0.85   % indigo
        0.60 0.00 0.60   % violet
        0.35 0.00 0.35   % dark purple (highest)
        ];

    colormap(cmap_I131);
    title(sprintf("Date :%s \n Location : %s, Maximum ground deposition : %d Bq/m^{2}",str_time(tt),Location(l),max(T.max_Cs137000)));
    hold on 
    geoscatter(Lat(l),Lon(l),"filled","r",Marker="diamond");
    hold on 
    % geobasemap colorterrain
    geolimits([-5,30],[90,122]);

    % Colorbar
    cb = colorbar;
    cb.Ticks = 1:4;

    cb.TickLabels = { ...
        ['≤ ' num2str(q(1),'%.2e')], ...
        ['≤ ' num2str(q(2),'%.2e')], ...
        ['≤ ' num2str(q(3),'%.2e')], ...
        ['> ' num2str(q(3),'%.2e')]};

    cb.Label.String = "Ground deposition (Bq/m^{2})";
    cb.Label.Rotation = 90;     % vertical label
    cb.Label.FontWeight = 'bold';

    % Ensure colorbar shows blocks properly
    caxis([1 4]);
    
    legend("NPP sites");


    end
end

%}

%% boxplot diagram 
figure()

Middledate = datetime(2024,1:12,15);
str_time ="";
for t = 1:numel(Middledate)
    str_time(t,:) ="T"+string(year(Middledate(t)))+sprintf("%02d",month(Middledate(t)))+sprintf("%02d",day(Middledate(t)))+sprintf("%02d",hour(Middledate(t)));
end

bx_datas=[];

for i = 1:numel(str_time)
    
    FieldN = str_time(i);
    bx_data = log(Data.(Location).(FieldN).ldata.mean_Cs137000);
    bx_datas = [bx_datas,bx_data];

end

bx = boxplot(bx_datas);

% Get the current axes
ax = gca;

% --- CHANGE 2: Format Y-Axis Labels to 10^x ---
% Get current Y-tick values (which are now log10 values)
current_yticks = ax.YTick; 
% Create new labels using LaTeX format for superscript
new_ylabels = arrayfun(@(y) sprintf('10^{%d}', round(y)), current_yticks, 'UniformOutput', false);
% Apply the new labels
set(ax, 'YTickLabel', new_ylabels);


% --- OPTIONAL: Format X-Axis Labels to Show Months ---
% Boxplot defaults to 1, 2, 3... This changes them to Jan, Feb, Mar...
set(ax, 'XTick', 1:numel(Middledate));
set(ax, 'YTickLabel', new_ylabels);
set(ax, 'TickLabelInterpreter', 'tex');
set(ax, 'XTickLabel', string(Middledate, 'MMM'));
xlabel('Date');
ylabel('Concentration Cs^{137}');
title('Monthly Cs-137 Distribution');