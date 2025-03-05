function [FirstReach,HighestReach,RecordReach,IdxAffectedregion,IdxUnaffectedregion] = ReachRegion(Data,Location,Time,all_geotable,SelectedCountry,Dose_Ratecoefficient)

arguments
    Data struct
    Location string
    Time datetime
    all_geotable table
    SelectedCountry string
    Dose_Ratecoefficient double = 3.89e-16
end

warning off

clc;
% Time = datetime('2023-01-01'):days(7):datetime('2023-12-31');
FirstReach = table;
HighestReach = table;
RecordReach = table;
for NPP = Location
    FirstReach.(NPP) = deal(cell(numel(Time),1));
    HighestReach.(NPP) = deal(cell(numel(Time),1));
    HighestReach.(NPP) = deal(cell(numel(Time),1));
    for Field = 1:numel(fieldnames(Data.(NPP)))
        FieldN = fieldnames(Data.(NPP));
        FieldN = FieldN{Field};
        SelectedData = Data.(NPP).(FieldN); % set the date
        Conc = SelectedData.Concentration;
        Variable = Conc.Properties.VariableNames(4:end-3);
        ldata = groupsummary(Conc,"UTC",["max","mean","std"],Variable);
        u = unique(Conc.UTC);
        [ConcCell,RegionCell,TrajCell,CentroidCell] = deal(cell(numel(u), 1));
        for d = 1:numel(u)
            ConcCell{d,:} = Conc(Conc.UTC==u(d),:);
            TrajCell{d,:} = SelectedData.Trajectory;
            RegionCell{d,:} = unique(ConcCell{d}.Region);
            RegionCell{d,:} = nonzeros(unique(cat(2,ConcCell{d}.Region)));
            [c1,c2] = centroid(polyshape(ConcCell{d}.Fun_LAT,ConcCell{d}.Fun_LON)); %,Trajectory = Traj
            CentroidCell{d,:} = [c1,c2];
        end
        ldata.Concentration = ConcCell;
        ldata.RegionCell = RegionCell;
        ldata.Trajectory = TrajCell;
        ldata.Centroid = CentroidCell;

        IdxAffectedregion = (unique(SelectedData.Concentration.Region,"rows","sorted"));
        IdxAffectedregion = IdxAffectedregion(2:end,:);
        Affectedregion = all_geotable(IdxAffectedregion,:);
        IdxUnaffectedregion = setdiff(1:height(all_geotable),IdxAffectedregion);
        UnAffectedregion = all_geotable(IdxUnaffectedregion,:);

        % First touch at the seletec country
        p = SelectedData.Concentration.Region;
        Countryidx = find(all_geotable.Country == SelectedCountry);

        FirstArrived = find(ismember(p,Countryidx));
        UTCFirstLand = SelectedData.Concentration(min(FirstArrived),:);
        UTCFirstLand = UTCFirstLand.UTC;

        % RecordReachTemp = unique(Conc(:,["UTC" "Region"]));
        RecordReachTemp = Conc(ismember(Conc.Region,Countryidx),:);
        RecordReachTemp = RecordReachTemp(:,["UTC" "Region" "Fun_C13700000" "Fun_C13700010" "Fun_C13700020" "Fun_C13700030"]);
        RecordReachTemp = RecordReachTemp(RecordReachTemp.Region ~= 0,:);
        RecordReachTemp.Location = all_geotable.District(RecordReachTemp.Region,:);
        RecordReachTemp.Shape = all_geotable.Shape(RecordReachTemp.Region,:);
        RecordReachTemp.Period = RecordReachTemp.UTC - Conc.UTC(1); 

        % Highest during the seventh day in the seletected country
        AllTouchedLand = SelectedData.Concentration(FirstArrived,:);
        [~,idxHighDose] = max(AllTouchedLand.Fun_C13700000);
        RegionHigh = AllTouchedLand.Region(idxHighDose);

        try
            FirstLand = SelectedData.Concentration(SelectedData.Concentration.UTC == UTCFirstLand & ismember(SelectedData.Concentration.Region,Countryidx) ,:);
            l = groupsummary(FirstLand,"Region","median","Fun_C13700000");
            FirstReach.(NPP){Field} = table(repelem(FirstLand.UTC(1),numel(all_geotable.District(unique(FirstLand.Region),:))).',string(all_geotable.District(unique(FirstLand.Region),:)),(l.median_Fun_C13700000)*(3.89E-16),all_geotable.Population(unique(FirstLand.Region),:),VariableNames=["Time","Location","Dose","Population"]);
            HighestReach.(NPP){Field} = table(AllTouchedLand.UTC(idxHighDose),all_geotable.District(RegionHigh,:),AllTouchedLand.Fun_C13700000(idxHighDose)*(Dose_Ratecoefficient),all_geotable.Population(RegionHigh,:),VariableNames={'Time','Location','Dose','Population'});
            RecordReach.(NPP){Field} = RecordReachTemp;
        catch
            FirstReach.(NPP){Field} = [];
            HighestReach.(NPP){Field} = [];
            RecordReach.(NPP){Field} = [];
        end

    end

    FirstReach.TimeStamp = Time.';
    FirstReach = movevars(FirstReach,"TimeStamp","Before",1);
    HighestReach.TimeStamp = Time.';
    HighestReach = movevars(HighestReach,"TimeStamp","Before",1);
    RecordReach.TimeStamp = Time.';
    RecordReach = movevars(RecordReach,"TimeStamp","Before",1);

end




end