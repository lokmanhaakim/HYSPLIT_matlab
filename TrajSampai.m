%%

TimeStamp_reach = [];
FieldN = string(FieldN);

figure()
for ff = FieldN.' %FieldN(contains(FieldN,"2023"),:).'

    if sum(DataTraj.(Location).(ff).Trajectory.Latitude <8) && sum(DataTraj.(Location).(ff).Trajectory.Longitude <110) && DataTraj.(Location).(ff).Trajectory.Longitude(end) >95
        geoplot(DataTraj.(Location).(ff).Trajectory,"Latitude","Longitude")
        hold on 
        TimeStamp_reach = [TimeStamp_reach;ff];
    end

end

% Remove the 'T' prefix
TimeStamp_reach = erase(TimeStamp_reach , "T");
legend(string(TimeStamp_reach))
%%
% Convert to datetime
TimeStamp_reach = datetime(TimeStamp_reach, "InputFormat", "yyyyMMddHH");
TimeStamp_reach = TimeStamp_reach(TimeStamp_reach.Hour == 0,:);
TimeStamp_reach = unique(TimeStamp_reach);
TimeStamp_reach(day(TimeStamp_reach)>26 & month(TimeStamp_reach) == 12,:)=[];
% TimeStamp_reach.Year = 2022;
Group_TimeStamp = groupsummary(table(TimeStamp_reach),"TimeStamp_reach");
figure()
stem(TimeStamp_reach, ones(size(TimeStamp_reach)), ...
    'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 6, ...
    'MarkerFaceColor', 'r', 'Color', 'b');

title(["Total frequency for "+Location+" : "+(sum(Group_TimeStamp.GroupCount)/(365+365+366)).*100+" %"]);
%%
% Time = RecordReach_Timestamp(~ismember(RecordReach_Timestamp,TimeStamp_reach),:);
% 
% str_time = "";
% 
% for t = 1:numel(Time)
%     str_time(t,:) ="T"+string(year(Time(t)))+sprintf("%02d",month(Time(t)))+sprintf("%02d",day(Time(t)))+sprintf("%02d",hour(Time(t)));
% end
% 
% FieldFailed = str_time;
% 
% figure()
% for ff = FieldFailed.'
% 
%     try
%     % if DataTraj.(Location).(ff).Trajectory.Latitude(end) <12 && DataTraj.(Location).(ff).Trajectory.Longitude(end) <110 && DataTraj.(Location).(ff).Trajectory.Longitude(end) >95
%         geoplot(DataTraj.(Location).(ff).Trajectory,"Latitude","Longitude")
%         hold on 
% 
%     catch
%         continue
%     end
%         % TimeStamp_reach = [TimeStamp_reach;ff];
%     % end
% 
% end
% 
% legend(str_time)