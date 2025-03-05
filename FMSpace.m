function [AreaA,AreaB,AreaIntersact,FigureofMerit] = FMSpace(Contour1,Contour2,time)

wgs84 = wgs84Ellipsoid("km");
u = unique(Contour1.UTC);
warning off

% figure

for d = time
    ConcCell{d,:} = Contour1(Contour1.UTC==u(d),:);
    [k,Area(d,:)] = convhull(ConcCell{d,:}.Fun_LON,ConcCell{d,:}.Fun_LAT);
    PointConvA = [ConcCell{d}.Fun_LAT(k,:),ConcCell{d}.Fun_LON(k,:)];
    PolyA = polyshape(PointConvA);
end

% geoplot(ConcCell{d,:}.Fun_LAT(k,:),ConcCell{d,:}.Fun_LON(k,:));
AreaA = areaint(PolyA.Vertices(:,1),PolyA.Vertices(:,2),wgs84);
% 
% hold on

u = unique(Contour2.UTC);

for d = time
    ConcCell{d,:} = Contour2(Contour2.UTC==u(d),:);
    [k,Area(d,:)] = convhull(ConcCell{d,:}.Fun_LON,ConcCell{d,:}.Fun_LAT);
    PointConvB = [ConcCell{d}.Fun_LAT(k,:),ConcCell{d}.Fun_LON(k,:)];
    PolyB = polyshape(PointConvB); 
end

% geoplot(ConcCell{d,:}.Fun_LAT(k,:),ConcCell{d,:}.Fun_LON(k,:));
AreaB = areaint(PolyB.Vertices(:,1),PolyB.Vertices(:,2),wgs84);
% 
% hold on 

IntersactPoly= intersect(PolyA,PolyB);

% geoplot(IntersactPoly.Vertices(:,1),IntersactPoly.Vertices(:,2),"k");
try
    AreaIntersact = areaint(IntersactPoly.Vertices(:,1),IntersactPoly.Vertices(:,2),wgs84);
    % text(mean(IntersactPoly.Vertices(:,1)), mean(IntersactPoly.Vertices(:,2)), ...
    %     ['Intersect Area: ', num2str(AreaIntersact),'km^2'], 'Color', 'k', 'FontSize', 12)

catch
    disp("Area doesn't overlap");

end

figure

plot(PolyA, 'FaceColor', 'g') 
hold on 
plot(PolyB, 'FaceColor', 'r')
hold on 

try
    plot(IntersactPoly, 'FaceColor', 'k')
    hold off

    FigureofMerit = (AreaIntersact./(AreaA + AreaB))*100*2;

    text(mean(IntersactPoly.Vertices(:,1)), mean(IntersactPoly.Vertices(:,2)), ...
        ['FMS: ', num2str(FigureofMerit),"%"], 'Color', 'k', 'FontSize', 12)

catch
    AreaIntersact = 0;
    FigureofMerit = 0;
end

end