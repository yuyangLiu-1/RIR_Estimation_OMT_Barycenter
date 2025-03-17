function Room_Visualization(roomsize,sourcepoint,meas_points)
%Visualize the room, source and the measurement poitns 
%It couldn't ensure the uniform distribution of the points.
%Which could be repaired.
X = [0;roomsize(1);roomsize(1);0;0];
Y = [0;0;roomsize(2);roomsize(2);0];
Z = [0;0;0;0;0];
figure;
hold on;
plot3(X,Y,Z,"k","LineWidth",1.5);   % draw a square in the xy plane with z = 0
plot3(X,Y,Z+roomsize(3),"k","LineWidth",1.5); % draw a square in the xy plane with z = 1
set(gca,"View",[-28,35]); % set the azimuth and elevation of the plot
for k=1:length(X)-1
    plot3([X(k);X(k)],[Y(k);Y(k)],[0;roomsize(3)],"k","LineWidth",1.5);
end
grid on
xlabel("X (m)")
ylabel("Y (m)")
zlabel("Z (m)")
plot3(sourcepoint(1),sourcepoint(2),sourcepoint(3),"bx","LineWidth",2)
plot3(meas_points(:,1),meas_points(:,2),meas_points(:,3),"ro","LineWidth",2)
end

