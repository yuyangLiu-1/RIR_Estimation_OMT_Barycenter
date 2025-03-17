function meas_cld = meas_cloud_u(center,r0,seed,n,shape)
%This function returns n measurement points coordinates surrounding the center.
%Where the points distributed in a circle with the radius r0
%There are n points in the cloud.
%seed is the fixed random root.

rng(seed)                                        % Fix the random values.
phi = 2*pi*linspace(0,1,n+1)+2*pi/n*rand(1);     % The atimuth phases are random.
phi = phi(1:n);
theta = pi*linspace(0,1,n+1)+pi/n*rand(1);       % Elevation phase
theta = theta(1:n);
meas_cld = zeros(n,3);                           % Empty list for the coordiantes

if shape == 'sphere'
    for i = 1:n 
        meas_cld(i,:)=[
            center(1)+r0*cos(phi(i))*cos(theta(i)) ...
            center(2)+r0*sin(phi(i))*cos(theta(i)) ...
            center(3)*sin(theta(i))];
    end
elseif shape == 'circle'
     for i = 1:n 
        meas_cld(i,:)=[
            center(1)+r0*cos(phi(i)) ...
            center(2)+r0*sin(phi(i)) ...
            center(3)];
    end
end

end
