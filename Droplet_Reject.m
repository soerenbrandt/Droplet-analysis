%Particle Rejection Script

%Intialize 'Objects' in order to save computation time
T = Test.Objects;
n = size(T);



%For loop to determine if a droplet should be ommitted from the analysis
omitions = [];
for i = 1:(n(1))
    %Remove all droplets that only have one recorded position and therefore
    %are not considered moving droplets
    [N_positions,pos] = size(T(i).position);
    if N_positions <= 2 
        omitions(end + 1) = i;
        continue
    end
    
    %Quantify the distance that a droplet has traveled using the standard
    %euclidian distance formula
    XDistance = (T(i).position(end,1) - T(i).position(1,1));
    YDistance = (T(i).position(end,2) - T(i).position(1,2));
    Total_Distance = sqrt((XDistance + YDistance)^2);
    
    %Determine the size of droplet to see if it has traveled outside its
    %radius
    meansize = mean(T(i).size);
    
    %Determine if a droplet has traveled outside its radius
    if Total_Distance < meansize
        
        omitions(end + 1) = i;
        
    end
 
end

%Remove the droplets that have not traveled over the treshold distance
T(omitions) = [];

%------------------------End of omission script ---------------------------


    