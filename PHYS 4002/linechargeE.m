function field = linechargeE(location, charge_location)
    % Calculate the Electric field at a location
    % From a point charge at charge_location
    distance_r = norm(location - charge_location);
    unit_r = (location - charge_location)/distance_r;
    % units where q/4pi epsilon_0 = 1
    field = unit_r*1/distance_r^2;
end

