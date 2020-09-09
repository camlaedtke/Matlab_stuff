function mydeq = springdde(t, y, w0, wball, n, beta)
    % sets up n-by-n matrix 
    
    y2bm0 = reshape(y(1:2*n), 2, n); % reshape and leave out last two
    
    ypos = y2bm0(1,:); % take position vectors (rod only)
    dydt = y2bm0(2,:); % take velocity vectors (rod only)
    ydiffpos = diff(ypos); % differentiate position vectors
    
    if(t == 0)
        vcm = 0;
    else
        vcm = ypos./t; % the difference between a given mass’s velocity 
                        % and the (current) COM velocity
    end
    
    
    d2ydt2 = -(w0^2)*([0, ydiffpos] - [ydiffpos,0]) - 2*beta.*(dydt-vcm);
    
    % if position of ball is greater than position of first mass
    if y(2*n+1) >= ypos(1)                
        % apply a force proportional to how far away it is 
        fb = (wball^2)*(y(2*n+1) - ypos(1)).^(3/2);              
        % add that force to the force already on the first mass        
        d2ydt2(1) = d2ydt2(1) + fb;
        % make the force on the ball *negative* that force 
        d2ydt2(n+1) = -fb; 
    % otherwise, make the force on ball zero
    else
        d2ydt2(n+1) = 0;
    end
    
    % set velocity of first ball equal to 
    dydt(n+1) = y(2*n+2);
    y2bm0final = [dydt; d2ydt2];
    mydeq = y2bm0final(:);

end