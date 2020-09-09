
function [t_sim, du_m0, v_m0, du_m1, v_m1, du_m100, v_m100, fb, idx_contact] = loadSol(S,N,wball)
    t_sim = S.x;
    du_m0 = S.y(2*N+1,:);
    v_m0 = S.y(2*N+2,:);
    du_m1 = S.y(1,:);
    v_m1 = S.y(2,:);
    du_m100 = S.y((2*N-1),:);
    v_m100 = S.y(2*N,:);
    
    % Get information during time where ball & rod are in contact
    idx_contact = find((du_m0 - du_m1) > 0);
    t_spr_contact = t_sim(idx_contact);
    du_ball_contact = du_m0(idx_contact);
    du_left_contact = du_m1(idx_contact);
    fb = (wball^2)*(du_ball_contact - du_left_contact).^(3/2);
    fprintf(' Contact time (model): t_c = %.0f microseconds\n',t_spr_contact(end)*10^6)
   
end
