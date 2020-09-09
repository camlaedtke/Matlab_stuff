function [idx_contact, t_contact, fb] = contactdata(du_ball, du_left, t_sim, wball)
    t_c = 1.2900e-04;
    idx_contact = find((du_ball - du_left) > 0);
    t_contact = t_sim(idx_contact);
    du_ball_contact = du_ball(idx_contact);
    du_left_contact = du_left(idx_contact);
    fb = (wball^2)*(du_ball_contact - du_left_contact).^(3/2);
    fprintf('Contact time (model): t_c = %.1f microseconds\n',t_contact(end)*10^6)
    fprintf('Contact time (data): t_c = %.1f microseconds\n\n',t_c*10^6)
end

