% PROBLEM 1
x = [14,15,16,16,16,22,22,24,24,25,25,25,25,25];
avg_of_squares = sum(x.^2)/length(x);
square_of_avg = (sum(x)/length(x))^2;
avg = sum(x)/length(x);
avgs = avg.*ones(1,length(x));
deltas = x - avgs;
avg_deltas = sum(deltas.^2)/length(deltas);
sqrt(avg_deltas);