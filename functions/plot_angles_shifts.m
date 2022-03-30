% plot_angles_shifts
interm_result_path = '../interm_res/';

figure
for j = 1:3
    [~, I] = sort(params.angles(:, j));
    subplot(2, 3, j)
    plot(params.angles(I, j), '.');hold on;plot(angles_init(I, j), 'x');
    subplot(2, 3, j + 3)
    tmp = params.angles(I, j) - angles_init(I, j);
    tmp = min(tmp, 2 * pi - tmp);
    tmp1 = 2 * pi + tmp;
    tmp = tmp .* (abs(tmp) < abs(tmp1)) + tmp1 .* (abs(tmp1) < abs(tmp));
    hist(tmp, 100)
end
saveas(gcf, [interm_result_path 'rot' num2str(iter) '.png'])

figure
for j = 1:2
    [~, I] = sort(params.shifts(:, j));
    subplot(2, 2, j)
    plot(params.shifts(I, j), '.');hold on;plot(shifts_init(I, j), 'x');
    subplot(2, 2, j + 2)
    tmp = params.shifts(I, j) - shifts_init(I, j);
    hist(tmp, 100)
end
saveas(gcf, [interm_result_path 'shifts' num2str(iter) '.png'])

close all
