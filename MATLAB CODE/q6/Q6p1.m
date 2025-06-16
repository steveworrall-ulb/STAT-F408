studentnumber = 123456;
for d = 2:30
    dimension = d;
    simulateGibbssamplerN2025    % This will display a new figure for each d
    % sgtitle(['Gibbs/MH output for d = ', num2str(d)]); % (if you want a supertitle)
    % disp(['Showing result for d = ', num2str(d)]);
    saveas(gcf, sprintf('myfigure_d%02d.jpg', d));
    % pause;   % Wait for you to press a key before showing the next (or use pause(1) for 1 sec)
    close(gcf);  % Close current figure to tidy up
end 
