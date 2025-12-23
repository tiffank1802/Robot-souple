% plot_heaviside.m
% Plot Heaviside, its integral and derivative (as implemented)

close all;

% Domain
t = linspace(-1,3,1000);

% Evaluate functions (they are now vectorized)
h = Heaviside(t);
int_h = intHeaviside(t);
dh = difHeaviside(t);

% Detect interactive mode; if not interactive, save figures
isInteractive = false;
try
    isInteractive = usejava('desktop');
catch
    isInteractive = false;
end

% If running headless, prefer the 'gnuplot' toolkit which supports
% offscreen rendering for PNG output in Octave.
if ~isInteractive
    try
        graphics_toolkit('gnuplot');
    catch
        % ignore if not available
    end
end

% Single figure with the three functions
f = figure('visible', isInteractive); hold on;
plot(t,h,'b-','LineWidth',2);
plot(t,int_h,'r-','LineWidth',2);
plot(t,dh,'g-','LineWidth',2);
legend('Heaviside','Integral','Derivative');
title('Heaviside and related functions');
xlabel('t'); ylabel('value');
set(gca,'FontSize',14);

if ~isInteractive
    % In headless environments printing via Octave can fail. As a
    % robust fallback, write the data to a file and use gnuplot to
    % render the PNG directly (gnuplot is typically available on
    % systems with Octave installed).
    datafile = 'heaviside_data.dat';
    D = [t(:), h(:), int_h(:), dh(:)];
    try
        save(datafile, 'D', '-ascii');
        gp_all = sprintf(['set terminal pngcairo size 800,400', char(10), ...
            'set output "heaviside_all.png"', char(10), ...
            'set key left top', char(10), ...
            'set xlabel "t"', char(10), ...
            'set ylabel "value"', char(10), ...
            'plot "heaviside_data.dat" using 1:2 with lines title "Heaviside", "heaviside_data.dat" using 1:3 with lines title "Integral", "heaviside_data.dat" using 1:4 with lines title "Derivative"']);
        fid = fopen('heaviside_all.plt','w'); fprintf(fid,'%s\n',gp_all); fclose(fid);
        system('gnuplot heaviside_all.plt');
    catch
        % If gnuplot isn't available, at least keep the data file for
        % external plotting.
    end
end

% Separate figure for clarity (optional)
f2 = figure('visible', isInteractive);
subplot(3,1,1); plot(t,h,'b-'); title('Heaviside'); ylabel('h(t)');
subplot(3,1,2); plot(t,int_h,'r-'); title('Integral of Heaviside'); ylabel('int_h(t)');
subplot(3,1,3); plot(t,dh,'g-'); title('Derivative (0)'); xlabel('t'); ylabel('dh(t)');
set(gcf,'Position',[100 100 600 800]);

if ~isInteractive
    % Separate PNG using gnuplot
    try
        gp_sep = sprintf(['set terminal pngcairo size 600,800', char(10), ...
            'set output "heaviside_sep.png"', char(10), ...
            'set multiplot layout 3,1 title "Heaviside and relatives"', char(10), ...
            'set xlabel ""', char(10), 'set ylabel "h(t)"', char(10), 'plot "heaviside_data.dat" using 1:2 with lines notitle', char(10), ...
            'set ylabel "int_h(t)"', char(10), 'plot "heaviside_data.dat" using 1:3 with lines notitle', char(10), ...
            'set xlabel "t"', char(10), 'set ylabel "dh(t)"', char(10), 'plot "heaviside_data.dat" using 1:4 with lines notitle', char(10), 'unset multiplot']);
        fid = fopen('heaviside_sep.plt','w'); fprintf(fid,'%s\n',gp_sep); fclose(fid);
        system('gnuplot heaviside_sep.plt');
    catch
    end
end
