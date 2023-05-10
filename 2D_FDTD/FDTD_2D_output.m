%% FDTD 2D Simulation Animation Generator
%  William Poland
clc
clear

%% Customizable Values
%   ****************** file name ******************************************
test_name = "indoor2"; % name generated by C++ code
figure_title = 'Indoor Simulation'; % plot title text
ofile_mod = "take2"; % modifier for output MP4 file
input_file_path = "/Users/williampoland/FDTD/2D/raw/";
output_file_path = "/Users/williampoland/FDTD/2D/";
%   ***********************************************************************

%   ********* simulation values *******************************************
num_field_data = 3; % the number of data points in the text file
num_geo_data = 2; % number of data points in geometry file

% - - - - - - - color axis - - - - - - - - - - - - - - - - - - - - - - - - 
% Note: caxis = color bar axis
use_lim = 1; % turn caxis limits on/off
plot_epsr = 0; % plot dielectric geometry features
plot_ec = 1; % plot conductive geometry features
plot_secondary = 1; % add extra conditions to plot another set of geometry
linear_cmap = 0; % restrict colormap to vary linearly

upp_lim = 0.35; % maximum caxis value 
low_lim = -0.35; % minimum caxis value

map_res = 16; % number of bits to make rows in our color map (resolution)
% tolerance around region we want to make white (blank)
tol_factor = 1e-4; % multiplies absolute average of upper and lower limits
nonlin_factor = 4; % scale intensity of nonlinearity
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   **** values for max field plot ****
max_plot_only = 1; % true to skip animations and only do max field plot
num_tiles = 2;
use_lim2 = 1; % turn caxis limits on/off
linear_cmap2 = 1; % restrict colormap to vary linearly
% its a dB scale
upp_lim2 = 0; % maximum caxis value 
low_lim2 = -70; % minimum caxis value
% - - - - - - - -
%   ***********************************************************************

%% Simulation
input_file = input_file_path + test_name + ".txt";
geometry_file = input_file_path + test_name + "_geometry.txt";

% read in geometry features here from a geometry file
file_id = fopen(geometry_file);
fgetl(file_id); % discard 1st header line

% now extract header values
first_line = fgetl(file_id);
header_info = str2double(split(first_line));
FREQ = header_info(1);
T_MAX = header_info(2);
NUM_X = header_info(3);
NUM_Y = header_info(4);
DEL_T = header_info(5);
DEL_X = header_info(6);
DEL_Y = header_info(7);

T_NEXT = NUM_X * NUM_Y; % number of lines before advancing to next time index

fgetl(file_id); % discard 2nd header line

% now read in geometry features
file_format = '%f %f';
file_size = [num_geo_data Inf];
disp("Reading in data from geometry file...")
geo_data = fscanf(file_id, file_format, file_size);
fclose(file_id);
disp("Finished reading in geometry file data")
geo_data = transpose(geo_data); % transpose array so orientation is correct

% - - - - - - - - - - - - - - - - - - -
% - variables to store geometry data -
% the length of these will be unkown until reading data
ec_xvals = zeros(1,0, 'single');
ec_yvals = zeros(1,0, 'single');
epsr_xvals = zeros(1,0, 'single');
epsr_yvals = zeros(1,0, 'single');
secondary_xvals = zeros(1,0, 'single');
secondary_yvals = zeros(1,0, 'single');
% - - - - - - - - - - - - - - - - - - -

disp("Analyzing geometry data")
for x = 1:NUM_X
    for y = 1:NUM_Y
        curr_row = (double(x)-1)*NUM_X + (double(y)-1) + 1;
        curr_epsr = geo_data(curr_row, 1);
        curr_ec = geo_data(curr_row, 2);

        if (plot_epsr && curr_epsr ~= 1) % dielectrics
            epsr_xvals(end+1) = x;
            epsr_yvals(end+1) = y;
        end
        if (plot_ec && curr_ec ~= 0) % conductors
            ec_xvals(end+1) = x;
            ec_yvals(end+1) = y;
        end
        if (plot_secondary)
            if(curr_epsr >= 2.12 && curr_epsr <= 2.14) % extra geometry
                secondary_xvals(end+1) = x;
                secondary_yvals(end+1) = y;
            end
        end
%         fprintf("x,y,row,eps,ec)=(%d,%d,%d,%d,%f,%f)",x-1,y-1,curr_row, ...
%             curr_epsr,curr_ec);
%         fprintf("\n")
    end
end


% M = repmat(struct([]),T_MAX,1); % doesn't work

% create a custom color map
disp("Creating color map")
my_map = turbo(2^map_res);
if (use_lim)
    step = (upp_lim - low_lim)/length(my_map); % step size for each color (row)
    zero_row = -low_lim / step + 1;
    des_tol = tol_factor * mean([upp_lim abs(low_lim)]);
    tol = ceil(des_tol / step);
    % tol = 0;
    for i = ceil(zero_row - tol):floor(zero_row + tol)
        if (i > 0 && i < length(my_map))
            my_map(i,:) = [1 1 1];
        end
    end
    if (~linear_cmap) % nonlinearize color map
        dist_center = 0; % center point for nonlinear distortion
        var = 1:length(my_map); 
        var = var - (dist_center - low_lim) * length(var) / (upp_lim - low_lim);
        var = nonlin_factor * var / max(abs(var));
        var = sign(var) .* exp(abs(var));
        var = var - min(var);
        var = var * (2^map_res - 1) / max(var) + 1;
        my_map = interp1(var, my_map, 1:(2^map_res));
    end
end

% now open data file
% each line: (Ex, Ey, Ez)
file_id = fopen(input_file);
fgetl(file_id); % discard 1st header line

file_format = '%f %f';
file_size = [num_field_data T_NEXT];
% file_size = [num_field_data Inf];
% disp("Reading in data from file...")
% data = fscanf(file_id, file_format, file_size);
% fclose(file_id);
% disp("Finished reading file data")
% data = transpose(data); % transpose array so orientation is correct
% data = single(data); % change precision to single floating point

output_file = output_file_path + test_name + ofile_mod + ".mp4";
v = VideoWriter(output_file, 'MPEG-4');
v.FrameRate = 8;
open(v);
disp("Opening output mp4 file")


% preallocate arrays/matrices for speed
E_vals_TM = zeros(NUM_X, NUM_Y, 'single');
E_vals_TE = zeros(NUM_X, NUM_Y, 'double');

% arrays to track maxmimum field values
max_vals_TM = zeros(NUM_X, NUM_Y, 'single');
max_vals_TE = zeros(NUM_X, NUM_Y, 'double');
% extract data and generate plots
disp("Beginning plotting animation")
for t = 1:T_MAX
%     fprintf("\nReading in data from file @ t=%d",t)
    data = fscanf(file_id, file_format, file_size);
    data = transpose(data); % transpose array so orientation is correct
    for x = 1:NUM_X
        for y = 1:NUM_Y

%             data_line = fgetl(file_id);
%             E_data = str2double(split(first_line));
%             curr_Ex = E_data(1);
%             curr_Ey = E_data(2);
%             curr_Exy = sqrt(curr_Ex*curr_Ex + curr_Ey*curr_Ey);
%             curr_Ez = E_data(3);


%             curr_row = (double(t)-1)*NUM_X*NUM_Y + (double(x)-1)*NUM_X + (double(y)-1) + 1;
            curr_row = (double(x)-1)*NUM_X + (double(y)-1) + 1;
            curr_Ex = data(curr_row, 1);
            curr_Ey = data(curr_row, 2);
            curr_Exy = sqrt(curr_Ex*curr_Ex + curr_Ey*curr_Ey);
            curr_Ez = data(curr_row, 3);

            if (~max_plot_only)
                E_vals_TM(x,y) = curr_Ez;
                E_vals_TE(x,y) = curr_Exy;
            end

           % test if current value is max value at that point
           % TM
           if (curr_Ez > max_vals_TM(x,y))
               max_vals_TM(x,y) = curr_Ez;
           end
           % TE
           if (curr_Exy > max_vals_TE(x,y))
               max_vals_TE(x,y) = curr_Exy;
           end
            
%             fprintf("(t,x,y,row,Ex,Ey,Ez)=(%d,%d,%d,%d,%f,%f,%f)\n",t-1,x-1,y-1,curr_row, ...
%             curr_Ex,curr_Ey, curr_Ez);
        end
    end

    % check to make sure we don't plot empty arrays (for dielectrics/conductors)
    if (isempty(epsr_xvals))
        plot_epsr = 0; % skip plotting if no dielectrics present
    end
    if (isempty(ec_xvals))
        plot_ec = 0; % skip plotting if no conductivity present
    end
    if (isempty(secondary_xvals))
        plot_secondary = 0;
    end
    
    if (~max_plot_only) % skip this if we are only doing max field plot
        fig1 = figure(1);
    %     tl = tiledlayout(1,num_tiles);
    %     title(tl, "Maximum Field Coverage for Indoor Area")
           
        % plot electric field
        X = DEL_X*(0:1:(NUM_X-1));
        Y = DEL_Y*(0:1:(NUM_Y-1));
    
        s = pcolor(X,Y,transpose(E_vals_TM)); % need to transpose again to reorient (x,y) axes
        
        xlim([X(1) X(end)]);
        ylim([Y(1) Y(end)]);
        colormap(my_map);
        if (use_lim)
            caxis([low_lim upp_lim]);
        end
        c = colorbar;
        c.Label.String = '|E_z| [V/m]';
        c.FontSize = 12;
        s.FaceColor = 'interp';
        s.EdgeColor = 'none';
        
        % plot geometry features
        hold on
        if (plot_epsr) % plot dielectric features
            scttr_X = DEL_X * epsr_xvals;
            scttr_Y = DEL_Y * epsr_yvals;
            scatter(scttr_X, scttr_Y, 5, '*', 'black');
            die_const = data(epsr_xvals(1)*NUM_X + epsr_yvals(1) + 1,5); % reread dielectric value from data
            text(scttr_X(1) - DEL_Y*3,scttr_Y(end) + DEL_X*3, "Dielectric of " + die_const);
        end
        
        if (plot_ec) % plot conductive features
            scttr_X = DEL_X * ec_xvals;
            scttr_Y = DEL_Y * ec_yvals;
            scatter(scttr_X, scttr_Y, 25, 'filled', 'square', ...
            'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[0 0 0]);
        end
    
        if (plot_secondary) % plot additional features
            scttr_X = DEL_X * secondary_xvals;
            scttr_Y = DEL_Y * secondary_yvals;
            scatter(scttr_X, scttr_Y, 25, 'filled', 'square', ...
            'MarkerEdgeColor',[200 120 35]./255, 'MarkerFaceColor',[200 120 35]./255);
        end
    
        hold off
        
        xlabel('x [m]','FontSize',12);
        ylabel('y [m]','FontSize',12);
        title(figure_title);
    
        % determine proper frequency and time units for subtitle
        if (FREQ > 1e9)
            freq_text = " @" + round(FREQ/1e9) + "GHz";
        elseif (FREQ > 1e6)
            freq_text = " @" + round(FREQ/1e6) + "MHz";
        elseif (FREQ > 1e3)
            freq_text = " @" +round(FREQ/1e3) + "kHz";
        else
            freq_text = " @" + FREQ + "Hz ";
        end
        
        if (t * DEL_T < 1e-9)
            time_text = " @t = " + (t * DEL_T * 1e12) + "ps";
        elseif (t * DEL_T < 1e-6)
            time_text = " @t = " + (t * DEL_T * 1e9) + "ns";
        elseif (t * DEL_T < 1e-3)
            time_text = " @t = " + (t * DEL_T * 1e6) + "µs";
        elseif (t * DEL_T < 1)
            time_text = " @t = " + (t * DEL_T * 1e3) + "ms";
        else
            time_text = " @t = " + (t * DEL_T) + "s";
        end
    
        subtitle("Magnitude of E_z(x,y)" + freq_text + time_text);
        M(t) = getframe(1);
        writeVideo(v,M(t));
    %     fprintf("Finished loop t=%d\n", i-1);
    end
fprintf("\nFinished t=%d",t);
end
fclose(file_id);

close(v);
fprintf("MP4 file created");

%% max field plot

% now plot maximum fields (convert to power)
max_vals_TM = abs(max_vals_TM).*abs(max_vals_TM)./2; % magnitude of field squared
max_vals_TM = 10 .* log10(1000.*max_vals_TM); % convert to dBm
max_vals_TE = abs(max_vals_TE).*abs(max_vals_TE)./2; % magnitude of field
max_vals_TE = 10 .* log10(1000.*max_vals_TE); %convert dBm

% create a custom color map for max field plot
my_map = turbo(2^map_res);
if (use_lim2)
    step = (upp_lim2 - low_lim2)/length(my_map); % step size for each color (row)
    zero_row = -low_lim2 / step + 1;
    des_tol = tol_factor * mean([upp_lim2 abs(low_lim2)]);
    tol = ceil(des_tol / step);
    % tol = 0;
    for i = ceil(zero_row - tol):floor(zero_row + tol)
        if (i > 0 && i < length(my_map))
            my_map(i,:) = [1 1 1];
        end
    end
    if (~linear_cmap2) % nonlinearize color map
        dist_center = 0; % center point for nonlinear distortion
        var = 1:length(my_map); 
        var = var - (dist_center - low_lim2) * length(var) / (upp_lim2 - low_lim2);
        var = nonlin_factor * var / max(abs(var));
        var = sign(var) .* exp(abs(var));
        var = var - min(var);
        var = var * (2^map_res - 1) / max(var) + 1;
        my_map = interp1(var, my_map, 1:(2^map_res));
    end
end

figure(2);
tl = tiledlayout(1,num_tiles);
title(tl, "Maximum Field Coverage for Indoor Area")

for i = 1:num_tiles
% plot electric field
    if (i == 1)
        plot_data = max_vals_TM;
    elseif (i == 2)
        plot_data = max_vals_TE;
    else
        disp("ERROR");
    end
    nexttile
    X = DEL_X*(0:1:(NUM_X-1));
    Y = DEL_Y*(0:1:(NUM_Y-1));
    s = pcolor(X,Y,transpose(plot_data)); % need to transpose again to reorient (x,y) axes
    xlim([X(1) X(end)]);
    ylim([Y(1) Y(end)]);
    colormap(my_map);
    if (use_lim2)
        caxis([low_lim2 upp_lim2]);
    end
    c = colorbar;
    c.Label.String = '|Power| [dBm]';
    c.FontSize = 12;
    s.FaceColor = 'interp';
    s.EdgeColor = 'none';
    
    % plot geometry features
    hold on
    if (plot_epsr) % plot dielectric features
        scttr_X = DEL_X * epsr_xvals;
        scttr_Y = DEL_Y * epsr_yvals;
        scatter(scttr_X, scttr_Y, 5, '*', 'black');
        die_const = data(epsr_xvals(1)*NUM_X + epsr_yvals(1) + 1,5); % reread dielectric value from data
        text(scttr_X(1) - DEL_Y*3,scttr_Y(end) + DEL_X*3, "Dielectric of " + die_const);
    end
    
    if (plot_ec) % plot conductive features
        scttr_X = DEL_X * ec_xvals;
        scttr_Y = DEL_Y * ec_yvals;
        scatter(scttr_X, scttr_Y, 25, 'filled', 'square', ...
        'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[0 0 0]);
    end

    if (plot_secondary) % plot additional features
        scttr_X = DEL_X * secondary_xvals;
        scttr_Y = DEL_Y * secondary_yvals;
        scatter(scttr_X, scttr_Y, 25, 'filled', 'square', ...
        'MarkerEdgeColor',[200 120 35]./255, 'MarkerFaceColor',[200 120 35]./255);
    end

    hold off
    
    xlabel('x [m]','FontSize',12);
    ylabel('y [m]','FontSize',12);
    title(figure_title);
    if (i == 1)
        title("TMz Polarization");
    elseif (i ==2)
        title("TEz Polarization");
    end

    % determine proper frequency and time units for subtitle
    if (FREQ > 1e9)
        freq_text = " @" + round(FREQ/1e9) + "GHz";
    elseif (FREQ > 1e6)
        freq_text = " @" + round(FREQ/1e6) + "MHz";
    elseif (FREQ > 1e3)
        freq_text = " @" + round(FREQ/1e3) + "kHz";
    else
        freq_text = " @" + FREQ + "Hz ";
    end
    if (i == 1)
        subtitle("Maximum Value of |E| [TMz]" + freq_text);
    elseif (i ==2)
        subtitle("Maximum Value of |E| [TEz]" + freq_text);
    end
end

% save graph    
pic_file = output_file_path + test_name + ofile_mod + ".png";
exportgraphics(tl, pic_file);




