%% FDTD 2D Simulation Animation Generator
%  William Poland
clc
clear

test_name = "basic_sine1Hz"; % name generated by C++ code

ofile_mod = "_lim"; % modifier for output MP4 file
% ofile_mod = "_lm"; % modifier for output MP4 file
max_lim = 0.35; % maximum amplitude (range of color bar values)
use_lim = 1; % turn caxis limits on/off

input_file = "/Users/williampoland/FDTD/2D/raw/" + test_name + ".txt";
% each line: (t, x, y, z, Ez, Hx, Hy, eps_r, mu_r, e_cond, m_cond)
file_id = fopen(input_file);
fgetl(file_id); % discard 1st header line
% now extract TMAX, NX, NY values
T_MAX = str2double(fgetl(file_id));
NUM_X = str2double(fgetl(file_id));
NUM_Y = str2double(fgetl(file_id));
T_NEXT = NUM_X * NUM_Y; % number of lines before advancing to next time index
fgetl(file_id); % discard 2nd header line

file_format = '%f %f';
file_size = [4 Inf];
data = fscanf(file_id, file_format, file_size);
data = data'; % transpose array so orientation is correct

output_file = "/Users/williampoland/FDTD/2D/" + test_name + ofile_mod + ".mp4";
v = VideoWriter(output_file, 'MPEG-4');
v.FrameRate = 8;
open(v);

% T_MAX = 20; %reduce this for debugging purposes
% M = repmat(struct([]),T_MAX,1); % doesn't work
for t = 1:T_MAX
    vals = zeros(NUM_X, NUM_Y);
    for x = 1:NUM_X
        for y = 1:NUM_Y
            curr_row = (double(t)-1)*NUM_X*NUM_Y + (double(x)-1)*NUM_X + (double(y)-1) + 1;
            curr_val = data(curr_row, 4);
            vals(x,y) = curr_val;
            fprintf("(t,x,y,row, Ez)=(%d,%d,%d,%d,%f)\n",t-1,x-1,y-1,curr_row,curr_val);
            if (curr_val ~= 0)
                fprintf("******val not 0 @ (t,x,y)=(%d,%d,%d)******\n",t-1,x-1,y-1);
            end
        end
    end
    figure(1);
    s = pcolor(vals);
    colormap(turbo(2*4096));
    if (use_lim)
        caxis([-max_lim max_lim]);
    end
    c = colorbar;
    c.Label.String = '|Ez|';
    s.FaceColor = 'interp';
    s.EdgeColor = 'none';
    xlabel('x');
    ylabel('y');
    title("Magnitude of Ez(x,y) @ t = " + t);
    M(t) = getframe(1);
    writeVideo(v,M(t));
%     fprintf("Finished loop t=%d\n", i-1);

end

close(v);

% figure(2)
% movie(M,1,80)

