clc
clear

%% Simulation
test_name = "readtest";
input_file = "/Users/williampoland/FDTD/2D/raw/" + test_name + ".txt";
% each line: (t, x, y, z, Ez, Hx, Hy, eps_r, mu_r, e_cond, m_cond)
file_id = fopen(input_file);
fgetl(file_id); % discard 1st header line (FREQ TMAX NX NY DT DX DY)
% now extract values
first_line = fgetl(file_id);
header_info = str2double(split(first_line));
% ii = 0;
% header_info = zeros(1,7);
% while (length(first_line) > 1)
%     header_info(1,ii) = split(
% end
FREQ = header_info(1);
% T_MAX = str2double(fgetl(file_id));
% NUM_X = str2double(fgetl(file_id));
% NUM_Y = str2double(fgetl(file_id));
% DEL_T = 0;
% DEL_X = 0;
% DEL_Y = 0;
% T_NEXT = NUM_X * NUM_Y; % number of lines before advancing to next time index
fgetl(file_id); % discard 2nd header line (t, x, y, Ez, ec)

file_format = '%f %f';
file_size = [5 Inf];
data = fscanf(file_id, file_format, file_size);