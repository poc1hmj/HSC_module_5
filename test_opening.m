file_names = {'M013_V02_H06_Achsialhub_C001H001S0001'};

jj=1;
%% Input parameter

sampling_frequency = 100000;
time_step = 1/sampling_frequency;
pump_speed=3500; %rpm
cam_lobes=4;
t_AS=60/pump_speed/cam_lobes; %s
OT=46.2;  %TDC

R_ball_mm = 3.969;
Lift_meas=0.710; % in mm, measured lift
Lift_h0=0.7; % in mm, nominal lift
% Lift_meas=0.610; % in mm, measured lift
% Lift_h0=0.6; % in mm, nominal lift
dLift_h1=0.0462; % in mm, early impact due to geometry at opening
Closing_h0=0;
dClosing_h1=0.0276; % in mm, early impact due to geometry at closing

figure(4500)
hold on, grid on
p_AS_ref=round(t_AS/time_step);

time_ref = 0:time_step:(p_AS_ref-1)*time_step;
for i=1:1
    output_directory = '\\bosch.com\dfsrb\DfsTH\LOC\Hmj\GS_ENGTH\100_Processes_Methods\01_Software\01_MATLAB\10_Project_request\30_HighSpeedCamera\Transfer_HSC\Re_module5\test\';
    
    my_video_name = [output_directory, file_names{i}];
    
    output_table = [my_video_name, '_output_01.xlsx'];
    info_prev = xlsread(output_table,'data_03', 'C10:F30000');
    
    [rows_excel, cols_excel] = size(info_prev);
    
    info_prev_2 = xlsread(output_table,'data_02', 'C12:C13');
    
    % info from the already existing table
    ball_radius_pixel = info_prev_2(1,1);
    pixel_ratio = info_prev_2(2,1);
    
    
    %% Get info from previous module
    
    x_disp_counter_bar = info_prev(:,1);
    y_disp_counter_bar = info_prev(:,2);
    x_disp_ball = info_prev(:,3);
    y_disp_ball = info_prev(:,4);
    
    % conversion to mm
    x_disp_counter_bar = x_disp_counter_bar*pixel_ratio;
    y_disp_counter_bar = y_disp_counter_bar*pixel_ratio;
    x_disp_ball = x_disp_ball*pixel_ratio;
    y_disp_ball = y_disp_ball*pixel_ratio;
    
    time_axis = time_step:time_step:time_step*rows_excel;
    time_axis = time_axis';
    
    % calculate the relative displacements
    x_ball_rel = x_disp_counter_bar - x_disp_ball;
    y_ball_rel = y_disp_ball - y_disp_counter_bar;
    
    %% Calculation
    
    % subtracting offset value
    sorted_x_ball_rel =sort(x_ball_rel);
    offset_val = sorted_x_ball_rel(round(0.01*rows_excel));
    x_ball_rel = x_ball_rel - offset_val;
    
    % subtracting offset value
    sorted_y_ball_rel =sort(y_ball_rel);
    offset_val_y = sorted_y_ball_rel(round(0.01*rows_excel));
    y_ball_rel = y_ball_rel - offset_val_y;
    
    
    % Calculate the velocities
    x_vel_ball = diff(x_ball_rel)/time_step;
    y_vel_ball = diff(x_ball_rel)/time_step;
    
    %% Find opening and closing flanks in hub signal
    
    detection_threshold = 0.13; % micrometer 0.05 by default;
    pos_flank = 0;
    neg_flank = 0;
    positive_flanks = zeros(30,1);
    negative_flanks = zeros(30,1);
    auxiliar_flanks_pos = zeros(rows_excel,1);
    auxiliar_flanks_neg = zeros(rows_excel,1);
    
    for ii=1:1:rows_excel-20
        if pos_flank==0
            if x_ball_rel(ii)<detection_threshold && x_ball_rel(ii+1)>=detection_threshold
                pos_flank = pos_flank + 1;
                positive_flanks(pos_flank) = ii;
                auxiliar_flanks_pos(ii) = detection_threshold;
            end
        else
            if x_ball_rel(ii)<detection_threshold && x_ball_rel(ii+1)>=detection_threshold && ii>(negative_flanks(neg_flank)+50)
                pos_flank = pos_flank + 1;
                positive_flanks(pos_flank) = ii;
                auxiliar_flanks_pos(ii) = detection_threshold;
            end
        end
        % ------------- negative site ------------------
        if neg_flank==0
            if x_ball_rel(ii)>=detection_threshold && x_ball_rel(ii+1)<detection_threshold
                neg_flank = neg_flank + 1;
                negative_flanks(neg_flank) = ii;
                auxiliar_flanks_neg(ii) = detection_threshold;
            end
        else
            if x_ball_rel(ii)>=detection_threshold && x_ball_rel(ii+1)<detection_threshold && ii>(negative_flanks(neg_flank)+50)%neg
                neg_flank = neg_flank + 1;
                negative_flanks(neg_flank) = ii;
                auxiliar_flanks_neg(ii) = detection_threshold;
            end
        end
    end
    
    %% Visualize the ball displacements
    
   
    
    
    %% split each lifts and plot
    
    pointer=find (negative_flanks(1:neg_flank)<positive_flanks(3)); %negative flank of second WC
    pointer=pointer(end);
    point_OT=round(t_AS/time_step/90*OT);
    ref_position_lift=negative_flanks(pointer)+round(t_AS/time_step)-point_OT; % AS END postion
    % ref_position_lift=negative_flanks(pointer)+round((positive_flanks(2)-negative_flanks(pointer))/2);
    % points per AS
    p_AS_ref=round(t_AS/time_step);
    
    time_ref = 0:time_step:(p_AS_ref-1)*time_step;
    
   
    data_hub = zeros(p_AS_ref, neg_flank*2);
    k=0;
    
%     % First lift
%     point_miss = round(2*t_AS/time_step)-ref_position_lift+1; %if fisrt lift is shorter than loop period
%     if point_miss>0
%         x_ball_rel_temp_1 = x_ball_rel(1:p_AS_ref-point_miss); %
%     else
%         x_ball_rel_temp_1 = x_ball_rel(abs(point_miss):p_AS_ref+abs(point_miss)-1);
%     end
%     x_ball_rel_temp = zeros(p_AS_ref,1);
%     x_ball_rel_temp((end-length(x_ball_rel_temp_1)+1):end,1) = x_ball_rel_temp_1; %if not enough point at missing point at start replac with 0
%     plot(time_ref, x_ball_rel_temp)
%     k=k+1;
%     data_hub(1:length(x_ball_rel_temp),k) = time_ref;
%     k=k+1;
%     data_hub(1:length(x_ball_rel_temp),k) = x_ball_rel_temp;
%     clear x_ball_rel_temp point_miss
    x_ball_rel_temp=x_ball_rel(round(ref_position_lift+(jj-3)*(t_AS/time_step)+1):round(ref_position_lift+(jj-2)*(t_AS/time_step))); %get point period of each lift
    x_ball_rel_temp = movmean(x_ball_rel_temp,5);
    x_ball_vel_temp=x_vel_ball(round(ref_position_lift+(jj-3)*(t_AS/time_step)+1):round(ref_position_lift+(jj-2)*(t_AS/time_step))); %get point period of each lift
    x_ball_vel_temp = movmean(x_ball_vel_temp,5);
%     subplot(2,1,1)
    plot(time_ref(1:length(x_ball_rel_temp)), x_ball_rel_temp)
%     hold on
%     subplot(2,1,2)
%     plot(time_ref(1:length(x_ball_vel_temp)), x_ball_vel_temp)
%     hold on
    k=k+1;
    data_hub(1:length(x_ball_rel_temp),k)=time_ref(1:length(x_ball_rel_temp));
    k=k+1;
    data_hub(1:length(x_ball_rel_temp),k)=x_ball_rel_temp;
    clear x_ball_rel_temp
    
end
legend()
