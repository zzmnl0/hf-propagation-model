function a4_c1p_pharlap_scan(input_mat, output_dir, summary_json)
% A4-C1P PHaRLAP fixed-angle numerical-parameter scan; no root solving.

pharlap_root = 'D:\pharlap\pharlap_4.7.4';
addpath(fullfile(pharlap_root, 'mex'), '-begin');
addpath(genpath(fullfile(pharlap_root, 'src', 'matlab')), '-end');

profile_ids = {'baseline', 'loose_tol', 'tight_tol', 'ultratight_tol', ...
    'max_step_5', 'max_step_2', 'max_step_1', 'min_step_0p1', ...
    'min_step_0p01', 'reference'};
tol_matrix = [ ...
    1e-8,  0.01,  10; ...
    1e-6,  0.01,  10; ...
    1e-10, 0.01,  10; ...
    1e-12, 0.01,  10; ...
    1e-10, 0.01,   5; ...
    1e-10, 0.01,   2; ...
    1e-10, 0.01,   1; ...
    1e-12, 0.1,    1; ...
    1e-12, 0.01,   1; ...
    1e-12, 0.001,  1];
case_ids = {'low_pharlap_anchor', 'low_ah3d_anchor', ...
    'high_pharlap_anchor', 'high_ah3d_anchor'};
branches = {'low_ray', 'low_ray', 'high_ray', 'high_ray'};
angles = [21.94290326504532, 21.94562872092976, ...
    39.54254824570964, 39.58701138438118];

result = struct('schema_version', 'a4-c1p-pharlap-scan-v1', ...
    'status', 'blocked_c1p1_pharlap_scan', ...
    'case_semantics', 'fixed_angles_parameter_scan_no_root_solving', ...
    'matlab_version', version, 'matlab_release', version('-release'), ...
    'pharlap_version', '4.7.4', 'raytrace_3d_sp_path', which('raytrace_3d_sp'), ...
    'profile_ids', {profile_ids}, 'tol_matrix', tol_matrix, ...
    'case_ids', {case_ids}, 'branches', {branches}, 'angles_deg', angles, ...
    'records', [], 'error_identifier', '', 'error_message', '');

try
    in = load(input_mat);
    assert(strcmp(string(in.schema_version), 'a4-c1-prep1-input-v1') ...
        && strcmp(string(in.grid_variant), 'refined_1km'), ...
        'A4C1P:WrongInput', 'C1P requires frozen refined_1km input');
    assert(in.frequency_MHz == 10 && in.OX_mode == 0 && in.nhops == 1, ...
        'A4C1P:WrongPhysics', 'C1P requires 10 MHz, B=0/no-field, one hop');
    assert(all(in.collision_freq_Hz(:) == 0) && in.earth_radius_m == 6371000, ...
        'A4C1P:WrongPhysics', 'C1P requires collision=0 and R=6371000 m');
    assert(endsWith(lower(result.raytrace_3d_sp_path), '.mexw64'), ...
        'A4C1P:MexNotSelected', 'raytrace_3d_sp MEX is not selected');
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end

    records = repmat(empty_record(), 1, numel(profile_ids) * numel(case_ids));
    record_idx = 0;
    started = tic;
    for profile_idx = 1:numel(profile_ids)
        repeat_count = 1 + double(profile_idx == 1 || profile_idx == numel(profile_ids));
        for case_idx = 1:numel(case_ids)
            record_idx = record_idx + 1;
            run_raw = repmat(struct('ray_data', [], 'ray_path_data', [], ...
                'ray_state_vec', []), 1, repeat_count);
            runs = repmat(empty_run(), 1, repeat_count);
            for repeat_idx = 1:repeat_count
                [ray_data, ray_path, state] = trace_one( ...
                    in, angles(case_idx), tol_matrix(profile_idx, :));
                runs(repeat_idx) = summarize( ...
                    ray_data, ray_path, angles(case_idx), in);
                run_raw(repeat_idx).ray_data = ray_data;
                run_raw(repeat_idx).ray_path_data = ray_path;
                run_raw(repeat_idx).ray_state_vec = state;
            end
            profile_id = profile_ids{profile_idx}; %#ok<NASGU>
            case_id = case_ids{case_idx}; %#ok<NASGU>
            branch_id = branches{case_idx}; %#ok<NASGU>
            elevation_deg = angles(case_idx); %#ok<NASGU>
            tol = tol_matrix(profile_idx, :); %#ok<NASGU>
            raw_mat = fullfile(output_dir, ...
                ['pharlap_' profile_id '_' case_id '_raw.mat']);
            save(raw_mat, 'run_raw', 'profile_id', 'case_id', 'branch_id', ...
                'elevation_deg', 'tol', '-v7.3');
            records(record_idx) = struct( ...
                'profile_id', profile_id, 'tol', tol, 'case_id', case_id, ...
                'branch_id', branch_id, 'elevation_deg', elevation_deg, ...
                'repeat_count', repeat_count, 'raw_mat', raw_mat, 'runs', runs);
        end
    end
    result.records = records;
    result.wall_time_s = toc(started);
    result.status = 'completed_c1p1_pharlap_fixed_angle_scan';
    result.checks = struct('no_root_solving', true, ...
        'ten_profiles', numel(profile_ids) == 10, ...
        'four_frozen_angles', numel(case_ids) == 4, ...
        'baseline_and_reference_repeated', true);
catch ME
    result.error_identifier = ME.identifier;
    result.error_message = ME.message;
end

fid = fopen(summary_json, 'w');
assert(fid >= 0, 'A4C1P:SummaryOpen', 'Cannot open summary JSON');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', jsonencode(result, PrettyPrint=true));
clear cleanup raytrace_3d_sp;
if ~strcmp(result.status, 'completed_c1p1_pharlap_fixed_angle_scan')
    error('A4C1P:ScanFailed', '%s: %s', ...
        result.error_identifier, result.error_message);
end
end


function [ray_data, ray_path, state] = trace_one(in, elevation, tol)
[ray_data, ray_path, state] = raytrace_3d_sp( ...
    in.origin_lat_deg, in.origin_lon_deg, in.origin_height_km, ...
    elevation, 0.0, in.frequency_MHz, in.OX_mode, in.nhops, tol, ...
    in.earth_radius_m, in.iono_en_grid_cm3, in.iono_en_grid_5_cm3, ...
    in.collision_freq_Hz, in.iono_grid_parms, in.Bx_T, in.By_T, in.Bz_T, ...
    in.geomag_grid_parms);
end


function run = summarize(data, path, elevation, in)
label = last_or(data.ray_label, -100);
lat = last_or(data.lat, NaN);
lon = last_or(data.lon, NaN);
group_path = last_or(data.group_range, NaN);
run = empty_run();
run.elevation_deg = elevation;
run.ray_label = label;
run.status = label_status(label);
run.ground_range_km = great_circle( ...
    in.origin_lat_deg, in.origin_lon_deg, lat, lon, in.earth_radius_m / 1000);
run.along_track_residual_km = run.ground_range_km - in.target_range_km;
run.apex_km = max(path.height);
run.geometric_path_km = last_or(path.geometric_distance, NaN);
run.group_path_km = group_path;
run.tau_ms = 1000 * group_path / 299792.458;
run.landing_lat_deg = lat;
run.landing_lon_deg = lon;
run.path_point_count = numel(path.lat);
run.turning_count = turning_count(path.height);
run.forward_solver = 'PHaRLAP_raytrace_3d_sp_fixed_angle_parameter_scan';
end


function value = last_or(values, fallback)
if isempty(values), value = fallback; else, value = values(end); end
end


function count = turning_count(height)
height = double(height(:));
if numel(height) < 3, count = 0; return; end
difference = diff(height);
scale = max(max(abs(height)), 1);
signs = sign(difference);
signs(abs(difference) <= 1e-10 * scale) = 0;
signs = signs(signs ~= 0);
count = sum(signs(1:end-1) > 0 & signs(2:end) < 0);
end


function distance = great_circle(lat1, lon1, lat2, lon2, radius)
lat1 = deg2rad(lat1); lon1 = deg2rad(lon1);
lat2 = deg2rad(lat2); lon2 = deg2rad(lon2);
dot = sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon2 - lon1);
distance = radius * acos(max(-1, min(1, dot)));
end


function status = label_status(label)
if label == 1
    status = 'ground';
elseif label == -2
    status = 'penetrated';
elseif label == -3
    status = 'grid_exit';
elseif label == -4
    status = 'max_points';
elseif label == -100
    status = 'catastrophic_error';
else
    status = 'other';
end
end


function run = empty_run()
run = struct('elevation_deg', NaN, 'ray_label', -100, ...
    'status', 'other', 'ground_range_km', NaN, ...
    'along_track_residual_km', NaN, 'apex_km', NaN, ...
    'geometric_path_km', NaN, 'group_path_km', NaN, 'tau_ms', NaN, ...
    'landing_lat_deg', NaN, 'landing_lon_deg', NaN, ...
    'path_point_count', 0, 'turning_count', 0, 'forward_solver', '');
end


function record = empty_record()
record = struct('profile_id', '', 'tol', [NaN, NaN, NaN], ...
    'case_id', '', 'branch_id', '', 'elevation_deg', NaN, ...
    'repeat_count', 0, 'raw_mat', '', 'runs', empty_run());
end
