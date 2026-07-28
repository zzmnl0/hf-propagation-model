function a4_c1_same_angle_pharlap(input_mat, output_dir, summary_json, angles)
% Fixed-angle PHaRLAP forward runs for A4-C1-1; no root solving.

pharlap_root = 'D:\pharlap\pharlap_4.7.4';
setenv('DIR_MODELS_REF_DAT', fullfile(pharlap_root, 'dat'));
setenv('NRT_NUM_THREADS', '1');
addpath(fullfile(pharlap_root, 'mex'), '-begin');
addpath(genpath(fullfile(pharlap_root, 'src', 'matlab')), '-end');

case_ids = {'low_pharlap_anchor', 'low_ah3d_anchor', ...
    'high_pharlap_anchor', 'high_ah3d_anchor'};
result = struct('schema_version', 'a4-c1-1-pharlap-fixed-angle-v1', ...
    'status', 'blocked_pharlap_fixed_angle', ...
    'case_semantics', 'same_launch_angle_forward_trace_no_root_solving', ...
    'matlab_version', version, 'matlab_release', version('-release'), ...
    'pharlap_version', '4.7.4', 'raytrace_3d_sp_path', which('raytrace_3d_sp'), ...
    'error_identifier', '', 'error_message', '');

try
    in = load(input_mat);
    assert(strcmp(string(in.schema_version), 'a4-c1-prep1-input-v1') ...
        && strcmp(string(in.grid_variant), 'refined_1km'), ...
        'A4C1C1:WrongInput', 'C1-1 requires frozen refined_1km input');
    assert(numel(angles) == 4 && all(isfinite(angles)), ...
        'A4C1C1:WrongAngles', 'Four finite frozen anchor angles are required');
    assert(endsWith(lower(result.raytrace_3d_sp_path), '.mexw64'), ...
        'A4C1C1:MexNotSelected', 'raytrace_3d_sp MEX is not selected');
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end

    started = tic;
    cases = repmat(struct('case_id', '', 'elevation_deg', NaN, ...
        'raw_mat', '', 'runs', []), 1, 4);
    for case_idx = 1:4
        run_raw = repmat(struct('ray_data', [], 'ray_path_data', [], ...
            'ray_state_vec', []), 1, 2);
        runs = repmat(empty_run(), 1, 2);
        for repeat_idx = 1:2
            [ray_data, ray_path, state] = trace_one(in, angles(case_idx));
            runs(repeat_idx) = summarize(ray_data, ray_path, angles(case_idx), in);
            run_raw(repeat_idx).ray_data = ray_data;
            run_raw(repeat_idx).ray_path_data = ray_path;
            run_raw(repeat_idx).ray_state_vec = state;
        end
        case_result = struct('case_id', case_ids{case_idx}, ...
            'elevation_deg', angles(case_idx), 'runs', runs);
        raw_mat = fullfile(output_dir, ['pharlap_' case_ids{case_idx} '_raw.mat']);
        save(raw_mat, 'run_raw', 'case_result', '-v7.3');
        cases(case_idx) = struct('case_id', case_ids{case_idx}, ...
            'elevation_deg', angles(case_idx), 'raw_mat', raw_mat, 'runs', runs);
    end
    result.wall_time_s = toc(started);
    result.angles_deg = angles;
    result.cases = cases;
    result.checks = struct('four_cases', numel(cases) == 4, ...
        'two_repeats_each', all(arrayfun(@(x) numel(x.runs) == 2, cases)), ...
        'all_ground', all(arrayfun(@(x) all([x.runs.ray_label] == 1), cases)), ...
        'all_single_turning', all(arrayfun(@(x) all([x.runs.turning_count] == 1), cases)), ...
        'raw_outputs_preserved', true, 'no_root_solving', true);
    result.status = 'passed_pharlap_fixed_angle_forward_runs';
catch exception
    result.error_identifier = exception.identifier;
    result.error_message = exception.message;
end

fid = fopen(summary_json, 'w');
assert(fid >= 0, 'A4C1C1:OutputOpenFailed', 'Cannot open summary JSON');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', jsonencode(result, PrettyPrint=true));
clear cleanup raytrace_3d_sp;
if ~strcmp(result.status, 'passed_pharlap_fixed_angle_forward_runs')
    error('A4C1C1:TraceFailed', '%s: %s', ...
        result.error_identifier, result.error_message);
end
end


function [ray_data, ray_path, state] = trace_one(in, elevation)
[ray_data, ray_path, state] = raytrace_3d_sp( ...
    in.origin_lat_deg, in.origin_lon_deg, in.origin_height_km, ...
    elevation, 0.0, in.frequency_MHz, in.OX_mode, in.nhops, in.tol, ...
    in.earth_radius_m, in.iono_en_grid_cm3, in.iono_en_grid_5_cm3, ...
    in.collision_freq_Hz, in.iono_grid_parms, in.Bx_T, in.By_T, in.Bz_T, ...
    in.geomag_grid_parms);
end


function run = summarize(data, path, elevation, in)
label = last_or(data.ray_label, -100);
lat = last_or(data.lat, NaN); lon = last_or(data.lon, NaN);
group_path = last_or(data.group_range, NaN);
run = empty_run();
run.elevation_deg = elevation;
run.ray_label = label;
run.status = label_status(label);
run.ground_range_km = last_or(data.ground_range, NaN);
run.along_track_residual_km = run.ground_range_km - in.target_range_km;
run.cross_track_km = cross_track(lat, lon, in.receiver_lon_deg, ...
    in.earth_radius_m / 1000.0);
run.apex_km = last_or(data.apogee, max(path.height, [], 'omitnan'));
run.geometric_path_km = last_or(data.geometric_path_length, NaN);
run.group_path_km = group_path;
run.tau_ms = group_path / 299792.458 * 1000.0;
run.landing_lat_deg = lat; run.landing_lon_deg = lon;
run.path_point_count = numel(path.lat);
run.turning_count = turning_count(path.height);
run.forward_solver = 'PHaRLAP_raytrace_3d_sp_fixed_angle_no_root_solving';
end


function count = turning_count(height)
delta = diff(double(height(:)));
scale = max([abs(delta); 1.0]);
signs = sign(delta(abs(delta) > 1.0e-9 * scale));
count = sum(signs(1:end-1) > 0 & signs(2:end) < 0);
end


function value = last_or(array, fallback)
if isempty(array), value = fallback; else, value = array(end); end
end


function status = label_status(label)
if label == 1
    status = 'ground';
elseif label == -2
    status = 'penetrated';
else
    status = 'non_ground';
end
end


function value = cross_track(lat, lon, rx_lon, radius)
if ~isfinite(lat) || ~isfinite(lon), value = NaN; return; end
p = radius * [cosd(lat)*cosd(lon), cosd(lat)*sind(lon), sind(lat)];
east = [-sind(rx_lon), cosd(rx_lon), 0.0];
value = dot(p, east);
end


function run = empty_run()
run = struct('elevation_deg', NaN, 'ray_label', -100, 'status', 'non_ground', ...
    'ground_range_km', NaN, 'along_track_residual_km', NaN, ...
    'cross_track_km', NaN, 'apex_km', NaN, 'geometric_path_km', NaN, ...
    'group_path_km', NaN, 'tau_ms', NaN, 'landing_lat_deg', NaN, ...
    'landing_lon_deg', NaN, 'path_point_count', 0, 'turning_count', 0, ...
    'forward_solver', '');
end
