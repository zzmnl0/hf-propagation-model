function a4_c1_prep5_pharlap(input_mat, raw_mat, summary_json, brackets)
% PHaRLAP 1-km local branch re-solve for A4-C1 Prep-5 convergence only.

pharlap_root = 'D:\pharlap\pharlap_4.7.4';
setenv('DIR_MODELS_REF_DAT', fullfile(pharlap_root, 'dat'));
setenv('NRT_NUM_THREADS', '1');
addpath(fullfile(pharlap_root, 'mex'), '-begin');
addpath(genpath(fullfile(pharlap_root, 'src', 'matlab')), '-end');

result = struct('schema_version', 'a4-c1-prep5-pharlap-v1', ...
    'status', 'blocked_pharlap_refined_roots', ...
    'case_semantics', 'pharlap_1km_branch_stability_not_a4_c1_validation', ...
    'matlab_version', version, 'matlab_release', version('-release'), ...
    'pharlap_version', '4.7.4', 'raytrace_3d_sp_path', which('raytrace_3d_sp'), ...
    'error_identifier', '', 'error_message', '');

try
    in = load(input_mat);
    assert(strcmp(string(in.schema_version), 'a4-c1-prep1-input-v1') ...
        && strcmp(string(in.grid_variant), 'refined_1km'), ...
        'A4C1P5:WrongInput', 'Prep-5 requires frozen refined_1km input');
    assert(isequal(size(brackets), [2 2]) && all(isfinite(brackets), 'all') ...
        && all(brackets(:, 1) < brackets(:, 2)), ...
        'A4C1P5:WrongBrackets', 'Two finite ordered branch brackets are required');
    assert(endsWith(lower(result.raytrace_3d_sp_path), '.mexw64'), ...
        'A4C1P5:MexNotSelected', 'raytrace_3d_sp MEX is not selected');

    started = tic;
    roots = repmat(empty_root(), 1, 2);
    root_raw = repmat(struct('ray_data', [], 'ray_path_data', [], ...
        'ray_state_vec', []), 1, 2);
    options = optimset('TolX', 1e-8, 'Display', 'off');
    for idx = 1:2
        bracket = brackets(idx, :);
        [elevation, ~, exitflag] = fzero(@root_residual, bracket, options);
        [ray_data, ray_path, state] = trace_one(in, elevation);
        roots(idx) = summarize(ray_data, ray_path, elevation, bracket, exitflag, in);
        root_raw(idx).ray_data = ray_data;
        root_raw(idx).ray_path_data = ray_path;
        root_raw(idx).ray_state_vec = state;
    end
    result.wall_time_s = toc(started);
    result.brackets_deg = brackets;
    result.roots = roots;
    result.checks = struct('two_ground_roots', all([roots.ray_label] == 1), ...
        'roots_inside_brackets', all([roots.elevation_deg]' >= brackets(:, 1) ...
        & [roots.elevation_deg]' <= brackets(:, 2)), ...
        'raw_outputs_preserved', true);
    result.status = 'passed_pharlap_1km_branch_resolve';
    save(raw_mat, 'root_raw', 'result', '-v7.3');
catch exception
    result.error_identifier = exception.identifier;
    result.error_message = exception.message;
end

fid = fopen(summary_json, 'w');
assert(fid >= 0, 'A4C1P5:OutputOpenFailed', 'Cannot open summary JSON');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', jsonencode(result, PrettyPrint=true));
clear cleanup raytrace_3d_sp;
if ~strcmp(result.status, 'passed_pharlap_1km_branch_resolve')
    error('A4C1P5:ResolveFailed', '%s: %s', ...
        result.error_identifier, result.error_message);
end

    function residual = root_residual(elevation)
        [ray_data, ray_path] = trace_one(in, elevation);
        record = summarize(ray_data, ray_path, elevation, [NaN NaN], 0, in);
        assert(record.ray_label == 1, 'A4C1P5:NonGroundRootIteration', ...
            'Bounded root iteration left the continuous ground branch');
        residual = record.along_track_residual_km;
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


function root = summarize(data, path, elevation, bracket, exitflag, in)
label = last_or(data.ray_label, -100);
lat = last_or(data.lat, NaN); lon = last_or(data.lon, NaN);
group_path = last_or(data.group_range, NaN);
root = empty_root();
root.elevation_deg = elevation;
root.ray_label = label;
root.status = label_status(label);
root.ground_range_km = last_or(data.ground_range, NaN);
root.along_track_residual_km = root.ground_range_km - in.target_range_km;
root.cross_track_km = cross_track(lat, lon, in.receiver_lon_deg, ...
    in.earth_radius_m / 1000.0);
root.apex_km = last_or(data.apogee, max(path.height, [], 'omitnan'));
root.geometric_path_km = last_or(data.geometric_path_length, NaN);
root.group_path_km = group_path;
root.tau_ms = group_path / 299792.458 * 1000.0;
root.landing_lat_deg = lat; root.landing_lon_deg = lon;
root.path_point_count = numel(path.lat);
root.bracket_low_deg = bracket(1); root.bracket_high_deg = bracket(2);
root.exitflag = exitflag;
root.root_solver = 'MATLAB_fzero_bracketed_Brent_Dekker_type_refined_1km';
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


function root = empty_root()
root = struct('elevation_deg', NaN, 'ray_label', -100, 'status', 'non_ground', ...
    'ground_range_km', NaN, 'along_track_residual_km', NaN, ...
    'cross_track_km', NaN, 'apex_km', NaN, 'geometric_path_km', NaN, ...
    'group_path_km', NaN, 'tau_ms', NaN, 'landing_lat_deg', NaN, ...
    'landing_lon_deg', NaN, 'path_point_count', 0, 'bracket_low_deg', NaN, ...
    'bracket_high_deg', NaN, 'exitflag', 0, 'root_solver', '');
end
