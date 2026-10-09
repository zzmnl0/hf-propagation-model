function a4_c1m_uniform_export(input_mat, output_dir, summary_json)
% C1M A1: explicit uniform-medium launch; immutable raw outputs, no P2P.
pharlap_root = 'D:\pharlap\pharlap_4.7.4';
addpath(fullfile(pharlap_root, 'mex'), '-begin');
addpath(genpath(fullfile(pharlap_root, 'src', 'matlab')), '-end');
header = load(input_mat);
if isfield(header, 'a2_jobs')
    a2_export(header, output_dir, summary_json);
    return
end
assert(~isfile(summary_json), 'C1M:Exists', 'Summary already exists');
if ~isfolder(output_dir), mkdir(output_dir); end
result = struct('schema_version', 'c1m-uniform-v1', 'status', 'running', ...
    'matlab_version', version, 'mex_path', which('raytrace_3d_sp'), ...
    'records', [], 'cache_records', [], 'error_identifier', '', 'error_message', '');
try
    in = load(input_mat);
    assert(endsWith(lower(result.mex_path), '.mexw64'), 'C1M:Mex', 'MEX not selected');
    assert(in.earth_radius_m == 6371000 && in.origin_height_km == 100, ...
        'C1M:Geometry', 'Unexpected launch geometry');
    assert(all(in.Bx_T(:) == 0) && all(in.By_T(:) == 0) && ...
        all(in.Bz_T(:) == 0) && all(in.collision_freq_Hz(:) == 0), ...
        'C1M:Physics', 'A1 requires B=0 and collision=0');
    assert(isa(in.Ne_m3, 'double') && ndims(in.Ne_m3) == 4, ...
        'C1M:Grid', 'Expected double latitude,longitude,height,density grid');
    assert(size(in.Ne_m3, 1) == numel(in.latitude_deg) && ...
        size(in.Ne_m3, 2) == numel(in.longitude_deg) && ...
        size(in.Ne_m3, 3) == numel(in.height_km), 'C1M:Axes', 'Axis mismatch');
    readback = fullfile(output_dir, 'input_readback.mat');
    assert(~isfile(readback), 'C1M:Exists', 'Readback already exists');
    save(readback, '-struct', 'in', '-v7');
    grid = grid_parameters(in.latitude_deg, in.longitude_deg, in.height_km);

    % Negative control: a jump into a constant slab does not provide a
    % resolved density gradient and cannot calibrate the dispersion constant.
    probe_grid = grid; probe_grid(7) = 60;
    probe_ne = ones(size(in.Bx_T)) * 1e5;  % electrons / cm^3
    clear raytrace_3d_sp;
    [ray_data, ray_path_data, ray_state_vec] = raytrace_3d_sp( ...
        in.origin_lat_deg, in.origin_lon_deg, 0, 22, 0, 10, 0, 1, ...
        in.tol_matrix(end,:), in.earth_radius_m, probe_ne, probe_ne, ...
        in.collision_freq_Hz, probe_grid, in.Bx_T, in.By_T, in.Bz_T, probe_grid);
    raw_mat = fullfile(output_dir, 'calibration_raw.mat');
    assert(~isfile(raw_mat), 'C1M:Exists', 'Calibration already exists');
    save(raw_mat, 'ray_data', 'ray_path_data', 'ray_state_vec', 'probe_grid', '-v7');
    density = ray_path_data.electron_density(:) * 1e6;
    phase = ray_path_data.refractive_index(:);
    valid = ray_path_data.height(:) >= 80 & abs(density / 1e11 - 1) < 1e-8 ...
        & isfinite(phase);
    assert(nnz(valid) >= 10, 'C1M:Calibration', 'Insufficient uniform probe samples');
    coefficients = (1 - phase(valid).^2) * 1e14 ./ density(valid);
    result.calibration_k2 = NaN;
    result.coefficient_identifiable = false;
    result.slab_state_derived_k2 = median(coefficients);
    result.launch_k2 = 80.6404; % Same controlled physical initial state as AH3D.
    result.launch_k2_semantics = 'controlled_common_initial_state_not_measured_PH_constant';
    result.calibration_k2_spread = max(abs(coefficients-median(coefficients)));
    result.calibration_sample_count = nnz(valid);
    fprintf('C1M uniform coefficient unidentifiable; common launch K2=%.15g\n', ...
        result.launch_k2);

    initialization_records = [];
    control_k2 = [0, result.launch_k2, 80.6163849431291];
    for ki = 1:numel(control_k2)
        for mode = [0, 1, -1]
            record = run_one(in, grid, in.Ne_m3(:,:,:,3)/1e6, 3, 2, 1, ...
                size(in.tol_matrix,1), control_k2(ki), output_dir, ...
                sprintf('init_n%d_mode%d',ki,mode), mode);
            initialization_records = [initialization_records; record]; %#ok<AGROW>
        end
    end
    result.initialization_records = initialization_records;

    records = [];
    for ni = 1:numel(in.densities_m3)
        ne = in.Ne_m3(:,:,:,ni) / 1e6;
        assert(all(ne(:) == in.densities_m3(ni)/1e6), ...
            'C1M:Uniform', 'Nonuniform A1 grid');
        for fi = 1:numel(in.frequencies_MHz)
            for ei = 1:numel(in.elevations_deg)
                for li = 1:size(in.tol_matrix,1)
                    case_id = sprintf('case_n%d_f%d_e%d_l%d', ni, fi, ei, li);
                    record = run_one(in, grid, ne, ni, fi, ei, li, ...
                        result.launch_k2, output_dir, case_id);
                    records = [records; record]; %#ok<AGROW>
                end
            end
        end
        fprintf('C1M PHaRLAP density %d/%d complete\n', ni, numel(in.densities_m3));
    end
    result.records = records;
    sequence = [2, 4, 2, 4]; cache_records = [];
    for ci = 1:numel(sequence)
        ni = sequence(ci);
        record = run_one(in, grid, in.Ne_m3(:,:,:,ni)/1e6, ni, 2, 1, ...
            size(in.tol_matrix,1), result.launch_k2, output_dir, ...
            sprintf('cache_%d',ci));
        cache_records = [cache_records; record]; %#ok<AGROW>
    end
    result.cache_records = cache_records;
    result.status = 'completed_uniform_export';
catch ME
    result.status = 'failed_uniform_export';
    result.error_identifier = ME.identifier;
    result.error_message = ME.message;
end

fid = fopen(summary_json, 'w');
assert(fid >= 0, 'C1M:Summary', 'Cannot open summary');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', jsonencode(result, PrettyPrint=true));
clear cleanup raytrace_3d_sp;
if ~strcmp(result.status, 'completed_uniform_export')
    error('C1M:ExportFailed', '%s: %s', result.error_identifier, result.error_message);
end
end

function grid = grid_parameters(lat, lon, height)
grid = [lat(1),lat(2)-lat(1),numel(lat),lon(1),lon(2)-lon(1), ...
    numel(lon),height(1),height(2)-height(1),numel(height)];
end

function record = run_one(in, grid, ne, ni, fi, ei, li, k2, output_dir, case_id, mode)
if nargin < 11, mode = 0; end
frequency = in.frequencies_MHz(fi);
index = sqrt(1-k2*in.densities_m3(ni)/(frequency*1e6)^2);
direction = in.launch_unit_ecef(ei,:) * index;
position = in.origin_ecef_m(:);
initial = struct('pos_x',position(1),'pos_y',position(2),'pos_z',position(3), ...
    'dir_x',direction(1),'dir_y',direction(2),'dir_z',direction(3), ...
    'group_path',0,'geom_path',0,'phase_path',0,'absorption',0, ...
    'indep_var',0,'ODE_step_size',in.tol_matrix(li,2)*1000);
raw_mat = fullfile(output_dir, [case_id '.mat']);
assert(~isfile(raw_mat), 'C1M:Exists', 'Raw case already exists');
% Full grids are supplied even after clear; never use the cached-grid API.
clear raytrace_3d_sp;
[ray_data, ray_path_data, ray_state_vec] = raytrace_3d_sp( ...
    in.origin_lat_deg,in.origin_lon_deg,in.origin_height_km, ...
    in.elevations_deg(ei),0,frequency,mode,1,in.tol_matrix(li,:), ...
    in.earth_radius_m,ne,ne,in.collision_freq_Hz,grid, ...
    in.Bx_T,in.By_T,in.Bz_T,grid,initial);
save(raw_mat,'ray_data','ray_path_data','ray_state_vec','initial','k2','-v7');
record = struct('case_id',case_id,'density_m3',in.densities_m3(ni), ...
    'frequency_MHz',frequency,'elevation_deg',in.elevations_deg(ei), ...
    'level',li,'raw_mat',raw_mat,'ray_label',ray_data.ray_label(end), ...
    'point_count',numel(ray_path_data.lat),'initial_index',index,'OX_mode',mode);
end

function a2_export(in, output_dir, summary_json)
% Independent immutable jobs; summaries are restartable checkpoints.
if ~isfolder(output_dir), mkdir(output_dir); end
jobs = in.a2_jobs;
result = struct('schema_version','c1m-a2-export-v1','status','running', ...
    'matlab_version',version,'mex_path',which('raytrace_3d_sp'),'records',{{}});
assert(endsWith(lower(result.mex_path),'.mexw64'),'C1MA2:Mex','MEX not selected');
failures = 0;
for ji = 1:numel(jobs)
    if iscell(jobs), job = jobs{ji}; else, job = jobs(ji); end
    record = struct('job_id',job.job_id,'stage',job.stage,'raw_mat','', ...
        'job_sha256',job.job_sha256,'grid_sha256',job.grid_sha256, ...
        'ray_label',NaN,'point_count',0,'elapsed_s',0,'reused',false, ...
        'error_identifier','','error_message','');
    started = tic;
    try
        assert(strcmp(file_sha256(job.grid_mat),job.grid_sha256), ...
            'C1MA2:GridHash','Input grid changed');
        raw_mat = fullfile(output_dir,[job.job_id '.mat']);
        exporter_sha256 = file_sha256([mfilename('fullpath') '.m']);
        if isfield(job,'exporter_sha256')
            assert(strcmp(job.exporter_sha256,exporter_sha256), ...
                'C1MP1:ExportSource','Frozen exporter changed');
        end
        if isfile(raw_mat)
            previous = load(raw_mat,'job_sha256','grid_sha256','ray_data','ray_path_data');
            assert(strcmp(previous.job_sha256,job.job_sha256) && ...
                strcmp(previous.grid_sha256,job.grid_sha256), ...
                'C1MA2:ResumeHash','Existing raw job identity differs');
            if isfield(job,'p1_case')
                identity = load(raw_mat,'exporter_sha256','mex_sha256');
                assert(isfield(identity,'exporter_sha256') && ...
                    strcmp(identity.exporter_sha256,exporter_sha256) && ...
                    strcmp(identity.mex_sha256,file_sha256(which('raytrace_3d_sp'))), ...
                    'C1MP1:ResumeSource','P1 execution identity differs');
            end
            ray_data = previous.ray_data; ray_path_data = previous.ray_path_data;
            record.reused = true;
        else
            medium = load(job.grid_mat);
            lat = medium.latitude_deg; lon = medium.longitude_deg; h = medium.height_km;
            assert(numel(h)<=401 && all(diff(h)>0),'C1MA2:Height','Invalid height grid');
            iono_grid = grid_parameters(lat,lon,h);
            bh = h;
            if isfield(medium,'geomag_height_km'), bh = medium.geomag_height_km; end
            assert(numel(bh)<=201,'C1MA2:BHeight','Too many magnetic heights');
            geomag_grid = grid_parameters(lat,lon,bh);
            ne_cm3 = double(medium.Ne_m3)/1e6;
            assert(all(ne_cm3(:)>=0) && all(isfinite(ne_cm3(:))), ...
                'C1MA2:Ne','Invalid density');
            assert(all(medium.Bx_T(:)==0) && all(medium.By_T(:)==0) && ...
                all(medium.Bz_T(:)==0) && all(medium.collision_freq_Hz(:)==0), ...
                'C1MA2:Physics','A2 must be collisionless and unmagnetized');
            readback = fullfile(output_dir,['grid_readback_' job.grid_sha256(1:16) '.mat']);
            if ~isfile(readback)
                save(readback,'-struct','medium','-v7');
            elseif isfield(job,'p1_case')
                assert(isequaln(load(readback),medium), ...
                    'C1MP1:Readback','Existing readback differs');
            end
            args = {30,120,job.origin_height_km,job.elevation_deg,0,job.frequency_MHz, ...
                0,1,job.tol,6371000,ne_cm3,ne_cm3,medium.collision_freq_Hz,iono_grid, ...
                medium.Bx_T,medium.By_T,medium.Bz_T,geomag_grid};
            initial = [];
            if isfield(job,'initial_xyz_km')
                args{1} = job.origin_lat_deg; args{2} = job.origin_lon_deg;
                position = job.initial_xyz_km(:)*1000;
                direction = job.initial_p(:);
                initial = struct('pos_x',position(1),'pos_y',position(2),'pos_z',position(3), ...
                    'dir_x',direction(1),'dir_y',direction(2),'dir_z',direction(3), ...
                    'group_path',0,'geom_path',0,'phase_path',0,'absorption',0, ...
                    'indep_var',0,'ODE_step_size',job.tol(2)*1000);
                args{end+1} = initial;
            elseif job.explicit_state
                ne0 = interp1(h,medium.ne_profile_m3,job.origin_height_km,'linear');
                n2 = 1-job.launch_k2*ne0/(job.frequency_MHz*1e6)^2;
                assert(isfinite(n2) && n2>0,'C1MA2:Launch','Invalid initial dispersion');
                up = [cosd(30)*cosd(120),cosd(30)*sind(120),sind(30)];
                north = [-sind(30)*cosd(120),-sind(30)*sind(120),cosd(30)];
                position = (6371000+job.origin_height_km*1000)*up;
                direction = sqrt(n2)*(cosd(job.elevation_deg)*north+sind(job.elevation_deg)*up);
                initial = struct('pos_x',position(1),'pos_y',position(2),'pos_z',position(3), ...
                    'dir_x',direction(1),'dir_y',direction(2),'dir_z',direction(3), ...
                    'group_path',0,'geom_path',0,'phase_path',0,'absorption',0, ...
                    'indep_var',0,'ODE_step_size',job.tol(2)*1000);
                args{end+1} = initial;
            end
            clear raytrace_3d_sp;
            [ray_data,ray_path_data,ray_state_vec] = raytrace_3d_sp(args{:});
            job_sha256 = job.job_sha256; grid_sha256 = job.grid_sha256;
            mex_sha256 = file_sha256(which('raytrace_3d_sp'));
            save(raw_mat,'ray_data','ray_path_data','ray_state_vec','initial', ...
                'job','job_sha256','grid_sha256','exporter_sha256','mex_sha256','-v7');
        end
        record.raw_mat = raw_mat;
        record.ray_label = ray_data.ray_label(end);
        record.point_count = numel(ray_path_data.lat);
    catch ME
        failures = failures+1;
        record.error_identifier = ME.identifier;
        record.error_message = ME.message;
    end
    record.elapsed_s = toc(started);
    result.records{end+1} = record;
    write_checkpoint(summary_json,result);
    fprintf('C1M A2 PH %d/%d %s label=%g error=%s\n',ji,numel(jobs), ...
        job.job_id,record.ray_label,record.error_identifier);
end
result.failures = failures;
if failures==0, result.status='completed_a2_export'; else, result.status='partial_a2_export'; end
write_checkpoint(summary_json,result);
clear raytrace_3d_sp;
if failures>0, error('C1MA2:JobsFailed','%d jobs failed; completed raw jobs preserved',failures); end
end

function digest = file_sha256(path)
fid = fopen(path,'rb'); assert(fid>=0,'C1MA2:Read','Cannot read input');
cleanup = onCleanup(@() fclose(fid)); bytes = fread(fid,Inf,'*uint8');
md = java.security.MessageDigest.getInstance('SHA-256'); md.update(bytes);
digest = lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end

function write_checkpoint(path,value)
temporary = [path '.tmp'];
fid = fopen(temporary,'w'); assert(fid>=0,'C1MA2:Write','Cannot write summary');
cleanup = onCleanup(@() fclose(fid)); fprintf(fid,'%s',jsonencode(value,PrettyPrint=true));
clear cleanup;
movefile(temporary,path,'f');
end
