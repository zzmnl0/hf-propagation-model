function a4_c1m_p2_exit_export(input_mat, output_dir, summary_json)
% Frozen explicit-state exit controls. The MEX and every raw array are intact.
assert(~isfolder(output_dir),'P2Exit:Exists','Use a fresh raw output directory');
mkdir(output_dir);
in = load(input_mat);
jobs = in.exit_jobs;
result = struct('schema','c1m-p2-exit-export-v1','status','running', ...
    'matlab_version',version,'records',{{}},'failures',0);
for ji = 1:numel(jobs)
    if iscell(jobs), job = jobs{ji}; else, job = jobs(ji); end
    record = struct('job_id',job.job_id,'raw_mat','','error_identifier','', ...
        'error_message','','point_count',0,'ray_label',NaN);
    try
        assert(strcmp(file_sha(job.grid_mat),job.grid_sha256),'P2Exit:GridHash','Grid changed');
        assert(strcmp(file_sha([mfilename('fullpath') '.m']),job.exporter_sha256), ...
            'P2Exit:SourceHash','Exporter changed');
        clear raytrace_3d_sp;
        mex_path = which('raytrace_3d_sp');
        assert(endsWith(lower(mex_path),'.mexw64'),'P2Exit:Mex','MEX not selected');
        assert(strcmp(file_sha(mex_path),job.mex_sha256),'P2Exit:MexHash','MEX changed');
        medium = load(job.grid_mat);
        assert(all(medium.Bx_T(:)==0) && all(medium.By_T(:)==0) && ...
            all(medium.Bz_T(:)==0) && all(medium.collision_freq_Hz(:)==0), ...
            'P2Exit:Physics','Require B0 and no collisions');
        h = double(medium.height_km(:));
        bh = double(medium.geomag_height_km(:));
        assert(numel(h)<=401 && numel(bh)<=201 && all(diff(h)>0));
        iono_grid = params(medium.latitude_deg,medium.longitude_deg,h);
        geomag_grid = params(medium.geomag_latitude_deg,medium.geomag_longitude_deg,bh);
        initial = job.initial;
        expected = {'pos_x';'pos_y';'pos_z';'dir_x';'dir_y';'dir_z'; ...
            'group_path';'geom_path';'phase_path';'absorption';'indep_var';'ODE_step_size'};
        assert(isequal(fieldnames(initial),expected),'P2Exit:FieldOrder','Explicit field order');
        [ray_data,ray_path_data,ray_state_vec] = raytrace_3d_sp( ...
            30,120,50,job.source_elevation_deg,0,10,0,1,job.tol,6371000, ...
            double(medium.Ne_m3)/1e6,double(medium.Ne_m3)/1e6, ...
            medium.collision_freq_Hz,iono_grid,medium.Bx_T,medium.By_T,medium.Bz_T, ...
            geomag_grid,initial);
        raw_mat = fullfile(output_dir,[job.job_id '.mat']);
        job_sha256 = job.job_sha256;
        grid_sha256 = job.grid_sha256;
        exporter_sha256 = job.exporter_sha256;
        mex_sha256 = job.mex_sha256;
        save(raw_mat,'ray_data','ray_path_data','ray_state_vec','initial','job', ...
            'job_sha256','grid_sha256','exporter_sha256','mex_sha256','-v7');
        if ji==1, save(fullfile(fileparts(output_dir),'grid_readback.mat'),'-struct','medium','-v7'); end
        record.raw_mat = raw_mat;
        record.point_count = numel(ray_state_vec.pos_x);
        record.ray_label = ray_data.ray_label(end);
    catch ME
        result.failures = result.failures+1;
        record.error_identifier = ME.identifier;
        record.error_message = ME.message;
    end
    result.records{end+1} = record;
    write_json(summary_json,result);
    fprintf('P2 EXIT %d/%d %s error=%s\n',ji,numel(jobs),job.job_id,record.error_identifier);
end
if result.failures==0, result.status='completed'; else, result.status='partial'; end
write_json(summary_json,result);
clear raytrace_3d_sp;
assert(result.failures==0,'P2Exit:Failed','Failed jobs preserved in checkpoint');
end

function grid = params(lat,lon,h)
grid = [lat(1),lat(2)-lat(1),numel(lat),lon(1),lon(2)-lon(1),numel(lon), ...
    h(1),h(2)-h(1),numel(h)];
end

function digest = file_sha(path)
fid = fopen(path,'rb');
assert(fid>=0,'P2Exit:Read','Cannot read source');
cleanup = onCleanup(@() fclose(fid));
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(fread(fid,Inf,'*uint8'));
digest = lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end

function write_json(path,value)
fid = fopen(path,'w');
assert(fid>=0,'P2Exit:Write','Cannot write checkpoint');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(value,PrettyPrint=true));
end
