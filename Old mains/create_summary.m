%% Create lightweight summary file (without heavy simulation data)
% This extracts just the key results, making a ~10 MB file instead of 270 GB

outputDir = 'sweep_results';
files = dir(fullfile(outputDir, 'run*.mat'));
nRuns = length(files);

fprintf('Creating summary from %d runs...\n', nRuns);

% Preallocate structure array
summary = struct();
summary.runs = repmat(struct(...
    'run_index', [], ...
    'chunk_id', [], ...
    'SF_flag', [], ...
    'SF_efficiency', [], ...
    'junction_depth', [], ...
    'sp_r', [], ...
    'nemitter_EF0', [], ...
    'rightelectrode', [], ...
    'N_d', [], ...
    'Voc', [], ...
    'PCE', [], ...
    'MPP_power', [], ...
    'Pin', [], ...
    'Vapp', [], ...
    'Jtot', [], ...
    'EQE_wavelengths', [], ...
    'EQE', [], ...
    'IQE', [], ...
    'Jsc_vs_lambda', [], ...
    'time', []), nRuns, 1);

% Load each run and extract key data
for i = 1:nRuns
    if mod(i, 100) == 0
        fprintf('  Processing %d/%d...\n', i, nRuns);
    end
    
    try
        data = load(fullfile(outputDir, files(i).name), 'thisRun');
        r = data.thisRun;
        
        % Copy lightweight fields only (no par, solEq, solCV - those are huge!)
        summary.runs(i).run_index = r.run_index;
        summary.runs(i).chunk_id = r.chunk_id;
        summary.runs(i).SF_flag = r.SF_flag;
        summary.runs(i).SF_efficiency = r.SF_efficiency;
        summary.runs(i).junction_depth = r.junction_depth;
        summary.runs(i).sp_r = r.sp_r;
        summary.runs(i).nemitter_EF0 = r.nemitter_EF0;
        summary.runs(i).rightelectrode = r.rightelectrode;
        summary.runs(i).N_d = r.N_d;
        summary.runs(i).Voc = r.Voc;
        summary.runs(i).PCE = r.PCE;
        summary.runs(i).MPP_power = r.MPP_power;
        summary.runs(i).Pin = r.Pin;
        summary.runs(i).Vapp = r.Vapp;
        summary.runs(i).Jtot = r.Jtot;
        summary.runs(i).EQE_wavelengths = r.EQE_wavelengths;
        summary.runs(i).EQE = r.EQE;
        summary.runs(i).IQE = r.IQE;
        summary.runs(i).Jsc_vs_lambda = r.Jsc_vs_lambda;
        summary.runs(i).time = r.time;
    catch ME
        fprintf('  Warning: Could not load run %d: %s\n', i, ME.message);
    end
end

% Add metadata
summary.total_runs = nRuns;
summary.created = datetime('now');
summary.description = 'Lightweight summary of parameter sweep results';

% Save lightweight summary
summaryFile = fullfile(outputDir, 'summary_lightweight.mat');
save(summaryFile, 'summary', '-v7.3');

fileInfo = dir(summaryFile);
fprintf('\n✓ Created summary file: %s\n', summaryFile);
fprintf('  Size: %.1f MB (vs %.1f GB for all runs)\n', ...
        fileInfo.bytes/1e6, nRuns*30/1000);
fprintf('\nDownload this file to your laptop for analysis!\n');