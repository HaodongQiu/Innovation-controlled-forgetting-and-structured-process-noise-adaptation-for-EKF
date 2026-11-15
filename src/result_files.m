function files = result_files(folder)
%RESULT_FILES Ignore temporary checkpoints and unrelated MAT files.
files=dir(fullfile(folder,'case_*.mat'));
keep=arrayfun(@(f)~isempty(regexp(f.name,'^case_\d{4}_[A-Za-z0-9-]+\.mat$','once')),files);
files=files(keep);
end
