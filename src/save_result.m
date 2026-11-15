function save_result(path,result)
%SAVE_RESULT Finalize checkpoints atomically in the same directory.
folder=fileparts(path);if ~isfolder(folder),mkdir(folder);end
temporary=[tempname(folder) '.tmp.mat'];
save(temporary,'result','-v7');
[ok,msg]=movefile(temporary,path,'f');assert(ok,'%s',msg);
end
