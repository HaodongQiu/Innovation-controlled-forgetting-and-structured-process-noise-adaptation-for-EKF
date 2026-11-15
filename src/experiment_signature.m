function signature = experiment_signature
%EXPERIMENT_SIGNATURE Hash all simulation code and input data for safe resume.
root=fileparts(fileparts(mfilename('fullpath')));
files=[dir(fullfile(root,'src','*.m'));dir(fullfile(root,'data','**','*.mat'))];
names=string(fullfile({files.folder},{files.name}));[~,order]=sort(names);files=files(order);
md=java.security.MessageDigest.getInstance('SHA-256');
for k=1:numel(files)
 path=fullfile(files(k).folder,files(k).name);
 relative=path(numel(root)+2:end);md.update(unicode2native(strrep(relative,'\','/'),'UTF-8'));
 f=fopen(path,'rb');assert(f>=0);c=onCleanup(@()fclose(f));
 while ~feof(f),md.update(fread(f,1024*1024,'*uint8'));end
 clear c
end
signature=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
