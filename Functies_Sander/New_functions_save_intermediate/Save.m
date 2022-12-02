function Save(alg)

tSave = tic;

fileName = alg.name; % alg has no name
temp = split(fileName, '/');
if numel(temp) > 1
    folder = sprintf('%s%s', char(join(temp(1:end-1), '/')), '/');
    if ~isfolder(folder), mkdir(folder); end
end

% construct data
fileData = struct;
switch alg.options.saveMethod
    case 'full'
        fileData.mps            = vumps.mps;
        fileData.environment    = vumps.environment;
        fileData.lambda         = vumps.lambda;
        fileData.error          = vumps.error;
        fileData.chi            = vumps.mps(1).bondDimension(1);
        fileData.config         = vumps.options;
        fileData.operator       = vumps.operator;
        fileData.history        = vumps.history;
        fileData.mode           = 'full';
    case 'minimal'
        fileData.A              = vumps.mps.AR;
        fileData.C              = vumps.mps.C;
        fileData.lambda         = vumps.lambda;
        fileData.error          = vumps.error;
        fileData.chi            = vumps.mps(1).bondDimension(1);
        fileData.mode           = 'minimal';
end

% save
if exist(fileName,'file')
    old_file=load(fileName);
    fileName_temp=[fileName(1:end-4),'_temp.mat'];
    save(fileName_temp, '-struct', 'old_file', '-v7.3');
    saved_temp=1;
else
    saved_temp=0;
end

save(fileName, '-struct', 'fileData', '-v7.3');

if saved_temp
    delete(fileName_temp);
end

vumps.Communication('save', toc(tSave), ByteSize(fileData, 'B'));

end
