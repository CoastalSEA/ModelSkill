function ms_open_manual()
%find the location of the asmita app and open the manual
appinfo = matlab.apputil.getInstalledAppInfo;
idx = find(strcmp({appinfo.name},'ModelSkill'));
fpath = [appinfo(idx(1)).location,[filesep,'ModelSkill',filesep,'doc',...
                                         filesep,'ModelSkill manual.pdf']];
open(fpath)
