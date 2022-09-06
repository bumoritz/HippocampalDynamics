function cmpr = repo2ws(path,files,skipIncompletelyProcessed)
%files = ["zdFF_beh"]; %files = cmpr_list_repo;

if nargin<3
    skipIncompletelyProcessed = false;
end

cmpr = {};
for i=1:length(files)
    
    % for trck structs
    if contains(files(i),'trck')
        if skipIncompletelyProcessed
            try
%                 cmpr.(char(files(i))) = load([path.filepart_in,char(files(i)),'.mat']);
%                 
%                 cmpr.(char(files(i))) = load([path.filepart_in,char(files(i)),'.mat']);
%                 cmpr.(char(files(i))) = cmpr.(char(files(i))).(cell2mat(fields(cmpr.(char(files(i)))))); 
%                 disp(['--- Loaded ',char(files(i)),' file from repo. - ',[path.filepart_in,char(files(i)),'.mat']])
            catch
            end
        else
            
        end

    % for all other structs
    else
        if skipIncompletelyProcessed
            try
                try
                    cmpr.(char(files(i))) = load([path.filepart_in,char(files(i)),'.mat']);
                    cmpr.(char(files(i))) = cmpr.(char(files(i))).(cell2mat(fields(cmpr.(char(files(i)))))); 
                    disp(['--- Loaded ',char(files(i)),' file from repo. - ',[path.filepart_in,char(files(i)),'.mat']])
                catch
                    try
                        cmpr.(char(files(i))) = load([path.filepart_in_repoX,char(files(i)),'.mat']);
                        cmpr.(char(files(i))) = cmpr.(char(files(i))).(cell2mat(fields(cmpr.(char(files(i)))))); 
                        disp(['--- Loaded ',char(files(i)),' file from repoX. - ',[path.filepart_in_repoX,char(files(i)),'.mat']])
                    catch
                        try
                            cmpr.(char(files(i))) = load([path.filepart_in_analysis,char(files(i)),'.mat']);
                            cmpr.(char(files(i))) = cmpr.(char(files(i))).(cell2mat(fields(cmpr.(char(files(i)))))); 
                            disp(['--- Loaded ',char(files(i)),' file from analysis. - ',[path.filepart_in_analysis,char(files(i)),'.mat']])
                        catch
                            cmpr.(char(files(i))) = load([path.filepart_in_analysisX,char(files(i)),'.mat']);
                            cmpr.(char(files(i))) = cmpr.(char(files(i))).(cell2mat(fields(cmpr.(char(files(i)))))); 
                            disp(['--- Loaded ',char(files(i)),' file from analysisX. - ',[path.filepart_in_analysisX,char(files(i)),'.mat']])
                        end
                    end
                end
            catch
            end

        else
            try
                cmpr.(char(files(i))) = load([path.filepart_in,char(files(i)),'.mat']);
                cmpr.(char(files(i))) = cmpr.(char(files(i))).(cell2mat(fields(cmpr.(char(files(i)))))); 
                disp(['--- Loaded ',char(files(i)),' file from repo. - ',[path.filepart_in,char(files(i)),'.mat']])
            catch
                try
                    cmpr.(char(files(i))) = load([path.filepart_in_repoX,char(files(i)),'.mat']);
                    cmpr.(char(files(i))) = cmpr.(char(files(i))).(cell2mat(fields(cmpr.(char(files(i)))))); 
                    disp(['--- Loaded ',char(files(i)),' file from repoX. - ',[path.filepart_in_repoX,char(files(i)),'.mat']])
                catch
                    try
                        cmpr.(char(files(i))) = load([path.filepart_in_analysis,char(files(i)),'.mat']);
                        cmpr.(char(files(i))) = cmpr.(char(files(i))).(cell2mat(fields(cmpr.(char(files(i)))))); 
                        disp(['--- Loaded ',char(files(i)),' file from analysis. - ',[path.filepart_in_analysis,char(files(i)),'.mat']])
                    catch
                        try
                            cmpr.(char(files(i))) = load([path.filepart_in_analysisX,char(files(i)),'.mat']);
                            cmpr.(char(files(i))) = cmpr.(char(files(i))).(cell2mat(fields(cmpr.(char(files(i)))))); 
                            disp(['--- Loaded ',char(files(i)),' file from analysisX. - ',[path.filepart_in_analysisX,char(files(i)),'.mat']])
                        catch
                            cmpr.(char(files(i))) = load([path.filepart_in_repoXX,char(files(i)),'.mat']);
                            cmpr.(char(files(i))) = cmpr.(char(files(i))).(cell2mat(fields(cmpr.(char(files(i)))))); 
                            disp(['--- Loaded ',char(files(i)),' file from repoXX. - ',[path.filepart_in_repoXX,char(files(i)),'.mat']])
                        end
                    end
                end
            end
        end
    end
end

end