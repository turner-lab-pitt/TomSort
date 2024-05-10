function output = run_spc(argstruct)

output=argstruct;
params=argstruct.params;

if size(output.spc.input_features,1)<100
    output.spc.status='Not enough input events to run SPC';
    return;
end

filename = argstruct.params.filename;
datafile = [filename '_spcinput'];
delete_executable = false;

tic;
fprintf('Preparing for SPC. Writing temporary files...');
write_temporary_files();
fprintf(' done in %0.1f seconds.\nRunning spc...',toc);
tic;
exec_spc();
output.timing.spc=toc;
fprintf(' done in %0.1f seconds. Cleaning up...',output.timing.spc);
cleanup();
fprintf(' done.\n');

    function write_temporary_files()
        if exist([filename '.dg_01.lab'],'file')
            delete([filename '.dg_01.lab']);
            delete([filename '.dg_01']);
        end
        inputfeatures=argstruct.spc.input_features; %#ok<NASGU>
        save(datafile,'inputfeatures','-ascii');
        clear inputfeatures;
        
        fid = fopen(sprintf('%s.run',filename),'wt');
        fprintf(fid,'NumberOfPoints: %s\n',num2str(argstruct.detection.num_spikes));
        fprintf(fid,'DataFile: %s\n',datafile);
        fprintf(fid,'OutFile: %s\n',filename);
        fprintf(fid,'Dimensions: %s\n',num2str(argstruct.spc.feature_count));
        fprintf(fid,'MinTemp: %s\n',num2str(params.spc.min_temp));
        fprintf(fid,'MaxTemp: %s\n',num2str(params.spc.max_temp));
        fprintf(fid,'TempStep: %s\n',num2str(params.spc.temp_step));
        fprintf(fid,'SWCycles: %s\n',num2str(params.spc.cycles));
        fprintf(fid,'KNearestNeighbours: %s\n',num2str(params.spc.k_nearest_neighbors));
        fprintf(fid,'MSTree|\n');
        fprintf(fid,'DirectedGrowth|\n');
        fprintf(fid,'SaveSuscept|\n');
        fprintf(fid,'WriteLables|\n');
        fprintf(fid,'WriteCorFile~\n');
        if isfield(params.spc, 'random_seed') && params.spc.random_seed ~= 0
            fprintf(fid,'ForceRandomSeed: %s\n',num2str(params.spc.random_seed));
        end    
        fclose(fid);
    end

    function exec_spc()
        system_type = computer;
        switch system_type
            case {'PCWIN'}    
                if exist([pwd '\cluster.exe'],'file')==0
                    directory = which('cluster.exe');
                    copyfile(directory,pwd);
                    delete_executable=[pwd '\cluster.exe'];
                end
                [output.spc.status,output.spc.result] = dos(sprintf('cluster.exe %s.run',filename));
            case {'PCWIN64'}    
                if exist([pwd '\cluster_64.exe'],'file')==0
                    directory = which('cluster_64.exe');
                    copyfile(directory,pwd);
                    delete_executable=[pwd '\cluster_64.exe'];
                end
                [output.spc.status,output.spc.result] = dos(sprintf('cluster_64.exe %s.run',filename));
            case {'MAC'}
                if exist([pwd '/cluster_mac.exe'],'file')==0
                    directory = which('cluster_mac.exe');
                    copyfile(directory,pwd);
                    delete_executable=[pwd '/cluster_mac.exe'];
                end
                run_mac = sprintf('./cluster_mac.exe %s.run',filename);
                [output.spc.status,output.spc.result] = unix(run_mac);
            case {'MACI','MACI64'}
                if exist([pwd '/cluster_maci.exe'],'file')==0
                    directory = which('cluster_maci.exe');
                    copyfile(directory,pwd);
                    delete_executable=[pwd '/cluster_maci.exe'];
                end
                run_maci = sprintf('./cluster_maci.exe %s.run',filename);
                [output.spc.status,output.spc.result] = unix(run_maci);
            case {'GLNX86'}      
                if exist([pwd '/cluster_linux.exe'],'file') == 0
                    directory = which('cluster_linux.exe');
                    copyfile(directory,pwd);
                    delete_executable=[pwd '/cluster_linux.exe'];
                end
                run_linux = sprintf('./cluster_linux.exe %s.run',filename);
                [output.spc.status,output.spc.result] = unix(run_linux);
            case {'GLNXA64', 'GLNXI64'}
                if exist([pwd '/cluster_linux64.exe'],'file') == 0
                    directory = which('cluster_linux64.exe');
                    copyfile(directory,pwd);
                    delete_executable=[pwd '/cluster_linux64.exe'];
                end
                run_linux = sprintf('./cluster_linux64.exe %s.run',filename);
                [output.spc.status,output.spc.result] = unix(run_linux);        
            otherwise 
            ME = MException('MyComponent:NotSupportedArq', '%s type of computer not supported.',com_type);
            throw(ME)
        end
        output.spc.results=load(sprintf('%s.dg_01.lab',filename));
    end

    function cleanup()
        if delete_executable~=false
            delete(delete_executable);
        end
        delete(datafile);
        delete(sprintf('%s.run',filename));
        delete([filename '.dg_01.lab']);
        delete([filename '.dg_01']);
        
        delete *.mag
        delete *.edges
        delete *.param

        if exist([filename '.knn'],'file')
            delete([filename '.knn']);
        end
    end

end

