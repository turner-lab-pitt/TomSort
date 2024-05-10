function output = spc_iteration(argstruct)

params=argstruct.params;

if size(argstruct.data,1)<params.spc.min_cluster
    output='Not enough input events to run SPC';
    return;
end

filename = strsplit(argstruct.params.filename,'\');
filename = filename{end};
datafile = [filename '_spcinput'];
delete_executable = false;


write_temporary_files();
exec_spc();
cleanup();

    function write_temporary_files()
        if exist([filename '.dg_01.lab'],'file')
            delete([filename '.dg_01.lab']);
            delete([filename '.dg_01']);
        end
        inputfeatures=argstruct.data;
        [ns,nf]=size(inputfeatures);
        
        save_tries = 25;
        save_try=0;
        while true
            save_try = save_try+1;
            try
                save(datafile,'inputfeatures','-ascii');
                clear inputfeatures;
                break;
            catch ME
                if save_try >= save_tries
                    throw(ME);
                end
                pause(.5);%wait half a second and try again
            end
        end
        
        fid = fopen(sprintf('%s.run',filename),'wt');
        fprintf(fid,'NumberOfPoints: %s\n',num2str(ns));
        fprintf(fid,'DataFile: %s\n',datafile);
        fprintf(fid,'OutFile: %s\n',filename);
        fprintf(fid,'Dimensions: %s\n',num2str(nf));
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
                [~,~] = dos(sprintf('cluster.exe %s.run',filename));
            case {'PCWIN64'}    
                if exist([pwd '\cluster_64.exe'],'file')==0
                    directory = which('cluster_64.exe');
                    copyfile(directory,pwd);
                    delete_executable=[pwd '\cluster_64.exe'];
                end
                [~,~] = dos(sprintf('cluster_64.exe %s.run',filename));
            case {'MAC'}
                if exist([pwd '/cluster_mac.exe'],'file')==0
                    directory = which('cluster_mac.exe');
                    copyfile(directory,pwd);
                    delete_executable=[pwd '/cluster_mac.exe'];
                end
                run_mac = sprintf('./cluster_mac.exe %s.run',filename);
                [~,~] = unix(run_mac);
            case {'MACI','MACI64'}
                if exist([pwd '/cluster_maci.exe'],'file')==0
                    directory = which('cluster_maci.exe');
                    copyfile(directory,pwd);
                    delete_executable=[pwd '/cluster_maci.exe'];
                end
                run_maci = sprintf('./cluster_maci.exe %s.run',filename);
                [~,~] = unix(run_maci);
            case {'GLNX86'}      
                if exist([pwd '/cluster_linux.exe'],'file') == 0
                    directory = which('cluster_linux.exe');
                    copyfile(directory,pwd);
                    delete_executable=[pwd '/cluster_linux.exe'];
                end
                run_linux = sprintf('./cluster_linux.exe %s.run',filename);
                [~,~] = unix(run_linux);
            case {'GLNXA64', 'GLNXI64'}
                if exist([pwd '/cluster_linux64.exe'],'file') == 0
                    directory = which('cluster_linux64.exe');
                    copyfile(directory,pwd);
                    delete_executable=[pwd '/cluster_linux64.exe'];
                end
                run_linux = sprintf('./cluster_linux64.exe %s.run',filename);
                [~,~] = unix(run_linux);        
            otherwise 
            ME = MException('MyComponent:NotSupportedArq', '%s type of computer not supported.',com_type);
            throw(ME)
        end
        
        load_tries = 5;
        load_try=0;
        while true
            load_try = load_try+1;
            try
                output=load(sprintf('%s.dg_01.lab',filename));
                
                break;
            catch ME
                if load_try >= load_tries
                    throw(ME);
                end
                pause(.5);%wait half a second and try again
            end
        end
        
        %output=load(sprintf('%s.dg_01.lab',filename));
    end

    function cleanup()
        if delete_executable~=false
            delete(delete_executable);
        end
        delete(datafile);
        delete(sprintf('%s.run',filename));
        delete([filename '.dg_01.lab']);
        delete([filename '.dg_01']);
        
        delete([filename '*.mag']);
        delete([filename '*.edges']);
        delete([filename '*.param']);

        if exist([filename '.knn'],'file')
            delete([filename '.knn']);
        end
    end

end

