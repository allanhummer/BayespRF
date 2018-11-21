clear all

ExceptionSubjects=[1,5];
ExceptionTasks={{'eightbars_tr15','eightbars_tr15_scotom5degree'},{'eightbars_tr15_scotom25degree','eightbars_tr15_scotom5degree'}};

for SubjectNumber=2:6
    
    for Task={'eightbars_tr15','eightbars_tr15_scotom125degree','eightbars_tr15_scotom25degree','eightbars_tr15_scotom5degree'}
        
        if ~isempty(find(ExceptionSubjects==SubjectNumber,1)) && sum(strcmp(ExceptionTasks{find(ExceptionSubjects==SubjectNumber,1)},Task{1}))~=0
            Runs=1;
        else
            Runs=2;
        end

        %No Smooth
        BAYESPRF_SINGLESUBJECT('ret15',sprintf('ret15-%02d',SubjectNumber),Runs,'StimFolderName',Task{1},'FolderNameComment','s0')
        %2mm Smooth
        BAYESPRF_SINGLESUBJECT('ret15',sprintf('ret15-%02d',SubjectNumber),Runs,'StimFolderName',Task{1},'NiftiName','s2ravols','FolderNameComment','s2')
        
    end
end