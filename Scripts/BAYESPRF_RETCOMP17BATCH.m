clear all

Study='retcomp17';
ExceptionSubjects=[];%[1,5];
ExceptionTasks=[];%{{'eightbars_tr15','eightbars_tr15_scotom5degree'},{'eightbars_tr15_scotom25degree','eightbars_tr15_scotom5degree'}};
Tasks={'eightbars_hr_tr1','eightbars_lr_tr1'};
TEs=[0.0234,0.0252]; % TE HR, TE LR

for SubjectNumber=1:10
    
    for Task=1:length(Tasks)
        
        if ~isempty(find(ExceptionSubjects==SubjectNumber,1)) && sum(strcmp(ExceptionTasks{find(ExceptionSubjects==SubjectNumber,1)},Task{1}))~=0
            Runs=1;
        else
            Runs=2;
        end

        %No Smooth
        BAYESPRF_SINGLESUBJECT(Study,sprintf('%s-%02d',Study,SubjectNumber),Runs,'Session','7t','StimFolderName',Tasks{Task},StimulusName,'eightbars_tr1','TR',1,'TE',TEs(Task),'MaxStimSize',7,'NiftiName','uravols','FolderNameComment','unwarped')
        
    end
end