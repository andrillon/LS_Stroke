function [real_out, cont_out, perm_out, cont_perm_out, out_pred_perm]=lme_perm_lsstroke(table,predictor,formula,totperm,all_pred_perm)
% table=GO_table;
%
% formula='SW~1+Group+(1|SubID)';
% permvars={'SubID'};
% totperm=100;
if nargin<5
    all_pred_perm=[];
end
out_pred_perm=[];
% run real model
if strcmp(predictor,'subGroupID')
    eval(sprintf('table.pred=table.%s;',predictor));
    model= fitlme(table,formula);
    ContNames=unique(table.pred);
    real_out=[];
    cont_out=[];
    for nCond=1:2
        real_out=[real_out ; double(model.Coefficients(2+nCond,2)) double(model.Coefficients(2+nCond,4)) double(model.Coefficients(2+nCond,6)) nCond];
        cont_out=[cont_out ; {sprintf('%s vs %s',char(ContNames(nCond+1)),char(ContNames(1)))}];
    end
    table2=table;
    table2.pred=reordercats(table2.pred,[2 1 3]);
    ContNames2=unique(table2.pred);
    model2= fitlme(table2,formula);
    real_out=[real_out ; double(model2.Coefficients(2+nCond,2)) double(model2.Coefficients(2+nCond,4)) double(model2.Coefficients(2+nCond,6)) nCond+1];
    cont_out=[cont_out ; {sprintf('%s vs %s',char(ContNames2(nCond+1)),char(ContNames2(1)))}];
    
    
    
    uniqueIDs=unique(table.SubID);
    for nSub=1:length(uniqueIDs)
        uniqueGroups(nSub)=unique(table.subGroupID(table.SubID==uniqueIDs(nSub)));
    end
    perm_out=[];
    cont_perm_out=[];
    fprintf('%4.0f/%4.0f\n',0,totperm)
    for np=1:totperm
        if size(all_pred_perm,1)~=totperm
            group_perm_idx=1:length(uniqueGroups);
            group_perm_idx=group_perm_idx(randperm(length(group_perm_idx)));
            out_pred_perm(np,:)=group_perm_idx;
        else
            group_perm_idx=all_pred_perm(np,:)';
            out_pred_perm(np,:)=group_perm_idx;
        end
        uniqueGroups_perm=uniqueGroups(group_perm_idx);
        table_perm=table;
        for nS=1:length(uniqueIDs)
            table_perm.pred(table.SubID==uniqueIDs(nS))=uniqueGroups_perm(nS);
        end
        model_perm= fitlme(table_perm,formula);
        ContNames=unique(table_perm.pred);
        for nCond=1:2
            perm_out=[perm_out ; double(model_perm.Coefficients(2+nCond,2)) double(model_perm.Coefficients(2+nCond,4)) double(model_perm.Coefficients(2+nCond,6)) nCond np];
            cont_perm_out=[cont_perm_out ; {sprintf('%s vs %s',char(ContNames(nCond+1)),char(ContNames(1)))}];
        end
        table2_perm=table_perm;
        table2_perm.pred=reordercats(table2_perm.pred,[2 1 3]);
        ContNames2=unique(table2_perm.pred);
        model2_perm= fitlme(table2_perm,formula);
        perm_out=[perm_out ; double(model2_perm.Coefficients(2+nCond,2)) double(model2_perm.Coefficients(2+nCond,4)) double(model2_perm.Coefficients(2+nCond,6)) nCond+1 np];
        cont_perm_out=[cont_perm_out ; {sprintf('%s vs %s',char(ContNames2(nCond+1)),char(ContNames2(1)))}];
        fprintf('\b\b\b\b\b\b\b\b\b\b%4.0f/%4.0f\n',np,totperm)
    end
    fprintf('\n');
    
else
    %     eval(sprintf('table.pred=table.%s;',predictor));
    %     model= fitlme(table,formula);
    %     real_out=[double(model.Coefficients(match_str(model.CoefficientNames,'pred'),2)) double(model.Coefficients(match_str(model.CoefficientNames,'pred'),4)) double(model.Coefficients(match_str(model.CoefficientNames,'pred'),6))];
    %
    %     uniqueIDs=unique(table.SubID);
    %     uniqueBlocks=unique(table.BlockN);
    %     perm_out=nan(totperm,4);
    %     fprintf('%4.0f/%4.0f\n',0,totperm)
    %     for np=1:totperm
    %         if size(all_pred_perm,1)~=totperm
    %             pred_perm_idx=1:length(table.pred);
    %             for nS=1:length(uniqueIDs)
    %                 idx=find(table.SubID==uniqueIDs(nS));
    %                 temp=pred_perm_idx(idx);
    %                 pred_perm_idx(idx)=temp(randperm(length(temp)));
    %             end
    %             out_pred_perm(np,:)=pred_perm_idx;
    %         else
    %             pred_perm_idx=all_pred_perm(np,:)';
    %         end
    %
    %         table2=table;
    %         table2.pred=table2.pred(pred_perm_idx);
    %         model= fitlme(table2,formula);
    %         perm_out(np,:)=[double(model.Coefficients(match_str(model.CoefficientNames,'pred'),2)) double(model.Coefficients(match_str(model.CoefficientNames,'pred'),4)) double(model.Coefficients(match_str(model.CoefficientNames,'pred'),6)) np];
    %         fprintf('\b\b\b\b\b\b\b\b\b\b%4.0f/%4.0f\n',np,totperm)
    %     end
    %     fprintf('\n');
    
end