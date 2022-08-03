function [sig_con, ccg_vec, Bounds, ccg_vec_inh, sig_con_inh,ccgR, Pval]=ccgMCS(spikes, ccg_params)

%TH network inference CCG calculation script 

binSize = ccg_params.binSize; %.5ms
duration = ccg_params.duration; %200ms
epoch = ccg_params.epoch; %whole session
conv_w = ccg_params.conv_w;  % 10ms window
alpha = ccg_params.alpha; %high frequency cut off, must be .001 for causal p-value matrix
Fs = ccg_params.Fs;

nCel=length(spikes);

% Create CCGs (including autoCG) for all cells
[ccgR1,tR] = CCG(spikes,[],'binSize',binSize,'duration',duration,'Fs',Fs);

ccgR = nan(size(ccgR1,1),nCel,nCel);
ccgR(:,1:size(ccgR1,2),1:size(ccgR1,2)) = ccgR1;


% get  CI for each CCG
Pval=nan(length(tR),nCel,nCel);
Pred=zeros(length(tR),nCel,nCel);
Bounds=zeros(size(ccgR,1),nCel,nCel);
sig_con = [];
sig_con_inh = [];

TruePositive = nan(nCel,nCel);
FalsePositive = nan(nCel,nCel);
Pcausal = nan(nCel,nCel);



for refcellID=1:max(nCel)
    for cell2ID=1:max(nCel)
        
        if(refcellID==cell2ID)
           continue; 
        end
        
        cch=ccgR(:,refcellID,cell2ID);			% extract corresponding cross-correlation histogram vector
        
%         refcellshank=completeIndex(completeIndex(:,3)==refcellID);
%         cell2shank=completeIndex(completeIndex(:,3)==cell2ID);
%         if refcellshank==cell2shank
%             
%             
%             % central 1.6 ms on same-shank = NaN due to limitations of
%             % extracting overlapping spikes
%             
%             centerbins = ceil(length(cch)/2);
%             sameshankcch=cch;
%             sameshankcch(centerbins)=[];
%             
%             [pvals,pred,qvals]=bz_cch_conv(sameshankcch,conv_w);
%             pred=[pred(1:(centerbins(1)-1));nan(length(centerbins),1);pred(centerbins(end)-length(centerbins)+1:end)];
%             
%             pvals=[pvals(1:(centerbins(1)-1));nan(length(centerbins),1);pvals(centerbins(end)-length(centerbins)+1:end)];
%         else
%             % calculate predictions using Eran's bz_cch_conv
            
            [pvals,pred,qvals]=bz_cch_conv(cch,conv_w);
%         end
        
        
        
        % Store predicted values and pvalues for subsequent plotting
        Pred(:,refcellID,cell2ID)=pred;
        Pval(:,refcellID,cell2ID)=pvals(:);
        Pred(:,cell2ID,refcellID)=flipud(pred(:));
        Pval(:,cell2ID,refcellID)=flipud(pvals(:));
        
        % Calculate upper and lower limits with bonferonni correction
        % monosynaptic connection will be +/- 4 ms
        
        nBonf = round(.005/binSize)*2;
        hiBound=poissinv(1-alpha/nBonf,pred);
        loBound=poissinv(alpha/nBonf, pred);
        Bounds(:,refcellID,cell2ID,1)=hiBound;
        Bounds(:,refcellID,cell2ID,2)=loBound;
        
        Bounds(:,cell2ID,refcellID,1)=flipud(hiBound(:));
        Bounds(:,cell2ID,refcellID,2)=flipud(loBound(:));
        
        sig = cch>hiBound; 
        sig_inh= cch < loBound;
        %sig = cch>hiBound;
        
        % Find if significant periods falls in monosynaptic window +/- 4ms
        prebins = round(length(cch)/2 - .0032/binSize):round(length(cch)/2);
        postbins = round(length(cch)/2 + .0008/binSize):round(length(cch)/2 + .004/binSize);
        cchud  = flipud(cch);
        sigud  = flipud(sig);
        sigud_inh=flipud(sig_inh);
        
        sigpost=max(cch(postbins))>poissinv(1-alpha,max(cch(prebins)));
        sigpre=max(cchud(postbins))>poissinv(1-alpha,max(cchud(prebins)));
        
        sigpost_inh=min(cch(postbins))<poissinv(alpha,mean(cch(prebins)));
        sigpre_inh=min(cchud(postbins))<poissinv(alpha,mean(cchud(prebins)));
%         
%         
%         %define likelihood of being a connection
%         pvals_causal = 1 - poisscdf( max(cch(postbins)) - 1, max(cch(prebins) )) - poisspdf( max(cch(postbins)), max(cch(prebins)  )) * 0.5;
%         pvals_causalud = 1 - poisscdf( max(cchud(postbins)) - 1, max(cchud(prebins) )) - poisspdf( max(cchud(postbins)), max(cchud(prebins)  )) * 0.5;
%         
%         %can go negative for very small p-val - beyond comp. sig. dig
%         
%         if pvals_causalud<0
%             pvals_causalud = 0;
%         end
%         
%         if pvals_causal<0
%             pvals_causal = 0;
%         end
        
%         
%         Pcausal(refcellID,cell2ID) = pvals_causal;
%         Pcausal(cell2ID,refcellID) = pvals_causalud;
%         
%             if any(Pval(postbins,cell2ID,refcellID)<.001)
%                 
%                 FP =  v.ProbSyn.FalsePositive((histc(pvals_causalud,v.ProbSyn.thres))>0);
%                 TP =  v.ProbSyn.TruePositive((histc(pvals_causalud,v.ProbSyn.thres))>0);
%                 TruePositive(cell2ID,refcellID) = TP;
%                 FalsePositive(cell2ID,refcellID) = FP;
%             end
%             
%             
%             if any(Pval(postbins,refcellID,cell2ID)<.001)
%                 
%                 FP =  v.ProbSyn.FalsePositive((histc(pvals_causal,v.ProbSyn.thres))>0);
%                 TP =  v.ProbSyn.TruePositive((histc(pvals_causal,v.ProbSyn.thres))>0);
%                 TruePositive(refcellID,cell2ID) = TP;
%                 FalsePositive(refcellID,cell2ID) = FP;
%             end
       
        
        %check which is bigger
        if (any(sigud(prebins)) && sigpre)
            
            %test if causal is bigger than anti causal
            
            
            sig_con = [sig_con;cell2ID refcellID];
            
        end
        
        if (any(sig(postbins)) && sigpost)
            
            sig_con = [sig_con;refcellID cell2ID];
        end
        
        
          if (any(sigud_inh(prebins)) && sigpre_inh)
            
            %test if causal is bigger than anti causal
            
            
            sig_con_inh = [sig_con_inh;cell2ID refcellID];
            
        end
        
        if (any(sig_inh(postbins)) && sigpost_inh)
            
            sig_con_inh = [sig_con_inh;refcellID cell2ID];
        end
        
        
        
        
        
    end
    
end

sig_con=unique(sig_con,'rows');
sig_con_inh=unique(sig_con_inh, 'rows');


if(~isempty(sig_con))

ccg_vec=[];
for jj=1:size(sig_con,1)

    ccg_vec=[ccg_vec ccgR(:,sig_con(jj,1),sig_con(jj,2))];
    
    
end


else
    ccg_vec=[];
end



if(~isempty(sig_con_inh))

ccg_vec_inh=[];
for jj=1:size(sig_con_inh,1)

    ccg_vec_inh=[ccg_vec_inh ccgR(:,sig_con_inh(jj,1),sig_con_inh(jj,2))];
    
    
end


else
    ccg_vec_inh=[];
end




end