function [blockOrder,attBlockOrder,targetBlockOrder, cueBlockOrder, exoBlockOrder] = block_gen_exo(blockNames, attBlockNames, targetBlockNames, cueBlockNames, exoBlockNames)
% Random block generator for makeTADetectStim: one run (64 trials in one 
% repetition = 4 trials for each target (4) x cue (2) condition ) with
% added blank trials every 4 target trials 
% Assume cue validity = 75% ( 3 trials valid and 1 trial invalid )
% One blockOrder condition: fast-left or slow-left
% One attBlockOrder condition: att-right


%% define indices
blank = find(ismember(blockNames,'blank'));
fast_side = find(ismember(blockNames,{'fast-left','slow-left'}));

no_att = find(ismember(attBlockNames,'no-att'));
att_right = find(ismember(attBlockNames,'att-right'));

nt = find(ismember(targetBlockNames,'no-targ'));
ll = find(ismember(targetBlockNames,'left-left'));
lr = find(ismember(targetBlockNames,'left-right'));
rl = find(ismember(targetBlockNames,'right-left'));
rr = find(ismember(targetBlockNames,'right-right'));

ne = find(ismember(exoBlockNames,'no-exo'));
pp = find(ismember(exoBlockNames,'pres-pres'));
pa = find(ismember(exoBlockNames,'pres-abs'));
ap = find(ismember(exoBlockNames,'abs-pres'));
aa = find(ismember(exoBlockNames,'abs-abs'));

nc   = find(ismember(cueBlockNames,'no-cue'));
c1c1 = find(ismember(cueBlockNames,'1-1'));
c1c2 = find(ismember(cueBlockNames,'1-2'));
c2c1 = find(ismember(cueBlockNames,'2-1'));
c2c2 = find(ismember(cueBlockNames,'2-2'));

%% cue order
% pre-cue = T1
A = [c1c1,c1c1,c1c1,c1c2]; % specify validity for pre-cue = T1
cueBlockOrder_cue1 = repmat(A,1,16); % for all target conditions
% 
% pre-cue = T2
B = [c2c2,c2c2,c2c2,c2c1]; % specify validity for cue = T2
cueBlockOrder_cue2 = repmat(B,1,16); % for all target conditions

cueBlockOrder = [cueBlockOrder_cue1 , cueBlockOrder_cue2];

% randomize cue order
indices = randperm(length(cueBlockOrder));
cueBlockOrder = cueBlockOrder(indices);

%% target order
dummy = repmat([ll,lr,rl,rr],[16,1]);
targetBlockOrder_cue1 = dummy(:)';
targetBlockOrder = repmat(targetBlockOrder_cue1,[1,2]);

% randomize target order
targetBlockOrder = targetBlockOrder(indices);

%% exo order
dummy = repmat([pp,pa,ap,aa],[4,1]);
exoBlockOrder_cue1 = dummy(:)';
exoBlockOrder = repmat(exoBlockOrder_cue1,[1,16]);

% randomize target order
exoBlockOrder = exoBlockOrder(indices);

%% insert blank trials 
% (every 4 target trials for target and cue block order)
blank_tar = repmat(nt,[1,(length(targetBlockOrder)/4) + 1]);
ind = zeros(1, length(targetBlockOrder)+ length(blank_tar));
ind (1:5:length(ind)) = blank_tar;
ind2 = ind;
ind2 (ind == nt ) = nc;
ind3 (ind == nt ) = ne;
ind(ind == 0) = targetBlockOrder;
targetBlockOrder = ind;
ind2(ind2 == 0) = cueBlockOrder;
cueBlockOrder = ind2;
ind3(ind3 == 0) = exoBlockOrder;
exoBlockOrder = ind3;

%% block order (one condition) and attention order
block = [blank,repmat(fast_side,1,4)];
blockOrder = [repmat(block,1,32),blank];

att = [no_att,repmat(att_right,1,4)];
attBlockOrder = [repmat(att,1,32),no_att];


end

