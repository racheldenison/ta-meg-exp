% rd_plotSeqAllTrials.m

load([pathToExpt('stimuli') '/taDetectDiscrim401.mat'])

blockLengths = stimulus.itiSeq + p.blockDur;
blockStarts = [0 cumsum(blockLengths(1:end-1))];
blockIdx = blockStarts*p.refrate + 1;
nBlocks = numel(blockIdx);

window = 250;
a = [];
for i = 1:nBlocks
    a(i,:) = stimulus.seq(blockIdx(i):blockIdx(i)+window-1);
end

figure
plot(a')
xlabel('time (frames)')
ylabel('image number in stimulus.seq')

offsets = repmat(1:nBlocks,window,1);
figure
plot(a' + offsets);
xlabel('time (frames)')
ylabel('image number in stimulus.seq')