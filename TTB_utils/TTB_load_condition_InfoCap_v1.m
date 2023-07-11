function [vcondition]=load_condition_Sleep_InfoCap(trial,cond)

for i=1:trial
    load(sprintf('Wtrials_%03d_%d.mat',i,cond));
    vcondition.infocap_all(i,:)=(infocapacity);
    vcondition.suscep_all(i,:)=(susceptibility);
    vcondition.turbu_all(i,:,:) = TurbulenceSim_sub;
end