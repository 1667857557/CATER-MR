from pathlib import Path

p = Path("R/CATER_MR_core_v05.R")
s = p.read_text()
old1 = 'if(min(ev)<-tol*max(1,max(abs(ev))))'
new1 = 'if(min(ev) < -tol*max(1,max(abs(ev))))'
old2 = 'allq<-do.call(rbind,lapply(names(candidate_by_exp),function(e){d<-candidate_by_exp[[e]];if(!nrow(d))return(NULL);d$exposure<-e;d}));'
new2 = 'core<-c("snp","chr","pos","a1","a2","beta","se","p","eaf","n");allq<-do.call(rbind,lapply(names(candidate_by_exp),function(e){d<-candidate_by_exp[[e]];if(!nrow(d))return(NULL);d<-d[,core,drop=FALSE];d$exposure<-e;d}));'
for old, new in ((old1,new1),(old2,new2)):
    if old in s:
        s = s.replace(old,new,1)
    elif new not in s:
        raise SystemExit(f"Expected patch target not found: {old[:80]}")
p.write_text(s)
