from pathlib import Path

p = Path("R/CATER_MR_core_v05.R")
s = p.read_text()
patches = [
    ('if(min(ev)<-tol*max(1,max(abs(ev))))',
     'if(min(ev) < -tol*max(1,max(abs(ev))))'),
    ('allq<-do.call(rbind,lapply(names(candidate_by_exp),function(e){d<-candidate_by_exp[[e]];if(!nrow(d))return(NULL);d$exposure<-e;d}));',
     'core<-c("snp","chr","pos","a1","a2","beta","se","p","eaf","n");allq<-do.call(rbind,lapply(names(candidate_by_exp),function(e){d<-candidate_by_exp[[e]];if(!nrow(d))return(NULL);d<-d[,core,drop=FALSE];d$exposure<-e;d}));'),
    ('if(min(ev)<-1e-8*max(1,max(abs(ev))))',
     'if(min(ev) < -1e-8*max(1,max(abs(ev))))'),
    ('if(!is.null(trans_hits)&&nrow(trans_hits)&&is.finite(trans_reporting_p)&&trans_reporting_p>instrument_p)warning("trans_reporting_p is less stringent than instrument_p; CATER-MR will apply instrument_p to reported hits")',
     'if(!is.null(trans_hits)&&nrow(trans_hits)&&is.finite(trans_reporting_p)){if(trans_reporting_p>instrument_p)warning("trans_reporting_p is less stringent than instrument_p; CATER-MR will apply instrument_p to reported hits");if(trans_reporting_p<instrument_p)warning("trans_reporting_p is more stringent than instrument_p; the trans instrument pool is availability-censored at the source reporting threshold")}' )
]
for old, new in patches:
    if old in s:
        s = s.replace(old, new, 1)
    elif new not in s:
        raise SystemExit(f"Expected patch target not found: {old[:100]}")
p.write_text(s)
