from pathlib import Path
p=Path('R/CATER_MR_core_v05.R')
s=p.read_text()
old='if(length(ii)) sumdf$primary_q<-p.adjust(sumdf$primary_p[ii],method="BH")'
new='if(length(ii)) sumdf$primary_q[ii]<-p.adjust(sumdf$primary_p[ii],method="BH")'
assert s.count(old)==1
p.write_text(s.replace(old,new,1))
