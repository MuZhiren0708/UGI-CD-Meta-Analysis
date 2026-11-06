# Reproduce corrected synthesis and all SVGs; optionally refit the frozen ML analysis.
options(stringsAsFactors=FALSE)
argv<-commandArgs(trailingOnly=TRUE)
script<-sub('^--file=','',commandArgs()[grepl('^--file=',commandArgs())])
here<-dirname(normalizePath(script,winslash='/',mustWork=TRUE))
root<-normalizePath(file.path(here,'..'),winslash='/',mustWork=TRUE)
out<-if(length(argv))argv[1]else file.path(root,'reproduced_results')
mode<-if(length(argv)>1)argv[2]else 'saved'
stopifnot(mode%in%c('saved','full'))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
out<-normalizePath(out,winslash='/',mustWork=TRUE)
rscript<-file.path(R.home('bin'),'Rscript')
if(.Platform$OS.type=='windows')rscript<-paste0(rscript,'.exe')
copy_dir<-function(from,to){
 dir.create(to,recursive=TRUE,showWarnings=FALSE)
 ff<-list.files(from,full.names=TRUE)
 stopifnot(all(file.copy(ff,to,recursive=TRUE,overwrite=TRUE)))
}
run<-function(name,args){
 status<-system2(rscript,c(shQuote(file.path(here,name)),vapply(args,shQuote,character(1))))
 if(status!=0)stop('R script failed: ',name,'; status ',status)
}
if(mode=='saved'){
 copy_dir(file.path(root,'analysis_results'),file.path(out,'analysis_results_final'))
 copy_dir(file.path(root,'ML_demographic_sensitivity'),file.path(out,'ML_demographic_sensitivity'))
}
# Keep SVG text editable and avoid optional Node/Python dependencies in the standard run.
Sys.setenv(UGI_FINISH_SVG='0')
run('UGI_CD_Meta_Analysis.R',c(out,if(mode=='full')'all'else'figures'))
run('Untitled1.R',out)
extra<-list.files(file.path(out,'SVG_extra'),pattern='\\.svg$',full.names=TRUE)
stopifnot(all(file.copy(extra,file.path(out,'SVG_figures'),overwrite=TRUE)))
copy_dir(file.path(root,'figure_data'),file.path(out,'figure_data'))
for(ff in list.files(file.path(out,'figure_data'))){
 fresh<-file.path(out,'analysis_results_final',ff)
 if(file.exists(fresh))file.copy(fresh,file.path(out,'figure_data',ff),overwrite=TRUE)
}
# These frozen input rows are unchanged by the quality-total correction.
run('Redraw_Figures_07_11.R',out)
run('UGI_CD_Meta_Analysis.R',c(out,'validate'))
stopifnot(length(list.files(file.path(out,'SVG_figures'),pattern='\\.svg$'))==22)
cat('Reproduction complete. Results:',out,'\n')
