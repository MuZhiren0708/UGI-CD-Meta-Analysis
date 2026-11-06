# Revised source decisions follow the original observational-study eligibility.
# All fifty reports / fifty-five original records remain visible in Fig. 2.
sensitivity_models<-function(d,fit,m,base){
 one_family<-d[!d$paper%in%c('R030','R047','R063','R037'),]
 collapsed<-read.csv(file.path(base,'source_audit_50.csv'),fileEncoding='UTF-8-BOM')
 collapsed<-collapsed[collapsed$primary==1,]
 ans<-list('Primary reported L4'=m,
  'Clear L4 definitions'=fit(d[d$definition_unclear==0,]),
  'Exclude case-control'=fit(d[!d$paper%in%c('R043','R050','R063'),]),
  'Exclude Guerra'=fit(d[d$paper!='R043',]),
  'Objective location'=fit(d[d$objective_protocol==1,]),
  'Exclude clinical-labelled'=fit(d[!trimws(d$modality_original)%in%'Clinical',]),
  'One report per source family'=fit(one_family),
  'Combine within each paper'=fit(collapsed))
 ans
}

extra_source_figures<-function(base,ml,out,fit,stat,forest_export,export,tags){
 a<-read.csv(file.path(base,'anatomic_phenotype_outcomes.csv'),fileEncoding='UTF-8-BOM')
 a<-a[a$synthesis_eligible==1,];a<-a[order(a$row),]
 definitions<-list('Strict EGD'=a[a$outcome=='EGD'&a$outcome_definition=='EGD',],
  'Isolated L4'=a[a$outcome=='Isolated broad L4',],
  'Concomitant L4'=a[a$outcome=='Concomitant broad L4',],
  'Proximal small bowel'=a[a$outcome=='Proximal small bowel',],
  'L4a / EGD proxy'=a[a$outcome=='EGD',])
 numbers<-c(6,7,8,9,12);summaries<-list();models<-list()
 for(i in seq_along(definitions)){
  label<-names(definitions)[i];z<-definitions[[i]]
  # Where two reports might concern the same source, select one per anatomical
  # outcome for this sparse secondary synthesis. Missing sites are never zero.
  z<-z[!duplicated(z$family),]
  mm<-fit(z);models[[label]]<-mm
  forest_export(mm,paste0('FigS',numbers[i]),max(4.2,2.8+.25*nrow(z)))
  ss<-stat(label,mm)
  hk<-try(suppressWarnings(update(mm,method.tau='REML',method.random.ci='HK')),silent=TRUE)
  ss$REML_HK_lower<-if(!inherits(hk,'try-error'))100*plogis(hk$lower.random)else NA_real_
  ss$REML_HK_upper<-if(!inherits(hk,'try-error'))100*plogis(hk$upper.random)else NA_real_
  summaries[[label]]<-ss
  write.csv(z,file.path(ml,paste0('anatomy_sources_',numbers[i],'.csv')),row.names=FALSE)
 }
 write.csv(do.call(rbind,summaries),file.path(ml,'anatomic_phenotype_meta.csv'),row.names=FALSE)
 saveRDS(models,file.path(ml,'anatomic_phenotype_models.rds'))

}

audit_validation<-function(base,out){
 d<-read.csv(file.path(base,'analysis_data.csv'),fileEncoding='UTF-8-BOM')
 p<-read.csv(file.path(out,'out_of_fold_predictions.csv'))
 a<-read.csv(file.path(out,'fold_membership.csv'))
 stopifnot(nrow(d)==41,length(unique(d$paper))==37,length(unique(d$family))==33,
  all(d$events>0&d$events<d$N),all(p$predicted>=0&p$predicted<=1))
 for(id in unique(a$repeat_id)){
  z<-a[a$repeat_id==id,]
  stopifnot(nrow(z)==nrow(d),!anyDuplicated(z$row),setequal(z$row,d$row),setequal(z$fold,1:5))
  for(f in 1:5){
   stopifnot(!length(intersect(z$family[z$fold==f],z$family[z$fold!=f])),
    !length(intersect(z$paper[z$fold==f],z$paper[z$fold!=f])))
  }
 }
 independent<-do.call(rbind,lapply(split(p,list(p$model,p$repeat_id),drop=TRUE),function(z){
  obs<-z$observed;pr<-z$predicted;w<-z$weight;mu<-weighted.mean(obs,w)
  data.frame(model=z$model[1],repeat_id=z$repeat_id[1],
   RMSE=100*sqrt(weighted.mean((obs-pr)^2,w)),MAE=100*weighted.mean(abs(obs-pr),w),
   R2=1-sum(w*(obs-pr)^2)/sum(w*(obs-mu)^2))}))
 reported<-read.csv(file.path(out,'repeat_metrics.csv'))
 q<-merge(independent,reported,by=c('model','repeat_id'),suffixes=c('_audit','_reported'))
 err<-max(abs(unlist(q[c('RMSE_audit','MAE_audit','R2_audit')])-unlist(q[c('RMSE_reported','MAE_reported','R2_reported')])))
 stopifnot(err<1e-9)
 write.csv(q,file.path(out,'independent_metric_audit.csv'),row.names=FALSE)
 jsonlite::write_json(list(original_reports=50,original_cohort_records=55,primary_cohorts=nrow(d),
  primary_reports=length(unique(d$paper)),possible_source_families=length(unique(d$family)),
  folds=5,repeats=length(unique(a$repeat_id)),maximum_metric_difference=err,
  family_and_publication_separation=TRUE,core_predictors=3,binary_metrics_used=FALSE,
  outcome_imputation=FALSE,negative_R2_retained=TRUE),file.path(out,'analysis_validation.json'),pretty=TRUE,auto_unbox=TRUE)
}
