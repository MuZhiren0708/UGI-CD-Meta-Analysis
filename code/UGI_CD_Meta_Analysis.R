# Source-verified revision of the ORIGINAL UGI_CD_Meta_Analysis.R plotting framework.
# Blue/red forest plots, original map, panel layouts and original palettes retained.
# Original manuscript eligibility retained; all50 reports /55 rows remain visible.
# Primary study-reported L4:41 cohorts /37 reports;33 conservative validation families.
# Valid original assessment times retained, including Fumery227/535 at follow-up.
# Original references, manuscript table diagnosis labels and figure numbers are protected.
# Run: Rscript UGI_CD_Meta_Analysis.R [output folder] [all|ml|figures|validate]
# All means, model scores and source outcomes are computed; none are assigned.
options(stringsAsFactors=FALSE,scipen=99)
if(.Platform$OS.type=='windows')suppressWarnings(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
extra<-Sys.getenv('UGI_R_LIBRARY','')
if(nzchar(extra)).libPaths(c(strsplit(extra,.Platform$path.sep,fixed=TRUE)[[1]],.libPaths()))
args<-commandArgs(trailingOnly=TRUE)
script<-sub('^--file=','',commandArgs()[grepl('^--file=',commandArgs())])
here<-if(length(script))dirname(normalizePath(script,winslash='/'))else getwd()
output_root<-if(length(args)>=1)args[1]else file.path(here,'results')
stage<-if(length(args)>=2)args[2]else 'all'
stopifnot(stage%in%c('all','ml','figures','validate','ml_sensitivity'))
require_pkgs<-function(pkgs){for(p in pkgs)if(!requireNamespace(p,quietly=TRUE))stop('Required package: ',p)}
data_dir<-file.path(output_root,'data');ml_dir<-file.path(output_root,'analysis_results_final');figure_dir<-file.path(output_root,'SVG_figures')
for(x in c(output_root,data_dir,ml_dir,figure_dir))dir.create(x,showWarnings=FALSE,recursive=TRUE)
for(ff in list.files(file.path(here,'..','data'),full.names=TRUE))file.copy(ff,data_dir,overwrite=TRUE)
source(file.path(here,'source_sensitivity.R'),encoding='UTF-8')
source(file.path(here,'figures_10_11.R'),encoding='UTF-8')
run_ml<-function(base,out,exact_demographics=FALSE){
 require_pkgs(c('glmnet','ranger','xgboost','metafor','jsonlite'))
 d<-read.csv(file.path(base,'analysis_data.csv'),fileEncoding='UTF-8-BOM',check.names=FALSE)
 stopifnot(nrow(d)>=10,all(d$events>0&d$events<d$N),nrow(d)==41)
 if(exact_demographics)d<-d[d$demographic_exact_denominator==1,]
 dir.create(out,showWarnings=FALSE,recursive=TRUE)
 d$p<-d$events/d$N;d$y<-qlogis(d$p)
 SEED<-20261028L;REPEATS<-as.integer(Sys.getenv('UGI_CV_REPEATS','10'))
 stopifnot(REPEATS>=1)
 lambda_grid<-c(.5,.3,.15,.08,.04,.02,.01,.005,.001)
 core_names<-c('Age','Male','LogN')
 models<-c('LASSO','Random forest','XGBoost','REML meta-regression','Mean-only baseline','LASSO + L3')
 family_w<-function(dat){tab<-table(dat$family);as.numeric(1/tab[dat$family])}
 folds_grouped<-function(dat,k=5,seed=SEED){
  groups<-unique(dat$family);paediatric<-sapply(groups,function(g)any(dat$age_group[dat$family==g]=='Paediatric'))
  set.seed(seed);ans<-setNames(integer(length(groups)),groups)
  for(flag in c(TRUE,FALSE)){
   gg<-sample(groups[paediatric==flag]);if(length(gg))ans[gg]<-rep(sample(seq_len(k)),length.out=length(gg))
  }
  as.integer(ans[dat$family])
 }
 raw_x<-function(dat,set='core',med_l3=NULL,med_core=NULL){
  x<-cbind(Age=dat$age,Male=dat$male,LogN=log(dat$N))
  if(!is.null(med_core))for(j in seq_len(ncol(x)))x[is.na(x[,j]),j]<-med_core[j]
  if(set=='l3'){
   stopifnot(is.finite(med_l3));x<-cbind(x,L3=ifelse(is.na(dat$l3),med_l3,dat$l3),L3_missing=as.integer(is.na(dat$l3)))
  }
  if(set=='categories'){
   dictionaries<-list(region=c('Asia','Europe','North America','South America','Oceania','Europe / Asia'),
    design=c('Retrospective cohort','Prospective cohort','Cross-sectional','Case-control'),
    modality_group=c('Endoscopy-based','Clinical assessment','Mixed / unspecified','Imaging-based','Pathology'))
   for(v in names(dictionaries))for(lev in dictionaries[[v]][-1])x<-cbind(x,setNames(data.frame(as.numeric(dat[[v]]==lev)),paste(v,lev,sep='_')))
  }
  as.matrix(x)
 }
 prep_fit<-function(train,set='core',w=family_w(train)){
  med<-if(set=='l3')median(train$l3,na.rm=TRUE)else 0
  med_core<-apply(raw_x(train),2,function(z){m<-median(z,na.rm=TRUE);if(is.finite(m))m else 0})
  x<-raw_x(train,set,med,med_core)
  mu<-colSums(x*w)/sum(w);ss<-sqrt(colSums(sweep(x,2,mu)^2*w)/sum(w));ss[ss<1e-10]<-1
  list(set=set,med_l3=med,med_core=med_core,mu=mu,sd=ss)
 }
 prep_apply<-function(dat,pp)sweep(sweep(raw_x(dat,pp$set,pp$med_l3,pp$med_core),2,pp$mu),2,pp$sd,'/')
 model_set<-function(name)if(name=='LASSO + L3')'l3'else if(name=='LASSO + reported categories')'categories'else'core'
 model_type<-function(name)if(grepl('^LASSO',name))'LASSO'else name
 fit_one<-function(name,train,param=NULL,seed=SEED,weights=NULL){
  w<-if(is.null(weights))family_w(train)else weights
  if(name=='Mean-only baseline')return(list(name=name,constant=weighted.mean(train$p,w)))
  pp<-prep_fit(train,model_set(name),w);x<-prep_apply(train,pp);type<-model_type(name)
  if(type=='LASSO'){
   fit<-glmnet::glmnet(x,train$y,alpha=1,lambda=param$lambda,weights=w,standardize=FALSE,intercept=TRUE)
  }else if(type=='Random forest'){
   fit<-ranger::ranger(x=as.data.frame(x),y=train$y,num.trees=500,mtry=param$mtry,
    min.node.size=param$min_node,case.weights=w,num.threads=1,seed=seed,write.forest=TRUE)
  }else if(type=='XGBoost'){
   fit<-xgboost::xgb.train(params=list(objective='reg:squarederror',tree_method='hist',max_depth=param$depth,
    eta=param$eta,min_child_weight=5,lambda=10,subsample=.8,colsample_bytree=1,nthread=1,seed=seed,
    base_score=weighted.mean(train$y,w)),data=xgboost::xgb.DMatrix(x,label=train$y,weight=w),nrounds=param$rounds,verbose=0)
  }else if(type=='REML meta-regression'){
   fit<-suppressWarnings(metafor::rma.uni(yi=train$y,vi=1/train$events+1/(train$N-train$events),
    mods=x,method='REML',control=list(maxiter=1000)))
  }else stop(name)
  list(name=name,type=type,fit=fit,pp=pp,param=param)
 }
 predict_logit<-function(fit,new){
  if(fit$name=='Mean-only baseline')return(rep(qlogis(fit$constant),nrow(new)))
  x<-prep_apply(new,fit$pp)
  switch(fit$type,LASSO=as.numeric(predict(fit$fit,x,s=fit$param$lambda)),
   'Random forest'=as.numeric(predict(fit$fit,data=as.data.frame(x),num.threads=1)$predictions),
   XGBoost=as.numeric(predict(fit$fit,x)),
   'REML meta-regression'=as.numeric(cbind(1,x)%*%coef(fit$fit)))
 }
 params_for<-function(name){
  switch(model_type(name),LASSO=lapply(lambda_grid,function(l)list(lambda=l)),
  'Random forest'=lapply(seq_len(6),function(i){g<-expand.grid(min_node=c(10,5,3),mtry=c(1,3));as.list(g[i,])}),
  XGBoost=lapply(seq_len(8),function(i){g<-expand.grid(depth=c(1,2),rounds=c(25,75),eta=c(.05,.1));as.list(g[i,])}),
  list(NULL))
 }
 score<-function(obs,pred,w){
  stopifnot(all(is.finite(pred)));mu<-weighted.mean(obs,w)
  c(RMSE=100*sqrt(weighted.mean((obs-pred)^2,w)),MAE=100*weighted.mean(abs(obs-pred),w),
  R2=1-sum(w*(obs-pred)^2)/sum(w*(obs-mu)^2))
 }
 tune<-function(name,train,seed=SEED){
  candidates<-params_for(name)
  if(length(candidates)==1)return(list(param=candidates[[1]],scores=NULL,index=1L))
  folds<-folds_grouped(train,5,seed);loss<-matrix(NA_real_,length(candidates),5)
  for(f in 1:5){
   tr<-train[folds!=f,,drop=FALSE];va<-train[folds==f,,drop=FALSE]
   stopifnot(length(intersect(tr$family,va$family))==0)
   for(j in seq_along(candidates)){
    fit<-fit_one(name,tr,candidates[[j]],seed+f*100+j)
    loss[j,f]<-score(va$p,plogis(predict_logit(fit,va)),family_w(va))['RMSE']
   }
  }
  means<-rowMeans(loss);ses<-apply(loss,1,sd)/sqrt(5);best<-which.min(means)
  eligible<-which(means<=means[best]+ses[best])
  # The grids above are ordered toward stronger penalties and larger leaves.
  # For XGBoost, choose the shallowest and then the fewest boosting rounds.
  if(model_type(name)=='XGBoost'){
   complexity<-sapply(candidates,function(x)x$depth*10000+x$rounds*10+x$eta)
   chosen<-eligible[which.min(complexity[eligible])]
  }else chosen<-eligible[1]
  list(param=candidates[[chosen]],index=chosen,scores=data.frame(candidate=seq_along(candidates),rmse=means,se=ses,selected=seq_along(candidates)==chosen))
 }
 predictions<-list();fold_audit<-list();tuning<-list();contributions<-list();fits_saved<-list();n<-0L
 start<-proc.time()[3]
 for(rep in seq_len(REPEATS)){
  outer<-folds_grouped(d,5,SEED+rep)
  for(f in 1:5){
   train<-d[outer!=f,,drop=FALSE];test<-d[outer==f,,drop=FALSE]
   stopifnot(!length(intersect(train$family,test$family)),!length(intersect(train$paper,test$paper)))
   fold_audit[[length(fold_audit)+1]]<-data.frame(repeat_id=rep,fold=f,row=test$row,paper=test$paper,family=test$family,
    training_cohorts=nrow(train),training_families=length(unique(train$family)),test_cohorts=nrow(test),test_families=length(unique(test$family)))
   for(name in models){
    seed<-SEED+1000*rep+100*f
    tu<-tune(name,train,seed);fit<-fit_one(name,train,tu$param,seed)
     fits_saved[[paste(name,rep,f,sep=':')]]<-fit
    lp<-predict_logit(fit,test);pr<-plogis(lp);n<-n+1L
    predictions[[n]]<-data.frame(model=name,repeat_id=rep,fold=f,row=test$row,paper=test$paper,family=test$family,
     observed=test$p,predicted=pr,predicted_logit=lp,weight=family_w(test))
    tuning[[n]]<-data.frame(model=name,repeat_id=rep,fold=f,parameters=jsonlite::toJSON(tu$param,auto_unbox=TRUE,null='null'),
     training_l3_median=if(name=='LASSO + L3')fit$pp$med_l3 else NA_real_)
    if(name=='LASSO'){
     x<-prep_apply(test,fit$pp);b<-as.numeric(coef(fit$fit,s=fit$param$lambda));phi<-sweep(x,2,b[-1],'*')
     stopifnot(max(abs(rowSums(phi)+b[1]-lp))<1e-7)
     for(j in seq_len(ncol(phi)))contributions[[length(contributions)+1]]<-data.frame(model='LASSO',repeat_id=rep,fold=f,row=test$row,
      feature=colnames(phi)[j],value=raw_x(test,med_core=fit$pp$med_core)[,j],contribution=phi[,j],baseline=b[1])
    }
    if(name=='XGBoost'){
     x<-prep_apply(test,fit$pp);sh<-predict(fit$fit,x,predcontrib=TRUE)
     stopifnot(is.matrix(sh),ncol(sh)==ncol(x)+1,max(abs(rowSums(sh)-lp))<1e-5)
     for(j in seq_len(ncol(x)))contributions[[length(contributions)+1]]<-data.frame(model='XGBoost',repeat_id=rep,fold=f,row=test$row,
      feature=colnames(x)[j],value=raw_x(test,med_core=fit$pp$med_core)[,j],contribution=sh[,j],baseline=sh[,ncol(sh)])
    }
   }
   cat(sprintf('Repeat %d/%d fold %d: %d train / %d test cohorts; elapsed %.1fs\n',rep,REPEATS,f,nrow(train),nrow(test),proc.time()[3]-start));flush.console()
  }
 }
 saveRDS(fits_saved,file.path(out,'outer_fold_fitted_models.rds'))
  pred<-do.call(rbind,predictions);fa<-do.call(rbind,fold_audit);tu<-do.call(rbind,tuning);ph<-do.call(rbind,contributions)
 write.csv(pred,file.path(out,'out_of_fold_predictions.csv'),row.names=FALSE,fileEncoding='UTF-8')
 write.csv(fa,file.path(out,'fold_membership.csv'),row.names=FALSE,fileEncoding='UTF-8')
 write.csv(tu,file.path(out,'selected_hyperparameters.csv'),row.names=FALSE,fileEncoding='UTF-8')
 write.csv(ph,file.path(out,'out_of_fold_SHAP.csv'),row.names=FALSE,fileEncoding='UTF-8')
 metrics<-do.call(rbind,lapply(models,function(name)do.call(rbind,lapply(seq_len(REPEATS),function(rep){
  z<-pred[pred$model==name&pred$repeat_id==rep,];sc<-score(z$observed,z$predicted,z$weight)
  data.frame(model=name,repeat_id=rep,RMSE=sc['RMSE'],MAE=sc['MAE'],R2=sc['R2'])}))))
 write.csv(metrics,file.path(out,'repeat_metrics.csv'),row.names=FALSE,fileEncoding='UTF-8')
 summary<-do.call(rbind,lapply(models,function(name){z<-metrics[metrics$model==name,];data.frame(model=name,
  RMSE_mean=mean(z$RMSE),RMSE_sd=sd(z$RMSE),MAE_mean=mean(z$MAE),MAE_sd=sd(z$MAE),R2_mean=mean(z$R2),R2_sd=sd(z$R2),
  RMSE_min=min(z$RMSE),RMSE_max=max(z$RMSE),R2_min=min(z$R2),R2_max=max(z$R2))}))
 write.csv(summary,file.path(out,'performance_summary.csv'),row.names=FALSE,fileEncoding='UTF-8')
 trainperf<-list()
 avg<-aggregate(cbind(predicted,observed)~model+row+paper+family,data=pred,mean)
 write.csv(avg,file.path(out,'mean_out_of_fold_predictions.csv'),row.names=FALSE,fileEncoding='UTF-8')
 cal<-do.call(rbind,lapply(models,function(name){
  z<-avg[avg$model==name,];w<-family_w(z);fit<-lm(observed~predicted,data=z,weights=w)
  data.frame(model=name,intercept_percentage_points=100*coef(fit)[1],slope=coef(fit)[2])
 }))
 write.csv(cal,file.path(out,'continuous_calibration.csv'),row.names=FALSE,fileEncoding='UTF-8')
 # Final all-data fits are for formulas and interpretation, not performance claims.
 full<-list();hyper<-list()
 for(name in models){t<-tune(name,d,SEED+99999);full[[name]]<-fit_one(name,d,t$param,SEED+99999);hyper[[name]]<-t$param}
 saveRDS(full,file.path(out,'final_models_for_interpretation.rds'))
 for(name in models){ff<-full[[name]];pr<-plogis(predict_logit(ff,d));sc<-score(d$p,pr,family_w(d));trainperf[[name]]<-data.frame(model=name,training_RMSE=sc['RMSE'],training_MAE=sc['MAE'],training_R2=sc['R2'])}
 write.csv(do.call(rbind,trainperf),file.path(out,'apparent_training_metrics.csv'),row.names=FALSE)
 writeLines(jsonlite::toJSON(hyper,auto_unbox=TRUE,pretty=TRUE,null='null'),file.path(out,'final_model_parameters.json'))
 preprocessing<-lapply(full,function(f){if(is.null(f$pp))return(NULL);list(feature_order=names(f$pp$mu),training_medians_core=as.list(f$pp$med_core),training_median_L3=f$pp$med_l3,training_weighted_means=as.list(f$pp$mu),training_weighted_SD=as.list(f$pp$sd),feature_set=f$pp$set)})
 writeLines(jsonlite::toJSON(preprocessing,pretty=TRUE,auto_unbox=TRUE,null='null'),file.path(out,'final_preprocessing.json'))
 lf<-full[['LASSO']];bb<-as.numeric(coef(lf$fit,s=lf$param$lambda))
 coeff<-data.frame(term=c('Intercept',core_names),standardized_coefficient=bb,training_mean=c(NA,lf$pp$mu),training_sd=c(NA,lf$pp$sd))
 write.csv(coeff,file.path(out,'lasso_prediction_formula.csv'),row.names=FALSE,fileEncoding='UTF-8')
 agegrid<-seq(min(d$age,na.rm=TRUE),max(d$age,na.rm=TRUE),length.out=80);malegrid<-seq(min(d$male,na.rm=TRUE),max(d$male,na.rm=TRUE),length.out=80)
 pdp<-function(fit,data,var,grid){sapply(grid,function(g){z<-data;z[[var]]<-g;weighted.mean(plogis(predict_logit(fit,z)),family_w(data))})}
 xf<-full[['XGBoost']]
 xgboost::xgb.save(xf$fit,file.path(out,'XGBoost_prediction_formula.json'))
 curves<-rbind(data.frame(feature='Reported age summary',x=agegrid,prediction=pdp(xf,d,'age',agegrid),lasso_prediction=pdp(lf,d,'age',agegrid)),
  data.frame(feature='Male proportion',x=malegrid,prediction=pdp(xf,d,'male',malegrid),lasso_prediction=pdp(lf,d,'male',malegrid)))
 # Conditional cluster bootstrap: fixed all-data penalty, resampling families.
 BOOT<-as.integer(Sys.getenv('UGI_PDP_BOOTSTRAPS','200'));set.seed(SEED+500000)
 curveboot<-matrix(NA_real_,BOOT,nrow(curves));families<-unique(d$family)
 for(bi in seq_len(BOOT)){
  counts<-table(factor(sample(families,length(families),replace=TRUE),levels=families));keep<-as.numeric(counts[d$family])>0
  bootdat<-d[keep,];weights<-as.numeric(counts[bootdat$family])*family_w(bootdat)
  bfit<-fit_one('XGBoost',bootdat,xf$param,SEED+bi,weights)
  curveboot[bi,]<-c(pdp(bfit,d,'age',agegrid),pdp(bfit,d,'male',malegrid))
 }
 curves$lower<-apply(curveboot,2,quantile,.025);curves$upper<-apply(curveboot,2,quantile,.975)
 write.csv(curves,file.path(out,'study_level_prediction_curves.csv'),row.names=FALSE,fileEncoding='UTF-8')
 writeLines(capture.output(sessionInfo()),file.path(out,'ML_sessionInfo.txt'))
 cat('\nActual repeated grouped five-fold results:\n');print(summary,row.names=FALSE)
 cat('\nPrimary LASSO formula coefficients:\n');print(coeff,row.names=FALSE)
 cat('No individual-level records, binary metrics, elderly extrapolation, or four-factor thresholds were generated.\n')

}
run_figures<-function(base,ml,out){
 require_pkgs(c('meta','metafor','ggplot2','patchwork','svglite'))
 suppressPackageStartupMessages({library(meta);library(ggplot2);library(patchwork);library(grid)})
# Original SCRIPT 1/6 forest templates; SCRIPT 5/9/10/8/11 composite templates.
# Numerical tables are updated; the original panel types and colours are retained.
d<-read.csv(file.path(base,'analysis_data.csv'),fileEncoding='UTF-8-BOM');d<-d[order(d$row),]
d$age_group<-ifelse(d$age_group=='Paediatric','Paediatric cohorts','Other / mixed cohorts')
d$region<-factor(d$region,levels=c('Asia','Europe','North America','South America','Oceania','Europe / Asia'))
d$modality_group<-factor(d$modality_group,levels=c('Endoscopy-based','Clinical assessment','Mixed / unspecified','Imaging-based','Pathology'))
fit<-function(z,tau='DL')metaprop(events,N,studlab=study,data=z,sm='PLOGIT',method='Inverse',method.tau=tau,
 method.random.ci='classic',method.ci='CP',common=FALSE,random=TRUE,prediction=TRUE,method.predict='V')
stat<-function(label,m){z<-data.frame(analysis=label,k=m$k,N=sum(m$n),events=sum(m$event),prevalence=100*plogis(m$TE.random),
 lower=100*plogis(m$lower.random),upper=100*plogis(m$upper.random),I2=100*m$I2,tau2=m$tau2,Q=m$Q)
 if(m$k==1){ci<-binom.test(m$event,m$n)$conf.int;z$lower<-100*ci[1];z$upper<-100*ci[2]};z}
m<-fit(d);age<-update(m,subgroup=d$age_group);reg<-update(m,subgroup=d$region);mod<-update(m,subgroup=d$modality_group)
group_stats<-function(v,mm){z<-do.call(rbind,lapply(unique(as.character(d[[v]])),function(g)stat(g,fit(d[as.character(d[[v]])==g,]))));tt<-metafor::rma.uni(yi=qlogis(d$events/d$N),vi=1/d$events+1/(d$N-d$events),mods=model.matrix(~factor(d[[v]]))[, -1,drop=FALSE],method='DL');z$p_between<-tt$QMp;z}
for(v in c('age_group','region','modality_group'))write.csv(group_stats(v,switch(v,age_group=age,region=reg,modality_group=mod)),file.path(ml,paste0(v,'_DL_results.csv')),row.names=FALSE)
d$design_group<-ifelse(grepl('Prospective',d$design),'Prospective',ifelse(grepl('Retrospective',d$design),'Retrospective','Other design'))
stopifnot(all(!is.na(d$nos)),d$nos[d$paper=='R065']==8)
d$quality_group<-ifelse(d$nos<7,'NOS <7','NOS ≥7')
for(v in c('design_group','quality_group'))write.csv(group_stats(v,NULL),file.path(ml,paste0(v,'_DL_results.csv')),row.names=FALSE)
models<-sensitivity_models(d,fit,m,base)
sr<-do.call(rbind,Map(stat,names(models),models));rownames(sr)<-NULL
standard_models<-models[c('Primary reported L4','Clear L4 definitions','Exclude case-control','Exclude Guerra','Objective location','Exclude clinical-labelled','One report per source family','Combine within each paper')]
ss<-do.call(rbind,Map(stat,names(standard_models),standard_models));rownames(ss)<-NULL
ld<-do.call(rbind,lapply(unique(d$paper),function(pp){mm<-fit(d[d$paper!=pp,]);data.frame(row=min(d$row[d$paper==pp]),paper=pp,study=paste(unique(d$study[d$paper==pp]),collapse=' / '),p=100*plogis(mm$TE.random),lower=100*plogis(mm$lower.random),upper=100*plogis(mm$upper.random))}))
eg<-metabias(m,method.bias='linreg',k.min=10);tf<-trimfill(m)
write.csv(sr,file.path(ml,'diagnostic_design_sensitivity.csv'),row.names=FALSE)
write.csv(ss,file.path(ml,'conventional_sensitivity.csv'),row.names=FALSE)
write.csv(rbind(stat('DL',m),stat('REML',fit(d,'REML'))),file.path(ml,'overall_results.csv'),row.names=FALSE)
write.csv(ld,file.path(ml,'leave_one_out_results.csv'),row.names=FALSE)
write.csv(data.frame(t=eg$statistic,df=eg$df,p=eg$p.value),file.path(ml,'egger_results.csv'),row.names=FALSE)
tf_summary<-stat('Trim and fill',tf);tf_summary$N<-NA_real_;tf_summary$events<-NA_real_
tf_summary$observed_N<-sum(d$N);tf_summary$observed_events<-sum(d$events)
tf_summary$note<-'k includes hypothetical imputed estimates; no additional patients or case counts are implied.'
write.csv(cbind(tf_summary,imputed=tf$k-m$k),file.path(ml,'trimfill_results.csv'),row.names=FALSE)
saveRDS(list(all=m,age=age,region=reg,modality=mod,sensitivity=models),file.path(ml,'classical_models.rds'))

TNR<-'Times New Roman';NAVY<-'#1D3557';BLUE<-'#00468B';RED<-'#ED0000'
original_theme<-theme_classic(base_size=11,base_family=TNR)+theme(
 plot.title=element_text(face='bold',size=11.5,color='black',hjust=.5,margin=margin(b=10)),
 axis.title=element_text(face='plain',size=10,color='black'),axis.text=element_text(face='plain',size=9,color='black'),
 plot.tag=element_text(face='bold',size=12,family=TNR),panel.grid=element_blank(),
 plot.margin=margin(10,16,12,12),legend.title=element_blank(),legend.text=element_text(size=9),legend.position='top')
theme_set(original_theme)
export<-function(p,name,w,h){
 # DDS: captions and figure titles belong in the manuscript, not inside artwork.
 if(inherits(p,'patchwork'))p<-p&theme(plot.title=element_blank())else p<-p+theme(plot.title=element_blank())
 ggsave(file.path(out,paste0(name,'.svg')),p,device=svglite::svglite,width=w,height=h,bg='white',limitsize=FALSE)
}
tags<-function(p)p+plot_annotation(tag_levels='a')&theme(plot.tag=element_text(face='bold',size=12,family=TNR))

# Original blue/red forest, same columns, RevMan5 subgroup template and line spacing.
forest_export<-function(mm,name,h,main=FALSE){
 svglite::svglite(file.path(out,paste0(name,'.svg')),width=if(main)11 else 8.5,height=h,pointsize=10,bg='white')
 aa<-list(x=mm,layout=if(main)'meta'else'RevMan5',prediction=FALSE,
 leftcols=c('studlab','event','n'),leftlabs=c('Study','Cases','Assessed CD'),
 rightcols=c('effect','ci','w.random'),rightlabs=c('Prevalence (%)','95% CI','Weight'),
 pscale=100,digits=2,smlab='Prevalence (%)',plotwidth=if(main)'4cm'else'6cm',
 subgroup.name='',print.byvar=FALSE,addrow=TRUE,addrow.subgroups=TRUE,addrow.overall=TRUE,addrows.below.overall=3,
 print.tau2=TRUE,print.I2=TRUE,print.pval.Q=TRUE,test.subgroup=FALSE,
 col.square=BLUE,col.square.lines=if(main)BLUE else'black',col.diamond=RED,col.diamond.lines=RED,
 col.diamond.random=RED,col.diamond.lines.random=RED,col.diamond.subgroup=RED,col.diamond.lines.subgroup=RED,
 fontfamily=TNR,fontsize=10,colgap='1.5mm',spacing=.9)
 if(!main){up<-100*max(plogis(c(mm$upper,mm$upper.random,mm$upper.random.w)),na.rm=TRUE);aa$xlim<-c(0,min(100,ceiling(up/5)*5))}
 if(main){aa$xlim<-c(0,55);aa$squaresize<-.6}
 do.call(forest,aa);dev.off()
}
allrows<-read.csv(file.path(base,'all_source_rows.csv'),fileEncoding='UTF-8-BOM');allrows<-allrows[order(allrows$row),]
allrows$study<-paste0(allrows$study,' [',allrows$scope_code,']')
mall<-metaprop(events,N,studlab=study,data=allrows,exclude=primary==0,sm='PLOGIT',method='Inverse',method.tau='DL',method.random.ci='classic',method.ci='CP',common=FALSE,random=TRUE)
stopifnot(length(mall$studlab)==55,abs(mall$TE.random-m$TE.random)<1e-12)
forest_export(mall,'Fig2',max(7,2.7+.23*nrow(allrows)),TRUE)
forest_export(age,'FigS1',max(8,4+.25*nrow(d)))
forest_export(reg,'FigS2',max(10,5+.25*nrow(d)))
forest_export(mod,'FigS3',max(10,5+.25*nrow(d)))


# Original Figure 3: a single four-colour quantile map; only separable locations.
# Multicountry totals are not allocated to each participating country.
location_key<-c('China'='China','Taiwan China'='Taiwan','Taiwan, China'='Taiwan','Denmark'='Denmark','Portugal'='Portugal','Brazil'='Brazil',
 'Italy'='Italy','UK'='United Kingdom','India'='India','Spain'='Spain','Turkey'='Turkey','France'='France',
 'Israel'='Israel','Sweden'='Sweden','Ireland'='Ireland','Chile'='Chile','Canada'='Canada','Croatia'='Croatia',
 'USA'='United States of America','Hungary'='Hungary','Netherlands'='Netherlands','New Zealand'='New Zealand')
d$map_location<-unname(location_key[trimws(d$country)])
loc<-do.call(rbind,lapply(sort(unique(na.omit(d$map_location))),function(g)stat(g,fit(d[d$map_location==g&!is.na(d$map_location),]))))
q<-as.numeric(quantile(loc$prevalence,c(0,.25,.5,.75,1)));stopifnot(length(unique(q))==5)
lev<-sprintf('%s percentile (%.2f%%-%.2f%%)',c('0-25th','26th-50th','51st-75th','76th-100th'),q[1:4],q[2:5])
loc$quartile<-as.character(cut(loc$prevalence,q,include.lowest=TRUE,labels=lev))
world<-read.csv(file.path(base,'world_polygons.csv'));world<-world[world$country!='Antarctica',]
world$quartile<-loc$quartile[match(world$country,loc$analysis)]
world$quartile[is.na(world$quartile)]<-'No separate estimate'
# The French cohorts were recruited in metropolitan France, not French Guiana.
world$quartile[world$country=='France'&world$lon< -20]<-'No separate estimate'
world$quartile<-factor(world$quartile,levels=c('No separate estimate',lev))
mapcols<-setNames(c('#EAEAEA','#C7E9B4','#7FCDBB','#41B6C4','#225EA8'),levels(world$quartile))
pmap<-ggplot(world,aes(lon,lat,group=group,fill=quartile))+geom_polygon(color='white',linewidth=.12)+
 coord_quickmap(xlim=c(-180,180),ylim=c(-56,84),expand=FALSE)+scale_fill_manual(values=mapcols,drop=FALSE,name=NULL)+
 theme_void(base_family=TNR)+theme(legend.position='inside',legend.position.inside=c(.005,.09),legend.justification=c(0,0),
 legend.background=element_rect(fill='white',color=NA),legend.text=element_text(size=9),legend.key.size=unit(4.5,'mm'),plot.margin=margin(5,5,5,5))
export(pmap,'Fig3',12,5.1)
write.csv(loc,file.path(ml,'map_single_location_DL_estimates.csv'),row.names=FALSE)
write.csv(d[is.na(d$map_location),c('row','study','country','N','events')],file.path(ml,'map_unallocated_multicountry_rows.csv'),row.names=FALSE)

# Original Figure 4: subgroup comparison and enhanced funnel, grey/navy colours.
# The unsupported "ML Corrected" estimates become a labelled DL/REML comparison.
subsets<-list('Overall'=d,'Paediatric cohorts'=d[d$age_group=='Paediatric cohorts',],'Other cohorts'=d[d$age_group!='Paediatric cohorts',],'Asia'=d[d$region=='Asia',])
df_forest<-do.call(rbind,lapply(names(subsets),function(g)do.call(rbind,lapply(c('DL','REML'),function(method){z<-stat(g,fit(subsets[[g]],method));z$Method<-method;z}))))
df_forest$analysis<-factor(df_forest$analysis,levels=rev(names(subsets)));df_forest$Method<-factor(df_forest$Method,levels=c('DL','REML'))
p4a<-ggplot(df_forest,aes(prevalence,analysis,color=Method))+
 geom_errorbar(aes(xmin=lower,xmax=upper),orientation='y',width=.2,linewidth=.8,position=position_dodge(width=.5))+
 geom_point(size=3.1,position=position_dodge(width=.5))+scale_color_manual(values=c('DL'='#999999','REML'=NAVY))+
 scale_x_continuous(limits=c(0,60),breaks=seq(0,60,10),labels=function(x)paste0(x,'%'))+
 labs(title='Subgroup Comparison',x='Reported L4 Prevalence (95% CI)',y=NULL)+original_theme+
 theme(panel.grid.major.x=element_line(color='#F0F0F0'),legend.position='top')
funnel_data<-data.frame(x=m$TE,se=m$seTE,N=d$N,
 modality=ifelse(grepl('endoscopy',d$modality_original,ignore.case=TRUE),'Endoscopy-labelled','Other / unspecified'))
funnel_data$modality<-factor(funnel_data$modality,levels=c('Other / unspecified','Endoscopy-labelled'))
se<-seq(0,max(m$seTE)*1.1,length.out=150)
band<-data.frame(se=se,lo=m$TE.random-1.96*se,hi=m$TE.random+1.96*se,lo99=m$TE.random-2.58*se,hi99=m$TE.random+2.58*se)
p4b<-ggplot()+geom_ribbon(data=band,aes(y=se,xmin=lo99,xmax=hi99),orientation='y',fill='#F2F2F2',alpha=.6)+
 geom_ribbon(data=band,aes(y=se,xmin=lo,xmax=hi),orientation='y',fill='#E5E5E5',alpha=.6)+
 geom_line(data=band,aes(lo,se),linetype=2,color=RED,linewidth=.45)+geom_line(data=band,aes(hi,se),linetype=2,color=RED,linewidth=.45)+
 geom_vline(xintercept=m$TE.random,linetype='dotted',color=NAVY,linewidth=.7)+
 geom_point(data=funnel_data,aes(x,se,size=N,color=modality),alpha=.8)+
 scale_color_manual(values=c('Other / unspecified'='#457B9D','Endoscopy-labelled'='#E63946'),name='Reported modality')+
 scale_size_continuous(name='Sample Size (CD)',range=c(1.5,5),breaks=c(1000,3000,5000))+
 scale_y_reverse(limits=c(max(se),0))+scale_x_continuous(breaks=qlogis(c(.01,.05,.1,.2,.4,.6,.8)),labels=c('1%','5%','10%','20%','40%','60%','80%'))+
 labs(title='Enhanced Funnel Plot',x='Reported L4 Prevalence (logit scale)',y='Standard Error of Logit Prevalence')+
 original_theme+theme(legend.position='right',legend.box='vertical',legend.title=element_text(size=9),legend.text=element_text(size=8))
export(tags(p4a+p4b+plot_layout(widths=c(1,1.2))),'Fig4',12,6.8)
write.csv(df_forest,file.path(ml,'Figure4_subgroup_DL_REML.csv'),row.names=FALSE)

# Original Figure 5: coloured sensitivity squares + the right validation panel.
# A continuous leave-one-out analysis replaces the invalid classification ROC.
senscols<-c('#8D99AE','#E9C46A','#2A9D8F','#457B9D','#DDA15E','#E76F51','#577590','#A2D2FF')
make_sens<-function(z,title){
 z$analysis<-factor(z$analysis,levels=rev(z$analysis));pal<-setNames(senscols[seq_len(nrow(z))],rev(levels(z$analysis)))
 text_x<-ceiling(max(z$upper,na.rm=TRUE)/5)*5+3;xmax<-text_x+19
 ggplot(z)+geom_vline(xintercept=100*plogis(m$TE.random),color='#E63946',linetype=2,linewidth=.8)+
 geom_errorbar(aes(y=as.numeric(analysis),xmin=lower,xmax=upper,color=analysis),orientation='y',width=.22,linewidth=.8)+
 geom_point(aes(y=as.numeric(analysis),x=prevalence,fill=analysis),shape=22,size=3.2,color='black',stroke=.55)+
 geom_text(aes(x=text_x,y=as.numeric(analysis),label=sprintf('%.1f%% (%.1f%%-%.1f%%)',prevalence,lower,upper)),hjust=0,size=3.0,family=TNR)+
 scale_color_manual(values=pal)+scale_fill_manual(values=pal)+scale_y_continuous(breaks=seq_len(nrow(z)),labels=levels(z$analysis))+
 scale_x_continuous(limits=c(0,xmax),breaks=seq(0,text_x-3,10),labels=function(x)paste0(x,'%'))+
 labs(title=title,x='Reported L4 Prevalence Rate (%)',y=NULL)+original_theme+theme(legend.position='none',
 panel.grid.major.x=element_line(color='#EAEAEA',linewidth=.4),plot.margin=margin(8,8,10,8))
}
p5a<-make_sens(ss,'Sensitivity Analysis')
p5b<-ggplot(ld,aes(p,row))+geom_vline(xintercept=100*plogis(m$TE.random),linetype=2,color='#777777',linewidth=.6)+
 geom_errorbar(aes(xmin=lower,xmax=upper),orientation='y',width=0,color='#A0BED7',linewidth=.35)+geom_point(color=NAVY,size=1.1)+
 labs(title='Leave-One-Paper-Out Analysis',x='Reported L4 Prevalence (95% CI) (%)',y='Omitted Table 1 Row')+original_theme
export(tags(p5a+p5b+plot_layout(widths=c(1.5,1))),'Fig5',13,5.5)


ps<-read.csv(file.path(ml,'performance_summary.csv'));av<-read.csv(file.path(ml,'mean_out_of_fold_predictions.csv'))
wfamily<-function(z)as.numeric(1/table(z$family)[z$family])
full<-readRDS(file.path(ml,'final_models_for_interpretation.rds'))
predict_xgb<-function(z){f<-full[['XGBoost']];xx<-cbind(Age=z$age,Male=z$male,LogN=log(z$N))
 for(j in seq_len(ncol(xx)))xx[is.na(xx[,j]),j]<-f$pp$med_core[j]
 xx<-sweep(sweep(xx,2,f$pp$mu),2,f$pp$sd,'/');plogis(as.numeric(predict(f$fit,xx)))}

# Original Figure 6: four blue bars and the grey continuous calibration curve.
four<-c('LASSO','Random forest','REML meta-regression','XGBoost');p6d<-ps[match(four,ps$model),]
p6d$model<-factor(p6d$model,levels=four)
fourcols<-setNames(c('#D1E5F0','#67A9CF','#3182BD','#08519C'),four)
baseline_r2<-ps$R2_mean[ps$model=='Mean-only baseline']
p6a<-ggplot(p6d,aes(model,R2_mean,fill=model))+geom_hline(yintercept=0,linewidth=.35)+
 geom_hline(yintercept=baseline_r2,color='#888888',linetype=2,linewidth=.5)+geom_col(width=.56,color='#444444',linewidth=.35)+
 geom_errorbar(aes(ymin=R2_mean-R2_sd,ymax=R2_mean+R2_sd),width=.15,linewidth=.4)+
 geom_text(aes(y=max(.15,max(p6d$R2_mean+p6d$R2_sd)+.08),label=sprintf('R² = %.3f\nRMSE = %.2f\nMAE = %.2f',R2_mean,RMSE_mean,MAE_mean)),size=3.0,family=TNR,lineheight=1.05)+
 scale_fill_manual(values=fourcols,guide='none')+scale_x_discrete(labels=c('LASSO','Random Forest','Meta-Regression','XGBoost'))+
 scale_y_continuous(limits=c(min(p6d$R2_mean-p6d$R2_sd)-.08,max(.3,max(p6d$R2_mean+p6d$R2_sd)+.15)),breaks=scales::breaks_pretty(5))+
 labs(title='Machine Learning Models Performance',x='Model Type',y='Predictive R²')+original_theme+
 theme(axis.text.x=element_text(angle=35,hjust=1))
cal<-av[av$model=='XGBoost',];cal$weight<-wfamily(cal);cal$obs<-100*cal$observed;cal$pred<-100*cal$predicted
cf<-lm(obs~pred,data=cal,weights=weight);cx<-seq(min(cal$pred),max(cal$pred),length.out=100)
cb<-matrix(NA_real_,200,length(cx));set.seed(20261027+602)
families<-unique(cal$family)
for(bi in 1:200){ct<-table(factor(sample(families,length(families),replace=TRUE),levels=families));ww<-cal$weight*as.numeric(ct[cal$family]);
 bf<-lm(obs~pred,data=cal,weights=ww);cb[bi,]<-predict(bf,newdata=data.frame(pred=cx))}
cd<-data.frame(pred=cx,observed=predict(cf,newdata=data.frame(pred=cx)),lower=apply(cb,2,quantile,.025),upper=apply(cb,2,quantile,.975))
p6b<-ggplot(cal,aes(pred,obs))+geom_abline(slope=1,intercept=0,linetype=2,color='black',linewidth=.8)+
 geom_ribbon(data=cd,aes(pred,ymin=lower,ymax=upper),inherit.aes=FALSE,fill='#DDDDDD',alpha=.5)+
 geom_point(color='#457B9D',alpha=.65,size=1.6)+geom_line(data=cd,aes(pred,observed),inherit.aes=FALSE,color='#666666',linewidth=1)+
 annotate('text',x=5,y=43,label=sprintf('Intercept = %.2f\nSlope = %.2f',coef(cf)[1],coef(cf)[2]),hjust=0,size=3.2,family=TNR)+
 coord_cartesian(xlim=c(0,55),ylim=c(0,55),expand=FALSE)+labs(title='XGBoost Model Calibration',x='Predicted Prevalence (%)',y='Observed Prevalence (%)')+original_theme
export(tags(p6a+p6b),'Fig6',12.6,7)
write.csv(cd,file.path(ml,'Figure6_calibration_curve.csv'),row.names=FALSE)

# Original Figure 7: signed blue/red absolute-beta bars and a real waterfall.
# This is an exploratory three-covariate REML model, not SHAP or a causal I2 partition.
rf<-full[['REML meta-regression']];regfit<-rf$fit
pretty_names<-c(Age='Reported Age Summary',Male='Male (%)',LogN='Log Sample Size')
coefdata<-data.frame(Covariate=unname(pretty_names),Beta=as.numeric(coef(regfit)[-1]))
coefdata$Effect<-abs(coefdata$Beta);coefdata$Impact<-ifelse(coefdata$Beta<0,'Negative','Positive')
coefdata$Covariate<-factor(coefdata$Covariate,levels=coefdata$Covariate[order(coefdata$Effect)])
p7a<-ggplot(coefdata,aes(Effect,Covariate,fill=Impact))+geom_col(color='black',linewidth=.3,width=.65)+
 scale_fill_manual(values=c('Negative'='#4289B9','Positive'='#DC5144'),
 labels=c('Negative'='Negative Association','Positive'='Positive Association'))+
 scale_x_continuous(expand=expansion(mult=c(0,.12)))+
 labs(title='Multivariate Standardized Effect Size',x='Standardized Effect Size\n(Absolute Beta)',y='Included Study Covariates')+
 original_theme+theme(legend.position='top',legend.direction='vertical',axis.title=element_text(face='bold',size=10))
rx<-cbind(Age=d$age,Male=d$male,LogN=log(d$N))
for(j in seq_len(ncol(rx)))rx[is.na(rx[,j]),j]<-rf$pp$med_core[j]
rx<-sweep(sweep(rx,2,rf$pp$mu),2,rf$pp$sd,'/')
rw<-wfamily(d);yi<-qlogis(d$events/d$N);vi<-1/d$events+1/(d$N-d$events)
null<-metafor::rma.uni(yi=yi,vi=vi,method='REML')
fits<-lapply(list(c('Age','Male'),c('Age','Male','LogN')),function(v)
 metafor::rma.uni(yi=yi,vi=vi,mods=rx[,v,drop=FALSE],method='REML'))
taus<-c(null$tau2,sapply(fits,function(f)f$tau2));pct<-100*taus/taus[1]
wf<-data.frame(Group=c('Demographics','Sample Size','Residual Heterogeneity'),y_pos=3:1,
 xmin=c(pmin(pct[-3],pct[-1]),0),xmax=c(pmax(pct[-3],pct[-1]),pct[3]),
 label=c(sprintf('%+.1f%%',pct[-3]-pct[-1]),sprintf('%.1f%%',pct[3])),
 fill_val=c('Red','LightRed','Grey'))
xhi<-max(100,wf$xmax)*1.2
p7b<-ggplot(wf)+geom_segment(aes(x=0,xend=xmin,y=y_pos,yend=y_pos),linetype=2,color='grey65')+
 geom_rect(aes(xmin=xmin,xmax=xmax,ymin=y_pos-.32,ymax=y_pos+.32,fill=fill_val))+
 geom_text(aes(x=(xmin+xmax)/2,y=y_pos,label=label),size=3.1,family=TNR)+
 geom_segment(x=100,xend=100,y=.5,yend=3.5,color='grey60',linewidth=1.3)+
 annotate('text',x=xhi-.5,y=3.68,label=sprintf('Baseline τ² = %.3f',taus[1]),hjust=1,size=3.4,family=TNR)+
 scale_fill_manual(values=c('Grey'='#808080','LightRed'='#F8B7B0','Red'='#E25D4F'),guide='none')+
 scale_y_continuous(breaks=1:3,labels=rev(wf$Group),limits=c(.5,3.9),expand=c(0,0))+
 scale_x_continuous(limits=c(0,xhi),expand=c(0,0))+
 labs(title='Heterogeneity Attribution',x='Sequential τ² Components\n(% of Baseline)',y='Covariate Groups')+original_theme+
 theme(axis.title=element_text(face='bold',size=10))
export(tags(p7a+p7b+plot_layout(widths=c(1,1.05))),'Fig7',13.2,7)
write.csv(coefdata,file.path(ml,'Figure7_standardized_REML_coefficients.csv'),row.names=FALSE)
write.csv(data.frame(step=c('Null','Age + Male','Add Log N'),tau2=taus,percent_of_baseline=pct),file.path(ml,'Figure7_sequential_tau2.csv'),row.names=FALSE)

# Original Figure 8: grey ICE curves + blue mean, and black age PDP on age zones.
# Sparse duration is replaced by the measured male proportion. No elderly extrapolation.
cu<-read.csv(file.path(ml,'study_level_prediction_curves.csv'))
malegrid<-cu$x[cu$feature=='Male proportion'];ice<-do.call(rbind,lapply(seq_len(nrow(d)),function(i){z<-d[rep(i,length(malegrid)),];z$male<-malegrid;
 data.frame(row=d$row[i],x=malegrid,prediction=100*predict_xgb(z))}))
male_pdp<-cu[cu$feature=='Male proportion',];age_pdp<-cu[cu$feature=='Reported age summary',]
p8a<-ggplot()+geom_line(data=ice,aes(x,prediction,group=row),color='grey80',alpha=.5,linewidth=.4)+
 geom_step(data=male_pdp,aes(x,100*prediction),color='#0C5F96',linewidth=1.2)+
 scale_y_continuous(expand=expansion(mult=c(.02,.1)))+
 labs(title='Male Proportion & UGI Involvement',x='Male Proportion (%)',y='Study-Level Predicted L4 Prevalence (%)')+original_theme
amin<-min(d$age,na.rm=TRUE);amax<-max(d$age,na.rm=TRUE)
p8b<-ggplot()+annotate('rect',xmin=amin,xmax=16,ymin=0,ymax=25,fill='#CCEBC5',alpha=.95)+
 annotate('rect',xmin=16,xmax=amax,ymin=0,ymax=25,fill='#B3CDE3',alpha=.95)+
 geom_vline(xintercept=16,linetype=2,color='grey30',linewidth=.8)+
 geom_step(data=age_pdp,aes(x,100*prediction),color='black',linewidth=1.2)+
 annotate('text',x=(amin+16)/2,y=22.5,label='≤16 Years',size=3.3,family=TNR)+
 annotate('text',x=(16+amax)/2,y=22.5,label='>16 Years',size=3.3,family=TNR)+
 scale_x_continuous(limits=c(amin,amax),breaks=c(10,16,25,35,45),expand=c(0,0))+
 scale_y_continuous(limits=c(0,25),breaks=seq(0,25,5),expand=c(0,0))+
 labs(title='Cohort Age & UGI Involvement',x='Reported Age Summary (Years)',y='Study-Level Predicted L4 Prevalence (%)')+original_theme
export(tags(p8a+p8b),'Fig8',12,7.1)
write.csv(ice,file.path(ml,'Figure8_actual_Male_ICE.csv'),row.names=FALSE)

# Original Figure 9: two tile maps above a horizontal effect panel.
# Unsupported patient phenotypes and triple interactions are replaced in-place
# by observed study-level cells and honest expanded-model validation.
hm<-d[!is.na(d$l3),];cuts<-as.numeric(quantile(hm$l3,c(0,1/3,2/3,1)))
hlabels<-sprintf('%s L3\n%s',c('Lower','Middle','Higher'),c(sprintf('≤%.1f%%',100*cuts[2]),sprintf('%.1f-%.1f%%',100*cuts[2],100*cuts[3]),sprintf('>%.1f%%',100*cuts[3])))
hm$l3_band<-cut(hm$l3,cuts,labels=hlabels,include.lowest=TRUE);hm$cohort_age<-ifelse(hm$age_group=='Paediatric cohorts','Paediatric cohorts','Other / mixed cohorts')
hp<-av[av$model=='LASSO + L3',];hm$pred<-hp$predicted[match(hm$row,hp$row)]
cells<-expand.grid(l3_band=hlabels,cohort_age=c('Other / mixed cohorts','Paediatric cohorts'),stringsAsFactors=FALSE)
cells$observed<-cells$predicted<-NA_real_;cells$k<-0L
for(i in seq_len(nrow(cells))){z<-hm[hm$l3_band==cells$l3_band[i]&hm$cohort_age==cells$cohort_age[i],];cells$k[i]<-nrow(z)
 if(nrow(z)){cells$observed[i]<-stat('',fit(z))$prevalence;cells$predicted[i]<-100*weighted.mean(z$pred,wfamily(z))}}
cells$l3_band<-factor(cells$l3_band,levels=hlabels);cells$cohort_age<-factor(cells$cohort_age,levels=c('Paediatric cohorts','Other / mixed cohorts'))
heatmap_panel<-function(v,title){z<-cells;z$value<-z[[v]];z$lab<-ifelse(is.na(z$value),'NR',sprintf('%.1f%%\n(k=%d)',z$value,z$k))
 ggplot(z,aes(l3_band,cohort_age,fill=value))+geom_tile(color='white',linewidth=1.2)+
 geom_text(aes(label=lab),size=3.6,family=TNR,lineheight=1.05)+
 scale_fill_gradient2(low='#B3CDE3',mid='white',high='#E56353',midpoint=20,limits=c(0,45),na.value='#EAEAEA',name='L4 Prev. (%)',guide=guide_colorbar(display='gradient',nbin=128))+
 labs(title=title,x='Study-Level L3 Proportion',y='Cohort Population')+theme_minimal(base_family=TNR,base_size=11)+
 theme(panel.grid=element_blank(),plot.title=element_text(face='bold',hjust=0,size=11.5,margin=margin(b=10)),
 axis.text=element_text(color='black',size=9),plot.margin=margin(10,14,10,10),plot.tag=element_text(face='bold',size=12,family=TNR),legend.title=element_text(size=9))}
p9a<-heatmap_panel('observed','Observed Pooled Prevalence');p9b<-heatmap_panel('predicted','Out-of-Fold Predicted Prevalence')
ext<-ps[match(c('LASSO','LASSO + L3','Mean-only baseline'),ps$model),];ext$model<-factor(ext$model,levels=rev(ext$model))
p9c<-ggplot(ext,aes(R2_mean,model,fill=model))+geom_vline(xintercept=0,color='#888888',linewidth=.4)+
 geom_col(width=.44)+geom_errorbar(aes(xmin=R2_min,xmax=R2_max),orientation='y',width=.1,linewidth=.6)+
 geom_text(aes(x=.13,label=sprintf('RMSE = %.2f; MAE = %.2f',RMSE_mean,MAE_mean)),hjust=0,size=3.0,family=TNR)+
 scale_fill_manual(values=setNames(c('#A9A9A9','#E86C59','#5C8FE4'),rev(levels(ext$model))),guide='none')+
 scale_x_continuous(limits=c(min(ext$R2_min)-.05,max(.43,max(ext$R2_max)+.1)),breaks=seq(-.3,.1,.1),labels=function(x)sprintf('%.1f',x))+
 labs(title='Expanded-Feature Validation',x='Predictive R² (Range Across Repeats)',y=NULL)+original_theme+
 theme(panel.grid.major.x=element_line(color='#EAEAEA',linetype=2),plot.title=element_text(hjust=0))
export(tags((p9a+p9b)/p9c+plot_layout(heights=c(1.5,1))),'Fig9',12.8,8.5)
write.csv(cells,file.path(ml,'Figure9_study_level_heatmap_cells.csv'),row.names=FALSE)
write.csv(data.frame(tertile_boundary=cuts),file.path(ml,'Figure9_L3_cutpoints.csv'),row.names=FALSE)

# Original supplemental Figure 4: exact blue/grey point and straight-red-line style.
pd<-data.frame(x=tf$TE,se=tf$seTE,imputed=tf$trimfill);se_max<-max(tf$seTE)*1.15;lim<-range(c(tf$TE,qlogis(c(.005,.8))))+c(-.2,.2);zcrit<-qnorm(.975)
svglite::svglite(file.path(out,'FigS5.svg'),width=8,height=7,pointsize=12,bg='white')
par(mar=c(5.8,5.2,3.8,1.3),family=TNR,xaxs='i',yaxs='i',mgp=c(3.1,.9,0),tcl=-.3)
plot(NA,xlim=lim,ylim=c(se_max,0),axes=FALSE,ann=FALSE)
ticks<-c(.5,1,2,5,10,20,40,60,80);axis(1,at=qlogis(ticks/100),labels=ticks,cex.axis=.9)
axis(2,at=seq(0,1,.2),labels=sprintf('%.1f',seq(0,1,.2)),cex.axis=.9);box(lwd=.8)
lines(c(tf$TE.random,tf$TE.random-zcrit*se_max),c(0,se_max),col='#EC0000',lty=2,lwd=1)
lines(c(tf$TE.random,tf$TE.random+zcrit*se_max),c(0,se_max),col='#EC0000',lty=2,lwd=1)
abline(v=m$TE.random,col='#00458B',lwd=1.5);abline(v=tf$TE.random,col='#EC0000',lty=2,lwd=1.5)
points(pd$x[!pd$imputed],pd$se[!pd$imputed],pch=19,col='#00458B',cex=1)
points(pd$x[pd$imputed],pd$se[pd$imputed],pch=1,col='#BEBEBE',cex=1,lwd=1.1)
title(xlab='Prevalence (%) on logit scale',ylab='Standard error of logit prevalence',
 main='',cex.lab=1.1,cex.main=1.12)
legend('topleft',inset=c(.015,.02),legend=c(sprintf('Observed estimates (n = %d)',m$k),sprintf('Imputed estimates (n = %d)',tf$k-m$k),'Original pooled estimate','Trim-fill adjusted estimate'),
 pch=c(19,1,NA,NA),lty=c(NA,NA,1,2),col=c('#00458B','#BEBEBE','#00458B','#EC0000'),lwd=c(NA,NA,1.5,1.5),bty='n',cex=.8,y.intersp=1.05)
dev.off()
# Original unadjusted blue/red funnel remains Supplemental Figure 4.
svglite::svglite(file.path(out,'FigS4.svg'),width=8,height=7,pointsize=12,bg='white')
par(mar=c(5.8,5.2,2,1.3),family=TNR,xaxs='i',yaxs='i',mgp=c(3.1,.9,0),tcl=-.3)
mx<-max(m$seTE)*1.15;ll<-range(c(m$TE,m$TE.random+c(-1,1)*1.96*mx))+c(-.2,.2)
plot(NA,xlim=ll,ylim=c(mx,0),axes=FALSE,ann=FALSE)
ticks<-c(.5,1,2,5,10,20,40,60,80);axis(1,at=qlogis(ticks/100),labels=ticks,cex.axis=.9)
axis(2,las=1);box(lwd=.8)
lines(m$TE.random-c(0,1.96*mx),c(0,mx),col=RED,lty=2)
lines(m$TE.random+c(0,1.96*mx),c(0,mx),col=RED,lty=2)
abline(v=m$TE.random,col=BLUE,lwd=1.5);points(m$TE,m$seTE,pch=19,col=BLUE)
title(xlab='Prevalence (%) on logit scale',ylab='Standard error of logit prevalence')
dev.off()
figures_10_11(d,av,ps,ml,out,export,tags,original_theme,TNR)
extra_source_figures(base,ml,out,fit,stat,forest_export,export,tags)
writeLines(capture.output(sessionInfo()),file.path(ml,'Figures_sessionInfo.txt'))
cat('Original figure numbers preserved: main 2–11; supplemental 1–12.\n')

}
finish_svgs<-function(figure_dir,output_root){
 editable<-file.path(output_root,'SVG_editable_text');dir.create(editable,showWarnings=FALSE)
 ff<-list.files(figure_dir,pattern='\\.svg$',full.names=TRUE)
 invisible(file.copy(ff,editable,overwrite=TRUE))
 if(Sys.getenv('UGI_FINISH_SVG','1')=='0')return(invisible(FALSE))
 find_runtime<-function(kind){
  override<-Sys.getenv(paste0('UGI_',toupper(kind)),'');if(nzchar(override)&&file.exists(override))return(override)
  suffix<-if(kind=='python')'python/python.exe'else'node/bin/node.exe'
  bundled<-file.path(Sys.getenv('USERPROFILE'),'.cache','codex-runtimes','codex-primary-runtime','dependencies',suffix)
  if(file.exists(bundled))return(bundled)
  unname(Sys.which(kind))
 }
 node<-find_runtime('node');python<-find_runtime('python')
 if(!nzchar(node)||!nzchar(python)){warning('Native R SVGs completed. Optional trimming/outlines need Node sharp and Python fontTools/lxml.');return(invisible(FALSE))}
 helper<-file.path(output_root,'svg_finish_helpers');dir.create(helper,showWarnings=FALSE)
js_source<-r"---(
const fs=require('fs'),path=require('path');
const candidates=[process.env.UGI_SHARP_MODULE,'sharp',path.join(process.env.USERPROFILE||'','.cache','codex-runtimes','codex-primary-runtime','dependencies','node','node_modules','sharp')].filter(Boolean);
let sharp;for(const c of candidates){try{sharp=require(c);break;}catch{}}
if(!sharp)throw Error('SVG trim requires Node package sharp; set UGI_SHARP_MODULE if needed.');
(async()=>{
 const dir=process.argv[2];
 for(const name of fs.readdirSync(dir).filter(x=>x.endsWith('.svg'))){
  const p=path.join(dir,name);let svg=fs.readFileSync(p,'utf8');
  const vb=svg.match(/viewBox=['"]([^'"]+)/)[1].split(/\s+/).map(Number);
  const {data,info}=await sharp(Buffer.from(svg),{density:144}).flatten({background:'#fff'}).removeAlpha().raw().toBuffer({resolveWithObject:true});
  let x0=info.width,y0=info.height,x1=0,y1=0;
  for(let y=0;y<info.height;y++)for(let x=0;x<info.width;x++){
   const i=(y*info.width+x)*3;if(Math.min(data[i],data[i+1],data[i+2])<235){x0=Math.min(x0,x);x1=Math.max(x1,x);y0=Math.min(y0,y);y1=Math.max(y1,y);}
  }
  if(x0<2||y0<2||x1>=info.width-2||y1>=info.height-2)throw Error('Content touches canvas '+name);
  const scale=vb[2]/info.width,pad=7;
  const v=[vb[0]+x0*scale-pad,vb[1]+y0*scale-pad,(x1-x0+1)*scale+2*pad,(y1-y0+1)*scale+2*pad];
  svg=svg.replace(/viewBox=['"][^'"]+['"]/,`viewBox="${v.join(' ')}"`).replace(/width=['"][^'"]+['"]/,`width="${v[2]}pt"`).replace(/height=['"][^'"]+['"]/,`height="${v[3]}pt"`);
  fs.writeFileSync(p,svg);console.log('Trimmed '+name);
 }
})().catch(e=>{process.stderr.write(e.stack);process.exit(1)});
)---"
py_source<-r"---(from pathlib import Path
import re,sys,os
from lxml import etree

from fontTools.ttLib import TTFont
from fontTools.pens.svgPathPen import SVGPathPen
B=Path(sys.argv[1])
NS='http://www.w3.org/2000/svg'
fonts={}
def props(el):
 return dict((k.strip(),v.strip())for k,v in [x.split(':',1)for x in el.get('style','').split(';')if ':'in x])
for p in B.glob('*.svg'):
 root=etree.parse(str(p)).getroot();count=0
 for el in list(root.iter()):
  tag=etree.QName(el).localname
  if tag in ['line','polyline','polygon','path','rect','circle']:
   defaults={'fill':'none','stroke':'#000000','stroke-linecap':'round','stroke-linejoin':'round','stroke-miterlimit':'10'}
   defaults.update(props(el))
   for k,v in defaults.items():el.set(k,v)
   el.attrib.pop('style',None)
 for el in list(root.xpath('//*[local-name()="text"]')):
  text=''.join(el.itertext());pr=props(el)
  if not text:continue
  assert not list(el),'Nested text unsupported'
  bold=pr.get('font-weight','') in ['bold','700'];italic=pr.get('font-style','')=='italic'
  name='timesbi.ttf' if bold and italic else 'timesbd.ttf' if bold else 'timesi.ttf' if italic else 'times.ttf'
  if name not in fonts:
   f=TTFont(str(Path(os.getenv('UGI_FONT_DIR',os.path.join(os.getenv('WINDIR','C:/Windows'),'Fonts')))/name));fonts[name]=(f,f.getGlyphSet(),f.getBestCmap())
  f,gs,cmap=fonts[name];units=f['head'].unitsPerEm
  fs=float(re.findall(r'[\d.]+',pr.get('font-size','9'))[0]);scale=fs/units
  glyphs=[cmap.get(ord(c),'.notdef')for c in text]
  widths=[f['hmtx'][g][0]for g in glyphs];advance=sum(widths)*scale
  target=el.get('textLength');target=float(re.findall(r'[\d.]+',target)[0])if target else advance
  sx=target/advance if advance else 1
  anchor=el.get('text-anchor',pr.get('text-anchor','start'));offset=-target if anchor=='end' else -target/2 if anchor=='middle'else 0
  x=float(el.get('x','0'))+offset;y=float(el.get('y','0'))
  group=etree.Element('{'+NS+'}g');group.set('data-font-family','Times New Roman');group.set('aria-label',text)
  if el.get('transform'):group.set('transform',el.get('transform'))
  inner=etree.SubElement(group,'{'+NS+'}g');inner.set('transform',f'translate({x},{y}) scale({scale*sx},{-scale})')
  inner.set('fill',pr.get('fill','#000000'));inner.set('stroke','none')
  advanceunits=0
  for g,w in zip(glyphs,widths):
   pen=SVGPathPen(gs);gs[g].draw(pen);commands=pen.getCommands()
   if commands:
    pa=etree.SubElement(inner,'{'+NS+'}path');pa.set('d',commands);pa.set('transform',f'translate({advanceunits},0)');pa.set('fill',pr.get('fill','#000000'));pa.set('stroke','none')
   advanceunits+=w
  el.getparent().replace(el,group);count+=1
 for st in root.xpath('//*[local-name()="style"]'):st.getparent().remove(st)
 p.write_bytes(etree.tostring(root,xml_declaration=True,encoding='UTF-8'))
 print(p.name,'outlined text',count)
)---"

 jsfile<-file.path(helper,'trim_svg.cjs');pyfile<-file.path(helper,'outline_svg.py')
 writeLines(js_source,jsfile,useBytes=TRUE);writeLines(py_source,pyfile,useBytes=TRUE)
 a<-system2(node,c(shQuote(jsfile),shQuote(normalizePath(figure_dir,winslash='/'))))
 if(a!=0){warning('SVG trimming failed; editable originals are in SVG_editable_text.');return(invisible(FALSE))}
 b<-system2(python,c(shQuote(pyfile),shQuote(normalizePath(figure_dir,winslash='/'))))
 if(b!=0){warning('Font outline conversion failed; editable originals remain available.');return(invisible(FALSE))}
 cat('SVG finishing complete: vectors retained, Times New Roman outlined.\n')
 invisible(TRUE)
}


if(stage%in%c('all','ml'))run_ml(data_dir,ml_dir)
if(stage%in%c('all','ml_sensitivity'))run_ml(data_dir,file.path(output_root,'ML_demographic_sensitivity'),TRUE)
if(stage%in%c('all','figures')){run_figures(data_dir,ml_dir,figure_dir);finish_svgs(figure_dir,output_root)}
if(stage%in%c('all','ml','validate'))audit_validation(data_dir,ml_dir)
cat('Completed source-verified stage:',stage,'\n')
