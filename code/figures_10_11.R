# Original curve/boxplot and PCA/radar layouts; actual study data only.
figures_10_11<-function(d,av,ps,ml,out,export,tags,original_theme,TNR){
 a<-av[av$model=='XGBoost',];a<-a[match(d$row,a$row),];a$observed<-100*a$observed;a$predicted<-100*a$predicted
 breaks<-unique(as.numeric(quantile(a$predicted,seq(0,1,length.out=5))))
 a$stratum<-cut(a$predicted,breaks,include.lowest=TRUE,labels=paste0('Q',seq_len(length(breaks)-1)))
 sm<-do.call(rbind,lapply(levels(a$stratum),function(g){z<-a[a$stratum==g,];w<-1/as.numeric(table(z$family)[z$family]);data.frame(stratum=g,k=nrow(z),observed=weighted.mean(z$observed,w),predicted=weighted.mean(z$predicted,w))}))
 long<-rbind(data.frame(sm[c('stratum','k')],value=sm$observed,type='Observed'),data.frame(sm[c('stratum','k')],value=sm$predicted,type='Out-of-fold prediction'))
 p10a<-ggplot2::ggplot(long,ggplot2::aes(stratum,value,color=type,group=type))+ggplot2::geom_line(linewidth=1)+ggplot2::geom_point(size=3)+
  ggplot2::scale_color_manual(values=c('Observed'='#E76F51','Out-of-fold prediction'='#264653'))+
  ggplot2::labs(x='Out-of-Fold Prediction Stratum',y='Mean Study Prevalence (%)')+original_theme
 p10b<-ggplot2::ggplot(a,ggplot2::aes(stratum,observed,fill=stratum))+ggplot2::geom_boxplot(width=.55,outlier.shape=NA,color='#444444')+
  ggplot2::geom_point(position=ggplot2::position_jitter(width=.1,seed=20261028),size=1.7,alpha=.7,color='#264653')+
  ggplot2::scale_fill_manual(values=c('#EAEAEA','#A8DADC','#8AB8C6','#E76F51')[seq_along(levels(a$stratum))],guide='none')+
  ggplot2::labs(x='Out-of-Fold Prediction Stratum',y='Observed Study Prevalence (%)')+original_theme
 export(tags(p10a+p10b),'Fig10',12,6.4)
 write.csv(a,file.path(ml,'Figure10_OOF_prediction_strata.csv'),row.names=FALSE)
 write.csv(sm,file.path(ml,'Figure10_stratum_summaries.csv'),row.names=FALSE)
 write.csv(data.frame(boundary=breaks),file.path(ml,'Figure10_prediction_cutpoints.csv'),row.names=FALSE)

 x<-cbind(Age=d$age,Male=d$male,LogN=log(d$N));med<-apply(x,2,median,na.rm=TRUE)
 for(j in seq_len(ncol(x)))x[is.na(x[,j]),j]<-med[j]
 xx<-scale(x);set.seed(20261028);km<-lapply(2:5,function(k)kmeans(xx,k,nstart=100));ds<-dist(xx)
 si<-sapply(km,function(z)mean(cluster::silhouette(z$cluster,ds)[,3]));best<-which.max(si);z<-km[[best]]
 ord<-order(z$centers[,'Age']);lab<-setNames(seq_along(ord),ord);grp<-as.integer(lab[as.character(z$cluster)])
 pc<-prcomp(xx,center=FALSE,scale.=FALSE);pv<-100*pc$sdev^2/sum(pc$sdev^2)
 zz<-data.frame(row=d$row,paper=d$paper,family=d$family,study=d$study,group=paste0('Group ',grp),PC1=pc$x[,1],PC2=pc$x[,2])
 hull<-do.call(rbind,lapply(split(zz,zz$group),function(v)if(nrow(v)>=3)v[chull(v$PC1,v$PC2),]else v[FALSE,]))
 pal<-setNames(c('#E76F51','#2A9D8F','#457B9D','#E9C46A','#9B5DE5')[seq_along(ord)],paste0('Group ',seq_along(ord)))
 p11a<-ggplot2::ggplot(zz,ggplot2::aes(PC1,PC2,color=group,fill=group))+
  ggplot2::geom_polygon(data=hull,ggplot2::aes(group=group),alpha=.12,linewidth=.4)+ggplot2::geom_point(size=2.4,alpha=.85)+
  ggplot2::scale_color_manual(values=pal)+ggplot2::scale_fill_manual(values=pal)+
  ggplot2::labs(x=sprintf('PC1 (%.1f%%)',pv[1]),y=sprintf('PC2 (%.1f%%)',pv[2]))+original_theme
 profile<-cbind(x,Observed=100*d$events/d$N,Predicted=a$predicted)
 scaled<-apply(profile,2,function(v)if(diff(range(v))>0)(v-min(v))/diff(range(v))else rep(0,length(v)))
 radar<-do.call(rbind,lapply(seq_along(ord),function(g){val<-colMeans(scaled[grp==g,,drop=FALSE]);data.frame(group=paste0('Group ',g),axis=c(seq_along(val),1),value=c(val,val[1]))}))
 axisnames<-c('Reported age','Male (%)','Log sample size','Observed L4*','OOF predicted L4*')
 grid<-do.call(rbind,lapply(c(.25,.5,.75,1),function(v)data.frame(axis=c(1:5,1),value=v,level=v)))
 radar$angle<-(radar$axis-1)*2*pi/5;radar$x<-radar$value*sin(radar$angle);radar$y<-radar$value*cos(radar$angle)
 grid$angle<-(grid$axis-1)*2*pi/5;grid$x<-grid$value*sin(grid$angle);grid$y<-grid$value*cos(grid$angle)
 labs<-data.frame(x=1.23*sin((0:4)*2*pi/5),y=1.23*cos((0:4)*2*pi/5),label=c('Reported age','Male (%)','Log sample size','Observed L4*','OOF predicted L4*'))
 p11b<-ggplot2::ggplot()+ggplot2::geom_polygon(data=grid,ggplot2::aes(x,y,group=level),fill=NA,color='#DDDDDD',linewidth=.35)+
  ggplot2::geom_segment(data=labs,ggplot2::aes(x=0,y=0,xend=x/1.23,yend=y/1.23),color='#EEEEEE',linewidth=.3)+
  ggplot2::geom_polygon(data=radar,ggplot2::aes(x,y,color=group,fill=group,group=group),alpha=.08,linewidth=.6)+
  ggplot2::geom_point(data=radar,ggplot2::aes(x,y,color=group),size=1.4)+
  ggplot2::geom_text(data=labs,ggplot2::aes(x,y,label=label),size=3,family=TNR)+
  ggplot2::scale_color_manual(values=pal)+ggplot2::scale_fill_manual(values=pal)+
  ggplot2::coord_fixed(xlim=c(-1.65,1.65),ylim=c(-1.4,1.4),expand=FALSE)+
  ggplot2::theme_void(base_family=TNR)+ggplot2::theme(plot.margin=ggplot2::margin(20,22,20,22),legend.position='none')
 export(tags(p11a+p11b+patchwork::plot_layout(widths=c(1.15,1))),'Fig11',12.8,6.8)
 write.csv(zz,file.path(ml,'Figure11_covariate_clusters.csv'),row.names=FALSE)
 write.csv(data.frame(k=2:5,silhouette=si,selected=(2:5)==best+1),file.path(ml,'Figure11_cluster_selection.csv'),row.names=FALSE)
 write.csv(pc$rotation,file.path(ml,'Figure11_PCA_loadings.csv'))
 write.csv(radar,file.path(ml,'Figure11_scaled_cluster_profiles.csv'),row.names=FALSE)
 write.csv(data.frame(feature=names(med),imputation_median=med),file.path(ml,'Figure11_descriptive_preprocessing.csv'),row.names=FALSE)
 saveRDS(list(kmeans=z,pca=pc,medians=med,scale=attributes(xx)),file.path(ml,'Figure11_descriptive_cluster_models.rds'))
}
