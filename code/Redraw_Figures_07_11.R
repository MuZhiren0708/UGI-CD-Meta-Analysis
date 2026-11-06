options(stringsAsFactors=FALSE)
if(.Platform$OS.type=='windows')suppressWarnings(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
extra<-Sys.getenv('UGI_R_LIBRARY','')
if(nzchar(extra)).libPaths(c(strsplit(extra,.Platform$path.sep,fixed=TRUE)[[1]],.libPaths()))
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(svglite)})
args<-commandArgs(trailingOnly=TRUE);base<-normalizePath(args[1],winslash='/',mustWork=TRUE)
dat<-file.path(base,'figure_data');out<-file.path(base,'SVG_figures');dir.create(out,showWarnings=FALSE,recursive=TRUE)
get<-function(n)read.csv(file.path(dat,n),check.names=FALSE,fileEncoding='UTF-8-BOM')
font<-'Times New Roman';navy<-'#244B63';teal<-'#238B8B';rust<-'#C66B4E';grey<-'#8797A1';light<-'#D9E6EB'
pal<-c('Group 1'=rust,'Group 2'=teal)
base_theme<-theme_classic(base_size=11,base_family=font)+theme(
 text=element_text(color='#25313A'),axis.text=element_text(color='#34434D',size=10),
 axis.title=element_text(size=11),axis.line=element_line(color='#697983',linewidth=.35),
 axis.ticks=element_line(color='#697983',linewidth=.3),axis.ticks.length=grid::unit(2,'mm'),
 plot.title=element_text(face='bold',size=12,hjust=0,margin=margin(b=11)),
 plot.tag=element_text(face='bold',size=13,family=font),plot.tag.position=c(0,.99),
 plot.margin=margin(12,16,10,12),legend.position='top',legend.title=element_blank(),
 legend.text=element_text(size=10),legend.key.width=grid::unit(13,'mm'),
 panel.grid.major.y=element_line(color='#EEF1F3',linewidth=.25))
theme_set(base_theme)
tag<-function(p,l)p+labs(tag=l)
save<-function(p,name,w,h){
 svglite(file.path(out,paste0(name,'.svg')),width=w,height=h,bg='white',pointsize=11)
 print(p);dev.off()
}

# Figure 7: preserve absolute beta values and sequential tau-squared components.
cf<-get('Figure7_standardized_REML_coefficients.csv');tt<-get('Figure7_sequential_tau2.csv')
cf$label<-c('Reported age','Male proportion','Log sample size');cf$y<-3:1
p7a<-ggplot(cf,aes(Effect,y,color=Impact))+geom_segment(aes(x=0,xend=Effect,yend=y),linewidth=1.4,alpha=.30)+
 geom_point(size=3.8)+geom_text(aes(label=sprintf('%.3f',Effect)),nudge_x=.025,hjust=0,size=3.6,family=font,show.legend=FALSE)+
 scale_color_manual(values=c(Negative=navy,Positive=rust),labels=c('Negative association','Positive association'))+
 scale_y_continuous(breaks=1:3,labels=rev(cf$label),limits=c(.5,3.5))+
 scale_x_continuous(limits=c(0,.48),breaks=c(0,.1,.2,.3,.4),expand=c(0,0))+
 labs(title='Standardised coefficients',x=expression('Absolute standardised coefficient ('*abs(beta)*')'),y=NULL)+
 theme(legend.position='bottom',panel.grid.major.y=element_blank())
pc<-tt$percent_of_baseline
wf<-data.frame(y=3:1,xmin=c(pc[2],pc[2],0),xmax=c(100,pc[3],pc[3]),type=c('reduction','increase','residual'),
 label=c(sprintf('%.1f%% lower',100-pc[2]),sprintf('+%.1f%% of null',pc[3]-pc[2]),sprintf('%.1f%% remaining',pc[3])))
p7b<-ggplot(wf)+geom_vline(xintercept=100,color='#B4BEC4',linetype=2,linewidth=.45)+
 geom_rect(aes(xmin=xmin,xmax=xmax,ymin=y-.24,ymax=y+.24,fill=type),color=NA)+
 geom_text(aes(x=105,y=y,label=label),hjust=0,size=3.35,family=font)+
 scale_fill_manual(values=c(reduction=teal,increase=rust,residual=navy),guide='none')+
 scale_y_continuous(breaks=1:3,labels=c('Residual','Add sample size','Age and sex'),limits=c(.5,3.5))+
 scale_x_continuous(limits=c(0,165),breaks=c(0,25,50,75,100),expand=c(0,0))+
 labs(title='Residual heterogeneity',x=expression(tau^2*' relative to null model (%)'),y=NULL)+
 annotate('text',x=0,y=.56,label=sprintf('Null model: τ² = %.3f',tt$tau2[1]),hjust=0,size=3.3,family=font,color=grey)+
 theme(panel.grid.major.y=element_blank())
save(tag(p7a,'a')+tag(p7b,'b')+plot_layout(widths=c(1,1.12)),'Fig7',10.8,4.7)

# Figure 8: true saved predictions, no smoothing or extrapolation.
cu<-get('study_level_prediction_curves.csv');ice<-get('Figure8_actual_Male_ICE.csv');d<-get('cohort_plot_data.csv')
male<-subset(cu,feature=='Male proportion');age<-subset(cu,feature=='Reported age summary')
line_theme<-theme(panel.grid.major.y=element_line(color='#E9EEF1',linewidth=.3),legend.position='none')
p8a<-ggplot()+geom_step(data=ice,aes(x,prediction,group=row),color='#B5CCD6',alpha=.48,linewidth=.35)+
 geom_step(data=male,aes(x,100*prediction),color=navy,linewidth=1.05)+
 geom_rug(data=d[!is.na(d$male),],aes(x=male),inherit.aes=FALSE,sides='b',color=grey,alpha=.6,linewidth=.3)+
 scale_y_continuous(limits=c(8.5,16),breaks=seq(9,15,2))+
 labs(title='Male proportion',x='Male proportion (%)',y='Predicted study-level L4 prevalence (%)')+line_theme
p8b<-ggplot(age,aes(x,100*prediction))+geom_step(color=teal,linewidth=1.05)+
 geom_vline(xintercept=16,linetype=3,color='#A7B2B9',linewidth=.45)+
 geom_rug(data=d[!is.na(d$age),],aes(x=age),inherit.aes=FALSE,sides='b',color=grey,alpha=.6,linewidth=.3)+
 annotate('text',x=16.5,y=15.6,label='16 years',hjust=0,size=3.2,family=font,color=grey)+
 scale_y_continuous(limits=c(8.5,16),breaks=seq(9,15,2))+scale_x_continuous(breaks=c(10,16,25,35,45))+
 labs(title='Reported age',x='Reported age summary (years)',y='Predicted study-level L4 prevalence (%)')+line_theme
save(tag(p8a,'a')+tag(p8b,'b'),'Fig8',10.8,4.6)

# Figure 9: common colour scale and the original validation ranges (not CIs).
cells<-get('Figure9_study_level_heatmap_cells.csv');ps<-get('performance_summary.csv')
cells$band<-factor(rep(1:3,2),levels=1:3,labels=c('Lower L3\n≤35.0%','Middle L3\n35.0–48.1%','Higher L3\n>48.1%'))
cells$age<-factor(cells$cohort_age,levels=c('Paediatric cohorts','Other / mixed cohorts'),labels=c('Paediatric','Other / mixed'))
heat<-function(v,title){z<-cells;z$value<-z[[v]]
 ggplot(z,aes(band,age,fill=value))+geom_tile(width=.96,height=.94,color='white',linewidth=1)+
 geom_text(aes(label=sprintf('%.1f%%\n(k = %d)',value,k),color=value>27),family=font,size=3.8,lineheight=1.15)+
 scale_color_manual(values=c('FALSE'='#24404B','TRUE'='white'),guide='none')+
 scale_fill_gradientn(colours=c('#F3F7F7','#B6DAD6','#56A9A6','#176B72'),limits=c(0,45),breaks=c(0,15,30,45),name='L4 prevalence (%)',guide=guide_colorbar(display='rectangles',nbin=128))+
 labs(title=title,x=NULL,y=NULL)+coord_fixed(ratio=.9)+
 theme(axis.line=element_blank(),axis.ticks=element_blank(),panel.grid=element_blank(),legend.position='bottom',
 legend.key.width=grid::unit(15,'mm'),legend.title=element_text(size=10),legend.title.position='top',plot.margin=margin(8,10,0,10))}
ext<-ps[match(c('LASSO','LASSO + L3','Mean-only baseline'),ps$model),];ext$y<-3:1
ext$display<-c('Core LASSO','LASSO + L3','Mean-only baseline')
p9c<-ggplot(ext,aes(R2_mean,y,color=model))+geom_vline(xintercept=0,color=grey,linetype=2,linewidth=.4)+
 geom_segment(aes(x=R2_min,xend=R2_max,yend=y),linewidth=.9)+geom_point(size=3.2)+
 scale_color_manual(values=setNames(c(navy,rust,grey),ext$model),guide='none')+
 scale_y_continuous(breaks=1:3,labels=rev(ext$display),limits=c(.5,3.65))+
 scale_x_continuous(limits=c(-.29,.065),breaks=c(-.25,-.15,-.05,0,.05))+
 labs(title='Validation with and without L3',x='Predictive R² (range across repeats)',y=NULL)
mt<-rbind(data.frame(y=ext$y,x=1,label=sprintf('%.2f',ext$RMSE_mean)),data.frame(y=ext$y,x=2,label=sprintf('%.2f',ext$MAE_mean)))
p9t<-ggplot(mt,aes(x,y,label=label))+geom_text(size=3.7,family=font)+
 annotate('text',x=c(1,2),y=3.55,label=c('RMSE','MAE'),fontface='bold',size=3.6,family=font)+
 coord_cartesian(xlim=c(.5,2.5),ylim=c(.5,3.65))+theme_void(base_family=font)+
 theme(plot.margin=margin(38,12,35,0))
upper<-(tag(heat('observed','Observed prevalence'),'a')+tag(heat('predicted','Out-of-fold prediction'),'b'))+plot_layout(guides='collect')&theme(legend.position='bottom')
lower<-tag(p9c,'c')+p9t+plot_layout(widths=c(3.3,1))
save(upper/lower+plot_layout(heights=c(1.5,1)),'Fig9',10.8,7.2)

# Figure 10: retain quartile strata and every observed cohort point.
a<-get('Figure10_OOF_prediction_strata.csv');sm<-get('Figure10_stratum_summaries.csv')
long<-rbind(data.frame(stratum=sm$stratum,value=sm$observed,type='Observed'),data.frame(stratum=sm$stratum,value=sm$predicted,type='Out-of-fold prediction'))
p10a<-ggplot(long,aes(stratum,value,color=type,shape=type,group=type))+geom_line(linewidth=.8)+geom_point(size=3.3)+
 scale_color_manual(values=c(Observed=rust,'Out-of-fold prediction'=navy))+
 scale_shape_manual(values=c(Observed=16,'Out-of-fold prediction'=17))+
 scale_x_discrete(labels=setNames(sprintf('%s\n(k = %d)',sm$stratum,sm$k),sm$stratum))+
 labs(title='Observed and predicted means',x='Out-of-fold prediction stratum',y='Mean study prevalence (%)')+
 theme(legend.position='bottom',panel.grid.major.y=element_line(color='#EEF1F3',linewidth=.3))
p10b<-ggplot(a,aes(stratum,observed))+geom_boxplot(width=.48,fill='#D5E7E5',color=navy,linewidth=.55,outlier.shape=NA)+
 geom_point(position=position_jitter(width=.10,seed=20261101),size=2,alpha=.8,color=teal)+
 scale_x_discrete(labels=setNames(sprintf('%s\n(k = %d)',sm$stratum,sm$k),sm$stratum))+
 labs(title='Observed cohort distributions',x='Out-of-fold prediction stratum',y='Observed study prevalence (%)')+
 theme(panel.grid.major.y=element_line(color='#EEF1F3',linewidth=.3))
save(tag(p10a,'a')+tag(p10b,'b'),'Fig10',10.8,4.9)

# Figure 11: retain the saved PCA coordinates, group memberships and profiles.
zz<-get('Figure11_covariate_clusters.csv');radar<-get('Figure11_scaled_cluster_profiles.csv')
mod<-readRDS(file.path(dat,'Figure11_descriptive_cluster_models.rds'));pv<-100*mod$pca$sdev^2/sum(mod$pca$sdev^2)
hull<-do.call(rbind,lapply(split(zz,zz$group),function(v)v[chull(v$PC1,v$PC2),]))
gl<-setNames(sprintf('%s (k = %d)',names(table(zz$group)),as.integer(table(zz$group))),names(table(zz$group)))
p11a<-ggplot(zz,aes(PC1,PC2,color=group,fill=group))+geom_hline(yintercept=0,color='#E3E8EB',linewidth=.35)+
 geom_vline(xintercept=0,color='#E3E8EB',linewidth=.35)+geom_polygon(data=hull,alpha=.075,linewidth=.5)+
 geom_point(size=2.8,alpha=.9,shape=21,stroke=.4)+
 scale_color_manual(values=pal,labels=gl)+scale_fill_manual(values=pal,labels=gl)+
 labs(title='Cohort groups',x=sprintf('PC1 (%.1f%%)',pv[1]),y=sprintf('PC2 (%.1f%%)',pv[2]))+
 theme(panel.grid=element_blank(),legend.position='bottom')
grid<-do.call(rbind,lapply(c(.25,.5,.75,1),function(v)data.frame(axis=c(1:5,1),value=v,level=v)))
grid$angle<-(grid$axis-1)*2*pi/5;grid$x<-grid$value*sin(grid$angle);grid$y<-grid$value*cos(grid$angle)
ax<-data.frame(angle=(0:4)*2*pi/5,label=c('Reported\nage','Male\nproportion','Log sample\nsize','Observed\nL4*','Predicted\nL4*'))
ax$x<-1.23*sin(ax$angle);ax$y<-1.23*cos(ax$angle)
p11b<-ggplot()+geom_polygon(data=grid,aes(x,y,group=level),fill=NA,color='#DDE4E7',linewidth=.4)+
 geom_segment(data=ax,aes(x=0,y=0,xend=x/1.23,yend=y/1.23),color='#E8EDEF',linewidth=.35)+
 geom_polygon(data=radar,aes(x,y,group=group,fill=group,color=group),alpha=.075,linewidth=.8)+
 geom_point(data=subset(radar,axis<=5),aes(x,y,color=group),size=2)+
 geom_text(data=ax,aes(x,y,label=label),size=3.7,family=font,lineheight=1.05)+
 annotate('text',x=.05,y=c(.25,.5,.75,1),label=c('.25','.50','.75','1.00'),size=2.9,hjust=0,color=grey,family=font)+
 scale_color_manual(values=pal,guide='none')+scale_fill_manual(values=pal,guide='none')+
 coord_fixed(xlim=c(-1.55,1.55),ylim=c(-1.4,1.5),clip='off')+
 labs(title='Scaled cohort profiles')+theme_void(base_family=font)+
 theme(text=element_text(color='#25313A'),plot.title=element_text(face='bold',size=12,hjust=0,margin=margin(b=10,l=24)),
 plot.tag=element_text(face='bold',size=13,family=font),plot.tag.position=c(0,.99),plot.margin=margin(12,20,10,20))
save(tag(p11a,'a')+tag(p11b,'b')+plot_layout(widths=c(1.06,1)),'Fig11',10.8,5.0)
capture.output(sessionInfo(),file=file.path(base,'R_sessionInfo.txt'))
cat('Five SVG figures generated from saved analysis results.\n')
