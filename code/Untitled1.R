# Untitled1.R — based on the supplied Untitled1.R plotting blocks.
# Original marginal scatter/residual and dual-axis SHAP appearance retained.
# Real OOF predictions and three source-aligned predictors replace the simulated inputs.
# First run UGI_CD_Meta_Analysis.R, then:
# Rscript Untitled1.R [the_same_output_folder]
options(stringsAsFactors=FALSE,scipen=99)
if(.Platform$OS.type=='windows')suppressWarnings(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
extra<-Sys.getenv('UGI_R_LIBRARY','')
if(nzchar(extra)).libPaths(c(strsplit(extra,.Platform$path.sep,fixed=TRUE)[[1]],.libPaths()))
for(p in c('ggplot2','patchwork','svglite'))if(!requireNamespace(p,quietly=TRUE))stop('Install required package: ',p)
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(grid)})
args<-commandArgs(trailingOnly=TRUE)
script<-sub('^--file=','',commandArgs()[grepl('^--file=',commandArgs())])
here<-if(length(script))dirname(normalizePath(script,winslash='/'))else getwd()
output_root<-if(length(args))args[1]else file.path(here,'DDS_results')
base<-file.path(output_root,'data');ml<-file.path(output_root,'analysis_results_final')
out<-file.path(output_root,'SVG_extra');dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(!file.exists(file.path(ml,'out_of_fold_SHAP.csv')))stop('Run UGI_CD_Meta_Analysis.R all first; use the same output directory.')
# Untitled1.R: keep the original marginal-scatter/residual and dual-axis SHAP templates.
# Replace actual-plus-random-noise predictions and simulated SHAP with real OOF values.
d<-read.csv(file.path(base,'analysis_data.csv'),fileEncoding='UTF-8-BOM')
av<-read.csv(file.path(ml,'mean_out_of_fold_predictions.csv'))
ps<-read.csv(file.path(ml,'performance_summary.csv'))
TNR<-'Times New Roman';col_train<-'#9e9e9e';col_test<-'#418ebd';col_darkred<-'#c21f1f'
colset<-c('Random Forest (OOF)'=col_train,'XGBoost (OOF)'=col_test)
cd_data<-av[av$model%in%c('Random forest','XGBoost'),]
cd_data$set<-factor(ifelse(cd_data$model=='XGBoost','XGBoost (OOF)','Random Forest (OOF)'),levels=names(colset))
cd_data$actual_prev<-100*cd_data$observed;cd_data$predicted_prev<-100*cd_data$predicted
cd_data$residual<-cd_data$predicted_prev-cd_data$actual_prev
cd_data$weight<-as.numeric(1/table(d$family)[cd_data$family])
xgb<-cd_data[cd_data$model=='XGBoost',];performance<-ps[ps$model=='XGBoost',]
fit_model<-lm(predicted_prev~actual_prev,data=xgb,weights=weight)
axis_min<--2;axis_max<-ceiling(max(c(cd_data$actual_prev,cd_data$predicted_prev))/5)*5
fitline<-data.frame(actual_prev=c(axis_min,axis_max));fitline$predicted_prev<-predict(fit_model,newdata=fitline)
box_theme<-theme_bw(base_family=TNR)+theme(panel.grid.major=element_line(linetype='dashed',color='gray85',linewidth=.5),
 panel.grid.minor=element_blank(),panel.border=element_rect(color='black',linewidth=.7,fill=NA),
 axis.title=element_text(size=12),axis.text=element_text(color='black',size=10),legend.key=element_blank(),
 legend.text=element_text(size=10,family=TNR),plot.margin=margin(0,0,0,0))
p_main<-ggplot(cd_data,aes(actual_prev,predicted_prev))+
 geom_abline(slope=1,intercept=0,linetype='dashed',color='black',linewidth=.65)+
 geom_line(data=fitline,color='black',linewidth=.9)+
 geom_point(aes(fill=set),shape=21,color='white',size=2.8,stroke=.4,alpha=.85)+
 scale_fill_manual(values=colset,name=NULL)+
 annotate('text',x=axis_min+2,y=axis_max-6,label=sprintf('XGBoost CV R² = %.3f\nRMSE = %.2f',performance$R2_mean,performance$RMSE_mean),
 color=col_darkred,size=4.0,hjust=0,family=TNR,lineheight=1.2)+
 coord_cartesian(xlim=c(axis_min,axis_max),ylim=c(axis_min,axis_max),expand=FALSE)+box_theme+
 theme(axis.title.x=element_blank(),axis.text.x=element_blank(),legend.position='inside',legend.position.inside=c(.75,.65),
 legend.background=element_rect(fill='white',color='black',linewidth=.4))+
 labs(y='Predicted Prevalence (%)')
p_top<-ggplot(cd_data,aes(actual_prev,fill=set,color=set))+
 geom_histogram(aes(y=after_stat(density)),position='identity',alpha=.4,bins=12,color='white',linewidth=.3)+
 geom_density(alpha=0,linewidth=.8)+scale_fill_manual(values=colset)+scale_color_manual(values=colset)+
 coord_cartesian(xlim=c(axis_min,axis_max),expand=FALSE)+theme_void(base_family=TNR)+theme(legend.position='none',plot.margin=margin(b=2))
p_right<-ggplot(cd_data,aes(predicted_prev,fill=set,color=set))+
 geom_histogram(aes(y=after_stat(density)),position='identity',alpha=.4,bins=12,color='white',linewidth=.3)+
 geom_density(alpha=0,linewidth=.8)+scale_fill_manual(values=colset)+scale_color_manual(values=colset)+
 coord_flip(xlim=c(axis_min,axis_max),expand=FALSE)+theme_void(base_family=TNR)+theme(legend.position='none',plot.margin=margin(l=2))
rfmae<-ps$MAE_mean[ps$model=='Random forest'];xgmae<-performance$MAE_mean
p_bottom<-ggplot(cd_data,aes(actual_prev,residual,fill=set))+geom_hline(yintercept=0,color='black',linewidth=.6)+
 geom_point(shape=21,color='white',size=2.5,stroke=.3,alpha=.85)+scale_fill_manual(values=colset)+
 annotate('text',x=axis_max-1,y=8,label=sprintf('MAE (RF CV) = %.2f\nMAE (XGB CV) = %.2f',rfmae,xgmae),hjust=1,size=3.7,family=TNR)+
 coord_cartesian(xlim=c(axis_min,axis_max),ylim=range(cd_data$residual)+c(-2,2),expand=FALSE)+box_theme+
 theme(legend.position='none',plot.margin=margin(t=5,r=0,b=5,l=0))+
 labs(x='Actual Observed Prevalence (%)',y='Residuals (pp)')
p_main<-p_main+labs(tag='a')+theme(plot.tag=element_text(face='bold',size=12,family=TNR),plot.tag.position=c(.035,.97))
p_bottom<-p_bottom+labs(tag='b')+theme(plot.tag=element_text(face='bold',size=12,family=TNR),plot.tag.position=c(.035,.93))
final_plot_1<-p_top+p_main+p_right+p_bottom+
 plot_layout(design='AAAA#\nBBBBC\nBBBBC\nBBBBC\nDDDD#',heights=c(.7,1,1,1,1.2),widths=c(1,1,1,1,.7))
ggsave(file.path(out,'FigS11.svg'),final_plot_1,width=8,height=9.5,device=svglite::svglite,bg='white')

# Original dual-axis SHAP code: same translucent bars, zero line, coloured dots.
sh<-read.csv(file.path(ml,'out_of_fold_SHAP.csv'));sh<-sh[sh$model=='XGBoost',]
shap_df<-aggregate(cbind(contribution,value)~row+feature,data=sh,mean)
shap_df$weight<-as.numeric(1/table(d$family)[d$family[match(shap_df$row,d$row)]])
pretty<-c(Age='Reported age summary',Male='Male proportion',LogN='Log sample size')
shap_summary<-do.call(rbind,lapply(unique(shap_df$feature),function(v){z<-shap_df[shap_df$feature==v,]
 data.frame(feature=v,mean_abs_shap=weighted.mean(abs(z$contribution),z$weight))}))
shap_summary<-shap_summary[order(shap_summary$mean_abs_shap),]
shap_summary$label_text<-sprintf('%s (%.4f)',pretty[shap_summary$feature],shap_summary$mean_abs_shap)
shap_df$label_text<-factor(shap_summary$label_text[match(shap_df$feature,shap_summary$feature)],levels=shap_summary$label_text)
shap_summary$label_text<-factor(shap_summary$label_text,levels=shap_summary$label_text)
shap_df$feature_value<-ave(shap_df$value,shap_df$feature,FUN=function(x)if(diff(range(x))>0)(x-min(x))/diff(range(x))else 0)
x_min<-min(-.08,min(shap_df$contribution)*1.15);x_max<-max(.08,max(shap_df$contribution)*1.12)
s_factor<-((x_max-x_min)*.70)/max(shap_summary$mean_abs_shap)
color_palette<-c('#143588','#6a92ce','#e8eaf1','#e97e6d','#b5182b')
p_shap<-ggplot()+geom_rect(data=shap_summary,
 aes(xmin=x_min,xmax=x_min+mean_abs_shap*s_factor,ymin=as.numeric(label_text)-.4,ymax=as.numeric(label_text)+.4,fill=as.numeric(label_text)),alpha=.25)+
 scale_fill_gradientn(colors=color_palette,guide='none')+geom_vline(xintercept=0,color='gray50',linewidth=.65)+
 geom_point(data=shap_df,aes(contribution,as.numeric(label_text),color=feature_value),
 position=position_jitter(height=.15,width=0,seed=999),size=1.8,alpha=.85,stroke=0)+
 scale_color_gradientn(colors=color_palette,name='Feature value',breaks=c(0,1),labels=c('Low','High'),
 guide=guide_colorbar(display='gradient',nbin=128,barheight=unit(3.4,'in'),barwidth=unit(.12,'in'),ticks=FALSE,title.position='right',
 title.theme=element_text(angle=-90,hjust=.5,family=TNR,size=10)))+
 scale_x_continuous(name='SHAP Value (Logit Units)',limits=c(x_min,x_max),expand=c(0,0),
 sec.axis=sec_axis(~(.-x_min)/s_factor,name='Mean(|SHAP Value|) (Logit Units)'))+
 scale_y_continuous(breaks=seq_len(nrow(shap_summary)),labels=levels(shap_summary$label_text),expand=c(.1,.1))+
 theme_bw(base_family=TNR)+theme(panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),
 panel.grid.major.y=element_line(color='gray92'),panel.border=element_rect(color='black',linewidth=.8,fill=NA),
 axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.text.y=element_text(size=11,color='black',hjust=1,margin=margin(r=5)),
 axis.text.x=element_text(size=10,color='black'),axis.title.x.bottom=element_text(margin=margin(t=10),size=12),
 axis.title.x.top=element_text(margin=margin(b=10),size=12),legend.position='right',
 legend.title=element_text(family=TNR),legend.text=element_text(family=TNR),plot.margin=margin(12,18,12,12))
ggsave(file.path(out,'FigS10.svg'),p_shap,width=9.5,height=6,device=svglite::svglite,bg='white')
write.csv(shap_summary,file.path(ml,'Untitled1_real_SHAP_summary.csv'),row.names=FALSE)
cat('Untitled1.R original templates updated with actual OOF data.\n')

# Optional deterministic SVG finishing: trimming retains vector geometry;
# font outlines retain Times New Roman in Word without font substitution.
# Pure-R editable SVGs are always saved in SVG_editable_text before finishing.
# Optional dependencies: Node sharp; Python fontTools and lxml; installed Times fonts.
# Set UGI_FINISH_SVG=0 to retain only native-R SVG export.
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


finish_svgs(out,output_root)
