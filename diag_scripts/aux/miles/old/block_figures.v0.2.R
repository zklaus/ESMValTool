######################################################
#--------Blocking routines plotting for ECE3---------#
#-------------P. Davini (Oct 2014)-------------------#
######################################################

miles.blockfigures<-function(exp,year1,year2,dataset_ref,year1_ref,year2_ref,season,FIGDIR,BLOCKDIR,cfg)
{
#R_LIBLOC<-Sys.getenv(c("R_LIBLOC"))
R_LIBLOC=.libPaths()[1]
source(cfg)

#correct folder to year and season dependence
REFDIR=paste(BLOCKDIR,"/",dataset_ref,"/",year1_ref,"_",year2_ref,"/",season,"/",sep="")
EXPDIR=paste(BLOCKDIR,"/",exp,"/",year1,"_",year2,"/",season,"/",sep="")
FIGDIR=paste(FIGDIR,"/",exp,"/",year1,"_",year2,"/",season,"/",sep="")
dir.create(FIGDIR,recursive=T)

#PROGDIR="/work/users/jost/blocking/MiLES"

#loadin packages
library("spam",lib.loc=R_LIBLOC)
library("maps",lib.loc=R_LIBLOC)
library("fields",lib.loc=R_LIBLOC)
library("ncdf4",lib.loc=R_LIBLOC)

#preparing routines
#source(paste(PROGDIR,"/script/basis_functions.R",sep=""))

#which fieds to plot
fieldlist=c("Block","Z500","BI","MGI")

#dataset_ref="ERAINTERIM"; year1_ref=1989; year2_ref=2010

#open reference fields
for (field in fieldlist)
		{
		nomefile=paste(REFDIR,"/BlockClim_",dataset_ref,"_",year1_ref,"_",year2_ref,"_",season,".nc",sep="")
		field_ref=ncdf.opener(nomefile,field,"Lon","Lat",rotate=F)
		assign(paste(field,"_ref",sep=""),field_ref)
		}

#loading Blocking file
#for (field in fieldlist)
#                {
#                nomefile=paste(BLOCKDIR,"/BlockClim_",exp,"_",year1,"_",year2,"_",season,".nc",sep="")
#                field_in=ncdf.opener(nomefile,field,"Lon","Lat",rotate=F)
#                assign(field,field_in)
#                }

outname=paste(EXPDIR,"/Block_",exp,"_",year1,"_",year2,"_",season,sep="")
load(outname)

#color palette to be used
palette0=tim.colors
palette1=colorRampPalette(c("white","orange","darkred"))
palette2=colorRampPalette(c("blue","white","red"))

#loop on fields
for (field in fieldlist)

	{
	#define field-dependent properties
	if (field=="Block")
	{
		field_exp=frequency; 
		color_field=palette1; color_diff=palette2
		lev_field=seq(0,35,3); lev_diff=seq(-10.5,10.5,1)
		legend_unit="Blocked Days (%)"; title_name="Instantaneous Blocking:"; legend_distance=3
	}

	if (field=="Z500")
        {
                field_exp=Z500mean;
                color_field=palette0; color_diff=palette2
                lev_field=seq(4800,6000,50); lev_diff=seq(-310,310,20)
                legend_unit="Geopotential Height (m)"; title_name="Z500:" ; legend_distance=4
        }

	if (field=="BI")
        {
                field_exp=BI;
                color_field=palette0; color_diff=palette2
                lev_field=seq(1,6,0.25); lev_diff=seq(-2.1,2.1,.2)
                legend_unit="BI index"; title_name="Blocking Intensity (BI):" ; legend_distance=3
        }

	if (field=="MGI")
        {
                field_exp=MGI;
                color_field=palette0; color_diff=palette2
                lev_field=seq(0,15,1); lev_diff=seq(-5.25,5.25,.5)
                legend_unit="MGI Index"; title_name="Meridional Gradient Inversion (MGI):" ; legend_distance=3
        }

	field_ref=get(paste(field,"_ref",sep=""))
	
	#secondary plot properties
	nlev_field=length(lev_field)-1
	nlev_diff=length(lev_diff)-1
	lat_lim=c(20,90)
	info_exp=paste(exp,year1,"-",year2,season)
	info_ref=paste(dataset_ref,year1_ref,"-",year2_ref,season)

	#final plot production
	figname=paste0(FIGDIR,"/",field,"_",exp,"_",year1,"_",year2,"_",season,".",output_file_type)
	print(figname)

        ## Chose output format for figure
        if (tolower(output_file_type) == "png") {
           png(filename = figname,
            width <- png_width,
            height <- png_height,
            units <- png_units,
            pointsize <- png_pointsize,
            bg <- png_bg)
        } else if (tolower(output_file_type) == "pdf") {
            pdf(file=figname,width=pdf_width,height=pdf_height,onefile=T)
        } else if (tolower(output_file_type) == "eps") {
            setEPS(width=pdf_width,height=pdf_height,onefile=T,paper="special")
            postscript(figname)
        }

	par(mfrow=c(3,1),cex.main=2,cex.axis=1.5,cex.lab=1.5,mar=c(5,5,4,8),oma=c(1,1,1,1))

	#main experiment plot
	filled.contour3(ics,ipsilon,field_exp,xlab="Longitude",ylab="Latitude",main=paste(title_name,info_exp),levels=lev_field,color.palette=color_field,ylim=lat_lim)
	map("world",regions=".",interior=F,exact=F,boundary=T,add=T)

	#reference field plot
	filled.contour3(ics,ipsilon,field_ref,xlab="Longitude",ylab="Latitude",main=paste(title_name,info_ref),levels=lev_field,color.palette=color_field,ylim=lat_lim)
	map("world",regions=".",interior=F,exact=F,boundary=T,add=T)
	image.plot(legend.only=TRUE,axis.args=list(cex.axis=1.5),legend.width=1,zlim=c(min(lev_field),max(lev_field)),col=color_field(nlev_field),legend.args=list(side=4,line=legend_distance,cex=1.2,text=legend_unit))

	#delta field plot
	filled.contour3(ics,ipsilon,field_exp-field_ref,xlab="Longitude",ylab="Latitude",main=paste(title_name,"Difference"),levels=lev_diff,color.palette=color_diff,ylim=lat_lim)
	map("world",regions=".",interior=F,exact=F,boundary=T,add=T)
	image.plot(legend.only=TRUE,axis.args=list(cex.axis=1.5),legend.width=1,zlim=c(min(lev_diff),max(lev_diff)),col=color_diff(nlev_diff),legend.args=list(side=4,line=legend_distance,cex=1.2,text=legend_unit))

	dev.off()
	}

}
