getauc<-function(x,y,ybase,splnmeth="fmm",nval=NA,lbl=NA,ylbl=NA,xlbl=NA,plot=T,fnprefix="",
                 replacefile=F,wd=1200,ht=6000,res=600,interpolmethod="spline",draw.polygon=F,
                 psymbol=1,psize=1,daysfrombasedate=0,ID_labels,ID_vals,dailyhlimit=3,
                 output_folder=NA,sensor_changes=NULL)
{
  library(rootSolve)
  library(MESS)
  #library(plyr)
  
  subarea<-function(x,y2,r1,r2)
  {
    r11<-r1
    r12<-r21<-r11+(r2-r1)/2
    r22<-r2
    cat(paste("      Splitting: ",r11," to ",r12," to ",r22,"\n",sep=""))  
    area1<-try(auc(x,y2,type=interpolmethod,from=r11,to=r12),silent=TRUE)
    area2<-try(auc(x,y2,type=interpolmethod,from=r21,to=r22),silent=TRUE)
    if (is.numeric(area1) & is.numeric(area2))
    {
      cat(paste("      1st Area + 2nd Area = ",area1+area2,"\n",sep=""))
      return(area1+area2)
    } else
    {
      # Break once more
      if (is.numeric(area1)==FALSE)
      {
        area1<-subarea(x,y2,r11,r12)
      } else
      {
        cat(paste("      1st Area = ",area1,"\n",sep=""))
      }
      if (is.numeric(area2)==FALSE)
      {
        area2<-subarea(x,y2,r21,r22)
      } else
      {
        cat(paste("      2nd Area = ",area2,"\n",sep=""))
      }
      if (is.na(area1) | is.na(area2)) return(NA) else return(area1+area2)
    }
  }
  
  if (anyNA(y)==F)
  {
    if (length(x)<5)
    {
      cat(paste("Calculating ",lbl," (< 5 points)\n",sep=""))
      return(list("ybase"=ybase,"ymin"=NA,"ymax"=NA,"xmin"=NA,"xmax"=NA,
                  "roots"=NA,"auc"=NA,"tuc"=NA,"aac"=NA,"tac"=NA,"freq"=NA))
    } 
    cat(paste("Calculating ",lbl,sep=""))
    y.min<-min(y)
    y.max<-max(y)
    x.min<-min(x)
    x.max<-max(x)
    y2<-y-ybase
    if (is.na(output_folder)) output_folder<-getwd()
    if (replacefile==T & plot==T) fn<-paste(output_folder,"/",fnprefix,"plot 1.tiff",sep="")
    if (replacefile==F & plot==T) fn<-funcSerialFilename(paste(output_folder,"/",fnprefix,"plot.tiff",sep=""))
    if (is.na(nval)) nval<-10*length(x)
    splfunc<-splinefun(x,y2,method=splnmeth,ties=mean)
    roots<-uniroot.all(splfunc,c(min(x),max(x)),n=nval)
    roots<-roots[order(roots)]
    # Add first and last X points as roots
    roots<-c(x[1],roots)
    roots<-c(roots,x[length(x)])
    if(length(roots)<2)
    {
      cat(" (< 2 roots)\n")
      if (plot==T)
      {
        tiff(fn,width=wd,height=ht,res=res,bg="white",compress="lzw")
        par(mfrow=c(2,1),bty="l")
        if (is.na(lbl)==T) lbl<-""
        if (is.na(ylbl)==T) ylbl1<-"Y values" else ylbl1<-ylbl
        if (is.na(xlbl)==T) xlbl1<-"X values" else xlbl1<-xlbl
        # Fit the linear model
        lm_res <- lm(y ~ x)
        # Extract the coefficients for the equation
        intercept <- coef(lm_res)[1]
        slope <- coef(lm_res)[2]
        # Round the coefficients to 3 significant digits
        intercept_rounded <- signif(intercept, digits = 3)
        slope_rounded <- signif(slope, digits = 3)
        # Format the rounded coefficients to fixed-point format
        intercept_formatted <- formatC(intercept_rounded, format = "f", digits = max(0, 3 - floor(log10(abs(intercept_rounded)) + 1)))
        slope_formatted <- formatC(slope_rounded, format = "f", digits = max(0, 3 - floor(log10(abs(slope_rounded)) + 1)))
        # Format the equation with the intercept first and then the slope
        if (slope_rounded < 0) {
          equation <- paste("y = ", intercept_formatted, " - ", formatC(abs(slope_rounded), format = "f", digits = max(0, 3 - floor(log10(abs(slope_rounded)) + 1))), "x")
        } else {
          equation <- paste("y = ", intercept_formatted, " + ", slope_formatted, "x")
        }
        plot(x,y,main=lbl,xlab=paste(xlbl1, equation, sep=" | "),ylab=ylbl1,type="n",
             sub=paste("Y: [min, max] = [",formatC(min(y),digits=3,width=7,format="f"),", ",formatC(max(y),digits=3,width=7,format="f"), "], ",
                       "X: [min, max] = [",formatC(min(x),digits=3,width=7,format="f"),", ",formatC(max(x),digits=3,width=7,format="f"), "]", sep=""))
        # Create color vector based on sensor changes
        point_colors <- rep("dodgerblue4", length(x))  # default color
        if (!is.null(sensor_changes) && nrow(sensor_changes) > 0) {
          # For each sensor change point
          for (i in seq_len(nrow(sensor_changes))) {
            # Alternate between green and blue
            current_color <- if (i %% 2 == 1) "darkgreen" else "dodgerblue4"
            point_colors[x > sensor_changes$t[i]] <- current_color
          }
        }
        # Replot the points
        points(x,y,col=point_colors,pch=psymbol,cex=psize)
        # Add the linear trendline
        abline(lm_res, col="red", lwd=2)
        # draw spline lines
        s<-spline(x,y,n=nval,ties=mean)
        lines(s,col="dodgerblue4")
        # Add vertical lines at sensor changes if they exist
        if (!is.null(sensor_changes) && nrow(sensor_changes) > 0) {
          abline(v=sensor_changes$t, col="darkgreen", lty=2, lwd=1)
        }
        # Plot 2
        # Difference plot
        if (is.na(ylbl)==T) ylbl2<-paste("Y values - ",ybase,sep="") else ylbl2<-paste(ylbl," - ",ybase,sep="")
        if (is.na(xlbl)==T) xlbl2<-"X values" else xlbl2<-xlbl
        plot(x,y2,xlab=xlbl2,ylab=ylbl2)
        abline(h=0)
        abline(v=roots,lty=3)
        lines(spline(x,y2,n=nval))
        # replot the points
        points(x,y2,col="firebrick",pch=psymbol,cex=psize)
        # Add vertical lines at sensor changes if they exist
        if (!is.null(sensor_changes) && nrow(sensor_changes) > 0) {
          abline(v=sensor_changes$t, col="darkgreen", lty=2, lwd=1)
        }
        dev.off()
      }
      return(list("ybase"=ybase,"ymin"=y.min,"ymax"=y.max,
                  "xmin"=x.min,"xmax"=x.max,"roots"=NA,
                  "auc"=NA,"tuc"=NA,"aac"=NA,"tac"=NA,"freq"=NA))        
    }
    # tb = time between
    below.areasum<-0
    below.tbsum<-0
    below.freqsum<-0
    below.df<-data.frame(t=numeric(),area=numeric(),tb=numeric())
    below.df[nrow(below.df)+1,]<-c(0,0,0)
    above.areasum<-0
    above.tbsum<-0
    above.freqsum<-0
    above.df<-data.frame(t=numeric(),area=numeric(),tb=numeric())
    above.df[nrow(above.df)+1,]<-c(0,0,0)
    df.inter<-data.frame(t1=as.double(),t2=as.double(),above.area=as.double(),above.time=as.double(),below.area=as.double(),below.time=as.double())
    if (length(ID_labels)>0) df.inter<-cbind(as.data.frame(replicate(length(ID_labels), numeric()), col.names = ID_labels), df.inter)
    areaOK<-TRUE
    # progress bar
    # https://ryouready.wordpress.com/2009/03/16/r-monitor-function-progress-with-a-progress-bar/
    pb.total<-length(roots)-1
    pb.title<-ifelse(is.na(lbl)==T,"Processing",paste("Processing ",lbl,sep=""))
    pb<-tkProgressBar(title=pb.title,min=0,max=pb.total,width=400)
    for (i in 1:(length(roots)-1))
    {
      setTkProgressBar(pb,i,label=paste("Calculating area (",round(i/pb.total*100, 0),"% done)",sep=""))
      tcl("update")
      r1<-roots[i]
      r2<-roots[i+1]
      if (r1!=0 | r2!=0)
      {
        above.area<-0
        above.time<-0
        below.area<-0
        below.time<-0
        if (r1<r2)
        {
          areacalc<-try(auc(x,y2,type=interpolmethod,from=r1,to=r2),silent=TRUE)
          if (!is.numeric(areacalc))
          {
            # Break the interval in half the try to solve it
            cat("\n")
            cat(paste("   * sub-splitting roots (",i," and ",i+1,") ", r1," and ", r2, "\n",sep=""))
            area1<-subarea(x,y2,r1,r2)
            if (is.na(area1))
            {
              areacalc<-0
              areaOK<-FALSE
            } else
            {
              cat(paste("      Total Area = ",area1,"\n",sep=""))  
              areacalc<-area1
            }
          }
          if (areacalc>0)
          {
            above.area<-above.area+areacalc
            above.time<-above.time+(r2-r1)
            above.areasum<-above.areasum+areacalc
            above.freqsum<-above.freqsum+1
            above.tbsum<-above.tbsum+(r2-r1)
            above.df[nrow(above.df)+1,]<-c(r2,above.areasum,above.tbsum)
          } else
          {
            below.area<-below.area+areacalc
            below.time<-below.time+(r2-r1)
            below.areasum<-below.areasum-areacalc   #areacalc is negative for below
            below.freqsum<-below.freqsum+1
            below.tbsum<-below.tbsum+(r2-r1)
            below.df[nrow(below.df)+1,]<-c(r2,below.areasum,below.tbsum)
          }
        }
        # Store current event
        if (length(ID_labels)>0)
        {
          df.inter[nrow(df.inter)+1,]<-c(c(ID_vals),c(r1,r2,above.area,above.time,-below.area,below.time))
        } else
        {
          df.inter[nrow(df.inter)+1,]<-c(r1,r2,above.area,above.time,-below.area,below.time)
        }
      }
    }
    df.inter<-lapply(df.inter, function(x) type.convert(as.character(x)))
    df.inter<-as.data.frame(df.inter)
    if (areaOK==FALSE)
    {
      below.areasum<--1
      below.tbsum<--1
      below.freqsum<--1
      above.areasum<--1
      above.tbsum<--1
      above.freqsum<--1
      cat(" (FAILED)\n")
    } else cat("   Finished calculation\n")
    if (plot==T)
    {
      setTkProgressBar(pb,0,label=paste("Creating graphics (",round(0, 0),"% done)",sep=""))
      tcl("update")
      tiff(fn,width=wd,height=ht,res=res,bg="white",compress="lzw")
      layout(matrix(c(1,1,2,2,3,4), 3, 2, byrow = TRUE))
      if (is.na(lbl)==T) lbl<-""
      if (is.na(ylbl)==T) ylbl1<-"Y values" else ylbl1<-ylbl
      if (is.na(xlbl)==T) xlbl1<-"X values" else xlbl1<-xlbl
      # plot 1
      # Fit the linear model
      lm_res <- lm(y ~ x)
      # Extract the coefficients for the equation
      intercept <- coef(lm_res)[1]
      slope <- coef(lm_res)[2]
      # Round the coefficients to 3 significant digits
      intercept_rounded <- signif(intercept, digits = 3)
      slope_rounded <- signif(slope, digits = 3)
      # Format the rounded coefficients to fixed-point format
      intercept_formatted <- formatC(intercept_rounded, format = "f", digits = max(0, 3 - floor(log10(abs(intercept_rounded)) + 1)))
      slope_formatted <- formatC(slope_rounded, format = "f", digits = max(0, 3 - floor(log10(abs(slope_rounded)) + 1)))
      # Format the equation with the intercept first and then the slope
      if (slope_rounded < 0) {
        equation <- paste("y = ", intercept_formatted, " - ", formatC(abs(slope_rounded), format = "f", digits = max(0, 3 - floor(log10(abs(slope_rounded)) + 1))), "x")
      } else {
        equation <- paste("y = ", intercept_formatted, " + ", slope_formatted, "x")
      }
      plot(x,y,main=lbl,xlab=paste(xlbl1, equation, sep=" | "),ylab=ylbl1,type="n",
           sub=paste("Y: [min, max] = [",formatC(min(y),digits=3,width=7,format="f"),", ",formatC(max(y),digits=3,width=7,format="f"), "], ",
                     "X: [min, max] = [",formatC(min(x),digits=3,width=7,format="f"),", ",formatC(max(x),digits=3,width=7,format="f"), "]", sep=""))
      # Create color vector based on sensor changes
      point_colors <- rep("dodgerblue4", length(x))  # default color
      if (!is.null(sensor_changes) && nrow(sensor_changes) > 0) {
        # For each sensor change point
        for (i in seq_len(nrow(sensor_changes))) {
          # Alternate between green and blue
          current_color <- if (i %% 2 == 1) "darkgreen" else "dodgerblue4"
          point_colors[x > sensor_changes$t[i]] <- current_color
        }
      }
      # Replot the points
      points(x,y,col=point_colors,pch=psymbol,cex=psize)
      # Add the linear trendline
      abline(lm_res, col="red", lwd=2)
      # draw spline lines
      s<-spline(x,y,n=nval,ties=mean)
      lines(s,col="dodgerblue4")
      # Add vertical lines at sensor changes if they exist
      if (!is.null(sensor_changes) && nrow(sensor_changes) > 0) {
        abline(v=sensor_changes$t, col="darkgreen", lty=2, lwd=1)
      }
      setTkProgressBar(pb,25,label=paste("Creating graphics (",round(25, 0),"% done)",sep=""))
      tcl("update")
      
      # plot 2
      ylbl1<-paste(ylbl1," (base = ", ybase, ")", sep="")
      plot(x,y,xlab=xlbl1,ylab=ylbl1,type="n",
           sub=paste("Y: [min, max] = [",formatC(min(y),digits=3,width=7,format="f"),", ",formatC(max(y),digits=3,width=7,format="f"), "], ",
                     "X: [min, max] = [",formatC(min(x),digits=3,width=7,format="f"),", ",formatC(max(x),digits=3,width=7,format="f"), "]", sep=""))
      # Create color vector based on sensor changes
      point_colors <- rep("dodgerblue4", length(x))  # default color
      if (!is.null(sensor_changes) && nrow(sensor_changes) > 0) {
        # For each sensor change point
        for (i in seq_len(nrow(sensor_changes))) {
          # Alternate between green and blue
          current_color <- if (i %% 2 == 1) "darkgreen" else "dodgerblue4"
          point_colors[x > sensor_changes$t[i]] <- current_color
        }
      }
      # Replot the points
      points(x,y,col=point_colors,pch=psymbol,cex=psize)
      # reference line
      abline(h=ybase)
      # draw spline lines
      s<-spline(x,y,n=nval,ties=mean)
      lines(s,col="dodgerblue4")
      # Add vertical lines at sensor changes if they exist
      if (!is.null(sensor_changes) && nrow(sensor_changes) > 0) {
        abline(v=sensor_changes$t, col="darkgreen", lty=2, lwd=1)
      }
      setTkProgressBar(pb,25,label=paste("Creating graphics (",round(50, 0),"% done)",sep=""))
      tcl("update")
      # draw polygon
      if (draw.polygon==T)
      {
        for (i in 1:(length(roots)-1))
        {
          setTkProgressBar(pb,i,label=paste("Drawing polygons (",round(i/pb.total*100, 0),"% done)",sep=""))
          tcl("update")
          nval2<-max(3,3*length(which(x>roots[i] & x<roots[i+1])))
          s2<-spline(x,y,n=nval2,xmin=roots[i],xmax=roots[i+1],ties=mean)
          if (sum(s2$y-ybase)<0) 
          {
            polygon(c(roots[i],s2$x,roots[i+1]),c(ybase,s2$y,ybase),col=rgb(1,0,0,0.1),border=NA)
          } else
          {
            polygon(c(roots[i],s2$x,roots[i+1]),c(ybase,s2$y,ybase),col=rgb(0,0,1,0.1),border=NA)
          }
        }
      }
      
      # plot 3
      # Cumulative area and tb for above and below the base
      if (is.na(ylbl)==T) ylbl2<-"Cummulative area" else ylbl2<-paste("Cumulative area of ",ylbl,sep="")
      if (is.na(xlbl)==T) xlbl2<-"X values" else xlbl2<-xlbl
      sub1<-paste("Roots = ",length(roots), ", ",
                  "AAC = ",formatC(above.areasum,digits=3,width=7,format="f"), ", ",
                  "AUC = ",formatC(below.areasum,digits=3,width=7,format="f"), sep="")
      plot(c(above.df$t,below.df$t),c(above.df$area,-below.df$area),type="n",xlab=xlbl2,ylab=ylbl2,sub=sub1)
      abline(h=0)
      lines(above.df$t,above.df$area,col=rgb(0, 0, 1, 0.7),lwd=2)
      lines(below.df$t,-below.df$area,col=rgb(1, 0, 0, 0.7),lwd=2)
      polygon(c(0,above.df$t,max(above.df$t)),c(0,above.df$area,0),col=rgb(0, 0, 1, 0.2),border=NA)
      polygon(c(0,below.df$t,max(below.df$t)),c(0,-below.df$area,0),col=rgb(1, 0, 0, 0.2),border=NA)
      setTkProgressBar(pb,75,label=paste("Creating graphics (",round(75, 0),"% done)",sep=""))
      tcl("update")
      
      # plot 4
      if (is.na(ylbl)==T) ylbl2<-"Cummulative time" else ylbl2<-paste("Cumulative time of ",ylbl,sep="")
      if (is.na(xlbl)==T) xlbl2<-"X values" else xlbl2<-xlbl
      sub1<-paste("TAC = ",formatC(above.tbsum,digits=3,width=7,format="f"), ", ",
                  "TUC = ",formatC(below.tbsum,digits=3,width=7,format="f"), ", ",
                  "Frequency = ",formatC(above.freqsum,digits=0,width=2,format="f"), sep="")      
      plot(c(above.df$t,below.df$t),c(above.df$tb,-below.df$tb),type="n",xlab=xlbl2,ylab=ylbl2,sub=sub1)
      abline(h=0)
      lines(above.df$t,above.df$tb,col=rgb(0, 0, 1, 0.7),lwd=2)
      lines(below.df$t,-below.df$tb,col=rgb(1, 0, 0, 0.7),lwd=2)
      polygon(c(0,above.df$t,max(above.df$t)),c(0,above.df$tb,0),col=rgb(0, 0, 1, 0.2),border=NA)
      polygon(c(0,below.df$t,max(below.df$t)),c(0,-below.df$tb,0),col=rgb(1, 0, 0, 0.2),border=NA)
      setTkProgressBar(pb,100,label=paste("Creating graphics (",round(100, 0),"% done)",sep=""))
      tcl("update")
      dev.off()
    }
    
    # Save interval duration file
    cat(paste("   Saving interval data on '",fn,"'\n",sep=""))
    fn<-paste(output_folder,"/",fnprefix,"inter.csv",sep="")
    df.inter$d1<-floor(as.numeric(daysfrombasedate)+df.inter$t1*60/1440)
    df.inter$d2<-floor(as.numeric(daysfrombasedate)+df.inter$t2*60/1440)
    write.csv(df.inter,fn)
    
    if (exists("df.inter.all", envir = .GlobalEnv) && !is.null(df.inter.all))  {
      df.inter.all <<- rbind(df.inter.all, df.inter)
    } else {
      df.inter.all <<- df.inter
    }
    
    # Calculate daily below and above
    cat(paste("   Saving interval DAILY data on '",fn,"'\n",sep=""))
    cat("\n\n")
    fn<-paste(output_folder,"/",fnprefix,"inter_daily.csv",sep="")
    #df.inter.daily<-ddply(df.inter,"d1",summarize,above.area=sum(above.area),above.time=sum(above.time),
    #                      below.area=sum(below.area),below.time=sum(below.time))
    df.inter.daily <- df.inter %>%
      group_by(d1) %>%
      summarize(above.area = sum(above.area), above.time = sum(above.time),
                below.area = sum(below.area), below.time = sum(below.time))
    write.csv(df.inter.daily,fn)
    
    if (plot==T)
    {
      setTkProgressBar(pb,0,label=paste("Creating graphics (",round(0, 0),"% done)",sep=""))
      tcl("update")
      df.inter.daily.limit<-subset(df.inter.daily,df.inter.daily$below.time>dailyhlimit)
      if (nrow(df.inter.daily.limit)>0)
      {
        dailyhlimit.events<-sum(df.inter.daily.limit$below.time>dailyhlimit)
        dailyhlimit.time<-sum(df.inter.daily.limit$below.time)
        if ((x.max-x.min)>0) dailyhlimit.prop<-dailyhlimit.time/(x.max-x.min) else dailyhlimit.prop<-0
        setTkProgressBar(pb,33,label=paste("Creating graphics (",round(25, 0),"% done)",sep=""))
        tcl("update")
        fn<-paste(output_folder,"/",fnprefix,"plot daily limit pH.tiff",sep="")
        tiff(fn,width=wd,height=ht,res=res,bg="white",compress="lzw")
        ylbl1<-paste("Time ruminal pH was below ",ybase," for more than ",dailyhlimit," consecutive hours",sep="")
        xlbl1<-"Time, days"
        plot(df.inter.daily.limit$d1,df.inter.daily.limit$below.time,main=lbl,xlab=xlbl1,ylab=ylbl1,type="h",
             sub=paste("Consecutive ruminal pH < ",ybase," for more than ",dailyhlimit," h: Events : ",formatC(dailyhlimit.events,digits=1,width=5,format="f"),"; ",
                       "Time : ",formatC(dailyhlimit.time,digits=1,width=5,format="f")," h, which represents ",
                       formatC(dailyhlimit.prop*100,digits=2,width=5,format="f"), "% of the total time (",formatC((x.max-x.min),digits=0,width=4,format="f")," h)", sep=""))
        setTkProgressBar(pb,33,label=paste("Creating graphics (",round(50, 0),"% done)",sep=""))
        tcl("update")
        # replot the points
        setTkProgressBar(pb,33,label=paste("Creating graphics (",round(75, 0),"% done)",sep=""))
        tcl("update")
        points(df.inter.daily.limit$d1,df.inter.daily.limit$below.time,col="dodgerblue4",pch=psymbol,cex=psize+2)
        # reference lines
        abline(h=24,lwd=2)
        abline(h=dailyhlimit,lwd=2)
        setTkProgressBar(pb,100,label=paste("Creating graphics (",round(100, 0),"% done)",sep=""))
        tcl("update")
        dev.off()
      }
    }
    close(pb)
    if (exists("df.inter.all")) df.inter.all<<-df.inter.all[complete.cases(df.inter.all),]
    output<-list("ybase"=ybase,"ymin"=y.min,"ymax"=y.max,
                 "xmin"=x.min,"xmax"=x.max,"roots"=length(roots),
                 "auc"=below.areasum,"tuc"=below.tbsum,
                 "aac"=above.areasum,"tac"=above.tbsum,
                 "freq"=above.freqsum)
    return(output)
  } else
  {
    cat(paste("Calculating ",lbl," (No valid points?)",sep=""))
    return(list("ybase"=ybase,"ymin"=NA,"ymax"=NA,
                "xmin"=NA,"xmax"=NA,"roots"=NA,
                "auc"=NA,"tuc"=NA,"aac"=NA,"tac"=NA,"freq"=NA))
  }
}
