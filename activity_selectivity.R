library(RCurl)

#reads names and ages of time intervals
time_url <- RCurl::getURL("https://paleobiodb.org/data1.2/intervals/list.txt?scale=1&scale_level=5&limit=all", ssl.verifypeer=F)
time_int <- read.csv(text=time_url)
stages <- subset(time_int, time_int$scale_level == 5)
bin_midpt <- colMeans(rbind(stages$max_ma, stages$min_ma))


taxa <- c("Brachiopoda,Bivalvia,Gastropoda,Porifera,Cnidaria,Echinodermata,Bryozoa,Trilobita,Malacostraca,Ostracoda")

#reads occurrences based on specified taxa and age range
occ_url <- RCurl::getURL(paste("https://paleobiodb.org/data1.2/occs/list.txt?base_name=",taxa,"&interval=Artinskian,Bathonian&envtype=marine&show=phylo,geo,ident&limit=all",sep=""), ssl.verifypeer=F, timeout=300)
occs <- read.csv(text=occ_url)

#FILE PREPARATION AND DATA CLEANING

#read stage match
time_url <- RCurl::getURL("https://docs.google.com/spreadsheets/d/1DPWsRF3DmGEo0eJXjvotuDJLphFaSYqjfCpYXWSgDvY/export?gid=0&format=csv",ssl.verifypeer = FALSE)
time_match <- read.csv(text=time_url)

#finds occurrences restricted to a single stage
occs$early_stage <- time_match$matched_stage[match(occs$early_interval, time_match$int_name)]
occs$late_stage <- time_match$matched_stage[match(occs$late_interval, time_match$int_name)]

single_int <- subset(occs, occs$late_interval=="" & is.na(occs$early_stage)==F)

multi_int <- subset(occs, occs$late_interval!="")
multi_int <- subset(multi_int, multi_int$early_stage==multi_int$late_stage)

resolved_occs <- rbind(single_int, multi_int)

#deletes occurrences not resolved to at least genus level using classified and unclassified species
cleaned_occs <- subset(resolved_occs, resolved_occs$accepted_rank %in% c("genus", "subgenus", "species"))

#deletes occurrences where genus is qualified with question mark, quotations, cf. or aff.
cleaned_occs <- subset(cleaned_occs, cleaned_occs$primary_reso=="" | cleaned_occs$primary_reso=="n. gen.")

cleaned_occs$primary_name <- factor(cleaned_occs$primary_name)
cleaned_occs$early_stage <- factor(cleaned_occs$early_stage)

stages_final <- subset(stages, stages$interval_name %in% levels(cleaned_occs$early_stage))

#reads activity data
activity_level_url <- RCurl::getURL("https://raw.githubusercontent.com/mclapham/grc2014/master/activity_level.csv",ssl.verifypeer=F)
activity_level <- read.csv(text=activity_level_url)

#finds row in activity data frame corresponding to each occurrence
activity_row <- cbind(match(cleaned_occs$class,activity_level$taxon),match(cleaned_occs$order,activity_level$taxon),match(cleaned_occs$family,activity_level$taxon))
activity_row <- apply(activity_row, 1, function(x) min(x, na.rm=T))

#adds activity quotient to occurrences
cleaned_occs$activity <- activity_level$activity_quotient[activity_row]


#ANALYSES
#mean activity per interval
#all taxa
all_occs <- cleaned_occs

stage_activity_all <- sapply(split(all_occs$activity, all_occs$early_stage), function(x) mean(x, na.rm=T))

activity_res_all <- data.frame(activity=stage_activity_all[match(stages_final$interval_name, names(stage_activity_all))], age=rowMeans(stages_final[,c("max_ma","min_ma")]))

#excluding Ostracoda
cleaned_occs <- subset(cleaned_occs, cleaned_occs$class!="Ostracoda")

stage_activity <- sapply(split(cleaned_occs$activity, cleaned_occs$early_stage), function(x) mean(x, na.rm=T))

activity_res <- data.frame(activity=stage_activity[match(stages_final$interval_name, names(stage_activity))], age=rowMeans(stages_final[,c("max_ma","min_ma")]))

perm <- subset(activity_res, activity_res$age>252)
trias <- subset(activity_res, activity_res$age<252 & activity_res$age>201)
jur <- subset(activity_res, activity_res$age<201)

perm_lm <- lm(activity~age, data=perm)
trias_lm <- lm(activity~age, data=trias)
jur_lm <- lm(activity~age, data=jur)

perm_prd <- predict(perm_lm, data.frame(age=seq(252.17,290,length.out=500)), interval="confidence")
trias_prd <- predict(trias_lm, data.frame(age=seq(201.3,252.17,length.out=500)),interval="confidence")
jur_prd <- predict(jur_lm, data.frame(age=seq(150,201.3,length.out=500)),interval="confidence")


#by genus
genus_activity <- sapply(split(all_occs, all_occs$early_stage), function(x) sapply(split(x$activity, x$primary_name), function(y) mean(y, na.rm=T)))

stage_activity_genus <- apply(genus_activity, 2, function(x) mean(x, na.rm=T))

activity_res_genus <- data.frame(activity=stage_activity_genus[match(stages_final$interval_name, names(stage_activity_genus))], age=rowMeans(stages_final[,c("max_ma","min_ma")]))

trias_genus <- subset(activity_res_genus, activity_res_genus$age<252 & activity_res_genus$age>201)
jur_genus <- subset(activity_res_genus, activity_res_genus$age<201)

trias_genus_lm <- lm(activity~age, data=trias_genus)
jur_genus_lm <- lm(activity~age, data=jur_genus)

trias_genus_prd <- predict(trias_genus_lm, data.frame(age=seq(201.3,252.17,length.out=500)),interval="confidence")
jur_genus_prd <- predict(jur_genus_lm, data.frame(age=seq(150,201.3,length.out=500)),interval="confidence")


#PLOT ACTIVITY TRENDS
pdf("activity_trends.pdf", width=7)
par(mgp=c(2,0.75,0))
par(mfrow=c(2,1))
par(mar=c(3, 3, 1, 1))
plot(activity_res$age, activity_res$activity, type="n", xlab="Age (Ma)", ylab="Mean activity level (by occurrence)", xlim=rev(range(activity_res$age)), ylim=c(1.5, 2.45))

polygon(c(seq(252.17,290,length.out=500),rev(seq(252.17,290,length.out=500))), c(perm_prd[,"lwr"],rev(perm_prd[,"upr"])), border=NA, col="gray")
polygon(c(seq(201.3,252.17,length.out=500),rev(seq(201.3,252.17,length.out=500))), c(trias_prd[,"lwr"],rev(trias_prd[,"upr"])), border=NA, col="gray")
polygon(c(seq(150,201.3,length.out=500),rev(seq(150,201.3,length.out=500))), c(jur_prd[,"lwr"],rev(jur_prd[,"upr"])), border=NA, col="gray")

segments(252.17, 1.59, 252.17, 3)
text(252.17, 1.57, "End-Permian", cex=0.6)
segments(201.3, 1.55, 201.3, 3)
text(201.3, 1.53, "End-Triassic", cex=0.6)
segments(259.9, 1.55, 259.9, 3, lty=2)
text(259.9, 1.53, "Guadalupian", cex=0.6)
segments(174.1, 1.55, 174.1, 3, lty=2)
text(174.1, 1.53, "Toarcian", cex=0.6)

segments(252.17, perm_prd[1,"fit"], 290, perm_prd[nrow(perm_prd),"fit"])
segments(201.3, trias_prd[1,"fit"], 252.17, trias_prd[nrow(trias_prd),"fit"])
segments(150, jur_prd[1,"fit"], 201.3, jur_prd[nrow(jur_prd),"fit"])

points(activity_res_all$age, activity_res_all$activity)

points(activity_res$age, activity_res$activity, pch=16)

legend(288, 2.47, c("All taxa", "Excluding Ostracoda"), cex=0.7, pch=c(1,16), bg="white", box.col="white")

#add timescale
rect(stages$max_ma, 1.47, stages$min_ma, 1.51, col=paste(stages$color))
text(bin_midpt, 1.49, strtrim(stages$interval_name,1),cex=0.6)

box()

mtext(font=2, "(a)", adj=0)

plot(activity_res_genus$age, activity_res_genus$activity, type="n", xlab="Age (Ma)", ylab="Mean activity level (by genus)", xlim=rev(range(activity_res$age)), ylim=c(1.7, 2.6))

polygon(c(seq(201.3,252.17,length.out=500),rev(seq(201.3,252.17,length.out=500))), c(trias_genus_prd[,"lwr"],rev(trias_genus_prd[,"upr"])), border=NA, col="gray")
polygon(c(seq(150,201.3,length.out=500),rev(seq(150,201.3,length.out=500))), c(jur_genus_prd[,"lwr"],rev(jur_genus_prd[,"upr"])), border=NA, col="gray")

segments(252.17, 1.79, 252.17, 3)
text(252.17, 1.77, "End-Permian", cex=0.6)
segments(201.3, 1.75, 201.3, 3)
text(201.3, 1.73, "End-Triassic", cex=0.6)
segments(259.9, 1.75, 259.9, 3, lty=2)
text(259.9, 1.73, "Guadalupian", cex=0.6)
segments(174.1, 1.75, 174.1, 3, lty=2)
text(174.1, 1.73, "Toarcian", cex=0.6)


segments(201.3, trias_genus_prd[1,"fit"], 252.17, trias_genus_prd[nrow(trias_genus_prd),"fit"])
segments(150, jur_genus_prd[1,"fit"], 201.3, jur_genus_prd[nrow(jur_genus_prd),"fit"])

points(activity_res_genus$age, activity_res_genus$activity, pch=16)

#add timescale
rect(stages$max_ma, 1.67, stages$min_ma, 1.71, col=paste(stages$color))
text(bin_midpt, 1.69, strtrim(stages$interval_name,1),cex=0.6)

box()

mtext(font=2, "(b)", adj=0)
dev.off()



#activity score for each genus
genus_activity <- sapply(split(cleaned_occs$activity,cleaned_occs$primary_name), mean)

genus_matrix_raw <- sapply(split(cleaned_occs, cleaned_occs$early_stage), function(x) table(x$primary_name))

#orders genus matrix so that time intervals are in chronologic order
genus_matrix <- genus_matrix_raw[,match(stages_final$interval_name, colnames(genus_matrix_raw))]

genus_matrix <- genus_matrix[,rev(seq(ncol(genus_matrix)))]

bin_top <- time_int$min_ma[match(colnames(genus_matrix),time_int$interval_name)]

#LOGISTIC REGRESSION FOR SURVIVAL AS FUNCTION OF ACTIVITY
log_odds_ratio<-numeric(0)
log_odds_error<-numeric(0)
log_odds_p<-numeric(0)

#logistic regression based on boundary-crosser extinction
for (i in 1:(ncol(genus_matrix)-2)) {
  j = i+1
  prior_int <- genus_matrix[,1:j]
  cohort <- subset(genus_matrix, apply(data.frame(prior_int[,seq(j-1)]), 1, sum) > 0 & prior_int[,j] > 0)
  extinct_table <- apply(cohort, 1, function(x) ifelse(x[seq((j+1),length(x))]>0,1,0))
  if (class(extinct_table)=="numeric") {
    extinct <- ifelse(extinct_table>0, 1, 0)
  } else {
    extinct <- ifelse(apply(extinct_table, 2, sum)>0, 1, 0)
  }
  extinct_select <- data.frame(extinct,activity=genus_activity[match(names(extinct),names(genus_activity))])
  activity_glm<-glm(extinct~activity,data=extinct_select)
  log_odds_ratio[i]<-summary(activity_glm)$coefficients[2]
  log_odds_error[i]<-summary(activity_glm)$coefficients[4]
  log_odds_p[i]<-summary(activity_glm)$coefficients[8]
}

activity_results_bc <- data.frame(age=bin_top[2:(length(bin_top)-1)], log_odds_ratio, log_odds_error, log_odds_p)


log_odds_ratio<-numeric(0)
log_odds_error<-numeric(0)
log_odds_p<-numeric(0)

#logistic regression based on modified three-timer extinction (three bin moving window)
for (i in 1:(ncol(genus_matrix)-2)) {
  moving_window<-genus_matrix[,seq(i,i+2)]
  moving_window<-subset(moving_window,moving_window[,1]>0 & moving_window[,2]>0) #examines cohort of taxa present in both time t-1 and time t
  extinct<-apply(moving_window,1,function(x) ifelse(x[3]>0,1,0)) #modified version of three-timer extinction; two-timers=extinct, three-timers=survive
  extinct_select<-data.frame(extinct,activity=genus_activity[match(names(extinct),names(genus_activity))])
  activity_glm<-glm(extinct~activity,data=extinct_select)
  log_odds_ratio[i]<-summary(activity_glm)$coefficients[2]
  log_odds_error[i]<-summary(activity_glm)$coefficients[4]
  log_odds_p[i]<-summary(activity_glm)$coefficients[8]
}

activity_results_3t <- data.frame(age=bin_top[2:(length(bin_top)-1)], log_odds_ratio, log_odds_error, log_odds_p)


#prepares info for plot (axis limits and colors)
pdf("selectivity.pdf",width=8,height=5.5,pointsize=12)
par(mgp=c(2,0.75,0))
par(mfrow=c(2,1))
par(mar=c(3, 3, 1, 1))

ylim_max <- max(activity_results_3t$log_odds_ratio + 1.96*activity_results_3t$log_odds_error)
ylim_min <- min(activity_results_3t$log_odds_ratio - 1.96*activity_results_3t$log_odds_error)
height_3t <- ylim_max - ylim_min

activity_results_3t$color <- ifelse(activity_results_3t$log_odds_p<0.05, ifelse(activity_results_3t$log_odds_ratio>0, "red", "blue"),"gray")

plot(activity_results_3t$age,activity_results_3t$log_odds_ratio,ylim=c(ylim_min*1.2,ylim_max),xlim=rev(range(activity_results_3t$age)),xlab="Age (Ma)",ylab="Log odds ratio",col=activity_results_3t$color,pch=16,bty="n")
abline(h=0,lty=3)
segments(activity_results_3t$age,activity_results_3t$log_odds_ratio-1.96*activity_results_3t$log_odds_error,activity_results_3t$age,activity_results_3t$log_odds_ratio+1.96*activity_results_3t$log_odds_error,col=activity_results_3t$color)

segments(252.17, 1.2*ylim_min+0.13, 252.17, -0.05)
text(252.17, 1.2*ylim_min+0.1, "End-Permian", cex=0.6)
segments(201.3, 1.2*ylim_min+0.1, 201.3, -0.05)
text(201.3, 1.2*ylim_min+0.07, "End-Triassic", cex=0.6)
segments(259.9, 1.2*ylim_min+0.1, 259.9, -0.05, lty=2)
text(259.9, 1.2*ylim_min+0.07, "Guadalupian", cex=0.6)
segments(174.1, 1.2*ylim_min+0.1, 174.1, -0.05, lty=2)
text(174.1, 1.2*ylim_min+0.07, "Toarcian", cex=0.6)

#adds geological timescale to plot
rect(stages$max_ma,ylim_min*1.2,stages$min_ma,1.2*ylim_min+0.05,col=paste(stages$color))
text(bin_midpt,ylim_min*1.2+.025,strtrim(stages$interval_name,1),cex=0.6)

mtext(font=2, "(a)", adj=0)


ylim_max <- max(activity_results_bc$log_odds_ratio + 1.96*activity_results_bc$log_odds_error)
ylim_min <- min(activity_results_bc$log_odds_ratio - 1.96*activity_results_bc$log_odds_error)
height_ratio <- (ylim_max-ylim_min)/height_3t

activity_results_bc$color <- ifelse(activity_results_bc$log_odds_p<0.05, ifelse(activity_results_bc$log_odds_ratio>0, "red", "blue"),"gray")

plot(activity_results_bc$age,activity_results_bc$log_odds_ratio,ylim=c(ylim_min*1.2,ylim_max),xlim=rev(range(activity_results_bc$age)),xlab="Age (Ma)",ylab="Log odds ratio",col=activity_results_bc$color,pch=16,bty="n")
abline(h=0,lty=3)
segments(activity_results_bc$age,activity_results_bc$log_odds_ratio-1.96*activity_results_bc$log_odds_error,activity_results_bc$age,activity_results_bc$log_odds_ratio+1.96*activity_results_bc$log_odds_error,col=activity_results_bc$color)

segments(252.17, 1.2*ylim_min+0.11, 252.17, 0)
text(252.17, 1.2*ylim_min+0.08, "End-Permian", cex=0.6)
segments(201.3, 1.2*ylim_min+0.08, 201.3, 0)
text(201.3, 1.2*ylim_min+0.05, "End-Triassic", cex=0.6)
segments(259.9, 1.2*ylim_min+0.08, 259.9, 0, lty=2)
text(259.9, 1.2*ylim_min+0.05, "Guadalupian", cex=0.6)
segments(174.1, 1.2*ylim_min+0.08, 174.1, 0, lty=2)
text(174.1, 1.2*ylim_min+0.05, "Toarcian", cex=0.6)

#adds geological timescale to plot
rect(stages$max_ma,ylim_min*1.2,stages$min_ma,1.2*ylim_min+0.05*height_ratio,col=paste(stages$color))
text(bin_midpt,ylim_min*1.2+0.025*height_ratio,strtrim(stages$interval_name,1),cex=0.6)

mtext(font=2, "(b)", adj=0)

dev.off()



