


##---------------------------------------------------------------------------##
## catch plots - projection period ####
## catch - loop through list of scenes and convert to long-form dataframe

if (exists("catch_long"))
  rm(catch_long)

catch_long <- catchExtract(scenes, catch)

if (exists("scene_order")) {
  catch_long <- catch_long %>% within({
    scene_f <- factor(scene, levels = scene_order)
  })
}

catchplotdat <- catch_long %>% filter(sau == SAUnames[1] & year %in% c(2025, 2030, 2050))  

catchplot <-  catchplotdat %>%
  ggplot(aes(x=scene_f, y=catch, fill = scene_f)) +
  geom_boxplot(outlier.colour = "orange", width = 0.5, color = "grey") +
  theme_bw() +
  facet_grid(~ year,  scales = "free") +
  scale_x_discrete(breaks = NULL) +
  guides(fill=guide_legend(title="Scenario")) +
  theme(legend.text = element_text(size = 10)) + 
  scale_fill_manual(name = "scene_f", values = myColors)

catchplots <- makelist(SAUnames)
for (i in 1:nSAU) {
  # i <- 1
  catchplot$data <-
    catch_long %>%
    filter(sau == SAUnames[i]  & year %in% c(2025, 2030, 2050))
  catchplots[[i]] <- catchplot
}

catchplot_no_outlier <-  catchplotdat %>%
  ggplot(aes(x=scene_f, y=catch, fill = scene_f)) +
  geom_boxplot(width = 0.5, color = "grey", outliers = FALSE) +
  theme_bw() +
  facet_grid(~ year,  scales = "free") +
  scale_x_discrete(breaks = NULL) +
  guides(fill=guide_legend(title="Scenario")) +
  theme(legend.text = element_text(size = 10)) + 
  scale_fill_manual(name = "scene_f", values = myColors)

catchplot_no_outliers <- makelist(SAUnames)
for (i in 1:nSAU) {
  # i <- 1
  catchplot_no_outlier$data <-
    catch_long %>%
    filter(sau == SAUnames[i]  & year %in% c(2025, 2030, 2050))
  catchplot_no_outliers[[i]] <- catchplot_no_outlier
}



##---------------------------------------------------------------------------##
## cpue - loop through ist of scenes and convert to long-form dataframe
## CPUE plots - projection period ####
if (exists("cpue_long"))
  rm(cpue_long)

cpue_long <- cpueExtract(scenes, cpue)

if (exists("scene_order")) {
  cpue_long <- cpue_long %>% within({
    scene_f <- factor(scene, levels = scene_order)
  })
}


##-------## cpue at 2025, 2030, 2050
cpueplotdat <- cpue_long %>% filter(sau == SAUnames[1] & year %in% c(2025, 2030, 2050))
#cpueplotdat$scene_f = factor(cpueplotdat$scene, levels=c("5years","10years"))

cpueplot <-  cpueplotdat %>%
  ggplot(aes(x=scene_f, y=cpue, fill = scene_f)) +
  geom_boxplot(outlier.colour = "orange", width = 0.5, color = "grey") +
  theme_bw() +
  facet_grid(~ year,  scales = "free") +
  scale_x_discrete(breaks = NULL) +
  guides(fill=guide_legend(title="Scenario")) +
  theme(legend.text = element_text(size = 10)) + 
  scale_fill_manual(name = "scene_f", values = myColors)

cpueplots <- makelist(SAUnames)
for (i in 1:nSAU) {
  # i <- 1
  cpueplot$data <-
    cpue_long %>%
    filter(sau == SAUnames[i]  & year %in% c(2025, 2030, 2050))
  cpueplots[[i]] <- cpueplot
}

cpueplot_no_outlier <-  cpueplotdat %>%
  ggplot(aes(x=scene_f, y=cpue, fill = scene_f)) +
  geom_boxplot(width = 0.5, color = "grey", outliers = FALSE) +
  theme_bw() +
  facet_grid(~ year,  scales = "free") +
  scale_x_discrete(breaks = NULL) +
  guides(fill=guide_legend(title="Scenario")) +
  theme(legend.text = element_text(size = 10)) + 
  scale_fill_manual(name = "scene_f", values = myColors)

cpueplot_no_outliers <- makelist(SAUnames)
for (i in 1:nSAU) {
  # i <- 1
  cpueplot_no_outlier$data <-
    cpue_long %>%
    filter(sau == SAUnames[i]  & year %in% c(2025, 2030, 2050))
  cpueplot_no_outliers[[i]] <- cpueplot_no_outlier
}


##---------------------------------------------------------------------------##
## Exploitable Biomass - loop through ist of scenes and convert to long-form dataframe
## expB plots -projection period ####
if (exists("expB_long"))
  rm(expB_long)
expB_long <- expBExtract(scenes, exploitB)

if (exists("scene_order")) {
  expB_long <- expB_long %>% within({
    scene_f <- factor(scene, levels = scene_order)
  })
}


##-------## cpue at 2025, 2030, 2050
expBplotdat <- expB_long %>% filter(sau == SAUnames[1] & year %in% c(2025, 2030, 2050))
#expBplotdat$scene_f = factor(expBplotdat$scene, levels=c("5years","10years"))

expBplot <-  expBplotdat %>%
  ggplot(aes(x=scene_f, y=expB, fill = scene_f)) +
  geom_boxplot(outlier.colour = "orange", width = 0.5, color = "grey") +
  theme_bw() +
  facet_grid(~ year,  scales = "free") +
  scale_x_discrete(breaks = NULL) +
  guides(fill=guide_legend(title="Scenario")) +
  theme(legend.text = element_text(size = 10)) + 
  scale_fill_manual(name = "scene_f", values = myColors)


expBplots <- makelist(SAUnames)
for (i in 1:nSAU) {
  # i <- 1
  expBplot$data <-
    expB_long %>%
    filter(sau == SAUnames[i]  & year %in% c(2025, 2030, 2050))
  expBplots[[i]] <- expBplot
}


expBplot_no_outlier <-  expBplotdat %>%
  ggplot(aes(x=scene_f, y=expB, fill = scene_f)) +
  geom_boxplot(width = 0.5, color = "grey", outliers = FALSE) +
  theme_bw() +
  facet_grid(~ year,  scales = "free") +
  scale_x_discrete(breaks = NULL) +
  guides(fill=guide_legend(title="Scenario")) +
  theme(legend.text = element_text(size = 10)) + 
  scale_fill_manual(name = "scene_f", values = myColors)


expBplot_no_outliers <- makelist(SAUnames)
for (i in 1:nSAU) {
  # i <- 1
  expBplot_no_outlier$data <-
    expB_long %>%
    filter(sau == SAUnames[i]  & year %in% c(2025, 2030, 2050))
  expBplot_no_outliers[[i]] <- expBplot_no_outlier
}



##---------------------------------------------------------------------------##
## Mature Biomass - loop through list of scenes and convert to long-form dataframe
## matB plots - projection period ####
if (exists("matB_long"))
  rm(matB_long)

matB_long <- matBExtract(scenes, matureB)

if (exists("scene_order")) {
  matB_long <- matB_long %>% within({
    scene_f <- factor(scene, levels = scene_order)
  })
}

##-------## cpue at 2025, 2030, 2050
matBplotdat <- matB_long %>% filter(sau == SAUnames[1] & year %in% c(2025, 2030, 2050))
#matBplotdat$scene_f <- factor(matBplotdat$scene, levels=c("5years","10years"))

matBplot <-  matBplotdat %>%
  ggplot(aes(x=scene_f, y=matB, fill = scene_f)) +
  geom_boxplot(outlier.colour = "orange", width = 0.5, color = "grey") +
  theme_bw() +
  facet_grid(~ year,  scales = "free") +
  scale_x_discrete(breaks = NULL) +
  guides(fill=guide_legend(title="Scenario")) +
  theme(legend.text = element_text(size = 10)) + 
  scale_fill_manual(name = "scene_f", values = myColors)


matBplots <- makelist(SAUnames)
for (i in 1:nSAU) {
  # i <- 1
  matBplot$data <-
    matB_long %>%
    filter(sau == SAUnames[i]  & year %in% c(2025, 2030, 2050))
  matBplots[[i]] <- matBplot
}

matBplot_no_outlier <-  matBplotdat %>%
  ggplot(aes(x=scene_f, y=matB, fill = scene_f)) +
  geom_boxplot(width = 0.5, color = "grey", outliers = FALSE) +
  theme_bw() +
  facet_grid(~ year,  scales = "free") +
  scale_x_discrete(breaks = NULL) +
  guides(fill=guide_legend(title="Scenario")) +
  theme(legend.text = element_text(size = 10)) + 
  scale_fill_manual(name = "scene_f", values = myColors)


matBplot_no_outliers <- makelist(SAUnames)
for (i in 1:nSAU) {
  # i <- 1
  matBplot_no_outlier$data <-
    matB_long %>%
    filter(sau == SAUnames[i]  & year %in% c(2025, 2030, 2050))
  matBplot_no_outliers[[i]] <- matBplot_no_outlier
}



##---------------------------------------------------------------------------##
## Exploitable Biomass depletion - loop through list of scenes and convert to long-form dataframe
## depleB plots - projection period ####
if (exists("depleB_long"))
  rm(depleB_long)

depleB_long <- depleBExtract(scenes, depleB)

if (exists("scene_order")) {
  depleB_long <- depleB_long %>% within({
    scene_f <- factor(scene, levels = scene_order)
  })
}

##-------## cpue at 2025, 2030, 2050
depleBplotdat <- depleB_long %>% filter(sau == SAUnames[1] & year %in% c(2025, 2030, 2050))
#depleBplotdat$scene_f = factor(depleBplotdat$scene, levels=c("5years","10years"))

depleBplot <-  depleBplotdat %>%
  ggplot(aes(x=scene_f, y=depleB, fill = scene_f)) +
  geom_boxplot(outlier.colour = "orange", width = 0.5, color = "grey") +
  theme_bw() +
  facet_grid(~ year,  scales = "free") +
  scale_x_discrete(breaks = NULL) +
  guides(fill=guide_legend(title="Scenario")) +
  theme(legend.text = element_text(size = 10)) + 
  scale_fill_manual(name = "scene_f", values = myColors)


depleBplots <- makelist(SAUnames)
for (i in 1:nSAU) {
  # i <- 1
  depleBplot$data <-
    depleB_long %>%
    filter(sau == SAUnames[i]  & year %in% c(2025, 2030, 2050))
  depleBplots[[i]] <- depleBplot
}

depleBplot_no_outlier <-  depleBplotdat %>%
  ggplot(aes(x=scene_f, y=depleB, fill = scene_f)) +
  geom_boxplot(width = 0.5, color = "grey", outliers = FALSE) +
  theme_bw() +
  facet_grid(~ year,  scales = "free") +
  scale_x_discrete(breaks = NULL) +
  guides(fill=guide_legend(title="Scenario")) +
  theme(legend.text = element_text(size = 10)) + 
  scale_fill_manual(name = "scene_f", values = myColors)


depleBplot_no_outliers <- makelist(SAUnames)
for (i in 1:nSAU) {
  # i <- 1
  depleBplot_no_outlier$data <-
    depleB_long %>%
    filter(sau == SAUnames[i]  & year %in% c(2025, 2030, 2050))
  depleBplot_no_outliers[[i]] <- depleBplot_no_outlier
}


##----------------------------------------------------------------##
## HarvestR - loop through ist of scenes and convert to long-form dataframe
## Harvest Rate plots - projection period ####

if (exists("hrate_long"))
  rm(hrate_long)
hrate_long <- hrExtract(scenes, harvestR)

if (exists("scene_order")) {
  hrate_long <- hrate_long %>% within({
    scene_f <- factor(scene, levels = scene_order)
  })
}

hrateplotdat <- hrate_long %>% filter(sau == SAUnames[1] & year %in% c(2025, 2030, 2050))

hrateplot <-  hrateplotdat %>%
  ggplot(aes(x=scene_f, y=hrate, fill = scene_f)) +
  geom_boxplot(outlier.colour = "orange", width = 0.5, color = "grey") +
  theme_bw() +
  facet_grid(~ year,  scales = "free") +
  scale_x_discrete(breaks = NULL) +
  guides(fill=guide_legend(title="Scenario")) +
  theme(legend.text = element_text(size = 10)) + 
  scale_fill_manual(name = "scene_f", values = myColors)

hrplots <- makelist(SAUnames)
for (i in 1:nSAU) {
  # i <- 1
  hrateplot$data <-
    hrate_long %>%
    filter(sau == SAUnames[i]  & year %in% c(2025, 2030, 2050))
  hrplots[[i]] <- hrateplot
}


hrateplot_no_outlier <-  hrateplotdat %>%
  ggplot(aes(x=scene_f, y=hrate, fill = scene_f)) +
  geom_boxplot(width = 0.5, color = "grey", outliers = FALSE) +
  theme_bw() +
  facet_grid(~ year,  scales = "free") +
  scale_x_discrete(breaks = NULL) +
  guides(fill=guide_legend(title="Scenario")) +
  theme(legend.text = element_text(size = 10)) + 
  scale_fill_manual(name = "scene_f", values = myColors)

hrateplot_no_outliers <- makelist(SAUnames)
for (i in 1:nSAU) {
  # i <- 1
  hrateplot_no_outlier$data <-
    hrate_long %>%
    filter(sau == SAUnames[i]  & year %in% c(2025, 2030, 2050))
  hrateplot_no_outliers[[i]] <- hrateplot_no_outlier
}



##----------------------------------------------------------------##
## hcrScore - loop through list of scenes and convert to long-form dataframe
## HCR Score plots - projection period ####

if (exists("hcrScore_long"))
  rm(hcrScore_long)

hcrScore_long <- hcrScoreExtract(scenes, outhcr)

if (exists("scene_order")) {
  hcrScore_long <- hcrScore_long %>% within({
    scene_f <- factor(scene, levels = scene_order)
  })
}


hscoreplotdat <- hcrScore_long %>% filter(sau == SAUnames[1] & year %in% c(seq(2023,2028,1)))

hcrscoreplot <-  hscoreplotdat %>%
  ggplot(aes(x=scene_f, y=catchmult, fill = scene_f)) +
  geom_boxplot(outlier.colour = "orange", width = 0.5, color = "grey") +
  theme_bw() +
  facet_grid(~ year,  scales = "free") +
  scale_x_discrete(breaks = NULL) +
  guides(fill=guide_legend(title="Scenario")) +
  theme(legend.text = element_text(size = 10)) + 
  scale_fill_manual(name = "scene_f", values = myColors)  +
  theme(legend.position="bottom")

hcrscoreplots <- makelist(SAUnames)
for (i in 1:nSAU) {
  # i <- 1
  hcrscoreplot$data <-
    hcrScore_long %>%
    filter(sau == SAUnames[i]  & year %in% c(seq(2023,2028,1)))
  hcrscoreplots[[i]] <- hcrscoreplot
}

hcrscoreplot_no_outlier <-  hscoreplotdat %>%
  ggplot(aes(x=scene_f, y=catchmult, fill = scene_f)) +
  geom_boxplot(width = 0.5, color = "grey", outliers = FALSE) +
  theme_bw() +
  facet_grid(~ year,  scales = "free") +
  scale_x_discrete(breaks = NULL) +
  guides(fill=guide_legend(title="Scenario")) +
  theme(legend.text = element_text(size = 10)) + 
  scale_fill_manual(name = "scene_f", values = myColors)  +
  theme(legend.position="bottom")

hcrscoreplot_no_outliers <- makelist(SAUnames)
for (i in 1:nSAU) {
  # i <- 1
  hcrscoreplot_no_outlier$data <-
    hcrScore_long %>%
    filter(sau == SAUnames[i]  & year %in% c(seq(2023,2028,1)))
  hcrscoreplot_no_outliers[[i]] <- hcrscoreplot_no_outlier
}




##---------------------------------------------------------------------------##
## year CPUE hits target ####
## Finds the first year in each run where CPUE exceeds the CPUE Target

if (exists("cpue_long_targ"))
  rm(cpue_long_targ)

cpue_targ <- ceTargExtract(scenes, outhcr)

if (exists("scene_order")) {
  cpue_long <- cpue_long %>% within({
    scene_f <- factor(scene, levels = scene_order)
  })
}

# cpue_targ <- data.frame(hcrout[[1]][,2])
# cpue_targ$sau <- row.names(cpue_targ)
# colnames(cpue_targ) <- c("target", "sau")

#cpue_long$sau <- as.character(cpue_long$sau)
cpue_long_targ <- filter(cpue_long, year >= 2020) %>% left_join(select(cpue_targ, -sau_f, -scene_f), by = c("scene", "sau",  "year", "run") ) 

# yrstoTarg <- cpue_long_targ %>% group_by(scene_f, sau, run) %>%
#   summarise(firstyr = min(which(cpue > target, arr.ind = TRUE), na.rm=T)) %>% 
#   filter(!is.infinite(firstyr))

dt <- data.table(cpue_long_targ)    
setkey(dt, scene_f, sau, run)
yrstoTarg <- dt[cpue > ceTarg, .SD[1], by=key(dt)] %>% as.data.frame()
rm(dt)


yrstoTarg <-
  yrstoTarg %>%   within({
    sau_f <-
      factor(
        sau,
        levels = c(
          "sau6","sau7","sau8","sau9","sau10","sau11","sau12","sau13")
      )
  })

##---------------------------------------------------------------------------##
## year of max CPUE ####
## Finds the first year in each run where catch exceeds MSY 


yrstoMaxCPUE <- cpue_long %>% 
  filter(year > 2020) %>% 
  group_by(scene_f, sau, run) %>%
  summarise(firstyr = min(which(cpue == max(cpue), arr.ind = TRUE), na.rm=T)) %>% 
  filter(!is.infinite(firstyr))

yrstoMaxCPUE <-
  yrstoMaxCPUE %>%   within({
    sau_f <-
      factor(
        sau,
        levels = c(
          "sau6",
          "sau7",
          "sau8",
          "sau9",
          "sau10",
          "sau11",
          "sau12",
          "sau13"
        )
      )
    
  })


##---------------------------------------------------------------------------##
## year catch exceeds MSY ####
## Finds the first year in each run where catch exceeds MSY 

#msy_sau <- data.frame(prods$BaseCase[3,])
msy_sau <- data.frame(prods[[1]][3,])
msy_sau$sau <- row.names(msy_sau)
colnames(msy_sau) <- c("msy", "sau")
msy_sau$sau <- as.character(msy_sau$sau)

catch_long$sau <- as.character(catch_long$sau)
catch_long_msy <- filter(catch_long, year > 2020) %>% left_join( msy_sau, by = c("sau" = "sau") ) 

dt1 <- data.table(catch_long_msy)
setkey(dt1, scene_f, sau, run)
yrstoMSY <- dt1[catch > msy, .SD[1], by=key(dt1)] %>% as.data.frame()
rm(dt1)

# yrstoMSY <- catch_long_msy %>% group_by(scene, sau, run) %>% 
#   arrange(scene, sau, run) %>% 
#   filter(catch >= msy) %>% 
#   filter(rank(year, ties.method="first")==1)


if (exists("scene_order")) {
  yrstoMSY <- yrstoMSY %>% within({
    scene_f <- factor(scene, levels = scene_order)
  })
}


## diagnostics
# pick <- which(tmp$scene == "BaseC_MR" & tmp$sau == "sau6" & as.numeric(as.character(tmp$run))==1)
# test <- tmp[pick,]

yrstoMSY <-
  yrstoMSY %>%   within({
    sau_f <-
      factor(
        sau,
        levels = c(
          "sau6",
          "sau7",
          "sau8",
          "sau9",
          "sau10",
          "sau11",
          "sau12",
          "sau13"
        )
      )
  })

##---------------------------------------------------------------------------##
## year catch exceeds MSY + Harvest Rate ####
## - didnt really  show anything different

msy_hr <- left_join(catch_long_msy, hrate_long) 


dt2 <- data.table(msy_hr)    
setkey(dt2,scene_f, sau, run)
yrstoMSY_hr <- dt2[catch > msy & hrate < 0.2, .SD[1], by=key(dt2)] %>% as.data.frame()
rm(dt2)

yrstoMSY_hr <-
  yrstoMSY_hr %>%   within({
    sau_f <-
      factor(
        sau,
        levels = c(
          "sau6",
          "sau7",
          "sau8",
          "sau9",
          "sau10",
          "sau11",
          "sau12",
          "sau13"
        )
      )
  })



##---------------------------------------------------------------------------##
## years to matB > B-MSY TRP ####
## Finds the first year in which B exceeds B-MSY 

#msy_sau <- data.frame(prods$BaseCase[3,])
Bmsy_sau <- data.frame(prods[[1]][2,])
Bmsy_sau$sau <- row.names(Bmsy_sau)
colnames(Bmsy_sau) <- c("Bmsy", "sau")
Bmsy_sau$sau <- as.character(Bmsy_sau$sau)

matB_long$sau <- as.character(matB_long$sau)
matB_long_Bmsy <- filter(matB_long, year > 2020) %>% left_join( Bmsy_sau, by = c("sau" = "sau") ) 

# yrstoMSY <- matB_long_msy %>% group_by(scene_f, sau, run) %>%
#   summarise(firstyr = min(which(catch > msy, arr.ind = TRUE), na.rm=T)) %>% 
#   filter(!is.infinite(firstyr))

dt1 <- data.table(matB_long_Bmsy)    
setkey(dt1,scene_f, sau, run)
yrstoBMSY <- dt1[matB > Bmsy, .SD[1], by=key(dt1)] %>% as.data.frame()
rm(dt1)


## diagnostics
# pick <- which(tmp$scene == "BaseC_MR" & tmp$sau == "sau6" & as.numeric(as.character(tmp$run))==1)
# test <- tmp[pick,]

yrstoBMSY <-
  yrstoBMSY %>%   within({
    sau_f <-
      factor(
        sau,
        levels = c(
          "sau6",
          "sau7",
          "sau8",
          "sau9",
          "sau10",
          "sau11",
          "sau12",
          "sau13"
        )
      )
  })


## years to matB > B-MSY  Conditional on Harvest Rate ####
## - didnt really  show anything different

Bmsy_hr <- left_join(matB_long_Bmsy, hrate_long) 


dt2 <- data.table(Bmsy_hr)    
setkey(dt2,scene_f, sau, run)
yrstoBMSY_hr <- dt2[matB > Bmsy & hrate < 0.2, .SD[1], by=key(dt2)] %>% as.data.frame()
rm(dt2)

yrstoBMSY_hr <-
  yrstoBMSY_hr %>%   within({
    sau_f <-
      factor(
        sau,
        levels = c(
          "sau6",
          "sau7",
          "sau8",
          "sau9",
          "sau10",
          "sau11",
          "sau12",
          "sau13"
        )
      )
  })




##---------------------------------------------------------------------------##
## years catch exceeds catch Target ####
## Finds the first year in each run where catch exceeds catch (MSY-Proxy)
## e.g. Mean catch of the Ref Period

#msy_sau <- data.frame(prods$BaseCase[3,])
msy_sau <- data.frame(prods[[1]][3,])
msy_sau$sau <- row.names(msy_sau)
colnames(msy_sau) <- c("msy", "sau")
msy_sau$sau <- as.character(msy_sau$sau)

catch_long$sau <- as.character(catch_long$sau)
catch_long_msy <- filter(catch_long, year > 2020) %>% left_join( msy_sau, by = c("sau" = "sau") ) 

# yrstoMSY <- catch_long_msy %>% group_by(scene_f, sau, run) %>%
#   summarise(firstyr = min(which(catch > msy, arr.ind = TRUE), na.rm=T)) %>% 
#   filter(!is.infinite(firstyr))

dt1 <- data.table(catch_long_msy)    
setkey(dt1,scene_f, sau, run)
yrstoMSY <- dt1[catch > msy, .SD[1], by=key(dt1)] %>% as.data.frame()
rm(dt1)




## diagnostics
# pick <- which(tmp$scene == "BaseC_MR" & tmp$sau == "sau6" & as.numeric(as.character(tmp$run))==1)
# test <- tmp[pick,]

yrstoMSY <-
  yrstoMSY %>%   within({
    sau_f <-
      factor(
        sau,
        levels = c(
          "sau6",
          "sau7",
          "sau8",
          "sau9",
          "sau10",
          "sau11",
          "sau12",
          "sau13"
        )
      )
  })

##---------------------------------------------------------------------------##
## YRS to MSY Conditional on Harvest Rate ####
## - didnt really  show anything different

msy_hr <- left_join(catch_long_msy, hrate_long) 


dt2 <- data.table(msy_hr)    
setkey(dt2,scene_f, sau, run)
yrstoMSY_hr <- dt2[catch > msy & hrate < 0.2, .SD[1], by=key(dt2)] %>% as.data.frame()
rm(dt2)

yrstoMSY_hr <-
  yrstoMSY_hr %>%   within({
    sau_f <-
      factor(
        sau,
        levels = c(
          "sau6",
          "sau7",
          "sau8",
          "sau9",
          "sau10",
          "sau11",
          "sau12",
          "sau13"
        )
      )
  })

##---------------------------------------------------------------------------##


##---------------------------------------------------------------------------##

