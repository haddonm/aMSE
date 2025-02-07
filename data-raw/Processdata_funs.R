
## outExtract ####
#' read a selection of aMSE out files, extract required details and return as a list object
#'
#' @param files  object containing directory listing
#' @param pickfiles vector of files to extract
#'
#' @return list of all objects required
#' @export
#'
#' @examples files=files; pickfiles=pickfiles[4]
outExtract <- function(files, pickfiles) {
  files2 <- files[pickfiles]
  nfile <- length(pickfiles)
  label <- vector(mode = "character", length = nfile)
  for (i in 1:nfile) label[i] <- unlist(strsplit(files2[i], ".", fixed = TRUE))[1]
  
  ans <- makelist(label) # vector(mode="list",length=nfile)
  dyn <- makelist(label) # vector(mode="list",length=nfile)
  glbc <- makelist(label) # vector(mode="list",length=nfile)
  prods <- makelist(label) # vector(mode="list",length=nfile)
  #hcrout <- makelist(label) # vector(mode="list",length=nfile)
  scenes <- vector(mode = "character", length = nfile)
  hsargs <- makelist(label)
  outhcr <- makelist(label) # vector(mode="list",length=nfile)
  scoremed <- makelist(label) # vector(mode="list",length=nfile)
  zone <- makelist(label)
  for (i in 1:nfile) {
    # i = 1
    filename <- paste0(outdir, files2[i])
    out <- NULL # so the function knows an 'out' exists
    load(filename)
    ans[[i]] <- out
    dyn[[i]] <- out$sauout
    glbc[[i]] <- out$glb
    prods[[i]] <- out$sauprod
    #hcrout[[i]] <- out$hcrout$refpts
    hsargs[[i]] <- out$hsargs
    scenes[i] <- out$ctrl$runlabel
    outhcr[[i]] <- out$outhcr
    scoremed[[i]] <- out$scoremed
    zone[[i]] <- out$outzone
  }
  nscenes <- length(scenes)
  catch <- makelist(scenes)
  acatch <- makelist(scenes)
  cpue <- makelist(scenes)
  matureB <- makelist(scenes)
  exploitB <- makelist(scenes)
  harvestR <- makelist(scenes)
  deplsB <- makelist(label)
  depleB <- makelist(label)
  medsc <- makelist(scenes)
  cetarg <- makelist(scenes)
  medtarg <- makelist(scenes)
    # j <- 1
    for (j in 1:nscenes) {
    catch[[j]] <- dyn[[j]]$catch
    acatch[[j]] <- dyn[[j]]$acatch
    cpue[[j]] <- dyn[[j]]$cpue
    matureB[[j]] <- dyn[[j]]$matureB
    exploitB[[j]] <- dyn[[j]]$exploitB
    harvestR[[j]] <- dyn[[j]]$harvestR
    deplsB[[j]] <- dyn[[j]]$deplsB
    depleB[[j]] <- dyn[[j]]$depleB
    medtarg[[j]] <- scoremed[[j]]   
    cetarg[[j]] <- outhcr[[j]]$cetarg  
  }
  
  extracted  <-
    list(
      ans = ans,
      dyn = dyn,
      glbc = glbc,
      prods = prods,
      #hcrout = hcrout,
      hsargs = hsargs,
      scenes = scenes,
      nscenes = nscenes,
      outhcr = outhcr,
      zone = zone,
      catch = catch,
      acatch = acatch,
      cpue = cpue,
      matureB = matureB,
      exploitB = exploitB,
      harvestR = harvestR,
      deplsB = deplsB,
      depleB = depleB,
      medsc = medsc,
      cetarg = cetarg,
      scoremed = scoremed
    )
  return(extracted)
}



## process list of out objects

## Catch ####
## - loop through list of scenes and convert to long-form dataframe

#' Title
#'
#' @param scenes 
#' @param basescen 
#' @param nscen 
#' @param catch 
#'
#' @return
#' @export
#'
#' @examples inscenes = scenes;  incatch = catch
catchExtract <- function(inscenes, incatch) {

  nscen <- length(inscenes)
  for (i in 1:nscen) {
    # i = 1
    
    df <- array2df(incatch[[i]])
    colnames(df) <- c("catch", "year", "sau", "run")
    df$scene <- inscenes[i]
    
    df <- df %>% within({
      year <- as.numeric(as.character(year))
    })
    
    if (exists("catch_long")) {
      catch_long <- rbind(catch_long, df)
    } else {
      catch_long <- df
    }
    
    
  } # End compile long form dataframe for MSE catch
  
  
  catch_long <-
    catch_long %>%   within({
      scene_f <- factor(scene, levels = inscenes)
      sau_f <-
        factor(sau,
               levels = c(
                 "sau6",
                 "sau7",
                 "sau8",
                 "sau9",
                 "sau10",
                 "sau11",
                 "sau12",
                 "sau13"
               ))
    })
  
  return(catch_long)
  
} # end function to extract catch in long format 



##---------------------------------------------------------------------------##
## CPUE ####

#'  Loop through list of scenes and convert to long-form dataframe
#'
#' @param inscenes 
#' @param inbasescen 
#' @param incpue 
#'
#' @return
#' @export
#'
#' @examples inscenes = scenes; incpue = cpue
cpueExtract <- function(inscenes, incpue) {
  nscen <- length(inscenes)
  
  for (i in 1:nscen) {
    # i = 1
    
    df <- array2df(incpue[[i]])
    colnames(df) <- c("cpue", "year", "sau", "run")
    df$scene <- inscenes[i]

    df <- df %>% within({
      year <- as.numeric(as.character(year))
    })

    
    if (exists("cpue_long")) {
      cpue_long <- rbind(cpue_long, df)
    } else {
      cpue_long <- df
    }

    
  } # End compile long form dataframe for MSE cpue
  
  cpue_long <-
    cpue_long %>%   within({
      scene_f <- factor(scene, levels = inscenes)
      sau_f <-
        factor(sau,
               levels = c(
                 "sau6",
                 "sau7",
                 "sau8",
                 "sau9",
                 "sau10",
                 "sau11",
                 "sau12",
                 "sau13"
               ))
    })
  
  return(cpue_long)
} # end function to extract cpue in long format 


##----------------------------------------------------------------##
## Exploitable Biomass ####
## loop through list of scenes and convert to long-form dataframe


#' Title
#'
#' @param inscenes 
#' @param inexploitB
#'
#' @return
#' @export
#'
#' @examples inscenes = scenes; inexploitB = exploitB
expBExtract <- function(inscenes, inexploitB) {
  nscen <- length(inscenes)
  for (i in 1:nscen) {
    # i = 1
    
    df <- array2df(inexploitB[[i]])
    colnames(df) <- c("expB", "year", "sau", "run")
    df$scene <- inscenes[i]
    df <- df %>% within({
      year <- as.numeric(as.character(year))
    })
    
    if (exists("explB_long")) {
      explB_long <- rbind(explB_long, df)
    } else {
      explB_long <- df
    }
    
    
  } # End compile long form dataframe for MSE expB
  
  explB_long <-
    explB_long %>%   within({
      scene_f <- factor(scene, levels = inscenes)
      sau_f <-
        factor(sau,
               levels = c(
                 "sau6",
                 "sau7",
                 "sau8",
                 "sau9",
                 "sau10",
                 "sau11",
                 "sau12",
                 "sau13"
               ))
    })
  
  return(explB_long)
} # end function to extract expB in long format 



##----------------------------------------------------------------##
## Mature Biomass ####
## loop through list of scenes and convert to long-form dataframe


#' Title
#'
#' @param inscenes 
#' @param inexploitB
#'
#' @return
#' @export
#'
#' @examples inscenes = scenes; inmatureB = matureB
matBExtract <- function(inscenes, inmatureB) {
  nscen <- length(inscenes)
  for (i in 1:nscen) {
    # i = 1
    
    df <- array2df(inmatureB[[i]])
    colnames(df) <- c("matB", "year", "sau", "run")
    df$scene <- inscenes[i]
    df <- df %>% within({
      year <- as.numeric(as.character(year))
    })
    
    if (exists("matB_long")) {
      matB_long <- rbind(matB_long, df)
    } else {
      matB_long <- df
    }
    
    
  } # End compile long form dataframe for MSE HarvestR
  
  matB_long <-
    matB_long %>%   within({
      scene_f <- factor(scene, levels = inscenes)
      sau_f <-
        factor(sau,
               levels = c(
                 "sau6",
                 "sau7",
                 "sau8",
                 "sau9",
                 "sau10",
                 "sau11",
                 "sau12",
                 "sau13"
               ))
    })
  
  return(matB_long)
} # end function to extract hr in long format 



##----------------------------------------------------------------##
## Depletion Exploitable Biomass ####
## loop through list of scenes and convert to long-form dataframe


#' Title
#'
#' @param inscenes 
#' @param indepleB
#'
#' @return
#' @export
#'
#' @examples inscenes = scenes; indepleB = depleB
depleBExtract <- function(inscenes, indepleB) {
  nscen <- length(inscenes)
  for (i in 1:nscen) {
    # i = 3
    
    df <- array2df(indepleB[[i]])
    colnames(df) <- c("depleB", "year", "sau", "run")
    df$scene <- inscenes[i]
    df <- df %>% within({
      year <- as.numeric(as.character(year))
    })
    
    if (exists("depleB_long")) {
      depleB_long <- rbind(depleB_long, df)
    } else {
      depleB_long <- df
    }
    
    
  } # End compile long form dataframe for MSE HarvestR
  
  depleB_long <-
    depleB_long %>%   within({
      scene_f <- factor(scene, levels = inscenes)
      sau_f <-
        factor(sau,
               levels = c(
                 "sau6",
                 "sau7",
                 "sau8",
                 "sau9",
                 "sau10",
                 "sau11",
                 "sau12",
                 "sau13"
               ))
    })
  
  return(depleB_long)
} # end function to extract hr in long format 



##----------------------------------------------------------------##
## Harvest rate ####
## loop through list of scenes and convert to long-form dataframe

#' Title
#'
#' @param inscenes 
#' @param inharvestR 
#'
#' @return
#' @export
#'
#' @examples inscenes = scenes; inharvestR = harvestR
hrExtract <- function(inscenes, inharvestR) {
  nscen <- length(inscenes)
  for (i in 1:nscen) {
    # i = 1
    
    df <- array2df(inharvestR[[i]])
    colnames(df) <- c("hrate", "year", "sau", "run")
    df$scene <- inscenes[i]
    df <- df %>% within({
      year <- as.numeric(as.character(year))
    })
    
    if (exists("hrate_long")) {
      hrate_long <- rbind(hrate_long, df)
    } else {
      hrate_long <- df
    }
    
    
  } # End compile long form dataframe for MSE HarvestR
  
  hrate_long <-
    hrate_long %>%   within({
      scene_f <- factor(scene, levels = inscenes)
      sau_f <-
        factor(sau,
               levels = c(
                 "sau6",
                 "sau7",
                 "sau8",
                 "sau9",
                 "sau10",
                 "sau11",
                 "sau12",
                 "sau13"
               ))
    })
  
  return(hrate_long)
} # end function to extract hr in long format 



##----------------------------------------------------------------##
## ceTarg  ####
## loop through list of scenes and convert to long-form dataframe


#' Title
#'
#' @param inscenes 
#' @param incatch 
#'
#' @return
#' @export
#'
#' @examples inscenes=scenes; inScores=outhcr
ceTargExtract <- function(inscenes, inceTarg) {
  if (exists("ceTarg_long"))
    rm(ceTarg_long)
  
  
  nscen <- length(inscenes)
  for (i in 1:nscen) {
    # i = 1
    
    df <- array2df(inceTarg[[i]]$cetarg)
    colnames(df) <- c("ceTarg", "year", "sau", "run")
    df$scene <- inscenes[i]
    
    df <- df %>% within({
      year <- as.numeric(as.character(year))
    })
    
    if (exists("ceTarg_long")) {
      ceTarg_long <- rbind(ceTarg_long, df)
    } else {
      ceTarg_long <- df
    }
    
    
  } # End compile long form dataframe for MSE catch
  
  
  ceTarg_long <-
    ceTarg_long %>% within({
      scene_f <- factor(scene, levels = inscenes)
      sau_f <-
        factor(sau,
               levels = c(
                 "sau6",
                 "sau7",
                 "sau8",
                 "sau9",
                 "sau10",
                 "sau11",
                 "sau12",
                 "sau13"
               ))
    })
  
  
  return(ceTarg_long)
  
} # end function to extract catch in long format 


#ceTarg_long <- ceTargExtract(scenes, outhcr)



##----------------------------------------------------------------##
## HCR scores  ####
## loop through list of scenes and convert to long-form dataframe

#' Title
#'
#' @param inscenes 
#' @param incatch 
#'
#' @return
#' @export
#'
#' @examples inscenes=scenes; inScores=outhcr
hcrScoreExtract <- function(inscenes, inScores) {
  if (exists("hcrScore_long"))
    rm(hcrScore_long)
  
  
  nscen <- length(inscenes)
  for (i in 1:nscen) {
    # i = 1
    
    df <- array2df(inScores[[i]]$catchmult)
    colnames(df) <- c("catchmult", "year", "sau", "run")
    df$scene <- inscenes[i]
    
    df <- df %>% within({
      year <- as.numeric(as.character(year))
    })
    
    if (exists("hcrScore_long")) {
      hcrScore_long <- rbind(hcrScore_long, df)
    } else {
      hcrScore_long <- df
    }
    
    
  } # End compile long form dataframe for MSE catch
  
  
  hcrScore_long <-
    hcrScore_long %>% within({
      scene_f <- factor(scene, levels = inscenes)
      sau_f <-
        factor(sau,
               levels = c(
                 "sau6",
                 "sau7",
                 "sau8",
                 "sau9",
                 "sau10",
                 "sau11",
                 "sau12",
                 "sau13"
               ))
    })
  
  
  return(hcrScore_long)
  
} # end function to extract HCR scores & mult in long format 




