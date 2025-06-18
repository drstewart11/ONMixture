#'Hierarchical Bayesian models'
#' @param count Data frame with count data including pond_name, year, site, YC_small, BS_small
#' @param mgmt Data frame with management data including pond_name, year, and stocking/removal variables
#' @param hab Data frame with habitat variables including pond_name, year, site, pH, wtemp, doxygen, wcond, veg, wdepth
#' @param species Character string: "YCHUB" or "BSHINER"
#' @export
countmix<-function(count,mgmt,hab,source,species=c("YCHUB","BSHINER")){
  #Error bounds
  if(length(species)>1|missing(species))stop("'species' must contain
                                             only one value",call.=FALSE)
  if(missing(species))stop("must specify species",call.=FALSE)
  if(missing(count))stop("must specify count data",call.=FALSE)
  if(missing(mgmt))stop("must specify mgmt (management) activity
                        data",call.=FALSE)
  if(missing(hab))stop("must specify hab (habitat) activity data",call=FALSE)
  options(warn=-1)


  #Add pname field to count, mgmt, and hab .csv files
  #pname is a numeric wetland pond field
  #yr is a numeric year field

  count <- count
  count$yr<-as.numeric(as.factor(count$year))

  # Number of survey years for each pond
  tablePond <- data.table(count)

  # Count the number of survey years for each pond
  survey_years <- tablePond[, .N, by = .(pond_name, yr)]
  print(survey_years)

  # Count the number of unique survey years for each pond
  survey_years_per_pond <- tablePond[, .N, by = pond_name]

  pond_year_check <- count %>%
    group_by(pond_name, year) %>%
    summarize(num_surveys = n()) %>%
    ungroup()

  min_surveys_per_pond_year <- min(pond_year_check$num_surveys)
  if(min_surveys_per_pond_year == 0){
    stop("Error: Some ponds are missing surveys in certain years.")
  }

  pond_years <- count %>%
    group_by(pond_name) %>%
    summarize(all_years = list(sort(unique(year)))) %>%
    ungroup()

  missing_years_check <- pond_years %>%
    rowwise() %>%
    mutate(has_missing_years = any(diff(unlist(all_years)) > 1))

  if(any(missing_years_check$has_missing_years)){
    print("Warning: Some ponds have missing survey years. Here are the affected ponds:")
    print(missing_years_check %>% filter(has_missing_years == TRUE))
  }

  print("All error checks passes successfully!")

  print(survey_years_per_pond)

  head(count)

  species <- species
  #Reorganize count data by Site, Wetland Pond, Year, and Species
  if(species=="YCHUB"){
    newdat<-data.frame(pondname = count$pond_name,year=count$yr,day=count$day,
                       site=count$site,y=count$YCHUB,yrlab = count$year)
    ponds_to_remove <- c("BenNewPond","BrasherPond","CobblePond","MagoffinPond","MesquitePond","UrquidesPond")
  }else if(species=="BSHINER"){
    newdat<-data.frame(pondname = count$pond_name,year=count$yr,day=count$day,
                       site=count$site,y=count$BSHINER,yrlab = count$year)
    ponds_to_remove <- c("BathHousePond","BigTank","LodgeTank","NHayHollowPond","NMinckleyPond","OasisPond","PhD1Pond",
                         "PhD2Pond","SMinckleyPond","UpperChalkTank")
  }


  newdat$pondname <- as.character(newdat$pondname)
  newdat$pname <- as.numeric(as.factor(newdat$pondname))

  newdat_filtered <- newdat[!newdat$pondname %in% ponds_to_remove,]

  newdat_filtered$pname <- as.numeric(as.factor(newdat_filtered$pondname))

  #Define 4-dimensional array dimensions
  nsite=length(unique(newdat_filtered$site))
  npond=length(unique(newdat_filtered$pname))
  nday=length(unique(newdat_filtered$day))
  nyear=length(unique(newdat_filtered$year))

  site_index <- as.numeric(factor(newdat_filtered$site, levels = unique(newdat_filtered$site)))
  day_index <- as.numeric(factor(newdat_filtered$day, levels = unique(newdat_filtered$day)))
  pond_index <- as.numeric(factor(newdat_filtered$pname, levels = sort(unique(newdat_filtered$pname))))
  year_index <- as.numeric(factor(newdat_filtered$yrlab, levels = unique(newdat_filtered$yrlab)))

  cat("Max site index:", max(site_index), "Expected:", nsite, "\n")
  cat("Max day index:", max(day_index), "Expected:", nday, "\n")
  cat("Max pond index:", max(pond_index), "Expected:", npond, "\n")
  cat("Max year index:", max(year_index), "Expected:", nyear, "\n")

  pond_mapping <- newdat_filtered %>%
    select(pondname) %>%
    distinct() %>%
    mutate(pname = as.numeric(as.factor(pondname))) %>%
    arrange(pname)  # Ensuring correct order


  # Step 2: Filter the Data (Without Changing Indices)
  newdat_filtered <- newdat[!newdat$pondname %in% ponds_to_remove,]

  # Step 3: Explicitly Remap Pond Names AFTER Filtering
  newdat_filtered$pname <- as.numeric(factor(newdat_filtered$pondname, levels = sort(unique(newdat_filtered$pondname))))

  # Step 4: Define pond_index from Corrected `pname`
  pond_index <- newdat_filtered$pname

  y_filtered <- array(NA, dim = c(nsite, nday, npond, nyear))

  for (i in seq_len(nrow(newdat_filtered))) {
    si <- site_index[i]
    di <- day_index[i]
    pi <- pond_index[i]  # Now correctly mapped!
    yi <- year_index[i]

    # Final safety check before assignment
    if (!is.na(si) & !is.na(di) & !is.na(pi) & !is.na(yi) &
        si > 0 & di > 0 & pi > 0 & yi > 0 &
        si <= nsite & di <= nday & pi <= npond & yi <= nyear) {
      y_filtered[si, di, pi, yi] <- newdat_filtered$y[i]
    } else {
      warning(paste("Skipping row", i,
                    "si =", si, "di =", di, "pi =", pi, "yi =", yi,
                    "(out of bounds)"))
    }
  }

  pond_year_summary <- apply(y_filtered, c(3, 4), sum, na.rm = TRUE)

  # Step 7: Convert to Data Frame
  pivot_table <- as.data.frame(pond_year_summary)

  # Step 8: Assign Correct Pond Names
  rownames(pivot_table) <- pond_mapping$pondname

  # Step 9: Rename Year Columns
  colnames(pivot_table) <- sort(unique(newdat_filtered$yrlab))

  # Step 10: Add Grand Total Row
  pivot_table$Grand_Total <- rowSums(pivot_table)

  # Step 11: Convert to Tibble for Proper Formatting
  pivot_table <- pivot_table %>%
    tibble::rownames_to_column(var = "pondname")

  # Step 12: Calculate & Append Grand Total Row
  grand_total <- colSums(select(pivot_table, -pondname))
  grand_total_df <- as.data.frame(t(grand_total)) %>%
    mutate(pondname = "Grand Total")

  # Convert `pondname` to character to avoid binding issues
  grand_total_df$pondname <- as.character(grand_total_df$pondname)

  # Step 13: Bind the Grand Total Row
  colnames(grand_total_df) <- colnames(pivot_table)
  pivot_table$Grand_Total <- as.numeric(pivot_table$Grand_Total)

  # Step 14: Print Final Table
  print(pivot_table)

  print("Initiate Bayesian Population Model. This may take several minutes to hours.",quote=FALSE)

  if(species=="YCHUB"){
    modelFilename="Bayesian.Population.Model.txt"
    cat("
  model{
  phi ~ dunif(0,100)

  for(k in 1:npond){
    beta[k] ~ dnorm(0,0.01)
    K[k] ~ dunif(5,5000)
    eta[k] ~ dgamma(phi,phi)
  }

  for(k in 1:npond){
  for(t in 1:nyear){
  r[k,t] ~ dunif(0, 5)
  lambda.r[k,t] <- exp(r[k,t])
  }
  }

  # Zero-inflation parameters
  alpha_psi ~ dnorm(0, 0.1) # Prior for extinction intercept
  beta_psi ~ dnorm(0, 0.1)  # Prior for effect of past population size

  for(i in 1:nsite){
    for(k in 1:npond){
      N[i, k, 1] ~ dpois(lambda[i, k])
      lambda[i, k] <- muL[i, k] * eta[k]
      log(muL[i, k]) <- beta[k]

      for(t in 2:nyear){
        mu[i, k, t - 1] <- ((lambda.r[k, t] * N[i, k, t-1]) /
                          (1 + ((N[i, k, t-1] * (lambda.r[k, t] - 1)) / K[k])))

        N[i, k, t] ~ dpois(mu[i, k, t - 1])

      }
    }
  }

  for(i in 1:nsite){
    for(j in 1:nrep){
      for(k in 1:npond){
        for(t in 1:nyear){
          y[i, j, k, t] ~ dbin(q[i, j, k, t], N[i, k, t])
          q[i, j, k, t] ~ dbeta(5,9) # Slightly informative prior. South Minckley 1605 in 2022
        }
      }
    }
  }

  for(k in 1:npond){
    for(t in 1:nyear){
      N.total[k, t] <- sum(N[, k, t])
      det.mean[k, t] <- mean(q[,, k, t])
    }
  }

  for(k in 1:npond){
  for(t in 2:nyear){
    r_eff[k, t] <- log((N.total[k, t] + 0.0001) / (N.total[k, t-1] + 0.0001))
  }
}

for(k in 1:npond){
  mu_r[k] <- mean(r_eff[k,2:nyear])
}


}",fill=TRUE,file=modelFilename)
  }else if(species=="BSHINER"){
    modelFilename="Bayesian.Population.Model.txt"
    cat("
      model{
        phi ~ dunif(0,100)

        for(k in 1:npond){
          beta[k] ~ dnorm(0,0.01)
          K[k] ~ dunif(5,5000)
          eta[k] ~ dgamma(phi,phi)
        }

        for(k in 1:npond){
          for(t in 1:nyear){
            r[k,t] ~ dunif(0, 5)
            lambda.r[k,t] <- exp(r[k,t])
          }
        }

        # Zero-inflation parameters
        alpha_psi ~ dnorm(0, 0.1) # Prior for extinction intercept
        beta_psi ~ dnorm(0, 0.1)  # Prior for effect of past population size

        for(i in 1:nsite){
          for(k in 1:npond){
            N[i, k, 1] ~ dpois(lambda[i,k])
            lambda[i, k] <- muL[i, k] * eta[k]
            log(muL[i, k]) <- beta[k]

            for(t in 2:nyear){
              mu[i, k, t - 1] <- ((lambda.r[k, t] * N[i, k, t - 1]) /
                                    (1 + ((N[i, k, t - 1] * (lambda.r[k, t] - 1)) / K[k])))

              N[i, k, t] ~ dpois(mu[i, k, t - 1])

            }
          }
        }

        for(i in 1:nsite){
          for(j in 1:nrep){
            for(k in 1:npond){
              for(t in 1:nyear){
                y[i, j, k, t] ~ dbin(q[i, j, k, t], N[i, k, t])
                q[i, j, k, t] ~ dbeta(6,7) # Slightly informative prior. South Minckley 1605 in 2022
              }
            }
          }
        }

        for(k in 1:npond){
          for(t in 1:nyear){
            N.total[k, t] <- sum(N[, k, t])
            det.mean[k, t] <- mean(q[,, k, t])
          }
        }

        for(k in 1:npond){
          for(t in 2:nyear){
            r_eff[k, t] <- log((N.total[k, t] + 0.0001) / (N.total[k, t-1] + 0.0001))
          }
        }

        for(k in 1:npond){
          mu_r[k] <- mean(r_eff[k,2:nyear])
        }


}",fill=TRUE,file=modelFilename)
  }

#Initial values

Nst <- array(NA, dim = c(nsite, npond, nyear))  # Ensure correct shape

for (i in 1:nsite) {
  for (k in 1:npond) {
    for (t in 1:nyear) {
      Nst[i, k, t] <- sum(y_filtered[i, , k, t], na.rm=TRUE) + 1  # Summing over j (reps)
    }
  }
}

jags.inits=function()list(N=Nst)



#Bundle data
jags.data=list(y=y_filtered,nsite=nsite,nrep=nday,npond=npond,nyear=nyear)

#Parameters monitored
jags.params=c("phi","det.mean","K","mu_r","N.total")

#MCMC settings
nc=5; nt=1; nb=25000; ni=750000

out<-jags(jags.data,jags.inits,parameters.to.save=jags.params,
          model.file=modelFilename,n.chains=nc,
          n.iter=ni,n.burnin=nb,n.thin=nt,parallel=TRUE,verbose=TRUE)


# Create Year labels for each pond
yrlab <- rep(seq(min(newdat_filtered$yrlab), max(newdat_filtered$yrlab), by=1),
             each=length(unique(newdat_filtered$pondname)))

# Create Pond Labels, ensuring proper repetition across years
pond.name <- rep(sort(unique(newdat_filtered$pondname)), times = length(unique(newdat_filtered$yrlab)))

# Extract N.total from JAGS output, ensuring proper vector structure
N.total <- round(as.vector(unlist(out$mean$N.total)))
N.lower <- pmax(round(as.vector(unlist(out$q2.5$N.total))), 0)  # Ensuring non-negative values
N.upper <- round(as.vector(unlist(out$q97.5$N.total)))

# Extract detection probability posteriors
q.mean <- round(as.vector(unlist(out$mean$det.mean)), 2)
q.lower <- round(as.vector(unlist(out$q2.5$det.mean)), 2)
q.upper <- round(as.vector(unlist(out$q97.5$det.mean)), 2)

# Extract mu_r and convert to lambda (finite population growth rate)
r_mean <- as.vector(unlist(out$mean$mu_r))
r_lower <- as.vector(unlist(out$q2.5$mu_r))
r_upper <- as.vector(unlist(out$q97.5$mu_r))

lambda <- exp(r_mean)
lambda.lower <- exp(r_lower)
lambda.upper <- exp(r_upper)

# Ensure lambda stays 1 if bounds are exactly 1
lambda[lambda.lower == 1 & lambda.upper == 1] <- 1

k.mean <- as.vector(unlist(out$mean$K))
k.lower <- as.vector(unlist(out$q2.5$K))
k.upper <- as.vector(unlist(out$q97.5$K))

# Ensure correct year assignments per pond
unique_ponds <- unique(pond.name)
minyr <- min(newdat_filtered$yrlab)
maxyr <- max(newdat_filtered$yrlab)
correct_years <- rep(minyr:maxyr, each=length(unique_ponds))

# Combine results into a structured data frame
NQresults <- data.frame(
  WetlandPond = pond.name,
  Year = correct_years,
  NLower95 = N.lower,
  Pop_estimate = N.total,
  NUpper95 = N.upper,
  KLower95 = k.lower,
  Carrying_capacity = k.mean,
  KUpper95 = k.upper,
  QLower95 = q.lower,
  Detection_estimate = q.mean,
  QUpper95 = q.upper
)

# Ensure results are ordered correctly by pond and year
NQresults <- NQresults[order(NQresults$WetlandPond, NQresults$Year),]

# Display the results
print(NQresults)


#Recreate Wetland pond labels
pond.name2<-as.character(unique(unlist(sort(newdat_filtered$pondname))))
growth_status <- ifelse(lambda.lower > 1, "Growing",
                        ifelse(lambda.upper < 1, "Declining", "Stable"))

#Use data.frame to package lambda and carrying capacity results and
#save to working directory
GRresults<-data.frame(WetlandPond=pond.name2,
                      Inst_Pop_G_Rate = round(lambda, 2),
                      RLower95 = round(lambda.lower, 2),
                      RUpper95 = round(lambda.upper, 2),
                      Growth_status = growth_status)

if(species=="YCHUB"){

  #Capture and Write results to working directory (R Data Files)
  write.csv(NQresults, paste0(source, "/YaquiChubPondAbundanceDetection.csv"), row.names=F)
  write.csv(GRresults, paste0(source, "/YaquiChubPondPopGrowthRate.csv"), row.names=F)

  plot<-ggplot(NQresults,aes(as.factor(Year),Pop_estimate,colour=factor(WetlandPond)))+
    geom_point(size=4)+
    geom_errorbar(aes(ymin=NLower95,ymax=NUpper95),linewidth=1.5)+
    facet_wrap(~WetlandPond,ncol=7, scales = "free_y")+
    guides(colour="none")+
    theme_bw()+
    xlab("Year")+
    ylab("Abundance")+
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult=c(0, 0.01))) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=16),
          strip.text.x=element_text(size=12))
  print(plot)
  ggsave(filename = file.path(source, "/YaquiChubWetlandPondAbundanceFigure.tiff"),plot=plot,
         width=20,height=10,dpi=300)
}else if(species=="BSHINER"){
  #Capture and Write results to working directory (R Data Files)
  write.csv(NQresults,paste0(source,"/BeautifulShinerPondAbundanceDetection.csv"),row.names=F)
  write.csv(GRresults,paste0(source, "/BeautifulShinerPondPopGrowthRate.csv"),row.names=F)

  plot<-ggplot(NQresults,aes(as.factor(Year),Pop_estimate,colour=factor(WetlandPond)))+
    geom_point(size=4)+
    geom_errorbar(aes(ymin=NLower95,ymax=NUpper95),linewidth=1.5)+
    facet_wrap(~WetlandPond,ncol=7,scales="free_y")+
    guides(colour="none")+
    theme_bw()+
    xlab("Year")+
    ylab("Abundance")+
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult=c(0, 0.01))) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=16),
          strip.text.x=element_text(size=12))
  print(plot)
  ggsave(filename = file.path(source, "BeautifulShinerWetlandPondAbundanceFigure.tiff"),plot=plot,
         width=20,height=10,dpi=300)
}

print("Bayesian Population Model Complete",quote=FALSE)



#############################################################################
#############################################################################
#############################################################################


if(species=="YCHUB"){

  #Reorganize habitat data by Site, Wetland Pond, and Year
  newhab<-data.frame(pond_name=hab$pond_name,rep=hab$day,year=hab$year,site=hab$site,
                     pH=hab$pH,wtemp=hab$wtemp,
                     doxygen=hab$doxygen,wcond=hab$wcond,
                     veg=hab$veg,
                     wdepth=hab$wdepth,include=hab$include,onrefuge=hab$onrefuge)
  #Filter


  newhab<-newhab %>% filter(include == 1 & onrefuge == 1)


  # Define dimensions
  sites <- unique(newhab$site)
  reps <- unique(newhab$rep)
  years <- unique(newhab$year)
  ponds <- sort(unique(newhab$pond_name))

  num_sites <- length(sites)
  num_reps <- length(reps)
  num_years <- length(years)
  num_ponds <- length(ponds)

  # Function to create a 4D array for a given variable
  create_4D_array <- function(var_name) {
    array_data <- array(NA, dim = c(num_sites, num_reps, num_years, num_ponds),
                        dimnames = list(site = sites, rep = reps, year = years, pond = ponds))

    for (i in 1:nrow(newhab)) {
      site_index <- which(sites == newhab$site[i])
      rep_index <- which(reps == newhab$rep[i])
      year_index <- which(years == newhab$year[i])
      pond_index <- which(ponds == newhab$pond_name[i])
      array_data[site_index, rep_index, year_index, pond_index] <- newhab[[var_name]][i]
    }

    return(array_data)
  }

  # Create 4D arrays for each habitat variable
  pH_array <- create_4D_array("pH")
  wtemp_array <- create_4D_array("wtemp")
  doxygen_array <- create_4D_array("doxygen")
  wcond_array <- create_4D_array("wcond")
  veg_array <- create_4D_array("veg")
  wdepth_array <- create_4D_array("wdepth")

  replace_dimnames_with_numbers <- function(array_4d) {
    dimnames(array_4d) <- list(
      site = as.character(seq_len(dim(array_4d)[1])),  # Numeric indices for sites
      rep = as.character(seq_len(dim(array_4d)[2])),
      year = as.character(seq_len(dim(array_4d)[3])),  # Numeric indices for years
      pond = as.character(seq_len(dim(array_4d)[4]))   # Numeric indices for ponds
    )
    return(array_4d)
  }

  pH_array <- replace_dimnames_with_numbers(pH_array)
  wtemp_array <- replace_dimnames_with_numbers(wtemp_array)
  doxygen_array <- replace_dimnames_with_numbers(doxygen_array)
  wcond_array <- replace_dimnames_with_numbers(wcond_array)
  veg_array <- replace_dimnames_with_numbers(veg_array)
  wdepth_array <- replace_dimnames_with_numbers(wdepth_array)

  remove_dimnames <- function(array_4d) {
    dimnames(array_4d) <- NULL  # Remove all dimension names
    return(array_4d)
  }

  pH_4D <- remove_dimnames(pH_array)
  wtemp_4D <- remove_dimnames(wtemp_array)
  doxygen_4D <- remove_dimnames(doxygen_array)
  wcond_4D <- remove_dimnames(wcond_array)
  veg_4D <- remove_dimnames(veg_array)
  wdepth_4D <- remove_dimnames(wdepth_array)

  center_4D_array <- function(array_4d) {
    mean_val <- mean(array_4d, na.rm = TRUE)  # Compute mean, ignoring NAs
    sd_val <- sd(array_4d, na.rm = TRUE)
    centered_array <- (array_4d - mean_val) / sd_val     # Subtract mean from all elements
    return(centered_array)
  }

  # Apply centering to all 4D arrays
  pH_centered <- center_4D_array(pH_4D)
  wtemp_centered <- center_4D_array(wtemp_4D)
  doxygen_centered <- center_4D_array(doxygen_4D)
  wcond_centered <- center_4D_array(wcond_4D)
  veg_centered <- center_4D_array(veg_4D)
  wdepth_centered <- center_4D_array(wdepth_4D)





  # Function to average over the 'rep' dimension while handling NAs
  average_over_rep <- function(array_4d) {
    apply(array_4d, c(1, 3, 4), function(x) {
      if (all(is.na(x))) {
        return(NA)  # If all values are NA, return NA
      } else {
        return(mean(x, na.rm = TRUE))  # Otherwise, compute the mean ignoring NAs
      }
    })
  }

  # Apply the function to create 3D arrays by averaging over rep
  pH_3D <- average_over_rep(pH_4D)
  wtemp_3D <- average_over_rep(wtemp_4D)
  doxygen_3D <- average_over_rep(doxygen_4D)
  wcond_3D <- average_over_rep(wcond_4D)
  veg_3D <- average_over_rep(veg_4D)
  wdepth_3D <- average_over_rep(wdepth_4D)

  # Function to replace dimnames with numeric indices
  replace_dimnames_with_numbers <- function(array_3d) {
    dimnames(array_3d) <- list(
      site = as.character(seq_len(dim(array_3d)[1])),  # Numeric indices for sites
      year = as.character(seq_len(dim(array_3d)[2])),  # Numeric indices for years
      pond = as.character(seq_len(dim(array_3d)[3]))   # Numeric indices for ponds
    )
    return(array_3d)
  }

  # Apply function to all 3D arrays
  pH_3D <- replace_dimnames_with_numbers(pH_3D)
  wtemp_3D <- replace_dimnames_with_numbers(wtemp_3D)
  doxygen_3D <- replace_dimnames_with_numbers(doxygen_3D)
  wcond_3D <- replace_dimnames_with_numbers(wcond_3D)
  veg_3D <- replace_dimnames_with_numbers(veg_3D)
  wdepth_3D <- replace_dimnames_with_numbers(wdepth_3D)

  # Function to remove dimnames entirely
  remove_dimnames <- function(array_3d) {
    dimnames(array_3d) <- NULL  # Remove all dimension names
    return(array_3d)
  }

  # Apply function to all 3D arrays
  pH_3D <- remove_dimnames(pH_3D)
  wtemp_3D <- remove_dimnames(wtemp_3D)
  doxygen_3D <- remove_dimnames(doxygen_3D)
  wcond_3D <- remove_dimnames(wcond_3D)
  veg_3D <- remove_dimnames(veg_3D)
  wdepth_3D <- remove_dimnames(wdepth_3D)

  pH_3D_centered <- center_4D_array(pH_3D)
  wtemp_3D_centered <- center_4D_array(wtemp_3D)
  doxygen_3D_centered <- center_4D_array(doxygen_3D)
  wcond_3D_centered <- center_4D_array(wcond_3D)
  veg_3D_centered <- center_4D_array(veg_3D)
  wdepth_3D_centered <- center_4D_array(wdepth_3D)



  # Format count data
  count$yr<-as.numeric(as.factor(count$year))


  #Reorganize count data by Site, Wetland Pond, Year, and Species
  newdat<-data.frame(year=count$yr,day=count$day,
                     site=count$site,y=count$YCHUB,include=count$include,
                     onrefuge=count$onrefuge,pond.name=count$pond_name)
  countdata<-newdat[order(newdat$day,newdat$pond.name,newdat$year,newdat$site),]

  #Filter
  countdat1<-countdata %>% filter(include == 1 & onrefuge == 1)
  countdat1$pname1<-as.numeric(as.factor(countdat1$pond.name))

  sites <- unique(countdat1$site)
  reps <- unique(countdat1$day)
  years <- unique(countdat1$year)
  ponds <- unique(countdat1$pname)

  num_sites <- length(sites)
  num_reps <- length(reps)
  num_years <- length(years)
  num_ponds <- length(ponds)

  y_array <- array(NA, dim = c(num_sites, num_reps, num_years, num_ponds))

  for (i in 1:nrow(countdat1)) {
    site_index <- which(sites == countdat1$site[i])
    rep_index <- which(reps == countdat1$day[i])
    year_index <- which(years == countdat1$year[i])
    pond_index <- which(ponds == countdat1$pname1[i])

    y_array[site_index, rep_index, year_index, pond_index] <- countdat1$y[i]
  }



  print("Executing JAGS model to assess relationships with select habitat variables. This may take several minutes to hours.",quote=FALSE)

  #Define model file name and create JAGS model to assess relationships with
  #known habitat variables
  modelFilename="Habitat.Model.txt"
  cat("
      model{
      phi~dunif(0.1,100)

      beta0 ~ dnorm(0, 0.01)
      beta.veg ~ dnorm(0, 0.01)
      beta.wtemp ~ dnorm(0, 0.01)
      #beta.wcond ~ dnorm(0, 0.01)
      #beta.doxygen ~ dnorm(0, 0.01)
      #beta.wdepth ~ dnorm(0, 0.01)

      alpha0 ~ dnorm(0, 1)I(-2,2)


      #alpha.veg ~ dnorm(0, 0.01)
      #alpha.pH ~ dnorm(0, 0.01)
      #alpha.wcond ~ dnorm(0, 0.01)
      #alpha.wdepth ~ dnorm(0, 0.01)
      #alpha.doxygen ~ dnorm(0, 0.01)

      for(i in 1:nsites){
        eta[i] ~ dgamma(phi, phi)
        for(t in 1:nyears){
          for(k in 1:nponds){
            N[i,t,k] ~ dpois(lambda[i,t,k])
            lambda[i,t,k] <- mu[i,t,k] * eta[i]
            log(mu[i,t,k]) <- beta0 + beta.veg * veg_3D[i,t,k] + beta.wtemp * wtemp_3D[i,t,k]
          }
        }
      }

      sigma_veg3D ~ dunif(0, 5)
      tau_veg3D <- pow(sigma_veg3D, -2)
      sigma_wtemp3D ~ dunif(0, 5)
      tau_wtemp3D <- pow(sigma_wtemp3D, -2)
      #sigma_wcond3D ~ dunif(0, 5)
      #tau_wcond3D <- pow(sigma_wcond3D, -2)
      #sigma_doxygen3D ~ dunif(0, 5)
      #tau_doxygen3D <- pow(sigma_doxygen3D, -2)
      #sigma_wdepth3D ~ dunif(0, 5)
      #tau_wdepth3D <- pow(sigma_wdepth3D,-2)

      for(i in 1:nsites){
        for(t in 1:nyears){
          for(k in 1:nponds){
            mu_veg3D[i,t,k] ~ dnorm(0, 0.01)
            mu_wtemp3D[i,t,k] ~ dnorm(0, 0.01)
            #mu_wcond3D[i,t,k] ~ dnorm(0, 0.01)
            #mu_wdepth3D[i,t,k] ~ dnorm(0, 0.01)
            #mu_doxygen3D[i,t,k] ~ dnorm(0, 0.01)
          }
        }
      }


      for(i in 1:nsites){
        for(t in 1:nyears){
          for(k in 1:nponds){
            veg_3D[i,t,k] ~ dnorm(mu_veg3D[i,t,k], tau_veg3D)
            wtemp_3D[i,t,k] ~ dnorm(mu_wtemp3D[i,t,k], tau_wtemp3D)
            #wcond_3D[i,t,k] ~ dnorm(mu_wcond3D[i,t,k], tau_wcond3D)
            #wdepth_3D[i,t,k] ~ dnorm(mu_wdepth3D[i,t,k], tau_wdepth3D)
            #doxygen_3D[i,t,k] ~ dnorm(mu_doxygen3D[i,t,k], tau_doxygen3D)
          }
        }
      }



      for(i in 1:nsites){
        for(j in 1:nreps){
          for(t in 1:nyears){
            for(k in 1:nponds){
              y[i,j,t,k] ~ dbin(q[i,j,t,k],N[i,t,k])
              logit(q[i,j,t,k]) <- alpha0
            }
          }
        }
      }

      #sigma_veg ~ dunif(0, 5)
      #tau_veg <- pow(sigma_veg,-2)
      #sigma_wcond ~ dunif(0, 5)
      #tau_wcond <- pow(sigma_wcond,-2)
      #sigma_wdepth ~ dunif(0, 5)
      #tau_wdepth <- pow(sigma_wdepth,-2)
      #sigma_doxygen ~ dunif(0, 5)
      #tau_doxygen <- pow(sigma_doxygen, -2)
      #sigma_pH ~ dunif(0, 5)
      #tau_pH <- pow(sigma_pH, -2)

      #for(i in 1:nsites){
      #  for(j in 1:nreps){
      #    for(t in 1:nyears){
      #      for(k in 1:nponds){
      #        mu_veg[i,j,t,k] ~ dnorm(0, 0.01)
              #mu_wcond[i,j,t,k] ~ dnorm(0, 0.01)
              #mu_wdepth[i,j,t,k] ~ dnorm(0, 0.01)
              #mu_doxygen[i,j,t,k] ~ dnorm(0, 0.01)
              #mu_pH[i,j,t,k] ~ dnorm(0, 0.01)
      #      }
      #    }
      #  }
      #}
      #for(i in 1:nsites){
      #  for(j in 1:nreps){
      #    for(t in 1:nyears){
      #      for(k in 1:nponds){
      #        veg_4D[i,j,t,k] ~ dnorm(mu_veg[i,j,t,k],tau_veg)
              #wcond_4D[i,j,t,k] ~ dnorm(mu_wcond[i,j,t,k],tau_wcond)
              #wdepth_4D[i,j,t,k] ~ dnorm(mu_wdepth[i,j,t,k], tau_wdepth)
              #doxygen_4D[i,j,t,k] ~ dnorm(mu_doxygen[i,j,t,k], tau_doxygen)
              #pH_4D[i,j,t,k] ~ dnorm(mu_pH[i,j,t,k], tau_pH)
      #      }
      #    }
      #  }
      #}



}",fill=TRUE,file=modelFilename)

  #Initial values
  Nst<-apply(y_array,c(1,3,4),function(x)sum(x,na.rm=TRUE))+10
  inits=function()list(N=Nst,phi=50)

  #ponds<-dplyr::dense_rank(countdat1$pname)
  nponds<-max(ponds)
  nsites <- max(sites)
  nreps <- max(reps)
  #years<-dplyr::dense_rank(countdat1$year)
  nyears<-max(years)

  #Bundle data
  ndata=list(y=y_array,nsites=nsites,nreps=nreps,nponds=nponds,nyears=nyears,
             veg_3D = veg_3D_centered, wtemp_3D = wtemp_3D_centered)

  #Parameters monitored
  params=c("phi","beta0","beta.veg","beta.wtemp","alpha0")
  #MCMC settings
  nc=3; nt=1; nb=15000; ni=350000

  #Call JAGS
  out2<-jags(ndata,inits,parameters.to.save=params,model.file=modelFilename,
             n.chains=nc,n.burnin=nb,n.thin=nt,n.iter=ni,parallel=TRUE,
             n.cores=nc,DIC=TRUE)

  #Create Wetland pond labels
  pond.name<-rep(as.character(unique(unlist(sort(countdata$pond.name)))))

  #Summarize posteriors for detection (alpha parameters)
  #Det.cond<-round(unlist(out2$mean$alpha.cond),3)
  #Det.cond<-as.vector(Det.cond)
  #Dcond.lower<-round(unlist(out2$q2.5$alpha.cond),3)
  #Dcond.lower<-as.vector(Dcond.lower)
  #Dcond.upper<-round(unlist(out2$q97.5$alpha.cond),3)
  #Dcond.upper<-as.vector(Dcond.upper)

  #Det.temp<-round(unlist(out2$mean$alpha.temp),2)
  #Det.temp<-as.vector(Det.temp)
  #DTemp.lower<-round(unlist(out2$q2.5$alpha.temp),2)
  #DTemp.lower<-as.vector(DTemp.lower)
  #DTemp.upper<-round(unlist(out2$q97.5$alpha.temp),2)
  #DTemp.upper<-as.vector(DTemp.upper)

  #Det.veg<-round(unlist(out2$mean$alpha.veg),3)
  #Det.veg<-as.vector(Det.veg)
  #DVeg.lower<-round(unlist(out2$q2.5$alpha.veg),3)
  #DVeg.lower<-as.vector(DVeg.lower)
  #DVeg.upper<-round(unlist(out2$q97.5$alpha.veg),3)
  #DVeg.upper<-as.vector(DVeg.upper)

  #Det.mean<-c(Det.cond,Det.veg)
  #Det.lower<-c(Dcond.lower,DVeg.lower)
  #Det.upper<-c(Dcond.upper,DVeg.upper)
  #Variable<-c("Water Conductivity","Vegetation")
  #Parameter<-rep(c("Detection"),2)

  #res2.D<-data.frame(Parameter=Parameter,Variable=Variable,Lower=Det.lower,
  #                   Mean=Det.mean,Upper=Det.upper)

  #Summarize posteriors for abundance (beta parameters)
  #N.depth<-round(unlist(out2$mean$beta.depth),3)
  #N.depth<-as.vector(N.depth)
  #NDepth.lower<-round(unlist(out2$q2.5$beta.depth),3)
  #NDepth.lower<-as.vector(NDepth.lower)
  #NDepth.upper<-round(unlist(out2$q97.5$beta.depth),3)
  #NDepth.upper<-as.vector(NDepth.upper)

  N.wtemp<-round(unlist(out2$mean$beta.wtemp),3)
  N.wtemp<-as.vector(N.wtemp)
  Nwtemp.lower<-round(unlist(out2$q2.5$beta.wtemp),3)
  Nwtemp.lower<-as.vector(Nwtemp.lower)
  Nwtemp.upper<-round(unlist(out2$q97.5$beta.wtemp),3)
  Nwtemp.upper<-as.vector(Nwtemp.upper)

  #N.cond<-round(unlist(out2$mean$beta.cond),2)
  #N.cond<-as.vector(N.cond)
  #Ncond.lower<-round(unlist(out2$q2.5$beta.cond),2)
  #Ncond.lower<-as.vector(Ncond.lower)
  #Ncond.upper<-round(unlist(out2$q97.5$beta.cond),2)
  #Ncond.upper<-as.vector(Ncond.upper)

  N.veg<-round(unlist(out2$mean$beta.veg),3)
  N.veg<-as.vector(N.veg)
  Nveg.lower<-round(unlist(out2$q2.5$beta.veg),3)
  Nveg.lower<-as.vector(Nveg.lower)
  Nveg.upper<-round(unlist(out2$q97.5$beta.veg),3)
  Nveg.upper<-as.vector(Nveg.upper)

  #N.algal<-round(unlist(out2$mean$beta.algal),2)
  #N.algal<-as.vector(N.algal)
  #Nalgal.lower<-round(unlist(out2$q2.5$beta.algal),2)
  #Nalgal.lower<-as.vector(Nalgal.lower)
  #Nalgal.upper<-round(unlist(out2$q97.5$beta.algal),2)
  #Nalgal.upper<-as.vector(Nalgal.upper)

  N.mean<-c(N.wtemp,N.veg)
  N.lower<-c(Nwtemp.lower,Nveg.lower)
  N.upper<-c(Nwtemp.upper,Nveg.upper)
  Variable<-c("Water Temperature",
              "Vegetation")
  Parameter<-rep(c("Abundance"),2)

  res2.N<-data.frame(Parameter=Parameter,Variable=Variable,Lower=N.lower,
                     Mean=N.mean,Upper=N.upper)

  #phi.mu<-round(unlist(out2$mean$phi),3)
  #phi.mu<-as.vector(phi.mu)
  #phi.lower<-round(unlist(out2$q2.5$phi),3)
  #phi.lower<-as.vector(phi.lower)
  #phi.upper<-round(unlist(out2$q97.5$phi),3)
  #phi.upper<-as.vector(phi.upper)

  #sdpond.mu<-round(unlist(out2$mean$sd_pond),3)
  #sdpond.mu<-as.vector(sdpond.mu)
  #sdpond.lower<-round(unlist(out2$q2.5$sd_pond),3)
  #sdpond.lower<-as.vector(sdpond.lower)
  #sdpond.upper<-round(unlist(out2$q97.5$sd_pond),3)
  #sdpond.upper<-as.vector(sdpond.upper)

  #N.lower=c(phi.lower,sdpond.lower)
  #N.upper=c(phi.upper,sdpond.upper)
  #N.mean<-c(phi.mu,sdpond.mu)

  #Variable<-c("phi (Overdispersion and site-level effect)","Random pond stdev")
  #Parameter<-rep(c("Random Effect"),2)
  #res2.R<-data.frame(Parameter=Parameter,Variable=Variable,Lower=N.lower,
  #Mean=N.mean,Upper=N.upper)

  res2<-rbind(res2.N)


  #Capture and Write results to working directory (R Data Files)
  write.csv(res2,"YaquiChubSHabitatModelParameters.csv",row.names=F)

  plot<-ggplot(res2,aes(Mean,Variable,colour=factor(Variable)))+
    geom_point(size=4)+
    xlim(-1,1) +
    geom_vline(aes(xintercept=0.0),color="black",linewidth=1)+
    geom_errorbarh(aes(xmin=Lower,xmax=Upper),height=.2,linewidth=1)+
    facet_wrap(~Parameter,ncol=2,scales="free_x")+
    guides(colour="none")+
    theme_bw()+
    xlab("Parameter Estimate")+
    ylab("")+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=16),
          strip.text.x=element_text(size=12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  print(plot)
  ggsave(paste0(source, "/YaquiChubWetlandPondModelParameterFigure.tiff"),plot=plot,
         width=12,height=7,dpi=300)

  print("Results for the Habitat Model Model are saved and stored in your working directory.",quote=FALSE)

}else if(species=="BSHINER"){


  #Reorganize habitat data by Site, Wetland Pond, and Year
  newhab<-data.frame(pond_name=hab$pond_name,rep=hab$day,year=hab$year,site=hab$site,
                     pH=hab$pH,wtemp=hab$wtemp,
                     doxygen=hab$doxygen,wcond=hab$wcond,
                     veg=hab$veg,
                     wdepth=hab$wdepth,include=hab$include,onrefuge=hab$onrefuge)
  #Filter
  library(tidyr)
  library(dplyr)

  newhab<-newhab %>% filter(include == 1 & onrefuge == 1)


  # Define dimensions
  sites <- unique(newhab$site)
  reps <- unique(newhab$rep)
  years <- unique(newhab$year)
  ponds <- sort(unique(newhab$pond_name))

  num_sites <- length(sites)
  num_reps <- length(reps)
  num_years <- length(years)
  num_ponds <- length(ponds)

  # Function to create a 4D array for a given variable
  create_4D_array <- function(var_name) {
    array_data <- array(NA, dim = c(num_sites, num_reps, num_ponds, num_years),
                        dimnames = list(site = sites, rep = reps, pond = ponds, year = years))

    for (i in 1:nrow(newhab)) {
      site_index <- which(sites == newhab$site[i])
      rep_index <- which(reps == newhab$rep[i])
      pond_index <- which(ponds == newhab$pond_name[i])
      year_index <- which(years == newhab$year[i])
      array_data[site_index, rep_index, pond_index, year_index] <- newhab[[var_name]][i]
    }

    return(array_data)
  }

  # Create 4D arrays for each habitat variable
  pH_array <- create_4D_array("pH")
  wtemp_array <- create_4D_array("wtemp")
  doxygen_array <- create_4D_array("doxygen")
  wcond_array <- create_4D_array("wcond")
  veg_array <- create_4D_array("veg")
  wdepth_array <- create_4D_array("wdepth")

  replace_dimnames_with_numbers <- function(array_4d) {
    dimnames(array_4d) <- list(
      site = as.character(seq_len(dim(array_4d)[1])),  # Numeric indices for sites
      rep = as.character(seq_len(dim(array_4d)[2])),
      pond = as.character(seq_len(dim(array_4d)[3])),   # Numeric indices for ponds
      year = as.character(seq_len(dim(array_4d)[4]))  # Numeric indices for years
    )
    return(array_4d)
  }

  pH_array <- replace_dimnames_with_numbers(pH_array)
  wtemp_array <- replace_dimnames_with_numbers(wtemp_array)
  doxygen_array <- replace_dimnames_with_numbers(doxygen_array)
  wcond_array <- replace_dimnames_with_numbers(wcond_array)
  veg_array <- replace_dimnames_with_numbers(veg_array)
  wdepth_array <- replace_dimnames_with_numbers(wdepth_array)

  remove_dimnames <- function(array_4d) {
    dimnames(array_4d) <- NULL  # Remove all dimension names
    return(array_4d)
  }

  pH_4D <- remove_dimnames(pH_array)
  wtemp_4D <- remove_dimnames(wtemp_array)
  doxygen_4D <- remove_dimnames(doxygen_array)
  wcond_4D <- remove_dimnames(wcond_array)
  veg_4D <- remove_dimnames(veg_array)
  wdepth_4D <- remove_dimnames(wdepth_array)

  center_4D_array <- function(array_4d) {
    mean_val <- mean(array_4d, na.rm = TRUE)  # Compute mean, ignoring NAs
    sd_val <- sd(array_4d, na.rm = TRUE)
    centered_array <- (array_4d - mean_val) / sd_val     # Subtract mean from all elements
    return(centered_array)
  }

  # Apply centering to all 4D arrays
  pH_centered <- center_4D_array(pH_4D)
  wtemp_centered <- center_4D_array(wtemp_4D)
  doxygen_centered <- center_4D_array(doxygen_4D)
  wcond_centered <- center_4D_array(wcond_4D)
  veg_centered <- center_4D_array(veg_4D)
  wdepth_centered <- center_4D_array(wdepth_4D)





  # Function to average over the 'rep' dimension while handling NAs
  average_over_rep <- function(array_4d) {
    apply(array_4d, c(1, 3, 4), function(x) {
      if (all(is.na(x))) {
        return(NA)  # If all values are NA, return NA
      } else {
        return(mean(x, na.rm = TRUE))  # Otherwise, compute the mean ignoring NAs
      }
    })
  }

  # Apply the function to create 3D arrays by averaging over rep
  pH_3D <- average_over_rep(pH_4D)
  wtemp_3D <- average_over_rep(wtemp_4D)
  doxygen_3D <- average_over_rep(doxygen_4D)
  wcond_3D <- average_over_rep(wcond_4D)
  veg_3D <- average_over_rep(veg_4D)
  wdepth_3D <- average_over_rep(wdepth_4D)

  # Function to replace dimnames with numeric indices
  replace_dimnames_with_numbers <- function(array_3d) {
    dimnames(array_3d) <- list(
      site = as.character(seq_len(dim(array_3d)[1])),  # Numeric indices for sites
      year = as.character(seq_len(dim(array_3d)[2])),  # Numeric indices for years
      pond = as.character(seq_len(dim(array_3d)[3]))   # Numeric indices for ponds
    )
    return(array_3d)
  }

  # Apply function to all 3D arrays
  pH_3D <- replace_dimnames_with_numbers(pH_3D)
  wtemp_3D <- replace_dimnames_with_numbers(wtemp_3D)
  doxygen_3D <- replace_dimnames_with_numbers(doxygen_3D)
  wcond_3D <- replace_dimnames_with_numbers(wcond_3D)
  veg_3D <- replace_dimnames_with_numbers(veg_3D)
  wdepth_3D <- replace_dimnames_with_numbers(wdepth_3D)

  # Function to remove dimnames entirely
  remove_dimnames <- function(array_3d) {
    dimnames(array_3d) <- NULL  # Remove all dimension names
    return(array_3d)
  }

  # Apply function to all 3D arrays
  pH_3D <- remove_dimnames(pH_3D)
  wtemp_3D <- remove_dimnames(wtemp_3D)
  doxygen_3D <- remove_dimnames(doxygen_3D)
  wcond_3D <- remove_dimnames(wcond_3D)
  veg_3D <- remove_dimnames(veg_3D)
  wdepth_3D <- remove_dimnames(wdepth_3D)

  pH_3D_centered <- center_4D_array(pH_3D)
  wtemp_3D_centered <- center_4D_array(wtemp_3D)
  doxygen_3D_centered <- center_4D_array(doxygen_3D)
  wcond_3D_centered <- center_4D_array(wcond_3D)
  veg_3D_centered <- center_4D_array(veg_3D)
  wdepth_3D_centered <- center_4D_array(wdepth_3D)



  # Format count data
  count$yr<-as.numeric(as.factor(count$year))


  #Reorganize count data by Site, Wetland Pond, Year, and Species
  newdat<-data.frame(year=count$yr,day=count$day,
                     site=count$site,y=count$BSHINER,include=count$include,
                     onrefuge=count$onrefuge,pond.name=count$pond_name)
  countdata<-newdat[order(newdat$day,newdat$pond.name,newdat$year,newdat$site),]

  #Filter
  countdat1<-countdata %>% filter(include == 1 & onrefuge == 1)
  countdat1$pname <- as.numeric(factor(countdat1$pond.name, levels = unique(countdat1$pond.name)))

  nsite=length(unique(countdat1$site))
  npond=length(unique(countdat1$pname))
  nday=length(unique(countdat1$day))
  nyear=length(unique(countdat1$year))

  site_index <- as.numeric(factor(countdat1$site, levels = unique(countdat1$site)))
  day_index <- as.numeric(factor(countdat1$day, levels = unique(countdat1$day)))
  pond_index <- as.numeric(factor(countdat1$pname, levels = sort(unique(countdat1$pname))))
  year_index <- as.numeric(factor(countdat1$year, levels = unique(countdat1$year)))

  cat("Max site index:", max(site_index), "Expected:", nsite, "\n")
  cat("Max day index:", max(day_index), "Expected:", nday, "\n")
  cat("Max pond index:", max(pond_index), "Expected:", npond, "\n")
  cat("Max year index:", max(year_index), "Expected:", nyear, "\n")

  y_filtered <- array(NA, dim = c(nsite, nday, npond, nyear))

  for (i in seq_len(nrow(countdat1))) {
    si <- site_index[i]
    di <- day_index[i]
    pi <- pond_index[i]  # Now correctly mapped!
    yi <- year_index[i]

    # Final safety check before assignment
    if (!is.na(si) & !is.na(di) & !is.na(pi) & !is.na(yi) &
        si > 0 & di > 0 & pi > 0 & yi > 0 &
        si <= nsite & di <= nday & pi <= npond & yi <= nyear) {
      y_filtered[si, di, pi, yi] <- countdat1$y[i]
    } else {
      warning(paste("Skipping row", i,
                    "si =", si, "di =", di, "pi =", pi, "yi =", yi,
                    "(out of bounds)"))
    }
  }


  modelFilename="Habitat.Model.txt"
  cat("
      model{
      phi~dunif(0.1,100)

      beta0 ~ dnorm(0, 0.01)
      beta.veg ~ dnorm(0, 0.01)I(-1.1,-0.80)
      beta.wtemp ~ dnorm(0, 0.01)I(-0.4,-0.05)
      #beta.doxygen ~ dnorm(0, 0.01)


      alpha0 ~ dnorm(0, 0.01)
      alpha.veg ~ dnorm(0, 0.01)I(-0.05,0.001)
      #alpha.wcond ~ dnorm(0, 0.01)
      alpha.wdepth ~ dnorm(0, 0.01)I(0.45,0.70)
      #alpha.doxygen ~ dnorm(0, 0.01)

      for(i in 1:nsites){
        eta[i] ~ dgamma(phi, phi)

          for(k in 1:nponds){
            for(t in 1:nyears){
              N[i,k,t] ~ dpois(lambda[i,k,t])
              lambda[i,k,t] <- mu[i,k,t] * eta[i]
              log(mu[i,k,t]) <- beta0 + beta.veg * veg_3D[i,k,t] + beta.wtemp * wtemp_3D[i,k,t]
          }
        }
      }

      sigma_veg3D ~ dunif(0, 5)
      tau_veg3D <- pow(sigma_veg3D, -2)
      sigma_wtemp3D ~ dunif(0, 5)
      tau_wtemp3D <- pow(sigma_wtemp3D, -2)
      #sigma_doxygen3D ~ dunif(0, 5)
      #tau_doxygen3D <- pow(sigma_doxygen3D, -2)

      for(i in 1:nsites){
          for(k in 1:nponds){
            for(t in 1:nyears){
            mu_veg3D[i,k,t] ~ dnorm(0, 0.01)
            mu_wtemp3D[i,k,t] ~ dnorm(0, 0.01)
            #mu_doxygen3D[i,k,t] ~ dnorm(0, 0.01)
          }
        }
      }


      for(i in 1:nsites){
          for(k in 1:nponds){
            for(t in 1:nyears){
            veg_3D[i,k,t] ~ dnorm(mu_veg3D[i,k,t], tau_veg3D)
            wtemp_3D[i,k,t] ~ dnorm(mu_wtemp3D[i,k,t], tau_wtemp3D)
            #doxygen_3D[i,k,t] ~ dnorm(mu_doxygen3D[i,k,t], tau_doxygen3D)
          }
        }
      }



      for(i in 1:nsites){
        for(j in 1:nreps){
          for(k in 1:nponds){
              for(t in 1:nyears){
              y[i,j,k,t] ~ dbin(q[i,j,k,t],N[i,k,t])
              logit(q[i,j,k,t]) <- alpha0 + alpha.veg * veg_4D[i,j,k,t] + alpha.wdepth * wdepth_4D[i,j,k,t]
            }
          }
        }
      }

      sigma_veg ~ dunif(0, 5)
      tau_veg <- pow(sigma_veg,-2)
      #sigma_wcond ~ dunif(0, 5)
      #tau_wcond <- pow(sigma_wcond,-2)
      sigma_wdepth ~ dunif(0, 5)
      tau_wdepth <- pow(sigma_wdepth,-2)
      #sigma_doxygen ~ dunif(0, 5)
      #tau_doxygen <- pow(sigma_doxygen, -2)

      for(i in 1:nsites){
        for(j in 1:nreps){
            for(k in 1:nponds){
              for(t in 1:nyears){
              mu_veg[i,j,k,t] ~ dnorm(0, 0.01)
              #mu_wcond[i,j,k,t] ~ dnorm(0, 0.01)
              mu_wdepth[i,j,k,t] ~ dnorm(0, 0.01)
              #mu_doxygen[i,j,k,t] ~ dnorm(0, 0.01)
            }
          }
        }
      }
      for(i in 1:nsites){
        for(j in 1:nreps){
            for(k in 1:nponds){
              for(t in 1:nyears){
              veg_4D[i,j,k,t] ~ dnorm(mu_veg[i,j,k,t],tau_veg)
              #wcond_4D[i,j,k,t] ~ dnorm(mu_wcond[i,j,k,t],tau_wcond)
              wdepth_4D[i,j,k,t] ~ dnorm(mu_wdepth[i,j,k,t], tau_wdepth)
              #doxygen_4D[i,j,k,t] ~ dnorm(mu_doxygen[i,j,k,t], tau_doxygen)
            }
          }
        }
      }



}",fill=TRUE,file=modelFilename)

  #Initial values
  Nst<-apply(y_filtered,c(1,3,4),function(x)sum(x,na.rm=TRUE))+1
  inits=function()list(N=Nst,phi=50)

  #Bundle data
  ndata=list(y=y_filtered,nsites=nsite,nreps=nday,nponds=npond,nyears=nyear,
             veg_4D = veg_centered, wdepth_4D = wdepth_centered,
             veg_3D = veg_3D_centered, wtemp_3D = wtemp_3D_centered)

  #Parameters monitored
  params=c("phi","beta0","beta.veg","beta.wtemp","alpha0","alpha.veg","alpha.wdepth")
  #MCMC settings
  nc=3; nt=1; nb=15000; ni=250000

  #Call JAGS
  out2<-jags(ndata,inits,parameters.to.save=params,model.file=modelFilename,
             n.chains=nc,n.burnin=nb,n.thin=nt,n.iter=ni,parallel=TRUE,
             n.cores=nc,DIC=TRUE)


  #Create Wetland pond labels
  pond.name<-rep(as.character(unique(unlist(sort(countdata$pond.name)))))

  #Summarize posteriors for detection (alpha parameters)

  Det.veg<-round(unlist(out2$mean$alpha.veg),3)
  Det.veg<-as.vector(Det.veg)
  DVeg.lower<-round(unlist(out2$q2.5$alpha.veg),3)
  DVeg.lower<-as.vector(DVeg.lower)
  DVeg.upper<-round(unlist(out2$q97.5$alpha.veg),3)
  DVeg.upper<-as.vector(DVeg.upper)

  Det.wdepth<-round(unlist(out2$mean$alpha.wdepth),3)
  Det.wdepth<-as.vector(Det.wdepth)
  Dwdepth.lower<-round(unlist(out2$q2.5$alpha.wdepth),3)
  Dwdepth.lower<-as.vector(Dwdepth.lower)
  Dwdepth.upper<-round(unlist(out2$q97.5$alpha.wdepth),3)
  Dwdepth.upper<-as.vector(Dwdepth.upper)

  Det.mean<-c(Det.veg,Det.wdepth)
  Det.lower<-c(DVeg.lower,Dwdepth.lower)
  Det.upper<-c(DVeg.upper,Dwdepth.upper)
  Variable<-c("Vegetation","Water Depth")
  Parameter<-rep(c("Detection"),2)

  res2.D<-data.frame(Parameter=Parameter,Variable=Variable,Lower=Det.lower,
                     Mean=Det.mean,Upper=Det.upper)

  #Summarize posteriors for abundance (beta parameters)
  N.wtemp<-round(unlist(out2$mean$beta.wtemp),3)
  N.wtemp<-as.vector(N.wtemp)
  Nwtemp.lower<-round(unlist(out2$q2.5$beta.wtemp),3)
  Nwtemp.lower<-as.vector(Nwtemp.lower)
  Nwtemp.upper<-round(unlist(out2$q97.5$beta.wtemp),3)
  Nwtemp.upper<-as.vector(Nwtemp.upper)

  N.veg<-round(unlist(out2$mean$beta.veg),3)
  N.veg<-as.vector(N.veg)
  Nveg.lower<-round(unlist(out2$q2.5$beta.veg),3)
  Nveg.lower<-as.vector(Nveg.lower)
  Nveg.upper<-round(unlist(out2$q97.5$beta.veg),3)
  Nveg.upper<-as.vector(Nveg.upper)

  N.mean<-c(N.wtemp,N.veg)
  N.lower<-c(Nwtemp.lower,Nveg.lower)
  N.upper<-c(Nwtemp.upper,Nveg.upper)
  Variable<-c("Water Temperature",
              "Vegetation")
  Parameter<-rep(c("Abundance"),2)

  res2.N<-data.frame(Parameter=Parameter,Variable=Variable,Lower=N.lower,
                     Mean=N.mean,Upper=N.upper)

  #phi.mu<-round(unlist(out2$mean$phi),3)
  #phi.mu<-as.vector(phi.mu)
  #phi.lower<-round(unlist(out2$q2.5$phi),3)
  #phi.lower<-as.vector(phi.lower)
  #phi.upper<-round(unlist(out2$q97.5$phi),3)
  #phi.upper<-as.vector(phi.upper)

  #sdpond.mu<-round(unlist(out2$mean$sd_pond),3)
  #sdpond.mu<-as.vector(sdpond.mu)
  #sdpond.lower<-round(unlist(out2$q2.5$sd_pond),3)
  #sdpond.lower<-as.vector(sdpond.lower)
  #sdpond.upper<-round(unlist(out2$q97.5$sd_pond),3)
  #sdpond.upper<-as.vector(sdpond.upper)

  #N.lower=c(phi.lower,sdpond.lower)
  #N.upper=c(phi.upper,sdpond.upper)
  #N.mean<-c(phi.mu,sdpond.mu)

  #Variable<-c("phi (Overdispersion and site-level effect)","Random pond stdev")
  #Parameter<-rep(c("Random Effect"),2)
  #res2.R<-data.frame(Parameter=Parameter,Variable=Variable,Lower=N.lower,
  #                   Mean=N.mean,Upper=N.upper)

  res2<-rbind(res2.N,res2.D)


  #Capture and Write results to working directory (R Data Files)
  write.csv(res2,paste0(source,"/BeautifulShinerHabitatModelParameters.csv"),row.names=F)

  plot<-ggplot(res2,aes(Mean,Variable,colour=factor(Variable)))+
    geom_point(size=4)+
    geom_vline(aes(xintercept=0.0),color="black",linewidth=1)+
    geom_errorbarh(aes(xmin=Lower,xmax=Upper),height=.2,linewidth=1)+
    facet_wrap(~Parameter,ncol=2,scales="free_x")+
    guides(colour="none")+
    theme_bw()+
    xlab("Parameter Estimate")+
    ylab("")+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=16),
          strip.text.x=element_text(size=12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  print(plot)
  ggsave(paste0(source,"/BeautifulShinerWetlandPondModelParameterFigure.tiff"),plot=plot,width=12,
         height=7,dpi=300)

  print("Results for the Habitat Model Model are saved and stored in your working directory.",quote=FALSE)
}

print("Executing JAGS model to assess relationships with select Management variables. This may take several minutes to hours.",quote=FALSE)





################################################################################
################################################################################
################################################################################


## Need to fix management model. The covariate dimensions need to be looked into to
## determine the best way to model them.

if(species=="YCHUB"){
  #Add pname field to count, mgmt, and hab .csv files
  #pname is a numeric wetland pond field
  #yr is a numeric year field
  newdat<-data.frame(year=count$yr,day=count$day,
                     site=count$site,y=count$YCHUB,include=count$include,
                     onrefuge=count$onrefuge,pond.name=count$pond_name)
  countdata<-newdat[order(newdat$day,newdat$pond.name,newdat$year,newdat$site),]

  #Filter
  #countdat1<-countdata %>% filter(include == 1)
  countdata$pname1<-as.numeric(as.factor(countdata$pond.name))

  sites <- unique(countdata$site)
  reps <- unique(countdata$day)
  years <- unique(countdata$year)
  ponds <- unique(countdata$pname)

  num_sites <- length(sites)
  num_reps <- length(reps)
  num_years <- length(years)
  num_ponds <- length(ponds)

  y_array <- array(NA, dim = c(num_sites, num_reps, num_ponds, num_years))

  for (i in 1:nrow(countdata)) {
    site_index <- which(sites == countdata$site[i])
    rep_index <- which(reps == countdata$day[i])
    year_index <- which(years == countdata$year[i])
    pond_index <- which(ponds == countdata$pname1[i])

    y_array[site_index, rep_index, pond_index, year_index] <- countdata$y[i]
  }

  #Reorganize management activity data by Site, Wetland Pond, and Year
  newmgmt<-data.frame(pond_name=mgmt$pond_name,year=mgmt$year,ychubrem=mgmt$YCHUB_removed,
                      bshinerrem=mgmt$BSHINER_removed,ytoprem=mgmt$YTOP_removed,
                      ychubstock=mgmt$YCHUB_stocked,
                      bshinerstock=mgmt$BSHINER_stocked,
                      ytopstock=mgmt$YTOP_stocked,include=mgmt$include,BShinerTimeSinceStocking = mgmt$BShinerTimeSinceStocking,
                      YChubTimeSinceStocking = mgmt$YChubTimeSinceStocking)

  #Filter
  mgmtdat1<-newmgmt %>% filter(include == 1)
  mgmtdat1$pname<-as.numeric(as.factor(mgmtdat1$pond_name))
  mgmtdat1$yr<-as.numeric(as.factor(mgmtdat1$year))

  #Reorder management activity data
  mgmtdata<-mgmtdat1[order(mgmtdat1$year,mgmtdat1$pname),]

  #Expand data.frame
  mgmtdat<-mgmtdata %>% group_by(pname,year,ychubrem,bshinerrem,ytoprem,
                                 ychubstock,bshinerstock,ytopstock, BShinerTimeSinceStocking, YChubTimeSinceStocking) %>% expand(sites = 1:10)

  mgmtdat<-mgmtdat[order(mgmtdat$pname),]

  ychubrem_array <- array(0, c(num_ponds, nyear))
  bshinerrem_array <- array(0, c(num_ponds, nyear))
  ytoprem_array <- array(0, c(num_ponds, nyear))
  ychubstock_array <- array(0, c(num_ponds, nyear))
  bshinerstock_array <- array(0, c(num_ponds, nyear))
  ytopstock_array <- array(0, c(num_ponds, nyear))

  pond_levels <- unique(as.character(mgmtdat$pname))
  year_levels <- unique(as.character(mgmtdat$year))

  for (i in seq_len(nrow(mgmtdat))) {
    pond_idx <- match(mgmtdat$pname[i], pond_levels)
    year_idx <- match(mgmtdat$year[i], year_levels)

    if (!is.na(pond_idx) && !is.na(year_idx) && pond_idx <= length(pond_levels) && year_idx <= length(year_levels)) {
      ychubrem_array[pond_idx, year_idx] <- ifelse(mgmtdat$ychubrem[i] >= 0, log(mgmtdat$ychubrem[i] + 1), 0)
      bshinerrem_array[pond_idx, year_idx] <- ifelse(mgmtdat$bshinerrem[i] >= 0, log(mgmtdat$bshinerrem[i] + 1), 0)
      ytoprem_array[pond_idx, year_idx] <- ifelse(mgmtdat$ytoprem[i] >= 0, log(mgmtdat$ytoprem[i] + 1), 0)
      ychubstock_array[pond_idx, year_idx] <- ifelse(mgmtdat$ychubstock[i] >= 0, log(mgmtdat$ychubstock[i] + 1), 0)
      bshinerstock_array[pond_idx, year_idx] <- ifelse(mgmtdat$bshinerstock[i] >= 0, log(mgmtdat$bshinerstock[i] + 1), 0)
      ytopstock_array[pond_idx, year_idx] <- ifelse(mgmtdat$ytopstock[i] >= 0, log(mgmtdat$ytopstock[i] + 1), 0)
    }
  }





  #Define model file name and create JAGS model to assess relationships with
  #known habitat variables
  modelFilename="Management.Model.txt"
  cat("
     model{
  phi ~ dunif(0,100)

  for(k in 1:npond){
    beta[k] ~ dnorm(0,0.01)
    eta[k] ~ dgamma(phi, phi)
  }


  beta0 ~ dnorm(0, 0.01)
  beta.ychubstock ~ dnorm(0, 0.1)


  for(i in 1:nsite){
    for(k in 1:npond){
      N[i, k, 1] ~ dpois(lambda[i, k])
      lambda[i, k] <- muL[i, k] * eta[k]
      log(muL[i, k]) <- beta[k]

      for(t in 2:nyear){
        log(mu[i, k, t - 1]) <- beta0 + beta.ychubstock * ychubstock[k,t]

        N[i, k, t] ~ dpois(mu[i, k, t - 1])

      }
    }
  }

  for(i in 1:nsite){
    for(j in 1:nrep){
      for(k in 1:npond){
        for(t in 1:nyear){
          y[i, j, k, t] ~ dbin(q[i, j, k, t], N[i, k, t])
          q[i, j, k, t] ~ dbeta(5,9) # Slightly informative prior. South Minckley 1605 in 2022
        }
      }
    }
  }




    }",fill=TRUE,file=modelFilename)


  #Initial values
  Nst<-apply(y_array,c(1,3,4),sum,na.rm=T)+1
  jags.inits=function()list(N=Nst)

  #Bundle data
  jags.data=list(y=y_array,nsite=num_sites,nrep=num_reps,npond=num_ponds,nyear=num_years,ychubstock = ychubstock_array)

  #Parameters monitored
  jags.params=c("phi","beta0", "beta.ychubstock")

  #MCMC settings
  nc=5; nt=1; nb=15000; ni=275000

  out<-jags(jags.data,jags.inits,parameters.to.save=jags.params,
            model.file=modelFilename,n.chains=nc,
            n.iter=ni,n.burnin=nb,n.thin=nt,parallel=TRUE,verbose=TRUE)

  #Summarize posteriors for abundance (beta parameters)
  N.ychubstock<-round(unlist(out$mean$beta.ychubstock),3)
  N.ychubstock<-as.vector(N.ychubstock)
  Nychubstock.lower<-round(unlist(out$q2.5$beta.ychubstock),3)
  Nychubstock.lower<-as.vector(Nychubstock.lower)
  Nychubstock.upper<-round(unlist(out$q97.5$beta.ychubstock),3)
  Nychubstock.upper<-as.vector(Nychubstock.upper)

  #N.ychubrem<-round(unlist(out$mean$beta.ychubrem),3)
  #N.ychubrem<-as.vector(N.ychubrem)
  #Nychubrem.lower<-round(unlist(out$q2.5$beta.ychubrem),3)
  #Nychubrem.lower<-as.vector(Nychubrem.lower)
  #Nychubrem.upper<-round(unlist(out$q97.5$beta.ychubrem),3)
  #Nychubrem.upper<-as.vector(Nychubrem.upper)

  #N.mean<-c(N.ytopstock,N.ytoprem,N.ychubstock,N.ychubrem,N.bshinestock,N.bshinerem)
  #N.lower<-c(Nytopstock.lower,Nytoprem.lower,Nychubstock.lower,Nychubrem.lower,Nbshinestock.lower,Nbshinerem.lower)
  #N.upper<-c(Nytopstock.upper,Nytoprem.upper,Nychubstock.upper,Nychubrem.upper,Nbshinestock.upper,Nbshinerem.upper)

  N.mean<-c(N.ychubstock)
  N.lower<-c(Nychubstock.lower)
  N.upper<-c(Nychubstock.upper)
  Variable<-c("Yaqui Chub Stocked")
  Parameter<-rep(c("Yaqui Chub Abundance"),1)

  res2.N<-data.frame(Parameter=Parameter,Variable=Variable,Lower=N.lower,
                     Mean=N.mean,Upper=N.upper)

  res2<-rbind(res2.N)


  #Capture and Write results to working directory (R Data Files)
  write.csv(res2,paste0(source, "/YaquiChubManagementStockingParameters.csv"),row.names=F)

  plot<-ggplot(res2,aes(Mean,Variable,colour=factor(Variable)))+
    geom_point(size=4)+
    geom_vline(aes(xintercept=0.0),color="black",linewidth=1)+
    geom_errorbarh(aes(xmin=Lower,xmax=Upper),height=.2,linewidth=1)+
    facet_wrap(~Parameter,ncol=2,scales="free_x")+
    guides(colour="none")+
    theme_bw()+
    xlab("Parameter Estimate")+
    ylab("")+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=16),
          strip.text.x=element_text(size=12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  print(plot)
  ggsave(paste0(source, "/YaquiChubWetlandManagementStockingParameterFigure.tiff"),plot=plot,
         width=12,height=7,dpi=300)

  print("Results of the Management Model (Removed) are saved and stored in your working directory.",quote=FALSE)

}else if(species=="BSHINER"){

  #Reorganize management activity data by Site, Wetland Pond, and Year
  #Reorganize management activity data by Site, Wetland Pond, and Year
  newmgmt<-data.frame(pond_name=mgmt$pond_name,
                      year=mgmt$year,
                      ychubrem=mgmt$YCHUB_removed,
                      bshinerrem=mgmt$BSHINER_removed,
                      ytoprem=mgmt$YTOP_removed,
                      ychubstock=mgmt$YCHUB_stocked,
                      bshinerstock=mgmt$BSHINER_stocked,
                      ytopstock=mgmt$YTOP_stocked,
                      onrefuge=mgmt$onrefuge,
                      include=mgmt$include,
                      BShinerTimeSinceStocking = mgmt$BShinerTimeSinceStocking,
                      YChubTimeSinceStocking = mgmt$YChubTimeSinceStocking)

  #Filter
  mgmtdat1<-newmgmt %>% filter(include == 1 & onrefuge == 1)
  mgmtdat1$pname<-as.numeric(as.factor(mgmtdat1$pond_name))
  mgmtdat1$yr<-as.numeric(as.factor(mgmtdat1$year))

  #Reorder management activity data
  mgmtdata<-mgmtdat1[order(mgmtdat1$year,mgmtdat1$pname),]

  #Expand data.frame
  mgmtdat<-mgmtdata %>% group_by(pname,year,ychubrem,bshinerrem,ytoprem,
                                 ychubstock,bshinerstock,ytopstock, BShinerTimeSinceStocking, YChubTimeSinceStocking) %>% expand(sites = 1:10)

  mgmtdat<-mgmtdat[order(mgmtdat$pname),]

  ychubrem_array <- array(0, c(num_ponds, nyear))
  bshinerrem_array <- array(0, c(num_ponds, nyear))
  ytoprem_array <- array(0, c(num_ponds, nyear))
  ychubstock_array <- array(0, c(num_ponds, nyear))
  bshinerstock_array <- array(0, c(num_ponds, nyear))
  ytopstock_array <- array(0, c(num_ponds, nyear))

  pond_levels <- unique(as.character(mgmtdat$pname))
  year_levels <- unique(as.character(mgmtdat$year))

  for (i in seq_len(nrow(mgmtdat))) {
    pond_idx <- match(mgmtdat$pname[i], pond_levels)
    year_idx <- match(mgmtdat$year[i], year_levels)

    if (!is.na(pond_idx) && !is.na(year_idx) && pond_idx <= length(pond_levels) && year_idx <= length(year_levels)) {
      ychubrem_array[pond_idx, year_idx] <- ifelse(mgmtdat$ychubrem[i] >= 0, log(mgmtdat$ychubrem[i] + 1), 0)
      bshinerrem_array[pond_idx, year_idx] <- ifelse(mgmtdat$bshinerrem[i] >= 0, log(mgmtdat$bshinerrem[i] + 1), 0)
      ytoprem_array[pond_idx, year_idx] <- ifelse(mgmtdat$ytoprem[i] >= 0, log(mgmtdat$ytoprem[i] + 1), 0)
      ychubstock_array[pond_idx, year_idx] <- ifelse(mgmtdat$ychubstock[i] >= 0, log(mgmtdat$ychubstock[i] + 1), 0)
      bshinerstock_array[pond_idx, year_idx] <- ifelse(mgmtdat$bshinerstock[i] >= 0, log(mgmtdat$bshinerstock[i] + 1), 0)
      ytopstock_array[pond_idx, year_idx] <- ifelse(mgmtdat$ytopstock[i] >= 0, log(mgmtdat$ytopstock[i] + 1), 0)
    }
  }





  #Define model file name and create JAGS model to assess relationships with
  #known habitat variables
  modelFilename="Management.Model.txt"
  cat("
     model{
  phi ~ dunif(0,100)

  for(k in 1:npond){
    beta[k] ~ dnorm(0,0.01)
    eta[k] ~ dgamma(phi, phi)
  }


  beta0 ~ dnorm(0, 0.01)
  beta.bshinerstock ~ dnorm(0, 0.1)
  beta.bshinerrem ~ dnorm(0, 0.1)


  for(i in 1:nsite){
    for(k in 1:npond){
      N[i, k, 1] ~ dpois(lambda[i, k])
      lambda[i, k] <- muL[i, k] * eta[k]
      log(muL[i, k]) <- beta[k]

      for(t in 2:nyear){
        log(mu[i, k, t - 1]) <- beta0 + beta.bshinerstock * bshinerstock[k,t] + beta.bshinerrem * bshinerrem[k,t]

        N[i, k, t] ~ dpois(mu[i, k, t - 1])

      }
    }
  }

  for(i in 1:nsite){
    for(j in 1:nrep){
      for(k in 1:npond){
        for(t in 1:nyear){
          y[i, j, k, t] ~ dbin(q[i, j, k, t], N[i, k, t])
          q[i, j, k, t] ~ dbeta(5,9) # Slightly informative prior. South Minckley 1605 in 2022
        }
      }
    }
  }




    }",fill=TRUE,file=modelFilename)


  #Initial values
  Nst<-apply(y_filtered,c(1,3,4),sum,na.rm=T)+1
  jags.inits=function()list(N=Nst)

  #Bundle data
  jags.data=list(y=y_filtered,nsite=num_sites,nrep=num_reps,npond=num_ponds,nyear=num_years,bshinerstock = bshinerstock_array, bshinerrem = bshinerrem_array)

  #Parameters monitored
  jags.params=c("phi","beta0", "beta.bshinerstock", "beta.bshinerrem")

  #MCMC settings
  nc=5; nt=1; nb=15000; ni=275000

  out2<-jags(jags.data,jags.inits,parameters.to.save=jags.params,
             model.file=modelFilename,n.chains=nc,
             n.iter=ni,n.burnin=nb,n.thin=nt,parallel=TRUE,verbose=TRUE)


  print(out2)

  #Create Wetland pond labels
  pond.name<-rep(as.character(unique(unlist(sort(countdat1$pond_name)))))

  #Summarize posteriors for abundance (beta parameters)
  N.bshinerstock<-round(unlist(out2$mean$beta.bshinerstock),2)
  N.bshinerstock<-as.vector(N.bshinerstock)
  Nbshinerstock.lower<-round(unlist(out2$q2.5$beta.bshinerstock),2)
  Nbshinerstock.lower<-as.vector(Nbshinerstock.lower)
  Nbshinerstock.upper<-round(unlist(out2$q97.5$beta.bshinerstock),2)
  Nbshinerstock.upper<-as.vector(Nbshinerstock.upper)

  N.bshinerrem<-round(unlist(out2$mean$beta.bshinerrem),2)
  N.bshinerrem<-as.vector(N.bshinerrem)
  Nbshinerrem.lower<-round(unlist(out2$q2.5$beta.bshinerrem),2)
  Nbshinerrem.lower<-as.vector(Nbshinerrem.lower)
  Nbshinerrem.upper<-round(unlist(out2$q97.5$beta.bshinerrem),2)
  Nbshinerrem.upper<-as.vector(Nbshinerrem.upper)

  #N.bshinerstock<-round(unlist(out2$mean$beta.bshinerstock),2)
  #N.bshinerstock<-as.vector(N.bshinerstock)
  #Nbshinerstock.lower<-round(unlist(out2$q2.5$beta.bshinerstock),2)
  #Nbshinerstock.lower<-as.vector(Nbshinerstock.lower)
  #Nbshinerstock.upper<-round(unlist(out2$q97.5$beta.bshinerstock),2)
  #Nbshinerstock.upper<-as.vector(Nbshinerstock.upper)

  #N.mean<-c(N.ytopstock,N.ytoprem,N.ychubstock,N.ychubrem,N.bshinestock,N.bshinerem)
  #N.lower<-c(Nytopstock.lower,Nytoprem.lower,Nychubstock.lower,Nychubrem.lower,Nbshinestock.lower,Nbshinerem.lower)
  #N.upper<-c(Nytopstock.upper,Nytoprem.upper,Nychubstock.upper,Nychubrem.upper,Nbshinestock.upper,Nbshinerem.upper)

  N.mean<-c(N.bshinerrem,N.bshinerstock)
  N.lower<-c(Nbshinerrem.lower,Nbshinerstock.lower)
  N.upper<-c(Nbshinerrem.upper,Nbshinerrem.upper)
  Variable<-c("Beautiful Shiner Removed","Beautiful Shiner Stocked")
  Parameter<-rep(c("Abundance","Abundance"))

  res2.N<-data.frame(Parameter=Parameter,Variable=Variable,Lower=N.lower,
                     Mean=N.mean,Upper=N.upper)

  res2<-rbind(res2.N)


  #Capture and Write results to working directory (R Data Files)
  write.csv(res2,paste0(source, "/BeautifulShinerManagementStockingParameters.csv"),row.names=F)

  plot<-ggplot(res2,aes(Mean,Variable,colour=factor(Variable)))+
    geom_point(size=4)+
    geom_vline(aes(xintercept=0.0),color="black",linewidth=1)+
    geom_errorbarh(aes(xmin=Lower,xmax=Upper),height=.2,linewidth=1)+
    facet_wrap(~Parameter,ncol=2,scales="free_x")+
    guides(colour="none")+
    theme_bw()+
    xlab("Parameter Estimate")+
    ylab("")+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=16),
          strip.text.x=element_text(size=12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  print(plot)
  ggsave(paste0(source, "/BeautifulShinerWetlandManagementStockingParameterFigure.tiff"),plot=plot,width=12,
         height=7,dpi=300)



}


print("If needed, please consult with the Regional Statistician (Dr. David R. Stewart) or the Regional Data Management Team once complete and if you have any concerns.")
}
