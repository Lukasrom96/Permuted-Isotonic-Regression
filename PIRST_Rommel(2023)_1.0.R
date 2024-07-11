###################################
####PIRST SIMULATION FUNCTIONS#####
###################################

# FINAL USED FUNCTION -> pirstsim()

#STEP 1
##Create Data for a single ST-Plot for one Person
pirstsim_sc <- function(npercon,      # Number of Points per Condition
                        ncons,        # Number of Conditions
                        overlap,      # Points overlapping
                        interactsize, # Single or multiple process model, vector
                        noise,        # Gaussian Noise
                        curve) {      # Linear, Concave, Sigmoid
  ##producing the linear funtion
  npoints <- npercon * ncons
  Condition <- c()
  for (i in 1:ncons) {
    Condition <- c(Condition, rep(i, npercon))
  }
  p <- overlap
  M <- matrix(nrow = ncons, ncol = npercon)
  # producing the points
  for (i in 1:ncons) {
    for (j in 1:npercon){
      M[i, j] <- 0.7 * (j / npercon) + 0.15
      M[i, j] <- 0.7 * (j / npercon) + 0.15
    }
  }
  X <- c()
  for (i in 1:nrow(M)){
    X <- c(X, M[i, ])
  }
  # Evenly strechtch points out
  if (ncons == 1) { 
    correct <- ncons+1
  } else {
    correct <- ncons
  }
  
  for (i in 1:npoints) {
    X[i] <- X[i] + (-0.35 + 0.7 * (Condition[i] - 1) / (correct - 1)) * ((p - 1 / ncons) / npercon)
  }
  ## Linear, curve, curve down, sigmoid or mirrored-sigmoid cases
  if (curve == 0) {
    # Linear
    Y <- X
    X <- X
  } else if (curve == 1) {
    # Transfroming to Concave Up
    Y <- sin((90 * X + 270) * pi / 180) + 1
    X <- cos((90 * X + 270) * pi / 180)
  } else if (curve == -1) {
    # Transforming to Concave Down
    Y <- cos((90 * X + 270) * pi / 180)
    X <- sin((90 * X + 270) * pi / 180) + 1
  } else if (curve == 2) {
    Y <- 0.5 * (1 + tanh(X * 10 - 5))
    X <- X
  }else if (curve == -2) {
    Y <- X
    X <- 0.5 * (1 + tanh(X * 10 - 5))
  } else {
    warning("invalid Curve input")
  }
  ## Adding the Interaction constant to Conditions
  for (i in 1:npoints) {
    Y[i] <- Y[i] + interactsize[Condition[i]]
  }
  ## Adding noise
  Y <- Y + rnorm(npoints, 0, noise)
  X <- X + rnorm(npoints, 0, noise)
  ## Combine to a dataframe
  Trace <- rep(seq(1, npercon, 1), ncons)
  simdata <- as.data.frame(cbind(X, Y, Condition, Trace))
  simdata
}  

#STEP 2
##Repeat pirstsim_sc for multiple measures on one Person
pirstsim_rep_sc <- function(npercon,
                            ncons,
                            overlap, 
                            interactsize, 
                            noise, 
                            curve,
                            nmeasures){
  measurements <- c()
  for (i in 1: nmeasures) { 
    m <- pirstsim_sc(npercon, ncons, overlap, interactsize, noise, curve)
    m$Measurement <- rep(i, nrow(m))
    measurements <- rbind(measurements, m)
  }
  measurements
}

#STEP 3
##Create a Full STA Data-Set of multiple persons with multiple measures
pirstsim <- function(npercon = 4,
                     ncons = 2,
                     overlap = 1,
                     interactsize = c(0, 0),
                     noise = 0,
                     curve = 0,
                     nmeasures = 10,
                     cases = 100) {
  
  if (length(interactsize) != ncons) {
    stop("Number of Conditions must match length(interactsize)")
  }
  exp_simdata <- c()
  for (i in 1:cases) {
    simdata_part <- pirstsim_rep_sc(npercon, ncons, overlap, interactsize, noise, curve, nmeasures)
    simdata_part$Person <- c(rep(i, npercon * ncons))
    exp_simdata <- rbind(exp_simdata, simdata_part)
  }
  exp_simdata
}

#ADDITIONAL TOOL
#set nmeasures and cases to 1 and use next function to inspect a single case
st_plot <- function(simdata, namex = "Variable A", namey = "Variable B") {
  names(simdata) <- c("X", "Y", "Condition")
  simdata$Condition <- as.factor(simdata$Condition)
  ggplot2::ggplot(simdata, ggplot2::aes(x = X, y = Y, xax)) +
    ggplot2::geom_point(ggplot2::aes(color = Condition), size = 2) +
    ggplot2::ggtitle("State-Trace-Plot") +
    ggplot2::xlab(namex) +
    ggplot2::ylab(namey) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}

################################################################################

###################################
####PIRST- Algorithm FUNCTIONS#####
###################################

# FINAL USED FUNCTION -> pirst()

#Aggregates X and Y over repeated measures of one Person
aggregate <- function(data) {
  newx <- c()
  newy <- c()
  newcondition <- c()
  newtrace <- c()
  newperson <- c()
  datanew <- c()
  for (i in 1:max(data$Person)) {
    for (j in 1:max(data$Condition)) {
      for (h in 1:max(data$Trace)) {
        newx <-  mean(data[data$Person == i &
                             data$Condition == j &
                             data$Trace == h, ]$X, na.rm = TRUE)
        newy <-  mean(data[data$Person == i &
                             data$Condition == j &
                             data$Trace == h, ]$Y, na.rm = TRUE)
        newcondition <- j
        newtrace <- h
        newperson <- i
        newrow <- cbind(newx, newy, newcondition, newtrace, newperson)
        datanew <- rbind.data.frame(datanew, newrow)
      }
    }
  }
  names(datanew) <- c("X", "Y", "Condition", "Trace", "Person")
  datanew
}

#Identifies eligible points
find_eligible_points <- function(data) {
  names(data) <- c("X", "Y", "Condition", "Trace")
  #Create a Dataframe of eligibilityvalues for the data
  eligibilityframe <- 
    as.data.frame(matrix(1, nrow = nrow(data), max(data$Condition)))
  eligibilitynames <- c()
  for (i in 1:(max(data$Condition))) {
    eligibilitynames[i] <- paste0("eligibility ", i)
  }
  names(eligibilityframe) <- eligibilitynames
  data <- cbind(data, eligibilityframe) # Add it to the data
  # in the next process the algorithm identifies all conditions  a point is
  # eligible with and adds that to the Dataframe
  data <-
    data[order(data$Condition, data$X), ] # needed for next function
  for (p in 1:max(data$Condition)) {
    for (q in 1:max(data$Condition)) {
      for (i in 1:(nrow(data))) {
        if ((data$Condition)[i] == p) {
          data[i, p + 4] <- 1
        } else if (data$Condition[i] == q) {
          if ((data$X[i] >= max(data$X[data$Condition == p]) &&
               data$Y[i] >= max(data$Y[data$Condition == p])) ||
              data$X[i] <= min(data$X[data$Condition == p]) &&
              data$Y[i] <= min(data$Y[data$Condition == p])) {
            data[i, p + 4] <- 0
          }
        }
      }
    }
  }
  data
}

#Identifies the Conditions a point is eligible with (for ncons > 2)
combfinder <- function(data) {
  eligibilityframe <- data[, 5:(4 + max(data$Condition))]
  eligibilitycombs <- dplyr::distinct(eligibilityframe)
  combination <- seq(1, nrow(eligibilitycombs), 1)
  eligibilitycombs$combination <- combination
  eligibility <- c(rep(0, nrow((data))))
  for (i in 1:nrow(eligibilitycombs)) {
    for (j in 1:(nrow(data))) {
      if (all(data[j, 5:(4 + max(data$Condition))] == eligibilitycombs[i, 1:max(data$Condition)]) == TRUE) {
        eligibility[j] <- eligibilitycombs$combination[i]
      }
    }
  }
  data$eligibility <- eligibility # Every Point gets a specific eligibility
  data                            # decoding-number
}

#Identifies if data is not completely overlapping(ncons >2) or not overlapping at all
find_cuts <- function(data) {
  cuts <- data
  result <- 0
  #starts at condition one and checks from conditions to condition if there is
  #any overlap to find out if there is a cut between conditions
  for (j in 1:max(cuts$Condition)) {
    for (i in 1:nrow(cuts)) {
      if ((cuts[i, 4 + 1] == 1) && (cuts$Condition[i] != 1)) {
        cuts[, 4 + j][cuts$Condition == cuts$Condition[i]] = 1
        cuts$Condition[i] = 1
      }
    }
  }
  if (all(cuts$Condition == data$Condition) == TRUE) {
    result <- 1 #no permutations possible
  } else if (all(cuts$Condition == rep(1, nrow(cuts))) == TRUE) {
    result <- 0 #no cuts, overlap of all Conditions
  } else {
    result <- 2 #at least one condition is isolated
  }
  result
}

#runs extremely fast isotonic regression
fastiso <- function (y,x) {
  if (anyDuplicated (x) == 0) {
    switchback <- seq_along(y) 
    isomat <- cbind(y,x,switchback)
    isomat <- isomat[Rfast::Order(isomat[, 2]),]
    isomat[, 1] <- monotone::monotone(isomat[, 1])
    isomat <- isomat[Rfast::Order(isomat[, 3]),]
    y <- isomat[, 1]
    y
  } else {#tie handling, because monotone does not have one (extremly few cases)
    isoreg <- isotone::gpava(x,y)
    y <- isoreg[["x"]]
    y
  }
}

#runs isotonic regression with x-axis a predictor and adds values to data
isoreg <- function(data) {
  Y_isoreg <- c()
  for (i in 1:max(data$Condition)) {
    isoreg <- fastiso(data$Y[data$Condition == i],
                      data$X[data$Condition == i])
    Y_isoreg <- c(Y_isoreg, isoreg)
  }
  data$Y_isoreg <- Y_isoreg
  data
}

#calculates SSE for the original Data
SSE_calc <- function(data) {
  SSE_gather <- c()
  for (i in 1:max(data$Condition)) {
    SSE_hold <- sum((data$Y_isoreg[data$Condition == i] -
                       data$Y[data$Condition == i]) ^ 2)
    SSE_gather[i] <- SSE_hold
  }
  SSE <- sum(SSE_gather)
  SSE
}

#Shuffles Conditionslabels from points of the same eligibilitycombination randomly
shuffle_more <- function(data) {
  new_perm <- data[c("Condition","eligibility")]
  for (i in 1:max(data$eligibility)) {
    if (length(new_perm$Condition[new_perm$eligibility == i]) > 1) {
      shuffle <- sample(new_perm$Condition[new_perm$eligibility == i])
      new_perm$Condition[data$eligibility == i] <- shuffle
    }
  }
  data$Condition <- new_perm$Condition
  data
}

#calculates SSE of a Permutation without adding any Values to a Dataframe
isoreg_ssepermcalc <- function(data) {
  SSE_gather <- c()
  for (i in 1:max(data$Condition)) {
    iso_reg <- fastiso(data$Y[data$Condition == i],data$X[data$Condition == i])
    SSE_gather[i] <- sum((iso_reg - data$Y[data$Condition == i]) ^ 2)
  }
  SSE <- sum(SSE_gather)
  SSE
}

#combines Shuffle_more and isoreg_ssepermcalc to permutate and get SSE in one step
permutate <- function(data, nperms) {
  cont_shuffle <-
    data$Condition #for shuffling again if permutation = original
  SSE_all <- c(rep(0, nperms))
  for (i in 1:nperms) {
    data <- shuffle_more(data)
    count <- 0
    while ((all(data$Condition == cont_shuffle) == TRUE) && (count < 1000) ) {
      data <- shuffle_more(data) #Bugfix for a very rare case of overlap but no
      count <- count + 1         # possible Permuations where the algorithm hangs up
      #...annoying
    }
    if (count < 1000){
      SSE_all[i] <- isoreg_ssepermcalc(data)
    } else {   #in this rare case, there are no possible permutations
      SSE_all <- c(rep(NaN, nperms))
      message("BUT no permutations are possible, this subject will not have impact on the Test-statistic")
      break
    }
  }
  SSE_all
} #nperms = number of Permutations

#HEART OF THE ALGORITHM
#combines all previous functions to a PIRST for a single case
pirst_sc <- function(data, nperms) {
  data <- find_eligible_points(data)
  data <- combfinder(data)
  data <- isoreg(data)
  SSE_compare <- SSE_calc(data)
  info <- find_cuts(data)
  if (info == 1) {
    message(
      "Because of no overlap, this subject will not have impact on the Test-statistic!"
    ) # No Overlap. Bad!
    SSE <- NaN
  } else if (info == 2) {
    message("Data is NOT completely overlapping!")
    SSE <- permutate(data, nperms) #just to inform Experimenter
  } else {
    message("Data is completely overlapping.")
    SSE <- permutate(data, nperms) #best case, everyting is fine
  }
  Proportion_greater <- sum(SSE > SSE_compare) / length(SSE)
  Proportion_equal <- sum(SSE == SSE_compare) / length(SSE)
  Proportion_smaller <- sum(SSE < SSE_compare) / length(SSE)
  Statistic <-
    list(Proportion_greater,
         Proportion_equal,
         Proportion_smaller,
         SSE,
         SSE_compare,
         data)
  names(Statistic) <-
    c(
      "Proportion greater",
      "Proportion equal",
      "Proportion smaller",
      "All SSE's",
      "SSE original",
      "Data"
    )
  Statistic
}

#applies pirst_sc to every person and calculates mean p-greater and other statistics
pirst <- function(data, nperms= 100) {
  names(data) <-
    c("X", "Y", "Condition", "Trace", "Measurement", "Person")
  SSE_greater <- c(rep(0, max(data$Person)))
  SSE_equal <- c(rep(0, max(data$Person)))
  SSE_smaller <- c(rep(0, max(data$Person)))
  SSE_original <- c(rep(0, max(data$Person)))
  output_data <- as.data.frame(c())
  data <- aggregate(data)
  for (i in 1:max(data$Person)) {
    calc_simdata <- data[data$Person == i, ]
    calc_simdata <- calc_simdata[, -5]
    message("For subject ", i)
    get_sc_data_greater <- pirst_sc(calc_simdata, nperms)
    SSE_greater[i] <- get_sc_data_greater[["Proportion greater"]]
    SSE_equal[i] <- get_sc_data_greater[["Proportion equal"]]
    SSE_smaller[i] <- get_sc_data_greater[["Proportion smaller"]]
    SSE_original[i] <- get_sc_data_greater[["SSE original"]]
    output_data <- rbind(output_data, get_sc_data_greater[["Data"]])
  }
  SSE_data <- as.data.frame(cbind(SSE_greater, SSE_equal, SSE_smaller, SSE_original))
  SSE_data <- na.omit(SSE_data)
  Prop_greater <- mean(SSE_data$SSE_greater)
  Prop_equal <- mean(SSE_data$SSE_equal)
  Prop_smaller <- mean(SSE_data$SSE_smaller)
  SSE_original_mean <- mean(SSE_data$SSE_original)
  Prop_greater_all <- SSE_data$SSE_greater
  SSE_original_all <- SSE_data$SSE_original
  Person <- data$Person
  output_data <- cbind(output_data, Person)
  
  Statistic <-
    list(
      Prop_greater,
      Prop_equal,
      Prop_smaller,
      SSE_original_mean,
      Prop_greater_all,
      SSE_original_all,
      output_data
    )
  names(Statistic) <-
    c(
      "Mean P greater",
      "Mean P equal",
      "Mean P smaller",
      "Mean of original SSE's",
      "Proportion of SSE's greater original SSE",
      "All original SSE's",
      "Data"
    )
  Statistic
}

################################################################################

##########################################
####PIRST- Significancetest FUNCTIONS#####
##########################################

# FINAL USED FUNCTION -> pirst_conf()

#calculates isotonic regressions for X and Y Values for a specific order
#Needed for easy_cmr
combined_isoreg_fit <- function(data) {
  isoY <- fastiso(data$Y, data$Order)
  isoX <- fastiso(data$X, data$Order)
  ErrorY <- sum((isoY - data$Y) ^ 2)
  ErrorX <- sum((isoX - data$X) ^ 2)
  fit <- ErrorY + ErrorX
  fit
}

#runs CMR for every condition seperately
#Needed for easy_cmr
preorder <- function(data) {
  for (i in 1:max(data$Condition)) {
    getpreordered_X <- fastiso(data$X[data$Condition == i],
                               data$Trace[data$Condition == i])
    getpreordered_Y <- fastiso(data$Y[data$Condition == i],
                               data$Trace[data$Condition == i])
    data$X[data$Condition == i] <- getpreordered_X
    data$Y[data$Condition == i] <- getpreordered_Y
  }
  data
}

#runs an easy cmr_Version with a sorting-algorithm
easy_cmr <- function(data) {
  data <- preorder(data)
  data$Order <- rep(0, nrow(data))
  cmr_dat <- data[data$Condition == 1, ]
  cmr_dat <- cmr_dat[order(cmr_dat$Condition, cmr_dat$Trace), ]
  cmr_dat$Order <- seq(1, nrow(cmr_dat))
  # Searches for the best fitting spot of a Point of an ordered Condition
  # in another ordered Condition and continues until the whole condition is
  # sorted in
  for (i in 2:max(data$Condition)) {
    shorten_search <- 0.5
    for (j in 1:max(data$Trace)) {
      best_order <- matrix(0, nrow = nrow(cmr_dat) + 1, ncol = 2)
      for (k in seq(from = shorten_search,
                    to = max(cmr_dat$Order) + 1,
                    by = 1)) {
        newrow <- data[data$Condition == i & data$Trace == j, ]
        newrow$Order[1] <- k
        cmr_dat_plus <- cmr_dat
        cmr_dat_plus <- dplyr::bind_rows(cmr_dat_plus, newrow)
        fit <- combined_isoreg_fit(cmr_dat_plus)
        best_order[k + 0.5, 1] <- fit
        best_order[k + 0.5, 2] <- k
      }
      if (shorten_search > 0.5) {
        best_order <- as.matrix(best_order[-c(1:(shorten_search - 0.5)), ])
      }
      #sometimes the best_order matrix get automatically transformed into a vector
      #this fixes indices and continues the algorithm
      if (ncol(best_order) > 1) {
        placement <- best_order[which.min(best_order[, 1]), 2]
      } else {
        placement <- best_order[2, 1]
      }
      newrow$Order <- placement
      cmr_dat <- rbind(cmr_dat, newrow)
      cmr_dat <- cmr_dat[order(cmr_dat$Order), ]
      cmr_dat$Order <- seq(1, nrow(cmr_dat))
      shorten_search <-
        cmr_dat$Order[cmr_dat$Condition == i & cmr_dat$Trace == j] + 0.5
    }
  }
  isoY <- fastiso(cmr_dat$Y, cmr_dat$Order)
  isoX <- fastiso(cmr_dat$X, cmr_dat$Order)
  cmr_dat$X <- isoX
  cmr_dat$Y <- isoY
  cmr_dat
}

#get the deviation-distribution of every point for every person
#uses the raw data and the pirst(data,...)-output
get_dist <- function(data, test) {
  datanew <- c()
  for (i in 1:max(data$Person)) {
    for (j in 1:max(data$Condition)) {
      for (h in 1:max(data$Trace)) {
        distx <- data[data$Person == i &
                        data$Condition == j &
                        data$Trace == h, ]$X -
          test[["Data"]][test[["Data"]]$Person == i &
                           test[["Data"]]$Condition == j &
                           test[["Data"]]$Trace == h, ]$X
        
        disty <-  data[data$Person == i &
                         data$Condition == j &
                         data$Trace == h, ]$Y -
          test[["Data"]][test[["Data"]]$Person == i &
                           test[["Data"]]$Condition == j &
                           test[["Data"]]$Trace == h, ]$Y
        
        newcondition <- rep(j, length(distx))
        newtrace <- rep(h, length(distx))
        newperson <- rep(i, length(distx))
        newrow <-
          cbind(distx, disty, newcondition, newtrace, newperson)
        datanew <- rbind.data.frame(datanew, newrow)
      }
    }
  }
  names(datanew) <-
    c("distX", "distY", "Condition", "Trace", "Person")
  datanew
}

# creates the nullmodel and the nulleffect-Dataframe 
create_null <- function(test) {
  #first create the nullmodel
  testdata <- test[["Data"]]
  nulldata <- c()
  for (i in 1:max(testdata$Person)) {
    nulldata_part <- easy_cmr(testdata[testdata$Person == i,])
    nulldata <- rbind(nulldata, nulldata_part)
  }
  nulldata <- cbind.data.frame(
    nulldata$X,
    nulldata$Y,
    nulldata$Condition,
    nulldata$Trace,
    rep(1, nrow(nulldata)),
    nulldata$Person
  )
  names(nulldata) <-
    c("X", "Y", "Condition", "Trace", "Measurement", "Person")
  nulldata
  #nullmodel is complete, now deviations can be added to points to generate
}
create_null_plus_dev <- function(nulldata, dist) {
  for (i in 1:max(dist$Person)) {
    for (j in 1:max(dist$Condition)) {
      for (h in 1:max(dist$Trace)) {
        randx <- dist[dist$Person == i &
                        dist$Condition == j &
                        dist$Trace == h,]$distX
        randy <- dist[dist$Person == i &
                        dist$Condition == j &
                        dist$Trace == h,]$distY
        randxmean <- 0
        randymean <- 0
  #####    UPDATE  
        for (k in 1:length(randx)) {
          randdev <- sample(1:length(randx), 1)
          randxmean <- randxmean + randx[randdev]
        }
        for (k in 1:length(randy)) {
          randdev <- sample(1:length(randy), 1)
          randymean <- randymean + randy[randdev]
        }
        randxmean <- randxmean / length(randx)
        randymean <- randymean / length(randy)
  #####    UPDATE
        nulldata[nulldata$Person == i &
                   nulldata$Condition == j &
                   nulldata$Trace == h,]$X <-
          nulldata[nulldata$Person == i &
                     nulldata$Condition == j &
                     nulldata$Trace == h, ]$X +
          randxmean
        nulldata[nulldata$Person == i &
                   nulldata$Condition == j &
                   nulldata$Trace == h,]$Y <-
          nulldata[nulldata$Person == i &
                     nulldata$Condition == j &
                     nulldata$Trace == h,]$Y +
          randymean
      }
    }
  }
  nulldata
}

# create null-distribution by drawing nulleffect-dataframes and collect mean p-greater
create_null_dist <- function(test, dist, ndraws, nperms = 100) {
  null_dist <- c()
  null_without_dev <- create_null(test)
  for (i in 1:ndraws) {
    nulldata <- create_null_plus_dev(null_without_dev, dist)
    null_test <- pirst(nulldata, nperms)
    null_dist[i] <- null_test[["Mean P greater"]]
  }
  null_dist
}


#Combines previous functions
#runs a significance Test for the found effect given the null-distribution
pirst_conf <- function(data, test, ndraws, nperms = 100, conf = 0.95) {
  dist <- get_dist(data, test)
  nulldist <- create_null_dist(test, dist, ndraws, nperms)
  upper <- round(conf * length(nulldist))
  nulldist <- sort(nulldist)
  upperbound <- nulldist[upper]
  p <- sum(nulldist > test[["Mean P greater"]]) / length(nulldist)
  significance <- FALSE
  if (upperbound < test[["Mean P greater"]]) {
    significance <- TRUE
  }
  Significancetest <-
    list(test[["Mean P greater"]], upperbound, significance, p, conf, nulldist)
  names(Significancetest) <-
    c("Mean P greater",
      "Upper",
      "Significance",
      "p-Value",
      "Confidence",
      "Null-Distribution")
  Significancetest
}

#All together- the final Test-function
pirst_test <- function(data, nperms = 100, nboot = 100, conf = 0.95){
  test <- pirst(data, nperms)
  Significancetest <- pirst_conf(data = data,
                                 test = test,
                                 ndraws = nboot,
                                 nperms = nperms,
                                 conf = conf)
  Significancetest
}

