
#' generate a meta-adjacency matrix with dispersal links
#' 
#' based on degree distribution from a metaweb, landscape spatial autocorrelation,
#' and dispersal kernel. In this version, these are input via their corresponding
#' categories, with values taken from the scripts that I first used to get them.
#'
#' @param param list of parameters, see ...
#' @param suit.df suitability of each sp given as a mean and sd
#' @param landscape.cat category of landscape spatial autocorrelation "saXX"
#' @param network.cat category of network degree distribution "ddXX"
#' @param dispersal.cat category of dispersal kernel distribution "dkXX"
#'
#' @return meta-adjacency matrix 
#' @export
#'
#' @examples
generate_landscape <- function(param,suit.df,
                               landscape.cat,
                               network.cat,
                               dispersal.cat){
  sp.names <- unique(suit.df$sp)
  num.sp <- length(sp.names)
  
  # initialize species info, to be filled below
  sp.traits <- data.frame(sp = sp.names, degree = 0, dispersal.distance = 0)
  
  richness <-  param$richness
  
  # to fully replicate the workflow of building networks, generate the set of
  # categories and take the exact value of the one chosen above.
  num.network.categories <- param$num.network.categories
  num.dispersal.categories <- param$num.dispersal.categories
  num.landscape.categories <- param$num.landscape.categories
  
  # landscape dimensions
  ncol <- param$ncol
  nrow <- param$nrow
  cells <- nrow*ncol
  
  # poisson mean
  min.lambda <- param$min.lambda
  max.lambda <- param$max.lambda # this should vary with richness. for S = 50, 15 gives connectance = 0.3
  
  # exponetial rate for dispersal
  min.rate <- param$min.rate
  max.rate <- param$max.rate
  
  # some constants for sampling interaction strengths
  int.mean <- param$int.mean
  int.sd <- param$int.sd
  tau <- param$tau
  diag.dom <- param$diag.dom
  
  # -------------------------------------------------------------------------
  # in this section, obtain the actual values of degree/dispersal/etc for the
  # network, based on the category specified above
  
  # network degree distribution
  degree.dist.gradient <- seq(from = min.lambda,
                              to = max.lambda, 
                              length.out = num.network.categories)
  # generative.models <- paste("dd",1:num.network.categories,sep="")
  generative.models <- paste("dd",sprintf("%02d", 1:num.network.categories),sep="")
  
  landscape.lambda <- degree.dist.gradient[which(generative.models == network.cat)]
  
  # landscape configuration
  landscape.gradient <- round(seq(0.1,0.9, length.out = num.landscape.categories),2)
  landscape.models <- paste("sa",sprintf("%02d", 1:num.landscape.categories),sep="")
  
  landscape.sa <- landscape.gradient[which(landscape.models == landscape.cat)]
  
  # dispersal kernel
  dispersal.rate.gradient <- seq(from = min.rate,
                                 to = max.rate, 
                                 length.out = num.dispersal.categories)
  
  dispersal.models <- paste("dk",sprintf("%02d", 1:num.dispersal.categories),sep="")
  
  landscape.dispersal <- dispersal.rate.gradient[which(dispersal.models == dispersal.cat)]
  
  # -------------------------------------------------------------------------
  # first, create landscape matrix
  
  my.landscape <- NLMR::nlm_fbm(ncol = ncol,nrow = nrow,fract_dim = landscape.sa)
  # show_landscape(my.landscape)
  my.landscape <- raster::as.matrix(my.landscape)
  
  # and the associated coordinates/distances
  land.df <- expand.grid(x = 1:ncol, y = 1:nrow)
  distance.matrix <- as.matrix(dist(land.df))
  
  land.df$cell <- 1:nrow(land.df)
  distance.df <- reshape2::melt(distance.matrix)
  names(distance.df) <- c("cell_from","cell_to","distance")
  
  # -------------------------------------------------------------------------
  # second, generate metaweb
  
  my.dist <- extraDistr::rtpois(n = richness,lambda = landscape.lambda,a = 0)
  # make sum even
  if (sum(my.dist) %% 2 != 0) { my.dist[1] <- my.dist[1] + 1 }
  my.net <- igraph::sample_degseq(my.dist,method = "vl")
  my.metaweb <- as.matrix(igraph::as_adjacency_matrix(my.net,type = "both"))
  # assign interaction strengths according to an "extended" normal dist
  weights <- abs(gamlss.dist::rSHASHo(sum(my.dist), mu = int.mean, 
                                      sigma = int.sd, nu = 0, tau = tau))
  my.metaweb[my.metaweb == 1] <- weights
  diag(my.metaweb) <- 1
  rownames(my.metaweb) <- sp.names
  colnames(my.metaweb) <- sp.names
  
  sp.traits$degree <- my.dist
  
  # -------------------------------------------------------------------------
  # third, assign species presences and absences based on landscape suitability
  landscape.n.rows <- richness*cells
  df.names <- expand.grid(sp.names,1:cells)
  landscape.names <- paste(df.names[,1],df.names[,2],sep="-")
  
  # presence-absence dataframe
  presence.df <- list()
  
  # create landscape matrix template
  landscape.template <- matrix(0,
                               nrow = landscape.n.rows,
                               ncol = landscape.n.rows,
                               dimnames = list(landscape.names,
                                               landscape.names))
  
  my.landscape.matrix <- landscape.template
  
  for(i.row in 1:nrow){
    for(i.col in 1:ncol){
      
      # this will be the local network, pruned from the metaweb
      my.network <- my.metaweb
      
      # cells are referred to by its row and col numbers, but I need
      # a cell number for the block matrix (see below)
      cell.id <- (nrow * (i.row - 1)) + i.col
      
      # suitability value of this cell
      landscape.value <- my.landscape[i.row,i.col]
      
      # presence dataframe for this cell
      cell.df <- tidyr::expand_grid(landscape.row = i.row,
                                    landscape.col = i.col,
                                    cell = cell.id,
                                    sp = sp.names,
                                    presence = FALSE)
      
      # go through all sp, checking if it should be present in this cell
      # and update the df and network matrix of the cell
      
      # a species is present if the suitability value of the cell
      # is within its mean +- sd 
      
      for(i.sp in 1:num.sp){
        my.suit <- which(suit.df$sp == sp.names[i.sp])
        min.suit <- suit.df$optimum[my.suit] - suit.df$sd[my.suit]
        max.suit <- suit.df$optimum[my.suit] + suit.df$sd[my.suit]
        
        if(landscape.value >= min.suit & landscape.value <= max.suit){
          cell.df$presence[i.sp] <- TRUE
        }else{
          # if species not present, "prune" the matrix, setting the species'
          # elements to 0
          my.network[sp.names[i.sp],] <- 0
          my.network[,sp.names[i.sp]] <- 0
        }# if-else species present
      }# for i.sp
      presence.df[[length(presence.df)+1]] <- cell.df
      
      # now, add the cell network to the landscape matrix
      init.row <- 1 + (num.sp * (cell.id - 1))
      init.col <- init.row
      end.row <- num.sp + (num.sp * (cell.id - 1))
      end.col <- end.row
      
      my.landscape.matrix[init.row:end.row,init.col:end.col] <- my.network
      
    }# i.col
  }# i.row
  
  presence.df <- bind_rows(presence.df)
  sp.presences <- presence.df %>% 
    group_by(sp) %>%
    summarise(presences = sum(presence))
  sp.traits <- left_join(sp.traits,sp.presences)
  
  # -------------------------------------------------------------------------
  # fourth, add dispersal links
  
  landscape <- my.landscape.matrix
  
  dispersal.values <- rexp(n = richness,rate = landscape.dispersal)
  sp.traits$dispersal.distance <- dispersal.values
  
  # expand the presence dataframe to see the cells to which
  # each sp can potentially disperse
  # by crossing presence information with cell distances
  # and with dispersal distances
  
  my.presence <- subset(presence.df,presence == TRUE) %>%
    dplyr::select(sp,cell) %>%
    rename(cell_from = cell)
  my.presence.full <- expand_grid(my.presence,cell_to = 1:cells)
  my.presence.full.2 <- left_join(my.presence.full,distance.df)
  my.presence.full.3 <- left_join(my.presence.full.2,sp.traits[,c("sp","dispersal.distance")])
  my.presence.full.3 <- subset(my.presence.full.3,cell_from != cell_to)
  
  my.presence.full.3$dispersal.potential <- ifelse(my.presence.full.3$distance <= 
                                                     my.presence.full.3$dispersal.distance,
                                                   TRUE,FALSE)
  my.dispersal.potential <- subset(my.presence.full.3[,c("sp","cell_from",
                                                         "cell_to","dispersal.potential")],
                                   dispersal.potential == TRUE)
  
  # now, check which of the target cells have populations of the species
  # because only those will have realized dispersal
  my.presence$presence <- TRUE
  my.dispersal.realized <- left_join(my.dispersal.potential,my.presence,
                                     by = c("cell_to" = "cell_from",
                                            "sp" = "sp")) %>%
    replace_na(list(presence = FALSE))
  my.dispersal.realized$dispersal <- as.logical(my.dispersal.realized$dispersal.potential * 
                                                  my.dispersal.realized$presence)
  my.dispersal.r2 <- subset(my.dispersal.realized[,c("sp","cell_from","cell_to","dispersal")],
                            dispersal == TRUE)
  
  # remove symmetrical information
  my.dispersal.r2$duplicated <- FALSE
  for(i in 1:nrow(my.dispersal.r2)){
    if(!my.dispersal.r2$duplicated[i]){
      dup <- which(my.dispersal.r2$cell_from == my.dispersal.r2$cell_to[i] &
                     my.dispersal.r2$cell_to == my.dispersal.r2$cell_from[i] &
                     my.dispersal.r2$sp == my.dispersal.r2$sp[i])
      my.dispersal.r2$duplicated[dup] <- TRUE
    }# if
  }# for i
  
  # this is the distilled dataframe
  # including the dispersal coefficient to assign
  # i.e. 1/number of dispersing cells (irrespective of distance)
  
  # this may change in the future
  
  my.dispersal <- subset(my.dispersal.r2,duplicated == FALSE) %>%
    group_by(sp,cell_from) %>%
    mutate(dispersal.coef = 1/n()) %>%
    dplyr::select(sp,cell_from,cell_to,dispersal.coef)
  
  # ----------------------------------------------------------------
  # update landscape matrix
  
  # I only need to go through the upper triangle, since 
  # dispersal is symmetric
  # hence the weird nested for loops
  for(i.row in 1:(nrow(landscape)-1)){
    
    # which sp?
    my.sp <- i.row %% richness
    if(my.sp == 0){my.sp <- richness}
    
    for(i.col in (i.row+1):ncol(landscape)){
      
      # if diagonal element and not main diagonal,
      # it is a dispersal coefficient
      # if(i.row %% richness == i.col %% richness & i.row != i.col){
      if(i.row %% richness == i.col %% richness){
        
        # source and dest cell
        source.cell <- ceiling(i.col/richness)
        dest.cell <- ceiling(i.row/richness)
        
        # double check
        if(source.cell != dest.cell){
          
          valid.dispersal <- which(my.dispersal$sp == sp.names[my.sp] &
                                     (my.dispersal$cell_from == source.cell &
                                        my.dispersal$cell_to == dest.cell | 
                                        my.dispersal$cell_from == dest.cell &
                                        my.dispersal$cell_to == source.cell) )
          if(length(valid.dispersal) == 1){
            
            disp.coef <- my.dispersal$dispersal.coef[valid.dispersal]
            
            # fill the symmetric positions
            landscape[i.row,i.col] <- disp.coef
            landscape[i.col,i.row] <- disp.coef
          }# if realized dispersal
          
        }# if different cell
      }# if dispersal cell
      
    }# for i.col
  }# for i.row
  
  
  # -----------------------------------------------------------------
  # remove absent species from every cell - this can be done simply
  # by keeping only columns and rows with 1 entries in their diagonal,
  # since each present species has 1 in the corresponding diagonal.
  
  present.sp <- diag(landscape) == 1
  landscape <- landscape[present.sp,present.sp] 
  
  return(list(landscape,sp.traits,distance.df))
  
}
