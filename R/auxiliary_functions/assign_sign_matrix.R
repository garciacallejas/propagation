
#' Generate a sign matrix considering all interaction types
#' This function is intended for the study on Spatial MetaEcosystems (SME).
#' Interactios can be assigned, for each type, randomly or following a structure.
#' For details of how are they structured, see the function InteractionTopology_SME
#' 
#' @param S number of speices
#' @param connectances vector with connectances of each type: amensalism,antagonism,commensalism,competition,mutualism
#' @param amensalism structured or random
#' @param antagonism structured or random
#' @param commensalism structured or random
#' @param competition structured or random
#' @param mutualism structured or random
#' @param connectance.type directed or undirected
#' @param niche numeric vector of length S representing a 1-dimensional niche axis used to generate network structure
#'
#' @return sign matrix
#' @export
#'
#' @examples
assign_sign_matrix <- function(S = 10,
                           connectances = rep(0.1,5),
                           amensalism = "random",
                           antagonism = "random",
                           commensalism = "random",
                           competition = "random",
                           mutualism = "random",
                           connectance.type = "undirected",
                           niche = 0){
  
  sign.matrix <- matrix(0,nr=S,nc=S)
  
  amensalism.matrix = matrix(0,nr=S,nc=S)
  antagonism.matrix = matrix(0,nr=S,nc=S)
  commensalism.matrix = matrix(0,nr=S,nc=S)
  competition.matrix = matrix(0,nr=S,nc=S)
  mutualism.matrix = matrix(0,nr=S,nc=S)
  
  # if any interaction type is structured, generate 1-dim niches and ranges
  # these niches will be assumed as determining trophic position
  # if(antagonism == "structured" |
  #    amensalism == "structured" |
  #    commensalism == "structured" |
  #    competition == "structured" |
  #    mutualism == "structured"){
  #   
  #   # Niche position
  #   niche = runif(S,0,1)
  #   names(niche) <- 1:length(niche)
  # }
  # if any interaction is structured, niche needs to be provided
  if(antagonism == "structured"){
    
    # Range
    # connectances is ordered alphabetically:
    # amensalism, antagonism, commensalism, competition, mutualism
    beta = 1/(2*connectances[2])-1
    range = rbeta(S,1,beta)*niche
    
    #Centroid
    centroid = runif(S,range/2,niche)
    
    # Make the first niche a producer
    range[which.min(niche)] = 0
    
    # Evaluate the matrix of links
    
    low = centroid - range/2
    high = centroid + range/2
    
    for(s1 in 1:S){
      for(s2 in 1:S){
        if(low[s1] < niche[s2] && high[s1] > niche[s2]){ 
          antagonism.matrix[s2,s1] = 1    	  
        }# if
      }# for
    }# for
    
  }else{
    rand <- matrix(runif(S^2,0,1),nr=S,nc=S)
    antagonism.matrix[rand<connectances[2]*2] = 1
  }# antagonism
  if(amensalism == "structured"){
    # amensalism.matrix <- InteractionTopology_SME(S = S,
    #                                           connectance = connectances[1],
    #                                           niches = niche,
    #                                           mode = "same",
    #                                           forbidden.positions = antagonism.matrix,
    #                                           connectance.type = connectance.type)
  }else{
    for(i.sp in 1:S){
      for(j.sp in 1:S){
        if(j.sp<i.sp){
          rand <- runif(1,0,1)
          if(rand<connectances[1]){
            amensalism.matrix[i.sp,j.sp] <- 1
            amensalism.matrix[j.sp,i.sp] <- 1
          }#if
        }#if
      }#for j
    }#for i
  }#amensalism
  if(commensalism == "structured"){
    # commensalism.matrix <- InteractionTopology_SME(S = S,
    #                                             connectance = connectances[3],
    #                                             niches = niche,
    #                                             mode = "same",
    #                                             forbidden.positions = antagonism.matrix,
    #                                             connectance.type = connectance.type)
  }else{
    for(i.sp in 1:S){
      for(j.sp in 1:S){
        if(j.sp<i.sp){
          rand <- runif(1,0,1)
          if(rand<connectances[3]){
            commensalism.matrix[i.sp,j.sp] <- 1
            commensalism.matrix[j.sp,i.sp] <- 1
          }#if
        }#if
      }#for j
    }#for i
  }#commensalism
  if(competition == "structured"){
    # competition.matrix <- InteractionTopology_SME(S = S,
    #                                            connectance = connectances[4],
    #                                            niches = niche,
    #                                            mode = "same",
    #                                            forbidden.positions = antagonism.matrix,
    #                                            connectance.type = connectance.type)
  }else{
    for(i.sp in 1:S){
      for(j.sp in 1:S){
        if(j.sp<i.sp){
          rand <- runif(1,0,1)
          if(rand<connectances[4]){
            competition.matrix[i.sp,j.sp] <- 1
            competition.matrix[j.sp,i.sp] <- 1
          }#if
        }#if
      }#for j
    }#for i
  }#competition
  if(mutualism == "structured"){
    # mutualism.matrix <- InteractionTopology_SME(S = S,
    #                                          connectance = connectances[5],
    #                                          niches = niche,
    #                                          mode = "adjacent",
    #                                          forbidden.positions = antagonism.matrix,
    #                                          connectance.type = connectance.type)
  }else{
    for(i.sp in 1:S){
      for(j.sp in 1:S){
        if(j.sp<i.sp){
          rand <- runif(1,0,1)
          if(rand<connectances[5]){
            mutualism.matrix[i.sp,j.sp] <- 1
            mutualism.matrix[j.sp,i.sp] <- 1
          }#if
        }#if
      }#for j
    }#for i
  }#mutualism
  
  # check if any position has two or more links, and assign randomly one of them to the final matrix. If it's not random,
  # priority would be given to the first interactions assigned
  # antagonism doesn't count, these do have priority
  
  for(i.sp in 1:S){
    for(j.sp in 1:S){
      existing.link <- rep(0,5)
      if(amensalism.matrix[i.sp,j.sp] == 1 | amensalism.matrix[j.sp,i.sp] == 1){
        existing.link[1] <- 1
      }
      if(antagonism.matrix[i.sp,j.sp] == 1 | antagonism.matrix[j.sp,i.sp] == 1){
        existing.link[2] <- 1
      }
      if(commensalism.matrix[i.sp,j.sp] == 1 | commensalism.matrix[j.sp,i.sp] == 1){
        existing.link[3] <- 1
      }
      if(competition.matrix[i.sp,j.sp] == 1 | competition.matrix[j.sp,i.sp] == 1){
        existing.link[4] <- 1
      }
      if(mutualism.matrix[i.sp,j.sp] == 1 | mutualism.matrix[j.sp,i.sp] == 1){
        existing.link[5] <- 1
      }
      
      # if any existing link 
      if(sum(existing.link == 1)>0){
        # choose its type:
        # if 1)antagonism, or 2)more than one link, or 3)just one
        if(existing.link[2] == 1){
          link.type <- 2
        }else if(sum(existing.link == 1)>1){
          link.type <- sample(which(existing.link == 1),1)
        }else{
          link.type <- which(existing.link == 1)
        }
        
        # add the link to the sign matrix
        # amensalist and commensalist links are randomly directed
        if(link.type == 1){#amensalist
          affected.sp <- sample(1:2,1)
          if(affected.sp == 1){
            sign.matrix[i.sp,j.sp] <- -1
            sign.matrix[j.sp,i.sp] <- 0
          }else{
            sign.matrix[i.sp,j.sp] <- 0
            sign.matrix[j.sp,i.sp] <- -1
          }
        }else if(link.type == 2){#antagonist
          if(antagonism.matrix[i.sp,j.sp] == 1){
            sign.matrix[i.sp,j.sp] <- -1
            sign.matrix[j.sp,i.sp] <- 1
          }else{
            sign.matrix[i.sp,j.sp] <- 1
            sign.matrix[j.sp,i.sp] <- -1
          }
        }else if(link.type == 3){#commensalist
          affected.sp <- sample(1:2,1)
          if(affected.sp == 1){
            sign.matrix[i.sp,j.sp] <- 1
            sign.matrix[j.sp,i.sp] <- 0
          }else{
            sign.matrix[i.sp,j.sp] <- 0
            sign.matrix[j.sp,i.sp] <- 1
          }
        }else if(link.type == 4){#competitive
          sign.matrix[i.sp,j.sp] <- -1
          sign.matrix[j.sp,i.sp] <- -1
        }else if(link.type == 5){#mutualistic
          sign.matrix[i.sp,j.sp] <- 1
          sign.matrix[j.sp,i.sp] <- 1
        }
      }# if any existing link
    }# for j
  }# for i
  
  # Intraspecific density-dependence
  diag(sign.matrix) <- -1
  
  return(sign.matrix)
  
}

