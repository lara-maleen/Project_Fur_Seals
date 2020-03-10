adjacent <- function(x,y){
  dists <- outer(1:length(x),1:length(x),FUN = function(i,j) sqrt((x[i]-x[j])^2+(y[i]-y[j])^2))
  apply(dists,1,FUN = function(x) any(x > 0 & x <= sqrt(2)) | sum(x==0) > 1)
}

adjacent_n <- function(x,y,n){
  dists <- outer(1:length(x),1:length(x),FUN = function(i,j) sqrt((x[i]-x[j])^2+(y[i]-y[j])^2))
  apply(dists,1,FUN = function(x) sum(x <= sqrt(2)) > n)
}

grow_islands <- function(world.mat,x0,y0,maxsurf=5e2){
  if(any(adjacent(x0,y0))){
    stop("initial points cannot be adjacent")
  }
  
  world.mat[] <- 0
  world.mat[cbind(x0,y0)] <- 1:length(x0)
  
  isles <- rep(1,length(x0))
  no_change <- 0
  while(no_change < 20){
    
    cands <- which(rbind(FALSE,(world.mat[-1,] == 0 & world.mat[-nrow(world.mat),] != 0)) |
                     rbind((world.mat[-nrow(world.mat),] == 0 & world.mat[-1,] != 0),FALSE) |
                     cbind(FALSE,(world.mat[,-1] == 0 & world.mat[,-ncol(world.mat)] != 0)) |
                     cbind((world.mat[,-ncol(world.mat)] == 0 & world.mat[,-1] != 0),FALSE))
    
    cand <- NA
    # pick candidate
    if(length(cands) > 1){
      cand <- sample(cands,1)
    }else if(length(cands == 1)){
      cand <- cands
    }
    
    id.neigh <- NA
    # cat(cand,"--")
    if(!is.na(cand)){
      # number of different islands it is adjacent too
      y.cand <- (cand-1) %% nrow(world.mat) + 1
      x.cand <- (cand - y.cand)/nrow(world.mat) + 1
      
      neighbors <- data.frame(x=c(x.cand + 1, x.cand - 1 , x.cand, x.cand),y=c(y.cand,y.cand,y.cand+1,y.cand-1))
      neighbors <- neighbors[neighbors$x >0 & neighbors$y > 0 & neighbors$x <=ncol(world.mat) & neighbors$y <= nrow(world.mat),]
      id.neigh <- as.numeric(world.mat[as.matrix(neighbors[,2:1])])
      id.neigh <- unique(id.neigh)
      id.neigh <- id.neigh[id.neigh!=0]
    }
    
    if(!is.na(cand) & length(id.neigh) == 1 & isles[id.neigh[1]] < maxsurf){
      world.mat[y.cand,x.cand] <- id.neigh
      isles[id.neigh] <- isles[id.neigh] + 1
      no_change <- 0
    }else{
      no_change <- no_change + 1
    }
    
  }  
  world.mat
}

draw <- function(animals,world.lan,xlim,ylim){
  par(bg=type_col[["water"]],mar=c(0,0,0,0))
  plot(0,0,type="n",xlim=xlim,ylim=ylim,axes=FALSE,xlab="",ylab="")    
  sapply(1:nrow(world.lan),FUN = function(x){
    with(world.lan[x,],polygon(c(x-0.5,x-0.5,x+0.5,x+0.5),c(y-0.5,y+0.5,y+0.5,y-0.5),col = type_col[[type]],border = NA))
  })
  
  return(NULL)
}

sketch <- function(edge,xlim,ylim){
  plot(0,0,type="n",axes=FALSE,xlim=xlim,ylim=ylim,xlab="",ylab="")
  lapply(edge,FUN = function(x) points(x[,1],x[,2],col="grey",pch=20))  
}

world.edge <- function(world.mat){
  isles <- 1:max(world.mat)
  out <- list()
  for(i in isles){
    # getting all boxeswith ID isles[i] that are next to a pieces of water
    cands <- which(rbind(FALSE,(world.mat[-1,] == 0 & world.mat[-nrow(world.mat),] == isles[i])) |
                     rbind((world.mat[-nrow(world.mat),] == 0 & world.mat[-1,] == isles[i]),FALSE) |
                     cbind(FALSE,(world.mat[,-1] == 0 & world.mat[,-ncol(world.mat)] == isles[i])) |
                     cbind((world.mat[,-ncol(world.mat)] == 0 & world.mat[,-1] == isles[i]),FALSE))
    
    y.cand <- (cands-1) %% nrow(world.mat) + 1
    x.cand <- (cands - y.cand)/nrow(world.mat) + 1
    for(j in 1:5){
      online <- adjacent_n(x.cand,y.cand,2)
      x.cand <- x.cand[online]
      y.cand <- y.cand[online]
    }
    # removing those pieces of water that form a block of smaller than 3 pieces of water
    out[[i]] <- cbind(x.cand,y.cand)
  }
  out
}

plot_animals <- function(animals,col){
  points(animals$x,animals$y,cex=animals$size/5,pch=18,col=col)  
}

move_moms <- function(moms,speed = 1){
  # for all 0 and negative state moms, move
  ind <- (moms$state == -1)
  dir.x <- moms$xof[ind] - moms$x[ind]
  dir.y <- moms$yof[ind] - moms$y[ind]
  nn <- sqrt(dir.x^2 + dir.y^2)
  dir.x <- dir.x/nn
  dir.y <- dir.y/nn
  moms$x[ind] <- moms$x[ind] + speed*dir.x + rnorm(sum(ind),0,0.7)
  moms$y[ind] <- moms$y[ind] + speed*dir.y + rnorm(sum(ind),0,0.7)
  
  ind <- (moms$state == 0)
  dir.x <- moms$xin[ind] - moms$x[ind]
  dir.y <- moms$yin[ind] - moms$y[ind]
  nn <- sqrt(dir.x^2 + dir.y^2)
  dir.x <- dir.x/nn#*sign(moms$state[ind]+0.1)
  dir.y <- dir.y/nn#*sign(moms$state[ind]+0.1)
  moms$x[ind] <- moms$x[ind] + speed*dir.x + rnorm(sum(ind),0,0.7)
  moms$y[ind] <- moms$y[ind] + speed*dir.y + rnorm(sum(ind),0,0.7)
  
  ind <- (moms$state == 2)
  moms$x[ind] <- moms$x[ind] + rnorm(sum(ind),0,0.7)
  moms$y[ind] <- moms$y[ind] + rnorm(sum(ind),0,0.7)
  
  moms
}

run_game <- function(init_pop,xsize,ysize,world,Nt=1,pause=1,live=FALSE){
  edge <- world[['edge']]
  moms <- init_pop[['moms']]
  offs <- init_pop[['offs']]
  
  for(t in 1:Nt){
    
    if(live){
      Sys.sleep(pause)
    }
    #draw(NA,world.lan,xlim=range(world.lan$x)+c(-3,3),ylim=range(world.lan$y)+c(-3,3))
    sketch(edge,xlim=c(1,xsize)+c(-3,3),ylim=c(1,ysize)+c(-3,3))
    
    
    plot_animals(moms[moms$state != 1,],col="red")
    plot_animals(offs,col="black")
    
    if(t<Nt){
      # moms move
      moms <- move_moms(moms)
      
      # mom state update
      # -1 :moving towards foraging, 0: moving back to young, 1: foraging at sea, 2: feeding pup
      # wait: number of time steps to wait in the current feeding or foraging state
      moms$wait[moms$state >= 1] <- moms$wait[moms$state >= 1] - 1
      
      # feeding to moving out
      moms$state[moms$state == 2 & moms$wait == 0] <- -1
      
      # moving out to foraging
      inds <- moms$state==-1 & sqrt((moms$x - moms$xof)^2+(moms$y - moms$yof)^2) < 1
      moms$state[inds] <- 1
      moms$wait[inds] <- 5
      
      # foraging to moving in
      moms$state[moms$state== 1 & moms$wait == 0] <- 0
      
      # moving in to feeding
      inds <- moms$state == 0 & sqrt((moms$x - moms$xin)^2+(moms$y - moms$yin)^2) < 1
      moms$state[inds] <- 2
      moms$wait[inds] <- 5
      
      
      # feeding takes place (detect range between pups and moms)
      growers <- (sqrt((offs$x - moms$x)^2 + (offs$y - moms$y)^2) < 1) & moms$state == 2 
      
      if(any(growers)){
        offs$size[growers] <- offs$size[growers] + 0.5*(20 - offs$size[growers]) 
      }
      
      # pups potentially get trampled (detect range between pups and non-moms)
      dist <- outer(1:nrow(moms),1:nrow(moms), FUN = function(i,j){
        sqrt((moms$x[j] - offs$x[i])^2 + (moms$y[j]-offs$y[i])^2)
      })
      
      diag(dist) <- 20 #hack, so the own mother cannot trample the pup
      tramples <- apply(dist,1,FUN = function(x) any(x < 3))
      if(any(tramples)){
        moms <- moms[!tramples,]
        offs <- offs[!tramples,]
      }
      # size is adapted for the pups that were not fed - small shrinking
      offs$size <- offs$size - 0.3
      
      # some pups die of hunger
      cand <- offs$size <= 0
      if(sum(cand)>0){
        moms <- moms[!cand,]
        offs <- offs[!cand,]
      }
    }
  }
  
  return(NA)
}


create_world <- function(xsize,ysize,maxsurf=5e2,x0,y0){
  
  world.mat <- matrix(NA,nrow=xsize,ncol=ysize)
  
  if(missing(x0) | missing(y0)){
    x0 <- c(round(0.2*xsize),round(0.75*ysize))
    y0 <- c(round(0.75*xsize),round(0.45*ysize))
  }
  
  world.mat <- grow_islands(world.mat,x0,y0,maxsurf = maxsurf)
  
  edge <- world.edge(world.mat)
  
  
  return(list(edge=edge,world.mat=world.mat))
}

create_init_pop <- function(world,x0,y0){
  edge <- world[['edge']]
  world.mat <- world[['world.mat']]
  
  if(missing(x0) | missing(y0)){
    offs <- data.frame(x=sample(ncol(world.mat),50),y=sample(nrow(world.mat),50),size=10)
  }else{
    offs <- data.frame(x=x0,y=y0,size=10)
  }
  # island check:
  offs <- offs[as.numeric(world.mat[as.matrix(offs[,2:1])]) != 0,]
  
  moms <- offs
  moms$size <- 20
  moms$xin <- moms$x
  moms$yin <- moms$y
  moms$wait <- 0
  moms$xof <- NA
  moms$yof <- NA
  moms$state <- -1 # -1 :moving towards foraging (other states: 0,1,2)
  
  for(i in 1:nrow(moms)){
    m.x <- moms$x[i]
    m.y <- moms$y[i]
    ed <- edge[[world.mat[m.y,m.x]]]
    point <- which.min((m.x-ed[,1])^2 + (m.y-ed[,2])^2)
    moms$xof[i] <- ed[point,1]
    moms$yof[i] <- ed[point,2]
  }
  
  return(list(offs=offs,moms=moms))
}

run_game_pdf <- function(...,filename){
  pdf(filename)
  run_game(...)
  dev.off()  
}

run_game_live <- function(...){
  run_game(...,live=TRUE)
}

