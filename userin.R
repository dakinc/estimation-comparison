###################################################
#The user is requried to have 3 packages installed#
###################################################

userin=function(pkg,x,y,p,nset,nlgt,nbeta){
  
  #############################################################
  #Purpose:
  #
  #An user interface to convert the data into desired format
  #
  #Note: The user should have bayesm,choicemodelr,RSGHB packages installed
  #
  #Author: 
  #   
  # Hanshui Zhang (Laste Modified on July 25, 2013)
  #
  #Arguments:
  #
  #   pkg: package indicator, the user could choose a package for the estimation by entering number 1, 2, or 3
  #       1 - bayesm package
  #       2 - ChoiceModelR package
  #       3 - RSGHB package
  #       4 - Mixed Probit
  #
  #   p: num of alternatives
  #   nset: num of choice sets
  #   nlgt: num of respodents
  #   nbeta: num of model parameters to be estimated
  #
  #x: design matrix which contains the attributes
  #   if the user is using bayesm/RSGHB, effect coding should be done before entering x
  #
  #y: (nset*nlgt) by 1 matrix, this matrix contains the responses made by the respondents
  #
  #######################################################################
  
  
  
  if(pkg==1){
    
    library(bayesm)
    print("The dataset is ready to be used for bayesm package")
    
    
    ###########################################
    ###an input procedure for bayesm package###
    ###########################################
    
    bayesmin=function(x,y,p,nset,nlgt){
      
      #transform y and x into matrix form
      y=t(t(y))
      x=t(t(x))
      
      lgtdata=NULL
      
      ni=rep(nset,nlgt)
      csa=nset*p
      for (i in 1:nlgt)
      { 
        #obtain y
        ychoice=NULL
        ybeg=nset*(i-1)+1
        yend=nset*i
        for(c in 1:nset){ychoice[1:nset]=y[ybeg:yend]}
        
        #transform x into dataframe
        xmat=NULL
        xbeg=csa*(i-1)+1
        xend=csa*i
        xmat[[i]]=x[xbeg:xend,]
        
        lgtdata[[i]]=list(y=ychoice,X=xmat[[i]])
      }
      #end of converting data
      
      #the bayesmin function returns a list of 2
      return(bayesmdata=list(p=p,lgtdata=lgtdata))
      
    }#end of bayesmin function
    
    #if user call pkg 1, the data will be converted to bayesm format
    return(bayesmin(x,y,p,nset,nlgt))
    
    
  }#end of calling bayesm package
  
  else if(pkg==2){
    
    library(ChoiceModelR)
    print("The dataset is ready to be used for ChoiceModelR package")
    
    #################################################
    ###an input procedure for ChoiceModelR package###
    #################################################
    
    choicemodelrin=function(x,y,p,nset,nlgt){
      
      #transform y and x into matrix form
      y=t(t(y))
      x=t(t(x))
      
      #input procedure for ChoiceModelR
      #setup choiceset column for each individual
      #choiceset row, total num of rows=p*nset*nlgt
      set=rep(1:nset, each = p, times=nlgt)
      
      #id row, total number of rows=p*nset*nlgt
      id=rep(1:nlgt,each=nset*p)
      
      #create alternative indicator, total number of rows=p*nset*nlgt
      alt=rep(c(1:p), nset*nlgt)
      
      #beginning of the matrix
      initialmat=t(rbind(id, set,alt))
      
      
      #combine xmat and attrmat
      xmat=cbind(initialmat,x)
      
      #make choice columns
      newchoice=y
      zeromat=matrix(0,nset*nlgt,p-1)
      choicemat=cbind(newchoice,zeromat)
      
      #This is the final y column representing choice
      choicecol=matrix(c(t(choicemat)))
      
      #bind everything together
      return(choicemodelrdata=cbind(xmat,choicecol))
      
      
    }#end of choicemodelrin function
    
    #if user call pkg 2, the data will be converted to ChoiceModelR format
    return(choicemodelrin(x,y,p,nset,nlgt))
    
  }#end of calling ChoiceModelR
  
  else if (pkg==3){
    library(RSGHB)
    print("The dataset is ready to be used for RSGHB package")
   
    ##########################################
    ###an input procedure for RSGHB package###
    ##########################################
    
    rsghbin=function(x,y,p,nset,nlgt,nbeta){
      
      #transform y and x into matrix form
      y=t(t(y))
      x=t(t(x))
      
      #nbeta is number of model parameters
      nbeta=ncol(x)
      
      #RSGHB id row, total number of rows=p*nset*nlgt
      rsghbid=rep(1:nlgt,each=nset)
      
      #define choiceset column
      ncs=rep(1:nset, times=nlgt)
      
      #first 2 columns of the RSGHB data structure 
      rsghbinitialmat=t(rbind(rsghbid,ncs))
      
      #attribute matrix
      rsghbattrmat=NULL
      indset=nset*nlgt
      for (cs in 1:indset){
        beg=p*(cs-1)+1
        end=p*cs
        
        xtemp=NULL
        for(col in 1:nbeta){
          xtemp=cbind(xtemp,t(x[beg:end,col]))
        }
        rsghbattrmat=rbind(rsghbattrmat,xtemp)
      }
      
      #choice column
      RSGHBchoice=y
      
      #combine columns to form the dataset
      RSGHBdata=data.frame(cbind(rsghbinitialmat,rsghbattrmat,RSGHBchoice))
      
      
      #first column is named as "ID", second column is called "choice set"
      colnames(RSGHBdata)[[1]]="ID"
      colnames(RSGHBdata)[[2]]="Choice Set"
      
      #name the last column (choice column) of the data structure
      cy=ncol(RSGHBdata)
      colnames(RSGHBdata)[[cy]]="Choice"
      
      #the RSGHB data is returned as a data frame
      return(RSGHBdata)
      
    }#end of rsghbin function
    
    #if user call pkg 3, the data will be converted to RSGHB format
    return(rsghbin(x,y,p,nset,nlgt,nbeta))
    
  }#end of calling RSGHB
  
  
  else if (pkg==4){
    print("The dataset is ready to be used for Mixed Probit Estimation")
    library(bayesm)
    
    ##########################################
    ###an input procedure for rmxp function###
    ##########################################
    mxpin=function(x,y,p,nset,nlgt,nbeta){
      
      x=t(t(x))
      y=t(t(y))
      
      #mxpin function transform the input data into the structure which can be used in rmxp function
      #code y matrix using indicator variables
      ynum=nrow(y)
      yind=NULL
      for (i in 1:ynum){
        zerotemp=matrix(0,p,1)
        index=y[i,]
        zerotemp[index]=1
        yind=rbind(yind,zerotemp)
      }#end for i
      
      #format y using the coded y matrix
      y = array(t(yind),dim=c(p,nset,nlgt))
      
      X = array(t(x),dim=c(nbeta,p,nset,nlgt))
      
      return(Data=list(y=y,X=X,nlgt=nlgt,nset=nset,p=p,nbeta=nbeta))
      
    }# end of mxpin
    
    return(mxpin(x,y,p,nset,nlgt,nbeta))
    
  }# end else if pkg==4
  
  
  else{
    
    return(print("please specify: 1-bayesm, 2-ChoiceModelR, 3-RSGHB, 4-Mixed Probit") )
    
  }#end else
  
  
}#end of userin function