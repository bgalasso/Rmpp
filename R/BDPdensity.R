### BDPdensity.R
### For density estimation using a Bernstein-Dirichlet prior
###
### Copyright: Alejandro Jara and Fernando Quintana, 2007-2012.
###
### Last modification: 30-08-2010.
###
### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation; either version 2 of the License, or (at
### your option) any later version.
###
### This program is distributed in the hope that it will be useful, but
### WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
### General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, write to the Free Software
### Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
###
### The authors' contact information:
###
###      Alejandro Jara
###      Department of Statistics
###      Facultad de Matematicas
###      Pontificia Universidad Catolica de Chile
###      Casilla 306, Correo 22
###      Santiago
###      Chile
###      Voice: +56-2-3544506  URL  : http://www.mat.puc.cl/~ajara
###      Fax  : +56-2-3547729  Email: atjara@uc.cl
###
###      Fernando Quintana
###      Departamento de Estadistica
###      Facultad de Matematicas
###      Pontificia Universidad Catolica de Chile
###      Casilla 306, Correo 22
###      Santiago
###      Voice: +56-2-3544464  URL  : http://www.mat.puc.cl/~quintana
###      Fax  : +56-2-3547229  Email: quintana@mat.puc.cl


BDPdensity<-function(y,support=3,ngrid=1000,grid=NULL,prior,mcmc,state,status,data=sys.frame(sys.parent()),na.action=na.fail)
UseMethod("BDPdensity")

BDPdensity.default<-function(y,support=3,ngrid=1000,grid=NULL,prior,mcmc,state,status,data,na.action=na.fail)
{
         #########################################################################################
         # call parameters
         #########################################################################################
           cl <- match.call()
		   resp <- na.action(as.matrix(y))
	       varnames <- all.vars(cl)[1]

         #########################################################################################
         # data structure
         #########################################################################################
     	   nrec <- nrow(resp)
		   nvar <- ncol(resp)

           if(nvar>1)
           {
             stop("This function can only be used for univariate density estimation.\n")
           }

           resp <- as.vector(resp)


           if(is.null(grid))
           {
			  grid <- seq(0,1,len=ngrid)
		   }
		   else
		   {
			  grid <- as.vector(grid)
			  ngrid <- length(grid)
		   }

           if(support==3)
           {
              left<-min(resp)-0.5*sqrt(var(resp))
              right<-max(resp)+0.5*sqrt(var(resp))
              x<-(resp-left)/(right-left)
              jacob<-1.0/(right-left)
              grids<- left+grid*(right-left)
           }
           if(support==1)
           {
              left<-0
              right<-1
              x<-resp
              jacob<-1
              grids <- grid
           }
           if(support==2)
           {
              left<-0
              right<-max(resp)+0.5*sqrt(var(resp))
              x<-(resp-left)/(right-left)
              jacob<-1.0/right
              grids <- grid*right
           }

         #########################################################################################
         # prior information
         #########################################################################################

		   if(is.null(prior$aa0))
		   {
				aa0<--1
				ab0<--1
				alpha<-prior$alpha
				alpharand<-0
		   }
           else
           {
				aa0<-prior$aa0
				ab0<-prior$ab0
				alpha<-1
				alpharand<-1
		   }
		   a0b0<-c(aa0,ab0)

           kmax<-prior$kmax
           a0<-prior$a0
           b0<-prior$b0


         #########################################################################################
         # mcmc specification
         #########################################################################################
           mcmcvec<-c(mcmc$nburn,mcmc$nskip,mcmc$ndisplay)
           nsave<-mcmc$nsave

         #########################################################################################
         # output
         #########################################################################################

           fun<-matrix(0,nrow=ngrid,ncol=nsave)
           fhat<-rep(0,ngrid)
           cpo<-matrix(0,nrow=nrec,ncol=2)
           thetasave<-matrix(0,nrow=nsave,ncol=3)
           randsave<-matrix(0,nrow=nsave,ncol=(nrec+1))

         #########################################################################################
         # parameters depending on status
         #########################################################################################

    	   if(status==TRUE)
		   {
                  yclus<-rbeta(nrec,a0,b0)
                  k<-kmax
                  ncluster<-nrec
                  ss<-seq(1,nrec)
		   }

      	   if(status==FALSE)
		   {
	          alpha<-state$alpha
                  k<-state$k
	          ncluster<-state$ncluster
	          yclus<-state$yclus
	          ss<-state$ss
		   }

         #########################################################################################
         # working space
         #########################################################################################
           cstrt<-matrix(0,nrow=nrec,ncol=nrec)
           ccluster<-rep(0,nrec)
           prob<-rep(0,nrec+1)
           probk<-rep(0,kmax)
           seed<-c(sample(1:29000,1),sample(1:29000,1))
           y<-rep(0,nrec)

         #########################################################################################
         # calling the fortran code
         #########################################################################################

          foo <- .Fortran("bdpdensity",
					nrec       =as.integer(nrec),
					jacob      =as.double(jacob),
					x          =as.double(x),
					a0b0       =as.double(a0b0),
					a0         =as.double(a0),
					b0         =as.double(b0),
					kmax       =as.integer(kmax),
					k          =as.integer(k),
					ncluster   =as.integer(ncluster),
					ss         =as.integer(ss),
					alpha      =as.double(alpha),
					yclus	     =as.double(yclus),
					mcmc       =as.integer(mcmcvec),
					nsave      =as.integer(nsave),
					cpo        =as.double(cpo),
					randsave   =as.double(randsave),
					thetasave  =as.double(thetasave),
					ngrid      =as.integer(ngrid),
					grid       =as.double(grid),
					fun        =as.double(fun),
					seed       =as.integer(seed),
					cstrt      =as.integer(cstrt),
					ccluster   =as.integer(ccluster),
					prob       =as.double(prob),
					probk      =as.double(probk),
					y          =as.double(y),
					PACKAGE    ="Rmpp")

         #########################################################################################
         # save state
         #########################################################################################

           cpom<-matrix(foo$cpo,nrow=nrec,ncol=2)
           cpo<-cpom[,1]
           fso<-cpom[,2]

           model.name<-"Bayesian semiparametric density estimation"

           state <- list(alpha=foo$alpha,
                         yclus=foo$yclus,
                         ncluster=foo$ncluster,
                         ss=foo$ss)

           randsave<-matrix(foo$randsave,nrow=nsave,ncol=(nrec+1))
           thetasave<-matrix(foo$thetasave,nrow=nsave,ncol=3)

           pnames<-c("k","ncluster","alpha")
           colnames(thetasave)<-pnames
           coeff<-apply(thetasave,2,mean)

           save.state <- list(thetasave=thetasave,randsave=randsave)

		   z<-list(	call=cl,
					y=resp,
					varnames=varnames,
					cpo=cpo,
					fso=fso,
					modelname=model.name,
					coefficients=coeff,
					prior=prior,
					mcmc=mcmc,
					state=state,
					save.state=save.state,
					nrec=foo$nrec,
					grid=grids,
					fun=foo$fun,
					fhat=foo$fhat)

          cat("\n\n")
		  class(z)<-"BDPdensity"
		  return(z)
}




