fuzzy.p.value <-
function(t, H0b, sig=0.05, p.value, knot.n=10, fig=c("1","2","3"), ...){

# options(warn = -1) 

  if(fig==1){
    P = f2apply(t, H0b, p.value, knot.n=knot.n, type="l", I.O.plot=FALSE, ...)
    } 
    else{
      if(fig==2){
        P = f2apply(t, H0b, p.value, knot.n=knot.n, type="l", I.O.plot=FALSE, ...)

        if(class(sig)=="numeric"){abline(v=sig, col=2, lty=3)}
           else{plot(sig, col=2, lty=3, add=TRUE)}
        legend("topright", c("Fuzzy p-value", "Significance level"), col = c(1,2), text.col = 1, lwd = c(1,1), lty = c(1,3), bty="n")
        }
        else{
            if(fig==3){
              P = f2apply(t, H0b, p.value, knot.n=knot.n, type="l", I.O.plot=TRUE, ...)
              x = t
              y = H0b
              } 

            }
       }


  if( class(sig) == "numeric" ){
         sig <- TriangularFuzzyNumber(sig, sig, sig)
         }


  P_L = P$cuts[,"L"]
  P_L = P_L[length(P_L):1]

  P_U = P$cuts[,"U"]
  P_U = P_U[length(P_U):1]

  S_L = alphacut(sig, round(seq(0, 1, len=knot.n), 5))[,"L"]

  S_U = alphacut(sig, round(seq(0, 1, len=knot.n), 5))[,"U"]


  Int1 = ( P_U - S_L ) * ( P_U > S_L )
  Int2 = ( P_L - S_U ) * ( P_L > S_U )

  Arz = 1 / (knot.n - 1)  #Arze Mostatilha baraye mohasebe-ye Integral

  Integral1 <- ( sum( Int1 ) - Int1[1]/2 - Int1[length(Int1)]/2 ) *Arz
  Integral2 <- ( sum( Int2 ) - Int2[1]/2 - Int2[length(Int2)]/2 ) *Arz

  Delta_PS = Integral1 + Integral2

  Int3 = ( S_U - P_L ) * ( S_U > P_L )
  Int4 = ( S_L - P_U ) * ( S_L > P_U )

  Integral3 <- ( sum( Int3 ) - Int3[1]/2 - Int3[length(Int3)]/2 ) *Arz
  Integral4 <- ( sum( Int4 ) - Int4[1]/2 - Int4[length(Int4)]/2 ) *Arz

  Delta_SP = Integral3 + Integral4


  Degree_P_biger_than_S = Delta_PS / (Delta_PS + Delta_SP)
  Degree_S_biger_than_P = 1- Degree_P_biger_than_S


  if (Degree_P_biger_than_S >= Degree_S_biger_than_P) 
    {
    result = noquote( paste( "The null hypothesis (H0) is accepted with degree D(P>S)=", round(Degree_P_biger_than_S, 4), ", at  the considered significance level." ) ) 
    accepted_hypothesis = noquote( "H0" )
    acceptance_degree = Degree_P_biger_than_S
    }
    else{if(Degree_P_biger_than_S < Degree_S_biger_than_P) 
           {
    result = noquote( paste( "The althernative hypothesis (H1) is accepted with degree D(S>P)=", round(Degree_S_biger_than_P, 4), ", at  the considered significance level." ) ) 
    accepted_hypothesis = noquote( "H1" )
    acceptance_degree = Degree_S_biger_than_P
           }
           else{print( noquote( paste0("Impossible case" ) ) )}
         }

  return( list( 
               result = result, 
               cuts = P$cuts, 
               core = P$core , 
               support = P$support , 
               Delta_PS = as.numeric( Delta_PS ) , 
               Delta_SP = as.numeric( Delta_SP ) , 
               Degree_P_biger_than_S = as.numeric( Degree_P_biger_than_S ) , 
               Degree_S_biger_than_P = as.numeric( Degree_S_biger_than_P ) ,
               accepted_hypothesis = accepted_hypothesis ,
               acceptance_degree = as.numeric( acceptance_degree )
               ) 
        )

}
