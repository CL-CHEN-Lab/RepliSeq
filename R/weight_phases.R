#' @title Weighting phases
#' @description Using a simulation of the S-phase (active forks in function of the percentage replicated genome), this function calculates a normalization factor for each phase based on the FACS gates and G1/G2 peaks
#'
#'
#' @param G1_peak aproximate position of the G1 peak 
#' @param G2_peak aproximate position of the G2 peak 
#' @param gates a data.frame containing start and end of each gate ordered as in rs_assay
#' @param rs_assay a Repli-seq assay (data.frame) loaded with readRS() and formatted as chr,start,stop,S1,...,Sx
#' @param SimSPhase a user constumise simulation of the S-phase progression (data.frame) formatted as Percentage of replicated DNA, Number of active forks
#' @param TotalReads Totala ammount of reads for normalization, by defaul 15000000 per phase
#' @return A  named vector containing the normalization factors for each phase to use with normalizeRS
#' @export
#'

weight_phases <- function(rs_assay, G1_peak,G2_peak,gates,SimSPhase=NA,TotalReads=NA) {
    
    #load default inputs if needed :
    if(is.na(SimSPhase)){
        load(system.file("extdata", "S-phase.RData", package = "RepliSeq"))
    }
    if(is.na(TotalReads)){
        TotalReads=length(gates[,1])*15000000
    }
    
    #convert gates positions into percent of DNA :
    gates=100*(gates-G1_peak)/(G2_peak-G1_peak)
    
    #calculate the mean number of forks inside each gate :
    gates$MeanForks=NA
    for (i in 1:length(gates[,1])){
        gates$MeanForks[i] = mean( SimSPhase$Forks[SimSPhase$percent > gates[i,1] & SimSPhase$percent < gates[i,2]])
    }
    
    MeanForks = TotalReads*gates$MeanForks/sum(gates$MeanForks)
    
    # calculate total reads for each fraction :
    Total_reads <- colSums(rs_assay[,!(names(rs_assay) %in% c("chr","start","stop"))])
    
    # return normalization factors if the number of phases and gates corresponds :
    if(length(Total_reads)==length(MeanForks)) {
        
        return(Total_reads/MeanForks)
    
    }else{
        
        stop('Gates and Number of Phases do not match')
    }
    

}