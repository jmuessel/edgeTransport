#' Harmonize the energy intensities / efficiencies to match the IEA energy balance from 2005 onward.
#' Thereby, we use the energy service demand trajectories from GCAM and the energy balances from IEA.
#' We scale the energy intensities of the GCAM energy services to match the IEA energy balances for the years
#' 2005, 2010, 2015, 2020 and after that we use scaling factor from 2020 until 2100
#'
#' @param demKm demand in million (p,t)km/yr from mrEDGET
#' @param int intensity MJ / (p,t)km (load factor is already applied) from mrEDGET
#' @param IEA IEA balances from mrCommons
#' @importFrom rmndt magpie2dt
#'
#'


toolIEAharmonization <- function(demKm, int, IEA, l_f) {
    te <- isbunk <- flow <- value <- `.` <- region <- EJ_Mpkm <- conv_pkm_MJ <- subsector_L1 <- 
    subsector_L3 <- technology <- EJ_Mpkm.mean <- vehicle_type <- sector <- EJ_Mpkm_ave_adjusted <- 
    Mpkm_tot <- factor_intensity <- EJ_Mpkm_ave <-EJ_Mpkm_final <- NULL
    
    #CHoose the right IEA data and convert to logit structure 
    #Select final energy demand for transportantion (fe) from diesel (die), petrol (pet), gases (gat), electricity (elt)
    IEA <- IEA[, , c("fedie", "fepet", "fegat", "feelt"), pmatch = TRUE]

    #se: seel, seliqfos, segafos, segabio, fe: fedie, fepet, fegat, feelt, mod: Additives, multiple Gases and ELECTR, flow: flow technology, such as PIPELINE or Rail
    #https://github.com/remindmodel/remind/blob/develop/core/sets.gms
    #TRNONSPE: Transport not elsewhere specified, NETRANS: Non-energy use in transport

    IEA_dt <- magpie2dt(IEA, datacols=c("se", "fe", "te", "mod", "flow"), regioncol="region", yearcol="year")

    #delete fedie.dot which is a Diesel Oil Turbine and not relevent for transport technologies
    IEA_dt <- IEA_dt[te != "dot"]

    #We have IEA values for 2005, 2010, 2015 and 2020
    IEA_dt = IEA_dt[year %in% c(2005, 2010, 2015, 2020)]
    
    #map IEA data to logit structure (Level0)
    #Differentiate between bunkers and other tranport modes, now described as "short-medium"
    IEA_dt[, isbunk := ifelse(grepl("BUNK", flow), flow, "short-medium")]
    IEA_dt[, c("se", "fe", "mod", "flow") := NULL]

    ## sum fossil liquids and biofuel to tdlit, and biogas and natural gas to tdgat
    IEA_dt[te %in% c("tdfospet", "tdfosdie", "tdbiopet", "tdbiodie"), te := "tdlit"]
    IEA_dt[te %in% c("tdfosgat", "tdbiogat"), te := "tdgat"]
    IEA_dt[, value := sum(value), by=.(region, year, te, isbunk)]

    IEA_dt <- unique(IEA_dt)

    ## load (p, v)km and intensity from GCAM database (mrEDGET)
    CONV_MJ_EJ <- 1e-12 #MJ->EJ
    CONV_millionkm_km <- 1e6 #million (p, v)km -> (p, v)km
    vehicle_intensity <- copy(int) #load factor is already applied
    vehicle_intensity[, EJ_Mpkm := conv_pkm_MJ * CONV_millionkm_km * CONV_MJ_EJ][,conv_pkm_MJ:=NULL] #EJ/(p, v)km this is the metric we want to return at the end this function

    ## Energy Service Demand is given in million (p, t)km / yr
    tech_output <- copy(demKm) 
    tech_output <- merge(tech_output, vehicle_intensity, all.x = TRUE,
                         by = intersect(colnames(tech_output), colnames(vehicle_intensity))) #merge ES demand with intensity
    tech_output <- tech_output[!subsector_L3 %in% c("Cycle", "Walk"), ]  #no non-motorized **this is not needed anymore, as we have hopefully already filtered out non-motorized in mrEDGET

    ## use only 2005, 2010, 2015 and 2020
    tech_output <- tech_output[year %in% c(2005, 2010, 2015, 2020)]

    #for merge to work, we need to rename the columns
    setkey(tech_output, "technology") 

    ## apply the IEA categories in (EJ/yr) (tdelt, tdlit, tdgat) to the GCAM categories (BEV, Electric, NG, ...) from tech_output (Energy Service Demand: (p, t)km/yr)
    elts <- c("Electric", "BEV") #Electric and tdel belong together
    tech_output[technology %in% elts, te := "tdelt"]
    gats <- "NG" #Natural Gas and tdgat belong together
    tech_output[technology %in% gats, te := "tdgat"]
    tech_output[is.na(te), te := "tdlit"] ## all other technologies belong to tdlit (all others are liquids, including some coal an Hybrid Electric)
    
    #** to be checked: why do we have to do this?
    tech_output <- tech_output[tech_output > 0]
    dups <- duplicated(tech_output, by=c("region", "technology", "vehicle_type"))
    if(any(dups)){
        warning("Duplicated techs found in supplied demand.")
        print(tech_output[dups])
        tech_output <- unique(tech_output, by=c("region", "technology", "vehicle_type"))
    }
    tech_output[, EJ_Mpkm.mean := mean(EJ_Mpkm, na.rm = T), by=.(year, technology, vehicle_type)
                ][is.na(EJ_Mpkm), EJ_Mpkm := EJ_Mpkm.mean
                  ][,EJ_Mpkm.mean := NULL]## if there is output but no intensity, we have to apply avg. intensity

    ## apply logit structure (Level0) to GCAM ES Demand data
    tech_output[, isbunk := ifelse(sector == "trn_aviation_intl", "AVBUNK", NA)]
    tech_output[, isbunk := ifelse(sector == "trn_shipping_intl", "MARBUNK", isbunk)]
    tech_output[, isbunk := ifelse(is.na(isbunk), "short-medium", isbunk)]

    ## aggregate to IEA demand data (EJ/yr)    
    tech_output_aggr <- tech_output[, .(Mpkm_tot = sum(tech_output),
                                        EJ_Mpkm_ave = sum(tech_output/sum(tech_output) * EJ_Mpkm)),
                                        by = c("region", "year", "te", "isbunk")]

    ## merge with IEA
    tech_output_iea <- merge(IEA_dt, tech_output_aggr, all.y = TRUE, by = c("year","region", "te", "isbunk"))

    ##** to be checked: why do we have to do this?
    ## inconsistencies such as IEA stating 0 feelt in AFG and GCAM saying that the
    ## total pkm there are1.791546e+07 are solved in the hard way, deleting the demand
    ## for the country that in theory is there according to GCAM
    tech_output_iea <- tech_output_iea[value > 0]

    ## calculate intensity factor
    tech_output_iea[, EJ_Mpkm_ave_adjusted := value/Mpkm_tot]

    ## calculate the ratio between the new and the old energy intensity to see if we have to adjust the GCAM data to hit the IEA data
    tech_output_iea[, factor_intensity := EJ_Mpkm_ave_adjusted/EJ_Mpkm_ave][, year := NULL]

    ## redistribute the energy intensity to the broader category they belong to, in the energy intensity dt
    tech_output <- merge(tech_output, tech_output_iea, all = FALSE,
                             by = c("region", "te", "isbunk"))

    # ** to be checked: why do we have to do this? And if it is wrong to drop the year
    ## remove all the columns that we don't need, including years (so to have it for all years)
    tech_output=tech_output[, c("region","technology","vehicle_type", "factor_intensity")]

    ## harmonize data
    merged_intensity <- tech_output[vehicle_intensity, on=c("region", "technology", "vehicle_type")]
    merged_intensity[, EJ_Mpkm_final := EJ_Mpkm * factor_intensity]

    ## if there is no harmonization data, lets use the existing one
    merged_intensity[is.na(EJ_Mpkm_final), EJ_Mpkm_final := EJ_Mpkm]

    ## delete columns not useful anymore c("EJ_Mpkm", "factor_intensity")
    merged_intensity[, c("EJ_Mpkm", "factor_intensity") := NULL]

    ##** to be checked: why do we have to do this?
    ## cut energy intensity to 3 for NG LDVs
    merged_intensity[subsector_L1 == "trn_pass_road_LDV_4W" & technology == "NG", EJ_Mpkm_final := min(EJ_Mpkm_final, 3.2e-06), by = c("region", "year", "vehicle_type")]

    return (merged_intensity)

}
