/***********************************************************************
 *
 *variable flow velocity, developed by Kerstin Verzano (former Schulze) and
 *Frank Voss as well as restructured routing scheme. This is one of the main
 *changes of WaterGAP 2.2 refvers compared to 2.1 versions. This is described
 *in the paper Verzano, K., I. Bärlund, M. Flörke, B. Lehner, E. Kynast, F.
 *Voß, and J. Alcamo (2012), Modeling variable river flow velocity on
 *continental scale: Current situation and climate change impacts in Europe,
 *Journal of Hydrology, 424-425, 238-251, doi:10.1016/j.jhydrol.2012.01.005.
 *[online] Available from:
 *http://linkinghub.elsevier.com/retrieve/pii/S0022169412000248 as well as in
 *the thesis Verzano, K. (2009), Climate Change Impacts on Flood Related
 *Hydrological Processes: Further Development and Application of a Global Scale
 *Hydrological Model, 166 pp., University of Kassel. Within that framework, the
 *bankfull flow concept was integrated by Cristof Schneider and described in the
 *paper Schneider, C., M. Flörke, S. Eisner, and F. Voss (2011), Large scale
 *modelling of bankfull flow: An example for Europe, Journal of Hydrology,
 *408(3-4), 235-245, doi:10.1016/j.jhydrol.2011.08.004. [online] Available
 *from: http://linkinghub.elsevier.com/retrieve/pii/S0022169411005300
 *
 *water use is now divided into water uptake from groundwater and surface
 *water. This was coded by Heike Hoffmann-Dobrev, a tool to calculate the net
 *uses of groundwater and surface water (GWSWUSE) was developed from Felix
 *Portmann and Heike Hoffman-Dobrev. The concept is described in the paper
 *Döll, P., Hoffmann-Dobrev, H., Portmann, F.T., Siebert, S., Eicker, A.,
 *Rodell, M., Strassberg, G., Scanlon, B.R. (2012): Impact of water withdrawals
 *from groundwater and surface water on continental water storage variations,
 *Journal of Geodynamics (in press).
 *
 *calculation of groundwater storages was moved from daily.cpp to routing.cpp
 *by Heike Hoffmann-Dobrev; storage otuputs were restrctured and enhanced by
 *Heike Hoffmann-Dobrev, Linda Adam, Stephanie Eisner and Hannes Müller
 *Schmied); several new outputs are integrated (by Ellen Kyast (former
 *Teichert) Frank Voss, Stephanie Eisner (e.g. max and min monhtly river
 *availbility), Hannes Müller Schmied (e.g. restructuring .31 and .365 output)
 *
 *new reservoir algorithm from Kristina Fiedler, described at the paper Döll,
 *P., K. Fiedler, and J. Zhang (2009), Global-scale analysis of river flow
 *alterations due to water withdrawals and reservoirs, Hydrology and Earth
 *System Sciences, 13(12), 2413-2432, doi:10.5194/hess-13-2413-2009. [online]
 *Available from: http://www.hydrol-earth-syst-sci.net/13/2413/2009/ as well as
 *work on that by Ellen Kynast (former Teichert)
 *
 *handling of .365 climate input by Heike Hoffmann-Dobrev
 *
 *litle inconsistencies due to numeric problems when reading in float values of
 *surface water body input files were catched throug a adapted storage number
 *of lake/wetland fraction error message (Linda Adam, Hannes Müller Schmied)
 *
 *all mean input (eg. mean temperature) were changed from 1961-1990 to
 *1971-2000 because this period is the focus period of calibration (Hannes
 *Müller Schmied)
 *
 *the allocated water use of surface water as preprocessing step for
 *calibration (where water use from neighbour cell is not possible at catchment
 *borders) is now G_ACTUAL_NAs (Stephanie Eisner, Hannes Müller Schmied)
 *
 *important: it was decided to do the calibration with station correction but
 *normal runs without station correction e.g. when the model time period
 *differs from calibration period. This can be switched on and of with the new
 *option 24 (Hannes Müller Schmied).
 *
 * see former changes at file routing.cpp.versioninfos.txt
 *
 ***********************************************************************/
//#define NDEBUG

// Felix Portmann (FP) 2015 & 2016 several changes
//  Saved in UTF-8 format
//  Inserted "// end"-statements
//      - for classes
//      - for/if/else statements (mostly class/method "routing")
//  Replacing all "++n" by "n++" (filling / plain treatment procedures)
//  Use of geo.G_contcell[n] (geo.cpp)
//      - With query (geo.G_landfreq[n] + geo.G_fwaterfreq[n]) < limit_contcell_pct)
//      instead of query ((geo.G_landfreq[n] + geo.G_fwaterfreq[n]) == 0) (or "== 0."
//      - For correctly identifying contintal cells vs. ocean cells of Caspian Sea
//  Moving some brackets "}" to correct position in class/method "routing"
//      for G_tobeCalculated[n]
//  Inserted checks functions for testing "less than zero" and "Not a Number" (NaN)
//  Using cell-specific Calibration Parameters
// FP20161010N001 Replacement: geo.G_contfreq[n] instead (geo.G_landfreq[n] + geo.G_fwaterfreq[n])
// FP20161018N002 Reservoir operation start years

#include <cstdio>
#include <cmath>
#include <cassert>
#include "timestring.h"
#include "def.h"

#include "common.h"
#include "globals.h"
#include "calib_param.h"

using namespace std;



// Conversion factors  // FP
#include "conversion.h"

#include "calcWaterTemp.h"

calcWaterTempClass cwt;


const int numberOfDaysInMonth[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

template<class T>
inline int round(const T value) {
    return (int) floor(value + 0.5);
}
// Definitions for global checks (all grid cells):
// use of extern variables (chk_year, chk_month, chk_day_in_month), classes/functions, exceptions
// Felix Portmann 2015
#ifdef _WATERGAP_CHECKS_GLOBAL_H
#include "routing_checks_global.h"
#endif
routingClass::routingClass() {
    defaultRiverVelocity = 86.4; // [km/d] = 1 m/s
    defaultTimeStepsPerDay = 1;
}

void routingClass::setParameters() {
    lakeOutflowExp = 1.5;
    wetlOutflowExp = 2.5;
    evapoReductionExp = 3.32193;
    evapoReductionExpReservoir = 2.81383;
}

void routingClass::init(const short nBasins, ConfigFile configFile, WghmStateFile &wghmState, AdditionalOutputInputFile &additionalOutIn) {
    // - reads the files with routing parameters and
    //   grids with lake and wetland locations
    // - calculates land area fraction
    // - reads files with flow direction, etc.
    // - dynamic memory allocation

    // requirements:
    // - G_sbasin has been read in 'main'
    char filename[250];
    // GFZ:GLOLAK begin - initialise and read variables of biglakefile (from Kristina)
    if (0 != ng_glolakcells) {
        for (int n = 0; n < ng_glolakcells; n++) {
            big_lakes_area_frac[n] = 0.0;
            big_lakes_cells(0,n) = -99;
            big_lakes_cells(1,n) = -99;
        }
        sprintf(filename, "%s/wg_big_lakes_area_frac.txt", options.input_dir.c_str());
        ifstream biglakefile(filename);
        if (biglakefile) {
            int n = 0;
            // float data; // FP20161018N002 Reservoir operation start years: uncomment unused variable
            while (biglakefile && (n < ng_glolakcells)) {
                biglakefile >> big_lakes_cells(0,n);
                biglakefile >> big_lakes_cells(1,n);
                biglakefile >> big_lakes_area_frac[n];
                n++;
            }

            biglakefile.close();
            cout << "Allocate GloLakeStorage for " << ng_glolakcells << " lakecells in wg_big_lakes_area_frac.txt."
                 << endl;
        } else {
            cout << "ERROR: wg_big_lakes_area_frac.txt not read." << endl;
            cout << "Using default values (no allocation of lakestorage on all lakecells)." << endl;
        }
    }
    // GFZ:GLOLAK end

    setParameters();

    if (options.calc_wtemp == 1 || options.calc_wtemp == 2) {
        // Initialize Temporary Storage to hold Temperatures and Inflowvolumes from previous day to starttemp. 15 degC
        G_RiverWTempPreStep.fill(15.);
        G_RiverAvailPreStep.fill(0.);
        G_locLakeTempPrevStep.fill(15.);
        G_locWetlTempPrevStep.fill(15.);
        G_gloLakeTempPrevStep.fill(15.);
        G_reservoirTempPrevStep.fill(15.);
        G_gloWetlandTempPrevStep.fill(15.);
        SumWaterTempRiver.fill(0.);
        NumWaterTempRiver.fill(0.);
        SumWaterTempLocLake.fill(0.);
        NumWaterTempLocLake.fill(0.);
        SumWaterTempLocWetl.fill(0.);
        NumWaterTempLocWetl.fill(0.);
        SumWaterTempGloLake.fill(0.);
        NumWaterTempGloLake.fill(0.);
        SumWaterTempReservoir.fill(0.);
        NumWaterTempReservoir.fill(0.);
        SumWaterTempGloWetl.fill(0.);
        NumWaterTempGloWetl.fill(0.);
        G_iceThicknesslocLake.fill(0.);
        G_iceThicknesslocWetland.fill(0.);
        G_iceThicknessGloLake.fill(0.);
        G_iceThicknessReservoir.fill(0.);
        G_iceThicknessGloWetland.fill(0.);
        G_iceThicknessRiver.fill(0.);
        // read in all inflow cells for every cell for temperature computation
        G_inflow_cells.read(options.routing_dir + "/G_INFLC.9.UNF4");
    }

    // read parameters from ROUTING.DAT
    double v = defaultRiverVelocity;
    short s = defaultTimeStepsPerDay;

    try {
        // sprintf(filename, "ROUTING.DAT");
        // AEMS: new code:
        sprintf(filename, configFile.routingoptionsfile.c_str());
        ifstream routingOptionFile(filename);

        if (!routingOptionFile) {
            throw routingOptionFile.rdstate();
        } else {
            char string[250];

            // read commentlines indicated by #
            while (!routingOptionFile.eof() && routingOptionFile.peek() == '#') {
                routingOptionFile.getline(string, sizeof(string));
            }

            routingOptionFile >> v; // river velocity
            if (!routingOptionFile.good())
                throw routingOptionFile.rdstate();
            else
                defaultRiverVelocity = v;

            routingOptionFile >> s;
            if (!routingOptionFile.good())
                throw routingOptionFile.rdstate();
            else
                defaultTimeStepsPerDay = s;

        }
    } catch (ios::iostate &status) {
        cout << filename << " not read. Using default values.\n";
        // cout << " Status: " << status << endl;
    }

    // defaults for velocity and timesteps
    riverVelocity = defaultRiverVelocity; // [km/day]
    timeStepsPerDay = defaultTimeStepsPerDay;

    cout << "Routing options: riverVelocity = " << riverVelocity << " km/day; "
         << timeStepsPerDay << " time steps per day." << endl;

    // read lakes, wetlands and reservoirs
    // and remove nodata values

    // fractions of grid cell area
    G_glo_lake.read(options.input_dir + "/G_GLOLAK.UNF0");
    G_loc_lake.read(options.input_dir + "/G_LOCLAK.UNF0");
    G_glo_wetland.read(options.input_dir + "/G_GLOWET.UNF0");
    G_loc_wetland.read(options.input_dir + "/G_LOCWET.UNF0");

    // FP20161018N002 Reservoir operation start years
    // read all areas in km2 into arrays of data type double
    G_lake_area.read(options.input_dir + "/G_LAKAREA.UNF0");

    // reservoir AND regulated lake areas when all reservoirs/reg.lakes are in operation
    G_reservoir_area_full.read(options.input_dir + "/G_RESAREA.UNF0");

    // fraction of grid cell area
    G_reg_lake.read(options.input_dir + "/G_REGLAKE.UNF0");

    //read regulated lake status (1 == lake is regulated lake), area included in G_RESAREA.UNF0
    G_reg_lake_status.read(options.input_dir + "/G_REG_LAKE.UNF1");

    // and remove nodata values
    for (int n = 0; n < ng; n++) {

        if (G_glo_lake[n] < 0.) G_glo_lake[n] = 0.;
        if (G_glo_wetland[n] < 0.) G_glo_wetland[n] = 0.;
        if (G_loc_lake[n] < 0.) G_loc_lake[n] = 0.;
        if (G_loc_wetland[n] < 0.) G_loc_wetland[n] = 0.;
        if (G_lake_area[n] < 0.) G_lake_area[n] = 0.;
        // FP20161018N002 Reservoir operation start years
        if (G_reservoir_area_full[n] < 0.) G_reservoir_area_full[n] = 0.;
        if (G_reg_lake[n] < 0.) G_reg_lake[n] = 0.;

        // initial values for fractions
        // (1) water areas (global lakes + local lakes + reservoirs (including regulated lakes)
        // (2) land areas (including wetlands)
        G_fwaterfreq[n] = 0.;
        G_landfreq[n] = 0.;
    }

    // GLWDunit pattern with index to outflow grid cell Image22-Number (n = Image22nr - 1)
    G_outflow_cell_assignment.read(options.input_dir + "/G_OUTFLOW_CELL_ASSIGNMENT.UNF4");


    // read cell number of downstream cell;
    // added for new use_allocation algorithm, M.Hunger 2/2006,
    // and also for the routing algorithm (F.Voss 4/2008)
    G_downstreamCell.read(options.routing_dir + "/G_OUTFLC.UNF4");

    //define status at model start / first occurence, e.g. when reading for annual reservoir fraction arrays
    // statusStarted_updateGloResPrevYear = 0; //obsolete as covered in initFractionStatus

    // Status of algorithm depends on
    // whether naturalized runs are performed or whether reservoirs (incl. regulated lakes) are to be used
    if (options.antNatOpt == 0) { // anthropogenic runs: consider reservoirs and regulated lakes as reservoirs
        cout << "routing.init: if options.antNatOpt == 0, anthropogenic run:"
             << " consider reservoirs and regulated lakes as reservoirs" << endl;

        //select reservoir algorithm yes/no

        // Adding to lake area depending on
        // (1) reservoir algorithm option
        // (2) status of the reservoir / regulated lake: type, mean outflow
        // Calculate here, otherwise, areas would be added each year.

        // Anthropogenic run
        // read files which do not change in time
        // Use reservoir algorithm for reservoirs/reg.lakes


        if (options.resOpt == 1) {
            cout << "routing.init: if options.resOpt == 1, initialize new reservoir algorithm" << endl;

            //reading in and adding local reservoirs to loc lakes(reservoirs below the
            // threshold of 0.5 km3 reservoir volume)
            G_loc_res.read(options.input_dir + "/G_LOCRES.UNF0");
            G_loc_lake+=G_loc_res;

            // read reservoir_type (char)
            G_res_type.read(options.input_dir + "/G_RES_TYPE.UNF1");

            // FP20161018N002 Reservoir operation start years
            // read reservoir operation start year (int)
            G_res_start_year.read(options.input_dir + "/G_START_YEAR.UNF4");

            // read file -> start month of operational year (char)
            G_start_month.read(options.routing_dir + "/G_START_MONTH.UNF1");

            // read mean annual outflow of reservoirs (long term mean [km3/month], calculated from G_RIVER_AVAIL),
            // was G_MEAN_INFLOW before but is in fact inflow to reservoir plus P-PET of reservoir
            G_mean_outflow.read(options.input_dir + "/G_MEAN_OUTFLOW.UNF0");

            // read file -> reservoir volume from published data (GLWD1), storage capacity [km3]
            // FP20161018N002 Reservoir operation start years
            G_stor_cap_full.read(options.input_dir + "/G_STORAGE_CAPACITY.UNF0");

            // read file -> mean annual net use of surface water per cell (e.g. mean value of 1971 to 2000)
            // net use is given in m3/yr
            // FP20161018N002 Reservoir operation start years:
            // Specify name with averaging period of full year start and end in options
            // e.g. G_NUs_7100.UNF0
            G_mean_NUs.read(options.input_dir + "/G_NUs_"
                + std::to_string(options.resNUsMeanYearFirst) + "_"
                + std::to_string(options.resNUsMeanYearLast) + ".UNF0");

            G_alloc_coeff.read(options.routing_dir + "/G_ALLOC_COEFF."
                + std::to_string(reservoir_dsc) + ".UNF0");

                // FP20161018N002 Reservoir operation start years:
                // Remark: Calculate this only once (1x)
                // Change in code: G_reservoir_area_full instead G_reservoir_area
                // Use all reservoirs and regulated lakes
            // calculate MEAN demand of downstream area
                for (int n = 0; n < ng; n++) {
                    G_mean_demand[n] = 0; //init
                    G_mean_demand[n] += G_mean_NUs[n]; //water demand of cell itself
                    G_glores_prevyear[n] = 0.;
                    G_glores_change[n] = 0.;
                    short i = 0;
                    int downstreamCell = G_downstreamCell[n];

                // max reservoir_dsc cells downstream or to the basin outlet or area up to next reservoir
                while (i < reservoir_dsc && downstreamCell > 0
                       && G_reservoir_area_full[downstreamCell] <= 0) {
                    //add demand of downstream cells, m3/yr
                    G_mean_demand[n] += G_mean_NUs[downstreamCell - 1] * G_alloc_coeff(n,i++);

                    // next downstream cell
                    downstreamCell = G_downstreamCell[downstreamCell - 1];
                }
            }

            // Remark: Calculate this only once (1x)
            // Change in code: G_stor_cap_full instead G_stor_cap
            // Calculate G_mean_outflow[n] & G_mean_demand[n]
            for (int n = 0; n < ng; n++) {
                if (G_stor_cap_full[n] < 0.)
                    G_stor_cap_full[n] = 0.;
                //AEMS
                if (additionalOutIn.additionalfilestatus == 0) {
                    K_release[n] = 0.1; //HMS 2015-03-02 Experience from ISIMIP, leads to a reduction of initial outflow before start month is reached
                }
                else { // when it should be read from file
                    K_release[n] = additionalOutIn.additionalOutputInput(n,5);
                }

                // FP20161018N002 Reservoir operation start years
                // No change in code, only 2 remarks:
                // (1) Calculate this only once (1x)
                // (2) Otherwise values are increasing until nan!
                // G_mean_outflow is long term mean annual value ([km3/month]); required number is mean annual value
                G_mean_outflow[n] = G_mean_outflow[n] * 12. * 1000000000. /
                                    31536000.; // 31536000=(365 * 24 * 60 * 60); //km3/month -> km3/yr -> m3/s
                G_mean_demand[n] = G_mean_demand[n] /
                                   31536000.; // 31536000=(365 * 24 * 60 * 60); //km3/yr -> m3/s // G_mean_demand is in m�/yr as G_NUs_6190 is in km�/yr
            }

        }
            // endif (options.resOpt == 1)
            // FP20161018N002 Reservoir operation start years
            // Anthropogenic run
            // Reservoirs are not considered, regulated lakes are handled as global lakes
        else if (options.resOpt == 0) {
            cout
                           << "routing.init: else if options.resOpt == 0, regulated lakes are handled as global lakes"
                           << endl;


            for (int n = 0; n < ng; n++) {
                // treat area
                // for regulated lakes (in G_reservoir_area_full[n]), add area to lake area
                if (G_reg_lake_status[n] == 1) {
                    if (G_reservoir_area_full[n] > 0.) {
                        G_lake_area[n] += G_reservoir_area_full[n]; //add reservoir area to lake area in case of regulated lake
                        G_reservoir_area_full[n] = 0.; //set area to zero (inhibits further accounting elsewhere)
                    }
                }
                // treat land cover fractions
                G_glo_lake[n] += G_reg_lake[n]; // add regulated lake fraction to glo_lake fraction
                G_reg_lake[n] = 0.;  //set fraction to zero (inhibits further accounting elsewhere)
                // set reservoir fraction to zero
                G_glo_res[n] = 0.;
                G_glores_prevyear[n] = 0.;
            }
        }
    }
    if ((options.antNatOpt == 1) && (options.resOpt ==1)){
        cout << "Since options.antNatOpt ==1, options.resOpt cannot be 1: using options.resOpt ==1" << endl;
        options.resOpt = 0;

    }
    if ((options.antNatOpt == 1) && (options.resOpt == 0)) { //naturalized run
        cout << "routing.init: if options.antNatOpt == 1, naturalized run: using regulated lakes as lakes" << endl;

        for (int n = 0; n < ng; n++) {
            // treat land cover fractions
            G_glo_res[n] = 0.; //no reservoir fractions are used at all  // set fraction to zero (inhibits further accounting elsewhere)
            G_glores_prevyear[n] = G_glo_res[n]; //initialize with percent fraction of reservoirs of start or reference year
            G_glo_lake[n] += G_reg_lake[n]; //global lake fraction is increased by regulated lake fraction only (NOT reservoirs)
            G_reg_lake[n] = 0.; // set fraction to zero (inhibits further accounting elsewhere)
            // treat area
            if (G_reg_lake_status[n] == 1) {
                G_lake_area[n] += G_reservoir_area_full[n]; //add reservoir area to lake area in case of regulated lakes
                G_reservoir_area_full[n] = 0.; //set area to zero (inhibits further accounting elsewhere)  // new FP vs. ISIMIP of HMS
            }
        }

    }

    // Reservoir operation start years
    // Counting number of reservoirs and regulated lakes
    long count = 0;
    for (int n = 0; n < ng; n++) {
        if (G_reservoir_area_full[n] > 0.) count++;
    }
    cout << "routing.init: nr. of reservoirs/regulated lakes detected (area > 0)  : " << count << endl;


    // set correction factor for upstream locations to 1 (= no correction) use statcorrOpt == 1 during calibration runs (here, station correction should occur). During normal runs, this line should be normally 0 (runs without station correction).
    if (1 == options.statcorrOpt) {
        G_statCorrFact.fill(1.0);
    }

    // before calculating landareaFrac, initial riverFrac needs to be calculated, based on bankfull flow and default flow velocity.

    // read stretch length for each river segment
    // (added for new routing-algorithm, F.Voss 4/2008)
    G_riverLength.read(options.routing_dir + "/G_RIVER_LENGTH.UNF0");
    G_riverLengthCheck.read(options.routing_dir + "/G_RIVER_LENGTH.UNF0");

    // Introduced in WG22b: Reduce river length in coastal cells in proportion to continental area fraction:
    for (int n = 0; n < ng; n++) {
        G_riverLength[n] *= geo.G_contfreq[n] / 100.;
    }

    // velocity -------------------------------------------------------

    G_RiverSlope.read(options.routing_dir + "/G_RIVERSLOPE.UNF0");

    G_Roughness.read(options.input_dir + "/G_ROUGHNESS.UNF0");

    G_bankfull_flow.read(options.input_dir + "/G_BANKFULL.UNF0");  //[m/sec]

    for (int n = 0; n < ng; n++) {
        // to avoid negative bottom width
        if (G_bankfull_flow[n] < 0.05)
            G_bankfull_flow[n] = 0.05;

        G_RiverWidth_bf[n] = 2.71 * pow(G_bankfull_flow[n], 0.557); //[m]
        G_RiverDepth_bf[n] = 0.349 * pow(G_bankfull_flow[n], 0.341); //[m]
        // calculate river bottom width asssuming a trapezoidal channes with 2/1 run to rise ratio; [m]
        G_riverBottomWidth[n] = G_RiverWidth_bf[n] - 2.0 * 2.0 * G_RiverDepth_bf[n];
        //CR: G_riverStorageMax[n] is required to compute riverEvapoReductionFactor introduced in WG22b.
        G_riverStorageMax[n] = G_riverLength[n] * 0.5 * G_RiverDepth_bf[n] / 1000. *
                               (G_riverBottomWidth[n] / 1000. + G_RiverWidth_bf[n] / 1000.); // [km3]
        //set initial value for G_riverAreaFrac to 0 (later updated)
        G_riverAreaFrac[n] = 0.;

    }

    // set initialization of G_landAreaFracNextTimestep to -99 for the first model day
    // to calculate new special value at the end of the first model day
    // Set to zero to prevent undesired results
    G_landAreaFracNextTimestep.fill(0);
    G_landAreaFrac.fill(0);
    G_landAreaFracPrevTimestep.fill(0);

    // set G_GloLakePrecip and G_ReservoirPrecip to 0
    G_GloLakePrecip.fill(0);
    G_ReservoirPrecip.fill(0);
    G_dailyRiverPrecip.fill(0);

    // read image numbers of neighbouring cells of each cell
    // (added for new use_allocation algorithm, M.Hunger 2/2006)
    if (options.use_alloc > 0) {
        G_neighbourCell.read(options.routing_dir + "/G_NEIGHBOUR_CELLS.8.UNF4");
    }

    // read ordered routing number of each cell
    // and assign array position for ordered routing scheme
    // (added for new routing-algorithm, F.Voss 4/2008)
    Grid<int> G_routOrder;
    G_routOrder.read(options.routing_dir + "/G_ROUT_ORDER.UNF4");

#pragma omp parallel for \
        shared(G_routingCell, G_routOrder)
    for (int i = 0; i < ng; ++i)
    {
        G_routingCell[G_routOrder[i] - 1] = i + 1;
    }

    // read frgi to adapt NAg if necessary
    if (options.subtract_use > 0) {
        G_fractreturngw_irrig.read(options.input_dir + "/G_FRACTRETURNGW_IRRIG.UNF0");
    }

    G_LDD.read(options.routing_dir + "/G_LDD_2.UNF1");

    // number of superbasins varies.
    // therefore size of the arrays is allocated dynamically.
    nSpecBasins = nBasins; // number of specified basins

    cout << " : ##### routing.init:  nSpecBasins= " << nSpecBasins << endl;

    dailyRiverDischarge.initialize(nSpecBasins);
    // initialize daily storages
    if (options.day_store == 1) {
        dailyLocLakeStorage.initialize(nSpecBasins);
        dailyLocWetlStorage.initialize(nSpecBasins);
        dailyGloLakeStorage.initialize(nSpecBasins);
        dailyGloWetlStorage.initialize(nSpecBasins);
        dailyRiverVelocity.initialize(nSpecBasins);
        if (options.resOpt == 1)
            dailyResStorage.initialize(nSpecBasins);
    }
    // endif options.day_store

    // initialize monthly storages
    if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {
        G_monthlyConsistentPrecip.initialize();
        G_monthlyCellRunoff.initialize();
        G_monthlyCellSurfaceRunoff.initialize();
        G_monthlyRiverAvail.initialize();
        G_monthlyPotCellRunoff.initialize();
        G_monthlyPotCellRunoffDeficit.initialize();
        G_monthlyGwrunSurfrun.initialize();

        if (options.outRiverInUpstream)
            G_monthlyRiverInUpstream.initialize(); // total inflow from upstream

        // cell AET calculation (since 2.1f)
        G_monthlyCellAET.initialize();
        G_monthlyCellAETWCa.initialize();
        G_monthlyOpenWaterEvap.initialize();

        G_monthlyVelocity.initialize(); // velocity grid
        G_monthlySurfStor.initialize(); // surface water storage
        G_monthlySurfStor_mm.initialize(); // surface water storage [mm]

        G_monthlyMinRiverAvail.initialize();
        G_monthlyMaxRiverAvail.initialize();

        if(options.outWaterTempMonthly) G_monthlyMeanWaterTempRiver.initialize();
        if(options.outWaterTempMonthlyAllSWB){
            G_monthlyMeanWaterTempLocLake.initialize();
            G_monthlyMeanWaterTempLocWetl.initialize();
            G_monthlyMeanWaterTempGloLake.initialize();
            G_monthlyMeanWaterTempReservoir.initialize();
            G_monthlyMeanWaterTempGloWetl.initialize();
        }

        G_monthlyGwrSwb.initialize();
        G_monthlyFswb.initialize();
        G_monthlyLandAreaFrac.initialize();
        G_monthlyRiverAreaFrac.initialize();
        G_monthlyRiverPET.initialize();
        G_monthlyGwRunoff.initialize();
        G_monthlyLocWetlExtent.initialize();
        G_monthlyGloWetlExtent.initialize();

        //define arrays for storage output
        if (options.resOpt == 1) {
            G_monthlyResStorage.initialize();
            G_monthlyResStorageRatio.initialize();
            G_monthlyResInflow.initialize();
            G_monthlyResOutflow.initialize();
        }

        if (options.glacierOpt == 1) {
            G_monthlyGlacierStorage.initialize();
            G_monthlyGlacierArea.initialize();
            G_monthlyGlacierAreaFrac.initialize();
            G_monthlyGlacierRunoff.initialize();
            G_monthlyGlacierPrecip.initialize();
        }

        G_monthlyRiverStorage.initialize();
        G_monthlyLocLakeStorage.initialize();
        G_monthlyLocWetlStorage.initialize();
        G_monthlyGloLakeStorage.initialize();
        G_monthlyGloWetlStorage.initialize();
        G_monthlySatisfiedUse.initialize();
        G_monthlyGwStorage.initialize();
        G_monthlyActualUse.initialize();
        G_monthlyNUg.initialize();
        G_monthlyAllocatedUse.initialize();
    }

    if ((3 == options.grid_store) || (4 == options.grid_store)) {
        if (options.outPrecDaily) G_daily31ConsistentPrecip.initialize();
        if (options.outCellRunoffDaily) G_daily31CellRunoff.initialize();
        if (options.outGWRunoffDaily) G_daily31GwRunoff.initialize();
        if (options.outCellSurfaceDaily) G_daily31CellSurfaceRunoff.initialize();
        if (options.outRiverAvailDaily) G_daily31RiverAvail.initialize();
        if (options.outRiverVeloDaily) G_daily31Velocity.initialize();
        if (options.outWaterTempDaily) G_daily31WaterTemp.initialize();
        if (options.outWaterTempDailyAllSWB) {
            G_daily31locLakeTemp.initialize();
            G_daily31locWetlTemp.initialize();
            G_daily31gloLakeTemp.initialize();
            G_daily31reservoirTemp.initialize();
            G_daily31gloWetlandTemp.initialize();
        }
        if (options.outCellAETDaily) G_daily31CellAET.initialize();
        if (options.outSurfStorDaily) G_daily31SurfStor.initialize();
        if (options.outSingleStoragesDaily) G_daily31LocLakeStor.initialize();
        if (options.outSingleStoragesDaily) G_daily31LocWetlStor.initialize();
        if (options.outSingleStoragesDaily) G_daily31GloLakeStor.initialize();
        if (options.outSingleStoragesDaily) G_daily31GloWetlStor.initialize();
        if (options.outSingleStoragesDaily) G_daily31RiverStor.initialize();
        if (options.outSingleStoragesDaily) G_daily31ResStor.initialize();
        if (options.outSingleStoragesDaily || options.outGWStorageDaily) G_daily31GwStor.initialize();
        if (options.outTotalWaterInStoragesDaily_km3 ||
            options.outTotalWaterInStoragesDaily_mm)
            G_daily31TotalWaterInStorages_km3.initialize();
        if (options.outTotalWaterInStoragesDaily_mm) G_daily31TotalWaterInStorages_mm.initialize();
        if (options.outGwrSwbDaily) G_daily31GwrSwb.initialize();
        if (options.outFswbDaily) G_daily31Fswb.initialize();
        if (options.outLandAreaFracDaily) G_daily31LandAreaFrac.initialize();
        if (options.outGwrunSurfrunDaily) G_daily31GwrunSurfrun.initialize();
        if (options.outCellAETWCaDaily) G_daily31CellAETWCa.initialize();

    }

    if (options.outTotalWaterInStoragesStartEndDaily_km3)
        G_startendTotalWaterInStorages_km3.initialize();

    if ((5 == options.grid_store) || (6 == options.grid_store)) {

        if ((options.outPrecDaily) || ((options.scoutcPrecip) && (2 == options.day_store)))
            G_daily365ConsistentPrecip.initialize();

        if ((options.outCellRunoffDaily) ||
            (((options.scoutCellRunoffkm3) || (options.scoutCellRunoffmm)) && (2 == options.day_store)))
            G_daily365CellRunoff.initialize();
        if ((options.outCellRunoffDaily) ||
            ((options.scoutCellRunoffmm) && (2 == options.day_store)))
            G_daily365CellRunoff_mm.initialize();
        if ((options.outCellSurfaceDaily) ||
            ((options.scoutCellSRunoff) && (2 == options.day_store)))
            G_daily365CellSurfaceRunoff.initialize(); // FP20160915N002
        if ((options.outRiverAvailDaily) || ((options.scoutQ) && (2 == options.day_store)))
            G_daily365RiverAvail.initialize();
        if ((options.outRiverVeloDaily) || ((options.scoutFlowVelo) && (2 == options.day_store)))
            G_daily365Velocity.initialize();
        if ((options.outGWRunoffDaily) || ((options.scoutGwRunoff) && (2 == options.day_store)))
            G_daily365GwRunoff.initialize();
        if ((options.outLandAreaFracDaily) ||
            ((options.scoutLandAreaFrac) && (2 == options.day_store)))
            G_daily365LandAreaFrac.initialize();
        if ((options.outCellAETDaily) || ((options.scoutCellAET) && (2 == options.day_store)))
            G_daily365CellAET.initialize();
        if ((options.outSurfStorDaily) || ((options.scoutSurfaceStor) && (2 == options.day_store)))
            G_daily365SurfStor.initialize();
        if ((options.outSingleStoragesDaily) ||
            ((options.scoutLocLake) && (2 == options.day_store)))
            G_daily365LocLakeStor.initialize();
        if ((options.outSingleStoragesDaily) ||
            ((options.scoutLocWet) && (2 == options.day_store)))
            G_daily365LocWetlStor.initialize();
        if ((options.outSingleStoragesDaily) ||
            ((options.scoutGloLake) && (2 == options.day_store)))
            G_daily365GloLakeStor.initialize();
        if ((options.outSingleStoragesDaily) ||
            ((options.scoutGloWet) && (2 == options.day_store)))
            G_daily365GloWetlStor.initialize();
        if ((options.outSingleStoragesDaily) || ((options.scoutRiver) && (2 == options.day_store)))
            G_daily365RiverStor.initialize();
        if ((options.outSingleStoragesDaily) || ((options.scoutReservoir) && (2 == options.day_store)))
            G_daily365ResStor.initialize();
        if ((options.outSingleStoragesDaily || options.outGWStorageDaily) ||
            ((options.scoutGwStor) && (2 == options.day_store)))
            G_daily365GwStor.initialize();
        if ((options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm) ||
            (((options.scoutTWSkm3) || (options.scoutTWSmm)) && (2 == options.day_store)))
            G_daily365TotalWaterInStorages_km3.initialize();
        if ((options.outTotalWaterInStoragesDaily_mm) ||
            ((options.scoutTWSmm) && (2 == options.day_store)))
            G_daily365TotalWaterInStorages_mm.initialize();
        if ((options.outGwrSwbDaily) || ((options.scoutGwrSwb) && (2 == options.day_store)))
            G_daily365GwrSwb.initialize();
        if ((options.outFswbDaily) || ((options.scoutFswb) && (2 == options.day_store)))
            G_daily365Fswb.initialize();
        if (options.outGwrunSurfrunDaily) G_daily365GwrunSurfrun.initialize();
        if (options.outCellAETWCaDaily) G_daily365CellAETWCa.initialize();
        if (options.outGlacierStorageDaily) G_daily365GlacierStorage.initialize();
        if (options.outGlacierStorageDaily_mm) G_daily365GlacierStorage_mm.initialize();
    }

    //setStoragesToZero(); // HMS old code; in principle this function is also called in setStorages.
    //AEMS
    //setStorages(wghmState, additionalOutIn); //HMS will be called twice otherwise

}
// end init

void routingClass::initFractionStatus() {
    // this subroutine is necessary to reset the fractions to its initial values
    // in each first year during a calibration loop
    // otherwise, fractions are used from the previous loop and calibration fails.

    for (int n = 0; n < ng; n++) {
        statusStarted_landfreq_fwaterfreq_landAreaFrac[n] = 0;
        statusStarted_landAreaFracNextTimestep[n] = 0;

        // initialize G_fswb and related fractions
        G_fLocLake[n] = G_loc_lake[n] / 100.;
        G_fLocWet[n] = G_loc_wetland[n] / 100.;
        G_fGloLake[n] = G_glo_lake[n] / 100.; // do not adapt at all, just for output purposes
        G_fGloWet[n] = G_glo_wetland[n] / 100.;
        G_fswbInit[n] = G_fLocLake[n] + G_fLocWet[n] +
                        G_fGloWet[n]; // for output purposes and for specifying catchment area (for fractional routing), outflow cells are not considered anymore
        G_fswbLandAreaFrac[n] = G_fswbInit[n]; // for changing LandAreaFrac only
        G_fswbLandAreaFracNextTimestep[n] = G_fswbLandAreaFrac[n];  // initialize for first time step
    }
    statusStarted_updateGloResPrevYear = 0;
}

void routingClass::initFractionStatusAdditionalOI(AdditionalOutputInputFile &additionalOutIn) {
    // this subroutine is necessary to read the fractions from additionalOI array for
    // the initial time step. If not done this way, the year break makes trouble
    // for the following fractions

    for (int n = 0; n < ng; n++) {
        statusStarted_landfreq_fwaterfreq_landAreaFrac[n] = 0;
        statusStarted_landAreaFracNextTimestep[n] = 0;

        // initialize G_fswb and related fractions

        G_fLocLake[n] = additionalOutIn.additionalOutputInput(n, 45);
        G_fLocWet[n] = additionalOutIn.additionalOutputInput(n, 47);
        G_fGloLake[n] = G_glo_lake[n] / 100.; // do not adapt at all, just for output purposes
        G_fGloWet[n] = additionalOutIn.additionalOutputInput(n, 46);
        G_fswbInit[n] = additionalOutIn.additionalOutputInput(n, 14);
        G_fswbLandAreaFrac[n] = G_fswbInit[n]; // for changing LandAreaFrac only
        G_fswbLandAreaFracNextTimestep[n] = G_fswbLandAreaFrac[n];  // initialize for first time step
    }
    statusStarted_updateGloResPrevYear = 0;
}

void routingClass::setStoragesToZero() {
    // this subroutine concerns those storages
    // which have to be set to zero at the beginning of
    // a normal simulation, or at the beginning of each
    // pass of a calibration or sensitivity loop
    //
    // additional storages have to be set to zero at the beginning
    // of each simulation year

    // set storage value of lakes, wetlands and rivers to zero
    for (int n = 0; n < ng; n++) {
        G_locLakeStorage[n] = 0.;
        G_locWetlStorage[n] = 0.;
        G_gloLakeStorage[n] = 0.;
        G_gloWetlStorage[n] = 0.;
        // FP20161018N002 Reservoir operation start years
        G_landStorageChange[n] = 0.;
        //initialize release coefficient for reservoirs
        if ((options.antNatOpt == 0) && (options.resOpt == 1)) {
            G_gloResStorage[n] = 0.; // for new reservoir algorithm
        }

        if (options.glacierOpt == 1)
            G_glacierStorage[n] = 0.;

        if ((options.riverveloOpt == 0) || (options.riverveloOpt == 1)) {
            G_riverStorage[n] = 0.;
        }
        else {
            G_riverStorage[n] = G_riverStorageMax[n];
        }

        G_groundwaterStorage[n] = 0.;
        // G_totalUnsatisfiedUse is copied to G_UnsatisfiedUsePrevYear
        // at the beginning of each year.
        // Therefore it has to be set to zero before doing this
        // for the first time.
        G_totalUnsatisfiedUse[n] = 0.;        // should be initialized in any case

        // Variables introduced for Option "aggrNUsGloLakResOpt":
        // The amount of unsatisfied use in riparian cells of global lakes or reservoirs
        // can be satisfied from river storage in the next time step.
        // This is done independently from Option "delayedUseSatisfaction".
        G_unsatUseRiparian[n] = 0.;

        // New approach in 2.2b for Option "use_alloc":
        // If the remainingUse cannot be satisfied in a cell, the value is added to the remainingUse of a second cell.
        // Depending on the routing order, the water balance of the second cell is computed in the same or the next timestep.
        // This way, NAs is consistently included in the storage equations of lakes and the river.
        G_AllocatedUse[n] = 0.;
        G_daily_allocatedUseNextDay[n] = 0.;
        G_dailyAllocatedUse[n] = 0.;
        G_UnsatAllocUse[n] = 0.;
        G_daily_UnsatAllocUseNextDay[n] = 0.;
        G_AllocUseToNeigborCell[n] = 0.;
        G_SecondCell[n] = 0;    // int

    }
}
// reading dailyNUs und dailyNUg once a year in m3/month

// AEMS: new subroutine to set storages according to the start values read.
void routingClass::setStorages(WghmStateFile &wghmState, AdditionalOutputInputFile &additionalOutIn)
{
    //setStoragesToZero(); not needed anymore as it gets overwritten partly hereafter.
    double cellArea = 0.;
    for (int n = 0; n < ng; n++)
    {
        cellArea = geo.areaOfCellByArrayPos(n);
        G_locLakeStorage[n] = wghmState.cell(n).locallake(0)*((cellArea * (geo.G_contfreq[n] / 100.)) / 1000000.);
        G_locWetlStorage[n] = wghmState.cell(n).localwetland(0)*((cellArea * (geo.G_contfreq[n] / 100.)) / 1000000.);
        G_gloLakeStorage[n] = wghmState.cell(n).globallake(0)*((cellArea * (geo.G_contfreq[n] / 100.)) / 1000000.);
        G_gloWetlStorage[n] = wghmState.cell(n).globalwetland(0)*((cellArea * (geo.G_contfreq[n] / 100.)) / 1000000.);
        G_riverStorage[n]   = wghmState.cell(n).river(0)*((cellArea * (geo.G_contfreq[n] / 100.)) / 1000000.);

        G_groundwaterStorage[n]        = wghmState.cell(n).groundwater(0) * ((cellArea * (geo.G_contfreq[n] / 100.)) / 1000000.); // HMS todo: move to routing.cpp
        if (options.resOpt == 1){ //KF *** for new reservoir algorithm
            G_gloResStorage[n]  = wghmState.cell(n).reservoir(0)*((cellArea * (geo.G_contfreq[n] / 100.)) / 1000000.);

        }else{
            G_gloResStorage[n] = 0.;	// to set/reset bc. setStoragesToZero() deactivated
        }

        G_totalUnsatisfiedUse[n] = additionalOutIn.additionalOutputInput(n, 3);
        G_UnsatisfiedUsePrevYear[n] = additionalOutIn.additionalOutputInput(n, 4);// AEMS: set G_totalUnsatisfiedUse values
        G_reducedReturnFlow[n] = additionalOutIn.additionalOutputInput(n, 38);
        G_reducedReturnFlowPrevYear[n] = additionalOutIn.additionalOutputInput(n, 50);
        G_unsatisfiedNAsFromIrrig[n] = additionalOutIn.additionalOutputInput(n, 39);
        G_unsatisfiedNAsFromIrrigPrevYear[n] =  additionalOutIn.additionalOutputInput(n, 51);
        G_unsatisfiedNAsFromOtherSectors[n] = additionalOutIn.additionalOutputInput(n, 49); // Grid with unsatisfied NAs from other sectors than irrig associated with G_reducedReturnFlow
        G_unsatisfiedNAsFromOtherSectorsPrevYear[n] = additionalOutIn.additionalOutputInput(n, 52);
        //TODO: check what happens with the other arrays of the function setStoragesToZero (e.g. G_totalDesiredUseCell[n].
    }
}

void routingClass::dailyNUInit(const std::string input_directory, const short new_year,calibParamClass &calParam)// OE: added calParam
            {
    short month;

    // if this is a calibration run, allocated surface water and adapted net abstraction from groundwater should be
    // used.
    if (3 == options.subtract_use) {

        G_dailyNUs.read(input_directory + "/G_ACTUAL_NAS_" + std::to_string(new_year) + ".12.UNF0");
        G_dailyNUg.read(input_directory + "/G_ACTUAL_NAG_" + std::to_string(new_year) + ".12.UNF0");

    } else { // normal runs

        G_dailyNUs.read(input_directory + "/G_NETUSE_SW_m3_" + std::to_string(new_year) + ".12.UNF0");
        G_dailyNUg.read(input_directory + "/G_NETUSE_GW_m3_" + std::to_string(new_year) + ".12.UNF0");

        if(options.calc_wtemp == 2){
            G_yearlyWW.read(input_directory + "/G_ELEC_WW_m3_" + std::to_string(new_year) + ".UNF0");
            G_yearlyWC.read(input_directory + "/G_ELEC_WC_m3_" + std::to_string(new_year) + ".UNF0");
        }
    }

    G_monthlyWUIrrigFromSwb.read(input_directory + "/G_IRRIG_WITHDRAWAL_USE_SW_m3_" + std::to_string(new_year) + ".12.UNF0");

    G_monthlyCUIrrigFromSwb.read(input_directory + "/G_IRRIG_CONS_USE_SW_m3_" + std::to_string(new_year) + ".12.UNF0");

    // Cell-specific calibration parameters - Apply multiplier  // FP
    for (int n_local = 0; n_local < ng; n_local++) {
        for (month = 0; month < 12; month++) {
            G_dailyNUs(n_local,month) = calParam.getValue(M_NETABSSW, n_local) * G_dailyNUs(n_local,month);
            G_dailyNUg(n_local,month) = calParam.getValue(M_NETABSGW, n_local) * G_dailyNUg(n_local,month);
        }
    }

    // read file -> read GLWD units of global lakes and reservoirs (Option 28)
    G_glwdUnit.read(std::string(options.input_dir) + "/GLWDunits.UNF4");

    // CR 2015: Calculation of array G_dailyNUsAggregated(n,month) for outflow and riparian cells of global lakes/reservoirs and all other cells.
    // outflow cells: all positive values of dailyNUs are aggregated over riparian cells and are added to the value of the outflow cell
    // negative values in riparian cells: G_dailyNUsAggregated(n,month) = dailyNUs(n,month) (initial value)
    // positive values in riparian cells (no outflow cell): G_dailyNUsAggregated(n,month) = 0 (since value is added to outflow cell)
    // all other cells (neither riparian nor outflow cell): G_dailyNUsAggregated(n,month) = dailyNUs(n,month) (initial value)
    if (1 == options.aggrNUsGloLakResOpt) {

        //  cout << "options.aggrNUsGloLakResOpt == 1: Aggregation of NUs over riparian cells of global lakes/reservoirs" << endl;

        for (int n = 0; n < ng; n++) {
            for (month = 0; month < 12; month++) {
                G_dailyNUsAggregated(n,month) = G_dailyNUs(n,month);
            }
        }

        // Aggregation of NUs values (except negative values) over GLWD units
        // (GLWD units include all riparian cells of a global lake or reservoir)

        for (month = 0; month < 12; month++) {

            for (int n = 0; n < ng; n++) {

                if (G_lake_area[n] > 0. || G_reservoir_area[n] > 0.) {        // if cell is an outflow cell

                    for (int i = 0; i < ng; i++) {

                        if (n == i) {        // to avoid counting the NUs value of the outflow cell twice
                            continue;
                        }

                        if (G_glwdUnit[n] == G_glwdUnit[i]) {

                            if (G_dailyNUs(i,month) > 0.) {

                                G_dailyNUsAggregated(n,month) += G_dailyNUsAggregated(i,month);
                                G_dailyNUsAggregated(i,month) = 0.;

                            } //else (negative value): G_dailyNUsAggregated remains G_dailyNUs
                        }
                    }
                }
            }
        }
    }

// calculate daily values, if this is a calibration run, G_dailyNUs is also in m3
    for (month = 0; month <= 11; month++) {
        for (int n = 0; n <= ng - 1; n++) {
            G_dailyNUs(n,month) = G_dailyNUs(n,month) / (1000000000. * (double) numberOfDaysInMonth[month]);
            G_dailyNUg(n,month) = G_dailyNUg(n,month) / (1000000000. * (double) numberOfDaysInMonth[month]);
            if (1 == options.aggrNUsGloLakResOpt) {
                G_dailyNUsAggregated(n,month) =
                               G_dailyNUsAggregated(n,month) / (1000000000. * (double) numberOfDaysInMonth[month]);
            }
        }
    }
}

void routingClass::annualInit(const short year, int start_month, AdditionalOutputInputFile &additionalOutIn) {
    // FP20161018N002: Do NOT introduce nSpecBasins as local variable nBasins, as superbasin area calculation probably wrong & obsolete
    short month;  // index for months
    int n;  // index for grid cells
    char filename[250];

    // Part 1: Initialize always needed arrays
    for (n = 0; n < ng; n++) {

        // store the values of the grid, so that we can compare
        // at the end of the year if unsatisfied use of the previous year
        // has been satisfied with water availability from the current year
        // however, this should only occur in std mode. If in AEMS mode
        // G_UnsatisfiedUsePrevYear[n] and G_totalUnsatisfiedUse[n] are set in setStorages
        // at beginning and only in the years after the start_year same behavior as std mode.
        if (additionalOutIn.additionalfilestatus == 0) {
            G_UnsatisfiedUsePrevYear[n] = G_totalUnsatisfiedUse[n];
            G_unsatisfiedNAsFromIrrigPrevYear[n] = G_unsatisfiedNAsFromIrrig[n];
            G_unsatisfiedNAsFromOtherSectorsPrevYear[n] = G_unsatisfiedNAsFromOtherSectors[n];
            G_reducedReturnFlowPrevYear[n] = G_reducedReturnFlow[n];
        }
        else if(start_month==1||year>options.start_year){    // 1 equals Jan. here bc. configfile->startMonth btw. 1-12
            G_UnsatisfiedUsePrevYear[n] = G_totalUnsatisfiedUse[n]; //must be done every beginning of year even with additionalOutIn.additionalfilestatus==1
            G_unsatisfiedNAsFromIrrigPrevYear[n] = G_unsatisfiedNAsFromIrrig[n];
            G_unsatisfiedNAsFromOtherSectorsPrevYear[n] = G_unsatisfiedNAsFromOtherSectors[n];
            G_reducedReturnFlowPrevYear[n] = G_reducedReturnFlow[n];
        }



        G_potCellRunoff[n] = 0.;

        G_actualUse[n] = 0.; // used for new use_allocation (M.Hunger 2/2006)
        // FP20161018N002 Reservoir operation start years
        // Reservoir and regulated lakes area and storage capacity of current year
        G_reservoir_area[n] = 0.;
        G_stor_cap[n] = 0.;

    }

    //
    // FP20161018N002 Reservoir operation start years
    //
    // Part 2: Fill/initialize arrays for the current year
    // e.g. G_reservoir_area, G_stor_cap
    //
    // Status of algorithm depends on
    // whether naturalized runs are performed or whether reservoirs (incl. regulated lakes) are to be used
    // Here, only anthrogenic runs are of interest, as then year-specific values have to be considered
    if (options.antNatOpt == 0) { // anthropogenic runs: consider reservoirs and regulated lakes as reservoirs

        // Anthropogenic run
        // Use reservoirs and regulated lakes
        if (options.resOpt == 1) {

            // Always for regulated lakes, reservoir_area and storage capacity needs to be considered
            for (n = 0; n < ng; n++) {

                // 1. Do not consider another time land cover fraction of regulated lakes
                //      because they are already in G_glo_res (that should be used in subsequent treatment, if concerned)
                //      But as G_reservoir_area is initialized with zero every year, the transfer has to be done every year
                //      from G_reservoir_area_full as read from the respective file (that always should contain regulated lakes, too)
                // 2. Treat area and storage capacity
                if (G_reg_lake_status[n] == 1) {

                    //reservoir area is needed anyhow from regulated lakes, otherwise they cannot be calculated
                    //in the reservoir algorithm
                    G_reservoir_area[n] = G_reservoir_area_full[n];

                    G_stor_cap[n] = G_stor_cap_full[n];
                }
            }

            // Define desired year of reservoir data (fraction, area, storage capacity)
            // Depending on reservoir start year option yes/no
            if (options.resYearOpt == 1) { // if reservoir year option is switched on

                // Limit years to valid domain specified by
                // options.resYearFirstToUse & options.resYearLastToUse
                // Maximum = last year
                if (year > options.resYearLastToUse) {
                    resYearCurrentToUse = options.resYearLastToUse;
                }
                    // Minimum = first year
                else if (year < options.resYearFirstToUse) {
                    resYearCurrentToUse = options.resYearFirstToUse;
                }
                    // Within domain: Current year
                else {
                    resYearCurrentToUse = year;
                }

            }
            else if (options.resYearOpt == 0) {

                // Using reference year
                resYearCurrentToUse = options.resYearReference;  // ISIMIP21a: G_RES_2000.UNF0

            }

            // Read desired reservoir and regulated lake fractions
            // within domain of years, e.g. current, min, max, reference year
            // ATTENTION: Fractions should be consistent to condition resYearCurrentToUse >= G_res_start_year[n]
            G_glo_res.read(std::string(options.input_dir) + "/G_RES/G_RES_" + std::to_string(resYearCurrentToUse) + ".UNF0");

            // Transfer the values at start of the first model year to array of previous
            if (0 == statusStarted_updateGloResPrevYear) {
                routingClass::updateGloResPrevYear_pct();  // copying values of current year
                statusStarted_updateGloResPrevYear = 1; // indicate for the next iteration that first year has been treated / model has started
            }

            // Treat area and storage capacity (outflow cells)
            for (n = 0; n < ng; n++) {
                // (1.1) update G_reservoir_area and G_stor_cap if reservoir is in operation, e.g. since reference year
                if (resYearCurrentToUse >= G_res_start_year[n]) {
                    G_reservoir_area[n] = G_reservoir_area_full[n];
                    G_stor_cap[n] = G_stor_cap_full[n];
                }
                // (1.2) Define current storage to be the maximum storage capacity of the selected year (current, max, min, reference)
                G_actualStorageCapacity[n] = G_stor_cap[n];
            }

            // In case of erroneous attributes of a reservoir, treat it as a global lake:
            // Unknown reservoir type & negative mean outflow
            // Only look at conditions of current year
            for (n = 0; n < ng; n++) {
                if ((G_reservoir_area[n] > 0.) && (resYearCurrentToUse >= G_res_start_year[n])) {

                    if (G_res_type[n] + 0 == 0) {//only IF reservoir type is UNKNOWN (should not be the case anyhow)
                        cerr << "ERROR: routing.annuaInit: Year " << year
                             << " - Reservoir management: G_res_type is unknown at n: " << n << endl;
                        cerr << "Reservoir will be treated as global lake!" << endl;
                        cerr << "Adding " << G_reservoir_area_full[n] << " km2 to global lake area " << G_lake_area[n]
                             << endl;
                        G_lake_area[n] += G_reservoir_area_full[n]; //use reservoir area as lake area, or sum up reservoir and lake area
                        G_reservoir_area[n] = 0.; //set reservoir area to zero (inhibits further accounting elsewhere)
                        G_reservoir_area_full[n] = 0.; //set full reservoir area to zero (inhibits further accounting in next year)
                        cerr << "Final global lake area " << G_lake_area[n] << endl;
                    } else if (G_mean_outflow[n] <= 0.) {
                        if  (!(G_LDD[n] == -1))
                        {
                            cerr << "ERROR: routing.annuaInit: Year " << year
                                 << " - Reservoir management: G_mean_outflow in reservoir is " << G_mean_outflow[n]
                                 << " i.e. less than or equal to zero at n: " << n << endl;
                            cerr << "Reservoir will be treated as global lake!" << endl;
                            cerr << "Adding " << G_reservoir_area_full[n] << " km2 to global lake area "
                                 << G_lake_area[n]
                                 << endl;
                        }
                        G_lake_area[n] += G_reservoir_area_full[n]; //use reservoir area as lake area, or sum up reservoir and lake area
                        G_reservoir_area[n] = 0.; //set reservoir area to zero (inhibits further accounting elsewhere)
                        G_reservoir_area_full[n] = 0.; //set full reservoir area to zero (inhibits further accounting in next year)
                        if (!(G_LDD[n] == -1))
                        {
                            cerr << "Final global lake area " << G_lake_area[n] << endl;
                        }

                    }

                }
                // end if option conditions
            }
            // end loop over grid cells

        }
        // end if options.resOpt == 1

    }
    //end calculate only for anthropogenic runs (options.antNatOpt == 0)

    // Only for runs which consider glaciers
    if (options.glacierOpt) {
        for (n = 0; n < ng; n++) {
            // fill G_glacAreaFracPrevYear at the beginning of each year, before calculation of glacAreaFrac of current year
            G_glacAreaFracPrevYear[n] = G_glacAreaFrac[n];
            if (glacierYear.G_glacier_area[n] > 0.) { // for glacierized cells only
                if (glacierYear.G_glacier_area[n] > (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n]
                - geo.G_fwaterfreq_const[n] - G_glo_wetland[n] - G_loc_wetland[n]) / 100.)) {
                    // if input glacier area is larger than cell continental area, then set it to be equal to cell continental area
                    G_glacAdaptedArea[n] = geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n]
                                   - geo.G_fwaterfreq_const[n] - G_glo_wetland[n] - G_loc_wetland[n]) / 100.;

                } else
                    G_glacAdaptedArea[n] = glacierYear.G_glacier_area[n];
                G_glacAreaFrac[n] = G_glacAdaptedArea[n] / geo.areaOfCellByArrayPos(n) * 100.; // calculate glacier area as a fraction of cell area (%)
            } else
                G_glacAdaptedArea[n] = 0.; // make sure that glacier area is equal to 0 in non-glacierized cells
        }
    } //end if glacierOpt


    // Count the number of cells with reservoir area
    short count = 0;
    for (n = 0; n < ng; n++) {
        if (G_reservoir_area[n] > 0.) count++;
    }
    cout << "routing.annualInit: " << count << " reservoirs (area > 0) are considered in this year " << year << endl;

    //
    // FP20161018N002 Reservoir operation start years
    //

    // Part 3: Always used arrays that depend on G_glo_res (Unit: percent fraction of grid cell area)
    for (n = 0; n < ng; n++) {
        // Now, calculate land area fraction (instead of once at model initialisation at routing.init because G_glo_res can be dynamic)

        // (1) Land area fraction excluding lakes and reservoirs/regulated lakes (G_glo_res includes G_reg_lake)
        // fixed G_fGloLake = G_glo_lake/100 and newly read-in G_glo_res
        G_landWaterExclGloLakAreaFrac[n] = (geo.G_contfreq[n] - G_glo_lake[n] - G_glo_res[n]);
        if (G_landWaterExclGloLakAreaFrac[n] < limit_errorstrange_landAreaFrac_pct) {
            cerr << "Warning: routing.annualInit: Year " << year
                 << " - Strange number for lake/wetland fraction 'G_landWaterExclGloLakAreaFrac' for cell " << n << ": "
                 << G_landWaterExclGloLakAreaFrac[n] << ", limit: " << limit_errorstrange_landAreaFrac_pct << endl;
            cerr << "Warning: routing.annualInit: Year " << year << " - Set values < 0 to zero" << endl;
        }
        if (G_landWaterExclGloLakAreaFrac[n] < 0.)
            G_landWaterExclGloLakAreaFrac[n] = 0.;

        // (2)
        // (2.1) Fraction of water and land: G_fwaterfreq and G_landfreq (model start and updatedevery year)
        // and associated changes
        //G_fwaterfreq is sum of G_glo_lake, G_loc_lake, and G_glo_res (includes G_reg_lake)
        //G_landfreq is G_cont - G_fwaterfreq, contains wetlands (is used to calculate PET in daily.cpp)
        // ATTENTION: Any change to previous year is assumed to be an increase in G_glo_res

        // (2.2) G_landAreaFrac for the first year = at model start
        // For first year

        if (0 == statusStarted_landfreq_fwaterfreq_landAreaFrac[n]) {

            // Water and area fractions (in percent)
            // fixed G_fLocLake from first reading in routing.init: G_loc_lake/100
            // fixed G_fGloLake = G_glo_lake/100
            // as read-in G_glo_res
            G_fwaterfreq[n] = G_glo_lake[n] + G_loc_lake[n] + G_glo_res[n];
            G_landfreq[n] = geo.G_contfreq[n] - G_fwaterfreq[n];


            G_fwaterfreq_change[n] = 0.;
            G_landfreq_change[n] = 0.;

            // Initial land area fraction (in percent)

            G_landAreaFrac[n] = /*(double)*/ (geo.G_contfreq[n]
                                          - (G_glo_lake[n] + G_glo_wetland[n] + G_loc_lake[n] + G_loc_wetland[n] + G_glo_res[n]));

            G_landAreaFracPrevTimestep[n] = 0.;

            // Only for runs which consider glaciers
            if (options.glacierOpt == 1) {
                G_landAreaFrac[n] = G_landAreaFrac[n] - G_glacAreaFrac[n]; // reduce land area fraction by glacier area fraction
            }

            // Treat undesired values
            if (G_landAreaFrac[n] < limit_errorstrange_landAreaFrac_pct) {
                cerr << "Warning: routing.annualInit: Year " << year
                     << " - Strange number for land area fraction 'G_landAreaFrac' for cell " << n << ": "
                     << G_landAreaFrac[n] << ", limit: " << limit_errorstrange_landAreaFrac_pct << endl;
                cerr << "Warning: routing.annualInit: Year " << year << " - Set values < 0 to zero" << endl;
            }
            if (G_landAreaFrac[n] < 0.)
                G_landAreaFrac[n] = 0.;

            // update indicator
            statusStarted_landfreq_fwaterfreq_landAreaFrac[n] = 1;

        }
            // For all years besides the first year, transfer data from previous year and recalculate
        else {

            G_fwaterfreq_prevyear[n] = G_fwaterfreq[n];
            G_landfreq_prevyear[n] = G_landfreq[n];

            // Consider fractional routing fractions of local lakes (as done in previous step, besides for the first year)
            // variable G_fLocLake from previous time step
            // fixed G_fGloLake = G_glo_lake/100
            // as read-in G_glo_res
            G_fwaterfreq[n] = G_glo_lake[n] + (G_fLocLake[n] * 100.) + G_glo_res[n];
            G_landfreq[n] = geo.G_contfreq[n] - G_fwaterfreq[n];

            G_fwaterfreq_change[n] = G_fwaterfreq[n] - G_fwaterfreq_prevyear[n];
            G_landfreq_change[n] = G_landfreq[n] - G_landfreq_prevyear[n];

            // Only for runs which consider glaciers
            if (options.glacierOpt == 1) {
                // Calculate glacier area fraction change for years after first year
                // If there was no glacier(s) on previous year, then G_glacAreaFrac_change is equal to G_glacAreaFrac
                G_glacAreaFrac_change[n] = G_glacAreaFrac[n] - G_glacAreaFracPrevYear[n];

                if (G_glacAreaFrac_change[n] < -0.00001 || G_glacAreaFrac_change[n] > 0.00001) {
                    // If glacier area fraction has increased or decreased with reference to the previous year,
                    // adapt land area fraction by glacier area fraction change
                    G_landAreaFrac[n] = G_landAreaFrac[n] - G_glacAreaFrac_change[n];

                    // Treat undesired values
                    if (G_landAreaFrac[n] < limit_errorstrange_landAreaFrac_pct) {
                        cerr << "Warning: routing.annualInit: (options.glacierOpt == 1) & (G_glacAreaFrac_change[n] != 0.) - Strange number for land area fraction 'G_landAreaFrac' for cell "
                        << n << ": " << G_landAreaFrac[n] << endl;
                    }
                    if (G_landAreaFrac[n] < 0.)
                        G_landAreaFrac[n] = 0.;

                    // Calculate land area fraction change
                    //G_landAreaFrac_change[n] = G_landAreaFrac[n] - G_landAreaFracPrevTimestep[n];

                    // Set land area fraction of next time step to adapted land area fraction
                    G_landAreaFracNextTimestep[n] = G_landAreaFrac[n];

                } // end if G_glacAreaFracChange != 0
            } //end if glacierOpt

        }
    }
    // end loop over grid cells

    // Part 4: Needed only in case of anthropogenic runs using reservoirs
    if (options.antNatOpt == 0) { // anthropogenic runs: consider reservoirs and regulated lakes as reservoirs
        //    cout << "routing.annualInit Part 4: if options.antNatOpt == 0, anthropogenic run: consider reservoirs and regulated lakes as reservoirs" << endl;

        // Anthropogenic run
        // Use reservoirs and regulated lakes
        if (options.resOpt == 1) {

            // HERE: Adapt code in order to add water from reduced land fraction to reservoir volume
            // include also G_riverAreaFrac[n]
            for (n = 0; n < ng; n++) {

                //check if change arrays etc. are already initialized
                //in case AEMS is running over year break, routing.updateGloResPrevYear_pct() destroys potentially added reservoirs
                if (additionalOutIn.additionalfilestatus == 1) {
                    G_glores_prevyear[n] = additionalOutIn.additionalOutputInput(n, 44);
                }
                // Calculate (possible) increase in reservoir area at start of the year
                // indicating operation start of a new reservoir
                if (start_month == 1 || year > options.start_year){

                    G_glores_change[n] = G_glo_res[n] - G_glores_prevyear[n];
                    // Perform changes in water storages only if a new reservoir started operation
                    if (G_glores_change[n] > 0.) {

                        // Subtract additional reservoir fraction from land area fraction
                        // (1) Land area fraction including wetlands, local lakes, rivers
                        // day 0 (current day, in current year: year / Jan / 01): G_glo_res[n]
                        // From routing.updateLandAreaFrac:
                        // day-1 (previous day, in previous year: year-1 / Dec / 31): G_landAreaFrac[n]
                        // day-2 (day before previous day, in previous year: year-1 / Dec / 30): G_landAreaFracPrevTimestep[n]

                        // 1 equals Jan. here bc. configfile->startMonth btw. 1-12
                        if (G_glores_prevyear[n] > 0.) { // in case there was already a reservoir fraction, adapt only the changes.
                            G_landAreaFrac[n] = G_landAreaFrac[n] - G_glores_change[n];
                        }
                        else { //reservoir fraction increased from 0.
                            G_landAreaFrac[n] = G_landAreaFrac[n] - G_glo_res[n];
                        }

                        G_glores_change[n]=0.;
                        // Treat undesired values
                        if (G_landAreaFrac[n] < limit_errorstrange_landAreaFrac_pct) {
                            cerr
                                           << "Warning: routing.annualInit: (options.antNatOpt == 0) & (options.resOpt == 1) & (G_glores_change[n] > 0.) - Strange number for land area fraction 'G_landAreaFrac' for cell "
                                           << n << ": " << G_landAreaFrac[n] << endl;
                        }
                        if (G_landAreaFrac[n] < 0.)
                            G_landAreaFrac[n] = 0.;

                        // START ISIMIP2b POSSIBLE uncomment 2016.
                        // ATTENTION: POSSIBLY For ISIMIP2B, GLWD units consistent to outflow cells are missing.
                        // Therefore, the calculation of net change and subsequent transfer of storage to reservoir storage is commented out
                        // Scaling of storage as usual in daily.calcNewDay
                        //

                        // Calculate net change
                        // Changes of day 0 from wetland / loc.lakes/rivers (fractional routing) + new reservoir
                        // with respect to day-1
                        G_landAreaFrac_change[n] = G_landAreaFrac[n] - G_landAreaFracPrevTimestep[n];

                        // In case that net change is not negative,
                        // assume that it is a result of fractional routing only, reset G_landAreaFrac to value of day-2
                        // and let daily.calcNewDay execute as usual, implicitly applying no changes (from fractional routing
                        if (G_landAreaFrac_change[n] >= 0.) {
                            G_landAreaFrac[n] = G_landAreaFracPrevTimestep[n];
                            G_landAreaFrac_change[n] = 0.;
                        }
                            // Net change is negative (reservoir and/or added: Adapt storage
                        else {

                            // Calculate absolute storage volume in km3 from reduced ("lost") land surface fraction
                            canopyWaterContent_change_km3 =
                                           dailyWaterBalance.G_canopyWaterContent[n] * geo.areaOfCellByArrayPos(n) /
                                           1000000.0 * ((-G_landAreaFrac_change[n]) / 100.0);
                            soilWaterContent_change_km3 =
                                           dailyWaterBalance.G_soilWaterContent[n] * geo.areaOfCellByArrayPos(n) /
                                           1000000.0 * ((-G_landAreaFrac_change[n]) / 100.0);
                            snowWaterContent_change_km3 =
                                           dailyWaterBalance.G_snow[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 *
                                           ((-G_landAreaFrac_change[n]) / 100.0);
                            SnowInElevation_n_change_km3 = 0.;
                            for (short elev = 0; elev < 101; elev++) {
                                SnowInElevation_elev_change_km3 = dailyWaterBalance.G_SnowInElevation(n,elev) / 100.0 *
                                                                  geo.areaOfCellByArrayPos(n) / 1000000.0 *
                                                                  ((-G_landAreaFrac_change[n]) / 100.0);
                                SnowInElevation_n_change_km3 += SnowInElevation_elev_change_km3;
                            }

                            if (((SnowInElevation_n_change_km3 - snowWaterContent_change_km3) >
                                 limit_errorstrange_snowWaterContent_change_km3)
                                || ((snowWaterContent_change_km3 - SnowInElevation_n_change_km3) >
                                    limit_errorstrange_snowWaterContent_change_km3)) {
                                cerr
                                               << "Warning: routing.annualInit: (options.antNatOpt == 0) & (options.resOpt == 1) & (G_glores_change[n] > 0.) - Exceeded limit of difference between 'SnowInElevation_n_change_km3' and 'snowWaterContent_change_km3' for cell "
                                               << n << ": " << SnowInElevation_n_change_km3 << ", "
                                               << snowWaterContent_change_km3 << ", limit: "
                                               << limit_errorstrange_snowWaterContent_change_km3 << endl;
                            }

                            // Add storage to reservoir storage volume located in outflow cell, as identified via GLWDunit grid

                            G_gloResStorage[G_outflow_cell_assignment[n] - 1] +=
                                           canopyWaterContent_change_km3 + soilWaterContent_change_km3 +
                                           snowWaterContent_change_km3; //

                            // ATTENTION: dailyWaterBalance.G_canopyWaterContent etc. are in units of Millimeter (mm)
                            // Through the reduction of land surface fraction, water storage is implicitly reduced
                            // To ensure that no further scaling is done within daily.calcNewDay, set equal previous value (day-1) to current newly determined uvalue (day 0)

                            G_landAreaFracNextTimestep[n] = G_landAreaFrac[n];

                        }
                        // end else
                    }
                    // end if (G_glores_change[n] > 0.)
                }
                // end if (start_month == 1 || year > options.start_year)
            }
            // end loop over grid cells

        }
        // end if options.resOpt == 1

    }
    //end calculate only for anthropogenic runs (options.antNatOpt == 0)

    // END ISIMIP21a insert ersion

    // Initialize arrays for writing output
    if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {
        for (n = 0; n < ng; n++)
            for (int month = 0; month < 12; month++) {
                G_monthlyConsistentPrecip(n,month) = 0.;

                G_monthlyRiverAvail(n,month) = 0.;
                G_monthlyPotCellRunoff(n,month) = 0.;
                G_monthlyPotCellRunoffDeficit(n,month) = 0.;
                G_monthlyCellRunoff(n,month) = 0.;
                G_monthlyCellSurfaceRunoff(n,month) = 0.;
                G_monthlyGwrunSurfrun(n,month) = 0.;

                G_monthlyCellAET(n,month) = 0.;
                G_monthlyCellAETWCa(n,month) = 0.;
                G_monthlyOpenWaterEvap(n,month) = 0.;
                G_monthlyVelocity(n,month) = 0.;
                G_monthlySurfStor(n,month) = 0.;
                G_monthlySurfStor_mm(n,month) = 0.;
                G_monthlyMinRiverAvail(n,month) = 0.;
                G_monthlyMaxRiverAvail(n,month) = 0.;
                if (options.outWaterTempMonthly) G_monthlyMeanWaterTempRiver(n,month) = 0.;
                if (options.outWaterTempMonthlyAllSWB) {
                    G_monthlyMeanWaterTempLocLake(n, month) = 0.;
                    G_monthlyMeanWaterTempLocWetl(n, month) = 0.;
                    G_monthlyMeanWaterTempGloLake(n, month) = 0.;
                    G_monthlyMeanWaterTempReservoir(n, month) = 0.;
                    G_monthlyMeanWaterTempGloWetl(n, month) = 0.;
                }
                G_monthlyGwrSwb(n,month) = 0.;
                G_monthlyLandAreaFrac(n,month) = 0.;
                G_monthlyFswb(n,month) = 0.;
                G_monthlyGwRunoff(n,month) = 0.;
                G_monthlyLocWetlExtent(n,month) = 0.;
                G_monthlyGloWetlExtent(n,month) = 0.;
                if (options.outRiverPET == 1) G_monthlyRiverAreaFrac(n,month) = 0.;
                if (options.outRiverPET == 1) G_monthlyRiverPET(n,month) = 0.;

                if (options.resOpt == 1) {
                    G_monthlyResStorage(n,month) = 0.;
                    G_monthlyResStorageRatio(n,month) = 0.;
                    G_monthlyResInflow(n,month) = 0.;
                    G_monthlyResOutflow(n,month) = 0.;
                }

                if (options.glacierOpt == 1) {
                    G_monthlyGlacierArea(n,month) = 0.;
                    G_monthlyGlacierAreaFrac(n,month) = 0.;
                    G_monthlyGlacierRunoff(n,month) = 0.;
                    G_monthlyGlacierPrecip(n,month) = 0.;
                    G_monthlyGlacierStorage(n,month) = 0.;
                }

                if (options.outRiverInUpstream)
                    G_monthlyRiverInUpstream(n,month) = 0.;
                G_monthlyRiverStorage(n,month) = 0.;
                G_monthlyGwStorage(n,month) = 0.;
                G_monthlyLocLakeStorage(n,month) = 0.;
                G_monthlyLocWetlStorage(n,month) = 0.;
                G_monthlyGloLakeStorage(n,month) = 0.;
                G_monthlyGloWetlStorage(n,month) = 0.;
                G_monthlySatisfiedUse(n,month) = 0.;
                G_monthlyActualUse(n,month) = 0.;
                G_monthlyNUg(n,month) = 0.;
                G_monthlyAllocatedUse(n,month) = 0.;
                G_monthlySatisAllocatedUseinSecondCell(n,month) = 0.;
                G_monthlyRedistributeAllocUse(n, month) = 0.;
                G_monthlyRedistributeGLWDUse(n, month) = 0.;
            }

    }

}

// HMS 2015-03-02 I do not understand why this is called monthlyInit - should it be named better daily31outInit? (same in daily.cpp)
// (FP 2015: perhaps naming is because the data arrays contain data from one month)
// daily output option 31 start
void routingClass::daily31outInit() {
    short d;
    int n;
    for (n = 0; n < ng; n++)
        for (d = 0; d < 31; d++) {
            if (options.outPrecDaily) G_daily31ConsistentPrecip(n,d) = 0.;

            if (options.outCellRunoffDaily) G_daily31CellRunoff(n,d) = 0.;
            if (options.outCellSurfaceDaily) G_daily31CellSurfaceRunoff(n,d) = 0.;
            if (options.outGWRunoffDaily) G_daily31GwRunoff(n,d) = 0.;
            if (options.outRiverAvailDaily) G_daily31RiverAvail(n,d) = 0.;
            if (options.outRiverVeloDaily) G_daily31Velocity(n,d) = 0.;
            if (options.outWaterTempDaily) G_daily31WaterTemp(n,d) = -9999.;
            if (options.outWaterTempDailyAllSWB){
                G_daily31locLakeTemp(n,d) = -9999.;
                G_daily31locWetlTemp(n,d) = -9999.;
                G_daily31gloLakeTemp(n,d) = -9999.;
                G_daily31reservoirTemp(n,d) = -9999.;
                G_daily31gloWetlandTemp(n,d) = -9999.;
            }
            if (options.outCellAETDaily) G_daily31CellAET(n,d) = 0.;
            if (options.outSurfStorDaily) G_daily31SurfStor(n,d) = 0.;
            if (options.outSingleStoragesDaily) G_daily31LocLakeStor(n,d) = 0.;
            if (options.outSingleStoragesDaily) G_daily31LocWetlStor(n,d) = 0.;
            if (options.outSingleStoragesDaily) G_daily31GloLakeStor(n,d) = 0.;
            if (options.outSingleStoragesDaily) G_daily31GloWetlStor(n,d) = 0.;
            if (options.outSingleStoragesDaily) G_daily31RiverStor(n,d) = 0.;
            if (options.outSingleStoragesDaily) G_daily31ResStor(n,d) = 0.;
            if (options.outSingleStoragesDaily || options.outGWStorageDaily) G_daily31GwStor(n,d) = 0.;
            if (options.outTotalWaterInStoragesDaily_km3 ||
                options.outTotalWaterInStoragesDaily_mm)
                G_daily31TotalWaterInStorages_km3(n,d) = 0.;
            if (options.outTotalWaterInStoragesDaily_mm) G_daily31TotalWaterInStorages_mm(n,d) = 0.;
            if (options.outGwrSwbDaily) G_daily31GwrSwb(n,d) = 0.;
            if (options.outFswbDaily) G_daily31Fswb(n,d) = 0.;
            if (options.outLandAreaFracDaily) G_daily31LandAreaFrac(n,d) = 0.;
            if (options.outGwrunSurfrunDaily) G_daily31GwrunSurfrun(n,d) = 0.;
            if (options.outCellAETWCaDaily) G_daily31CellAETWCa(n,d) = 0.;
        }
}

// to store daily TWS value for first and last day in year.
void routingClass::startendOutInit() {
    short d;
    int n;

    for (n = 0; n < ng; n++)
        for (d = 0; d < 1; d++) {
            G_startendTotalWaterInStorages_km3(n,d) = 0.;
        }
}

// daily output option 365 start
void routingClass::daily365outInit() {  // HMS 2015-03-02 where is this class called (not in watergap.cpp)
    short d;
    int n;

    for (n = 0; n < ng; n++)
        for (d = 0; d < 365; d++) {

            if ((options.outPrecDaily) ||
                ((options.scoutcPrecip) && (2 == options.day_store)))
                G_daily365ConsistentPrecip(n,d) = 0.;
            if ((options.outCellRunoffDaily) ||
                (((options.scoutCellRunoffkm3) || (options.scoutCellRunoffmm)) && (2 == options.day_store)))
                G_daily365CellRunoff(n,d) = 0.;
            if ((options.outCellRunoffDaily) ||
                ((options.scoutCellRunoffmm) && (2 == options.day_store)))
                G_daily365CellRunoff_mm(n,d) = 0.;
            if ((options.outCellSurfaceDaily) ||
                ((options.scoutCellSRunoff) && (2 == options.day_store)))
                G_daily365CellSurfaceRunoff(n,d) = 0.;
            if ((options.outLandAreaFracDaily) ||
                ((options.scoutLandAreaFrac) && (2 == options.day_store)))
                G_daily365LandAreaFrac(n,d) = 0.;
            if ((options.outGWRunoffDaily) || ((options.scoutGwRunoff) && (2 == options.day_store)))
                G_daily365GwRunoff(n,d) = 0.;
            if ((options.outRiverAvailDaily) || ((options.scoutQ) && (2 == options.day_store)))
                G_daily365RiverAvail(n,d) = 0.;
            if ((options.outRiverVeloDaily) || ((options.scoutFlowVelo) && (2 == options.day_store)))
                G_daily365Velocity(n,d) = 0.;
            if ((options.outCellAETDaily) || ((options.scoutCellAET) && (2 == options.day_store)))
                G_daily365CellAET(n,d) = 0.;
            if ((options.outSurfStorDaily) ||
                ((options.scoutSurfaceStor) && (2 == options.day_store)))
                G_daily365SurfStor(n,d) = 0.;
            if ((options.outSingleStoragesDaily) ||
                ((options.scoutLocLake) && (2 == options.day_store)))
                G_daily365LocLakeStor(n,d) = 0.;
            if ((options.outSingleStoragesDaily) ||
                ((options.scoutLocWet) && (2 == options.day_store)))
                G_daily365LocWetlStor(n,d) = 0.;
            if ((options.outSingleStoragesDaily) ||
                ((options.scoutGloLake) && (2 == options.day_store)))
                G_daily365GloLakeStor(n,d) = 0.;
            if ((options.outSingleStoragesDaily) ||
                ((options.scoutGloWet) && (2 == options.day_store)))
                G_daily365GloWetlStor(n,d) = 0.;
            if ((options.outSingleStoragesDaily) ||
                ((options.scoutRiver) && (2 == options.day_store)))
                G_daily365RiverStor(n,d) = 0.;
            if ((options.outSingleStoragesDaily) ||
                ((options.scoutReservoir) && (2 == options.day_store)))
                G_daily365ResStor(n,d) = 0.;
            if ((options.outSingleStoragesDaily || options.outGWStorageDaily) ||
                ((options.scoutGwStor) && (2 == options.day_store)))
                G_daily365GwStor(n,d) = 0.;
            if ((options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm) ||
                (((options.scoutTWSkm3) || (options.scoutTWSmm)) && (2 == options.day_store)))
                G_daily365TotalWaterInStorages_km3(n,d) = 0.;
            if ((options.outTotalWaterInStoragesDaily_mm) ||
                ((options.scoutTWSmm) && (2 == options.day_store)))
                G_daily365TotalWaterInStorages_mm(n,d) = 0.;
            if ((options.outGwrSwbDaily) || ((options.scoutGwrSwb) && (2 == options.day_store)))
                G_daily365GwrSwb(n,d) = 0.;
            if ((options.outFswbDaily) || ((options.scoutFswb) && (2 == options.day_store))) G_daily365Fswb(n,d) = 0.;
            if (options.outGwrunSurfrunDaily) G_daily365GwrunSurfrun(n,d) = 0.;
            if (options.outCellAETWCaDaily) G_daily365CellAETWCa(n,d) = 0.;
            if (options.outGlacierStorageDaily) G_daily365GlacierStorage(n,d) = 0.;
            if (options.outGlacierStorageDaily_mm) G_daily365GlacierStorage_mm(n,d) = 0.;
        }
    // endfor loop over 365 days


}
//}
// end daily365outInit
// daily output option 365 end

void routingClass::routing(const short year, short day, short month, const short day_in_month,
                           const short last_day_in_month, WghmStateFile &wghmState,
                           AdditionalOutputInputFile &additionalOutIn, const short readinstatus,calibParamClass &calParam)// OE: added calParam
            {
    // routing through river network


    const short int first_day_in_month[12] = {1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335};

    Grid<double> G_localRunoff;
    Grid<double> G_localRunoffIntoRiver; // for options.fractionalRouting
    Grid<double> G_localGWRunoff; // needed for different values for sub-/surface runoff
    Grid<double> G_localGWRunoffIntoRiver; // for options.fractionalRouting
    Grid<double> G_glacRunoff; // for options.glacierOpt

    // initialize
    G_localRunoff.fill(0.0);
    G_localRunoffIntoRiver.fill(0.0);
    G_localGWRunoff.fill(0.0);
    G_localGWRunoffIntoRiver.fill(0.0);
    G_glacRunoff.fill(0.0);

    double CellRunoff = 0.; // actual cell runoff generated within the cell
    double cellArea = 0., transportedVolume = 0., dailyRiverEvapo = 0.;
    double inflow = 0., inflowUpstream = 0., totalInflow = 0., outflow = 0., maxStorage = 0.;
    double K = 0.; // for new routing-equation
    short b = 0; // counters for nSpecBasins (avoid dynamic memory allocation)
    double potCellRunoff = 0.;

    //for water temp calculation
    double outflowlocLake = 0.;
    double outflowlocWetl = 0.;
    double outflowGloLake = 0.;
    double outflowGloWetland = 0.;

    double Kswbgw = 10.; // constant groundwater recharge below surface water bodies in mm/d --> 10 mm per day as gwr

    double fswbFracCorr = 0.; // correction for swb (loclak, locwet, glowet) if riverAreaFrac can not be satisfied
    double riverAreaFracDeficit = 0.; //to correct river area fraction if necessary


    // open water PET reduction (2.1f)
    double locLakeAreaReductionFactor = 0., gloLakeEvapoReductionFactor = 0.,
                   locWetlAreaReductionFactor = 0., gloWetlAreaReductionFactor = 0.;
    double gloResEvapoReductionFactor = 0., riverAreaReductionFactor = 0.;
    //CR: introduced in 2.2b; openwaterPET, gwr, and RemainingUse are included in one storage equation.
    //If they cannot be fully satisfied, they are reduced by the ReductionFactor = PETgwrMax/PETgwr.
    double PETgwr = 0., PETgwrMax = 0., PETgwrRedFactor = 0., PETgwrRemUseRedFactor = 0.;
    double groundwater_runoff_km3 = 0.;
    double G_groundwaterStoragePrevStep = 0.;
    // WG22b: GWR beneath SWB [km3] is included in the GW storage equation of the next timestep
    double netGWin = 0.;
    double dailyLocalSurfaceRunoff = 0.;    // in km3

    // set arrays with daily values to zero
    // this has to be done already here because
    // these arrays are summed up with finer time steps later
    // Reservoir operation start years
    // Homogenized index n to general b for "nSpecBasins"
    for (b = 0; b < nSpecBasins; b++) {
        dailyRiverDischarge(b,day - 1) = 0.;
        if (options.day_store == 1) {
            dailyLocLakeStorage(b,day - 1) = 0.;
            dailyLocWetlStorage(b,day - 1) = 0.;
            dailyGloLakeStorage(b,day - 1) = 0.;
            dailyGloWetlStorage(b,day - 1) = 0.;
            dailyRiverVelocity(b,day - 1) = 0.;    // velocity
            if (options.resOpt == 1)
                dailyResStorage(b,day - 1) = 0.;
        }
        // end if options.day_store == 1
    }
    // end loop over nSpecBasins

    //AEMS:
    //#pragma omp parallel for shared(G_localRunoff, G_localGWRunoff, dailyWaterBalance, geo, options, month)

    if ((additionalOutIn.additionalfilestatus == 1) && (1 == readinstatus)) {

        for (int n = 0; n < ng; n++) {

            // G_locWetlAreaReductionFactor[n] = additionalOutIn.additionalOutputInput(n, 8);
            // G_gloLakeEvapoReductionFactor[n] = additionalOutIn.additionalOutputInput(n, 9);
            // G_groundwaterStorage[n] = additionalOutIn.additionalOutputInput(n, 10);
            // G_locLakeAreaReductionFactor[n] = additionalOutIn.additionalOutputInput(n, 11);
            // G_gloWetlAreaReductionFactor[n] = additionalOutIn.additionalOutputInput(n, 12);
            // G_gloResEvapoReductionFactor[n] = additionalOutIn.additionalOutputInput(n, 13);
            //G_fswbInit[n] = additionalOutIn.additionalOutputInput(n, 14);
            // G_gloWetlStorage[n] = additionalOutIn.additionalOutputInput(n, 15);
            // G_fswbLandAreaFrac[n] = additionalOutIn.additionalOutputInput(n, 16);
            // G_locWetlStorage[n] = additionalOutIn.additionalOutputInput(n, 17);
            // G_locLakeStorage[n] = additionalOutIn.additionalOutputInput(n, 18);
            // G_riverStorage[n] = additionalOutIn.additionalOutputInput(n, 19);
            // G_gloLakeStorage[n] = additionalOutIn.additionalOutputInput(n, 21);
            // G_gloResStorage[n] = additionalOutIn.additionalOutputInput(n, 23);
            G_UnsatAllocUse[n] = additionalOutIn.additionalOutputInput(n, 27);
            G_AllocatedUse[n] = additionalOutIn.additionalOutputInput(n, 28);
            G_SecondCell[n] = additionalOutIn.additionalOutputInput(n, 29);
            G_daily_allocatedUseNextDay[n] = additionalOutIn.additionalOutputInput(n, 31);
            G_daily_UnsatAllocUseNextDay[n] = additionalOutIn.additionalOutputInput(n, 30);
            G_AllocUseToNeigborCell[n] = additionalOutIn.additionalOutputInput(n, 32);
            // G_fswbLandAreaFracNextTimestep[n] = additionalOutIn.additionalOutputInput(n, 33);
            G_PrevUnsatAllocUse[n] = additionalOutIn.additionalOutputInput(n, 34);
            G_dailyRemainingUse[n] = additionalOutIn.additionalOutputInput(n, 35);
            G_fGloLake[n] = additionalOutIn.additionalOutputInput(n, 36);
            G_PrevTotalUnstatisfiedUse[n] = additionalOutIn.additionalOutputInput(n, 37);
            G_dailySatisAllocatedUseInSecondCell[n] = additionalOutIn.additionalOutputInput(n, 40);
            G_dailyAllocatedUse[n] = additionalOutIn.additionalOutputInput(n, 41);
            G_unsatUseRiparian[n] = additionalOutIn.additionalOutputInput(n, 42);
            G_withdrawalIrrigFromSwb[n] = additionalOutIn.additionalOutputInput(n, 25) ;
            G_consumptiveUseIrrigFromSwb[n] = additionalOutIn.additionalOutputInput(n, 26) ;

        }

    }
    // Only for runs which consider glaciers
    if (options.glacierOpt == 1) {
        for (int n = 0; n < ng; n++) {

            if (glacierYear.G_glacier_mass_change_d365(n,day - 1) > glacierYear.G_precipitation_on_glacier_d365(n,day - 1)) {
                // If input glacier mass change is larger than input precipitation on glacier area,
                // increase precipitation on glacier area by the difference (necessary to avoid
                // negative values when calculating glacier runoff and to close the water budget)
                G_glacPrecip[n] = glacierYear.G_precipitation_on_glacier_d365(n,day - 1) -
                               (glacierYear.G_precipitation_on_glacier_d365(n,day - 1) - glacierYear.G_glacier_mass_change_d365(n,day - 1));
            } else
                G_glacPrecip[n] = glacierYear.G_precipitation_on_glacier_d365(n,day - 1);

            if (G_glacAdaptedArea[n] > 0.) { // for glacierized cells only
                G_glacRunoff[n] = G_glacPrecip[n] - glacierYear.G_glacier_mass_change_d365(n,day - 1); // calculate glacier runoff
            } else
                G_glacRunoff[n] = 0.;

            G_glacierStorage[n] += glacierYear.G_glacier_mass_change_d365(n,day - 1);  // G_glacierStorage is equal to cumulated input glacier mass change

            G_monthlyGlacierArea(n,month) = G_glacAdaptedArea[n];
            G_monthlyGlacierAreaFrac(n,month) = G_glacAreaFrac[n];
            G_monthlyGlacierRunoff(n,month) += G_glacRunoff[n];
            G_monthlyGlacierPrecip(n,month) += G_glacPrecip[n];
        } // end loop over grid cells
    } //end if glacierOpt


    for (int n = 0; n < ng; n++) {
        cellArea = geo.areaOfCellByArrayPos(n);

        // added for cell AET calculation (2.1f)
        // set daily variables to zero
        G_dailyLocLakeEvapo[n] = 0.;
        G_dailyGloLakeEvapo[n] = 0.;
        G_dailyResEvapo[n] = 0.;
        G_dailyLocWetlEvapo[n] = 0.;
        G_dailyGloWetlEvapo[n] = 0.;
        G_dailyRiverEvapo[n] = 0.;
        // added for new routing algorithm
        G_riverInflow[n] = 0.;
        G_daily_res_in[n] = 0.;
        G_daily_res_out[n] = 0.;
        G_locwetlextent[n] = 0.;
        G_glowetlextent[n] = 0.;
        G_transportedVolume_for_AET[n] = 0.;


        // WG22b (option 'use_alloc'): MaxStorage of local lakes and global lakes
        // required to identify the second cell with the maximum 'storageSum'
        if (options.use_alloc > 0) {

            G_locLakeMaxStorage[n] = 0.;
            G_gloLakeMaxStorage[n] = 0.;

            if (G_loc_lake[n] > 0.) {
                // Cell-specific active lake depth - Use array (initialized from cell-specific parameter)
                G_locLakeMaxStorage[n] = ((G_loc_lake[n]) / 100.) * cellArea * G_lakeDepthActive[n];
            }

            if (G_lake_area[n] > 0.) {
                // Cell-specific active lake depth - Use array (initialized from cell-specific parameter)
                G_gloLakeMaxStorage[n] = ((double) G_lake_area[n]) * G_lakeDepthActive[n];
            }
        }
        //	end if (options.use_alloc > 0)
    }
    // end loop over grid cells (initialize)


    double daily_totalUnsatisfiedUse = 0.;

    double remainingUse = 0., dailyUse = 0., totalDesiredUse = 0., totalNeighbourStorage = 0., totalDesiredUseCell = 0.;
    double daily_allocatedUseNextDay = 0., daily_UnsatAllocUseNextDay = 0., dailyActualUse = 0., demandWithoutAllocUse = 0., unsatAllocUse = 0., unsatUse = 0.;
    double unsatGLWDUseOfRiparians = 0.; // unsatisfied NAs which steem of riparian cells and is reallocated to riparian (not outflow) cells of global lakes or reservoirs
    //introduced in WG2.2b to implement NUs into storage equations:
    double PETgwrRemUseMax = 0., PETgwrRemUse = 0., RiverEvapoRemUse = 0., RivEvapoRemUseMax = 0.;
    //introduced for the new approach in WG2.2b to abstract 50% of NUs from each SWB, if both SWB exist:
    double remainingUseGloLake = 0., remainingUseGloLakeStart = 0., remainingUseRes = 0., remainingUseResStart = 0., remainingUseRiverStart = 0.;


#pragma omp parallel for  \
            private(dailyUse, totalDesiredUse, remainingUse) \
            shared(G_toBeCalculated, options)

    // loop over timeStepsPerDay (currently fixed to 1)
    for (short st = 1; st <= timeStepsPerDay; st++) {

        /* routing order from cell to cell: no parallelization allowed */

        // routingCell: (WG2.2) added for ordered routing scheme
        for (int routingCell = 0; routingCell < ng; routingCell++) {

            // added for new routing order in Wg2.2 // moved some lines up
            int n = (G_routingCell[routingCell] - 1);
            // Efficiency enhancement through local copies instead of often used method getValue
            double P_GWOUTF_C__at__n_routing = calParam.getValue(P_GWOUTF_C, n);
            double M_EVAREDEX__at__n_routing = calParam.getValue(M_EVAREDEX, n);

            if (options.use_alloc > 0)
                G_CellPositionRoutOrder[n] = routingCell;


            // initialize values
            remainingUse = 0.;
            dailyUse = 0.;
            totalDesiredUse = 0.;
            totalNeighbourStorage = 0.;
            totalDesiredUseCell = 0.;
            PETgwrRemUseMax = 0.;
            PETgwrRemUse = 0.;
            RiverEvapoRemUse = 0.;
            RivEvapoRemUseMax = 0.;
            remainingUseGloLake = 0.;
            remainingUseGloLakeStart = 0.;
            remainingUseRes = 0.;
            remainingUseResStart = 0.;
            remainingUseRiverStart = 0.;
            dailyActualUse = 0.;
            daily_allocatedUseNextDay = 0.;
            daily_UnsatAllocUseNextDay = 0.;
            demandWithoutAllocUse = 0.;
            unsatAllocUse = 0.;
            unsatUse = 0.;

            daily_totalUnsatisfiedUse = 0.;

            //for water temp calculation
            outflowlocLake = 0.;
            outflowlocWetl = 0.;
            outflowGloLake = 0.;
            outflowGloWetland = 0.;

            // Only execute for continental grid cells
            if (geo.G_contcell[n]) {

                cellArea = geo.areaOfCellByArrayPos(n);


                // If landAreaFraction == 0., the remaining canopy, snow, and soil storage from the previous timestep
                // becomes surface runoff (see daily.cpp).
                if (G_landAreaFrac[n] <= 0.) {
                    dailyLocalSurfaceRunoff = dailyWaterBalance.G_dailyStorageTransfer[n] * cellArea / 1000000. *
                                              G_landAreaFracPrevTimestep[n] / 100.;
                } else {
                    dailyLocalSurfaceRunoff = dailyWaterBalance.G_dailyLocalSurfaceRunoff[n] * cellArea / 1000000. *
                                              G_landAreaFrac[n] / 100.;
                }

                // If GWR under SWB / fractionalRouting is allowed,
                // GW storage is computed before local lakes in humid cells (baseflow directed into local lake)
                // and before the river in (semi-)arid cells (baseflow directed into river)

                // Calculate fswb*Rs and (a-fswb)*Rs for arid cells for Options fractionalRouting and gwr_swb:
                if (options.fractionalRoutingOpt && G_LDD[n] >= 0){
                    G_fswb_catchment[n] = G_fswbInit[n] *
                                          20.;    // G_fswb_catchment[n] = G_fswb[n] * 5.; // add a catchment area to surface water bodies (specified as x times the size of swb)...

                    if (G_fswb_catchment[n] > 1.) // ... and prevent that G_fswb is greater then 1
                        G_fswb_catchment[n] = 1.;

                    // Firstly, amount of local runoff flowing directly into the river is calculated, then local gw runoff is added for the rest term.
                    G_localRunoffIntoRiver[n] = (1. - G_fswb_catchment[n]) *
                                                dailyLocalSurfaceRunoff;    // dailyLocalSurfaceRunoff in km3
                }

                if ((1 == options.fractionalRoutingOpt) && (1 == G_aindex[n]) && (G_LDD[n] >= 0)) {    // (semi-)arid cells except inland sinks; fractionalRoutingOpt automatically set to 1 if aridareaOpt == 1


                    G_localGWRunoff[n] = 0.; // in arid regions, no GW flow into surface water bodies

                    G_localRunoff[n] = G_fswb_catchment[n] * dailyLocalSurfaceRunoff;

                    // Only for runs which consider glaciers
                    if ((options.glacierOpt == 1) && (G_glacRunoff[n] > 0.)) { // for glacierized cells only
                        // Add a fraction of glacier runoff to local runoff flowing directly into the river and route the rest through SWB
                        G_localRunoffIntoRiver[n] += (1. - G_fswb_catchment[n]) * G_glacRunoff[n];
                        G_localRunoff[n] += G_fswb_catchment[n] * G_glacRunoff[n];
                    }

                    if (options.aridareaOpt == 0){
                        G_groundwaterStoragePrevStep = G_groundwaterStorage[n];
                        netGWin = dailyWaterBalance.G_dailyGwRecharge[n] * cellArea * (G_landAreaFrac[n] / 100.) /
                                  1000000.;        // mm -> km3
                        // G_dailydailyNUg is adapted in this timestep based on the remaining use of last time step
                        if (options.subtract_use > 0) {
                            if (G_dailyRemainingUse[n] != 0.) {
                                G_dailydailyNUg[n] = updateNetAbstractionGW(n);
                            }
                            netGWin -= G_dailydailyNUg[n];
                        }

                        // Analytical solution of dS/dt = netGWR - k*S
                        // k_g in [1/d]
                        // Cell-specific calibration parameters - Use parameter  // FP

                        G_groundwaterStorage[n] = G_groundwaterStoragePrevStep * exp(-1. * P_GWOUTF_C__at__n_routing) +
                                                  (1. / P_GWOUTF_C__at__n_routing) * netGWin *
                                                  (1. - exp(-1. * P_GWOUTF_C__at__n_routing));

                        if(abs(G_groundwaterStorage[n])<=minStorVol)    //to counter numerical inaccuracies
                            G_groundwaterStorage[n]=0.;

                        groundwater_runoff_km3 = G_groundwaterStoragePrevStep - G_groundwaterStorage[n] + netGWin;
                        if (groundwater_runoff_km3 <=
                            0.) {    // different differential equation (ODE) applies: dS/dt = GWR - NAg (without k*S) -> S(t) = S(t-1) + netGWin

                            groundwater_runoff_km3 = 0.;

                            G_groundwaterStorage[n] = G_groundwaterStoragePrevStep + netGWin;

                            if(abs(G_groundwaterStorage[n])<=minStorVol)    //to counter numerical inaccuracies
                                G_groundwaterStorage[n]=0.;

                        }
                        if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store))
                            G_monthlyGwRunoff(n,month) += groundwater_runoff_km3;
                        G_localGWRunoffIntoRiver[n] = groundwater_runoff_km3;
                        potCellRunoff = (G_localRunoff[n] + G_localRunoffIntoRiver[n] +
                                         G_localGWRunoffIntoRiver[n]    // is in km3
                                         + ((dailyWaterBalance.G_lakeBalance[n]    // is in mm
                                             / 1000000.0)    // mm -> km3
                                            * (((double) G_lake_area[n] + (double) G_reservoir_area[n])    // km2
                                               + ((geo.areaOfCellByArrayPos(n) / 100.0)    // % -> km2
                                                  * (G_loc_lake[n] + G_loc_wetland[n] + G_glo_wetland[n])))));    // [%]

                        G_potCellRunoff[n] += potCellRunoff;
                        if (2 == options.grid_store) {
                            G_monthlyPotCellRunoff(n,month) += potCellRunoff;
                        }

                    }
                }

                // Calculate baseflow, runoff into local lakes and direct runoff into river for humid cells:
                if ((1 == options.fractionalRoutingOpt) && (0 == G_aindex[n]) && (G_LDD[n] >= 0)) {    // humid cells only (except for inland sinks)

                    G_groundwaterStoragePrevStep = G_groundwaterStorage[n];

                    // WG22b: ODE (ordinary differential equation) is applied for groundwater storage:
                    // dS/dt = GWR - NAg - k*S  (k*S = groundwater_runoff_km3)

                    // humid cells -> gwr_swb = 0
                    netGWin = (double) dailyWaterBalance.G_dailyGwRecharge[n] * cellArea *
                              ((double) G_landAreaFrac[n] / 100.) / 1000000.;        // mm -> km3

                    // G_dailydailyNUg is adapted in this time step based on the remaining use of last time step
                    if (options.subtract_use > 0) {
                        if (G_dailyRemainingUse[n] != 0.) {
                            G_dailydailyNUg[n] = updateNetAbstractionGW(n);
                        }
                        netGWin -= G_dailydailyNUg[n];
                    }

                    // Analytical solution of dS/dt = netGWR - k*S
                    // k_g in [1/d]
                    // Cell-specific calibration parameters - Use parameter  // FP

                    G_groundwaterStorage[n] = G_groundwaterStoragePrevStep * exp(-1. * P_GWOUTF_C__at__n_routing) +
                                              (1. / P_GWOUTF_C__at__n_routing) * netGWin *
                                              (1. - exp(-1. * P_GWOUTF_C__at__n_routing));

                    if(abs(G_groundwaterStorage[n])<=minStorVol)    //to counter numerical inaccuracies
                        G_groundwaterStorage[n]=0.;

                    groundwater_runoff_km3 = G_groundwaterStoragePrevStep - G_groundwaterStorage[n] + netGWin;

                    if (groundwater_runoff_km3 <=
                        0.) {    // different differential equation (ODE) applies:
                        // dS/dt = GWR - NAg (without k*S) -> S(t) = S(t-1) + netGWin

                        groundwater_runoff_km3 = 0.;

                        G_groundwaterStorage[n] = G_groundwaterStoragePrevStep + netGWin;

                        if(abs(G_groundwaterStorage[n])<=minStorVol)    //to counter numerical inaccuracies
                            G_groundwaterStorage[n]=0.;

                    }
                    if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store))
                        G_monthlyGwRunoff(n,month) += groundwater_runoff_km3;

                    // in humid areas, groundwater is also routed through surface water bodies
                    G_localGWRunoffIntoRiver[n] = (1. - G_fswb_catchment[n]) * groundwater_runoff_km3;

                    // route groundwater through surface water bodies in humid regions
                    G_localGWRunoff[n] = G_fswb_catchment[n] * groundwater_runoff_km3;


                    G_localRunoff[n] = (G_fswb_catchment[n] * dailyLocalSurfaceRunoff) + G_localGWRunoff[n];

                    // Only for runs which consider glaciers
                    if ((options.glacierOpt == 1) && (G_glacRunoff[n] > 0.)) { // for glacierized cells only
                        // Add a fraction of glacier runoff to local runoff flowing directly into the river and route the rest through SWB
                        G_localRunoffIntoRiver[n] += (1. - G_fswb_catchment[n]) * G_glacRunoff[n];
                        G_localRunoff[n] += G_fswb_catchment[n] * G_glacRunoff[n];
                    }


                    // G_localRGWRunoffIntoRiver
                    // reintroduced for CFA calculation
                    potCellRunoff = (G_localRunoff[n] + G_localRunoffIntoRiver[n] +
                                     G_localGWRunoffIntoRiver[n]    // is in km3
                                     + ((dailyWaterBalance.G_lakeBalance[n]    // is in mm
                                         / 1000000.0)    // mm -> km3
                                        * (((double) G_lake_area[n] + (double) G_reservoir_area[n])    // km2
                                           + ((geo.areaOfCellByArrayPos(n) / 100.0)    // % -> km2
                                              * (G_loc_lake[n] + G_loc_wetland[n] + G_glo_wetland[n])))));    // [%]

                    G_potCellRunoff[n] += potCellRunoff;

                    if (2 == options.grid_store) {
                        G_monthlyPotCellRunoff(n,month) += potCellRunoff;
                    }

                }    // end of computation of GW storage and baseflow for options.aridareaOpt == 1 in humid cells

                if ((0 == options.fractionalRoutingOpt) && (0 == options.aridareaOpt) && (G_LDD[n] >= 0)) {


                    G_localRunoff[n] = dailyLocalSurfaceRunoff;
                    // Only for runs which consider glaciers
                    if ((options.glacierOpt == 1) && (G_glacRunoff[n] > 0.)) { // for glacierized cells only
                        // Add a fraction of glacier runoff to local runoff flowing directly into the river and route the rest through SWB
                        G_localRunoff[n] += G_glacRunoff[n];
                    }

                    G_groundwaterStoragePrevStep = G_groundwaterStorage[n];
                    netGWin = dailyWaterBalance.G_dailyGwRecharge[n] * cellArea * (G_landAreaFrac[n] / 100.) /
                                   1000000.;        // mm -> km3

                                   // G_dailydailyNUg is adapted in this timestep based on the remaining use of last time step
                                   if (options.subtract_use > 0) {
                                       if (G_dailyRemainingUse[n] != 0) {
                                           G_dailydailyNUg[n] = updateNetAbstractionGW(n);
                                       }
                                       netGWin -= G_dailydailyNUg[n];
                                   }

                                   // Analytical solution of dS/dt = netGWR - k*S
                                   // k_g in [1/d]
                                   // Cell-specific calibration parameters - Use parameter  // FP

                                   G_groundwaterStorage[n] = G_groundwaterStoragePrevStep * exp(-1. * P_GWOUTF_C__at__n_routing) +
                                                  (1. / P_GWOUTF_C__at__n_routing) * netGWin *
                                                  (1. - exp(-1. * P_GWOUTF_C__at__n_routing));

                    if(abs(G_groundwaterStorage[n])<=minStorVol)    //to counter numerical inaccuracies
                        G_groundwaterStorage[n]=0.;

                                   groundwater_runoff_km3 = G_groundwaterStoragePrevStep - G_groundwaterStorage[n] + netGWin;
                                   if (groundwater_runoff_km3 <=
                                   0.) {    // different differential equation (ODE) applies: dS/dt = GWR - NAg (without k*S) -> S(t) = S(t-1) + netGWin

                                       groundwater_runoff_km3 = 0.;

                                       G_groundwaterStorage[n] = G_groundwaterStoragePrevStep + netGWin;

                                       if(abs(G_groundwaterStorage[n])<=minStorVol)    //to counter numerical inaccuracies
                                           G_groundwaterStorage[n]=0.;

                                   }
                                   if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store))
                                       G_monthlyGwRunoff(n,month) += groundwater_runoff_km3;
                                   G_localGWRunoff[n] = groundwater_runoff_km3;
                                   G_localRunoff[n] += G_localGWRunoff[n];
                                   potCellRunoff = (G_localRunoff[n] + G_localGWRunoff[n]    // is in km3
                                                  + ((dailyWaterBalance.G_lakeBalance[n]    // is in mm
                                                  / 1000000.0)    // mm -> km3
                                                  * (((double) G_lake_area[n] + (double) G_reservoir_area[n])    // km2
                                                  + ((geo.areaOfCellByArrayPos(n) / 100.0)    // % -> km2
                                                  * (G_loc_lake[n] + G_loc_wetland[n] + G_glo_wetland[n])))));    // [%]

                                                  G_potCellRunoff[n] += potCellRunoff;
                                                  if (2 == options.grid_store) {
                                                      G_monthlyPotCellRunoff(n,month) += potCellRunoff;
                                                  }

                }
                if (G_LDD[n] < 0) { // in all inland sinks

                    G_groundwaterStoragePrevStep = G_groundwaterStorage[n];

                    // WG22b: ODE (ordinary differential equation) is applied for groundwater storage: dS/dt = GWR - NAg - k*S  (k*S = groundwater_runoff_km3)

                    netGWin = dailyWaterBalance.G_dailyGwRecharge[n] * cellArea * (G_landAreaFrac[n] / 100.) /
                              1000000.;        // mm -> km3
                    // G_dailydailyNUg is adapted in this timestep based on the remaining use of last time step
                    if (options.subtract_use > 0) {
                        if (G_dailyRemainingUse[n] != 0.) {
                            G_dailydailyNUg[n] = updateNetAbstractionGW(n);
                        }
                        netGWin -= G_dailydailyNUg[n];
                    }

                    // Analytical solution of dS/dt = netGWR - k*S
                    // k_g in [1/d]
                    // Cell-specific calibration parameters - Use parameter  // FP

                    G_groundwaterStorage[n] = G_groundwaterStoragePrevStep * exp(-1. * P_GWOUTF_C__at__n_routing) +
                                              (1. / P_GWOUTF_C__at__n_routing) * netGWin *
                                              (1. - exp(-1. * P_GWOUTF_C__at__n_routing));

                    if(abs(G_groundwaterStorage[n])<=minStorVol)    //to counter numerical inaccuracies
                        G_groundwaterStorage[n]=0.;

                    groundwater_runoff_km3 = G_groundwaterStoragePrevStep - G_groundwaterStorage[n] + netGWin;

                    if (groundwater_runoff_km3 <=
                        0.) {    // different differential equation (ODE) applies: dS/dt = GWR - NAg (without k*S) -> S(t) = S(t-1) + netGWin

                        groundwater_runoff_km3 = 0.;

                        G_groundwaterStorage[n] = G_groundwaterStoragePrevStep + netGWin;

                        if(abs(G_groundwaterStorage[n])<=minStorVol)    //to counter numerical inaccuracies
                            G_groundwaterStorage[n]=0.;

                    }
                    if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store))
                        G_monthlyGwRunoff(n,month) += groundwater_runoff_km3;

                    // CR 2015-08: If landAreaFraction == 0., the remaining canopy, snow, and soil storage from the previous timestep becomes surface runoff (see daily.cpp).
                    if (G_landAreaFrac[n] == 0.) {
                        dailyLocalSurfaceRunoff = dailyWaterBalance.G_dailyStorageTransfer[n] * cellArea / 1000000. *
                                                  G_landAreaFracPrevTimestep[n] / 100.;
                    } else {
                        dailyLocalSurfaceRunoff = dailyWaterBalance.G_dailyLocalSurfaceRunoff[n] * cellArea / 1000000. *
                                                  G_landAreaFrac[n] / 100.;
                    }

                    G_localGWRunoff[n] = groundwater_runoff_km3;
                    G_localRunoff[n] = dailyLocalSurfaceRunoff + G_localGWRunoff[n];
                    if ((options.glacierOpt == 1) && (G_glacRunoff[n] > 0.)) // for glacierized cells only
                        G_localRunoff[n] += G_glacRunoff[n]; // add glacier runoff to local runoff

                    potCellRunoff = (G_localRunoff[n]     // is in km3
                                     + ((dailyWaterBalance.G_lakeBalance[n]    // is in mm
                                         / 1000000.0)    // mm -> km3
                                        * (((double) G_lake_area[n] + (double) G_reservoir_area[n])    // km2
                                           + ((geo.areaOfCellByArrayPos(n) / 100.0)    // % -> km2
                                              * (G_loc_lake[n] + G_loc_wetland[n] + G_glo_wetland[n])))));    // [%]

                    G_potCellRunoff[n] += potCellRunoff;
                    if (2 == options.grid_store) {
                        G_monthlyPotCellRunoff(n, month) += potCellRunoff;
                    }
                }

                if (options.subtract_use > 0) {

                    if (1 == G_toBeCalculated[n]) {

                        // estimation of dailyUse and totalDesiredUseCell
                        // option for satisfaction of NUs of riparian cells of global lakes and reservoirs
                        if (options.aggrNUsGloLakResOpt == 1) {

                            if (G_reservoir_area[n] > 0. || G_lake_area[n] > 0.) {    // outflow cell of global lake and/or reservoir

                                dailyUse = G_dailydailyNUsAggregated[n];    // [km3/day]

                                // totalDesiredUseCell (introduced for option 'aggrNUsGloLakResOpt'):
                                // required to compute G_satisfiedUse and G_actualUse of outflow and riparian cells of global lakes/reservoirs
                                // based on the initial G_dailydailyNUs[n] (not aggregated over riparian cells)
                                totalDesiredUseCell = G_dailydailyNUs[n];

                            } else {    // no outflow cell (either riparian cell or other cell)

                                if (G_glwdUnit[n] > 0) {    // riparian cell, but no outflow cell

                                    // G_unsatUseRiparian[n] is set to 0 for outflow cells (see computation of water balance for global lakes and reservoirs)
                                    // For outflow cells,
                                    dailyUse = G_dailydailyNUsAggregated[n] + G_unsatUseRiparian[n];
                                    totalDesiredUseCell = G_dailydailyNUs[n] +
                                                          G_unsatUseRiparian[n]; // only required to compute G_satisfiedUse and G_actualUse for option 'aggrNUsGloLakResOpt' (see above)
                                } else {    // neither outflow nor riparian cell
                                    dailyUse = G_dailydailyNUs[n];
                                    totalDesiredUseCell = G_dailydailyNUs[n];  // only required to redistribute remainingUse over riparian cells for option 'aggrNUsGloLakResOpt' (see above)
                                }
                            }

                        } else {

                            // G_dailydailyNUs (surface water abstractions) can be negative (i.e. return flow to surface water bodies)
                            dailyUse = G_dailydailyNUs[n]; // [km3/day]
                            totalDesiredUseCell = G_dailydailyNUs[n]; // only required to redistribute remainingUse over riparian cells for option 'aggrNUsGloLakResOpt' (see above)

                        }

                        totalDesiredUse = dailyUse;


                        if (options.delayedUseSatisfaction == 1) {
                            totalDesiredUse += G_totalUnsatisfiedUse[n];
                            totalDesiredUseCell += G_totalUnsatisfiedUse[n];// only required to redistribute remainingUse over riparian cells for option 'aggrNUsGloLakResOpt' (see above)
                            // G_PrevTotalUnstatisifiedUse[n] is needed for comparing it to G_totalUnstatisfiedUse to
                            // adapt NAg based on this difference
                            G_PrevTotalUnstatisfiedUse[n] = G_totalUnsatisfiedUse[n];
                            G_totalUnsatisfiedUse[n] = 0;

                        }

                        if (options.use_alloc > 0) {

                            // New approach in WG22b:
                            // If a second cell exists, all unsatisfied use of cell n is added to the remainingUse
                            // in the second cell (G_AllocatedUse[secondCell]) and is satisfied in the same
                            // or the next time step, depending on the routing order.
                            // If there is no second cell and "options.delayedUseSatisfaction == 1":
                            // The remainingUse is satisfied in the next timestep (G_totalUnsatisfiedUse[n]).
                            // Thus, G_AllocatedUse[n] and G_totalUnsatisfiedUse[n] are included in the storage
                            // equations of surface water bodies like dailydailyNUs.
                            // G_UnsatAllocUse[n] is the fraction of allocated use that could not be satisfied in a nbc
                            // in the last timestep. Thus it is reallocated to cell n.

                            totalDesiredUse += G_UnsatAllocUse[n];
                            totalDesiredUseCell += G_UnsatAllocUse[n];
                            // G_PrevUnsatAllocUse is needed for comparing it to G_UnsatAllocUse[n] to adapt NAg based
                            // on this difference
                            G_PrevUnsatAllocUse[n] = G_UnsatAllocUse[n];
                            G_UnsatAllocUse[n] = 0.;

                            demandWithoutAllocUse = totalDesiredUse;

                            totalDesiredUse += G_AllocatedUse[n];
                            totalDesiredUseCell += G_AllocatedUse[n]; // only required to compute G_satisfiedUse and G_actualUse (see above)

                            // If n is identified as a "second cell" again during this time step,
                            // G_AllocatedUse[n] will be satisfied in the next timestep.
                            // G_dailyAllocatedUse[n] only used as output
                            G_dailyAllocatedUse[n] = G_AllocatedUse[n];
                            G_AllocatedUse[n] = 0.;
                            // daily_allocatedUseNextDay only used as output
                            daily_allocatedUseNextDay = G_daily_allocatedUseNextDay[n];
                            G_daily_allocatedUseNextDay[n] = 0.;
                            // daily_UnsatAllocUseNextDay only used as output
                            daily_UnsatAllocUseNextDay = G_daily_UnsatAllocUseNextDay[n];
                            G_daily_UnsatAllocUseNextDay[n] = 0.;

                        }
                        // end if (options.use_alloc > 0)

                        // 'remainingUse' contains always the amount of water that
                        // has not been satisfied (if "options.use_alloc" > 0 and/or "options.delayedUseSatisfaction" > 0)
                        // At this point, remainingUse can get negative!
                        // G_dailydailyNUs (surface water abstractions) can be negative (i.e. return flow to surface water bodies)
                        remainingUse = totalDesiredUse;
                    }
                    // end if (1 == G_toBeCalculated[n])

                } else {    // if (options.subtract_use == 0)
                    remainingUse = 0.;
                }

                // added for cell AET calculation (2.2)
                // added for new routing algorithm WG2.2
                // store value for surface water storages in previous time step
                double G_riverStoragePrevStep = 0.;
                double G_locLakeStoragePrevStep = 0.;
                double G_locWetlStoragePrevStep = 0.;
                double G_gloLakeStoragePrevStep = 0.;
                double G_gloResStoragePrevStep = 0.;
                double G_gloWetlStoragePrevStep = 0.;

                if (0 != G_toBeCalculated[n]) {

                    // G_localRunoffPerStep contributes to the inflow of lake/wetland system
                    inflow = G_localRunoff[n];
                    cellArea = geo.areaOfCellByArrayPos(n);

                    // Efficiency enhancement through local copies instead of often used method getValue // FP
                    double P_SWOUTF_C__at__n = calParam.getValue(P_SWOUTF_C, n);

                    // -------------- Routing through local lake ----------------------------------------------------
                    if (G_loc_lake[n] > 0.) {

                        // assign correct storage value which was calculated in the last timestep or maxStorage at model day one
                        // Note: G_locLakeStorage[n] at this point can differ from the output written at the end of the last time step,
                        // since water abstractions from a second cell (Option "use_alloc") are carried out after writing grids.
                        G_locLakeStoragePrevStep = G_locLakeStorage[n];

                        // maximum storage capacity
                        // Cell-specific active lake depth - Use array (initialized from cell-specific parameter)  // FP
                        maxStorage = ((G_loc_lake[n]) / 100.) * cellArea * G_lakeDepthActive[n]; // [km3]
                        G_locLakeMaxStorage[n] = maxStorage; // G_locLakeMaxStorage is used when NUs is abstracted from locLakeStorage at the end of the routing loop (in cell n and in second cell)

                        // This factor was added to reduce PET as a function of actual lake storage (2.1f).
                        // Without reduction PET would lead to a continuous decline of lake level in some cases ((semi)arid regions)
                        // use reduction factor from previous time step (or the initial value)
                        locLakeAreaReductionFactor = G_locLakeAreaReductionFactor[n];


                        //reintroduced CFA

                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        // CFA approach for surface water bodies from version 2.1f
                        //
                        // Water balance of lake / reservoir:
                        // dS/dt = CFA*(P - PET) + Qin - NAs - gwr_swb
                        //
                        // (P - PET) is considered as the runoff from the respective SWB. Therefore, only (P - PET) is multiplied by CFA
                        // and only P and PET are taken into account to compute PET_corrected.
                        // (different from soil water balance, where R equals (P - AET - dS/dt), and dS/dt is taken into account to compute PET_corrected.
                        //
                        // SWB: CFA affects R, AET and dS/dt.
                        // Land water balance (daily.cpp): CFA affects only R and AET.
                        // (as in version 2.1f, but code implementation changed in 2.2b)
                        //
                        // Approach to compute PET_corrected:
                        //
                        //	equation 1) P - PET = R
                        // 	equation 2) P - PET_corr = R * CFA
                        // 						   R = (P - PET_corr)/CFA
                        //
                        // 	equation 2) in equation 1):
                        // 	P - PET = (P - PET_corr)/CFA			| *CFA
                        // 	(P - PET) * CFA = P - PET_corr
                        //
                        // 	CFA*P - CFA * PET = P - PET_corr
                        // 	PET_corr = P - CFA*P + CFA * PET
                        //
                        // Final equation in routing.cpp:
                        //
                        // 	PET_corr = (1 - CFA) * P + CFA * PET
                        //
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        //
                        // compare daily.cpp, line 1217:
                        // 	AET_corr = dS/dt*(CFA - 1) - P*(CFA - 1) + CFA * AET_uncorr
                        //
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

                        G_dailyLocLakeEvapo[n] =
                                       ((1.0 - dailyWaterBalance.G_cellCorrFact[n])                        // [mm]
                                        * dailyWaterBalance.G_openWaterPrec[n] * locLakeAreaReductionFactor)
                                       + (dailyWaterBalance.G_cellCorrFact[n]
                                          * (dailyWaterBalance.G_openWaterPET[n] * locLakeAreaReductionFactor));

                        if (G_dailyLocLakeEvapo[n] < 0.) { // to avoid unphysical negative aet
                            G_dailyLocLakeEvapo[n] = 0.;
                        }


                        totalInflow = inflow + (dailyWaterBalance.G_openWaterPrec[n] * locLakeAreaReductionFactor)
                                               * (cellArea / 1000000.) * (G_loc_lake[n] / 100.);    //[mm]->[km3]


                        // gwr below surface water bodies (except in arid inland sinks)
                        if ((1 == options.aridareaOpt) && (1 == G_aindex[n]) && (G_LDD[n] >= 0)){
                            gwr_loclak[n] = Kswbgw * locLakeAreaReductionFactor * G_loc_lake[n] / 100. /
                                            (geo.G_contfreq[n] /
                                             100.); // gwr below local lakes in mm/d with respect to continental area

                        } else
                            gwr_loclak[n] = 0.;    // -> to avoid repeating the following equations for options.aridareaOpt==0

                        PETgwr = G_dailyLocLakeEvapo[n] * (cellArea / 1000000.) *
                                 (G_loc_lake[n] / 100.)                            //[mm]->[km3]
                                 + gwr_loclak[n] * cellArea * (geo.G_contfreq[n] / 100.) / 1000000.;    //[mm -> km3]

                        // Maximum amount of (openwaterPET + gwr_loclak) to ensure that G_locLakeStorage is always larger than -maxStorage:
                        PETgwrMax = G_locLakeStoragePrevStep + maxStorage + totalInflow;

                        if (PETgwrMax < 0.) // avoids numeric problems when reading in and prevents that value gets negatively
                            PETgwrMax = 0.;

                        if (PETgwr > PETgwrMax) {    // reduction of PET and gwr_locLak required

                            G_locLakeStorage[n] = (-1.) * maxStorage;
                            gwr_loclak[n] *= PETgwrMax / PETgwr;    //[mm]
                            G_dailyLocLakeEvapo[n] *= PETgwrMax / PETgwr;    //[mm]

                        } else { // no reduction of PET and gwr_loclak required

                            G_locLakeStorage[n] = G_locLakeStoragePrevStep + totalInflow - PETgwr;
                        }

                        // WG22b: outflow is a function of storage from previous time step (outflow = f(St-1)). Until WG22: outflow = f(S(t))
                        if (G_locLakeStoragePrevStep > 0.) {
                            // Cell-specific calibration parameters - Use parameter  // FP
                            // loc_storageFactor (now: only one factor for surface water P_SWOUTF_C) in [1/d]
                            outflow = P_SWOUTF_C__at__n * G_locLakeStoragePrevStep
                                      * pow((G_locLakeStoragePrevStep / maxStorage), lakeOutflowExp);

                            if (G_locLakeStorage[n] <=
                                0.) {    // Even if S(t-1) is positive, G_locLakeStorage[n] can become negative, if PET+gwr > totalInflow.
                                outflow = 0;                    // In this case, the outflow is set to zero.
                            } else {
                                if (outflow > G_locLakeStorage[n])
                                    outflow = G_locLakeStorage[n];
                            }

                        } else // if G_locLakeStoragePrevStep <= 0, no outflow should occur
                            outflow = 0.;

                        G_locLakeStorage[n] -= outflow;

                        if(abs(G_locLakeStorage[n])<=minStorVol) //to counter numerical inaccuracies
                            G_locLakeStorage[n]=0.;

                        // reduce G_locLakeStorage[n] to maximum storage capacity
                        if (G_locLakeStorage[n] > maxStorage) {
                            outflow += (G_locLakeStorage[n] - maxStorage);
                            G_locLakeStorage[n] = maxStorage;
                        }

                        // If (PETgwr > PETgwrMax), PET is reduced by PETgwrRedFactor. If a reduction is not required, PETgwrRedFactor is set to 1.
                        // HMS 2014-09-26 if PETgwr is 0 (because gwr_swb is 0 because reductionfactor is 0 because storage is empty),
                        // division through 0 --> has to be catched, because otherwise G_CELL_RUNOFF_TOTAL and G_CELL_AET gets nan.
                        if (PETgwr == 0.)
                            PETgwrRedFactor = 1.;
                        else
                            PETgwrRedFactor = PETgwrMax / PETgwr;

                        if (PETgwrMax > PETgwr) {// -> reduction of PET not required
                            PETgwrRedFactor = 1.;
                        }

                        if (2 == options.grid_store)
                            // sum up the negative fraction of the potential cell runoff which has not been realized
                            // due to reduction of potential evaporation
                            G_monthlyPotCellRunoffDeficit(n,month)
                                           += (0.0 - dailyWaterBalance.G_openWaterPET[n] * PETgwrRedFactor
                                                     * (1.0 - locLakeAreaReductionFactor)
                                                     * dailyWaterBalance.G_cellCorrFact[n]
                                                     * (cellArea / 1000000.0)
                                                     * (G_loc_lake[n] / 100.0));


                        inflow = outflow;
                        outflowlocLake = outflow;

                        // Area reduction factor is now calculated for the next day.
                        // Cell-specific calibration parameters - Apply multiplier
                        G_locLakeAreaReductionFactor[n] = 1. -
                                                          pow(fabs(G_locLakeStorage[n] - maxStorage)
                                                              / (2. * maxStorage),
                                                              (M_EVAREDEX__at__n_routing * evapoReductionExp));

                        if (G_locLakeAreaReductionFactor[n] <
                            0.)
                            G_locLakeAreaReductionFactor[n] = 0.;

                        if (G_locLakeAreaReductionFactor[n] > 1.)  // (could occur due to numerical inaccuracies)
                            G_locLakeAreaReductionFactor[n] = 1.;

                    }
                    // outflow of local lakes is inflow to local wetlands


                    // transport through local wetlands
                    if (G_loc_wetland[n] > 0.) {

                        // assign correct storage value which was calculated in the last timestep or maxStorage at model day one
                        G_locWetlStoragePrevStep = G_locWetlStorage[n];

                        // maximum storage capacity
                        // Cell-specific active wetland depth - Use array (initialized from cell-specific parameter)
                        maxStorage = ((G_loc_wetland[n]) / 100.) * cellArea * G_wetlDepthActive[n];

                        // use reduction factor from previous time step (or the initial value)
                        locWetlAreaReductionFactor = G_locWetlAreaReductionFactor[n];

                        // G_dailyLocWetlEvapo[n]: see local lakes for explanation

                        G_dailyLocWetlEvapo[n] =
                                       ((1.0 - dailyWaterBalance.G_cellCorrFact[n])                        // [mm]
                                        * dailyWaterBalance.G_openWaterPrec[n] * locWetlAreaReductionFactor)
                                       + (dailyWaterBalance.G_cellCorrFact[n]
                                          * (dailyWaterBalance.G_openWaterPET[n] * locWetlAreaReductionFactor));

                        if (G_dailyLocWetlEvapo[n] < 0.) { // avoid unphysical negative aet values
                            G_dailyLocWetlEvapo[n] = 0.;
                        }


                        // calculate inflow only in dependence of reduction factor to avoid inconsistent counting of precipitation
                        totalInflow = inflow + (dailyWaterBalance.G_openWaterPrec[n] * locWetlAreaReductionFactor
                                                * (cellArea / 1000000.) * (G_loc_wetland[n] / 100.));

                        G_locwetlextent[n] = locWetlAreaReductionFactor * cellArea * (G_loc_wetland[n] / 100.);
                        // gwr below surface water bodies (except in arid inland sinks)
                        if ((1 == options.aridareaOpt) && (1 == G_aindex[n]) && (G_LDD[n] >= 0)) {
                            gwr_locwet[n] = Kswbgw * locWetlAreaReductionFactor * G_loc_wetland[n] / 100. /
                                            (geo.G_contfreq[n] /
                                             100.); // gwr below local wetlands in mm/d with respect to continental area

                        } else
                            gwr_locwet[n] = 0.; // -> to avoid repeating the following equations for options.aridareaOpt==0

                        PETgwr = G_dailyLocWetlEvapo[n] * (cellArea / 1000000.) *
                                 (G_loc_wetland[n] / 100.)                        //[mm -> km3]
                                 + gwr_locwet[n] * cellArea * (geo.G_contfreq[n] / 100.) / 1000000.;    //[mm -> km3]

                        // Maximum amount of (openwaterPET + gwr_locwet) to ensure that G_locWetlStorage does not fall below zero:
                        PETgwrMax = G_locWetlStoragePrevStep + totalInflow;

                        if (PETgwr > PETgwrMax) {

                            G_locWetlStorage[n] = 0.;
                            gwr_locwet[n] *= PETgwrMax / PETgwr;    // [mm]
                            G_dailyLocWetlEvapo[n] *= PETgwrMax / PETgwr;    // [mm]

                        } else {
                            G_locWetlStorage[n] = G_locWetlStoragePrevStep + totalInflow - PETgwr;
                        }

                        if(abs(G_locWetlStorage[n])<=minStorVol)    //to counter numerical inaccuracies
                            G_locWetlStorage[n]=0.;

                        if (G_locWetlStorage[n] > 0.) {
                            // Cell-specific calibration parameters - Use parameter  // FP
                            // loc_storageFactor (now: only one factor for surface water P_SWOUTF_C) in [1/d]
                            outflow = P_SWOUTF_C__at__n * G_locWetlStorage[n]
                                      * pow((G_locWetlStorage[n] / maxStorage), wetlOutflowExp);

                            if (outflow > G_locWetlStorage[n]) {    // this could happen if S(t-1) >> S(t)
                                outflow = G_locWetlStorage[n];
                            }

                        } else
                            outflow = 0.;

                        // substract outflow from storage
                        G_locWetlStorage[n] -= outflow;

                        // reduce G_locWetlStorage to maximum storage capacity

                        if (G_locWetlStorage[n] > maxStorage) {
                            outflow += (G_locWetlStorage[n] - maxStorage);
                            G_locWetlStorage[n] = maxStorage;
                        }

                        // If (PETgwr > PETgwrMax), PET is reduced by PETgwrRedFactor. If a reduction is not required, PETgwrRedFactor is set to 1.
                        // if PETgwr is 0 (because gwr_swb is 0 because reductionfactor is 0 because storage is empty),
                        // division through 0 --> has to be catched, because otherwise G_CELL_RUNOFF_TOTAL and G_CELL_AET gets nan.
                        if (PETgwr ==
                            0.)
                            PETgwrRedFactor = 1.;
                        else
                            PETgwrRedFactor = PETgwrMax / PETgwr;

                        if (PETgwrMax > PETgwr)    // -> reduction of PET not required
                            PETgwrRedFactor = 1.;

                        // reintroduced CFA
                        if (2 == options.grid_store)
                            // sum up the negative fraction of the potential cell runoff which has not been realized
                            // due to reduction of potential evaporation
                            G_monthlyPotCellRunoffDeficit(n,month)
                                           += (0.0 - dailyWaterBalance.G_openWaterPET[n] * PETgwrRedFactor
                                                     * (1.0 - locWetlAreaReductionFactor)
                                                     * dailyWaterBalance.G_cellCorrFact[n]
                                                     * (cellArea / 1000000.0)
                                                     * (G_loc_wetland[n] / 100.0));

                        inflow = outflow;
                        outflowlocWetl = outflow;

                        // Area reduction factor is now calculated for the next day.
                        // Cell-specific calibration parameters - Apply multiplier
                        G_locWetlAreaReductionFactor[n] = 1. -
                                                          pow(fabs(G_locWetlStorage[n] - maxStorage)
                                                              / (maxStorage),
                                                              (M_EVAREDEX__at__n_routing * evapoReductionExp));


                        if (G_locWetlAreaReductionFactor[n] < 0.)    // (could occur due to numerical inaccuracies)
                            G_locWetlAreaReductionFactor[n] = 0.;

                        if (G_locWetlAreaReductionFactor[n] > 1.)  // (could occur due to numerical inaccuracies)
                            G_locWetlAreaReductionFactor[n] = 1.;

                    }
                    // end if (G_loc_wetland[n] > 0.)

                    // water that comes out of the system of local lakes/wetlands
                    // is inflow into the river and is being routed through global lakes, reservoirs and global wetlands

                    inflow += G_riverInflow[n];
                    inflowUpstream = G_riverInflow[n];

                    // transport through global lakes
                    // input is defined as the inflow from the local wetlands in addition to the upstream flow

                    // new reservoir algorithm
                    if (G_lake_area[n] > 0.) //cell is a natural lake, or the reservoir type is unknown
                    {
                        //assign correct storage value which was calculated in the last timestep or maxStorage at model day one
                        G_gloLakeStoragePrevStep = G_gloLakeStorage[n];

                        // maximum storage capacity
                        // Cell-specific active lake depth - Use array (initialized from cell-specific parameter)  // FP
                        maxStorage = ((double) G_lake_area[n]) * G_lakeDepthActive[n];
                        G_gloLakeMaxStorage[n] = maxStorage; // CR: Introduced in WG22b: G_gloLakeMaxStorage is used when NUs is abstracted from G_gloLakeMaxStorage of a second (neighboring) cell.
                        //  maxStG_locWetlStorage[n]orage = (((double) G_glo_lake[n]) / 100.) * cellArea * lakeDepth;	// [km3]

                        // This factor was added to reduce PET as a function of actual lake storage (2.1f).
                        // Without reduction PET would lead to a continuous decline of lake level in some cases ((semi)arid regions)
                        // HMS use reduction factor from previous time step (or the initial value)
                        gloLakeEvapoReductionFactor = G_gloLakeEvapoReductionFactor[n];

                        // G_dailyGloLakeEvapo[n]: see local lakes for explanation

                        G_dailyGloLakeEvapo[n] =
                                       ((1.0 - dailyWaterBalance.G_cellCorrFact[n])                            // [mm]
                                        * dailyWaterBalance.G_openWaterPrec[n])
                                       + (dailyWaterBalance.G_cellCorrFact[n]
                                          * dailyWaterBalance.G_openWaterPET[n] * gloLakeEvapoReductionFactor);

                        if (G_dailyGloLakeEvapo[n] < 0.) { // avoid unphysical negative aet values
                            G_dailyGloLakeEvapo[n] =0.;
                        }


                        // reintroduced CFA in calculation (!no distinguishing for aridareaOpt!)
                        totalInflow = inflow + (dailyWaterBalance.G_openWaterPrec[n]
                                                * ((double) G_lake_area[n] / 1000000.));

                        G_GloLakePrecip[n] = dailyWaterBalance.G_openWaterPrec[n] * ((double) G_lake_area[n]) /
                                             1000000.; // km3

                        // gwr below surface water bodies (except in arid inland sinks)
                        if ((1 == options.aridareaOpt) && (1 == G_aindex[n]) && (G_LDD[n] >= 0)) {

                            gwr_glolak[n] = Kswbgw * gloLakeEvapoReductionFactor * ((double) G_lake_area[n] /
                                                                                    (cellArea * (geo.G_contfreq[n] /
                                                                                                 100.))); // gwr below global lakes in mm/d

                        } else {
                            gwr_glolak[n] = 0.;    // to avoid repeating the following equations for options.aridareaOpt==0
                        }


                        // remainingUse is evenly distributed between global reservoirs and global lakes (if existing) within a grid cell
                        // If remainingUse is negative (i.e. return flow to surface water bodies, from negative surface water abstractions G_dailydailyNUs)
                        // then it is added to the storage
                        if (G_reservoir_area[n] >
                            0.) {    // If the cell also contains a reservoir, 50% of NAs is withdrawn from each of the SWB. If "remainingUseGloLake" cannot be fully satisfied, it is added to the 50% NAs of the reservoir.
                            remainingUseGloLake = 0.5 * remainingUse;
                        } else {
                            remainingUseGloLake = remainingUse;
                        }

                        // Required in 22b to compute actual use for Option 28 "aggrNUsGloLakResOpt"
                        remainingUseGloLakeStart = remainingUseGloLake;


                        PETgwrRemUse = G_dailyGloLakeEvapo[n] * ((double) G_lake_area[n] /
                                                                 1000000.)                                    //[mm]->[km3]
                                       + gwr_glolak[n] * cellArea * (geo.G_contfreq[n] / 100.) / 1000000.    //[mm]->[km3]
                                       + remainingUseGloLake;

                        // Maximum amount of water available to satisfy PET, gwr and NAs,
                        // until St = -Smax (outflow = 0 if storage <= 0).
                        PETgwrRemUseMax = totalInflow + maxStorage + G_gloLakeStoragePrevStep;

                        if (PETgwrRemUse > PETgwrRemUseMax) {
                            // PETgwrRemUse could be zero -> division through zero.

                            G_gloLakeStorage[n] = (-1.) * maxStorage;
                            outflow = 0.;
                            G_dailyGloLakeEvapo[n] *=
                                           PETgwrRemUseMax / PETgwrRemUse;                                // [mm]
                            gwr_glolak[n] *= PETgwrRemUseMax / PETgwrRemUse;                                // [mm]

                            if (remainingUseGloLake > 0.) {

                                remainingUseGloLake -= remainingUseGloLake * PETgwrRemUseMax / PETgwrRemUse;    //[km3]

                            } else {
                                remainingUseGloLake = 0.;
                            }


                        }
                        else {    // PETgwrRemUse <= PETgwrRemUseMax

                            // This ODE (ordinary differential equation) applies for S >= 0 only:
                            // dS/dt = totalInflow - PET - gwr - NAs - k*S
                            // If the outflow becomes zero, S is computed again using a different ODE.

                            // Cell-specific calibration parameters - Use parameter
                            // glo_storageFactor (now: only one factor for surface water P_SWOUTF_C) in [1/d]
                            G_gloLakeStorage[n] = G_gloLakeStoragePrevStep * exp(-1. * P_SWOUTF_C__at__n)
                                                  + (1. / P_SWOUTF_C__at__n) * (totalInflow - PETgwrRemUse)
                                                  * (1. - exp(-1. * P_SWOUTF_C__at__n));

                            outflow = totalInflow + G_gloLakeStoragePrevStep - G_gloLakeStorage[n] - PETgwrRemUse;

                            // reduce G_gloLakeStorage[n] to maximum storage capacity
                            if (G_gloLakeStorage[n] > maxStorage) {
                                outflow += (G_gloLakeStorage[n] - maxStorage);
                                G_gloLakeStorage[n] = maxStorage;
                            }

                            if (outflow < 0.) {

                                // Different ODE (ordinary differential equation) applies for S < 0:
                                // dS/dt = totalInflow - PET - gwr - NAs  (without -k*S)
                                // S(t) = S(t-1) + inflow - PET - gwr - NAs

                                outflow = 0.;

                                // Storage equation without outflow:
                                G_gloLakeStorage[n] = G_gloLakeStoragePrevStep + totalInflow - PETgwrRemUse;
                            }

                            remainingUseGloLake = 0.;
                            // Reduction of gwr_glolak[n] and evaporation not required in this case.

                        }

                        if(abs(G_gloLakeStorage[n])<=minStorVol)    //to counter numerical inaccuracies
                            G_gloLakeStorage[n]=0.;

                        // If (PETgwrRemUse > PETgwrRemUseMax), PET is reduced by PETgwrRemUseRedFactor.
                        // If a reduction is not required, PETgwrRemUseRedFactor is set to 1.
                        if (PETgwrRemUse == 0.) {
                            PETgwrRemUseRedFactor = 1.;
                        } else {
                            PETgwrRemUseRedFactor = PETgwrRemUseMax / PETgwrRemUse;
                        }

                        if (PETgwrRemUseMax > PETgwrRemUse) {    // -> reduction of PET not required
                            PETgwrRemUseRedFactor = 1.;
                        }


                        if (2 == options.grid_store) {
                            // sum up the negative fraction of the potential cell runoff which has not been realized
                            // due to reduction of potential evaporation
                            G_monthlyPotCellRunoffDeficit(n,month)
                                           += (0.0 -
                                               (dailyWaterBalance.G_openWaterPET[n] * PETgwrRemUseRedFactor
                                                * (1.0 - gloLakeEvapoReductionFactor)
                                                * dailyWaterBalance.G_cellCorrFact[n]
                                                * G_lake_area[n]
                                                / 1000000.0));
                        }

                        inflow = outflow;
                        outflowGloLake = outflow;

                        dailyActualUse = remainingUseGloLakeStart - remainingUseGloLake;

                        //Evaporation reduction factor is now calculated for the next day.
                        // Cell-specific calibration parameters - Apply multiplier  // FP
                        G_gloLakeEvapoReductionFactor[n] = 1. -
                                                           pow(fabs(G_gloLakeStorage[n] - maxStorage)
                                                               / (2. * maxStorage),
                                                               (M_EVAREDEX__at__n_routing * evapoReductionExp));


                        if (G_gloLakeEvapoReductionFactor[n] < 0.)    // (could occur due to numerical inaccuracies)
                            G_gloLakeEvapoReductionFactor[n] = 0.;

                        if (G_gloLakeEvapoReductionFactor[n] > 1.)  // (could occur due to numerical inaccuracies)
                            G_gloLakeEvapoReductionFactor[n] = 1.;

                    }

                    // transport through reservoir if (options.resOpt == 1)
                    if (G_reservoir_area[n] > 0.) { //cell is a reservoir cell with a well known type

                        double c_ratio;
                        double prov_rel; //provisional release
                        double release;
                        double dailyUse;
                        double monthlyUse;

                        // define storage capacity to mean annual outflow ratio (c_ratio)
                        // G_mean_outflow:  m3/s,  G_stor_cap:  km3
                        // 31536000=(365 * 24 * 60 * 60)
                        c_ratio = G_stor_cap[n] / (G_mean_outflow[n] * 31536000. / 1000000000.);
                        // Reservoir operation start years
                        // maxStorage = G_reservoir_area[n] * 0.01; //10 m depth for active storage volume
                        // not calculated here because for global lake handling of regulated lakes before operation year we must seperate calculation.
                        maxStorage = G_stor_cap[n];

                        G_gloResMaxStorage[n] = maxStorage;    // required for Option "use_alloc" (water abstraction from second cell)

                        G_daily_res_in[n] = inflow; // [km3/day]

                        //assign correct storage value which was calculated in the last timestep or maxStorage at model day one
                        G_gloResStoragePrevStep = G_gloResStorage[n];

                        //same as gloLakeEvapo
                        // use reduction factor from previous time step (or the initial value)
                        gloResEvapoReductionFactor = G_gloResEvapoReductionFactor[n];

                        // G_dailyResEvapo[n]: see local lakes for explanation

                        G_dailyResEvapo[n] =
                                       ((1.0 - dailyWaterBalance.G_cellCorrFact[n])                            // [mm]
                                        * dailyWaterBalance.G_openWaterPrec[n])
                                       + (dailyWaterBalance.G_cellCorrFact[n]
                                          * (dailyWaterBalance.G_openWaterPET[n] * gloResEvapoReductionFactor));

                        if (G_dailyResEvapo[n] < 0.) { // avoid unphysical negative aet values
                            G_dailyResEvapo[n] =0.;
                        }


                        //inflow from upstream PLUS lake water balance
                        // reintroduced CFA in calculation (!no distinguishing for aridareaOpt!)
                        totalInflow = inflow + (dailyWaterBalance.G_openWaterPrec[n]
                                                * ((double) G_reservoir_area[n] / 1000000.));

                        G_ReservoirPrecip[n] = dailyWaterBalance.G_openWaterPrec[n] * ((double) G_reservoir_area[n]) /
                                               1000000.; // km3

                        // gwr below surface water bodies (except in arid inland sinks)
                        if ((1 == options.aridareaOpt) && (1 == G_aindex[n]) && (G_LDD[n] >= 0)) {
                            gwr_res[n] = Kswbgw * gloResEvapoReductionFactor * ((double) G_reservoir_area[n] /
                                                                                (cellArea * (geo.G_contfreq[n] /
                                                                                             100.))); // gwr below reservoirs in mm/d per continental area

                        } else {
                            gwr_res[n] = 0.;    // -> to avoid repeating the following equations for options.aridareaOpt==0
                        }

                        // remainingUse is evenly distributed between global reservoirs and global lakes (if existing) within a grid cell
                        // If remainingUse is negative (i.e. return flow to surface water bodies, from negative surface water abstractions G_dailydailyNUs)
                        //      then it is added to the storage
                        if (G_lake_area[n] >
                            0.) // If the cell also contains a global lake, 50% of NAs has already been subtracted from global lake storage.
                            remainingUseRes = 0.5 * remainingUse + remainingUseGloLake;

                        else
                            remainingUseRes = remainingUse;

                        // Required in 22b to compute actual use for Option 28 "aggrNUsGloLakResOpt"
                        remainingUseResStart = remainingUseRes;

                        // PET and gwr can occur until St = 0, while NAs is only withdrawn until St = 0.1*stor_capacity.
                        // Therefore, in case of reservoirs, the three variables cannot be equally satisfied.
                        // First, PET and gwr are satisfied. Then, if St > 0.1*stor_capacity, 'remainingUseRes' is subtracted from storage.

                        PETgwr = G_dailyResEvapo[n] *
                                 ((double) G_reservoir_area[n] / 1000000.) //[mm]->[km3]
                                 + gwr_res[n] * cellArea * (geo.G_contfreq[n] / 100.) / 1000000.;    //[mm]->[km3]

                        //WG22b: Water abstractions from reservoirs are reduced to 'remainingUseMax' to ensure that G_gloResStorage does not fall below 10% of the storage capacity.
                        PETgwrMax = G_gloResStoragePrevStep + totalInflow;

                        if (PETgwr > PETgwrMax) {

                            G_gloResStorage[n] = G_gloResStoragePrevStep + totalInflow - PETgwrMax;
                            gwr_res[n] *= PETgwrMax / PETgwr;
                            G_dailyResEvapo[n] *= PETgwrMax / PETgwr;

                        } else {
                            G_gloResStorage[n] = G_gloResStoragePrevStep + totalInflow - PETgwr;
                        }

                        // Treatment of NAs (net abstraction of surface water): add or substract from available storage
                        // ---
                        // Rationale behind net abstraction satisfaction, documented by Tim Trautmann Jan. 2022
                        // ---
                        // There are three processes to redistribute net abstraction from surface water bodies spatially
                        // and temporally a) aggregation of positive NAs over GLWD units, b) allocation of use to
                        // neigbhouring cell and c) delayed satisfaction.
                        // a) (options.aggrNUsGloLakResOpt): the water balance of global lakes, regulated lakes and
                        //   reservoirs are calculated in the outflow cell. Thus positive NAs is aggregated over
                        //   riparian cells of a global lake and is tried to be satisfied in outflow cell with global
                        //   lake and or reservoir storage. If not possible the part of this unsatisfied aggregated
                        //   NAs is reassign to originating riparian cell. This reassigning happens after the reservoir
                        //   resp. global lake storage. This share of NAs is tried to be satisfied first.
                        // b) (options.use_alloc): if a cell n cannot satisfy its own desired NAs it is tried to assign
                        //   this remaining use to the neighbouring cell (secondCell) if possible. If it cannot be
                        //   satisfied also in secondCell it is reassign to cell n in next time step. This reassigned
                        //   unsatisfied may be piled up with delayed satisfaction. The remaining use of a cell n is
                        //   satisfied first before the allocated use of another cell is satisfied.
                        // c) (options.delayedUseSatisfaction): the unsatisfied use of cell n can be used to be
                        //   satisfied for up to one year respecting the upper redistribution mechanism.
                        // The return flows into groundwater are adapted on daily basis by calculating a daily change
                        // of remaining use of cell n. For more details see updateNetAbstractionGW.
                        if (options.subtract_use > 0) {

                            // Negative values of NAs (return flow to surface water) are added (signs "- plus - equals +") to the storage:
                            if (remainingUseRes < 0.) {
                                G_gloResStorage[n] -= remainingUseRes;
                                remainingUseRes = 0.;
                            }
                                // Subtract remainingUseRes only if G_gloResStorage[n] > (0.1 * G_stor_cap[n]):
                            else {    // remainingUseRes >=0.
                                if (G_gloResStorage[n] > (0.1 * G_stor_cap[n])) {

                                    if (remainingUseRes < (G_gloResStorage[n] - (G_stor_cap[n] * 0.1))) {
                                        G_gloResStorage[n] -= remainingUseRes;
                                        remainingUseRes = 0.;
                                    }
                                        // When remainingUseRes is larger or equal to available storage
                                        // then subtract the available storage amount until 10% stor_cap is reached
                                    else {
                                        remainingUseRes -= (G_gloResStorage[n] - (G_stor_cap[n] * 0.1));
                                        G_gloResStorage[n] = G_stor_cap[n] * 0.1;
                                    }

                                }
                            }
                        }

                        if(abs(G_gloResStorage[n])<=minStorVol) //to counter numerical inaccuracies
                            G_gloResStorage[n]=0.;

                        // Different from local lakes and local wetlands, the outflow (release) of reservoirs is
                        // computed based on S(t) and not S(t-1).

                        // Calculate release (outflow) of reservoir:
                        // set rules at the beginning of the operational year! (only once!)
                        if (month == G_start_month[n] - 1) {
                            if (day == first_day_in_month[month]) {
                                    //1. calculate release coefficient for the actual year:
                                    //reduce release coefficent in this year to refill storage volume in reservoir
                                    if (G_gloResStorage[n] < (G_stor_cap[n] * 0.1))
                                        K_release[n] = 0.1;
                                    else
                                        K_release[n] = G_gloResStorage[n] / (maxStorage * 0.85);

                                //2 calculate annual release
                                annual_release[n] = K_release[n] * G_mean_outflow[n];
                            }
                        }

                        if ((G_res_type[n] + 0) == 1) {
                            // (irrigation reservoir)
                            //calculate monthly demand of downstream area
                            dailyUse = 0.0;

                            // read in net surface water use
                            dailyUse = G_dailydailyNUs[n];// [km3/day]

                            short i = 0;
                            int downstreamCell = G_downstreamCell[n];
                            // max reservoir_dsc cells downstream or to the basin outlet or area up to next reservoir
                            while (i < reservoir_dsc && downstreamCell > 0 && G_reservoir_area[downstreamCell] <= 0) {
                                dailyUse += G_dailydailyNUs[downstreamCell - 1] * G_alloc_coeff(n,i++);

                                // next downstream cell
                                downstreamCell = G_downstreamCell[downstreamCell - 1];
                            }

                            monthlyUse = dailyUse * numberOfDaysInMonth[month]; // [km3/month]
                            monthlyUse = monthlyUse * 1000000000. / (numberOfDaysInMonth[month] * 86400.); //[m3/s]

                            if (G_mean_demand[n] >= 0.5 * G_mean_outflow[n])
                                //provisional monthly release = i_mean /2. * (1. + monthly_demand_sum / d_mean_annual)
                                prov_rel = G_mean_outflow[n] / 2. * (1. + monthlyUse / G_mean_demand[n]); //[m3/s]
                            else
                                //provisional monthly release = i_mean + monthly_demand_sum - d_mean)
                                prov_rel = G_mean_outflow[n] + monthlyUse - G_mean_demand[n]; //[m3/s]
                        } else if ((G_res_type[n] + 0) == 2) {//(non-irrigation)
                            prov_rel = G_mean_outflow[n]; //[m3/s]
                        } else
                            cerr << "unknown reservoir type in gcrcNumber " << n + 1 << '\n';

                        // calculate release [m3/s]
                        if (c_ratio >= 0.5) {
                            release = K_release[n] * prov_rel;
                        } else { //monthly_release = (c_ratio/0,5)^2 * provisional monthly release + (1-(c_ratio/0,5)^2) * monthly_inflow (just daily inflow?)
                            release = ((4. * c_ratio * c_ratio) * K_release[n] * prov_rel)
                                      + ((1.0 - ((4. * c_ratio * c_ratio))) * inflow * 1000000000. / (24. * 3600.));
                        }

                        //new outflow
                        //reservoir storage volume should not be less than 10% of maximum!)...
                        if (G_gloResStorage[n] >= (G_stor_cap[n] * 0.1))
                            outflow = release * (24. * 3600.) / 1000000000.; //m3/s  ->  km3/routing time step
                            //...otherwise outflow will be reduced!
                            //(outflow should not be stopped, because we should serve ecosystem demands)
                        else
                            outflow = 0.1 * release * (24. * 3600.) / 1000000000.;

                        // to prevent that outflow (= river availability) gets negative in rare cases when net uses are negative downstream of a irrigation reservoir
                        if (outflow < 0.)
                            outflow = 0.;

                        // substract outflow from storage
                        G_gloResStorage[n] -= outflow;

                        // reduce G_gloResStorage to maximum storage capacity
                        if (G_gloResStorage[n] > maxStorage) {
                            outflow += (G_gloResStorage[n] - maxStorage);
                            G_gloResStorage[n] = maxStorage;
                        }

                        // if G_gloResStorage is below '0' there is no outflow anymore
                        // 'gloResEvapoReductionFactor' should prevent that G_gloResStorage falls below zero.
                        if (G_gloResStorage[n] < 0.) {
                            outflow += G_gloResStorage[n];
                            G_gloResStorage[n] = 0.;
                        }

                        //daily outflow from reservoirs
                        G_daily_res_out[n] = outflow; // [km3/day]

                        // If (PETgwr > PETgwrMax), PET is reduced by PETgwrRedFactor.
                        // If a reduction is not required, PETgwrRedFactor is set to 1.
                        if (PETgwr == 0.)
                            PETgwrRemUseRedFactor = 1.;
                        else
                            PETgwrRemUseRedFactor = PETgwrMax / PETgwr;

                        if (PETgwrMax > PETgwr) {    // -> reduction of PET not required
                            PETgwrRemUseRedFactor = 1.;
                        }


                        if (2 == options.grid_store)
                            // sum up the negative fraction of the potential cell runoff which has not been realized
                            // due to reduction of potential evaporation
                            G_monthlyPotCellRunoffDeficit(n,month)
                                           += (0.0 - (dailyWaterBalance.G_openWaterPET[n] * PETgwrRemUseRedFactor // ???
                                                      * (1.0 - gloResEvapoReductionFactor)
                                                      * dailyWaterBalance.G_cellCorrFact[n]
                                                      * G_reservoir_area[n]
                                                      / 1000000.0));

                        // reintroduced CFA
                        inflow = outflow; // outflow of reservoir is inflow to global wetlands

                        // introduced for 22b, Option 28 "aggrNUsGloLakResOpt"
                        dailyActualUse += remainingUseResStart - remainingUseRes;

                        //Evaporation reduction factor is now calculated for the next day.
                        G_gloResEvapoReductionFactor[n] = 1. - pow(fabs(G_gloResStorage[n] - maxStorage)
                                                        / maxStorage, evapoReductionExpReservoir);

                        if (G_gloResEvapoReductionFactor[n] < 0.)    // (could occur due to numerical inaccuracies)
                            G_gloResEvapoReductionFactor[n] = 0.;

                        if (G_gloResEvapoReductionFactor[n] > 1.)  // (could occur due to numerical inaccuracies)
                            G_gloResEvapoReductionFactor[n] = 1.;

                        additionalOutIn.additionalOutputInput(n, 5) = K_release[n]; // AEMS: save K_release

                    }
                    if (options.aggrNUsGloLakResOpt == 1 && (G_reservoir_area[n] > 0. || G_lake_area[n] > 0.)) {
                        double remainingUseAfterGloLakandRes = 0.; // is remainingUse (all allocated, unsatisfied etc. remaining use) after treating reservoirs and global lakes and regulated lakes
                        double rUseAllocUnsat = (totalDesiredUse - G_dailydailyNUsAggregated[n]); // variable needed to distinguish GLWD aggregated remaining Use and other remaining Use, (totalDesiredUse - G_dailydailyNUsAggregated[n]) = unsatisfied + allocated + unsatalloc + unsatUseRiparian use, which can only positive or zero
                        unsatGLWDUseOfRiparians = 0.;

                        if (G_reservoir_area[n] > 0.) {
                            remainingUseAfterGloLakandRes = remainingUseRes; // is correct if solely res in cell or both glo lak and res
                        } else if (G_lake_area[n] > 0.){
                            remainingUseAfterGloLakandRes = remainingUseGloLake; // only if solely glo lak in cell
                        }

                        for (int i = 0; i < ng; i++) {

                            if (G_glwdUnit[n] == G_glwdUnit[i]) {

                                if ((remainingUseAfterGloLakandRes - rUseAllocUnsat )> 0.) { // is the case if remainingUse of NAs aggregated is left

                                    if (G_dailydailyNUs[n] <= 0.) { // different calculation if own NAs is negative see documentation of 22b
                                        if (n == i ){ // i == n we are calculating the outflow cell
                                            if (G_dailydailyNUsAggregated[n] < 0) // shouldn't be reached as remainingUseAfterGloLakandRes - rUseAllocUnsat > 0.
                                                remainingUse = remainingUseAfterGloLakandRes;
                                            else
                                                remainingUse = rUseAllocUnsat;

                                            G_unsatUseRiparian[n] = 0.;
                                        } else { // riparian cell
                                        if (G_dailydailyNUs[i] <= 0.) {

                                            G_unsatUseRiparian[i] = 0.;

                                        } else {

                                            G_unsatUseRiparian[i] = (G_dailydailyNUs[i] /
                                                                     (G_dailydailyNUsAggregated[n] -
                                                                      G_dailydailyNUs[n])) *
                                                                    (remainingUseAfterGloLakandRes - rUseAllocUnsat); // (G_remainingUse[n] - G_dailydailyNUs[n]) always larger than zero in this loop
                                        }
                                        }
                                        //AEMS:
                                        if (day == last_day_in_month + 1) {
                                            additionalOutIn.additionalOutputInput(n, 42) = G_unsatUseRiparian[n];
                                            additionalOutIn.additionalOutputInput(i, 42) = G_unsatUseRiparian[i];
                                        }

                                    } else {  // if daily NAs of outflow cell is positive
                                        if (n == i) { // i == outflow cell (n)
                                            remainingUse = ( G_dailydailyNUs[n] / G_dailydailyNUsAggregated[n]) * (remainingUseAfterGloLakandRes - rUseAllocUnsat) + rUseAllocUnsat;
                                            G_unsatUseRiparian[n] = 0.;
                                        } else { // i == riparian cell
                                            if (G_dailydailyNUs[i] < 0.)
                                                G_unsatUseRiparian[i] = 0.;
                                            else
                                                G_unsatUseRiparian[i] =
                                                               G_dailydailyNUs[i] / G_dailydailyNUsAggregated[n] *
                                                               (remainingUseAfterGloLakandRes -
                                                                        rUseAllocUnsat);
                                        }
                                        //AEMS
                                        if (day == last_day_in_month + 1) {
                                            additionalOutIn.additionalOutputInput(n, 42) = G_unsatUseRiparian[n];
                                            additionalOutIn.additionalOutputInput(i, 42) = G_unsatUseRiparian[i];
                                        }
                                    }
                                } else {
                                    G_unsatUseRiparian[i] = 0.;
                                    G_unsatUseRiparian[n] = 0.;
                                    remainingUse = remainingUseAfterGloLakandRes;
                                    //AEMS
                                    if (day == last_day_in_month + 1) {
                                        additionalOutIn.additionalOutputInput(i, 42) = G_unsatUseRiparian[i];
                                        additionalOutIn.additionalOutputInput(n, 42) = G_unsatUseRiparian[n];   //before n,48 but 48 never read in
                                    }
                                }
                                if ((n != i) && (G_dailydailyNUs[i] > 0)){
                                    G_monthlyRedistributeGLWDUse(i,month) += (G_dailydailyNUs[i] -
                                                   G_unsatUseRiparian[i]);

                                    G_monthlyRedistributeGLWDUse(n,month) -= (G_dailydailyNUs[i] -
                                                   G_unsatUseRiparian[i]);
                                }
                            }
                        }
                        unsatGLWDUseOfRiparians = remainingUseAfterGloLakandRes - remainingUse;
                    } else if (options.aggrNUsGloLakResOpt == 0 && (G_reservoir_area[n] > 0. || G_lake_area[n] > 0.)) {
                        if (G_reservoir_area[n] > 0.) {
                            remainingUse = remainingUseRes; // is correct if solely res in cell or both glo lak and res
                        } else if (G_lake_area[n] > 0.){
                            remainingUse = remainingUseGloLake; // only if solely glo lak in cell
                        }
                    }

                    // outflow of a global lake is inflow to global wetlands


                    // ------------------  Routing through global wetlands --------------------------------------------
                    if (G_glo_wetland[n] > 0) {

                        // assign correct storage value (for t) which was calculated in the last timestep or maxStorage at model day one
                        G_gloWetlStoragePrevStep = G_gloWetlStorage[n];

                        // maximum storage capacity
                        // Cell-specific active wetland depth - Use array (initialized from cell-specific parameter)
                        maxStorage = ((G_glo_wetland[n]) / 100.) * cellArea * G_wetlDepthActive[n];

                        // This factor was added to reduce PET as a function of actual wetland storage (2.1f).
                        // use reduction factor from previous time step (or the initial value)
                        gloWetlAreaReductionFactor = G_gloWetlAreaReductionFactor[n];


                        // G_dailyGloWetlEvapo[n]: see local lakes for explanation

                        G_dailyGloWetlEvapo[n] =
                                       ((1.0 - dailyWaterBalance.G_cellCorrFact[n])                            // [mm]
                                        * (dailyWaterBalance.G_openWaterPrec[n] * gloWetlAreaReductionFactor))
                                       + (dailyWaterBalance.G_cellCorrFact[n]
                                          * (dailyWaterBalance.G_openWaterPET[n] * gloWetlAreaReductionFactor));

                        if (G_dailyGloWetlEvapo[n] < 0.) { // avoid unphysical negative aet values
                            G_dailyGloWetlEvapo[n] =0.;
                        }



                        // calculate inflow only in dependence of reduction factor to avoid inconsistent counting of precipitation
                        totalInflow = inflow + (dailyWaterBalance.G_openWaterPrec[n] * gloWetlAreaReductionFactor
                                                * (cellArea / 1000000.) * (G_glo_wetland[n] / 100.));    //[mm]->[km3]

                        G_glowetlextent[n] = gloWetlAreaReductionFactor * cellArea * (G_glo_wetland[n] / 100.);
                        // gwr below surface water bodies (except in arid inland sinks)
                        if ((1 == options.aridareaOpt) && (1 == G_aindex[n]) && (G_LDD[n] >= 0)) {
                            gwr_glowet[n] = Kswbgw * gloWetlAreaReductionFactor * G_glo_wetland[n] / 100. /
                                            (geo.G_contfreq[n] / 100.); // gwr below global wetlands in mm/d

                        } else
                            gwr_glowet[n] = 0.; // -> to avoid repeating the following equations for options.aridareaOpt==0

                        PETgwr = G_dailyGloWetlEvapo[n] * (cellArea / 1000000.) *
                                 ((G_glo_wetland[n]) / 100.)                        //[mm]->[km3]
                                 + gwr_glowet[n] * cellArea * (geo.G_contfreq[n] / 100.) / 1000000.;    //[mm -> km3]

                        // There are two options to compute PETgwrMax: Based on the analytical solution of 1) dS/dt=In-Out-k*S; 2) dS/dt=In-Out
                        // In the first equation, the outflow is appr. 0.6% of S(t-1), although the storage is empty at the end of the timestep.
                        // In this case, it seems more realistic if the remaining water evaporates or infiltrates (second equation).

                        // PETgwrMax = totalInflow + (glo_storageFactor * G_gloWetlStoragePrevStep * exp(-1. * glo_storageFactor)) / (1. - exp(-1 * glo_storageFactor));

                        PETgwrMax = totalInflow + G_gloWetlStoragePrevStep;

                        if (PETgwr > PETgwrMax) {    // reduction of PET and gwr_glowet required

                            G_gloWetlStorage[n] = 0.;
                            outflow = 0.;
                            gwr_glowet[n] *= PETgwrMax / PETgwr;        // [mm]
                            G_dailyGloWetlEvapo[n] *= PETgwrMax / PETgwr;        // [mm]

                        } else {

                            // Cell-specific calibration parameters - Use parameter
                            // glo_storageFactor (now: only one factor for surface water P_SWOUTF_C) in [1/d]
                            G_gloWetlStorage[n] = G_gloWetlStoragePrevStep * exp(-1. * P_SWOUTF_C__at__n)
                                           + (1. / P_SWOUTF_C__at__n) * (totalInflow - PETgwr)
                                           * (1. - exp(-1. * P_SWOUTF_C__at__n));

                            outflow = totalInflow + G_gloWetlStoragePrevStep - G_gloWetlStorage[n] - PETgwr;
                        }

                        if (G_gloWetlStorage[n] > maxStorage) {
                            outflow += (G_gloWetlStorage[n] - maxStorage);
                            G_gloWetlStorage[n] = maxStorage;
                        }

                        if(abs(G_gloWetlStorage[n])<=minStorVol)    //to counter numerical inaccuracies
                            G_gloWetlStorage[n]=0.;

                        // If (PETgwr > PETgwrMax), PET is reduced by PETgwrRedFactor. If a reduction is not required,
                        // PETgwrRedFactor is set to 1.
                        // if PETgwr is 0 (because gwr_swb is 0 because reductionfactor is 0 because storage is empty),
                        // division through 0 --> has to be catched, because otherwise G_CELL_RUNOFF_TOTAL
                        // and G_CELL_AET gets nan.
                        if (PETgwr == 0.)
                            PETgwrRedFactor = 1.;
                        else
                            PETgwrRedFactor = PETgwrMax / PETgwr;

                        if (PETgwrMax > PETgwr)    // -> reduction of PET not required
                            PETgwrRedFactor = 1.;

                        if (2 == options.grid_store)
                            // sum up the negative fraction of the potential cell runoff which has not been realized
                            // due to reduction of potential evaporation
                            G_monthlyPotCellRunoffDeficit(n,month)
                                           += (0.0 - dailyWaterBalance.G_openWaterPET[n] * PETgwrRedFactor
                                                     * (1.0 - gloWetlAreaReductionFactor)
                                                     * dailyWaterBalance.G_cellCorrFact[n]
                                                     * (cellArea / 1000000.0)
                                                     * (G_glo_wetland[n] / 100.0));


                        inflow = outflow;
                        outflowGloWetland = outflow;

                        // Area reduction factor is now calculated for the next day.
                        // Cell-specific calibration parameters - Apply multiplier  // FP
                        G_gloWetlAreaReductionFactor[n] = 1. -
                                                          pow(fabs(G_gloWetlStorage[n] - maxStorage)
                                                              / maxStorage,
                                                              (M_EVAREDEX__at__n_routing * evapoReductionExp));

                        if (G_gloWetlAreaReductionFactor[n] < 0.)  // (could occur due to numerical inaccuracies)
                            G_gloWetlAreaReductionFactor[n] = 0.;

                        if (G_gloWetlAreaReductionFactor[n] > 1.)  // (could occur due to numerical inaccuracies)
                            G_gloWetlAreaReductionFactor[n] = 1.;

                    }    // end if (G_glo_wetland[n] > 0)



                    // If GWR under SWB / fractionalRouting is allowed,
                    // GW storage is computed before local lakes in humid cells (baseflow directed into local lake)
                    // and before the river in (semi-)arid cells (baseflow directed into river)

                    if ((1 == options.aridareaOpt) && (1 == G_aindex[n]) && (G_LDD[n] >= 0)) {        // (semi-)arid cells; fractionalRoutingOpt automatically set to 1, if aridareaOpt==1

                        // sum up total groundwater recharge [mm/d] from surface water bodies and calculate fraction
                        // of actual surface waterbody area
                        // gwr_swb is added to the groundwater storage in the next timestep

                        G_gwr_swb[n] = gwr_loclak[n] + gwr_glolak[n] + gwr_locwet[n] + gwr_glowet[n] +
                                       gwr_res[n]; // sum up gwr but do not use gwr_rivers because rivers are the final storage [mm/d]
                        G_groundwaterStoragePrevStep = G_groundwaterStorage[n];

                        // WG22b: ODE (ordinary differential equation) is applied for groundwater storage:
                        // dS/dt = GWR - NAg - k*S  (k*S = groundwater_runoff_km3)

                        // Groundwater recharge below surface water bodies from previous timestep:
                        netGWin = G_gwr_swb[n] * cellArea * (geo.G_contfreq[n] / 100.) /
                                  1000000.                        // mm -> km3
                                  + (double) dailyWaterBalance.G_dailyGwRecharge[n] * cellArea *
                                    ((double) G_landAreaFrac[n] / 100.) / 1000000.;    // mm -> km3

                        // G_dailydailyNUg is adapted in this timestep based on the remaining use of last time step
                        if (options.subtract_use > 0) {
                            if (G_dailyRemainingUse[n] != 0.) {
                                G_dailydailyNUg[n] = updateNetAbstractionGW(n);
                            }
                            netGWin -= G_dailydailyNUg[n];
                        }

                        // Analytical solution of dS/dt = netGWR - k*S
                        // k_g (now: P_GWOUTF_C) in [1/d]
                        // Cell-specific calibration parameters - Use parameter

                        G_groundwaterStorage[n] = G_groundwaterStoragePrevStep * exp(-1. * P_GWOUTF_C__at__n_routing) +
                                                  (1. / P_GWOUTF_C__at__n_routing) * netGWin *
                                                  (1. - exp(-1. * P_GWOUTF_C__at__n_routing));

                        if(abs(G_groundwaterStorage[n])<=minStorVol)    //to counter numerical inaccuracies
                            G_groundwaterStorage[n]=0.;


                        groundwater_runoff_km3 = G_groundwaterStoragePrevStep - G_groundwaterStorage[n] + netGWin;



                        if (groundwater_runoff_km3 <= 0.) {
                            // different differential equation (ODE) applies:
                            // dS/dt = GWR - NAg (without k*S) -> S(t) = S(t-1) + netGWin

                            groundwater_runoff_km3 = 0.;

                            G_groundwaterStorage[n] = G_groundwaterStoragePrevStep + netGWin;

                            if(abs(G_groundwaterStorage[n])<=minStorVol)    //to counter numerical inaccuracies
                                G_groundwaterStorage[n]=0.;

                        }
                        if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store))
                            G_monthlyGwRunoff(n, month) += groundwater_runoff_km3;
                        // if daily 365 output and gwout (same for 31; same for land area fraction)
                        // G_daily365GwRunoff(n,day - 1) = groundwater_runoff_km3;

                        // in semi-arid/arid areas, all groundwater reaches the river directly
                        G_localGWRunoffIntoRiver[n] = groundwater_runoff_km3;

                        // Firstly, amount of local runoff flowing directly into the river is calculated, then local gw runoff is added for the rest term.

                        // G_localRunoff[n] and G_localRunoffIntoRiver[n] in arid cells calculated before local lakes.
                        // reintroduced for CFA calculation
                        potCellRunoff = (G_localRunoff[n] + G_localRunoffIntoRiver[n] +
                                         G_localGWRunoffIntoRiver[n]    // is in km3
                                         + ((dailyWaterBalance.G_lakeBalance[n]    // is in mm
                                             / 1000000.0)    // mm ->*/ km
                                            * (((double) G_lake_area[n] + (double) G_reservoir_area[n])    // km2
                                               + ((geo.areaOfCellByArrayPos(n) / 100.0)    // % -> km2
                                                  * (G_loc_lake[n] + G_loc_wetland[n] + G_glo_wetland[n])))));    // [%]

                        G_potCellRunoff[n] += potCellRunoff;

                        if (2 == options.grid_store) {
                            G_monthlyPotCellRunoff(n, month) += potCellRunoff;
                        }

                    }
                    // end of computation of GW storage and baseflow for options.aridareaOpt == 1 && (semi-)arid cells

                    // calculate the amount of water which flows out of local/global lakes/wetlands
                    // and contributes to the actual cell runoff within the cell
                    // Different from WG22, cell runoff is computed as (transportedVolume - inflow from upstream cells) in WG2.2b .
                    // if (1 == options.fractionalRoutingOpt) // reduce G_riverInflow by the amount of water which flows not through the storage cascade
                    // CellRunoff = (inflow - G_riverInflow[n] + G_localRunoffIntoRiver[n] + G_localGWRunoffIntoRiver[n]);

                    // else  // normal version (without fractional routing)
                    // CellRunoff = (inflow - G_riverInflow[n]);

                    // if ((0 == options.fractionalRoutingOpt) && (G_aindex[n] == 1) && (options.aridareaOpt == 1)) // in arid regions with gwr_swb groundwater should run directly into the river
                    // CellRunoff = (inflow - G_riverInflow[n] + G_localGWRunoffIntoRiver[n]);

                    // outflow from global wetlands is inflow to the river system
                    G_riverInflow[n] = inflow;



                    // At this point G_riverInflow is the total inflow to the river system:
                    // runoff->loca_lakes->local_wetlands->global_lakes->global_wetlands->river_network

                    if ((1 ==
                        options.fractionalRoutingOpt) && (G_LDD[n] >= 0)) { // if option is activated, add the water which is flowing directly to the river
                        G_riverInflow[n] += G_localRunoffIntoRiver[n];
                        G_riverInflow[n] += G_localGWRunoffIntoRiver[n];
                    }

                    // riverVelocity in [km/d]

                    if (options.riverveloOpt == 0) {
                        riverVelocity = defaultRiverVelocity;
                    }
                    else if (options.riverveloOpt == 1) {
                        // Cell-specific calibration parameters - Apply multiplier M_RIVRGH_C within method getRiverVelocity // FP
                        riverVelocity = getRiverVelocity(G_RiverSlope[n], G_riverBottomWidth[n], G_Roughness[n],
                                                         G_riverInflow[n], n, calParam);// OE: added calParam
                    }
                    else if (options.riverveloOpt == 2) {
                        riverVelocity = getNewRiverVelocity(G_RiverSlope[n], G_riverBottomWidth[n], G_Roughness[n],
                                                            G_riverStorage[n], G_riverLength[n], n, calParam);// OE: added calParam
                    }

                    // WG22b: To be consistent with global lakes/wetlands/groundwater (storages with ord. diff.
                    // equation), K is expressed in 1/d.
                    K = riverVelocity / G_riverLength[n]; // [(km/d)/km] = [1/d]

                    //assign correct storage value which was calculated in the last timestep or maxStorage at model day one

                    G_riverStoragePrevStep = G_riverStorage[n];

                    // Computation of river evaporation prior to river storage equation to implement G_dailyRiverEvapo
                    // into storage equation.
                    // River evaporation start
                    if (options.riverEvapoOpt == 1) {

                        // use reduction factor from previous time step (or the initial value)
                        //riverAreaReductionFactor = G_riverAreaReductionFactor[n];	// riverAreaReductionFactor is not used anymore

                        G_riverAreaFrac[n] = G_riverAreaFracNextTimestep_Frac[n];

                        //HMS both calculated later if riverareafrac is finally calculated
                        G_dailyRiverEvapo[n] = ((1.0 - dailyWaterBalance.G_cellCorrFact[n]) *
                                                (dailyWaterBalance.G_openWaterPrec[n])
                                                + (dailyWaterBalance.G_cellCorrFact[n] *
                                                   dailyWaterBalance.G_openWaterPET[n]))
                                               * G_riverAreaFrac[n] / 100. * cellArea / 1000000.;

                        G_dailyRiverPrecip[n] =
                                       dailyWaterBalance.G_openWaterPrec[n] * G_riverAreaFrac[n] / 100. * cellArea /
                                       1000000.; // precip on river fraction in km3

                    } else {// (no evaporation from rivers)
                        G_dailyRiverPrecip[n] = 0.;
                        G_dailyRiverEvapo[n] = 0.;
                    }

                    G_riverInflow[n] += G_dailyRiverPrecip[n];

                    // At this point, RiverEvapoRemUse can get negative because remainingUse can be negative!
                    RiverEvapoRemUse = remainingUse +
                                       G_dailyRiverEvapo[n]; // [km3] // remainingUse = 0, if options.subtract_use == 0

                    // introduced for 22b, Option 28 "aggrNUsGloLakResOpt"
                    remainingUseRiverStart = remainingUse;

                    // WG22b: Unit of K changed from [d] to [1/d]
                    // Different from global wetlands, RivEvapoRemUseMax is computed based on the ODE dS/dt=In-Out-k*S.
                    // (outflow is assumed to occur even if S = 0 at the end of the time step.
                    // See comment in global wetland loop)
                    RivEvapoRemUseMax = G_riverInflow[n] + (K * G_riverStoragePrevStep * exp(-1. * K))
                                                           / (1. - exp(-1. * K));

                    if (RiverEvapoRemUse > RivEvapoRemUseMax) {

                        G_riverStorage[n] = 0.;

                        transportedVolume = G_riverInflow[n] + G_riverStoragePrevStep - RivEvapoRemUseMax;

                        if (transportedVolume < 0.)   // to avoid the negative transportedvolume error
                            transportedVolume = 0.;

                        if (remainingUse > 0.)
                            remainingUse -= remainingUse * RivEvapoRemUseMax / RiverEvapoRemUse;
                        else
                            remainingUse = 0.;

                        G_dailyRiverEvapo[n] *= RivEvapoRemUseMax / RiverEvapoRemUse; // [km3]

                    } else {


                        // WG22b: Unit of K changed from [d] to [1/d]
                        G_riverStorage[n] = G_riverStoragePrevStep * exp(-1. * K)
                                            + (1. / K) * (G_riverInflow[n] - RiverEvapoRemUse)
                                              * (1. - exp(-1. * K));

                        if(abs(G_riverStorage[n])<=minStorVol)  //to counter numerical inaccuracies
                            G_riverStorage[n]=0.;

                        transportedVolume =
                                       G_riverInflow[n] + G_riverStoragePrevStep - G_riverStorage[n] - RiverEvapoRemUse;

                        if (transportedVolume < 0.)       // to avoid the negative transportedvolume error
                            transportedVolume = 0.;
                        remainingUse = 0.;

                    }


                    if (transportedVolume <
                        0.) {    // (due to numerical inaccuracies)
                        cout << "ERROR transported volume < 0 in cell: " << n << "routing.routing: year-month-day "
                             << year << "-" << month << "-" << day_in_month << " - julian day: " << day << " with: "
                             << transportedVolume << endl;

                        cout << "transportedVolume = G_riverInflow[n]: " << G_riverInflow[n]
                             << " + G_riverStoragePrevStep: " << G_riverStoragePrevStep << endl;
                        cout << " - G_riverStorage[n]: " << G_riverStorage[n] << " - RiverEvapoRemUse: "
                             << RiverEvapoRemUse << endl;
                        cout << "RivEvapoRemUseMax: " << RivEvapoRemUseMax << ", remainingUse (set to zero): "
                             << remainingUse << endl;
                        cout << "K: " << K << ", riverVelocity: " << riverVelocity << ", G_riverLength[n]: "
                             << G_riverLength[n] << endl;
                        cout << "G_riverAreaFrac[n]: " << G_riverAreaFrac[n] << ", G_dailyRiverEvapo[n]: "
                             << G_dailyRiverEvapo[n] << ", G_dailyRiverPrecip[n]: " << G_dailyRiverPrecip[n] << endl;

                    }


                    // Calculate river width (based on storage at the end of current time step) according to new river
                    // flow velocity approach (not implemented for river velocity yet,
                    // but documented in documentation for WG22b)

                    // In WG22, G_RiverWidth was a function of "transportedVolume" and could therefore drop to zero.
                    // G_RiverWidth[n] = 2.71 * pow(transportedVolume * 1/0.0000864, 0.557); // calculate G_RiverWidth in m after Allen et al. 1994; transportedVolume has to be converted from km3/d into m3/s (1/((60*60*24)/1000000000)
                    // In WG22b, the minimum river width is the river bottom width meaning that river evaporation could occur even if the river depth is zero.
                    // Therefore, G_riverAreaFrac[n] is reduced by the G_riverAreaReductionFactor.

                    if (options.riverEvapoOpt == 1) {

                        crossSectionalArea = G_riverStorage[n] / G_riverLength[n];    // [km3/km = km2]
                        // G_riverBottomWidth[n]: m -> km
                        riverDepth = -G_riverBottomWidth[n] / (4. * 1000.) +
                                     sqrt(G_riverBottomWidth[n] / 1000. * G_riverBottomWidth[n] / (16. * 1000.) +
                                          0.5 * crossSectionalArea); // [km]
                        G_RiverWidth[n] = G_riverBottomWidth[n] / 1000. + 4. * riverDepth; //[km]

                        if (G_RiverWidth[n] > G_RiverWidth_bf[n] / 1000.) //G_RiverWidth_bf[n] in m
                            G_RiverWidth[n] = G_RiverWidth_bf[n] / 1000.;


                        // G_riverAreaReductionFactor and G_riverAreaFracNextTimestep are now calculated for the next day.
                        // Cell-specific calibration parameters - Apply multiplier  // FP
                        G_riverAreaReductionFactor[n] = 1. -
                                                        pow(fabs(G_riverStorage[n] - G_riverStorageMax[n])
                                                            / G_riverStorageMax[n],
                                                            (M_EVAREDEX__at__n_routing * evapoReductionExp));

                        if (G_riverAreaReductionFactor[n] < 0.)  // (could occur due to numerical inaccuracies)
                            G_riverAreaReductionFactor[n] = 0.;

                        if (G_riverAreaReductionFactor[n] > 1.)  // (could occur due to numerical inaccuracies)
                            G_riverAreaReductionFactor[n] = 1.;

                        // Only continental cells
                        if (geo.G_contcell[n]) {
                            G_riverAreaFracNextTimestep_Frac[n] =
                                           G_riverAreaReductionFactor[n] * G_riverLength[n] * G_RiverWidth[n] * 100. /
                                           cellArea;
                        }
                            // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                            // to avoid not defined G_riverAreaFrac at Caspian Sea
                        else {
                            G_riverAreaFracNextTimestep_Frac[n] = 0.;
                        }

                        // Calculate the necessary percent change toward the next time step
                        G_riverAreaFrac_ChangeFrac[n] = G_riverAreaFracNextTimestep_Frac[n] - G_riverAreaFrac[n];
                    }

                    // WG22b: Different from global lakes and reservoirs, remaining use (if > 0) is subtracted
                    // from local lake storage at the end of routing.
                    if (1 <= options.subtract_use) {

                        if ((G_loc_lake[n] > 0.) && (remainingUse > 0.)) {

                            // WG22b: Water can be abstracted as long as locLakeStorage exceeds -maxStorage
                            // (and not zero, as it was the case in all previous versions).
                            if (G_locLakeStorage[n] > (-1.) * G_locLakeMaxStorage[n]) {
                                if (remainingUse < (G_locLakeMaxStorage[n] + G_locLakeStorage[n])) {
                                    G_locLakeStorage[n] -= remainingUse;
                                    remainingUse = 0.;
                                } else {
                                    remainingUse -= (G_locLakeMaxStorage[n] + G_locLakeStorage[n]);
                                    G_locLakeStorage[n] = (-1.) * G_locLakeMaxStorage[n];
                                }
                            }

                            // Area reduction factor for the next day is calculated again based on the reduced local
                            // lake storage.
                            // Cell-specific calibration parameters - Apply multiplier  // FP
                            G_locLakeAreaReductionFactor[n] = 1. -
                                                              pow(fabs(G_locLakeStorage[n] - G_locLakeMaxStorage[n])
                                                                  / (2. * G_locLakeMaxStorage[n]),
                                                                  (M_EVAREDEX__at__n_routing * evapoReductionExp));

                            if (G_locLakeAreaReductionFactor[n] <
                                0.)
                                G_locLakeAreaReductionFactor[n] = 0.;

                            if (G_locLakeAreaReductionFactor[n] > 1.)  // (could occur due to numerical inaccuracies)
                                G_locLakeAreaReductionFactor[n] = 1.;

                        }
                        // end of abstraction from local lakes ((G_loc_lake[n] > 0.) && (remainingUse > 0.))

                        // introduced for 22b, Option 28 "aggrNUsGloLakResOpt"
                        dailyActualUse += remainingUseRiverStart - remainingUse;



                        // satisfied use of the current day summed up for the whole year
                        // G_satisfiedUse[n] += totalDesiredUse - remainingUse;	// WG22
                        // for new Use allocation (M.Hunger 2/2006)
                        // all actual satisfied water uses
                        // G_actualUse[n] += totalDesiredUse - remainingUse;	// WG22

                        // WG22b: REDISTRIBUTE UNSATISFIED ALLOCATED USE OVER RESPECTIVE NEIGHBORING CELLS:
                        // Rationale: First priority: Satisfy water demand of cell n (from water storage in cell n).
                        // Second priority: Satisfy water demand allocated from neighboring cell(s) (from water storage
                        // in cell n).

                        if ((options.use_alloc > 0) && (G_dailyAllocatedUse[n] > 0.)) {

                            if (remainingUse > 0.) {

                                if (demandWithoutAllocUse > 0.) {

                                    unsatUse = remainingUse; // remainingUse of cell n + allocated use

                                    remainingUse = demandWithoutAllocUse - dailyActualUse - unsatGLWDUseOfRiparians;    // remainingUse of cell n

                                    if (remainingUse < 0.)  // since "dailyActualUse" can exceed "demandWithoutAllocUse"
                                        remainingUse = 0.;

                                    unsatAllocUse = unsatUse -
                                                    remainingUse;    // total unsatisfiedUse (including allocated use)
                                                                     // minus remainingUse of cell n

                                } else {    // if (demandWithoutAllocUse <= 0.) -> negative G_dailydailyNUs

                                    // The initial water demand of cell n (demandWithoutAllocUse) is <= 0.
                                    // In this case, the remainingUse at this point only stems from the respective
                                    // neighboring cells.
                                    // The remainingUse of cell n is set to zero.
                                    unsatAllocUse = remainingUse;    // 	remainingUse to be redistributed over respective nbc

                                    remainingUse = 0.;    // remainingUse of cell n

                                }

                            }
                                // end (remainingUse > 0.)

                            else {    // if (remainingUse == 0)

                                unsatAllocUse = 0.;
                            }


                            for (i_nbc = 0; i_nbc < ng; i_nbc++) {

                                if (G_SecondCell[i_nbc] == n) {
                                    // Scale G_UnsatAllocUse[i_nbc] to just the allocated use coming from cell i_nbc
                                    G_UnsatAllocUse[i_nbc] = G_AllocUseToNeigborCell[i_nbc] * unsatAllocUse /
                                                             G_dailyAllocatedUse[n];    // At this point, G_dailyAllocatedUse[n] is always > 0.
                                    if (options.delayedUseSatisfaction == 0){
                                        G_dailyRemainingUse[i_nbc] = G_UnsatAllocUse[i_nbc];
                                    }else{
                                        // The adaption of NAg shall be done on a daily basis, hence the difference from
                                        // prev timestep is calculated. It can happen that secondCell is used in one
                                        // timestep and not in another thus both has to be subtracted.
                                        G_dailyRemainingUse[i_nbc] = (G_UnsatAllocUse[i_nbc] -
                                                                      G_PrevUnsatAllocUse[i_nbc] -
                                                                      G_PrevTotalUnstatisfiedUse[i_nbc]);
                                    }

                                    // In order to have information about the total satisfied water use in the
                                    // cell where it was demanded an array is created which contains information
                                    // about which saturated allocated use should be subtracted and which added where
                                    // to have this as output option.
                                    double sat_alloc_use = (G_AllocUseToNeigborCell[i_nbc] -
                                                                         G_UnsatAllocUse[i_nbc]);
                                    G_monthlyRedistributeAllocUse(i_nbc,month) += sat_alloc_use;
                                    G_monthlyRedistributeAllocUse(n,month) -= sat_alloc_use;

                                    // If "G_UnsatAllocUse[i_nbc]" will be added to the "remainingUse" of "i_nbc" the
                                    // next day, this amount of unsatisfied use must be taken into account to check the
                                    // following condition:
                                    // sum(actual use) + sum(unsat. use) = sum(NAs)

                                    if (options.delayedUseSatisfaction == 0) {

                                        G_totalUnsatisfiedUse[i_nbc] += G_UnsatAllocUse[i_nbc];
                                        G_UnsatAllocUse[i_nbc] = 0.;

                                    } else {    // options.delayedUseSatisfaction == 1:

                                        // used as output only:
                                        if (G_CellPositionRoutOrder[i_nbc] < G_CellPositionRoutOrder[n])
                                            G_daily_UnsatAllocUseNextDay[i_nbc] = G_UnsatAllocUse[i_nbc];

                                        else
                                            G_daily_UnsatAllocUseNextDay[i_nbc] = 0.;

                                        // WG22b: As in previous versions, subtraction of NUs can be delayed for up to one year.
                                        if (day == 365 &&
                                            (G_CellPositionRoutOrder[i_nbc] < G_CellPositionRoutOrder[n])) {
                                            G_totalUnsatisfiedUse[i_nbc] = G_UnsatAllocUse[i_nbc];
                                            G_UnsatAllocUse[i_nbc] = 0.;
                                            G_dailyRemainingUse[i_nbc] = 0.;
                                        }

                                    }
                                }
                            }
                        }

                        G_actualUse[n] += dailyActualUse;    // required to write annual value of G_actualUse[n]

                        // calculate satisfied allocated use in second cell for output
                        G_dailySatisAllocatedUseInSecondCell[n] = G_dailyAllocatedUse[n] - unsatAllocUse;


                        // monthly satisfied use: satisfied uses from this cell only
                        // monthly actual use: all satisfied uses (see below)
                        if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {
                            //    G_monthlySatisfiedUse(n,month) += totalDesiredUseCell - remainingUse;	// CR 2015-09: A distinction between satisfied use and actual use is not possible in 22b.
                            G_monthlyActualUse(n,month) += dailyActualUse;
                            G_monthlyNUg(n,month) += G_dailydailyNUg[n];
                            // Different from WG22, G_monthlyAllocatedUse does not necessarily have to be fully satisfied in the neighboring cell (nbc).
                            // It is only the water demand allocated to a nbc (without taking into account the current river and SWB storage in the nbc).
                            G_monthlyAllocatedUse(n,month) += G_dailyAllocatedUse[n];
                            G_monthlySatisAllocatedUseinSecondCell(n,month) += G_dailySatisAllocatedUseInSecondCell[n];
                        }


                        // ALLOCATION OF REMAINING USE TO SECOND CELL
                        // At this point, 'remainingUse' only includes unsatisfied water demand of the cell itself.
                        // Unsatisfied allocated use has been redistributed over the respective neighboring cells before.

                        if (options.use_alloc > 0) {

                            if (remainingUse > 0.) {    // Identification of secondCell

                                if (G_toBeCalculated[n] == 1) {
                                    // a second cell is allowed to satisfy the remaining demand
                                    short i;
                                    int nbc;
                                    double storageSum;

                                    // choose the neighbouring cell with the highest
                                    // amount of river and lake storage
                                    // which is not upstream of cell n
                                    totalNeighbourStorage = 0.;
                                    secondCell = -1; // standard setting "no downstream cell"

                                    for (i = 0; i < 8; i++) {

                                        nbc = G_neighbourCell(n,i) - 1; // GCRC id -1 of neighbour cell

                                        // For non-continent cells (Caspian Sea in SLM) G_contcell[nbc] == 0,
                                        // set nbc to unaccepted value "-1", to exclude this cell from further calculations here
                                        if (geo.G_contcell[nbc] == 0) {
                                            nbc = -1;
                                        }

                                        if (nbc < 0) {
                                            continue; //  neighbouring cell does not exist
                                        } else {
                                            if (G_toBeCalculated[nbc] == 1) {
                                                // " + G_locLakeMaxStorage[nbc]" and " + G_gloLakeMaxStorage[nbc]" inserted,
                                                // since in WG22b local/global lake storage can decrease to -maxStorage.
                                                // Otherwise, 'storageSum' could become negative, and the initial value of secondCell = -1 would not be updated (secondCell = nbc).

                                                storageSum = G_riverStorage[nbc]
                                                                + G_locLakeStorage[nbc] + G_locLakeMaxStorage[nbc]
                                                                + G_gloLakeStorage[nbc] + G_gloLakeMaxStorage[nbc];

                                                if (options.resOpt == 1) {

                                                    storageSum += G_gloResStorage[nbc];

                                                }

                                              if(storageSum < 1.e-12) //to counter numerical inaccuracies
                                                  storageSum = 0.;

                                                if (storageSum > totalNeighbourStorage) {

                                                    // check if neighbour cell is a direct upstream cell
                                                    if ((G_downstreamCell[nbc] - 1) != n) {
                                                        // set second cell only if it is not a direct
                                                        // upstream cell of current cell n
                                                        totalNeighbourStorage = storageSum;
                                                        secondCell = nbc;

                                                    }
                                                }

                                            }

                                        }

                                    }

                                }

                                // WG22b: As in previous versions, subtraction of NUs can be delayed for up to one year.
                                if (day == 365 && (G_CellPositionRoutOrder[secondCell] < G_CellPositionRoutOrder[n]))
                                    secondCell = -1;

                            }
                            else {
                                secondCell = -1;
                            }

                            if (secondCell < 0) {    // remainingUse can be >= 0 at this point

                                G_AllocUseToNeigborCell[n] = 0.;

                                if (options.delayedUseSatisfaction == 1) {
                                    // The adaption of NAg shall be done on a daily basis, hence the difference from
                                    // prev timestep is calculated. It can happen that secondCell is used in one
                                    // timestep and not in another thus both has to be subtracted.
                                    G_dailyRemainingUse[n] = (remainingUse -
                                                              G_PrevTotalUnstatisfiedUse[n] -
                                                              G_PrevUnsatAllocUse[n]);
                                    G_totalUnsatisfiedUse[n] = remainingUse;

                                } else {    // if options.delayedUseSatisfaction == 0: G_totalUnsatisfiedUse[n] is only stored as output variable
                                    G_totalUnsatisfiedUse[n] += remainingUse;
                                    // alloc usage but no second cell identified
                                    // G_dailyRemainingUse is assigned to  adapt NAg in next timestep
                                    G_dailyRemainingUse[n] = remainingUse;
                                }


                            }

                            else {    // secondCell > 0: remainingUse is always > 0. in this case

                                // G_AllocatedUse[secondCell] is increased (+= instead of =),
                                // as it may be identified as a "second cell" for other cells in the same timestep.
                                // When the remainingUse of the second cell is computed at the beginning of each
                                // timestep, G_AllocatedUse is set to zero directly after adding the value to the
                                // remainingUse.

                                G_AllocatedUse[secondCell] += remainingUse;

                                // WG22b: G_AllocUseToNeigborCell and G_SecondCell used to redistribute the unsatisfied
                                // fractions of allocated use to the original cells:

                                // If G_AllocatedUse cannot be fully satisfied in the second cell, the remainin
                                // allocated use is redistributed over the neighboring cells using
                                // G_AllocUseToNeigborCell[n] and G_SecondCell[n].
                                G_AllocUseToNeigborCell[n] = remainingUse;

                                // WG22b: G_daily_allocatedUseNextDay only used for output G_ALLOC_USE_NEXT_DAY_%d.365.UNF0:
                                //	If G_routingCell[n] > G_routingCell[secondCell],
                                //	remainingUse allocated to a nbc will be satisfied in the next timestep.
                                //	To check if water use components are written correctly, this unsatisfied allocated use
                                //	must be taken into account.
                                if (G_CellPositionRoutOrder[secondCell] <
                                    G_CellPositionRoutOrder[n])    // (Only G_AllocatedUse[secondCell] of the next time step is taken into account)
                                    G_daily_allocatedUseNextDay[secondCell] += remainingUse;

                                if (options.delayedUseSatisfaction == 1) {
                                    G_totalUnsatisfiedUse[n] = 0.;
                                    // Otherwise, remainingUse of cell n would be increased by G_totalUnsatisfiedUse[n] from the previous timestep again.
                                }
                                // else (options.delayedUseSatisfaction == 0):
                                // G_totalUnsatisfiedUse[n] is not changed in this timestep.

                            }

                            G_SecondCell[n] = secondCell;

                        }

                        else {    // (options.use_alloc == 0)	(remainingUse can be >= 0 at this point)

                            if (options.delayedUseSatisfaction == 1) {
                                // NAg adaption happens on a daily basis thus it is compared to the
                                // G_totalUnsatisfiedUse from last time step before assigning remainingUse.
                                G_dailyRemainingUse[n] = remainingUse - G_PrevTotalUnstatisfiedUse[n];
                                G_totalUnsatisfiedUse[n] = remainingUse;
                            }
                            else {
                                // In this case, 'G_totalUnsatisfiedUse' is not added to 'remainingUse' at the beginning,
                                // but the sum of 'G_totalUnsatisfiedUse' is stored as an output variable.
                                G_totalUnsatisfiedUse[n] += remainingUse;
                                //NAg adaption based on remainingUse
                                G_dailyRemainingUse[n] = remainingUse;
                            }
                        }
                        // together with G_dailyRemainingUse G_withdrawalIrrigFromSWB and
                        // G_consumptiveUseIrrigFromSWB are used to calculate correct returnflow reduction in
                        // case of allocation this information must origin from the timestep, when it is
                        // allocated, thus it is assigned here
                        G_withdrawalIrrigFromSwb[n] = G_monthlyWUIrrigFromSwb(n, month) / (1000000000. * (double) numberOfDaysInMonth[month]);
                        G_consumptiveUseIrrigFromSwb[n] = G_monthlyCUIrrigFromSwb(n, month) / (1000000000. * (double) numberOfDaysInMonth[month]);
                    }

                    // Computation of SWB and land area fractions at the end of 'routing'.

                    //                 riverAvail = transportedVolume; // CR 2015-08-04 not used in WG22b

                    // WG22b: The river availability represents the water resource of a cell (including inland sinks).
                    // Thus, the outflow of SWB should be taken into account in inland sinks as well.

                    if (1 == options.statcorrOpt)
                        // if station correction is used, multiply discharge with station correction (e.g. for calibration).
                        transportedVolume *= G_statCorrFact[n];

                    // WG22b: New approach to compute CellRunoff:
                    // in inland sinks, transportedVolume will be added to AET later, hence CellRunoff gets negative.
                    if (G_LDD[n] < 0) {
                        CellRunoff = 0. - inflowUpstream;
                        G_transportedVolume_for_AET[n] = transportedVolume;
                    } else {
                        CellRunoff = (transportedVolume - inflowUpstream);
                    }


                    //-------------------------------------for daily velocity file----------------------
                    if (options.day_store == 1) {
                        for (b = 0; b <= nSpecBasins - 1; b++) {
                            if (n == cbasin.cellNum[b] - 1) {
                                dailyRiverVelocity(b,day - 1) += (riverVelocity /
                                                                   86.4); // //km/d --> m/s ((60*60*24)/1000)
                                break;
                            }
                        }
                    }

                    // ... and put the water into the downstream cell
                    if ((G_downstreamCell[n] - 1) >= 0) {
                        if (0 != G_toBeCalculated[(G_downstreamCell[n] - 1)])
                            G_riverInflow[(G_downstreamCell[n] - 1)] += transportedVolume;
                    }

                    // annualCellRunoff over the whole year
                    G_AnnualCellRunoff[n] += CellRunoff;

                    // monthly output (possibly combined with daily output) start
                    if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {
                        // G_monthlyCellRunoff should be defined as the inflow from global lakes,
                        // when calculating total inflow to reservoirs.

                        // use precip on land, wetlands and local lakes and global lakes and reservoirs if they occur ans all only on continental area in km3
                        // Only continental cells
                        // each other cell that contains continental area
                        // ATTENTION/POTENTIAL ERROR: G_dailyRiverPrecip[n] is not taken into account in calculationg G_landWaterExclGloLakAreaFrac[n]

                        if (geo.G_contcell[n]) {
                            G_monthlyConsistentPrecip(n,month) += (
                                           (dailyWaterBalance.G_openWaterPrec[n] * G_landWaterExclGloLakAreaFrac[n] /
                                            geo.G_contfreq[n]) *
                                           ((geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) / 1000000.) +
                                           G_GloLakePrecip[n] + G_ReservoirPrecip[n] + G_dailyRiverPrecip[n]);

                            // Only for runs which consider glaciers
                            if ((options.glacierOpt == 1) && (G_glacAdaptedArea[n] > 0.)) // for glacierized cells only
                                // Adapt consistent precipitation by assigning precipitation from GGM forcing
                                // over glacier area fraction instead of precipitation from WGHM forcing
                                G_monthlyConsistentPrecip(n,month) =
                                               G_monthlyConsistentPrecip(n,month) - (dailyWaterBalance.G_openWaterPrec[n] * G_glacAdaptedArea[n] / 1000000.) + G_glacPrecip[n];
                        }
                            // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                            // inner Caspian lake cell (not belonging to a continent)
                        else {
                            G_monthlyConsistentPrecip(n,month) = 0.;
                        }

                        // add cell runoff to monthly cell runoff
                        G_monthlyCellRunoff(n,month) += CellRunoff;

                        // Only continental cells
                        // each other cell that contains continental area
                        if (geo.G_contcell[n]) {
                            G_monthlyCellSurfaceRunoff(n,month) += ((CellRunoff - G_localGWRunoff[n])
                                                                     / (cellArea * (geo.G_contfreq[n] / 100.)) *
                                                                     1000000.);
                        }
                            // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                            // inner Caspian lake cell (not belonging to a continent)
                        else {
                            G_monthlyCellSurfaceRunoff(n,month) = 0.;
                        }

                        // add transportedVolume to monthly cell runoff (only outside of inland sinks)
                        if (G_LDD[n] >= 0)
                            G_monthlyRiverAvail(n,month) += transportedVolume;

                        if (options.outRiverInUpstream)
                            G_monthlyRiverInUpstream(n,month) += inflowUpstream;


                        // add riverVelocity to monthlyVelocity
                        G_monthlyVelocity(n,month) += (riverVelocity / 86.4); // velocity grid [m/s]
                        // for those cells which are defined as inland sinks and do have considerable lake/wetland area (>'0'),
                        // river availability should be updated, otherwise river availability on basin scale could be '0'

                        // WG22b: The river availability represents the water resource of a cell.
                        // Thus, the outflow of SWB should be taken into account in inland sinks as well.

                        // monthly min and max river availability output start
                        //most certainly to be rationalised
                        if (options.outMinMaxRiverAvail) {
                            if (G_LDD[n] >= 0) {
                                if (day_in_month == 1) {
                                    G_monthlyMinRiverAvail(n,month) = transportedVolume;
                                    G_monthlyMaxRiverAvail(n,month) = transportedVolume;
                                } else {
                                    if (transportedVolume > G_monthlyMaxRiverAvail(n,month))
                                        G_monthlyMaxRiverAvail(n,month) = transportedVolume;
                                    if (transportedVolume < G_monthlyMinRiverAvail(n,month))
                                        G_monthlyMinRiverAvail(n,month) = transportedVolume;
                                }
                            }
                        }
                        // monthly min and max river availability output end
                        if (options.outRiverPET == 1) G_monthlyRiverAreaFrac(n,month) += G_riverAreaFrac[n];
                        if (options.outRiverPET == 1)
                            G_monthlyRiverPET(n,month) += G_dailyRiverEvapo[n] * 1000000. /
                                                           (G_riverAreaFrac[n] / 100. * cellArea); //km3 > mm
                        //output for reservoir and storages start
                        if ((options.antNatOpt == 0) && (options.resOpt == 1)) {
                            G_monthlyResInflow(n,month) += G_daily_res_in[n]; // total inflow in month [m3/timeStepsPerDay] -> [m3/month]
                            G_monthlyResOutflow(n,month) += G_daily_res_out[n]; // total outflow in month [m3/timeStepsPerDay] -> [m3/month]

                            if (options.grid_store_TypeForStorages == 1)
                                G_monthlySurfStor(n,month) +=
                                               (G_locLakeStorage[n]
                                                + G_locWetlStorage[n]
                                                + G_gloLakeStorage[n]
                                                + G_gloWetlStorage[n]
                                                + G_riverStorage[n]
                                                + G_gloResStorage[n]);

                        }
                        else {
                            // mean value of total surface water storage in grid cell (lakes/wetlands/rivers) [km3] (converted to mm for optional output in a later step)
                            if (options.grid_store_TypeForStorages == 1)
                                G_monthlySurfStor(n,month) +=
                                               (G_locLakeStorage[n]
                                                + G_locWetlStorage[n]
                                                + G_gloLakeStorage[n]
                                                + G_gloWetlStorage[n]
                                                + G_riverStorage[n]);
                        }

                        if (options.outGwrSwb == 1)
                            G_monthlyGwrSwb(n,month) += G_gwr_swb[n];

                        G_monthlyLocWetlExtent(n,month) += G_locwetlextent[n];
                        G_monthlyGloWetlExtent(n,month) += G_glowetlextent[n];

                        //output of monthly groundwater recharge below surface water bodies in a gridcell
                        if (options.outGwrunSurfrun == 1)
                            G_monthlyGwrunSurfrun(n,month) += groundwater_runoff_km3 + dailyLocalSurfaceRunoff; // in km3 day-1

                    }

                    // daily output option 31 start
                    if ((3 == options.grid_store) || (4 == options.grid_store)) {

                        if (options.outPrecDaily) {
                            // Only continental cells
                            if (geo.G_contcell[n]) {
                                G_daily31ConsistentPrecip(n,day_in_month - 1) = (
                                               (dailyWaterBalance.G_openWaterPrec[n] *
                                                G_landWaterExclGloLakAreaFrac[n] / geo.G_contfreq[n]) *
                                               ((geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) / 1000000.) +
                                               G_GloLakePrecip[n] + G_ReservoirPrecip[n] + G_dailyRiverPrecip[n]);

                                // Only for runs which consider glaciers
                                if ((options.glacierOpt == 1) && (G_glacAdaptedArea[n] > 0.)) // for glacierized cells only
                                    // Adapt consistent precipitation by assigning precipitation from GGM forcing
                                    // over glacier area fraction instead of precipitation from WGHM forcing
                                    G_daily31ConsistentPrecip(n,day_in_month - 1) =
                                                   G_daily31ConsistentPrecip(n,day_in_month - 1) - (dailyWaterBalance.G_openWaterPrec[n] * G_glacAdaptedArea[n] / 1000000.) + G_glacPrecip[n];
                            }
                                // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                            else {
                                G_daily31ConsistentPrecip(n,day_in_month - 1) = 0.;
                            }
                        }

                        if (options.outCellRunoffDaily)
                            G_daily31CellRunoff(n,day_in_month - 1) = CellRunoff;

                        if (options.outCellSurfaceDaily) {
                            // Only continental cells
                            if (geo.G_contcell[n]) {
                                G_daily31CellSurfaceRunoff(n,day_in_month - 1) =
                                               (CellRunoff - G_localGWRunoff[n])
                                               / (cellArea * (geo.G_contfreq[n] / 100.)) * 1000000.;
                            }
                                // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                            else {
                                G_daily31CellSurfaceRunoff(n,day_in_month - 1) = 0.;
                            }
                        }

                        if ((options.outRiverAvailDaily) && (G_LDD[n] >= 0)){
                            G_daily31RiverAvail(n,day_in_month - 1) = transportedVolume;

                            // for those cells which are defined as inland sinks and do have considerable lake/wetland area (>'0'),
                            // river availability should be updated, otherwise river availability on basin scale could be '0'
                        }

                        if (options.outRiverVeloDaily)
                            G_daily31Velocity(n,day_in_month - 1) = (riverVelocity /
                                                                      86.4); // velocity grid [m/s]

                        if (options.outGwrSwbDaily)
                            G_daily31GwrSwb(n,day_in_month - 1) = G_gwr_swb[n]; // daily gwr below surface water bodies [mm]
                        if (options.outGWRunoffDaily)
                            // daily fraction of surface water bodies [mm]
                            G_daily31GwRunoff(n,day_in_month -1) = groundwater_runoff_km3; // daily fraction of surface water bodies [mm]
                        if (options.outGwrunSurfrunDaily)
                            G_daily31GwrunSurfrun(n,day_in_month - 1) =
                                           groundwater_runoff_km3 + dailyLocalSurfaceRunoff; // in km3 day-1
                        if (options.outLandAreaFracDaily)
                            G_daily31LandAreaFrac(n,day_in_month - 1) = G_landAreaFrac[n];
                        if (options.outCellAETWCaDaily) // add NAg and actual use SW and later AET
                            G_daily31CellAETWCa(n,day_in_month - 1) =
                                           dailyActualUse + G_dailydailyNUg[n]; //already in km3 day-1
                    }

                    // daily output option 365 start
                    if ((5 == options.grid_store) || (6 == options.grid_store)) {

                        if (options.outPrecDaily || ((options.scoutcPrecip) && (2 == options.day_store))) {
                            // Only continental cells
                            if (geo.G_contcell[n]) {
                                G_daily365ConsistentPrecip(n,day - 1) = ((dailyWaterBalance.G_openWaterPrec[n] *
                                                                           G_landWaterExclGloLakAreaFrac[n] /
                                                                           geo.G_contfreq[n]) *
                                                                          ((geo.areaOfCellByArrayPos(n) *
                                                                            (geo.G_contfreq[n] / 100.)) / 1000000.) +
                                                                          G_GloLakePrecip[n] + G_ReservoirPrecip[n] +
                                                                          G_dailyRiverPrecip[n]);

                                // Only for runs which consider glaciers
                                if ((options.glacierOpt == 1) && (G_glacAdaptedArea[n] > 0.)) // for glacierized cells only
                                    // Adapt consistent precipitation by assigning precipitation from GGM forcing
                                    // over glacier area fraction instead of precipitation from WGHM forcing
                                    G_daily365ConsistentPrecip(n,day - 1) =
                                                   G_daily365ConsistentPrecip(n,day - 1) - (dailyWaterBalance.G_openWaterPrec[n] * G_glacAdaptedArea[n] / 1000000.) + G_glacPrecip[n];

                            }
                                // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                            else {
                                G_daily365ConsistentPrecip(n,day - 1) = 0.;
                            }
                        }

                        if ((options.outGwrSwbDaily) || ((options.scoutGwrSwb) && (2 == options.day_store)))
                            G_daily365GwrSwb(n,day - 1) = G_gwr_swb[n];
                        if ((options.outGWRunoffDaily) || ((options.scoutGwRunoff) && (2 == options.day_store)))
                            //HMS 2017-02-03 the following line was included but makes no sense, or?
                            G_daily365GwRunoff(n,day - 1) = groundwater_runoff_km3;
                        if (options.outGwrunSurfrunDaily)
                            G_daily365GwrunSurfrun(n,day - 1) = groundwater_runoff_km3 + dailyLocalSurfaceRunoff;
                        if ((options.outLandAreaFracDaily) || ((options.scoutLandAreaFrac) && (2 == options.day_store)))
                            G_daily365LandAreaFrac(n,day - 1) = G_landAreaFrac[n];

                        if (options.outCellRunoffDaily ||
                            (((options.scoutCellRunoffkm3) || (options.scoutCellRunoffmm)) &&
                             (2 == options.day_store))) {
                            G_daily365CellRunoff(n,day - 1) = CellRunoff;
                            // Only continental cells
                            if (geo.G_contcell[n]) {
                                G_daily365CellRunoff_mm(n,day - 1) =
                                               CellRunoff / (cellArea * (geo.G_contfreq[n] / 100.)) *
                                               1000000.;
                            }
                                // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                            else {
                                G_daily365CellRunoff_mm(n,day - 1) = 0.;
                            }
                        }

                        if (options.outCellSurfaceDaily ||
                            ((options.scoutCellSRunoff) && (2 == options.day_store))) {
                            // Only continental cells
                            // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                            if (geo.G_contcell[n]) {
                                G_daily365CellSurfaceRunoff(n,day - 1) = (CellRunoff - G_localGWRunoff[n]) /
                                                                          (cellArea * (geo.G_contfreq[n] / 100.)) *
                                                                          1000000.;
                            }
                            // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                            else {
                                G_daily365CellSurfaceRunoff(n,day - 1) = 0.;
                            }
                        }

                        if ((options.outRiverAvailDaily || ((options.scoutQ) && (2 == options.day_store))) && (G_LDD[n] >= 0)) {
                            G_daily365RiverAvail(n,day - 1) = transportedVolume;
                        }
                        if (options.outCellAETWCaDaily)
                            G_daily365CellAETWCa(n,day - 1) =
                                           dailyActualUse + G_dailydailyNUg[n]; //already in km3 day-1

                        if ((options.outRiverVeloDaily) || ((options.scoutFlowVelo) && (2 == options.day_store)))
                            G_daily365Velocity(n,day - 1) = (riverVelocity /
                                                              86.4); // velocity grid [m/s]

                    }

                    for (b = 0; b <= nSpecBasins - 1; b++) {
                        if (n == cbasin.cellNum[b] - 1) {
                            if (G_LDD[n] >= 0)
                                dailyRiverDischarge(b,day - 1) += transportedVolume;
                            break;
                        }
                    }
                }
                //--------------------------------------------------------------------------
                // calculation of water temperature 1 without power plant 2 with power plant
                if (options.calc_wtemp == 1 || options.calc_wtemp == 2) {

                    if(options.calc_wtemp == 2) {
                        cwt.tsf.calcTempPP = true;
                        cwt.tsf.powerPlantWW = G_yearlyWW[n] / 1.e9 / 365.; // [m3/a] --> [km3/d]
                        cwt.tsf.powerPlantWC = G_yearlyWC[n] / 1.e9 / 365.; // [m3/a] --> [km3/d]
                    }
                    else
                        cwt.tsf.calcTempPP = false;

                    for (int i = 0; i < 9; i++) {
                        if (G_inflow_cells(n, i) != 0) {
                            cwt.tsf.cellInflows.push_back(G_RiverAvailPreStep[G_inflow_cells(n,i)-1]);       // Inflow previous step
                            cwt.tsf.cellInflowsTemp.push_back(G_RiverWTempPreStep[G_inflow_cells(n,i)-1]);   // temperature previous step
                        }
                        else if(i==4) { // current cell
                            cwt.tsf.cellInflows.push_back(0.);
                            cwt.tsf.cellInflowsTemp.push_back(G_RiverWTempPreStep[n]);    // for river calculation temperature prev day needed
                        }
                        else { // no inflow from this direction
                            cwt.tsf.cellInflows.push_back(0.);
                            cwt.tsf.cellInflowsTemp.push_back(15.);    // temperature irrelevant bc. no associated volume; only to write something into this place
                        }
                    }

                    // ice thickness from prev day
                    cwt.tsf.iceThickness[0] = G_iceThicknessRiver[n];
                    cwt.tsf.iceThickness[1] = G_iceThicknesslocLake[n];
                    cwt.tsf.iceThickness[2] = G_iceThicknesslocWetland[n];
                    cwt.tsf.iceThickness[3] = G_iceThicknessGloLake[n];
                    cwt.tsf.iceThickness[4] = G_iceThicknessReservoir[n];
                    cwt.tsf.iceThickness[5] = G_iceThicknessGloWetland[n];

                    cwt.tsf.ShortRad = climate.G_shortwave_d(n, day_in_month - 1);  // [W/m2]
                    cwt.tsf.LongRad = climate.G_longwave_d(n, day_in_month - 1);    // [W/m2]
                    cwt.tsf.AirTempC = climate.G_temperature_d(n, day_in_month - 1); // [degC]
                    if (climate.G_GW_Temperature_y[n] < 0.) // groundwater cannot be below 0 degC
                        cwt.tsf.GWTemp = 0.;    // [degC]
                    else
                        cwt.tsf.GWTemp = climate.G_GW_Temperature_y[n];     // [degC]

                    cwt.tsf.CellCorrFact = dailyWaterBalance.G_cellCorrFact[n];
                    cwt.tsf.openWaterPET = dailyWaterBalance.G_openWaterPET[n];     // [mm]
                    cwt.tsf.openWaterPrecip = dailyWaterBalance.G_openWaterPrec[n]; // [mm]
                    cwt.tsf.CellArea = cellArea;    // [km2]

                    // locLake
                    cwt.tsf.locLakeArea = G_loc_lake[n];    //[%]
                    cwt.tsf.locLakeStorPrevStep = G_locLakeStoragePrevStep + (G_loc_lake[n] / 100.) * cellArea * G_lakeDepthActive[n]; // + maxStorage so no neg. vol. occurs
                    cwt.tsf.locLakeTempPrevStep = G_locLakeTempPrevStep[n];
                    cwt.tsf.dailyLocLakeEvapo = G_dailyLocLakeEvapo[n];
                    cwt.tsf.locLakeAreaRedFactor = locLakeAreaReductionFactor;
                    cwt.tsf.localRunoff = G_localRunoff[n];
                    cwt.tsf.localGWRunoff = G_localGWRunoff[n];

                    // locWetland
                    cwt.tsf.locWetlandArea = G_loc_wetland[n];  //[%]
                    cwt.tsf.locWetlandStorPrevStep = G_locWetlStoragePrevStep;
                    cwt.tsf.locWetlandTempPrevStep = G_locWetlTempPrevStep[n];
                    cwt.tsf.dailyLocWetlandEvapo = G_dailyLocWetlEvapo[n];
                    cwt.tsf.locWetlandAreaRedFactor = locWetlAreaReductionFactor;
                    cwt.tsf.inflowLocWetland = outflowlocLake;
                    cwt.tsf.outflowLocWetland = outflowlocWetl;

                    // gloLake
                    cwt.tsf.gloLakeArea = G_lake_area[n];   //[km2]
                    cwt.tsf.gloLakeStorPrevStep = G_gloLakeStoragePrevStep + ((double) G_lake_area[n]) * G_lakeDepthActive[n];  // + maxStorage so no neg. volume occurs
                    cwt.tsf.gloLakeTempPrevStep = G_gloLakeTempPrevStep[n];
                    cwt.tsf.dailyGloLakeEvapo = G_dailyGloLakeEvapo[n];
                    // total river inflows from different cells
                    cwt.tsf.inflowGloLake = 0.;
                    for (int i = 0; i < 9; i++) {
                        cwt.tsf.inflowGloLake += cwt.tsf.cellInflows[i];
                    }
                    if (G_loc_wetland[n] > 0.)                      // add outflow from locWetland or locLake if exists
                        cwt.tsf.inflowGloLake += outflowlocWetl;
                    else if (G_loc_lake[n] > 0.)
                        cwt.tsf.inflowGloLake += outflowlocLake;

                    // reservoir
                    if (G_lake_area[n] > 0.)                        // for correct inflow volume if no gloLake
                        cwt.tsf.inflowReservoir = outflowGloLake;
                    else
                        cwt.tsf.inflowReservoir = cwt.tsf.inflowGloLake;
                    cwt.tsf.reservoirArea = G_reservoir_area[n];    //[km2]
                    cwt.tsf.reservoirStorPrevStep = G_gloResStoragePrevStep;
                    cwt.tsf.reservoirTempPrevStep = G_reservoirTempPrevStep[n];
                    cwt.tsf.dailyReservoirEvapo = G_dailyResEvapo[n];

                    // gloWetlands
                    if (G_reservoir_area[n] > 0.)              // for correct inflow for different conditions of existence of lakes,... in cell
                        cwt.tsf.gloWetlandInflow = G_daily_res_out[n];
                    else if (G_lake_area[n] > 0.)
                        cwt.tsf.gloWetlandInflow = outflowGloLake;
                    else
                        cwt.tsf.gloWetlandInflow = cwt.tsf.inflowGloLake;
                    cwt.tsf.gloWetlandArea = G_glo_wetland[n];  //[%]
                    cwt.tsf.gloWetlandStorPrevStep = G_gloWetlStoragePrevStep;
                    cwt.tsf.gloWetlandTempPrevStep = G_gloWetlandTempPrevStep[n];
                    cwt.tsf.dailyGloWetlandEvapo = G_dailyGloWetlEvapo[n];
                    cwt.tsf.gloWetlandAreaReductionFactor = gloWetlAreaReductionFactor;

                    // River (if River Evapo is implemented in WaterGAP (currently not as of May 2020) handover here for temperature computation)
                    if (G_glo_wetland[n] > 0.)                 // for correct inflow for different conditions of existence of lakes,... in cell
                        cwt.tsf.RiverInflow = outflowGloWetland;
                    else if (G_reservoir_area[n] > 0.)
                        cwt.tsf.RiverInflow  = G_daily_res_out[n];
                    else if (G_lake_area[n] > 0.)
                        cwt.tsf.RiverInflow  = outflowGloLake;
                    else
                        cwt.tsf.RiverInflow  = cwt.tsf.inflowGloLake;
                    cwt.tsf.RiverStor = G_riverStorage[n];
                    cwt.tsf.localGWRunoffIntoRiver = G_localGWRunoffIntoRiver[n];
                    cwt.tsf.localRunoffIntoRiver = G_localRunoffIntoRiver[n];
                    cwt.tsf.RiverBottomWidth = G_riverBottomWidth[n];
                    cwt.tsf.RiverLength = G_riverLength[n];
                    cwt.tsf.RiverLengthCheck = G_riverLengthCheck[n];       // needed bc. G_riverLength is changed after read in

                    if (G_downstreamCell[n] > 0)                // check if cell has outflow;
                        cwt.tsf.hasOutflow = true;
                    else
                        cwt.tsf.hasOutflow = false;

                    // [0] outflow/river, [1] locLake, [2] locWetland, [3] gloLake, [4] reservoir, [5] gloWetland; if non existent in cell -9999. as temp
                    // ice thickness [6] river, [7] locLake, [8] locWetland, [9] gloLake, [10] reservoir, [11] gloWetland
                    calcWTempReturns = cwt.calcWTemp(cwt.tsf);

                    G_iceThicknessRiver[n] = calcWTempReturns[6];
                    G_iceThicknesslocLake[n] = calcWTempReturns[7];
                    G_iceThicknesslocWetland[n] = calcWTempReturns[8];
                    G_iceThicknessGloLake[n] = calcWTempReturns[9];
                    G_iceThicknessReservoir[n] = calcWTempReturns[10];
                    G_iceThicknessGloWetland[n] = calcWTempReturns[11];

                    if(options.outWaterTempDaily)
                        G_daily31WaterTemp(n, day_in_month - 1) = calcWTempReturns[0];      // daily31 output

                    if (options.outWaterTempMonthly) {
                        if (calcWTempReturns[0] >= 0.) {
                            SumWaterTempRiver[n] += calcWTempReturns[0];                              // for calculation of monthly mean Water Temp
                            NumWaterTempRiver[n] += 1.;                                         // needed bc. of possibility of -9999. (no calculated temperature)
                        }
                        if (day_in_month == numberOfDaysInMonth[month] && NumWaterTempRiver[n] > 0.) {
                            G_monthlyMeanWaterTempRiver(n, month) = SumWaterTempRiver[n] / NumWaterTempRiver[n];    // monthly output
                            SumWaterTempRiver[n] = 0.;
                            NumWaterTempRiver[n] = 0.;
                        } else
                            G_monthlyMeanWaterTempRiver(n, month) = -9999.;
                    }

                    // output of all SWB separately as daily31
                    if (options.outWaterTempDailyAllSWB) {
                        G_daily31locLakeTemp(n, day_in_month - 1) = calcWTempReturns[1];
                        G_daily31locWetlTemp(n, day_in_month - 1) = calcWTempReturns[2];
                        G_daily31gloLakeTemp(n, day_in_month - 1) = calcWTempReturns[3];
                        G_daily31reservoirTemp(n, day_in_month - 1) = calcWTempReturns[4];
                        G_daily31gloWetlandTemp(n, day_in_month - 1) = calcWTempReturns[5];
                    }

                    // output of all SWB separately as monthly
                    if(options.outWaterTempMonthlyAllSWB){
                        // locLake
                        if (calcWTempReturns[1] >= 0.) {
                            SumWaterTempLocLake[n] += calcWTempReturns[1];                              // for calculation of monthly mean Water Temp locLake
                            NumWaterTempLocLake[n] += 1.;                                         // needed bc. of possibility of -9999. (no calculated temperature)
                        }
                        if (day_in_month == numberOfDaysInMonth[month] && NumWaterTempLocLake[n] > 0.) {
                            G_monthlyMeanWaterTempLocLake(n, month) = SumWaterTempLocLake[n] / NumWaterTempLocLake[n];    // monthly output
                            SumWaterTempLocLake[n] = 0.;
                            NumWaterTempLocLake[n] = 0.;
                        } else
                            G_monthlyMeanWaterTempLocLake(n, month) = -9999.;
                        // locWetland
                        if (calcWTempReturns[2] >= 0.) {
                            SumWaterTempLocWetl[n] += calcWTempReturns[2];                              // for calculation of monthly mean Water Temp locWetland
                            NumWaterTempLocWetl[n] += 1.;                                         // needed bc. of possibility of -9999. (no calculated temperature)
                        }
                        if (day_in_month == numberOfDaysInMonth[month] && NumWaterTempLocWetl[n] > 0.) {
                            G_monthlyMeanWaterTempLocWetl(n, month) = SumWaterTempLocWetl[n] / NumWaterTempLocWetl[n];    // monthly output
                            SumWaterTempLocWetl[n] = 0.;
                            NumWaterTempLocWetl[n] = 0.;
                        } else
                            G_monthlyMeanWaterTempLocWetl(n, month) = -9999.;
                        // gloLake
                        if (calcWTempReturns[3] >= 0.) {
                            SumWaterTempGloLake[n] += calcWTempReturns[3];                              // for calculation of monthly mean Water Temp gloLake
                            NumWaterTempGloLake[n] += 1.;                                         // needed bc. of possibility of -9999. (no calculated temperature)
                        }
                        if (day_in_month == numberOfDaysInMonth[month] && NumWaterTempGloLake[n] > 0.) {
                            G_monthlyMeanWaterTempGloLake(n, month) = SumWaterTempGloLake[n] / NumWaterTempGloLake[n];    // monthly output
                            SumWaterTempGloLake[n] = 0.;
                            NumWaterTempGloLake[n] = 0.;
                        } else
                            G_monthlyMeanWaterTempGloLake(n, month) = -9999.;
                        // reservoir
                        if (calcWTempReturns[4] >= 0.) {
                            SumWaterTempReservoir[n] += calcWTempReturns[4];                              // for calculation of monthly mean Water Temp reservoir
                            NumWaterTempReservoir[n] += 1.;                                         // needed bc. of possibility of -9999. (no calculated temperature)
                        }
                        if (day_in_month == numberOfDaysInMonth[month] && NumWaterTempReservoir[n] > 0.) {
                            G_monthlyMeanWaterTempReservoir(n, month) = SumWaterTempReservoir[n] / NumWaterTempReservoir[n];    // monthly output
                            SumWaterTempReservoir[n] = 0.;
                            NumWaterTempReservoir[n] = 0.;
                        } else
                            G_monthlyMeanWaterTempReservoir(n, month) = -9999.;
                        // gloWetland
                        if (calcWTempReturns[5] >= 0.) {
                            SumWaterTempGloWetl[n] += calcWTempReturns[5];                              // for calculation of monthly mean Water Temp gloWetland
                            NumWaterTempGloWetl[n] += 1.;                                         // needed bc. of possibility of -9999. (no calculated temperature)
                        }
                        if (day_in_month == numberOfDaysInMonth[month] && NumWaterTempGloWetl[n] > 0.) {
                            G_monthlyMeanWaterTempGloWetl(n, month) = SumWaterTempGloWetl[n] / NumWaterTempGloWetl[n];    // monthly output
                            SumWaterTempGloWetl[n] = 0.;
                            NumWaterTempGloWetl[n] = 0.;
                        } else
                            G_monthlyMeanWaterTempGloWetl(n, month) = -9999.;
                    }

                    // water temperature save for next day, if not calculated (-9999.) then reset to air temp
                    if (calcWTempReturns[0] == -9999.)
                        if (cwt.tsf.AirTempC < 0.)
                            G_RiverWTempPreStep[n] = 0.;
                        else
                            G_RiverWTempPreStep[n] = cwt.tsf.AirTempC;                     // river
                    else
                        G_RiverWTempPreStep[n] = calcWTempReturns[0];

                    if (calcWTempReturns[1] == -9999.)
                        if (cwt.tsf.AirTempC < 0.)
                            G_locLakeTempPrevStep[n] = 0.;
                        else
                            G_locLakeTempPrevStep[n] = cwt.tsf.AirTempC;                  // locLake
                    else
                        G_locLakeTempPrevStep[n] = calcWTempReturns[1];

                    if (calcWTempReturns[2] == -9999.)
                        if (cwt.tsf.AirTempC < 0.)
                            G_locWetlTempPrevStep[n] = 0.;
                        else
                            G_locWetlTempPrevStep[n] = cwt.tsf.AirTempC;                  // locWetland
                    else
                        G_locWetlTempPrevStep[n] = calcWTempReturns[2];

                    if (calcWTempReturns[3] == -9999.)
                        if (cwt.tsf.AirTempC < 0.)
                            G_gloLakeTempPrevStep[n] = 0.;
                        else
                            G_gloLakeTempPrevStep[n] = cwt.tsf.AirTempC;                  // gloLake
                    else
                        G_gloLakeTempPrevStep[n] = calcWTempReturns[3];

                    if (calcWTempReturns[4] == -9999.)
                        if (cwt.tsf.AirTempC < 0.)
                            G_reservoirTempPrevStep[n] = 0.;
                        else
                            G_reservoirTempPrevStep[n] = cwt.tsf.AirTempC;                // reservoir
                    else
                        G_reservoirTempPrevStep[n] = calcWTempReturns[4];

                    if (calcWTempReturns[5] == -9999.)
                        if (cwt.tsf.AirTempC < 0.)
                            G_gloWetlandTempPrevStep[n] = 0.;
                        else
                            G_gloWetlandTempPrevStep[n] = cwt.tsf.AirTempC;               // gloWetland
                    else
                        G_gloWetlandTempPrevStep[n] = calcWTempReturns[5];

                    G_RiverAvailPreStep[n] = transportedVolume;                // cell outflow [km3]

                    cwt.tsf.cellInflows.clear();        // reset of vectors
                    cwt.tsf.cellInflowsTemp.clear();
                }
                // end of water temperature calculation
                //-------------------------------------

                G_riverInflow[n] = 0.;
            }
        }
    }


    // Calculate total actual Evapotranspiration (land and open water areas); consistent with corrected runoff (2.1f)
    // Store storage values and end of timestep
    if (options.grid_store > 1) {

        double dailyTotalAET;
        double dailyTotalAET_km3;
        double dailyOpenWaterEvap;
        double dailyOpenWaterEvap_km3;

        for (int n = 0; n < ng; n++) {
            if (0 != G_toBeCalculated[n]) {
                dailyTotalAET = 0.;
                dailyTotalAET_km3 = 0.;
                dailyOpenWaterEvap = 0.;
                dailyOpenWaterEvap_km3 = 0.;
                cellArea = geo.areaOfCellByArrayPos(n);


                // first store only land area evapotranspiration
                dailyTotalAET = dailyWaterBalance.getDailyCellLandAET(n);
                // then sum up land area evapotranspiration and evaporation
                // from open water areas
                // all areas in [km2]
                // resulting data unit: mm * km2

                // G_dailyLocLakeEvapo[n] etc. = openWaterPET * area reduction factor -> SWB area already reduced;
                // therefore, G_loc_lake[n] etc. is used in the following equation:
                dailyOpenWaterEvap_km3 =    // in km3
                               (G_dailyLocLakeEvapo[n] / 1000000. * cellArea * G_loc_lake[n] / 100.)
                               + (G_dailyLocWetlEvapo[n] / 1000000. * cellArea * G_loc_wetland[n] / 100.)
                               + (G_dailyGloLakeEvapo[n] / 1000000. * G_lake_area[n])
                               + (G_dailyResEvapo[n] / 1000000. * G_reservoir_area[n])
                               + (G_dailyGloWetlEvapo[n] / 1000000. * cellArea * G_glo_wetland[n] / 100.);

                //add transported volume to AET in inland sinks to close the water balance
                if (G_LDD[n] < 0) {
                    dailyOpenWaterEvap_km3 += G_transportedVolume_for_AET[n];
                }

                if (options.riverEvapoOpt == 1) {
                    // POSSIBLE ERROR data unit of dailyOpenWaterEvap: mm * km2 does not correspond to km3 of G_dailyRiverEvapo[n];
                    dailyOpenWaterEvap_km3 += G_dailyRiverEvapo[n];        // G_dailyRiverEvapo[n] in km3
                }

                dailyTotalAET_km3 = (dailyTotalAET / 1000000. * cellArea * G_landAreaFrac[n] /
                                     100.)        // dailyTotalAET: mm -> km3
                                    + dailyOpenWaterEvap_km3;

                // monthly output
                if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {
                    // Only continental cells
                    if (geo.G_contcell[n]) {

                        // AET in mm
                        G_monthlyCellAET(n,month) += (dailyTotalAET_km3 * 1000000. /
                                                       ((cellArea * geo.G_contfreq[n] / 100.)));    // [mm]
                        G_monthlyCellAETWCa(n,month) += dailyTotalAET_km3;    // [km3]

                        // open water evaporation from surface storages (without rivers)
                        G_monthlyOpenWaterEvap(n,month) += (dailyOpenWaterEvap_km3 * 1000000. /
                                                             ((cellArea * geo.G_contfreq[n] / 100.)));    // [mm]

                    }
                        // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                    else {
                        G_monthlyCellAET(n,month) = 0.;
                        G_monthlyCellAETWCa(n,month) = 0.;
                        G_monthlyOpenWaterEvap(n,month) = 0.;
                    }
                }

                // daily output option 31
                if ((3 == options.grid_store) || (4 == options.grid_store)) {
                    if (options.outCellAETDaily) {
                        // Only continental cells
                        if (geo.G_contcell[n]) {
                            // AET in mm
                            G_daily31CellAET(n,day_in_month - 1) = (dailyTotalAET_km3 * 1000000. /
                                                                     ((cellArea * geo.G_contfreq[n] / 100.)));
                        }
                            // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                        else {
                            G_daily31CellAET(n,day_in_month - 1) = 0.;
                        }
                    }
                    if (options.outCellAETWCaDaily) {
                        if (geo.G_contcell[n]) {
                            // AETWCa in km3
                            G_daily31CellAETWCa(n,day_in_month - 1) += dailyTotalAET_km3;
                        }
                            // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                        else {
                            G_daily31CellAETWCa(n,day_in_month - 1) = 0.;
                        }
                    }
                }

                // daily output option 365
                if ((5 == options.grid_store) || (6 == options.grid_store)) {

                    if ((options.outCellAETDaily) || ((options.scoutCellAET) && (2 == options.day_store))) {
                        // Only continental cells
                        if (geo.G_contcell[n]) {
                            // AET in mm
                            G_daily365CellAET(n,day - 1) = (dailyTotalAET_km3 * 1000000. /
                                                             ((cellArea * geo.G_contfreq[n] / 100.)));
                        }
                            // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                        else {
                            G_daily365CellAET(n,day - 1) = 0.;
                        }

                    }
                    if (options.outCellAETWCaDaily) {

                        if (geo.G_contcell[n]) {
                            // AET in km3
                            G_daily365CellAETWCa(n,day - 1) += dailyTotalAET_km3;
                        }
                            // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                        else {
                            G_daily365CellAETWCa(n,day - 1) = 0.;
                        }
                    }
                }
            }
        }
    }


    //calculate total monthly storage per grid cell
    // here: from lateral water balance (second part)
    // in km3
    if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {
        for (int n = 0; n < ng; n++) {
            if (0 != G_toBeCalculated[n]) {

                //output for reservoir and storages start
                if (options.resOpt == 1) {
                    if (G_reservoir_area[n] > 0.) {
                        G_monthlyResInflow(n,month) += G_daily_res_in[n]; // total inflow in month [km3/day] -> [km3/month]
                        G_monthlyResOutflow(n,month) += G_daily_res_out[n]; // total outflow in month [km3/day] -> [km3/month]
                    } else {
                        G_monthlyResInflow(n,month) = 0.; // no reservoir area --> no Inflow (anyhow, output is not yet implemented)
                        G_monthlyResOutflow(n,month) = 0.; // no reservoir area --> no Outflow
                    }
                }

                // Only for runs which consider glaciers
                if (options.glacierOpt == 1) {
                    if (options.grid_store_TypeForStorages == 0)  // values at end of time step
                        G_monthlyGlacierStorage(n,month) = G_glacierStorage[n];
                    else  // mean values of time step
                        G_monthlyGlacierStorage(n,month) += G_glacierStorage[n];
                }

                // output if
                // G_monthlyResStorage, G_monthlyLocLakeStorage, G_monthlyGloLakeStorage, G_monthlyRiverStorage,
                // G_monthlySurfStor, G_monthlyLocWetlStorage, G_monthlyGloWetlStorage are values at end of month
                if (options.grid_store_TypeForStorages == 0) {
                    if (options.resOpt == 1) {

                        G_monthlySurfStor(n,month) =
                                       (G_locLakeStorage[n]
                                        + G_locWetlStorage[n]
                                        + G_gloLakeStorage[n]
                                        + G_gloWetlStorage[n]
                                        + G_riverStorage[n]
                                        + G_gloResStorage[n]);
                    }
                    else {
                        // total surface water storage in grid cell (lakes/wetlands/rivers) // [km3]
                        G_monthlySurfStor(n,month) =
                                       (G_locLakeStorage[n]
                                        + G_locWetlStorage[n]
                                        + G_gloLakeStorage[n]
                                        + G_gloWetlStorage[n]
                                        + G_riverStorage[n]);
                    }
                }

                // Probably double-counting here, as same procedure already present in loop over time steps
                // if mean monthly surface storage value is of interest
                if (options.grid_store_TypeForStorages == 1) {
                    if ((options.antNatOpt == 0) && (options.resOpt == 1)) {
                        G_monthlySurfStor(n,month) +=
                                       (G_locLakeStorage[n]
                                        + G_locWetlStorage[n]
                                        + G_gloLakeStorage[n]
                                        + G_gloWetlStorage[n]
                                        + G_riverStorage[n]
                                        + G_gloResStorage[n]);


                    } else {
                        // mean value of total surface water storage in grid cell (lakes/wetlands/rivers) [km3] (converted to mm for optional output in a later step)
                            G_monthlySurfStor(n,month) +=
                                           (G_locLakeStorage[n]
                                            + G_locWetlStorage[n]
                                            + G_gloLakeStorage[n]
                                            + G_gloWetlStorage[n]
                                            + G_riverStorage[n]);
                    }
                }
                //output for natural lakes and wetlands
                if (options.grid_store_TypeForStorages == 0) {
                    G_monthlyLocLakeStorage(n,month) = G_locLakeStorage[n];
                    G_monthlyGloLakeStorage(n,month) = G_gloLakeStorage[n];
                    G_monthlyRiverStorage(n,month) = G_riverStorage[n];
                    G_monthlyLocWetlStorage(n,month) = G_locWetlStorage[n];
                    G_monthlyGloWetlStorage(n,month) = G_gloWetlStorage[n];
                    G_monthlyGwStorage(n,month) = G_groundwaterStorage[n];

                    if ((options.antNatOpt == 0) && (options.resOpt == 1))
                        G_monthlyResStorage(n,month) = G_gloResStorage[n];

                } else if (options.grid_store_TypeForStorages == 1) {
                    G_monthlyLocLakeStorage(n,month) += G_locLakeStorage[n];
                    G_monthlyGloLakeStorage(n,month) += G_gloLakeStorage[n];
                    G_monthlyRiverStorage(n,month) += G_riverStorage[n];
                    G_monthlyLocWetlStorage(n,month) += G_locWetlStorage[n];
                    G_monthlyGloWetlStorage(n,month) += G_gloWetlStorage[n];
                    G_monthlyGwStorage(n,month) += G_groundwaterStorage[n];

                    if ((options.antNatOpt == 0) && (options.resOpt == 1))
                        G_monthlyResStorage(n,month) += G_gloResStorage[n];
                }

//FP20161213N001 Item 2 end
            }
        }
    }


    if ((3 == options.grid_store) || (4 == options.grid_store)) {
        // .31 gw storage output calculation in km3
        // calculation G_monthlyGwStorage formerly in daily.cpp, has to be calculated after NUg is subtracted
        for (int n = 0; n < ng; n++) {
            if (0 != G_toBeCalculated[n]) {
                if (options.outGWStorageDaily || options.outSingleStoragesDaily)
                    G_daily31GwStor(n,day_in_month - 1) += G_groundwaterStorage[n];
            }
        }
    }


    if ((5 == options.grid_store) || (6 == options.grid_store)) {
        // .365 gw storage output - calculation in km3
        for (int n = 0; n < ng; n++) {
            if (0 != G_toBeCalculated[n]) {
                if ((options.outGWStorageDaily || options.outSingleStoragesDaily) ||
                    ((options.scoutGwStor) && (2 == options.day_store)))
                    G_daily365GwStor(n,day - 1) = G_groundwaterStorage[n];
            }
        }
    }

    // daily 31 storage calculation
    if (((3 == options.grid_store) || (4 == options.grid_store)) &&
        ((options.outTotalWaterInStoragesDaily_km3) || (options.outTotalWaterInStoragesDaily_mm))) {
        for (int n = 0; n < ng; n++) {
            if (0 != G_toBeCalculated[n]) {

                if (options.outSingleStoragesDaily) {
                    G_daily31LocLakeStor(n,day_in_month - 1) = G_locLakeStorage[n];
                    G_daily31LocWetlStor(n,day_in_month - 1) = G_locWetlStorage[n];
                    G_daily31GloLakeStor(n,day_in_month - 1) = G_gloLakeStorage[n];
                    G_daily31GloWetlStor(n,day_in_month - 1) = G_gloWetlStorage[n];
                    G_daily31RiverStor(n,day_in_month - 1) = G_riverStorage[n];
                    if (options.resOpt == 1)
                        G_daily31ResStor(n,day_in_month - 1) = G_gloResStorage[n];
                }

                // daily total water storage output (calculate here after subtracting water use from groundwater)
                // calculate daily 31 storage per grid cell
                // here: get from vertical water balance (first part)
                // in km³
                double storage_daily31 = 0.;
                storage_daily31 = dailyWaterBalance.G_daily31CanopyStorage(n,day_in_month - 1) // [km3]
                                  + dailyWaterBalance.G_daily31SnowStorage(n,day_in_month - 1)
                                  + dailyWaterBalance.G_daily31SoilStorage(n,day_in_month - 1);

                // here: add the storages from routing (and later groundwater) (second part)
                G_daily31TotalWaterInStorages_km3(n,day_in_month - 1) =
                               storage_daily31
                               + G_locLakeStorage[n]
                               + G_locWetlStorage[n]
                               + G_gloLakeStorage[n]
                               + G_gloWetlStorage[n]
                               + G_riverStorage[n]
                               + G_groundwaterStorage[n];
                if (options.resOpt) {
                    G_daily31TotalWaterInStorages_km3(n,day_in_month - 1) += G_gloResStorage[n];
                }

                // Only continental cells
                if (geo.G_contcell[n]) {
                    // convert total storage from km3 to mm
                    G_daily31TotalWaterInStorages_mm(n,day_in_month - 1) =
                                   G_daily31TotalWaterInStorages_km3(n,day_in_month - 1)
                                   / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000.;
                }
                    // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                else {
                    G_daily31TotalWaterInStorages_mm(n,day_in_month - 1) = 0.;
                }

                // daily surface water storage km3
                if (options.outSurfStorDaily) {
                    if (options.resOpt == 1)
                        G_daily31SurfStor(n,day_in_month - 1) +=
                                       (G_locLakeStorage[n]
                                        + G_locWetlStorage[n]
                                        + G_gloLakeStorage[n]
                                        + G_gloWetlStorage[n]
                                        + G_riverStorage[n]
                                        + G_gloResStorage[n]);
                    else
                        G_daily31SurfStor(n,day_in_month - 1) += // [km3]
                                       (G_locLakeStorage[n]
                                        + G_locWetlStorage[n]
                                        + G_gloLakeStorage[n]
                                        + G_gloWetlStorage[n]
                                        + G_riverStorage[n]);
                }
            }
        }
    }


    if (((5 == options.grid_store) || (6 == options.grid_store)) &&
        ((options.outTotalWaterInStoragesDaily_km3) || (options.outTotalWaterInStoragesDaily_mm) ||
         (((options.scoutTWSkm3) || (options.scoutTWSmm)) && (2 == options.day_store)))) {
        for (int n = 0; n < ng; n++) {
            if (0 != G_toBeCalculated[n]) {

                if ((options.outSingleStoragesDaily) || ((options.scoutLocLake) && (2 == options.day_store)))
                    G_daily365LocLakeStor(n,day - 1) = G_locLakeStorage[n];
                if ((options.outSingleStoragesDaily) || ((options.scoutLocWet) && (2 == options.day_store)))
                    G_daily365LocWetlStor(n,day - 1) = G_locWetlStorage[n];
                if ((options.outSingleStoragesDaily) || ((options.scoutGloLake) && (2 == options.day_store)))
                    G_daily365GloLakeStor(n,day - 1) = G_gloLakeStorage[n];
                if ((options.outSingleStoragesDaily) || ((options.scoutGloWet) && (2 == options.day_store)))
                    G_daily365GloWetlStor(n,day - 1) = G_gloWetlStorage[n];
                if ((options.outSingleStoragesDaily) || ((options.scoutRiver) && (2 == options.day_store)))
                    G_daily365RiverStor(n,day - 1) = G_riverStorage[n];
                if ((options.outSingleStoragesDaily) || ((options.scoutReservoir) && (2 == options.day_store))) {
                    if (options.resOpt == 1)
                        G_daily365ResStor(n,day - 1) = G_gloResStorage[n];
                }
                if (((options.outSingleStoragesDaily) && (options.glacierOpt == 1)) || ((options.outGlacierStorageDaily) && (options.glacierOpt == 1)))
                    G_daily365GlacierStorage(n,day - 1) = G_glacierStorage[n];
                if ((options.outGlacierStorageDaily_mm) && (options.glacierOpt == 1))
                    G_daily365GlacierStorage_mm(n,day - 1) = G_glacierStorage[n] / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000.0;

                // daily total water storage output (calculate here after subtracting water use from groundwater)
                // calculate daily 365 storage per grid cell
                // here: get from vertical water balance (first part)
                // in km³

                double storage_daily365 = 0.;
                storage_daily365 = dailyWaterBalance.G_daily365CanopyStorage(n,day - 1) // [km3]
                                   + dailyWaterBalance.G_daily365SnowStorage(n,day - 1)
                                   + dailyWaterBalance.G_daily365SoilStorage(n,day - 1);

                // here: add the storages from routing (and later groundwater) (second part)
                G_daily365TotalWaterInStorages_km3(n,day - 1) =
                               storage_daily365
                               + G_locLakeStorage[n]
                               + G_locWetlStorage[n]
                               + G_gloLakeStorage[n]
                               + G_gloWetlStorage[n]
                               + G_riverStorage[n]
                               + G_groundwaterStorage[n];
                if (options.resOpt || (((options.scoutTWSkm3) || (options.scoutTWSmm)) && (2 == options.day_store)))
                    G_daily365TotalWaterInStorages_km3(n,day - 1) += G_gloResStorage[n];
                if (options.glacierOpt == 1)
                    G_daily365TotalWaterInStorages_km3(n,day - 1) += G_daily365GlacierStorage(n,day - 1);


                if (options.outTotalWaterInStoragesDaily_mm || options.scoutTWSmm) {
                    // Only continental cells
                    if (geo.G_contcell[n]) {
                        // convert total storages from km3 to mm
                        G_daily365TotalWaterInStorages_mm(n,day - 1) =
                                       G_daily365TotalWaterInStorages_km3(n,day - 1)
                                       / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000.;
                    }
                        // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                    else {
                        G_daily365TotalWaterInStorages_mm(n,day - 1) = 0.;
                    }
                }

                // daily surface water storage km3
                if ((options.outSurfStorDaily) || ((options.scoutSurfaceStor) && (2 == options.day_store))) {
                    if (options.resOpt == 1)
                        G_daily365SurfStor(n,day - 1) =
                                       (G_locLakeStorage[n]
                                        + G_locWetlStorage[n]
                                        + G_gloLakeStorage[n]
                                        + G_gloWetlStorage[n]
                                        + G_riverStorage[n]
                                        + G_gloResStorage[n]);

                    else
                        G_daily365SurfStor(n,day - 1) = // [km3]
                                       (G_locLakeStorage[n]
                                        + G_locWetlStorage[n]
                                        + G_gloLakeStorage[n]
                                        + G_gloWetlStorage[n]
                                        + G_riverStorage[n]);
                }
            }
        }
    }

    //HMS startend TWS output
    if (options.outTotalWaterInStoragesStartEndDaily_km3) {
        for (int n = 0; n < ng; n++) {
            if (0 != G_toBeCalculated[n]) {
                if (day == 1) { // calculate for first day in year
                    G_startendTotalWaterInStorages_km3(n,0) =
                                   dailyWaterBalance.G_canopyWaterContent[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 *
                                   G_landAreaFrac[n] / 100.0 // [km3]
                                   + dailyWaterBalance.G_snow[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 *
                                     G_landAreaFrac[n] / 100.0
                                   + dailyWaterBalance.G_soilWaterContent[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 *
                                     G_landAreaFrac[n] / 100.0
                                   + G_locLakeStorage[n]
                                   + G_locWetlStorage[n]
                                   + G_gloLakeStorage[n]
                                   + G_gloWetlStorage[n]
                                   + G_riverStorage[n]
                                   + G_groundwaterStorage[n];
                    if (options.resOpt)
                        G_startendTotalWaterInStorages_km3(n,0) += G_gloResStorage[n];
                    if (options.glacierOpt)
                        G_startendTotalWaterInStorages_km3(n,0) += G_glacierStorage[n];
                }
                if (day == 365) { // calculate for last day in year
                    G_startendTotalWaterInStorages_km3(n,1) =
                                   dailyWaterBalance.G_canopyWaterContent[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 *
                                   G_landAreaFrac[n] / 100.0 // [km3]
                                   + dailyWaterBalance.G_snow[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 *
                                     G_landAreaFrac[n] / 100.0
                                   + dailyWaterBalance.G_soilWaterContent[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 *
                                     G_landAreaFrac[n] / 100.0
                                   + G_locLakeStorage[n]
                                   + G_locWetlStorage[n]
                                   + G_gloLakeStorage[n]
                                   + G_gloWetlStorage[n]
                                   + G_riverStorage[n]
                                   + G_groundwaterStorage[n];
                    if (options.resOpt)
                        G_startendTotalWaterInStorages_km3(n,1) += G_gloResStorage[n];
                    if (options.glacierOpt)
                        G_startendTotalWaterInStorages_km3(n,1) += G_glacierStorage[n];
                }
            }
        }
    }


    if (options.day_store == 1) {
        for (short b = 0; b < nSpecBasins; b++) {
            dailyLocLakeStorage(b,day - 1) = getActualStorageRatio(cbasin.cellNum[b],
                                                                    1);  // ATTENTION: method using only cell-specific calibration parameter P_LAK_D  // FP
            dailyLocWetlStorage(b,day - 1) = getActualStorageRatio(cbasin.cellNum[b], 2);
            dailyGloLakeStorage(b,day - 1) = getActualStorageRatio(cbasin.cellNum[b],
                                                                    3);  // ATTENTION: method using only cell-specific calibration parameter P_LAK_D  // FP
            dailyGloWetlStorage(b,day - 1) = getActualStorageRatio(cbasin.cellNum[b], 4);
            if (options.resOpt == 1)
                dailyResStorage(b,day - 1) = getActualStorageRatio(cbasin.cellNum[b], 5);
        }
    }


    // AEMS: Store daily values in wghmStateFile:
    for (int n = 0; n < ng; n++) {
        //HMS I don't understand why .locallake(day_in_month-1) is written, and not (0) - the wghmState.cell-Array is in my eyes not multidimensional? Check with Maike / Annette needed.
        //
        cellArea = geo.areaOfCellByArrayPos(n);
        wghmState.cell(n).locallake(day_in_month - 1) =
                       G_locLakeStorage[n] / ((cellArea * (geo.G_contfreq[n] / 100.)) / 1000000.);
        wghmState.cell(n).localwetland(day_in_month - 1) =
                       G_locWetlStorage[n] / ((cellArea * (geo.G_contfreq[n] / 100.)) / 1000000.);
        wghmState.cell(n).globallake(day_in_month - 1) =
                       G_gloLakeStorage[n] / ((cellArea * (geo.G_contfreq[n] / 100.)) / 1000000.);
        wghmState.cell(n).globalwetland(day_in_month - 1) =
                       G_gloWetlStorage[n] / ((cellArea * (geo.G_contfreq[n] / 100.)) / 1000000.);
        wghmState.cell(n).reservoir(day_in_month - 1) =
                       G_gloResStorage[n] / ((cellArea * (geo.G_contfreq[n] / 100.)) / 1000000.);
        wghmState.cell(n).river(day_in_month - 1) =
                       G_riverStorage[n] / ((cellArea * (geo.G_contfreq[n] / 100.)) / 1000000.);
        wghmState.cell(n).groundwater(day_in_month - 1) =
                       G_groundwaterStorage[n] / ((cellArea * (geo.G_contfreq[n] / 100.)) / 1000000.);
    }
    // AEMS: store G_totalUnsatisfiedUse of last day of a month in additionalOutputInputFile
    if (day == last_day_in_month + 1) // for the last day of a month
    {
        for (int n = 0; n < ng; n++) {
            additionalOutIn.additionalOutputInput(n, 3) = G_totalUnsatisfiedUse[n];
            additionalOutIn.additionalOutputInput(n, 4) = G_UnsatisfiedUsePrevYear[n];
            additionalOutIn.additionalOutputInput(n, 50) = G_reducedReturnFlowPrevYear[n];
            additionalOutIn.additionalOutputInput(n, 51) = G_unsatisfiedNAsFromIrrigPrevYear[n];
            additionalOutIn.additionalOutputInput(n, 52) = G_unsatisfiedNAsFromOtherSectorsPrevYear[n];

        }
    }
    // Computation of SWB and land area fractions for the next time step
    for (int n = 0; n < ng; n++) {

        // a second cell is allowed to satisfy the remaining demand

        // multiply actual size of each swb with specific area in % and divide through continental area in % to get actual fraction of surface water bodies for the next time step
        // G_lake_area and G_reservoir_area are in km2 and have therefore considered differently to get the fraction of the continental area of those
        // to avoid numerical problems (sometimes the area in % gets nan), each swb type is summed up seperately.
        // for arid and for humid cells

        // Reservoir operation start years
        // Initial values are now assigned in routing.init
        if ((G_loc_lake[n] > 0.) && (G_locLakeAreaReductionFactor[n] > 0.)) {
            G_fLocLake[n] = (G_locLakeAreaReductionFactor[n] * G_loc_lake[n] / 100.);
        } else {
            G_fLocLake[n] = 0.; // to prevent that a single fraction values of a previous grid cell is used (if the conditions of the if clause is not true), the fractions will be 0 if the if clause is not true.
        }
        if (day == last_day_in_month + 1) {
            additionalOutIn.additionalOutputInput(n, 45) = G_fLocLake[n];
        }
        if ((G_loc_wetland[n] > 0.) && (G_locWetlAreaReductionFactor[n] > 0.)) {
            G_fLocWet[n] = (G_locWetlAreaReductionFactor[n] * G_loc_wetland[n] / 100.);
        } else {
            G_fLocWet[n] = 0.;
        }

        if (day == last_day_in_month + 1) {
            additionalOutIn.additionalOutputInput(n, 47) = G_fLocWet[n];
        }

        if ((G_glo_wetland[n] > 0.) && (G_gloWetlAreaReductionFactor[n] > 0.)) {
            G_fGloWet[n] = (G_gloWetlAreaReductionFactor[n] * G_glo_wetland[n] / 100.);
        } else {
            G_fGloWet[n] = 0.;
        }

        if (day == last_day_in_month + 1) {
            additionalOutIn.additionalOutputInput(n, 46) = G_fGloWet[n];
        }

        // Reservoir operation start years
        // calculate changes in G_fswb (initialization is already done in routing.init)

        // Old fraction from previously calculated time step
        G_fswbLandAreaFrac[n] = G_fswbLandAreaFracNextTimestep[n];  // In first time step identical values from routing.init
        // New current fraction
        G_fswbLandAreaFracNextTimestep[n] = G_fLocLake[n] + G_fLocWet[n] + G_fGloWet[n];
        G_fswbLandAreaFrac_ChangePct[n] = G_fswbLandAreaFracNextTimestep[n] * 100. - G_fswbLandAreaFrac[n] * 100.0;

        // adapt changes in landareafrac due to changes in river area fraction
        if (options.riverEvapoOpt == 1) {

            G_maxRiverAreaFrac[n] = geo.G_contfreq[n] / 100. -
                                    G_fGloLake[n]; //river area fraction can not be larger than "available" land area fraction

            // FP comment: inland sink?? and 100% glo lake
            if (((G_lake_area[n] > 0.) || (G_reservoir_area[n] > 0.)) && (G_fGloLake[n] == 1.)) {
                // for outlet cells containing 100% global lakes, no river evapo at all
                G_riverAreaFracNextTimestep_Frac[n] = 0.;
                G_riverAreaFrac_ChangeFrac[n] = 0.;
                G_riverAreaReductionFactor[n] = 0.;
                G_dailyRiverEvapo[n] = 0.;
            } else {
                // no inland sink and no 100% glo lake
                if (G_riverAreaFracNextTimestep_Frac[n] <= G_maxRiverAreaFrac[n]) {

                    //river fraction can reach max value
                    if (G_fswbLandAreaFracNextTimestep[n] >
                        (G_maxRiverAreaFrac[n] - G_riverAreaFracNextTimestep_Frac[n])) {

                        fswbFracCorr = (G_maxRiverAreaFrac[n] - G_riverAreaFracNextTimestep_Frac[n]) /
                                       G_fswbLandAreaFracNextTimestep[n];
                        G_locLakeAreaReductionFactor[n] *= fswbFracCorr;
                        G_locWetlAreaReductionFactor[n] *= fswbFracCorr;
                        G_gloWetlAreaReductionFactor[n] *= fswbFracCorr;

                        if ((G_fLocLake[n] > 0.) && (G_locLakeAreaReductionFactor[n] > 0.)) {
                            G_fLocLake[n] = (G_locLakeAreaReductionFactor[n] * G_loc_lake[n] /
                                             100.); // and finalize swb fraction
                        } else {
                            G_locLakeAreaReductionFactor[n] = 0.;
                            G_fLocLake[n] = 0.;
                        }
                        if ((G_fLocWet[n] > 0.) && (G_locWetlAreaReductionFactor[n] > 0.)) {
                            G_fLocWet[n] = (G_locWetlAreaReductionFactor[n] * G_loc_wetland[n] /
                                            100.); // and finalize swb fraction
                        } else {
                            G_locWetlAreaReductionFactor[n] = 0.;
                            G_fLocWet[n] = 0.;
                        }

                        if ((G_fGloWet[n] > 0.) && (G_gloWetlAreaReductionFactor[n] > 0.)) {
                            G_fGloWet[n] = (G_gloWetlAreaReductionFactor[n] * G_glo_wetland[n] /
                                            100.); // and finalize swb fraction
                        } else {
                            G_gloWetlAreaReductionFactor[n] = 0.;
                            G_fGloWet[n] = 0.;
                        }

                    }
                }
                else {
                    //river fraction cannot reach max value, thus has to be limited as well as river evaporation
                    riverAreaFracDeficit = G_riverAreaFracNextTimestep_Frac[n] -
                                           G_maxRiverAreaFrac[n]; //calculate missing fraction
                    G_riverAreaFrac_ChangeFrac[n] -= riverAreaFracDeficit; //reduce fraction change
                    G_riverAreaReductionFactor[n] *= G_maxRiverAreaFrac[n] /
                                                     G_riverAreaFracNextTimestep_Frac[n];//adapt (reduce) riverAreaReductionFactor
                    G_riverAreaFracNextTimestep_Frac[n] = G_maxRiverAreaFrac[n]; //set fraction to max fraction

                    G_fLocLake[n] = 0.;
                    G_fLocWet[n] = 0.;
                    G_fGloWet[n] = 0.;
                }

            }

            G_fswbLandAreaFracNextTimestep[n] = G_fLocLake[n] + G_fLocWet[n] + G_fGloWet[n];
            G_fswbLandAreaFrac_ChangePct[n] = G_fswbLandAreaFracNextTimestep[n] * 100. - G_fswbLandAreaFrac[n] *
                                                                                         100.; // G_fwsb is in fraction, G_landAreaFrac (where G_fswbChange is calculated) is in %

        }
        else {
            //no evaporation from rivers
            G_riverAreaFracNextTimestep_Frac[n] = 0.;
            G_riverAreaFrac_ChangeFrac[n] = 0.;
        }

        //calculate final G_fswb (should be now max. continental frac)
        G_fswb[n] = G_fLocLake[n] + G_fLocWet[n] + G_fGloLake[n] + G_fGloWet[n] + G_glo_res[n] +
                    G_riverAreaFracNextTimestep_Frac[n]; // for output purposes only!

        if (day == last_day_in_month + 1) {
            if (month==11)  // month == 11 is december see routing call in integrateWGHM
                additionalOutIn.additionalOutputInput(n, 44) = G_glo_res[n];    // next month is next year so prevYear update
            else
                additionalOutIn.additionalOutputInput(n, 44) = G_glores_prevyear[n];    // other months prevYear stays the same
        }

        //now, adapt land area fraction for changes in fswb and riverareafrac
        statusStarted_landAreaFracNextTimestep[n] = 1;
        G_landAreaFracNextTimestep[n] = G_landAreaFrac[n] -
                                            (G_fswbLandAreaFrac_ChangePct[n] + (G_riverAreaFrac_ChangeFrac[n] * 100.));
        // Control not allowed range of values: Set to zero if calculated land fraction is less than a limit
        // limit currently set to zero
        // e.g. also for Caspian Sea cells
        if (G_landAreaFracNextTimestep[n] < 0.) {
            G_landAreaFracNextTimestep[n] = 0.;
        }

        // Sum values for monthly averaging, only with specific writing options
        if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {
            G_monthlyLandAreaFrac(n,month) += G_landAreaFrac[n];
        }

    }
    // end of computation of SWB and land area fractions for the next day

    // now write out final G_fswb
    for (int n = 0; n < ng; n++) {
        if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {
            if (options.outFswb == 1)
                G_monthlyFswb(n,month) += G_fswb[n];
        }
        if ((3 == options.grid_store) || (4 == options.grid_store)) {
            if (options.outFswbDaily)
                G_daily31Fswb(n,day_in_month - 1) = G_fswb[n]; // daily fraction of surface water bodies [mm]
        }
        if ((5 == options.grid_store) || (6 == options.grid_store)) {
            if ((options.outFswbDaily) || ((options.scoutFswb) && (2 == options.day_store)))
                G_daily365Fswb(n,day - 1) = G_fswb[n];
        }
    }
    //AEMS:
    if (day == last_day_in_month + 1){
         for (int n = 0; n < ng; n++){
            additionalOutIn.additionalOutputInput(n, 8) = G_locWetlAreaReductionFactor[n];
            additionalOutIn.additionalOutputInput(n, 9) = G_gloLakeEvapoReductionFactor[n];
            additionalOutIn.additionalOutputInput(n, 10) = G_groundwaterStorage[n];
            additionalOutIn.additionalOutputInput(n, 11) = G_locLakeAreaReductionFactor[n];
            additionalOutIn.additionalOutputInput(n, 12) = G_gloWetlAreaReductionFactor[n];
            additionalOutIn.additionalOutputInput(n, 13) = G_gloResEvapoReductionFactor[n];
            additionalOutIn.additionalOutputInput(n, 14) = G_fswbInit[n];
            additionalOutIn.additionalOutputInput(n, 15) = G_gloWetlStorage[n];
//            additionalOutIn.additionalOutputInput(n, 16) = G_fswbLandAreaFrac[n];
            additionalOutIn.additionalOutputInput(n, 17) = G_locWetlStorage[n];
            additionalOutIn.additionalOutputInput(n, 18) = G_locLakeStorage[n];
            additionalOutIn.additionalOutputInput(n, 19) = G_riverStorage[n];
            additionalOutIn.additionalOutputInput(n, 21) = G_gloLakeStorage[n];
            additionalOutIn.additionalOutputInput(n, 23) = G_gloResStorage[n];
            additionalOutIn.additionalOutputInput(n, 28) = G_AllocatedUse[n];
            additionalOutIn.additionalOutputInput(n, 29) = G_SecondCell[n];
            additionalOutIn.additionalOutputInput(n, 32) = G_AllocUseToNeigborCell[n];
            additionalOutIn.additionalOutputInput(n, 27) = G_UnsatAllocUse[n];
            additionalOutIn.additionalOutputInput(n, 30) = G_daily_UnsatAllocUseNextDay[n];
            additionalOutIn.additionalOutputInput(n, 31) = G_daily_allocatedUseNextDay[n];
            additionalOutIn.additionalOutputInput(n, 33) = G_fswbLandAreaFracNextTimestep[n];
            additionalOutIn.additionalOutputInput(n, 34) = G_PrevUnsatAllocUse[n];
            additionalOutIn.additionalOutputInput(n, 35) = G_dailyRemainingUse[n];
            additionalOutIn.additionalOutputInput(n, 36) = G_fGloLake[n];
            additionalOutIn.additionalOutputInput(n, 37) = G_PrevTotalUnstatisfiedUse[n];
            additionalOutIn.additionalOutputInput(n, 40) = G_dailySatisAllocatedUseInSecondCell[n];
            additionalOutIn.additionalOutputInput(n, 41) = G_dailyAllocatedUse[n];
            additionalOutIn.additionalOutputInput(n, 38) = G_reducedReturnFlow[n];
            additionalOutIn.additionalOutputInput(n, 39) = G_unsatisfiedNAsFromIrrig[n];
            additionalOutIn.additionalOutputInput(n, 49) = G_unsatisfiedNAsFromOtherSectors[n]; // Grid with unsatisfied NAs from other sectors than irrig associated with G_reducedReturnFlow
            additionalOutIn.additionalOutputInput(n, 25) = G_withdrawalIrrigFromSwb[n] ;
            additionalOutIn.additionalOutputInput(n, 26) = G_consumptiveUseIrrigFromSwb[n];

         }
    }
}

void routingClass::annualWaterUsePostProcessing(short year, AdditionalOutputInputFile &additionalOutIn) {
    if (options.subtract_use >= 1) {
        char filename[250];
        int n;

        // Has there been a satisfaction of the water use of the previous year?
        // We want to store a grid with the amount of use which is still
        // unsatisfied from the previous year for control purposes.
        // After this year the unsatisfied use of the previous year will
        // not be considered any more. Therefore it is removed from
        // the grid 'G_totalUnsatisfiedUse'.
        for (n = 0; n <= ng - 1; n++) {
            if (G_totalUnsatisfiedUse[n] < G_UnsatisfiedUsePrevYear[n]) {
                // This means: all the desired uses of the current year have been
                // satisfied and a fraction of the previous year also.
                G_UnsatisfiedUsePrevYear[n] = G_totalUnsatisfiedUse[n];
                G_unsatisfiedNAsFromIrrigPrevYear[n] = G_unsatisfiedNAsFromIrrig[n];
                G_unsatisfiedNAsFromOtherSectorsPrevYear[n] = G_unsatisfiedNAsFromOtherSectors[n];
                G_reducedReturnFlowPrevYear[n] = G_reducedReturnFlow[n];
                G_unsatisfiedNAsFromIrrig[n] = 0.;
                G_unsatisfiedNAsFromOtherSectors[n] = 0.;
                G_reducedReturnFlow[n] = 0.;
                G_totalUnsatisfiedUse[n] = 0.;
            } else {
                // This means: we have even more unsatisfied use than at the end
                // of the previous year.
                // The grid of the year before has not to be changed,
                // but the grid of the current year should only contain the values
                // of the current year.
                G_totalUnsatisfiedUse[n] -= G_UnsatisfiedUsePrevYear[n];
                G_unsatisfiedNAsFromIrrig[n] -= G_unsatisfiedNAsFromIrrigPrevYear[n];
                G_unsatisfiedNAsFromOtherSectors[n] -= G_unsatisfiedNAsFromOtherSectorsPrevYear[n];
                G_reducedReturnFlow[n] -= G_reducedReturnFlowPrevYear[n];
                if ((G_unsatisfiedNAsFromIrrig[n] + G_unsatisfiedNAsFromOtherSectors[n]) < 0){
                    G_unsatisfiedNAsFromIrrig[n] = 0;
                    G_unsatisfiedNAsFromOtherSectors[n] = 0;
                    G_reducedReturnFlow[n] = 0;
                } else if (G_unsatisfiedNAsFromOtherSectors[n] < 0)
                {
                    G_unsatisfiedNAsFromIrrig[n] += G_unsatisfiedNAsFromOtherSectors[n];
                    G_unsatisfiedNAsFromOtherSectors[n] = 0;
                } else if (G_unsatisfiedNAsFromIrrig[n] < 0) {
                    G_unsatisfiedNAsFromOtherSectors[n] += G_unsatisfiedNAsFromIrrig[n];
                    G_unsatisfiedNAsFromIrrig[n] = 0;
                    G_reducedReturnFlow[n]=0;
                }
                if (G_reducedReturnFlow[n] > 0) {
                    G_reducedReturnFlow[n] = 0;
                    G_unsatisfiedNAsFromIrrig[n] = 0;
                }
            }
            // AEMS: store G_totalUnsatisfiedUse and G_UnsatisfiedUsePrevYear of last day of the year
            additionalOutIn.additionalOutputInput(n, 3) = G_totalUnsatisfiedUse[n];
            additionalOutIn.additionalOutputInput(n, 4) = G_UnsatisfiedUsePrevYear[n];
            additionalOutIn.additionalOutputInput(n, 50) = G_reducedReturnFlowPrevYear[n];
            additionalOutIn.additionalOutputInput(n, 51) = G_unsatisfiedNAsFromIrrigPrevYear[n];
            additionalOutIn.additionalOutputInput(n, 52) = G_unsatisfiedNAsFromOtherSectorsPrevYear[n];
            additionalOutIn.additionalOutputInput(n, 38) = G_reducedReturnFlow[n];
            additionalOutIn.additionalOutputInput(n, 39) = G_unsatisfiedNAsFromIrrig[n];
            additionalOutIn.additionalOutputInput(n, 49) = G_unsatisfiedNAsFromOtherSectors[n]; // Grid with unsatisfied NAs from other sectors than irrig associated with G_reducedReturnFlow
        }
        if ((options.grid_store > 0) && (year >= options.evalStartYear)) {
            Grid<> G_array;

            if (options.outUnsatUseSWPrev) {
                G_UnsatisfiedUsePrevYear.write(options.output_dir + "/G_UNSAT_USE_SW_PREV_" + to_string(year) + ".UNF0");
            }
            if (options.outUnsatUseSW) {
                G_totalUnsatisfiedUse.write(options.output_dir + "/G_UNSAT_USE_SW_" + to_string(year) + ".UNF0");
            }
            if (options.outActualNAS) { // new output options 2.1f
                for (n = 0; n < ng; n++) {
                    G_actualUse[n] = G_actualUse[n] *
                                     1000000000.; // conversion km3/month --> m3/month to avoid loss of small numbers during casting operation.
                    G_array[n] = G_actualUse[n]; // conversion double -> float
                }
                G_array.write(options.output_dir + "/G_ACTUAL_NAS_" + to_string(year) + ".UNF0");
            }

        }
        // Note: ACTUAL_USE_[YEAR] + UNSAT_USE_[YEAR] does not equal CONS_WATER_USE_[YEAR]
        // as unsatisfied use from the previous year (YEAR-1) might have partially
        // been satisfied in YEAR
        // ACTUAL_USE_[YEAR] + UNSAT_USE_[YEAR] = CONS_WATER_USE[YEAR] + (UNSAT_USE_[YEAR-1] - UNSAT_USE_PREV_[YEAR])

    }
}

void routingClass::checkTimeStep() {
    // this is the public function
    // should be called after G_toBeCaluclated has been generated
    // and after initialization of this class
    defaultTimeStepsPerDay = checkTimeStep(defaultRiverVelocity, defaultTimeStepsPerDay);
    // velocity
    riverVelocity = defaultRiverVelocity; // [km/day]
}

void routingClass::updateLandAreaFrac(AdditionalOutputInputFile &additionalOutIn) {
    //LandAreaFracNextTimestep is already calculated, now G_LandAreaFrac has to be updated to be used consistently in daily and routing for the next timestep
    for (int n = 0; n < ng; n++) {
        G_landAreaFracPrevTimestep[n] = G_landAreaFrac[n];
        G_landAreaFrac[n] = G_landAreaFracNextTimestep[n];
        //AEMS:
        additionalOutIn.additionalOutputInput(n,6) = G_landAreaFrac[n]; // G_landAreaFrac[n]
        additionalOutIn.additionalOutputInput(n,7) = G_landAreaFracPrevTimestep[n]; // G_landAreaFracPrevTimestep[n]
    }
}

void routingClass::update_landarea_red_fac_PDAF(calibParamClass &calParam, AdditionalOutputInputFile &additionalOutIn) {
    double cellArea = 0., maxStorage = 0.;
    evapoReductionExp = 3.32193;

    // routingCell: (WG2.2) added for ordered routing scheme
    for (int routingCell = 0; routingCell < ng; routingCell++) {

        int n = (G_routingCell[routingCell] - 1);
        double M_EVAREDEX__at__n_routing = calParam.getValue(M_EVAREDEX, n);

        cellArea = geo.areaOfCellByArrayPos(n);
        // local Lake
        if (G_loc_lake[n] > 0.) {
            maxStorage = ((G_loc_lake[n]) / 100.) * cellArea * G_lakeDepthActive[n]; // [km3]
            G_locLakeAreaReductionFactor[n] = 1. -
                                              pow(fabs(G_locLakeStorage[n] - maxStorage)
                                                  / (2. * maxStorage),
                                                  (M_EVAREDEX__at__n_routing * evapoReductionExp));

            if (G_locLakeAreaReductionFactor[n] < 0.)
                G_locLakeAreaReductionFactor[n] = 0.;

            if (G_locLakeAreaReductionFactor[n] > 1.)  // (could occur due to numerical inaccuracies)
                G_locLakeAreaReductionFactor[n] = 1.;
        }

        // local wetland
        if (G_loc_wetland[n] > 0.) {
            maxStorage = ((G_loc_wetland[n]) / 100.) * cellArea * G_wetlDepthActive[n];

            // Area reduction factor is now calculated for the next day.
            // Cell-specific calibration parameters - Apply multiplier
            G_locWetlAreaReductionFactor[n] = 1. -
                                              pow(fabs(G_locWetlStorage[n] - maxStorage)
                                                  / (maxStorage),
                                                  (M_EVAREDEX__at__n_routing * evapoReductionExp));


            if (G_locWetlAreaReductionFactor[n] < 0.)    // (could occur due to numerical inaccuracies)
                G_locWetlAreaReductionFactor[n] = 0.;

            if (G_locWetlAreaReductionFactor[n] > 1.)  // (could occur due to numerical inaccuracies)
                G_locWetlAreaReductionFactor[n] = 1.;
        }
        // global wetlands
        if (G_glo_wetland[n] > 0.) {
            maxStorage = ((G_glo_wetland[n]) / 100.) * cellArea * G_wetlDepthActive[n];

            G_gloWetlAreaReductionFactor[n] = 1. -
                                              pow(fabs(G_gloWetlStorage[n] - maxStorage)
                                                  / maxStorage,
                                                  (M_EVAREDEX__at__n_routing * evapoReductionExp));

            if (G_gloWetlAreaReductionFactor[n] < 0.)  // (could occur due to numerical inaccuracies)
                G_gloWetlAreaReductionFactor[n] = 0.;

            if (G_gloWetlAreaReductionFactor[n] > 1.)  // (could occur due to numerical inaccuracies)
                G_gloWetlAreaReductionFactor[n] = 1.;
        }
        // global lakes
        if (G_lake_area[n] > 0.) {
            maxStorage = ((double) G_lake_area[n]) * G_lakeDepthActive[n];

            G_gloLakeEvapoReductionFactor[n] = 1. -
                                               pow(fabs(G_gloLakeStorage[n] - maxStorage)
                                                   / (2. * maxStorage),
                                                   (M_EVAREDEX__at__n_routing * evapoReductionExp));


            if (G_gloLakeEvapoReductionFactor[n] < 0.)    // (could occur due to numerical inaccuracies)
                G_gloLakeEvapoReductionFactor[n] = 0.;

            if (G_gloLakeEvapoReductionFactor[n] > 1.)  // (could occur due to numerical inaccuracies)
                G_gloLakeEvapoReductionFactor[n] = 1.;
        }
        // reservoirs
        if (G_reservoir_area[n] > 0.) { //cell is a reservoir cell with a well known type
            maxStorage = G_stor_cap[n];

            //Evaporation reduction factor is now calculated for the next day.
            G_gloResEvapoReductionFactor[n] = 1. - pow(fabs(G_gloResStorage[n] - maxStorage)
                                                       / maxStorage, evapoReductionExpReservoir);

            if (G_gloResEvapoReductionFactor[n] < 0.)    // (could occur due to numerical inaccuracies)
                G_gloResEvapoReductionFactor[n] = 0.;

            if (G_gloResEvapoReductionFactor[n] > 1.)  // (could occur due to numerical inaccuracies)
                G_gloResEvapoReductionFactor[n] = 1.;
        }

        if ((G_loc_lake[n] > 0.) && (G_locLakeAreaReductionFactor[n] > 0.)) {
            G_fLocLake[n] = (G_locLakeAreaReductionFactor[n] * G_loc_lake[n] / 100.);
        } else {
            G_fLocLake[n] = 0.; // to prevent that a single fraction values of a previous grid cell is used (if the conditions of the if clause is not true), the fractions will be 0 if the if clause is not true.
        }

        if ((G_loc_wetland[n] > 0.) && (G_locWetlAreaReductionFactor[n] > 0.)) {
            G_fLocWet[n] = (G_locWetlAreaReductionFactor[n] * G_loc_wetland[n] / 100.);
        } else {
            G_fLocWet[n] = 0.;
        }

        if ((G_glo_wetland[n] > 0.) && (G_gloWetlAreaReductionFactor[n] > 0.)) {
            G_fGloWet[n] = (G_gloWetlAreaReductionFactor[n] * G_glo_wetland[n] / 100.);
        } else {
            G_fGloWet[n] = 0.;
        }
        // New current fraction

        G_fswbLandAreaFracNextTimestep[n] = G_fLocLake[n] + G_fLocWet[n] + G_fGloWet[n]; // nextTimestep bc. in std code in function routing nextTimestep is copied into fswbLandAreaFrac
        G_fswbLandAreaFrac[n]=additionalOutIn.additionalOutputInput(n, 33); // G_fswbLandAreaFracNextTimestep[n] originally calculated a day before, put into G_fswbLandAreaFrac[n] for convenience, will be overwritten anyway
        G_fswbLandAreaFrac_ChangePct[n]=G_fswbLandAreaFracNextTimestep[n] * 100. - G_fswbLandAreaFrac[n] * 100.;

        G_landAreaFrac[n] = additionalOutIn.additionalOutputInput(n,6);
        G_landAreaFrac[n] =  G_landAreaFrac[n] - (G_fswbLandAreaFrac_ChangePct[n]);

        if (G_landAreaFrac[n] < 0.)
            G_landAreaFrac[n] = 0.;
        G_landAreaFracPrevTimestep[n] = additionalOutIn.additionalOutputInput(n,7); //G_landAreaFracPrevTimestep[n];
        statusStarted_landfreq_fwaterfreq_landAreaFrac[n] = 1;
    }
}

// Reservoir operation start years
// Transfer land area fraction of reservoirs (in percent) of current year in array of previous year for use in the next year
void routingClass::updateGloResPrevYear_pct() {
    for (int n = 0; n < ng; n++) {
        G_glores_prevyear[n] = G_glo_res[n];
    }
}

/**
 *
 * @brief This method updates the groundwater net abstraction due to not statisfied SW usage
 *
 * This method corrects the net abstraction from groundwater in case of not fullfilled net abstraction from surface
 * water bodies(swb). GWSWUSE is assuming that potential water abstractions from swb are fullfilled and accordingly
 * calculates return flows to groundwater resulting from irrrigation with surface water, which impact net abstractions
 * from groundwater NAg. Thus return flows to groundwater (which can only result from irrigation in WGHM) can be
 * overestimated.  For this adaption it is assumed that irrigation is satisfied with the lowest priority,
 * i.e. the unsatisfied NAs is assumed to be caused by less NAs for irrigation, up to setting irrigation water use
 * to zero. Based on this assumption, a new value for water abstractions from surface water for irrigation (wusi_new)
 * is computed and NAg is then reduced to NAgnew following the equations provided in Döll et al. (2012).
 * Furthermore due to delayed use satisfaction and use satisfaction through neighbouring cells a helper variable is
 * used (G_dailyRemainingUse[n]) which is calculated as RemainingUse - RemainingUse from last time step. By this return
 * flows can be readded to NAg, when potential NAs is satisfied delayed.
 * @param n
 * @param month
 */
double routingClass::updateNetAbstractionGW(int n) {
    double NAgnew = 0.; // modified NAg (return flows reduced or reintroduced)
    double returnflowChange = 0.; // amount by which original NAg is changed could be negative (return flows are reduced) or positive (return flows are reintroduced)
    double eff = 0.; // efficency of irrig water usage from surface water bodies
    double frgi = 0.; // fraction of return flows which flow into groundwater
    double factor = 0.; // factor to calculate new waterwithdrawal for irrigation from swb
    double WUsiNew = 0.; // waterwithdrawal for irrig from swb calculated based on unsatisfied NAs
    double dailyRemainingUseFromIrrig = 0.; // dailyremaininguse change which originates from irrigation
    double reintroducedRatio = 0.; // negative ratio of reduced return flow which should be reintroduced due to delayed satisfaction

    // G_dailyRemainingUse is > 0 when returnflows need to be reduced, < 0 when returnflows need to be reintroduced due
    // to delayed satisfaction, due to numerical issues the threshold is set to 1.e-12 resp -1.e-12
    if (G_dailyRemainingUse[n] > 1.e-12) {
        // Only update NAg if WUsi is existent
        if (G_withdrawalIrrigFromSwb[n] > 0.) {
            eff = G_consumptiveUseIrrigFromSwb[n] / G_withdrawalIrrigFromSwb[n];
            frgi = G_fractreturngw_irrig[n];
            factor = 1 - (1 - frgi) * (1 - eff);
            WUsiNew = 1 / factor * (G_withdrawalIrrigFromSwb[n] * factor - G_dailyRemainingUse[n]);
            if (WUsiNew < 0.) {
                // WUsiNew is below zero, when dailyRemainingUse is bigger than the estimated partion
                // of NAs from Irrigation
                WUsiNew = 0.;
                // unsatisfiedNAsfromirrig is needed to scale the reintroduction of returnflows in NAg
                G_unsatisfiedNAsFromIrrig[n] += G_withdrawalIrrigFromSwb[n] * factor;
                G_unsatisfiedNAsFromOtherSectors[n] += G_dailyRemainingUse[n] - (G_withdrawalIrrigFromSwb[n] * factor);
            } else { // unsatisfied use stems solely from irrigation
                G_unsatisfiedNAsFromIrrig[n] += G_dailyRemainingUse[n];
            };
            returnflowChange = (frgi * (1 - eff) * (WUsiNew - G_withdrawalIrrigFromSwb[n]));
            G_reducedReturnFlow[n] += returnflowChange;
            NAgnew = G_dailydailyNUg[n] - returnflowChange;
            G_dailyRemainingUse[n] = 0.;
            return NAgnew;
        } else {
            G_unsatisfiedNAsFromOtherSectors[n] += G_dailyRemainingUse[n];
            return G_dailydailyNUg[n];
        }
    } else if (G_dailyRemainingUse[n] < -1.e-12 ) {
        dailyRemainingUseFromIrrig = G_dailyRemainingUse[n] + G_unsatisfiedNAsFromOtherSectors[n];
        if (dailyRemainingUseFromIrrig < 0.) {
            G_unsatisfiedNAsFromOtherSectors[n] = 0.;
        } else {
            G_unsatisfiedNAsFromOtherSectors[n] += G_dailyRemainingUse[n];
            G_dailyRemainingUse[n] = 0.;
            return G_dailydailyNUg[n];
        }
        if (G_unsatisfiedNAsFromIrrig[n] == 0.){
            // this if clause is catching cases where reduced return flows are tried to be reintroduced, which aren't
            // existing. This can happen in rare cases of numerical inaccuracies.
            G_dailyRemainingUse[n] = 0.;
            return G_dailydailyNUg[n];
        }
        reintroducedRatio = (dailyRemainingUseFromIrrig / G_unsatisfiedNAsFromIrrig[n]);
        if (reintroducedRatio < -1.){
            reintroducedRatio = -1.;
            dailyRemainingUseFromIrrig = G_unsatisfiedNAsFromIrrig[n] * -1.;
        }
        returnflowChange =  (reintroducedRatio * G_reducedReturnFlow[n]);

        G_unsatisfiedNAsFromIrrig[n] += dailyRemainingUseFromIrrig;
        G_reducedReturnFlow[n] += returnflowChange;
        NAgnew = G_dailydailyNUg[n] - returnflowChange;
        G_dailyRemainingUse[n] = 0.;
        return NAgnew;
    }
    else {
        return G_dailydailyNUg[n];
    }
}

short routingClass::checkTimeStep(double riverVelocity, short numberOfTimeSteps) {


    // find most nothern and southern row
    short minRow = 1000;
    short maxRow = 0;

    for (int n = 0; n < ng; n++) {
        if (G_toBeCalculated[n] != 0) {
            if (geo.G_row[n] < minRow)
                minRow = geo.G_row[n];
            if (geo.G_row[n] > maxRow)
                maxRow = geo.G_row[n];
        }
    }

    // transfer to northern hemisphere
    if (minRow > 180)
        minRow = 360 - minRow + 1;
    if (maxRow > 180)
        maxRow = 360 - maxRow + 1;

    // find most northern row
    if (maxRow < minRow)
        minRow = maxRow;

    // calculate horizontal extend of cell
    const double pi = 3.141592653589793;
    const double earthRadius = 6371.211; // [km]
    double lat, horizDist;

    lat = 0.25 + minRow * 0.5;
    horizDist = 2 * pi * earthRadius * sin(lat * pi / 180.0) * 0.5 / 360.0;

    short newNumberOfTimeSteps;

    newNumberOfTimeSteps = numberOfTimeSteps;
    if (riverVelocity / (horizDist * numberOfTimeSteps) > 0.75) {
        cout << "Time step to great: " << numberOfTimeSteps << " steps per day\n";
        newNumberOfTimeSteps = (short) ceil(riverVelocity / (0.6 * horizDist));
        cout << "Adjusted to: " << newNumberOfTimeSteps << endl;
    }
    return newNumberOfTimeSteps;
}

// Set (initialize) cell-specific active lake depth (in km) as array to enable individual time-series
// Start with calibration parameter
void routingClass::initLakeDepthActive(calibParamClass &calParam)// OE: added calParam
            {

    // Initialize array with calibration parameter values
    // Data unit of depth has to be converted (m to km)
    for (int n = 0; n < ng; n++) {
        G_lakeDepthActive[n] = (calParam.getValue(P_LAK_D, n) *
                                CONV_FACTOR__LENGTH__M_TO_KM); // [km]  // OBSOLETE now read from calibration parameter file // FP
    }

}

// Set (initialize) cell-specific active wetland depth (in km) as array to enable individual time-series
// Start with calibration parameter
void routingClass::initWetlDepthActive(calibParamClass &calParam) // OE: added calParam
            {

    // Initialize array with calibration parameter values
    // Data unit of depth has to be converted (m to km)
    for (int n = 0; n < ng; n++) {
        G_wetlDepthActive[n] = (calParam.getValue(P_WET_D, n) *
                                CONV_FACTOR__LENGTH__M_TO_KM); // [km]  // OBSOLETE now read from calibration parameter file // FP
    }
}


void routingClass::setLakeWetlToMaximum(const short start_year) {
    // execute at start of each calibration year
    // or at first model year // FP20161018N002: any possible specified year
    for (int n = 0; n <= ng - 1; n++) {

        // Cell-specific active lake depth - Use array (initialized from cell-specific parameter)
        // Cell-specific active wetland depth - Use array (initialized from cell-specific parameter)
        G_locLakeStorage[n] = (G_loc_lake[n] / 100.)
                              * geo.areaOfCellByArrayPos(n) * G_lakeDepthActive[n];
        G_locWetlStorage[n] = (G_loc_wetland[n] / 100.)
                              * geo.areaOfCellByArrayPos(n) * G_wetlDepthActive[n];

        G_gloLakeStorage[n] = G_lake_area[n] *
                              G_lakeDepthActive[n]; // added for new reservoir algorithm  (remark: int G_lake_area, in full km2)

        G_gloWetlStorage[n] = (G_glo_wetland[n] / 100.)
                              * geo.areaOfCellByArrayPos(n) * G_wetlDepthActive[n];

        // Use reservoir algorithm with reservoirs and regulated lakes
        if ((options.antNatOpt == 0) && (options.resOpt == 1)) {

            //set to storage capacity if no operational year and start year not later than reference year (e.g. 2000)
            // for both, reg. lakes and reservoirs
            if ((options.resYearOpt == 0) && (options.resYearReference >= G_res_start_year[n]) &&
                (G_stor_cap_full[n] > -99)) {
                G_gloResEvapoReductionFactor[n] = 1.;
                G_gloResStorage[n] = G_stor_cap_full[n]; //for new reservoir algorithm (vol. from published data)
            }
            else
                G_gloResEvapoReductionFactor[n] = 0.;
            //set reservoir storage to capacity for reservoirs and regulated lakes which are already operating at model start year
            if ((options.resYearOpt == 1) && (start_year >= G_res_start_year[n]) && (G_stor_cap_full[n] > -99)) {
                G_gloResEvapoReductionFactor[n] = 1.;
                G_gloResStorage[n] = G_stor_cap_full[n]; //for new reservoir algorithm (vol. from published data)
            } else
                G_gloResEvapoReductionFactor[n] = 0.;

            // HERE CHECK IF CORRECT / NO DOUBLE COUNTING
            // In case of regulated lake: storage = product of area * lake_depth
            //increase reservoir storage if regulated lake is not yet operating but handle it like if it would be a global lake
            if (((options.resYearOpt == 1) && (G_reg_lake_status[n] == 1) && (start_year < G_res_start_year[n]))
                || ((options.resYearOpt == 0) && (G_reg_lake_status[n] == 1) &&
                    (options.resYearReference < G_res_start_year[n]))) {
                G_gloResStorage[n] += G_reservoir_area_full[n] * G_lakeDepthActive[n];
                // G_gloResEvapoReductionFactor[n] = 0. due to if else statements above
            }
        }

        //initialize area reduction factors
        // initial value at start of first model year or at each calibration year)
        G_locLakeAreaReductionFactor[n] = 1.;
        G_locWetlAreaReductionFactor[n] = 1.;
        G_gloLakeEvapoReductionFactor[n] = 1.;
        G_gloWetlAreaReductionFactor[n] = 1.;

        if (options.riverEvapoOpt == 1) {
            G_riverAreaReductionFactor[n] = 0.5; //corresponds to around 20% of bankfull storage for the first model day
            // river area fraction (as fraction of cell area)
            // with river width of bankfull flow [unit: m]
            // ATTENTION:
            // Area adds to water area!
            // Make sure that total area does not exceed cell area!
            // E.g. subtract from surface water bodies (e.g. lakes, wetlands)
            G_riverAreaFracNextTimestep_Frac[n] =
                           G_riverAreaReductionFactor[n] * G_riverLength[n] * G_RiverWidth_bf[n] / 1000. /
                           geo.areaOfCellByArrayPos(n);
        } else {
            //no evaporation from rivers
            G_riverAreaReductionFactor[n] = 0.;
            G_riverAreaFracNextTimestep_Frac[n] = 0.;
            G_riverAreaFrac_ChangeFrac[n] = 0.;
        }
    }
}

double routingClass::getActualStorage(int cellNumber, short selection) {
    switch (selection) {
        case '1':
            return G_locLakeStorage[cellNumber - 1];
        case '2':
            return G_locWetlStorage[cellNumber - 1];
        case '3':
            return G_gloLakeStorage[cellNumber - 1];
        case '4':
            return G_gloWetlStorage[cellNumber - 1];
        case '5':
            return G_gloResStorage[cellNumber - 1];
        default:
            cerr << "Your selection is not available\n";
            return -9999;
    }
}

double routingClass::getActualStorageRatio(int cellNumber, short selection) {
    double maxStorage;

    switch (selection) {
        case 1: {
            if (G_loc_lake[cellNumber - 1] <= 0)
                return -9999;
            else {
                // Cell-specific active lake depth - Use array (initialized from cell-specific parameter)  // FP
                maxStorage = (G_loc_lake[cellNumber - 1] / 100.0)
                             * geo.areaOfCell(cellNumber) * G_lakeDepthActive[cellNumber - 1];
                return G_locLakeStorage[cellNumber - 1] / maxStorage;
            }
        }
        case 2: {
            if (G_loc_wetland[cellNumber - 1] <= 0)
                return -9999;
            else {
                // Cell-specific active wetland depth - Use array (initialized from cell-specific parameter)
                maxStorage = (G_loc_wetland[cellNumber - 1] / 100.0)
                             * geo.areaOfCell(cellNumber) * G_wetlDepthActive[cellNumber - 1];
                return G_locWetlStorage[cellNumber - 1] / maxStorage;
            }
        }
        case 3: {
            if (G_lake_area[cellNumber - 1] <= 0)
                return -9999;
            else {
                // Cell-specific active lake depth - Use array (initialized from cell-specific parameter)
                maxStorage = G_lake_area[cellNumber - 1] * G_lakeDepthActive[cellNumber - 1];

                return G_gloLakeStorage[cellNumber - 1] / maxStorage;
            }
        }
        case 4: {
            if (G_glo_wetland[cellNumber - 1] <= 0)
                return -9999;
            else {
                // Cell-specific active wetland depth - Use array (initialized from cell-specific parameter)
                maxStorage = (G_glo_wetland[cellNumber - 1] / 100.0)
                             * geo.areaOfCell(cellNumber) * G_wetlDepthActive[cellNumber - 1];
                return G_gloWetlStorage[cellNumber - 1] / maxStorage;
            }
        }
        case 5: {
            if (G_reservoir_area[cellNumber - 1] <= 0)
                return -9999;
            else {
                maxStorage = G_stor_cap[cellNumber - 1];

                return G_gloResStorage[cellNumber - 1] / maxStorage;
            }
        }

        default:
            cerr << "Selection not available\n";
            return -9999;
    }
}


// Reservoir operation start years
void routingClass::writeAnnualReservoirCapacityandLandStorageChange(const std::string outputDirectory, const short year) {



    // Reservoir capacity of current year
    if (options.outResCapacity) {
        G_actualStorageCapacity.write(outputDirectory + "/G_RESERVOIR_STORAGE_CAPACITY_km3_" + std::to_string(year) + ".UNF0");
    }

    // Land storage change
    if (options.outLandStorageChange) {

        Grid<float> G_array;
        G_landStorageChange.write(outputDirectory + "/G_LAND_STORAGE_CHANGE_km3_" + std::to_string(year) + ".UNF0");

        for (int n = 0; n < ng; n++) {
            // for output in mm
            if (geo.G_contfreq[n] == 0) {
                G_array[n] = 0.;
            } else {
                G_array[n] = (float) G_landStorageChange[n]// conversion double -> float
                             / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000;
            }
        }
        G_array.write(outputDirectory + "/G_LAND_STORAGE_CHANGE_mm_" + std::to_string(year) + ".UNF0");

        //Now, set back G_landStorageChange to 0 (to avoid that this value occurs now in each year)
        G_landStorageChange.fill(0.);
    }
}

void routingClass::storeSingleCellDailyValuesToFile(const short year)
{
    // Storage of single cell values

    // store daily values of single cells defined in SINGLECELLS.DAT
    FILE *file_ptr;
    char filename[250];
//    extern optionClass options;
//    extern dailyWaterBalanceClass dailyWaterBalance;
    int actualSingleCellNumber;
    short calibcounter = calibGamma.getCallCounter();


    for (short i = 0; i < options.numberOfSingleCells; ++i) {
        actualSingleCellNumber = options.cellnr[i] - 1;
        sprintf(filename, "%s/SINGLE_CELL_DAILY_VALUES_%d_%d.%s", options.output_dir.c_str(), year, calibcounter,
                options.cellname[i]);
        file_ptr = fopen(filename, "w");
        fprintf(file_ptr, "# Daily values of single cell '%s':\n", options.cellname[i]);
        fprintf(file_ptr, "Day\t");
        if (options.scoutTemp) fprintf(file_ptr, "Temperature\t");
        if (options.scoutExtRad) fprintf(file_ptr, "ExtRad\t");
        if (options.scoutShortDownRad) fprintf(file_ptr, "ShortDownRad\t");
        if (options.scoutAlbedo) fprintf(file_ptr, "Albedo\t");
        if (options.scoutShortUpRad) fprintf(file_ptr, "ShortUpRad\t");
        if (options.scoutNetShortRad) fprintf(file_ptr, "NetShortRad\t");
        if (options.scoutLongDownRad) fprintf(file_ptr, "LongDownRad\t");
        if (options.scoutLongUpRad) fprintf(file_ptr, "LongUpRad\t");
        if (options.scoutNetLongRad) fprintf(file_ptr, "NetLongRad\t");
        if (options.scoutNetRad) fprintf(file_ptr, "NetRad\t");
        if (options.scoutLAI) fprintf(file_ptr, "LAI\t");
        if (options.scoutKc) fprintf(file_ptr, "Kc\t");
        if (options.scoutLandPET) fprintf(file_ptr, "LandPET\t");
        if (options.scoutCellPET) fprintf(file_ptr, "CellPET\t");
        if (options.scoutPrecip) fprintf(file_ptr, "Precip\t");
        if (options.scoutcPrecip) fprintf(file_ptr, "cPrecip\t");
        if (options.scoutCanopyWater) fprintf(file_ptr, "CanopyWater\t");
        if (options.scoutmaxCanopyWater) fprintf(file_ptr, "maxCanopyWater\t");
        if (options.scoutInterception) fprintf(file_ptr, "Interception\t");
        if (options.scoutSnowfall) fprintf(file_ptr, "Snowfall\t");
        if (options.scoutSnowCov) fprintf(file_ptr, "SnowCov\t");
        if (options.scoutSnowmelt) fprintf(file_ptr, "Snowmelt\t");
        if (options.scoutSnowWater) fprintf(file_ptr, "SnowWater\t");
        if (options.scoutSoilWater) fprintf(file_ptr, "SoilWater\t");
        if (options.scoutSurfaceRunoff) fprintf(file_ptr, "SurfaceRunoff\t");
        if (options.scoutGwRunoff) fprintf(file_ptr, "GwRunoff\t");
        if (options.scoutGwRecharge) fprintf(file_ptr, "GwRecharge\t");
        if (options.scoutCellAET) fprintf(file_ptr, "CellAET\t");
        if (options.scoutCellRunoffkm3) fprintf(file_ptr, "CellRunoffkm3\t");
        if (options.scoutCellRunoffmm) fprintf(file_ptr, "CellRunoffmm\t");
        if (options.scoutCellSRunoff) fprintf(file_ptr, "CellSRunoff\t");
        if (options.scoutQ) fprintf(file_ptr, "Q\t");
        if (options.scoutFlowVelo) fprintf(file_ptr, "FlowVelo\t");
        if (options.scoutLocLake) fprintf(file_ptr, "LocLake\t");
        if (options.scoutLocWet) fprintf(file_ptr, "LocWet\t");
        if (options.scoutGloLake) fprintf(file_ptr, "GloLake\t");
        if (options.scoutReservoir) fprintf(file_ptr, "Reservoir\t");
        if (options.scoutGloWet) fprintf(file_ptr, "GloWet\t");
        if (options.scoutRiver) fprintf(file_ptr, "River\t");
        if (options.scoutSurfaceStor) fprintf(file_ptr, "SurfaceStor\t");
        if (options.scoutTWSkm3) fprintf(file_ptr, "TWSkm3\t");
        if (options.scoutTWSmm) fprintf(file_ptr, "TWSmm\t");
        if (options.scoutGwStor) fprintf(file_ptr, "GwStor\t");
        if (options.scoutCanopyStor) fprintf(file_ptr, "CanopyStor\t");
        if (options.scoutSnowStor) fprintf(file_ptr, "SnowStor\t");
        if (options.scoutSoilStor) fprintf(file_ptr, "SoilStor\t");
        if (options.scoutGwrSwb) fprintf(file_ptr, "GwrSwb\t");
        if (options.scoutFswb) fprintf(file_ptr, "Fswb\t");
        if (options.scoutLandAreaFrac) fprintf(file_ptr, "LandAreaFrac\t");
        fprintf(file_ptr, "Day\n");

        for (short day = 0; day <= 364; day++) {
            fprintf(file_ptr, "%3d\t", day + 1);

            if (options.scoutTemp)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365TempC(actualSingleCellNumber,day));
            if (options.scoutExtRad)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365ExtRad(actualSingleCellNumber,day));
            if (options.scoutShortDownRad)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365ShortDownRad(actualSingleCellNumber,day));
            if (options.scoutAlbedo)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365Albedo(actualSingleCellNumber,day));
            if (options.scoutShortUpRad)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365ShortUpRad(actualSingleCellNumber,day));
            if (options.scoutNetShortRad)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365NetShortRad(actualSingleCellNumber,day));
            if (options.scoutLongDownRad)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365LongDownRad(actualSingleCellNumber,day));
            if (options.scoutLongUpRad)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365LongUpRad(actualSingleCellNumber,day));
            if (options.scoutNetLongRad)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365NetLongRad(actualSingleCellNumber,day));
            if (options.scoutNetRad)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365NetRadiation(actualSingleCellNumber,day));
            if (options.scoutLAI)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365LAI(actualSingleCellNumber,day));
            if (options.scoutKc)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365Kc(actualSingleCellNumber,day));
            if (options.scoutLandPET)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365PET(actualSingleCellNumber,day));
            if (options.scoutCellPET)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365TotalPET(actualSingleCellNumber,day));
            if (options.scoutPrecip)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365Precip(actualSingleCellNumber,day));
            if (options.scoutcPrecip)
                fprintf(file_ptr, "%10.6f\t", G_daily365ConsistentPrecip(actualSingleCellNumber,day));
            if (options.scoutCanopyWater)
                fprintf(file_ptr, "%10.6f\t",
                        dailyWaterBalance.G_daily365canopyWaterContent(actualSingleCellNumber,day));
            if (options.scoutmaxCanopyWater)
                fprintf(file_ptr, "%10.6f\t",
                        dailyWaterBalance.G_daily365maxcanopyWaterContent(actualSingleCellNumber,day));
            if (options.scoutInterception)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365Interception(actualSingleCellNumber,day));
            if (options.scoutSnowfall)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365SnowFall(actualSingleCellNumber,day));
            if (options.scoutSnowCov)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365SnowCoverAvg(actualSingleCellNumber,day));
            if (options.scoutSnowmelt)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365SnowMelt(actualSingleCellNumber,day));
            if (options.scoutSnowWater)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365SWE(actualSingleCellNumber,day));
            if (options.scoutSoilWater)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365SoilWaterAvg(actualSingleCellNumber,day));
            if (options.scoutSurfaceRunoff)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365SurfaceRunoff(actualSingleCellNumber,day));
            if (options.scoutGwRunoff) fprintf(file_ptr, "%10.6f\t", G_daily365GwRunoff(actualSingleCellNumber,day));
            if (options.scoutGwRecharge)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365GwRecharge(actualSingleCellNumber,day));
            if (options.scoutCellAET) fprintf(file_ptr, "%10.6f\t", G_daily365CellAET(actualSingleCellNumber,day));
            if (options.scoutCellRunoffkm3)
                fprintf(file_ptr, "%10.6f\t", G_daily365CellRunoff(actualSingleCellNumber,day));
            if (options.scoutCellRunoffmm)
                fprintf(file_ptr, "%10.6f\t", G_daily365CellRunoff_mm(actualSingleCellNumber,day));
            if (options.scoutCellSRunoff)
                fprintf(file_ptr, "%10.6f\t", G_daily365CellSurfaceRunoff(actualSingleCellNumber,day));
            if (options.scoutQ) fprintf(file_ptr, "%10.6f\t", G_daily365RiverAvail(actualSingleCellNumber,day));
            if (options.scoutFlowVelo) fprintf(file_ptr, "%10.6f\t", G_daily365Velocity(actualSingleCellNumber,day));
            if (options.scoutLocLake) fprintf(file_ptr, "%10.6f\t", G_daily365LocLakeStor(actualSingleCellNumber,day));
            if (options.scoutLocWet) fprintf(file_ptr, "%10.6f\t", G_daily365LocWetlStor(actualSingleCellNumber,day));
            if (options.scoutGloLake) fprintf(file_ptr, "%10.6f\t", G_daily365GloLakeStor(actualSingleCellNumber,day));
            if (options.scoutReservoir) fprintf(file_ptr, "%10.6f\t", G_daily365ResStor(actualSingleCellNumber,day));
            if (options.scoutGloWet) fprintf(file_ptr, "%10.6f\t", G_daily365GloWetlStor(actualSingleCellNumber,day));
            if (options.scoutRiver) fprintf(file_ptr, "%10.6f\t", G_daily365RiverStor(actualSingleCellNumber,day));
            if (options.scoutSurfaceStor)
                fprintf(file_ptr, "%10.6f\t", G_daily365SurfStor(actualSingleCellNumber,day));
            if (options.scoutTWSkm3)
                fprintf(file_ptr, "%10.6f\t", G_daily365TotalWaterInStorages_km3(actualSingleCellNumber,day));
            if (options.scoutTWSmm)
                fprintf(file_ptr, "%10.6f\t", G_daily365TotalWaterInStorages_mm(actualSingleCellNumber,day));
            if (options.scoutGwStor) fprintf(file_ptr, "%10.6f\t", G_daily365GwStor(actualSingleCellNumber,day));
            if (options.scoutCanopyStor)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365CanopyStorage(actualSingleCellNumber,day));
            if (options.scoutSnowStor)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365SnowStorage(actualSingleCellNumber,day));
            if (options.scoutSoilStor)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365SoilStorage(actualSingleCellNumber,day));
            if (options.scoutGwrSwb) fprintf(file_ptr, "%10.6f\t", G_daily365GwrSwb(actualSingleCellNumber,day));
            if (options.scoutFswb) fprintf(file_ptr, "%10.6f\t", G_daily365Fswb(actualSingleCellNumber,day));
            if (options.scoutLandAreaFrac)
                fprintf(file_ptr, "%10.6f\t", G_daily365LandAreaFrac(actualSingleCellNumber,day));
            fprintf(file_ptr, "%3d\n", day + 1);
        }
        fclose(file_ptr);
    }
}

void routingClass::writeLakeStorageGrids(const std::string outputDirectory, int year) {

    // values in [km3]  // yearly outputs (1 == options.grid_store)

    G_locLakeStorage.write(outputDirectory + "/G_LOC_LAKE_STORAGE_km3_" + std::to_string(year) + ".UNF0");

    G_locWetlStorage.write(outputDirectory + "/G_LOC_WETL_STORAGE_km3_" + std::to_string(year) + ".UNF0");

    G_gloLakeStorage.write(outputDirectory + "/G_GLO_LAKE_STORAGE_km3_" + std::to_string(year) + ".UNF0");

    G_gloWetlStorage.write(outputDirectory + "/G_GLO_WETL_STORAGE_km3_" + std::to_string(year) + ".UNF0");

    if (options.resOpt == 1) {
        G_gloResStorage.write(outputDirectory + "/G_RES_STORAGE_km3_" + std::to_string(year) + ".UNF0");
    }
}

void routingClass::writeLakeStorRatioGrids(const std::string outputDirectory, int year) {
    // range of values: 1 is maximum, yearly outputs (1 == options.grid_store)

    Grid<float> G_storageRatio;
    int n;

    // local lakes
    for (n = 0; n <= ng - 1; n++) {
        G_storageRatio[n] = G_locLakeStorage[n] / ((G_loc_lake[n] / 100.0)
                                                   * geo.areaOfCellByArrayPos(n) * G_lakeDepthActive[n]);
        // term in brackets is equal to maxStorage
    }
    G_storageRatio.write(outputDirectory + "/G_LOC_LAKE_STOR_RATIO_" + std::to_string(year) + ".UNF0");

    // local wetlands
    for (int n = 0; n <= ng - 1; n++) {
        // Cell-specific active wetland depth - Use array (initialized from cell-specific parameter)
        G_storageRatio[n] = G_locWetlStorage[n] / ((G_loc_wetland[n] / 100.0)
                                                   * geo.areaOfCellByArrayPos(n) * G_wetlDepthActive[n]);
    }
    G_storageRatio.write(outputDirectory + "/G_LOC_WETL_STOR_RATIO_" + std::to_string(year) + ".UNF0");

    // global lakes
    for (int n = 0; n <= ng - 1; n++) {
        G_storageRatio[n] = G_gloLakeStorage[n]
                            // Cell-specific active lake depth - Use array (initialized from cell-specific parameter)
                            / ((double) G_lake_area[n] * G_lakeDepthActive[n]);


    }
    G_storageRatio.write(outputDirectory + "/G_GLO_LAKE_STOR_RATIO_" + std::to_string(year) + ".UNF0");

    // global wetlands
    for (int n = 0; n <= ng - 1; n++) {
        // Cell-specific active wetland depth - Use array (initialized from cell-specific parameter)
        G_storageRatio[n] = G_gloWetlStorage[n] / ((G_glo_wetland[n] / 100.0)
                                                   * geo.areaOfCellByArrayPos(n) * G_wetlDepthActive[n]);
    }
    G_storageRatio.write(outputDirectory + "/G_GLO_WETL_STOR_RATIO_" + std::to_string(year) + ".UNF0");

    // reservoirs
    if (options.resOpt == 1) {
        for (int n = 0; n <= ng - 1; n++) {
            G_storageRatio[n] = G_gloResStorage[n] / G_stor_cap[n];
        }
        G_storageRatio.write(outputDirectory + "/G_RES_STOR_RATIO_" + std::to_string(year) + ".UNF0");
    }
}

void routingClass::writeMonthlyGrids(const std::string outputDirectory, int year, int month2store) {
    if (2 != options.grid_store && 4 != options.grid_store && 6 != options.grid_store)
        return;

    int n;
    short month;
    MonthlyGrid<> G_monthlyArray;
    int numberOfDaysInMonth[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    if (options.outRiverAvail) {
        G_monthlyRiverAvail.write(outputDirectory + "/G_RIVER_AVAIL_" + std::to_string(year) + ".12.UNF0");
    }

    if (options.outPrec) {
        //  Precip which is really used by WaterGAP
        G_monthlyConsistentPrecip.write(outputDirectory + "/G_CONSISTENT_PRECIPITATION_km3_" + std::to_string(year) + ".12.UNF0");
    }

    // output option min/max daily discharge
    if (options.outMinMaxRiverAvail) {
        G_monthlyMinRiverAvail.write(outputDirectory + "/G_MIN_RIVER_AVAIL_" + std::to_string(year) + ".12.UNF0");

        G_monthlyMaxRiverAvail.write(outputDirectory + "/G_MAX_RIVER_AVAIL_" + std::to_string(year) + ".12.UNF0");
    }

    if (options.outRiverPET) {
        for (n = 0; n < ng; n++)
            for (month = 0; month < 12; month++)
                G_monthlyArray(n,month) = G_monthlyRiverAreaFrac(n,month) / numberOfDaysInMonth[month];

        G_monthlyArray.write(outputDirectory + "/G_RIVER_FRACTIONS_" + std::to_string(year) + ".12.UNF0");

        G_monthlyRiverPET.write(outputDirectory + "/G_RIVER_PET_" + std::to_string(year) + ".12.UNF0");

    }

    if (options.outRiverInUpstream) {
        G_monthlyRiverInUpstream.write(outputDirectory + "/G_RIVER_IN_UPSTREAM_" + std::to_string(year) + ".12.UNF0");
    }

    if (options.outRiverVelo) {
        for (n = 0; n < ng; n++)
            for (month = 0; month < 12; month++)
                G_monthlyArray(n,month) = G_monthlyVelocity(n,month) / numberOfDaysInMonth[month];

        G_monthlyArray.write(outputDirectory + "/G_RIVER_VELOCITY_" + std::to_string(year) + ".12.UNF0");
    }

    if (options.outLocWetlExt) {
        for (n = 0; n < ng; n++)
            for (month = 0; month < 12; month++)
                G_monthlyArray(n,month) = G_monthlyLocWetlExtent(n,month) / numberOfDaysInMonth[month];

        G_monthlyArray.write(outputDirectory + "/G_LOC_WETL_EXTENT_km2_" + std::to_string(year) + ".12.UNF0");
    }
    if (options.outGloWetlExt) {
        for (n = 0; n < ng; n++)
            for (month = 0; month < 12; month++)
                G_monthlyArray(n,month) = G_monthlyGloWetlExtent(n,month) / numberOfDaysInMonth[month];

        G_monthlyArray.write(outputDirectory + "/G_GLO_WETL_EXTENT_km2_" + std::to_string(year) + ".12.UNF0");
    }
    if (options.outCellRunoff) {

        G_monthlyCellRunoff.write(outputDirectory + "/G_CELL_RUNOFF_" + std::to_string(year) + ".12.UNF0");

        for (n = 0; n < ng; n++) { // total water resources in mm/month
            for (month = 0; month < 12; month++) {
                // Only continental cells
                if (geo.G_contcell[n]) {
                    G_monthlyArray(n,month) = G_monthlyCellRunoff(n,month)
                                           / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000;
                }
                // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0.
                // The water balance is calculated in the outlet cell of large lakes; therefore a calculation of
                // e.g. precipitation for this cells leads to inconsistent quantification of global fluxes.
                // A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                else {
                    G_monthlyArray(n,month) = 0.;
                }
            }
        }

        G_monthlyArray.write(outputDirectory + "/G_CELL_RUNOFF_mm_" + std::to_string(year) + ".12.UNF0");

        for (n = 0; n < ng; n++)
            for (month = 0; month < 12; month++)
                G_monthlyPotCellRunoff(n,month) -= G_monthlyPotCellRunoffDeficit(n,month);    // contains negative values

        G_monthlyPotCellRunoff.write(outputDirectory + "/G_CELL_RUNOFF_TOTAL_" + std::to_string(year) + ".12.UNF0");
    }

    if (options.outGwrunSurfrun) {
        G_monthlyGwrunSurfrun.write(outputDirectory + "/G_GWRUN_SURFRUN_km3_" + to_string(year) + ".12.UNF0");
    }

    if (options.outCellSurface) {
        G_monthlyCellSurfaceRunoff.write(outputDirectory + "/G_CELL_SURFACE_RUNOFF_" + std::to_string(year) + ".12.UNF0");
    }

    if (options.outSurfStor) {

        if (options.grid_store_TypeForStorages == 0) {

            G_monthlySurfStor.write(outputDirectory + "/G_SURFACE_WATER_STORAGE_km3_" + std::to_string(year) + ".12.UNF0");

            for (n = 0; n < ng; n++) {
                for (month = 0; month < 12; month++) {
                    // Only continental cells
                    if (geo.G_contcell[n]) {
                        G_monthlyArray(n,month) =
                                       (G_monthlySurfStor(n,month)
                                        / (geo.areaOfCellByArrayPos(n) * geo.G_contfreq[n] / 100.)) * 1000000;
                    }
                        // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0.
                        // The water balance is calculated in the outlet cell of large lakes; therefore a calculation of
                        // e.g. precipitation for this cells leads to inconsistent quantification of global fluxes.
                        // A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                    else {
                        G_monthlyArray(n,month) = 0.;
                    }
                }
            }
            G_monthlyArray.write(outputDirectory + "/G_SURFACE_WATER_STORAGE_mm_" + std::to_string(year) + ".12.UNF0");
        }
        else {
            for (n = 0; n < ng; n++)
                for (month = 0; month < 12; month++)
                    G_monthlyArray(n,month) = G_monthlySurfStor(n,month) / numberOfDaysInMonth[month];

            G_monthlyArray.write(outputDirectory + "/G_SURFACE_WATER_STORAGE_MEAN_km3_" + std::to_string(year) + ".12.UNF0");

            for (n = 0; n < ng; n++) {
                for (month = 0; month < 12; month++) {
                    // Only continental cells
                    if (geo.G_contcell[n]) {
                        G_monthlyArray(n,month) =
                                       (G_monthlySurfStor(n,month) / numberOfDaysInMonth[month]
                                        / (geo.areaOfCellByArrayPos(n) * geo.G_contfreq[n] / 100.)) * 1000000;
                    }
                    // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0.
                    // The water balance is calculated in the outlet cell of large lakes; therefore a calculation of
                    // e.g. precipitation for this cells leads to inconsistent quantification of global fluxes.
                    // A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                    else {
                        G_monthlyArray(n,month) = 0.;
                    }
                }
            }
            G_monthlyArray.write(outputDirectory + "/G_SURFACE_WATER_STORAGE_MEAN_mm_" + std::to_string(year) + ".12.UNF0");
        }

    }

    if (options.outCellAET) {
        G_monthlyCellAET.write(outputDirectory + "/G_CELL_AET_" + std::to_string(year) + ".12.UNF0");
    }

    if (options.outOpenWaterEvap) {
        G_monthlyOpenWaterEvap.write(outputDirectory + "/G_OPEN_WATER_EVAP_" + std::to_string(year) + ".12.UNF0");
    }

    if (options.outGwrSwb) {

        // monthly sum of groundwater recharge below surface water bodies (with respect to continental area) [mm/month]
        for (n = 0; n < ng; n++)
            for (month = 0; month < 12; month++)
                G_monthlyArray(n,month) = G_monthlyGwrSwb(n,month);

        G_monthlyArray.write(outputDirectory + "/G_GWR_SURFACE_WATER_BODIES_mm_" + std::to_string(year) + ".12.UNF0");

        for (n = 0; n < ng; n++)
            for (month = 0; month < 12; month++)
                G_monthlyArray(n,month) = G_monthlyGwrSwb(n,month)
                                       * (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) / 1000000.;

        G_monthlyArray.write(outputDirectory + "/G_GWR_SURFACE_WATER_BODIES_km3_" + std::to_string(year) + ".12.UNF0");
    }

    if (options.outTotalGWR) {
        for (n = 0; n < ng; n++)
            for (month = 0; month < 12; month++)
                G_monthlyArray(n,month) = G_monthlyGwrSwb(n,month) + dailyWaterBalance.G_monthlyGwRecharge(n,month);

        G_monthlyArray.write(outputDirectory + "/G_TOTAL_GW_RECHARGE_mm_" + std::to_string(year) + ".12.UNF0");

        for (n = 0; n < ng; n++)
            for (month = 0; month < 12; month++)
                G_monthlyArray(n,month) = (G_monthlyGwrSwb(n,month) + dailyWaterBalance.G_monthlyGwRecharge(n,month))
                                       * (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) / 1000000.;

        G_monthlyArray.write(outputDirectory + "/G_TOTAL_GW_RECHARGE_km3_" + std::to_string(year) + ".12.UNF0");
    }

    if (options.outFswb) {

        // mean monthly fraction of surface water bodies
        for (n = 0; n < ng; n++)
            for (month = 0; month < 12; month++)
                G_monthlyArray(n,month) = G_monthlyFswb(n,month) / numberOfDaysInMonth[month];

        G_monthlyArray.write(outputDirectory + "/G_FRACTION_SURFACE_WATER_BODIES_" + std::to_string(year) + ".12.UNF0");
    }

    if (options.outLandAreaFrac) {

        // Write land area fractions
        for (n = 0; n < ng; n++)
            for (month = 0; month < 12; month++)
                G_monthlyArray(n,month) = G_monthlyLandAreaFrac(n,month) / numberOfDaysInMonth[month];

        G_monthlyArray.write(outputDirectory + "/G_LAND_AREA_FRACTIONS_" + std::to_string(year) + ".12.UNF0");
    }

    // GFZ:GLOLAK begin - Allocation of GloLakeStorage of Outflowcell to all Lakecells
    //for lakes and reservoirs
    if (0 != ng_glolakcells) {

        Grid<float,ng_glolakcells> G_tempGloLakeStorage;

        long ofcell, acell;
        for (month = 0; month < 12; month++) {
            //Initialisation of the temporary Glolakstorage-array
            for (n = 0; n < ng_glolakcells; n++) {
                G_tempGloLakeStorage[n] = 0.0;
            }
            //computation of glolakstorage for each lakecell and saving in temporary array
            for (n = 0; n < ng_glolakcells; n++) {
                ofcell = big_lakes_cells(0,n);
                acell = big_lakes_cells(1,n);
                if (ofcell != -99)
                    G_tempGloLakeStorage[n] += G_monthlyGloLakeStorage(acell - 1,month) * big_lakes_area_frac[n];
            }
            //Rewriting of global array for Glolakstorage
            for (n = 0; n < ng_glolakcells; n++) {
                ofcell = big_lakes_cells(0,n);
                G_monthlyGloLakeStorage(ofcell - 1,month) = G_tempGloLakeStorage[n];
            }
        }
    }

    // monthly storage grids
    if (options.grid_store_TypeForStorages == 1) {
        for (n = 0; n <= ng - 1; n++) {
            for (month = 0; month < 12; month++) {
                G_monthlyGwStorage(n,month) /= numberOfDaysInMonth[month];
                G_monthlyLocLakeStorage(n,month) /= numberOfDaysInMonth[month];
                G_monthlyLocWetlStorage(n,month) /= numberOfDaysInMonth[month];
                G_monthlyGloLakeStorage(n,month) /= numberOfDaysInMonth[month];
                G_monthlyGloWetlStorage(n,month) /= numberOfDaysInMonth[month];
                G_monthlyRiverStorage(n,month) /= numberOfDaysInMonth[month];
                if (options.resOpt == 1)
                    G_monthlyResStorage(n,month) /= numberOfDaysInMonth[month];
                if (options.glacierOpt == 1)
                    G_monthlyGlacierStorage(n,month) /= numberOfDaysInMonth[month];
            }
        }
    }

    if (options.outGWRunoff) {

        G_monthlyGwRunoff.write(outputDirectory + "/G_GW_RUNOFF_km3_" + std::to_string(year) + ".12.UNF0");

        for (n = 0; n < ng; n++) { // for output in mm
            for (month = 0; month < 12; month++) {
                // Only continental cells
                if (geo.G_contcell[n]) {
                    G_monthlyArray(n,month) = G_monthlyGwRunoff(n,month)
                                           / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000.0;
                }
                // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0.
                // The water balance is calculated in the outlet cell of large lakes; therefore a calculation of
                // e.g. precipitation for this cells leads to inconsistent quantification of global fluxes.
                // A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                else {
                    G_monthlyArray(n,month) = 0.;
                }
            }
        }
        G_monthlyArray.write(outputDirectory + "/G_GW_RUNOFF_mm_" + std::to_string(year) + ".12.UNF0");
    }

    if (options.outGlacArea)
        G_monthlyGlacierArea.write(outputDirectory + "/G_GLACIER_AREA_km2_" + std::to_string(year) + ".12.UNF0");

    if (options.outGlacAreaFrac)
        G_monthlyGlacierAreaFrac.write(outputDirectory + "/G_GLACIER_AREA_FRACTION_" + std::to_string(year) + ".12.UNF0");

    if (options.outGlacPrecip) {

        G_monthlyGlacierPrecip.write(outputDirectory + "/G_PRECIPITATION_ON_GLACIER_km3_" + std::to_string(year) + ".12.UNF0");

        for (n = 0; n < ng; n++) {
            for (month = 0; month < 12; month++) {
                if (geo.G_contcell[n])
                    G_monthlyArray(n,month) = G_monthlyGlacierPrecip(n,month) / G_glacAdaptedArea[n] * 1000000.0;  // for output in mm
                else
                    G_monthlyArray(n,month) = 0.;
            }
        }
        G_monthlyArray.write(outputDirectory + "/G_PRECIPITATION_ON_GLACIER_mm_" + std::to_string(year) + ".12.UNF0");
    }

    if (options.outGlacRunoff) {

        G_monthlyGlacierRunoff.write(outputDirectory + "/G_GLACIER_RUNOFF_km3_" + std::to_string(year) + ".12.UNF0");

        for (n = 0; n < ng; n++) {
            for (month = 0; month < 12; month++) {
                if (geo.G_contcell[n])
                    G_monthlyArray(n,month) = G_monthlyGlacierRunoff(n,month)
                                           / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000.0;  // for output in mm
                else
                    G_monthlyArray(n,month) = 0.;
            }
        }
        G_monthlyArray.write(outputDirectory + "/G_GLACIER_RUNOFF_mm_" + std::to_string(year) + ".12.UNF0");
    }

    //write output of monthly storage (components)

    if (options.outLocLakeStorage || options.outSingleStorages) {

        G_monthlyGwStorage.write(outputDirectory + "/G_GROUND_WATER_STORAGE_km3_" + std::to_string(year) + ".12.UNF0");

        for (n = 0; n < ng; n++) {
            // for output in mm
            for (month = 0; month < 12; month++) {
                // Only continental cells
                if (geo.G_contcell[n]) {
                    G_monthlyArray(n,month) =
                                   G_monthlyGwStorage(n,month)
                                   / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000.0;
                }
                // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0.
                // The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation
                // for this cells leads to inconsistent quantification of global fluxes.
                // A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                else {
                    G_monthlyArray(n,month) = 0.;
                }
            }
        }

        G_monthlyArray.write(outputDirectory + "/G_GROUND_WATER_STORAGE_mm_" + std::to_string(year) + ".12.UNF0");
    }
    if (options.outLocLakeStorage || options.outSingleStorages) {

        // LOCAL LAKES
        G_monthlyLocLakeStorage.write(outputDirectory + "/G_LOC_LAKE_STORAGE_km3_" + std::to_string(year) + ".12.UNF0");

        for (n = 0; n < ng; n++) {
            // for output in mm
            for (month = 0; month < 12; month++) {
                // Only continental cells
                if (geo.G_contcell[n]) {
                    G_monthlyArray(n,month) = G_monthlyLocLakeStorage(n,month)
                                           / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000;
                }
                // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0.
                // The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation
                // for this cells leads to inconsistent quantification of global fluxes.
                // A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                else {
                    G_monthlyArray(n,month) = 0.;
                }
            }
        }
        G_monthlyArray.write(outputDirectory + "/G_LOC_LAKE_STORAGE_mm_" + std::to_string(year) + ".12.UNF0");
    }

    if (options.outLocWetlStorage || options.outSingleStorages) {

        // LOCAL WETLANDS
        G_monthlyLocWetlStorage.write(outputDirectory + "/G_LOC_WETL_STORAGE_km3_" + std::to_string(year) + ".12.UNF0");

        for (n = 0; n < ng; n++) {
            // for output in mm
            for (month = 0; month < 12; month++) {
                // Only continental cells
                if (geo.G_contcell[n]) {
                    G_monthlyArray(n,month) = G_monthlyLocWetlStorage(n,month)
                                           / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000;
                }
                // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0.
                // The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation
                // for this cells leads to inconsistent quantification of global fluxes.
                // A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                else {
                    G_monthlyArray(n,month) = 0.;
                }
            }
        }
        G_monthlyArray.write(outputDirectory + "/G_LOC_WETL_STORAGE_mm_" + std::to_string(year) + ".12.UNF0");
    }

    if (options.outGloLakeStorage || options.outSingleStorages) {

        // GLOBAL LAKES
        G_monthlyGloLakeStorage.write(outputDirectory + "/G_GLO_LAKE_STORAGE_km3_" + std::to_string(year) + ".12.UNF0");

        for (n = 0; n < ng; n++) {
            // for output in mm
            for (month = 0; month < 12; month++) {
                // Only continental cells
                if (geo.G_contcell[n]) {
                    G_monthlyArray(n,month) = G_monthlyGloLakeStorage(n,month)
                                           / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000;
                }
                // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0.
                // The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation
                // for this cells leads to inconsistent quantification of global fluxes.
                // A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                else {
                    G_monthlyArray(n,month) = 0.;
                }
            }
        }
        G_monthlyArray.write(outputDirectory + "/G_GLO_LAKE_STORAGE_mm_" + std::to_string(year) + ".12.UNF0");
    }

    if ((options.outResStorage || options.outSingleStorages) && options.resOpt == 1) {

        // RESERVOIRS
        std::string filename;
        if (options.grid_store_TypeForStorages == 0)
            filename = outputDirectory + "/G_RES_STORAGE_km3_" + std::to_string(year) + ".12.UNF0";
        else
            filename = outputDirectory + "/G_RES_STORAGE_MEAN_km3_" + std::to_string(year) + ".12.UNF0";

        G_monthlyResStorage.write(filename);

        for (n = 0; n < ng; n++) {
            // for output in mm
            for (month = 0; month < 12; month++) {
                // Only continental cells
                if (geo.G_contcell[n]) {
                    G_monthlyArray(n,month) = G_monthlyResStorage(n,month)
                                           / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000;
                }
                // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0.
                // The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation
                // for this cells leads to inconsistent quantification of global fluxes.
                // A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                else {
                    G_monthlyArray(n,month) = 0.;
                }
            }
        }

        if (options.grid_store_TypeForStorages == 0)
            filename = outputDirectory + "/G_RES_STORAGE_mm_" + std::to_string(year) + ".12.UNF0";
        else
            filename = outputDirectory + "/G_RES_STORAGE_MEAN_mm_" + std::to_string(year) + ".12.UNF0";

        G_monthlyArray.write(filename);
    }

    if ((options.outResStorage || options.outSingleStorages) && options.resOpt == 1) {

        // RESERVOIRS
        for (n = 0; n < ng; n++)
            for (month = 0; month < 12; month++)
                G_monthlyArray(n,month) = G_monthlyResOutflow(n,month) /
                                       1000000000.; // [m3/month] -> [km3/month]
        G_monthlyArray.write(outputDirectory + "/G_RES_OUT_km3_" + std::to_string(year) + ".12.UNF0");
    }

    if (options.outGloWetlStorage || options.outSingleStorages) {

        // GLOBAL WETLANDS
        G_monthlyGloWetlStorage.write(outputDirectory + "/G_GLO_WETL_STORAGE_km3_" + std::to_string(year) + ".12.UNF0");

        for (n = 0; n < ng; n++) {
            // for output in mm
            for (month = 0; month < 12; month++) {
                // Only continental cells
                if (geo.G_contcell[n]) {
                    G_monthlyArray(n,month) = G_monthlyGloWetlStorage(n,month)
                                           / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000;
                }
                // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                else {
                    G_monthlyArray(n,month) = 0.;
                }
            }
        }
        G_monthlyArray.write(outputDirectory + "/G_GLO_WETL_STORAGE_mm_" + std::to_string(year) + ".12.UNF0");
    }

    if (options.outRiverStorage || options.outSingleStorages) {

        // RIVER STORAGE
        G_monthlyRiverStorage.write(outputDirectory + "/G_RIVER_STORAGE_km3_" + std::to_string(year) + ".12.UNF0");

        for (n = 0; n < ng; n++) {
            // for output in mm
            for (month = 0; month < 12; month++) {
                // Only continental cells
                if (geo.G_contcell[n]) {
                    G_monthlyArray(n,month) = G_monthlyRiverStorage(n,month)
                                           / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000;
                }
                // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                else {
                    G_monthlyArray(n,month) = 0.;
                }
            }
        }
        G_monthlyArray.write(outputDirectory + "/G_RIVER_STORAGE_mm_" + std::to_string(year) + ".12.UNF0");
    }

    if ((options.outGlacierStorage || options.outSingleStorages) && (options.glacierOpt == 1)) {

        G_monthlyGlacierStorage.write(outputDirectory + "/G_GLACIER_STORAGE_km3_" + std::to_string(year) + ".12.UNF0");

        for (n = 0; n < ng; n++) {
            for (month = 0; month < 12; month++) {
                if (geo.G_contcell[n])
                    G_monthlyArray(n,month) = G_monthlyGlacierStorage(n,month) / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000.0;  // for output in mm
                else
                    G_monthlyArray(n,month) = 0.;
            }
        }
        G_monthlyArray.write(outputDirectory + "/G_GLACIER_STORAGE_mm_" + std::to_string(year) + ".12.UNF0");
    }

    if (options.outActualNAS) {

        // ACTUAL USE
        for (n = 0; n < ng; n++)
            for (month = 0; month < 12; month++) {
                G_monthlyArray(n,month) = G_monthlyActualUse(n,month) * 1000000000;
            }
        G_monthlyArray.write(outputDirectory + "/G_ACTUAL_NAS_" + std::to_string(year) + ".12.UNF0");
    }

    if (options.outActualNAG) {

        // ACTUAL USE
        for (n = 0; n < ng; n++)
            for (month = 0; month < 12; month++) {
                G_monthlyArray(n,month) = G_monthlyNUg(n,month) * 1000000000;
            }
        G_monthlyArray.write(outputDirectory + "/G_ACTUAL_NAG_" + std::to_string(year) + ".12.UNF0");
    }

    if (options.outSatisAllocUsein2ndCell) {
        G_monthlySatisAllocatedUseinSecondCell.write(outputDirectory + "/G_SATIS_ALLOC_USE_IN_2NDCELL_" + std::to_string(year) + ".12.UNF0");
    }

    if (options.outWCa) {
        for (n = 0; n < ng; n++) {
            for (month = 0; month < 12; month++) {
                G_monthlyArray(n,month) = (G_monthlyActualUse(n,month) + G_monthlyNUg(n,month)) * 1000000000;
            }
        }
        G_monthlyArray.write(outputDirectory + "/G_ACTUAL_WATER_CONSUMPTION_" + std::to_string(year) + ".12.UNF0");
    }

    if (options.outWCaInCellUsed) {
        for (n = 0; n < ng; n++) {
            for (month = 0; month < 12; month++) {
                G_monthlyArray(n,month) = (G_monthlyRedistributeAllocUse(n,month) + G_monthlyRedistributeGLWDUse(n, month) + G_monthlyActualUse(n,month) + G_monthlyNUg(n,month)) * 1000000000;
            }
        }

        G_monthlyArray.write(outputDirectory + "/G_ACTUAL_WATER_CONSUMPTION_INCELLUSED_" + std::to_string(year) + ".12.UNF0");
    }

    if (options.outCellAETWCa) { // new output options 2.2b [km3/month]
        for (n = 0; n < ng; n++) {
            for (month = 0; month < 12; month++) {
                G_monthlyArray(n,month) = G_monthlyCellAETWCa(n,month) + G_monthlyActualUse(n,month) +
                                       G_monthlyNUg(n,month);
            }
        }
        G_monthlyArray.write(outputDirectory + "/G_CELLAET_CONSUSE_km3_" + std::to_string(year) + ".12.UNF0");
    }

    // G_AllocatedUse = test output: demand allocated from neighboring cell(s)
    if (options.outAllocUsein2ndCell) {
        G_monthlyAllocatedUse.write(outputDirectory + "/G_ALLOC_USE_IN_2NDCELL_" + std::to_string(year) + ".12.UNF0");
    }

    if (options.outWaterTempMonthly) {
        G_monthlyMeanWaterTempRiver.write(outputDirectory + "/G_RIVER_WATER_TEMP_" + std::to_string(year) + ".12.UNF0");
    }

    if (options.outWaterTempMonthlyAllSWB){
        G_monthlyMeanWaterTempLocLake.write(outputDirectory + "/G_LOC_LAKE_WATER_TEMP_" + std::to_string(year) + ".12.UNF0");
        G_monthlyMeanWaterTempLocWetl.write(outputDirectory + "/G_LOC_WETLAND_WATER_TEMP_" + std::to_string(year) + ".12.UNF0");
        G_monthlyMeanWaterTempGloLake.write(outputDirectory + "/G_GLO_LAKE_WATER_TEMP_" + std::to_string(year) + ".12.UNF0");
        G_monthlyMeanWaterTempReservoir.write(outputDirectory + "/G_RESERVOIR_WATER_TEMP_" + std::to_string(year) + ".12.UNF0");
        G_monthlyMeanWaterTempGloWetl.write(outputDirectory + "/G_GLO_WETLAND_WATER_TEMP_" + std::to_string(year) + ".12.UNF0");
    }

}

void routingClass::writeDaily31Grids(const std::string outputDirectory, const int year, const int month) {
    if (3 != options.grid_store && 4 != options.grid_store)
        return;

    if (options.outCellRunoffDaily) {
        G_daily31CellRunoff.write(outputDirectory + "/G_CELL_RUNOFF_" + std::to_string(year) + "_" + std::to_string(month+1) + ".31.UNF0");
    }
    if (options.outGWRunoffDaily) {
        G_daily31GwRunoff.write(outputDirectory + "/G_GW_RUNOFF_" + std::to_string(year) + "_" + std::to_string(month+1) + ".31.UNF0");
    }
    if (options.outGwrunSurfrunDaily) {
        G_daily31GwrunSurfrun.write(outputDirectory + "/G_GWRUN_SURFRUN_" + std::to_string(year) + "_" + std::to_string(month+1) + ".31.UNF0");
    }
    if (options.outPrecDaily) {
        G_daily31ConsistentPrecip.write(outputDirectory + "/G_CONSISTENT_PRECIPITATION_km3_" + std::to_string(year) + "_" + std::to_string(month+1) + ".31.UNF0");
    }
    if (options.outCellSurfaceDaily) {
        G_daily31CellSurfaceRunoff.write(outputDirectory + "/G_CELL_SURFACE_RUNOFF_" + std::to_string(year) + "_" + std::to_string(month+1) + ".31.UNF0");
    }
    if (options.outRiverAvailDaily) {
        G_daily31RiverAvail.write(outputDirectory + "/G_RIVER_AVAIL_" + std::to_string(year) + "_" + std::to_string(month+1) + ".31.UNF0");
    }
    if (options.outRiverVeloDaily) {
        G_daily31Velocity.write(outputDirectory + "/G_RIVER_VELOCITY_" + std::to_string(year) + "_" + std::to_string(month+1) + ".31.UNF0");
    }
    if (options.outWaterTempDaily) {
        G_daily31WaterTemp.write(outputDirectory + "/G_RIVER_WATER_TEMP_" + std::to_string(year) + "_" + std::to_string(month + 1) + ".31.UNF0");
    }
    if (options.outWaterTempDailyAllSWB){
        G_daily31locLakeTemp.write(outputDirectory + "/G_LOC_LAKE_WATER_TEMP_" + std::to_string(year) + "_" + std::to_string(month+1) + ".31.UNF0");
        G_daily31locWetlTemp.write(outputDirectory + "/G_LOC_WETLAND_WATER_TEMP_" + std::to_string(year) + "_" + std::to_string(month+1) + ".31.UNF0");
        G_daily31gloLakeTemp.write(outputDirectory + "/G_GLO_LAKE_WATER_TEMP_" + std::to_string(year) + "_" + std::to_string(month+1) + ".31.UNF0");
        G_daily31reservoirTemp.write(outputDirectory + "/G_RESERVOIR_WATER_TEMP_" + std::to_string(year) + "_" + std::to_string(month+1) + ".31.UNF0");
        G_daily31gloWetlandTemp.write(outputDirectory + "/G_GLO_WETLAND_WATER_TEMP_" + std::to_string(year) + "_" + std::to_string(month+1) + ".31.UNF0");
    }
    if (options.outCellAETDaily) {
        G_daily31CellAET.write(outputDirectory + "/G_CELL_AET_" + std::to_string(year) + "_" + std::to_string(month+1) + ".31.UNF0");
    }
    if (options.outCellAETWCaDaily) {
        G_daily31CellAETWCa.write(outputDirectory + "/G_CELLAET_CONSUSE_km3_" + std::to_string(year) + "_" + std::to_string(month+1) + ".31.UNF0");
    }
    if (options.outSurfStorDaily) {
        G_daily31SurfStor.write(outputDirectory + "/G_SURFACE_WATER_STORAGE_km3_" + std::to_string(year) + "_" + std::to_string(month+1) + ".31.UNF0");
    }
    if (options.outSingleStoragesDaily) {
        G_daily31LocLakeStor.write(outputDirectory + "/G_LOC_LAKE_STORAGE_km3_" + std::to_string(year) + "_" + std::to_string(month+1) + ".31.UNF0");
    }
    if (options.outSingleStoragesDaily) {
        G_daily31LocWetlStor.write(outputDirectory + "/G_LOC_WETL_STORAGE_km3_" + std::to_string(year) + "_" + std::to_string(month+1) + ".31.UNF0");
    }
    if (options.outSingleStoragesDaily) {
        G_daily31GloLakeStor.write(outputDirectory + "/G_GLO_LAKE_STORAGE_km3_" + std::to_string(year) + "_" + std::to_string(month+1) + ".31.UNF0");
    }
    if (options.outSingleStoragesDaily) {
        G_daily31GloWetlStor.write(outputDirectory + "/G_GLO_WETL_STORAGE_km3_" + std::to_string(year) + "_" + std::to_string(month+1) + ".31.UNF0");
    }
    if (options.outSingleStoragesDaily) {
        G_daily31RiverStor.write(outputDirectory + "/G_RIVER_STORAGE_km3_" + std::to_string(year) + "_" + std::to_string(month+1) + ".31.UNF0");
    }
    if (options.outSingleStoragesDaily) {
        G_daily31ResStor.write(outputDirectory + "/G_RESERVOIR_STORAGE_km3_" + std::to_string(year) + "_" + std::to_string(month+1) + ".31.UNF0");
    }
    if (options.outSingleStoragesDaily || options.outGWStorageDaily) {
        G_daily31GwStor.write(outputDirectory + "/G_GROUND_WATER_STORAGE_km3_" + std::to_string(year) + "_" + std::to_string(month+1) + ".31.UNF0");
    }
    if (options.outTotalWaterInStoragesDaily_km3) {
        G_daily31TotalWaterInStorages_km3.write(outputDirectory + "/G_TOTAL_STORAGES_km3_" + std::to_string(year) + "_" + std::to_string(month+1) + ".31.UNF0");
    }
    if (options.outTotalWaterInStoragesDaily_mm) {
        G_daily31TotalWaterInStorages_mm.write(outputDirectory + "/G_TOTAL_STORAGES_mm_" + std::to_string(year) + "_" + std::to_string(month+1) + ".31.UNF0");
    }
    if (options.outGwrSwbDaily) {
        // groundwater recharge below surface waterbodies (with respect to continental area) [mm/month]
        G_daily31GwrSwb.write(outputDirectory + "/G_GWR_SURFACE_WATER_BODIES_mm_" + std::to_string(year) + "_" + std::to_string(month+1) + ".31.UNF0");
    }
    if (options.outFswbDaily) {
        // fraction of surface water bodies
        G_daily31Fswb.write(outputDirectory + "/G_FRACTION_SURFACE_WATER_BODIES_" + std::to_string(year) + "_" + std::to_string(month+1) + ".31.UNF0");
    }
    if (options.outLandAreaFracDaily) {
        // fraction of surface water bodies
        G_daily31LandAreaFrac.write(outputDirectory + "/G_LAND_AREA_FRACTION_" + std::to_string(year) + "_" + std::to_string(month+1) + ".31.UNF0");
    }
}

// daily output option 365 start
void routingClass::writeDaily365Grids(const std::string outputDirectory, const int year) {
    if (5 != options.grid_store && 6 != options.grid_store)
        return;

    if (options.outCellRunoffDaily) {
        G_daily365CellRunoff.write(outputDirectory + "/G_CELL_RUNOFF_" + std::to_string(year) + ".365.UNF0");
    }
    if (options.outCellRunoffDaily) {
        G_daily365CellRunoff_mm.write(outputDirectory + "/G_CELL_RUNOFF_mm_" + std::to_string(year) + ".365.UNF0");
    }
    if (options.outLandAreaFracDaily) {
        G_daily365LandAreaFrac.write(outputDirectory + "/G_LAND_AREA_FRACTION_" + std::to_string(year) + ".365.UNF0");
    }
    if (options.outGWRunoffDaily) {
        G_daily365GwRunoff.write(outputDirectory + "/G_GW_RUNOFF_" + std::to_string(year) + ".365.UNF0");
    }
    if (options.outGwrunSurfrunDaily) {
        G_daily365GwrunSurfrun.write(outputDirectory + "/G_GWRUN_SURFRUN_" + std::to_string(year) + ".365.UNF0");
    }
    if (options.outCellAETWCaDaily) {
        G_daily365CellAETWCa.write(outputDirectory + "/G_CELLAET_CONSUSE_km3_" + std::to_string(year) + ".365.UNF0");
    }
    if (options.outPrecDaily) {
        G_daily365ConsistentPrecip.write(outputDirectory + "/G_CONSISTENT_PRECIPITATION_km3_" + std::to_string(year) + ".365.UNF0");
    }
    if (options.outCellSurfaceDaily) {
        G_daily365CellSurfaceRunoff.write(outputDirectory + "/G_CELL_SURFACE_RUNOFF_" + std::to_string(year) + ".365.UNF0");
    }
    if (options.outRiverAvailDaily) {
        G_daily365RiverAvail.write(outputDirectory + "/G_RIVER_AVAIL_" + std::to_string(year) + ".365.UNF0");
    }
    if (options.outRiverVeloDaily) {
        G_daily365Velocity.write(outputDirectory + "/G_RIVER_VELOCITY_" + std::to_string(year) + ".365.UNF0");
    }
    if (options.outCellAETDaily) {
        G_daily365CellAET.write(outputDirectory + "/G_CELL_AET_" + std::to_string(year) + ".365.UNF0");
    }
    if (options.outSurfStorDaily) {
        G_daily365SurfStor.write(outputDirectory + "/G_SURFACE_WATER_STORAGE_km3_" + std::to_string(year) + ".365.UNF0");
    }
    if (options.outSingleStoragesDaily) {
        G_daily365LocLakeStor.write(outputDirectory + "/G_LOC_LAKE_STORAGE_km3_" + std::to_string(year) + ".365.UNF0");
    }
    if (options.outSingleStoragesDaily) {
        G_daily365LocWetlStor.write(outputDirectory + "/G_LOC_WETL_STORAGE_km3_" + std::to_string(year) + ".365.UNF0");
    }
    if (options.outSingleStoragesDaily) {
        G_daily365GloLakeStor.write(outputDirectory + "/G_GLO_LAKE_STORAGE_km3_" + std::to_string(year) + ".365.UNF0");
    }
    if (options.outSingleStoragesDaily) {
        G_daily365GloWetlStor.write(outputDirectory + "/G_GLO_WETL_STORAGE_km3_" + std::to_string(year) + ".365.UNF0");
    }
    if (options.outSingleStoragesDaily) {
        G_daily365RiverStor.write(outputDirectory + "/G_RIVER_STORAGE_km3_" + std::to_string(year) + ".365.UNF0");
    }
    if (options.outSingleStoragesDaily) {
        G_daily365ResStor.write(outputDirectory + "/G_RESERVOIR_STORAGE_km3_" + std::to_string(year) + ".365.UNF0");
    }
    if (options.outSingleStoragesDaily || options.outGWStorageDaily) {
        G_daily365GwStor.write(outputDirectory + "/G_GROUND_WATER_STORAGE_km3_" + std::to_string(year) + ".365.UNF0");
    }
    if ((options.outGlacierStorageDaily) || (options.outSingleStoragesDaily)) {
        G_daily365GlacierStorage.write(outputDirectory + "/G_GLACIER_STORAGE_km3_" + std::to_string(year) + ".365.UNF0");
    }
    if (options.outGlacierStorageDaily_mm) {
        G_daily365GlacierStorage_mm.write(outputDirectory + "/G_GLACIER_STORAGE_mm_" + std::to_string(year) + ".365.UNF0");
    }
    if (options.outTotalWaterInStoragesDaily_km3) {
        G_daily365TotalWaterInStorages_km3.write(outputDirectory + "/G_TOTAL_STORAGES_km3_" + std::to_string(year) + ".365.UNF0");
    }
    if (options.outTotalWaterInStoragesDaily_mm) {
        G_daily365TotalWaterInStorages_mm.write(outputDirectory + "/G_TOTAL_STORAGES_mm_" + std::to_string(year) + ".365.UNF0");
    }
    if (options.outGwrSwbDaily) {
        // groundwater recharge below surface waterbodies (with respect to continental area) [mm/month]
        G_daily365GwrSwb.write(outputDirectory + "/G_GWR_SURFACE_WATER_BODIES_mm_" + std::to_string(year) + ".365.UNF0");
    }
    if (options.outFswbDaily) {
        // fraction of surface water bodies
        G_daily365Fswb.write(outputDirectory + "/G_FRACTION_SURFACE_WATER_BODIES_" + std::to_string(year) + ".365.UNF0");
    }
}

void routingClass::writeStartendGrids(const std::string outputDirectory, const int year) {
    G_startendTotalWaterInStorages_km3.write(outputDirectory + "/G_TOTAL_STORAGES_STARTEND_km3_" + std::to_string(year) + ".2.UNF0");
}

void routingClass::writeTotalCellRunoffGrid(const std::string outputDirectory, int year) {
    // used for calculation of correction factor during calibration

    if (options.outPotCellRunoffAnnual) {
        G_potCellRunoff.write(outputDirectory + "/G_POT_CELL_RUNOFF_" + std::to_string(year) + ".UNF0");
    }
}

// Reservoir operation start years
void routingClass::writeAnnualFreqGrids(const std::string outputDirectory, int year) {
    // yearly output for GFREQ, GFREQW and its changes (unit: percent of grid cell area)

    if (options.outLandWaterFractions) {

        G_fwaterfreq.write(outputDirectory + "/GFREQW_" + std::to_string(year) + ".UNF0");

        G_fwaterfreq_change.write(outputDirectory + "/GFREQW_change_" + std::to_string(year) + ".UNF0");

        G_landfreq.write(outputDirectory + "/GFREQ_" + std::to_string(year) + ".UNF0");

        G_landfreq_change.write(outputDirectory + "/GFREQ_change_" + std::to_string(year) + ".UNF0");

        G_landAreaFrac.write(outputDirectory + "/GLAND_AREA_FRAC_" + std::to_string(year) + ".UNF0");

        G_landAreaFrac_change.write(outputDirectory + "/GLAND_AREA_FRAC_change_" + std::to_string(year) + ".UNF0");
    }
}

void createAsciiFile(const std::string filename, ofstream &streamName) {
    streamName.open(filename);
    if (!streamName) {
        cerr << "Can not open file " << filename << " for writing" << endl;
        exit(-1);
    }
    streamName << "# " << getTimeString();
}

void routingClass::createDailyFiles(const std::string outputDir) {

//    extern cbasinClass cbasin;

    if (options.outStationDischargeDaily) {
        createAsciiFile(outputDir + "/STATION_DISCHARGE_DAILY.OUT", dailyDischargeFile);
        dailyDischargeFile << "# Year\tDay\t";
        for (int n = 0; n <= cbasin.numberOfCalibBasins - 1; n++)
            dailyDischargeFile << cbasin.name[n] << '\t';
        dailyDischargeFile << endl;
    }
    if (options.outStationVelocityDaily) {
        createAsciiFile(outputDir + "/STATION_VELOCITY_DAILY.OUT", dailyVelocityFile);
        dailyVelocityFile << "# Year\tDay\t";
        for (int n = 0; n <= cbasin.numberOfCalibBasins - 1; n++)
            dailyVelocityFile << cbasin.name[n] << '\t';
        dailyVelocityFile << endl;
    }
    if (options.outLocLakeStorageDaily) {
        createAsciiFile(outputDir + "/LOC_LAKE_STORAGE.OUT", locLakeFile);
        for (int n = 0; n <= cbasin.numberOfCalibBasins - 1; n++)
            locLakeFile << cbasin.name[n] << '\t';
        locLakeFile << endl;
    }
    if (options.outLocWetlStorageDaily) {
        createAsciiFile(outputDir + "/LOC_WETL_STORAGE.OUT", locWetlFile);
        for (int n = 0; n <= cbasin.numberOfCalibBasins - 1; n++)
            locWetlFile << cbasin.name[n] << '\t';
        locWetlFile << endl;
    }
    if (options.outGloLakeStorageDaily) {
        createAsciiFile(outputDir + "/GLO_LAKE_STORAGE.OUT", gloLakeFile);
        for (int n = 0; n <= cbasin.numberOfCalibBasins - 1; n++)
            gloLakeFile << cbasin.name[n] << '\t';
        gloLakeFile << endl;
    }

    if (options.outResStorageDaily && (options.resOpt == 1)) {
        createAsciiFile(outputDir + "/RESERVOIR_STORAGE.OUT", gloResFile);
        for (int n = 0; n <= cbasin.numberOfCalibBasins - 1; n++)
            gloResFile << cbasin.name[n] << '\t';
        gloResFile << endl;

    }
    if (options.outGloWetlStorage) {
        createAsciiFile(outputDir + "/GLO_WETL_STORAGE.OUT", gloWetlFile);
        for (int n = 0; n <= cbasin.numberOfCalibBasins - 1; n++)
            gloWetlFile << cbasin.name[n] << '\t';
        gloWetlFile << endl;
    }
}

void routingClass::closeDailyFiles() {
    if (options.outStationDischargeDaily) dailyDischargeFile.close();
    if (options.outStationVelocityDaily) dailyVelocityFile.close();
    if (options.outLocLakeStorageDaily) locLakeFile.close();
    if (options.outLocWetlStorageDaily) locWetlFile.close();
    if (options.outGloLakeStorageDaily) gloLakeFile.close();
    if (options.outResStorageDaily && (options.resOpt == 1)) gloResFile.close();
    if (options.outGloWetlStorageDaily) gloWetlFile.close();
}

void routingClass::appendDailyVelocity(const short actualYear, const short remainingInitYears) {
    // write daily river velocity information into file
    if (options.outStationVelocityDaily) {
        for (short d = 0; d <= 364; d++) {
            dailyVelocityFile << setw(4) << actualYear << '\t' << setw(3) << d + 1 << '\t';
            for (short b = 0; b < nSpecBasins; ++b)
                dailyVelocityFile << setiosflags(ios::fixed | ios::showpoint)
                                  << setprecision(6) << setw(13)
                                  << dailyRiverVelocity(b,d) << '\t';
            dailyVelocityFile << endl;
        }
    }
}


void routingClass::appendDailyDischarge(const short actualYear, const short remainingInitYears) {
    // write daily river discharge information into file
    if (options.outStationDischargeDaily) {
        for (short d = 0; d <= 364; d++) {
            dailyDischargeFile << setw(4) << actualYear << '\t' << setw(3) << d + 1 << '\t';
            for (short b = 0; b < nSpecBasins; ++b)
                dailyDischargeFile << setiosflags(ios::fixed | ios::showpoint)
                                   << setprecision(6) << setw(13)
                                   << dailyRiverDischarge(b,d) << '\t';
            dailyDischargeFile << endl;
        }
    }
}

void routingClass::appendDailyLakeWetl(const short actualYear, const short remainingInitYears) {
    if (options.outLocLakeStorageDaily) {
        // store daily values of lake and wetland water content
        for (short d = 0; d <= 364; d++) {
            locLakeFile << setw(4) << actualYear << '\t' << setw(3) << d + 1 << '\t';
            for (short b = 0; b < nSpecBasins; ++b)
                locLakeFile << setiosflags(ios::fixed | ios::showpoint)
                            << setprecision(6) << setw(13)
                            << dailyLocLakeStorage(b,d) << '\t';
            locLakeFile << endl;
        }
    }
    if (options.outLocWetlStorageDaily) {
        for (short d = 0; d <= 364; d++) {
            locWetlFile << setw(4) << actualYear << '\t' << setw(3) << d + 1 << '\t';
            for (short b = 0; b < nSpecBasins; ++b)
                locWetlFile << setiosflags(ios::fixed | ios::showpoint)
                            << setprecision(6) << setw(13)
                            << dailyLocWetlStorage(b,d) << '\t';
            locWetlFile << endl;
        }
    }
    if (options.outGloLakeStorageDaily) {
        for (short d = 0; d <= 364; d++) {
            gloLakeFile << setw(4) << actualYear << '\t' << setw(3) << d + 1 << '\t';
            for (short b = 0; b < nSpecBasins; ++b)
                gloLakeFile << setiosflags(ios::fixed | ios::showpoint)
                            << setprecision(6) << setw(13)
                            << dailyGloLakeStorage(b,d) << '\t';
            gloLakeFile << endl;
        }
    }

    if (options.outResStorageDaily && (options.resOpt == 1)) {
        for (short d = 0; d <= 364; d++) {
            gloResFile << setw(4) << actualYear << '\t' << setw(3) << d + 1 << '\t';
            for (short b = 0; b < nSpecBasins; ++b)
                gloResFile << setiosflags(ios::fixed | ios::showpoint)
                           << setprecision(6) << setw(13) << dailyResStorage(b,d)
                           << '\t';
            gloResFile << endl;
        }
    }

    if (options.outGloWetlStorageDaily) {
        for (short d = 0; d <= 364; d++) {
            gloWetlFile << setw(4) << actualYear << '\t' << setw(3) << d + 1 << '\t';
            for (short b = 0; b < nSpecBasins; ++b)
                gloWetlFile << setiosflags(ios::fixed | ios::showpoint)
                            << setprecision(6) << setw(13)
                            << dailyGloWetlStorage(b,d) << '\t';
            gloWetlFile << endl;
        }
    }
}

void routingClass::createMonthlyDischargeFile(const string outputDir) {
    char filename[250];
    if (options.outStationDischargeMonthly) {

        sprintf(filename, "%s/STATION_DISCHARGE_MONTHLY.OUT", outputDir.c_str());
        createAsciiFile(filename, monthlyDischargeFile);
        monthlyDischargeFile
                       << "Stationname\t"; //HMS 2017-03-16: changed from "# Station name" to one single word in order to easy handling of output file e.g. in R-scripts.
        for (int n = 0; n <= cbasin.numberOfCalibBasins - 1; n++)
            monthlyDischargeFile << cbasin.name[n] << '\t';
        monthlyDischargeFile << endl;
    }

    if (options.outResvoirMonthly && // output option for this 4 files
        options.resOpt == 1 && // New reservoir algorithm (N. Hanasaki) will be used (*)
        ((2 == options.grid_store) || (4 == options.grid_store) ||
         (6 == options.grid_store)) // Store grid files, including monthly values
                   ) {
        // create reservoir storage file (ascii)
        sprintf(filename, "%s/RESERVOIR_STORAGE_MONTHLY.OUT", outputDir.c_str());
        createAsciiFile(filename, monthlyReservoirStorageFile);
        monthlyReservoirStorageFile << "# Year" << '\t' << "Month" << '\t';
        for (int n = 0; n < ng; n++) {
            if (G_reservoir_area[n] > 0.0) {
                monthlyReservoirStorageFile << n + 1 << '\t';
            }
        }
        monthlyReservoirStorageFile << endl;

        // create reservoir storage ratio file (ascii)
        sprintf(filename, "%s/RESERVOIR_STORAGE_RATIO_MONTHLY.OUT", outputDir.c_str());
        createAsciiFile(filename, monthlyReservoirStorageRatioFile);
        monthlyReservoirStorageRatioFile << "# Year" << '\t' << "Month" << '\t';
        for (int n = 0; n < ng; n++) {
            if (G_reservoir_area[n] > 0.0) {
                monthlyReservoirStorageRatioFile << n + 1 << '\t';
            }
        }
        monthlyReservoirStorageRatioFile << endl;

        // create reservoir inflow file (ascii)
        sprintf(filename, "%s/RESERVOIR_INFLOW_MONTHLY.OUT", outputDir.c_str());
        createAsciiFile(filename, monthlyReservoirInflowFile);
        monthlyReservoirInflowFile << "# Year" << '\t' << "Month" << '\t';
        for (int n = 0; n < ng; n++) {
            if (G_reservoir_area[n] > 0.0) {
                monthlyReservoirInflowFile << n + 1 << '\t';
            }
        }
        monthlyReservoirInflowFile << endl;

        // create reservoir outflow file (ascii)
        sprintf(filename, "%s/RESERVOIR_OUTFLOW_MONTHLY.OUT", outputDir.c_str());
        createAsciiFile(filename, monthlyReservoirOutflowFile);
        monthlyReservoirOutflowFile << "# Year" << '\t' << "Month" << '\t';
        for (int n = 0; n < ng; n++) {
            if (G_reservoir_area[n] > 0.0) {
                monthlyReservoirOutflowFile << n + 1 << '\t';
            }
        }
        monthlyReservoirOutflowFile << endl;
    }
}

void routingClass::closeMonthlyDischargeFile() {
    if (options.outStationDischargeMonthly)
        monthlyDischargeFile.close();

    if (options.outResvoirMonthly && // output option for this 4 files
        options.resOpt == 1 && // reservoir algorithm (N. Hanasaki) will be used (*)
        ((2 == options.grid_store) || (4 == options.grid_store) ||
         (6 == options.grid_store)) // Store grid files, including monthly values
                   ) {
        monthlyReservoirStorageFile.close();
        monthlyReservoirStorageRatioFile.close();
        monthlyReservoirInflowFile.close();
        monthlyReservoirOutflowFile.close();
    }
}

void routingClass::appendMonthlyDischarge(const short actualYear, const short remainingInitYears) {
    char number_of_days_in_month[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    if (options.outStationDischargeMonthly) {
        short day = 0;

        for (short i = 0; i <= 11; i++) {
            if (remainingInitYears != 0)
                monthlyDischargeFile << "# ";
            monthlyDischargeFile << setw(4) << (actualYear - 1901) * 12 + i + 1 << '\t';
            for (short b = 0; b < nSpecBasins; ++b) {
                double monthlyDischarge = 0;
                for (short d = 0; d <= number_of_days_in_month[i] - 1; d++) {
                    monthlyDischarge += dailyRiverDischarge(b,day + d);
                }
                monthlyDischargeFile
                               << setiosflags(ios::fixed | ios::showpoint)
                               << setprecision(6) << setw(13) << monthlyDischarge << '\t';
            }
            monthlyDischargeFile << endl;
            day += number_of_days_in_month[i];
        }
    }

    if (options.outResvoirMonthly && // output option for this 4 files
        options.resOpt == 1 && // New reservoir algorithm (N. Hanasaki) will be used
        ((2 == options.grid_store) || (4 == options.grid_store) ||
         (6 == options.grid_store)) // Store grid files, including monthly values
                   ) {
        for (short month = 0; month < 12; month++) {
            if (remainingInitYears != 0) {
                monthlyReservoirStorageFile << "# ";
                monthlyReservoirStorageRatioFile << "# ";
                monthlyReservoirInflowFile << "# ";
                monthlyReservoirOutflowFile << "# ";
            }
            monthlyReservoirStorageFile << setw(4) << actualYear << '\t';
            monthlyReservoirStorageFile << setw(4) << month + 1 << '\t';

            monthlyReservoirStorageRatioFile << setw(4) << actualYear << '\t';
            monthlyReservoirStorageRatioFile << setw(4) << month + 1 << '\t';

            monthlyReservoirInflowFile << setw(4) << actualYear << '\t';
            monthlyReservoirInflowFile << setw(4) << month + 1 << '\t';

            monthlyReservoirOutflowFile << setw(4) << actualYear << '\t';
            monthlyReservoirOutflowFile << setw(4) << month + 1 << '\t';

            for (int n = 0; n < ng; n++) {

                if (G_reservoir_area[n] > 0.0) {
                    if (options.grid_store_TypeForStorages == 0)
                        G_monthlyResStorageRatio(n,month) = G_monthlyResStorage(n,month) / G_stor_cap[n];
                    else
                        G_monthlyResStorageRatio(n,month) =
                                       G_monthlyResStorage(n,month) / number_of_days_in_month[month] / G_stor_cap[n];

                    // (in routing, timestep flow values [m3/timestep] are summed up to calculate total monthly values here [m3/month]!)
                    // [m3/s]
                    G_monthlyResInflow(n,month) = G_monthlyResInflow(n,month) / (number_of_days_in_month[month] * 86400.);
                    G_monthlyResOutflow(n,month) = G_monthlyResOutflow(n,month) / (number_of_days_in_month[month] * 86400.);

                    double tmp_G_monthlyResStorage;
                    if (options.grid_store_TypeForStorages == 0) tmp_G_monthlyResStorage = G_monthlyResStorage(n,month);
                    else tmp_G_monthlyResStorage = G_monthlyResStorage(n,month) / number_of_days_in_month[month];

                    monthlyReservoirStorageFile << setiosflags(ios::fixed | ios::showpoint)
                                                << setprecision(6) << setw(13) << tmp_G_monthlyResStorage << '\t';
                    monthlyReservoirStorageRatioFile << setiosflags(ios::fixed | ios::showpoint)
                                                     << setprecision(6) << setw(13) << G_monthlyResStorageRatio(n,month)
                                                     << '\t';
                    monthlyReservoirInflowFile << setiosflags(ios::fixed | ios::showpoint)
                                               << setprecision(6) << setw(13) << (G_monthlyResInflow(n,month)) << '\t';
                    monthlyReservoirOutflowFile << setiosflags(ios::fixed | ios::showpoint)
                                                << setprecision(6) << setw(13) << (G_monthlyResOutflow(n,month)) << '\t';
                } else
                    G_monthlyResStorageRatio(n,month) = -99.99;

            }
            monthlyReservoirStorageFile << endl;
            monthlyReservoirStorageRatioFile << endl;
            monthlyReservoirInflowFile << endl;
            monthlyReservoirOutflowFile << endl;
        }
    }
}

void routingClass::createAnnualDischargeFile(const string outputDir) {
    if (options.outStationDischargeAnnual) {
        char filename[250];

        sprintf(filename, "%s/STATION_DISCHARGE_ANNUAL.OUT", outputDir.c_str());
        createAsciiFile(filename, annualDischargeFile);
        annualDischargeFile
                       << "Stationname\t"; // single word in order to easy handling of output file e.g. in R-scripts.
        for (int n = 0; n <= cbasin.numberOfCalibBasins - 1; n++)
            annualDischargeFile << cbasin.name[n] << '\t';
        annualDischargeFile << endl;
    }
}

void routingClass::closeAnnualDischargeFile() {
    if (options.outStationDischargeAnnual) annualDischargeFile.close();
}

void routingClass::appendAnnualDischarge(const short actualYear, const short remainingInitYears) {
    if (options.outStationDischargeAnnual) {
        if (remainingInitYears != 0)
            annualDischargeFile << "# ";
        annualDischargeFile << setw(4) << actualYear << '\t';

        for (short b = 0; b < nSpecBasins; ++b) {
            double annualDischarge = 0;
            for (short i = 0; i <= 364; i++) {
                annualDischarge += dailyRiverDischarge(b,i);
            }
            annualDischargeFile << setiosflags(ios::fixed | ios::showpoint)
                                << setprecision(6) << setw(13) << annualDischarge << '\t';
        }
        annualDischargeFile << endl;
    }
}

void routingClass::appendCalibInfoToAnnualDischargeFile(const double gamma) {
    if (options.outStationDischargeAnnual) {
        annualDischargeFile << "# Next step of calibration." << endl;
        annualDischargeFile << "# Gamma = " << gamma << endl;
    }
}

double routingClass::getAnnualDischarge(const int stationNumber) {
    double annualDischarge = 0;

    for (short i = 0; i <= 364; i++) {
        annualDischarge += dailyRiverDischarge(stationNumber - 1,i);
    }
    return annualDischarge;
}

double routingClass::getAnnualUpstreamInflow(const int stationNumber) {
    // calculates annual inflow from upstream basins

//    extern upstreamStationClass directUpstSt;
    int upSt;
    short i;
    double stationAnnualDischarge = 0, totalAnnualDischarge = 0;

    for (int n = 1; n <= directUpstSt.getNumberOfUpstreamStations(stationNumber); n++) {
        stationAnnualDischarge = 0;
        upSt = directUpstSt.getUpstreamStation(stationNumber, n);
        for (i = 0; i <= 364; i++) {
            stationAnnualDischarge += dailyRiverDischarge(upSt - 1,i);
        }
        totalAnnualDischarge += stationAnnualDischarge;
    }
    return totalAnnualDischarge;
}

double routingClass::getAnnualSatisfWaterUse(const int stationNumber) {

    double satisfiedBasinUse = 0;

    for (int n = 0; n <= ng - 1; n++) {
        if (stationNumber == G_sbasin[n]) {
            // A distinction between satisfied use and actual use is not possible in 22b.
            satisfiedBasinUse += G_actualUse[n];
        }
    }

    return satisfiedBasinUse;
}

double routingClass::getRiverVelocity(const double RiverSlope, const double G_riverBottomWidth,
                                      const double Roughness, double G_riverInflow, const int n, calibParamClass &calParam) // OE: added calParam
            {
    double river_velocity;

    double riverDepth;
    double wettedPerimeter;
    double crossSectionalArea;
    double hydraulicRad;
    double incoming_discharge;

    incoming_discharge = (G_riverInflow * 1000. * 1000. * 1000.) / (60. * 60. * 24.); //[km3/day]  --> [m3/sec]

    riverDepth = 0.349 * pow(incoming_discharge, 0.341); //[m]

    // trapezoidal channel shape with channel sides 2:1 run to rise ratio
    crossSectionalArea = riverDepth * (2.0 * riverDepth + G_riverBottomWidth); // trapezoidal river shape
    wettedPerimeter = G_riverBottomWidth + 2.0 * riverDepth * sqrt(5.0); // sqrt(1+2^2)

    hydraulicRad = crossSectionalArea / wettedPerimeter;

    // calculate riverVelocity
    // Cell-specific calibration parameters - Apply multiplier
    river_velocity = 1. / (calParam.getValue(M_RIVRGH_C, n) * Roughness) * pow(hydraulicRad, (2. / 3.)) *
                     pow(RiverSlope, 0.5); //[m/sec]
    river_velocity = river_velocity * 86.4; //m/sec -->km/day (60*60*24)/1000

    // return value in [km/day]; if value is below 1cm/day then set limit

    if (river_velocity < 0.00001)
        return 0.00001;
    else
        return river_velocity;
}

double routingClass::getNewRiverVelocity(const double RiverSlope, const double G_riverBottomWidth,
                                         const double Roughness, double G_riverStorage, double G_riverLength, const int n,calibParamClass &calParam)// OE: added calParam
            {
    double river_velocity;

    double riverDepth;
    double wettedPerimeter;
    double crossSectionalArea;
    double hydraulicRad;
    double incoming_discharge;

    // trapezoidal channel shape with channel sides 2:1 run to rise ratio
    crossSectionalArea = (G_riverStorage / G_riverLength) * 1000000.;	// km3/km *10^6 -> m2

    // G_riverBottomWidth in m
    riverDepth = - G_riverBottomWidth/4.0 + sqrt(G_riverBottomWidth * G_riverBottomWidth/16.0 + 0.5 * crossSectionalArea);

    wettedPerimeter = G_riverBottomWidth + 2.0 * riverDepth * sqrt(5.0); // sqrt(1+2^2)

    hydraulicRad = crossSectionalArea / wettedPerimeter;

    // calculate riverVelocity
    // Cell-specific calibration parameters - Apply multiplier
    river_velocity = 1. / (calParam.getValue(M_RIVRGH_C,n) * Roughness) * pow(hydraulicRad, (2. / 3.)) * pow(RiverSlope, 0.5); //[m/sec]
    river_velocity = river_velocity * 86.4; //m/sec -->km/day (60*60*24)/1000

    // return value in [km/day]; if value is below 1cm/day then set limit

    if (river_velocity < 0.00001)
        return 0.00001;
    else
        return river_velocity;
}

void routingClass::getAnnualDischargeVector(vector<float>&annualDischargeVector) {
    for (short n = 0; n < nSpecBasins; n++) {
        float annualDischarge = 0;
        for (short i = 0; i <= 364; i++) {
            annualDischarge += dailyRiverDischarge(n,i);
        }
        annualDischarge *= 1000000000 / (double) (365 * 24 * 60 * 60); // km3/a -> m3/s
        annualDischargeVector.push_back(annualDischarge);
    }
}

void routingClass::getMonthlyDischargeVector(short month, vector<float> &monthlyDischargeVector) {
    short day = 0;
    char number_of_days_in_month[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    for (short i = 0; i <= 11; i++) {
        if (i == month) {
            for (short b = 0; b < nSpecBasins; ++b) {
                float monthlyDischarge = 0;
                for (short d = 0; d < number_of_days_in_month[i]; d++) {
                    monthlyDischarge += dailyRiverDischarge(b,day + d);
                }
                monthlyDischarge *= 1000000000 / (float) (24 * 60 * 60 * number_of_days_in_month[i]);
                monthlyDischargeVector.push_back(monthlyDischarge);
            }
            return;
        }
        day += number_of_days_in_month[i];
    }
}

void routingClass::getMinMonthlyDischargeVector(vector<float> &monthlyDischargeVector) {
    char number_of_days_in_month[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    for (short b = 0; b < nSpecBasins; ++b) {
        float minValue = 99E29;
        short day = 0;
        for (short i = 0; i <= 11; i++) {
            float monthlyDischarge = 0;
            for (short d = 0; d <= number_of_days_in_month[i] - 1; d++) {
                monthlyDischarge += dailyRiverDischarge(b,day + d);
            }
            monthlyDischarge *= 1000000000 / (float) (24 * 60 * 60 * number_of_days_in_month[i]);
            if (monthlyDischarge < minValue)
                minValue = monthlyDischarge;
            day += number_of_days_in_month[i];
        }
        monthlyDischargeVector.push_back(minValue);
    }
}

void routingClass::getMinToMaxDischargeVector(vector<float> &monthlyDischargeVector) {
    char number_of_days_in_month[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    for (short b = 0; b < nSpecBasins; ++b) {
        float minValue = 99E29;
        float maxValue = -99E29;
        short day = 0;
        for (short i = 0; i <= 11; i++) {
            float monthlyDischarge = 0;
            for (short d = 0; d <= number_of_days_in_month[i] - 1; d++) {
                monthlyDischarge += dailyRiverDischarge(b,day + d);
            }
            monthlyDischarge *= 1000000000 / (float) (24 * 60 * 60 * number_of_days_in_month[i]);
            if (monthlyDischarge < minValue)
                minValue = monthlyDischarge;
            if (monthlyDischarge > maxValue)
                maxValue = monthlyDischarge;
            day += number_of_days_in_month[i];
        }
        monthlyDischargeVector.push_back(minValue / maxValue);
    }
}

// Implementation of class method calcNextDay_M(...)
void routingClass::calcNextDay_M(const short month) {
    /*
    static short month = 11;
    static short old_day;
    short old_month;

    const short int last_day_per_month[12]
                   = {31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365};
    old_month = month;

    // does the day beint to the next year already ?
    if (new_day == 1) {
        month = 0;
        old_day = 0;
    }
    // does the day beint to the next month ?
    while (new_day > last_day_per_month[month])
        month++;

    if (new_day != old_day + 1) {
        cerr << "Day counter mismatch!\n";
        exit(-1);
    }
*/
    for (int n = 0; n < ng ; n++) {
        G_dailydailyNUs[n] = G_dailyNUs(n,month);
        G_dailydailyNUg[n] = G_dailyNUg(n,month);

        if (1 == options.aggrNUsGloLakResOpt) {
            G_dailydailyNUsAggregated[n] = G_dailyNUsAggregated(n,month);
        }
    }
    //old_day = new_day;
}

routingClass::~routingClass() {

    monthlyDischargeFile.close();
    annualDischargeFile.close();
    if (options.outStationDischargeDaily) dailyDischargeFile.close();
    if (options.outStationVelocityDaily) dailyVelocityFile.close();
    if (options.outLocLakeStorageDaily) locLakeFile.close();
    if (options.outLocWetlStorageDaily) locWetlFile.close();
    if (options.outGloLakeStorageDaily) gloLakeFile.close();
    if (options.outResStorageDaily && (options.resOpt == 1))
        gloResFile.close();
    if (options.outGloWetlStorageDaily) gloWetlFile.close();
}

