/*
 * GrayFox
 * V1.0, June 2014
 * 
 * C. Ian Short
 * Saint Mary's University
 * Department of Astronomy and Physics
 * Halifax, NS, Canada
 *
 * 1D, static, plane-parallel, LTE, gray stellar atmospheric model
 * core + wing approximation to Voigt spectral line profile
 *
 * Suitable for pedagogical purposes only
 * 
 * Logic written in Java SE 8.0, JDK 1.8
 * GUI written with JavaFX 8.0
 *
 * System requirements: Java run-time environment (JRE)
 *
 * Code provided "as is" - there is no formal support 
 *
 * Java default license applies:
 * End users may adapt, modify, develop, debug, and deploy at will
 * for academic and othe non-profit uses, but are asked to leave this
 * header text in place (although they may add to the header text).
 *
 */
package grayfox;
// Returns depth distribution of occupation numbers in lower level of b-b transition,
// and in ground states of neutral and singly ionized stages for reference

// Input parameters:
// lam0 - line centre wavelength in nm
// logNl - log_10 column density of absorbers in lower E-level, l (cm^-2)
// logFlu - log_10 oscillator strength (unitless)
// chiL - energy of lower atomic E-level of b-b transition in eV
// chiI - ground state ionization energy to niext higher stage in (ev)
//   - we are assuming this is the neutral stage
// Also needs atsmopheric structure information:
// numDeps
// tauRos structure
// temp structure 
// rho structure
public class LevelPops {

    //public static double[][] levelPops(double lam0In, double logNlIn, double[][] Ne, boolean ionized, double chiI, double chiL, double[][] linePoints, double[][] lineProf,
    //        int numDeps, double kappaScale, double[][] tauRos, double[][] temp, double[][] rho) {
    public static double[][] levelPops(double lam0In, double logNlIn, double[][] Ne, boolean ionized, double chiI1, double chiI2, double chiL, double gw1, double gw2, double gwL,
            int numDeps, double kappaScale, double[][] tauRos, double[][] temp, double[][] rho) {

        double c = Useful.c;
        double logC = Useful.logC();
        double k = Useful.k;
        double logK = Useful.logK();
        double logH = Useful.logH();
        double logEe = Useful.logEe();
        double logMe = Useful.logMe();

        double ln10 = Math.log(10.0);
        double logE = Math.log10(Math.E); // for debug output
        double log2pi = Math.log(2.0 * Math.PI);
        double log2 = Math.log(2.0);

        double logNl = logNlIn * ln10;  // Convert to base e
        double logKScale = Math.log(kappaScale);

        //Assume ground state statistical weight (or partition fn) of Stage III is 1.0;
        double logGw3 = 0.0;

        //For now:
        //double gw1 = 1.0;  //stat weight ground state of Stage I
        //double gw2 = 1.0;  //stat weight ground state of Stage II
        //double gwL = 1.0;  //stat weight lower E-level
        double gwU = 1.0;  //stat weight upper E-level 

        double logGw1 = Math.log(gw1);
        double logGw2 = Math.log(gw2);
        double logGwL = Math.log(gwL);
        double logGwU = Math.log(gwU);

        //System.out.println("chiL before: " + chiL);
        // If we need to subtract chiI from chiL, do so *before* converting to tiny numbers in ergs!
        if (ionized) {
            //System.out.println("ionized, doing chiL - chiI: " + ionized);
//            chiL = chiL - chiI;
            chiL = chiL - chiI1;
        }
        //System.out.println("chiL after: " + chiL);

//        chiI = chiI * Useful.eV;  // Convert lower E-level from eV to ergs
        chiI1 = chiI1 * Useful.eV;  // Convert lower E-level from eV to ergs
        chiI2 = chiI2 * Useful.eV;  // Convert lower E-level from eV to ergs

//        double boltzFacI = chiI / k; // Pre-factor for exponent of ionization Boltzmann factor
        double boltzFacI1 = chiI1 / k; // Pre-factor for exponent of ionization Boltzmann factor for ion stage I
        double boltzFacI2 = chiI2 / k; // Pre-factor for exponent of ionization Boltzmann factor for ion stage I
        //System.out.println("boltzFacI1 " + boltzFacI1 + " boltzFacI2 " + boltzFacI2 + " chiI1 " + chiI1 + " chiI2 " + chiI2);

        double logSahaFac = log2 + (3.0 / 2.0) * (log2pi + logMe + logK - 2.0 * logH);

        chiL = chiL * Useful.eV;  // Convert lower E-level from eV to ergs

        //Log of line-center wavelength in cm
        double logLam0 = Math.log(lam0In); // * 1.0e-7);

        // energy of b-b transition
        double logTransE = logH + logC - logLam0; //ergs
        // Energy of upper E-level of b-b transition
        double chiU = chiL + Math.exp(logTransE);  //ergs

        double boltzFacL = chiL / k; // Pre-factor for exponent of excitation Boltzmann factor
        double boltzFacU = chiU / k; // Pre-factor for exponent of excitation Boltzmann factor

        double boltzFacGround = 0.0 / k; //I know - its zero, but let's do it this way anyway'

        int refRhoIndx = TauPoint.tauPoint(numDeps, tauRos, 1.0);
        double refLogRho = rho[1][refRhoIndx];
        //System.out.println("LINEKAPPA: refRhoIndx, refRho " + refRhoIndx + " " + logE*refRho);

        // return a 2D 3 x numDeps array of logarithmic number densities
        // Row 0: neutral stage ground state population
        // Row 1: ionized stage ground state population
        // Row 2: level population of lower level of bb transition (could be in either stage I or II!) 
        // Row 3: level population of upper level of bb transition (could be in either stage I or II!) 
        double[][] logNums = new double[4][numDeps];

        double num, logNum, expFac, logSaha, saha, logIonFracI, logIonFracII, logNumI, logNumII;
        double saha21, logSaha21, saha32, logSaha32;
        logNumI = 0.0;
        logNumII = 0.0;
        double logNe;

        //Get electron densities from class State
        //No! double[][] mmwNe = State.mmwNeFn(numDeps, temp, kappaScale);
        for (int id = 0; id < numDeps; id++) {

            // reduce or enhance number density by over-all Rosseland opcity scale parameter
            logNum = logNl + logKScale;

            // scale numer density by relative depth variation of mass density
            logNum = logNum + rho[1][id] - refLogRho;

            //// reduce number density by temperature-dependent factor of Saha equation:
            // Normalize wrt to solar Teff so we don't have to evaluate all the other stuff
            //Row 1 of Ne is log_e Ne in cm^-3
            logNe = Ne[1][id];

            /*
             //
             // ********** Accounting for only TWO ionization stages (I & II):
             //
             // This assumes partition fns of unity - not any more:
             logSaha = logSahaFac - logNe - boltzFacI / temp[0][id] + (3.0) / (2.0) * temp[1][id]; // log(RHS) of standard Saha equation
             saha = Math.exp(logSaha);   //RHS of standard Saha equation
             //System.out.println("logSahaFac, logNe, logSaha= " + logE*logSahaFac + " " + logE*logNe + " " + logE*logSaha);

             logIonFracII = logSaha - Math.log(1.0 + saha); // log ionization fraction in stage II
             logIonFracI = -1.0 * Math.log(1.0 + saha);     // log ionization fraction in stage I
             */
            //
            // ********** Accounting for THREE ionization stages (I, II, III):
            //
            logSaha21 = logSahaFac - logNe - boltzFacI1 / temp[0][id] + (3.0 * temp[1][id] / 2.0) + logGw2 - logGw1; // log(RHS) of standard Saha equation
            saha21 = Math.exp(logSaha21);   //RHS of standard Saha equation
            logSaha32 = logSahaFac - logNe - boltzFacI2 / temp[0][id] + (3.0 * temp[1][id] / 2.0) + logGw3 - logGw2; // log(RHS) of standard Saha equation
            saha32 = Math.exp(logSaha32);   //RHS of standard Saha equation
            //System.out.println("logSahaFac, logNe, logSaha= " + logE*logSahaFac + " " + logE*logNe + " " + logE*logSaha);

            logIonFracII = logSaha21 - Math.log(1.0 + saha21 + saha32 * saha21); // log ionization fraction in stage II
            logIonFracI = -1.0 * Math.log(1.0 + saha21 + saha32 * saha21);     // log ionization fraction in stage I
            //if (id == 36) {
            //    System.out.println("logSaha21 " + logE*logSaha21 + " logSaha32 " + logE*logSaha32);
            //    System.out.println("IonFracII " + Math.exp(logIonFracII) + " IonFracI " + Math.exp(logIonFracI) + " logNe " + logE*logNe);
            //}

            //System.out.println("LevelPops: id, ionFracI, ionFracII: " + id + " " + Math.exp(logIonFracI) + " " + Math.exp(logIonFracII) );
            if (ionized) {
                //System.out.println("LevPops: ionized branch taken, ionized =  " + ionized);

                logNums[0][id] = logNum + logIonFracI; // Ascribe entire neutral stage pop to its ground level
                logNumII = logNum + logIonFracII;
                logNums[1][id] = logNumII - boltzFacGround / temp[0][id] + logGw2; // ground level of ionized stage
                logNums[2][id] = logNumII - boltzFacL / temp[0][id] + logGwL; // lower level of b-b transition
                logNums[3][id] = logNumII - boltzFacU / temp[0][id] + logGwU; // upper level of b-b transition

            } else {
                //System.out.println("LevPops: neutral branch taken, ionized =  " + ionized);

                logNumI = logNum + logIonFracI;
                logNums[0][id] = logNumI - boltzFacGround / temp[0][id] + logGw1;  // ground level of neutral stage
                logNums[1][id] = logNum + logIonFracII; // Ascribe entire ionized stage pop to its ground level
                logNums[2][id] = logNumI - boltzFacL / temp[0][id] + logGwL; // lower level of b-b transition

                logNums[3][id] = logNumI - boltzFacU / temp[0][id] + logGwU; // upper level of b-b transition

            }

            // System.out.println("LevelPops: id, logNums[0][id], logNums[1][id], logNums[2][id], logNums[3][id]: " + id + " "
            //          + Math.exp(logNums[0][id]) + " "
            //         + Math.exp(logNums[1][id]) + " "
            //          + Math.exp(logNums[2][id]) + " "
            //        + Math.exp(logNums[3][id]));
            //System.out.println("LevelPops: id, logNums[0][id], logNums[1][id], logNums[2][id], logNums[3][id]: " + id + " "
            //       + logE * (logNums[0][id]) + " "
            //        + logE * (logNums[1][id]) + " "
            //        + logE * (logNums[2][id]) + " "
            //       + logE * (logNums[3][id]) );
            //System.out.println("LevelPops: id, logIonFracI, logIonFracII: " + id + " " + logE*logIonFracI + " " + logE*logIonFracII
            //        + "logNum, logNumI, logNums[0][id], logNums[1][id] "
            //        + logE*logNum + " " + logE*logNumI + " " + logE*logNums[0][id] + " " + logE*logNums[1][id]);
            //System.out.println("LevelPops: id, logIonFracI: " + id + " " + logE*logIonFracI
            //        + "logNums[0][id], boltzFacL/temp[0][id], logNums[2][id]: " 
            //        + logNums[0][id] + " " + boltzFacL/temp[0][id] + " " + logNums[2][id]);
        } //id loop

        return logNums;
    }
}
