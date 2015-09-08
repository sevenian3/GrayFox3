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

public class LineKappa {

// Assumes CRD, LTE, ???
// Input parameters:
// lam0 - line centre wavelength in nm
// logNl - log_10 column density of absorbers in lower E-level, l (cm^-2)
// logFlu - log_10 oscillator strength (unitless)
// chiL - energy of lower atomic E-level of b-b transition in eV
// chiI - ground state ionization energy to niext higher stage in (ev)
    //   
//     * PROBLEM: line kappaL values converted to mass extinction by division by rho() are 
// * not consistent with fake Kramer's Law based scaling of kappa_Ros with g.
    //*   Try leaving kappaLs as linear extinctions and converting the scaled kappa_Ros back to linear units
// * with solar rho() in LineTau2
    //
// Also needs atsmopheric structure information:
// numDeps
// tauRos structure
// temp structure 
// rho structure
//    public static double[][] lineKap(double lam0In, double logNlIn, double logFluIn, boolean ionized, double chiI, double chiL, double[][] linePoints, double[][] lineProf,
//            int numDeps, double kappaScale, double[][] tauRos, double[][] temp, double[][] rho) {
// Level population now computed in LevelPops.levelPops()
    public static double[][] lineKap(double lam0In, double[][] logNums, double logFluIn, double[][] linePoints, double[][] lineProf,
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

        double lam0 = lam0In; // * 1.0E-7; //nm to cm
        double logLam0 = Math.log(lam0);
        //double logNl = logNlIn * ln10;  // Convert to base e
        double logFlu = logFluIn * ln10; // Convert to base e
        double logKScale = Math.log10(kappaScale);

        //chiI = chiI * Useful.eV;  // Convert lower E-level from eV to ergs
        //double boltzFacI = chiI / k; // Pre-factor for exponent of excitation Boltzmann factor
        //double logSahaFac = log2 + (3.0/2.0) * ( log2pi + logMe + logK - 2.0*logH);
        //chiL = chiL * Useful.eV;  // Convert lower E-level from eV to ergs
        //double boltzFac = chiL / k; // Pre-factor for exponent of excitation Boltzmann factor
        int numPoints = linePoints[0].length;
        //System.out.println("LineKappa: numPoints: " + numPoints);

        double logPreFac;
        //This converts f_lu to a volume extinction coefficient per particle - Rutten, p. 23
        logPreFac = logFlu + Math.log(Math.PI) + 2.0 * logEe - logMe - logC;
        //System.out.println("LINEKAPPA: logPreFac " + logPreFac);

        //Assume wavelength, lambda, is constant throughout line profile for purpose
        // of computing the stimulated emission correction
        double logExpFac;
        logExpFac = logH + logC - logK - logLam0;
        //System.out.println("LINEKAPPA: logExpFac " + logExpFac);

        // int refRhoIndx = TauPoint.tauPoint(numDeps, tauRos, 1.0);
        // double refLogRho = rho[1][refRhoIndx];
        //System.out.println("LINEKAPPA: refRhoIndx, refRho " + refRhoIndx + " " + logE*refRho);
        // return a 2D numPoints x numDeps array of monochromatic *LINE* extinction line profiles
        double[][] logKappaL = new double[numPoints][numDeps];

        double num, logNum, logExpFac2, expFac, stimEm, logStimEm, logSaha, saha, logIonFrac;
        double logNe;

        for (int id = 0; id < numDeps; id++) {

            logExpFac2 = logExpFac - temp[1][id];
            expFac = -1.0 * Math.exp(logExpFac2);

            stimEm = 1.0 - Math.exp(expFac);
            logStimEm = Math.log(stimEm);

            // logNums is a 2D 3 x numDeps array of logarithmic number densities
            // Row 0: neutral stage ground state population
            // Row 1: ionized stage ground state population
            // Row 2: level population of lower level of bb transition (could be in either stage I or II!) 
            logNum = logNums[2][id];

            //if (id == refRhoIndx) {
            //    System.out.println("LINEKAPPA: logStimEm " + logE*logStimEm);
            //}
            for (int il = 0; il < numPoints; il++) {

                // From Radiative Transfer in Stellar Atmospheres (Rutten), p.31
                // This is a *volume* co-efficient ("alpha_lambda") in cm^-1:
                logKappaL[il][id] = logPreFac + logStimEm + logNum + Math.log(lineProf[il][id]);
                //if (id == 36) {
                //    System.out.println("il " + il + " logNum " + logE*logNum + " Math.log(lineProf[il][id]) " + logE*Math.log(lineProf[il][id]));
                ////    //System.out.println("logPreFac " + logPreFac + " logStimEm " + logStimEm);
                //}
                //System.out.println("LINEKAPPA: id, il " + id + " " + il + " logKappaL " + logE * logKappaL[il][id]);

                //Convert to mass co-efficient in g/cm^2:                
                // This direct approach won't work - is not consistent with fake Kramer's law scaling of Kapp_Ros with g instead of rho
                logKappaL[il][id] = logKappaL[il][id] - rho[1][id];

                // if (id == refRhoIndx) {
                //   System.out.println("LINEKAPPA: id, il " + id + " " + il + " logKappaL " + logE * logKappaL[il][id]);
                //    + " logPreFac " + logE*logPreFac + " logStimEm " + logE*logStimEm + " logNum " + logE*logNum 
                //   + " log(lineProf[1]) " + logE*Math.log(lineProf[1][il]) );
                //}
                //if (id == refRhoIndx-45) {
                //    System.out.println("LINEKAPPA: id, il " + id + " " + il + " logKappaL " + logE*logKappaL[il][id]
                //    + " logPreFac " + logE*logPreFac + " logStimEm " + logE*logStimEm + " logNum " + logE*logNum + " logRho " + logE*rho[1][id] 
                //    + " log(lineProf[1]) " + logE*Math.log(lineProf[1][il]) );
                //}
            } // il - lambda loop

        } // id - depth loop

        return logKappaL;

    }

    //Create total extinction throughout line profile:
    public static double[][] lineTotalKap(double[][] linePoints, double[][] logKappaL,
            int numDeps, double teff, double lam0, double kappaScale, double[][] kappa) {

        double logE = Math.log10(Math.E); // for debug output
        int numPoints = linePoints[0].length;

        // return a 2D numPoints x numDeps array of monochromatic *TOTAL* extinction line profiles
        double[][] logTotKappa = new double[numPoints][numDeps];

        double kappaL, logKappaC;

        // Set up multi-Gray continuum info:
        double isCool = 7300.0;  //Class A0

        //Set up multi-gray opacity:
        // lambda break-points and gray levels:
        // No. multi-gray bins = num lambda breakpoints +1
        double minLambda = 30.0;  //nm
        double maxLambda = 1.0e6;  //nm
        int maxNumBins = 11;
        double[][] grayLevelsEpsilons = MulGrayTCorr.grayLevEps(maxNumBins, minLambda, maxLambda, teff, isCool);
        //Find actual number of multi-gray bins:
        int numBins = 0; //initialization
        for (int i = 0; i < maxNumBins; i++) {
            if (grayLevelsEpsilons[0][i] < maxLambda) {
                numBins++;
            }
        }
        //put in multi-Gray opacities here:
        //Find which gray level bin the line falls in - assume that the first gray-level bin
        // is always at a shorter wavelength than the line!:
        int whichBin = 0;  //initialization
        for (int iB = 0; iB < numBins; iB++) {
            if (grayLevelsEpsilons[0][iB] >= lam0) {
                whichBin = iB;  //found it!
                break;
            }
        }

        for (int id = 0; id < numDeps; id++) {
            for (int il = 0; il < numPoints; il++) {
                //Both kappaL and kappa (continuum) are *mass* extinction (cm^2/g) at thsi point: 
                logKappaC = kappa[1][id]; // + Math.log(grayLevelsEpsilons[1][whichBin]); //kappaL = Math.exp(logKappaL[il][id]) + kappa[0][id];
                kappaL = Math.exp(logKappaL[il][id]) + Math.exp(logKappaC);
                logTotKappa[il][id] = Math.log(kappaL);
                //logTotKappa[il][id] = kappa[1][id];   //test - no line opacity
                //if (id == 36) {
                //    System.out.println("il " + il + " linePoints[0][il] " + 1.0e7*linePoints[0][il] + " logKappaL[il][id] " + logE*logKappaL[il][id] + " kappa[1][id] " + logE*kappa[1][id]);
                //}
            }
        }

        return logTotKappa;

    }

}
