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

/**
 * Mthod 1 This is dangerous - depends on depth scale calculations which is
 * flaky because of dependence on rho at upper boundary. Use METHOD 2 in
 * LineTau2 instead.
 *
 * Compute the monochromatic optical depth scale, tau_lambda, at each lambda
 * across the line profile Also computes line centre Continuum monochromatic
 * optical dpeth scale for continuum rectification And stores it in array
 * element numPoints+1 This may differ from the prescribed Tau_Ross because
 * kappa_Ross was NOT computed consistently with it Inputs: lineProf - only
 * neded for row 0 - wavelengths logKappaL - monochromatic line extinction
 * co-efficients kappa - we need the background continuum extinction
 * co-efficients, too!
 */
public class LineTau {

    /**
     * 
     * @param numDeps
     * @param lineProf
     * @param logKappaL
     * @param kappa
     * @param tauRos
     * @param rho
     * @param depths
     * @return 
     */
    public static double[][] tauLambda(int numDeps, double[][][] lineProf, double[][] logKappaL,
            double[][] kappa, double[][] tauRos, double[][] rho, double[] depths) {

        //No monochromatic optical depth can be less than the Rosseland optical depth,
        // so prevent zero tau_lambda values by setting each tau_lambda(lambda) at the 
        //top of the atmosphere to the tau_Ross value at the top 
        // - prevents trying to take a log of zero!
        double logE = Math.log10(Math.E); // for debug output
        double minTauL = tauRos[0][0];
        double minLogTauL = tauRos[1][0];

        int numPoints = lineProf[0].length;

        // returns numPoints+1 x numDeps array: the numPoints+1st row holds the line centre continuum tau scale
        double[][] logTauL = new double[numPoints + 1][numDeps];

        double kapTot1, logKapTot1, kapTot2, logKapTot2, kapTot4, logKapTot4,
                tau1, tau2, tau3, delta, tauL,
                logInteg1, integ1, logInteg2, logInteg4,
                h, k1, k2, k3, k4, logK1, logK2, logK3, logK4;

        for (int il = 0; il < numPoints; il++) {

            tau1 = minTauL; //initialize accumulator
            logTauL[il][0] = minLogTauL; // Set upper boundary TauL

            //System.out.println("LineTau: minTauL: " + minTauL);
            // Euler's method to get Tau at 2nd depth point in:
            //total extinction co-efficient
            kapTot1 = kappa[0][0] + Math.exp(logKappaL[il][0]);
            logKapTot1 = Math.log(kapTot1);

            delta = depths[1] - depths[0];

            //line Tau scale:
            logInteg1 = rho[1][0] + logKapTot1 + Math.log(delta);
            integ1 = Math.exp(logInteg1);

            //line Tau scale:
            tau2 = tau1 + integ1;
            logTauL[il][1] = Math.log(tau2);
            //tau1 = tau2;

            //if (il == 7) {
            //        //System.out.println("LineTau: id " + id + " k1 " + k1 + " k2 " + k2 + " k3 " + k3 + " k4 " + k4 + " log tau_l " + logE*logTauL[il][id]);
            //        System.out.println("LineTau: id 0 " + " integ1 " + integ1 + " log tau_l " + logE * logTauL[il][0]);
            //        System.out.println("LineTau: id 1 " + " integ1 " + integ1 + " log tau_l " + logE * logTauL[il][1]);  
            //    }
            for (int id = 2; id < numDeps; id++) {

                ////Euler's method:
                ////total extinction co-efficient
                //kapTot1 = kappa[0][id] + Math.exp(logKappaL[il][id]);
                //logKapTot1 = Math.log(kapTot1);
                //
                //delta = depths[id] - depths[id-1];
                //
                ////line Tau scale:
                //logInteg1 = rho[1][id] + logKapTot1 + Math.log(delta);
                //integ1 = Math.exp(logInteg1);
                //
                ////line Tau scale:
                //tau2 = tau1 + integ1;
                //logTauL[il][id] = Math.log(tau2);
                //tau1 = tau2;
                //
                // 4th order Runge-Kutte (mid-point), p. 705, Numerical Recipes in F77, 2nd Ed.
                h = depths[id] - depths[id - 2];

                kapTot1 = kappa[0][id - 2] + Math.exp(logKappaL[il][id - 2]);
                logKapTot1 = Math.log(kapTot1);
                kapTot2 = kappa[0][id - 1] + Math.exp(logKappaL[il][id - 1]);
                logKapTot2 = Math.log(kapTot2);
                kapTot4 = kappa[0][id] + Math.exp(logKappaL[il][id]);
                logKapTot4 = Math.log(kapTot4);

                //line Tau scale:
                logInteg1 = rho[1][id - 2] + logKapTot1 + Math.log(h);
                k1 = Math.exp(logInteg1);
                logInteg2 = rho[1][id - 1] + logKapTot2 + Math.log(h);
                k2 = Math.exp(logInteg2);
                k3 = k2;
                logInteg4 = rho[1][id] + logKapTot4 + Math.log(h);
                k4 = Math.exp(logInteg4);

                //line Tau scale:
                tau3 = tau1 + (k1 / 6.0) + (k2 / 3.0) + (k3 / 3.0) + (k4 / 6.0);
                if (tau3 <= tau2) {
                    tau3 = tau2 + minTauL;
                }
                logTauL[il][id] = Math.log(tau3);
                tau1 = tau2;
                tau2 = tau3;

                //logKapTot = Math.log(kapTot);
                //if ( il == (numPoints-1) ){
                //if (il == 7) {
                //    //System.out.println("LineTau: id " + id + " k1 " + k1 + " k2 " + k2 + " k3 " + k3 + " k4 " + k4 + " log tau_l " + logE*logTauL[il][id]);
                //    System.out.println("LineTau: id " + id + " integ1 " + integ1 + " log tau_l " + logE * logTauL[il][id]);
                //}
            } //id loop - depths

        } //il loop - lambdas

//Now compute the monochromatic line centre continuum optical depth scale and store it in an numPoints+1st column of
// logTauL array:
// This may differ from the prescribed Tau_Ross because kappa_Ross was NOT computed consistently with it 
        tau1 = minTauL; //initialize accumulator
        logTauL[numPoints][0] = minLogTauL; // Set upper boundary TauL

        // Euler's method to get Tau at 2nd depth point in:
        //total extinction co-efficient
        delta = depths[1] - depths[0];

        //line Tau scale:
        logInteg1 = rho[1][0] + kappa[1][0] + Math.log(delta);
        integ1 = Math.exp(logInteg1);

        //line Tau scale:
        tau2 = tau1 + integ1;
        logTauL[numPoints][1] = Math.log(tau2);
                //tau1 = tau2;

        //System.out.println("LineTau: minTauL: " + minTauL);
        for (int id = 2; id < numDeps; id++) {

            // Euler's method:
            //delta = depths[id] - depths[id - 1];
//
            //continuum tau scale:
            //          logInteg1 = rho[1][id] + kappa[1][id] + Math.log(delta);
            //        integ1 = Math.exp(logInteg1);
//
            //continuum tau scale:
            //          tau2 = tau1 + integ1;
            //        logTauL[numPoints][id] = Math.log(tau2);
            //      tau1 = tau2;
            // 4th order Runge-Kutte (mid-point), p. 705, Numerical Methods in F77, 2nd Ed.
            h = depths[id] - depths[id - 2];

            logInteg1 = rho[1][id - 2] + kappa[1][id - 2] + Math.log(h);
            k1 = Math.exp(logInteg1);
            logInteg2 = rho[1][id - 1] + kappa[1][id - 1] + Math.log(h);
            k2 = Math.exp(logInteg2);
            k3 = k2;
            logInteg4 = rho[1][id] + kappa[1][id] + Math.log(h);
            k4 = Math.exp(logInteg4);

            //line Tau scale:
            tau3 = tau1 + (k1 / 6.0) + (k2 / 3.0) + (k3 / 3.0) + (k4 / 6.0);
            if (tau3 <= tau2) {
                tau3 = tau2 + minTauL;
            }
            logTauL[numPoints][id] = Math.log(tau3);
            tau1 = tau2;
            tau2 = tau3;

            //if ( il == (numPoints-1) ){
            ////if ( il == 0 ){    
            //    System.out.println("LineTau: id " + id + " integ " + integ + " log tau_l " + logE*logTauL[il][id]);
            //}
        } //id loop - depths

        return logTauL;

    }

}
