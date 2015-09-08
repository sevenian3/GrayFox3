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
 * Returns vector of numDep linear geometric DEPTHS below top of atmosphere - in
 * cm (cgs) for consistency with log(g) units
 *
 * *May not be useful - scale depends on log(g)
 */
public class DepthScale {

    /**
     *
     * @param numDeps
     * @param tauRos
     * @param kappa
     * @param rho
     * @return
     */
    public static double[] depthScale(int numDeps, double[][] tauRos, double[][] kappa, double[][] rho) {

        double logE = Math.log10(Math.E); // for debug output

        //double ln10 = Math.log(10.0); //handy wee quantity 
        //log_10 Rosseland optical depth scale  
        double depths[] = new double[numDeps];

        // Upper bounday condition: 
        // Zero point at top of atmosphere - this can be shifted later?
        // log(z) cannot really correspond to zero 
        //double logZ1 = -10.0;  // log(cm)
        //depths[0] = Math.pow(10.0, logZ1);  //cm
        //Start at this depth index - the topmost layers have such low rhos that they correspond to HUUUGE geometric depths!
        int iStart = 10;
        double z1 = 0;   //cm
        //double z1 = -500.0 * 1.0e5; // FOr comparison to O&ASP 3rd Ed. (D.F. Gray), Table 9.2
        for (int i = 0; i <= iStart; i++) {
            depths[i] = z1;
        }

        //double minZ = 1.0E5; // = 1km - Minimum increase in depth from one point to the next
        // declare scratch variables
        //double deltaX, deltaZ, logZ2;
        double deltaX, deltaZ, z2, z3, help, logHelp, helpNext;
        //        h, k1, k2, k3, k4, logH, logK1, logK2, logK3, logK4;

        //Trapezoid method for depth at 2nd point in
        // Need to avoid using rho at upper boundary, so rho value must be taken at y_n+2 on all RHSs
        /*
         deltaX = tauRos[1][1] - tauRos[1][0];
         logHelp = tauRos[1][0] - kappa[1][0] - rho[1][2];
         System.out.format("%12.8f   %12.8f   %12.8f%n", logE*tauRos[1][0], logE*kappa[1][0], logE*rho[1][2]);
         //help = ( tauRos[0][0] / kappa[0][0] ) / rho[0][1];
         help = Math.exp(logHelp);
         */
        //First integrand:
        //deltaX = tauRos[1][iStart+1] - tauRos[1][iStart];
        logHelp = tauRos[1][iStart] - kappa[1][iStart] - rho[1][iStart];
        helpNext = Math.exp(logHelp);

//        deltaZ = (deltaX) * (0.5 * (help + helpNext));
        //       z2 = z1 + deltaZ;
//        depths[1] = z2;
        help = helpNext;

        //z1 =z2;
        for (int i = iStart + 1; i < numDeps; i++) {

            //Trapezoid method:
            deltaX = tauRos[1][i] - tauRos[1][i - 1];
            logHelp = tauRos[1][i] - kappa[1][i] - rho[1][i];
            helpNext = Math.exp(logHelp);
            //System.out.format("%12.8f   %12.8f   %12.8f%n", logE*tauRos[1][i], logE*kappa[1][i], logE*rho[1][i]);
            deltaZ = deltaX * (0.5 * (help + helpNext));
            //System.out.println("i " + i + " tauRos[1] " + logE*tauRos[1][i] + " kappa[1] " + logE*kappa[1][i] + " rho[1] " + logE*rho[1][i] + " deltaX " + deltaX + " deltaZ " + deltaZ);
            z2 = z1 + deltaZ;

            depths[i] = z2;
            z1 = z2;
            help = helpNext;

            //System.out.format("%12.8f   %12.8f%n", logE*tauRos[1][i], z2);

        }

        return depths;

    }

}
