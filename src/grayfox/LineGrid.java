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
 * Line profile, phi_lambda(lambda): Assume Voigt function profile - need H(a,v)
 * Assumes CRD, LTE, ??? Input parameters: lam0 - line center wavelength in nm
 * mass - mass of absorbing particle (amu) logGammaCol - log_10(gamma) - base 10
 * logarithmic collisional (pressure) damping co-efficient (s^-1) epsilon -
 * convective microturbulence- non-thermal broadening parameter (km/s) Also
 * needs atmospheric structure information: numDeps WON'T WORK - need observer's
 * frame fixed lambda at all depths: temp structure for depth-dependent thermal
 * line broadening Teff as typical temp instead of above pressure structure,
 * pGas, if scaling gamma
 */
public class LineGrid {

    /**
     *
     * @param lam0In
     * @param massIn
     * @param xiTIn
     * @param numDeps
     * @param teff
     * @return
     */
    public static double[][] lineGrid(double lam0In, double massIn, double xiTIn,
            int numDeps, double teff, int numCore, int numWing) {

        double c = Useful.c;
        double logC = Useful.logC();
        //double k = Useful.k;
        double logK = Useful.logK();
        //double e = Useful.e;
        //double mE = Useful.mE;
        double amu = Useful.amu;

        double ln10 = Math.log(10.0);
        double ln2 = Math.log(2.0);

        double logE = Math.log10(Math.E); // for debug output

        //Put input parameters into linear cgs units:
        //double gammaCol = Math.pow(10.0, logGammaCol);
        double logTeff = Math.log(teff);

        double xiT = xiTIn * 1.0E5; //km/s to cm/s
        double lam0 = lam0In; // * 1.0E-7; //nm to cm
        double logLam0 = Math.log(lam0);
        double logMass = Math.log(massIn * amu);  //amu to g 

        // Compute depth-independent Doppler width, Delta_lambda_D:
        double doppler, logDopp;
        double logHelp, help; //scratch

        logHelp = ln2 + logK + logTeff - logMass; // M-B dist, square of v_mode
        help = Math.exp(logHelp) + xiT * xiT; // quadratic sum of thermal v and turbulent v
        logHelp = 0.5 * Math.log(help);
        logDopp = logHelp + logLam0 - logC;

        doppler = Math.exp(logDopp);  // cm

        //System.out.println("LineGrid: doppler, logDopp: " + doppler + " " + logE*logDopp);
        //Set up a half-profile Delta_lambda grid in Doppler width units 
        //from line centre to wing
        //int numCore = 5;
        //int numWing = 5;
        //int numWing = 0;  //debug
        int numPoints = numCore + numWing;

        // a 2D 2 X numPoints array of Delta Lambdas 
        // Row 0 : Delta lambdas in cm - will need to be in nm for Planck and Rad Trans?
        // Row 1 : Delta lambdas in Doppler widths
        double[][] linePoints = new double[2][numPoints];

        // Line profiel points in Doppler widths - needed for Voigt function, H(a,v):
        double[] v = new double[numPoints];

        double maxCoreV = 3.5; //core half-width ~ in Doppler widths
        //double maxWingDeltaLogV = 1.5 * ln10; //maximum base e logarithmic shift from line centre in Doppler widths
        double minWingDeltaLogV = Math.log(maxCoreV + 1.5);
        double maxWingDeltaLogV = 9.0 + minWingDeltaLogV;

        double logV, ii, jj;

        for (int il = 0; il < numPoints; il++) {

            ii = (double) il;

            if (il < numCore) {

                // In core, space v points linearly:
                // Voigt "v" parameter
                // v > 0 --> This is the *red* wing:
                v[il] = ii * maxCoreV / (numCore - 1);
                linePoints[0][il] = doppler * v[il];
                linePoints[1][il] = v[il];

            } else {

                //Space v points logarithmically in wing
                jj = ii - numCore;
                logV = (jj * (maxWingDeltaLogV - minWingDeltaLogV) / (numPoints - 1)) + minWingDeltaLogV;
                v[il] = Math.exp(logV);
                linePoints[0][il] = doppler * v[il];
                linePoints[1][il] = v[il];

            } // end else

            //System.out.println("LineGrid: il, lam, v: " + il + " " + 
            //        linePoints[0][il] + " " + linePoints[1][il]);
        } // il lambda loop

        // Add the negative DeltaLambda half of the line:
        int numPoints2 = (2 * numPoints) - 1;
        //System.out.println("LineGrid: numpoints2: " + numPoints2);

        // Return a 2D 2 X (2xnumPoints-1) array of Delta Lambdas 
        // Row 0 : Delta lambdas in cm - will need to be in nm for Planck and Rad Trans?
        // Row 1 : Delta lambdas in Doppler widths
        double[][] linePoints2 = new double[2][numPoints2];

        //wavelengths are depth-independent - just put them in the 0th depth slot:
        for (int il2 = 0; il2 < numPoints2; il2++) {

            if (il2 < numPoints - 1) {

                int il = (numPoints - 1) - il2;
                linePoints2[0][il2] = -1.0 * linePoints[0][il];
                linePoints2[1][il2] = -1.0 * linePoints[1][il];

            } else {

                //Positive DelataLambda half:   
                int il = il2 - (numPoints - 1);
                linePoints2[0][il2] = linePoints[0][il];
                linePoints2[1][il2] = linePoints[1][il];

            }

            //System.out.println("LineGrid: il2, lam, v: " + il2 + " " + 
            //        linePoints2[0][il2] + " " + linePoints2[1][il2]);
        } //il2 loop

        return linePoints2;

    }

}
