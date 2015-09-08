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

import java.text.DecimalFormat;

public class Flux {

    //public static double[][] flux(double[][] intens, double[][] cosTheta, int numLams, double[] lambdas) {
    /**
     * @param intens
     * @param cosTheta
     * @return
     */
    public static double[] flux(double[] intens, double[][] cosTheta) {

        //double[][] fluxSurfSpec = new double[2][numLams];
        double[] fluxSurfSpec = new double[2];
        //double fluxSurfBol, logFluxSurfBol, lambda2, lambda1;  // Bolometric quantities for reality check

        // returns surface flux as a 2XnumLams vector
        //  - Row 0 linear flux (cgs units)
        //  - Row 1 log_e flux
        //  cosTheta is a 2xnumThetas array:
 
        // Gaussian quadrature
 
        int numThetas = cosTheta[0].length;

        //fluxSurfBol = 0;
        //for (int il = 0; il < numLams; il++) {
        double flx = 0.0;

        for (int it = 0; it < numThetas; it++) {

            //flx = flx + intens[il][it] * cosTheta[1][it] * cosTheta[0][it];
            flx = flx + intens[it] * cosTheta[1][it] * cosTheta[0][it];

        }  // it - theta loop

            //fluxSurfSpec[0][il] = 2.0 * Math.PI * flx;
        //fluxSurfSpec[1][il] = Math.log(fluxSurfSpec[0][il]);
        fluxSurfSpec[0] = 2.0 * Math.PI * flx;
        fluxSurfSpec[1] = Math.log(fluxSurfSpec[0]);

        /* Can no longer do this test here:
         if (il > 1) {
         lambda2 = lambdas[il]; // * 1.0E-7;  // convert nm to cm
         lambda1 = lambdas[il - 1]; // * 1.0E-7;  // convert nm to cm
         fluxSurfBol = fluxSurfBol
         + fluxSurfSpec[0][il] * (lambda2 - lambda1);
         }
         */
        //} //il - lambda loop

        /* Can no longer do this test here:
         logFluxSurfBol = Math.log(fluxSurfBol);
         double logTeff = (logFluxSurfBol - Useful.logSigma()) / 4.0;
         double teff = Math.exp(logTeff);

         String pattern = "0000.00";
         //String pattern = "#####.##";
         DecimalFormat myFormatter = new DecimalFormat(pattern);

         System.out.println("FLUX: Recovered Teff = " + myFormatter.format(teff));
         */
        return fluxSurfSpec;

    }

}
