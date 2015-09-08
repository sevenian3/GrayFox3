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
 * Compute the equivalent width of the Voigt line in pm - picometers NOTE: The
 * input parameter 'flux' should be a 2 x (numPoints+1) array where the
 * numPoints+1st value is the line centre monochromatic Continuum flux
 */
public class EqWidth {

    /**
     *
     * @param flux
     * @param linePoints
     * @param lam0
     * @return
     */
    public static double eqWidth(double[][] flux, double[][] linePoints, double lam0) {

        double Wlambda = 0.0; // Equivalent width in pm - picometers

        int numPoints = linePoints[0].length;
        //System.out.println("EqWidth: numPoints: " + numPoints);

        double delta, logDelta, term, normFlux, logNormFlux, integ, integ2, logInteg, lastInteg, lastTerm, term2;

        // Spectrum not normalized - try this instead (redefines input parameter fluxCont):
        double logFluxCont = Math.log((flux[0][0] + flux[0][numPoints - 1]) / 2.0);

        int iCount = (int) Math.floor(numPoints / 2) - 1; //initialize
        //console.log("numPoints " + numPoints + " iCount " + iCount);
        for (int il = (int) Math.floor(numPoints / 2); il < numPoints; il++) {
            //System.out.println("il " + il + " flux[0][il]/flux[0][numPoints - 1] " + flux[0][il]/flux[0][numPoints - 1]);
            if (flux[0][il] < 0.99 * flux[0][numPoints - 1]) {
                //console.log("Condition triggered");
                iCount++;
                //System.out.println("iCount " + iCount);
            }
        }
        //System.out.println("iCount " + iCount);
        //One more or two more if we can accomodate them:
        if (iCount < numPoints - 1) {
            iCount++;
        }
        if (iCount < numPoints - 1) {
            iCount++;
        }
        int iStart = numPoints - iCount;
        int iStop = iCount;
        //System.out.println("eqwidth: numPoints " + numPoints + " iStart " + iStart + " iStop " + iStop);
        //Trapezoid rule:       
        // First integrand:
        // Single-point normalization to line-centre flux suitable for narrow lines:
        //normFlux = flux[0][il] / flux[0][numPoints];
        //logNormFlux = flux[1][0] - fluxCont[1];
        logNormFlux = flux[1][iStart] - logFluxCont;
        //normFlux = flux[0][il] / fluxCont[0];
        //System.out.println("flux[0][il] " + flux[0][il] + " fluxCont[0] " + fluxCont[0]);
        if (logNormFlux >= -0.01) {
            lastInteg = 1.0e-99;
        } else {
            lastInteg = 1.0 - Math.exp(logNormFlux);
        }
        lastTerm = lastInteg; //initialization

        for (int il = iStart + 1; il < iStop; il++) {

            delta = linePoints[0][il] - linePoints[0][il - 1];
            delta = delta * 1.0E+7;  // cm to nm - W_lambda in pm
            logDelta = Math.log(delta);

            // Single-point normalization to line-centre flux suitable for narrow lines:
            //normFlux = flux[0][il] / flux[0][numPoints];
            //logNormFlux = flux[1][il] - fluxCont[1];
            logNormFlux = flux[1][il] - logFluxCont;
            //normFlux = flux[0][il] / fluxCont[0];
            //System.out.println("flux[0][il] " + flux[0][il] + " fluxCont[0] " + fluxCont[0]);
            // flux should be less than 0.99 of continuum flux:
            if (logNormFlux >= -0.01) {
                integ = 1.0e-99;
            } else {
                integ = 1.0 - Math.exp(logNormFlux);
            }

            //Trapezoid rule:
            integ2 = 0.5 * (lastInteg + integ);
            logInteg = Math.log(integ2);
            term = Math.exp(logInteg + logDelta);

            //Make sure weird features near the red edge don't pollute our Wlambda:
            // for lambda > line centre, area sould be monotically *decreasing*
            //console.log("il " + il + " linePoints[0][il] " + linePoints[0][il] + " lam0 " + lam0 + " integ " + integ + " lastInteg " + lastInteg);
            if ((linePoints[0][il] > 0.0) && (term > lastTerm)) {
                //console.log("Condition triggered");
                //term2 = lastTerm / 2.0;
                term2 = term; //above condition gives too small EWs
            } else {
                term2 = term;
            }

            //Wlambda = Wlambda + (term * delta);
            Wlambda = Wlambda + term2;

            lastTerm = term; //For catching problems
            lastInteg = integ;

            //System.out.println("EqWidth: il " + il + " delta " + delta + " term " + term + " normFlux " + normFlux );
            //System.out.println("EqWidth: Wlambda: " + Wlambda);
        }

        // Convert area in nm to pm - picometers
        Wlambda = Wlambda * 1.0E3;

        return Wlambda;

    }

}
