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
 * Return the array index of the wavelength array (lambdas) closest to a desired
 * value of wavelength (lam)
 */
public class LamPoint {

    /**
     * 
     * @param numLams
     * @param lambdas
     * @param lam
     * @return 
     */
    public static int lamPoint(int numLams, double[] lambdas, double lam) {

        int index;

        double[] help = new double[numLams];

        for (int i = 0; i < numLams; i++) {

            help[i] = lambdas[i] - lam;
            help[i] = Math.abs(help[i]);

        }
        index = 0;
        double min = help[index];

        for (int i = 1; i < numLams; i++) {

            if (help[i] < min) {
                min = help[i];
                index = i;
            }

        }

        return index;

    }

}
