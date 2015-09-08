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
 * Find the wavelength where a spectral line profile is half of its minimum,
 * line-centre brightness Assume symmetric (pure Voigt) profile?
 *
 * Returns the integer indices of the red half-power and line-centre wavelengths
 */
public class HalfPower {

    /**
     * 
     * @param numPoints
     * @param lineFlux
     * @return 
     */
    public static int[] halfPower(int numPoints, double[][] lineFlux) {

        int[] keyLambdas = new int[2];

    // CAUTION; numPoints-1st value holds the line centre monochromatic *continuum* flux for normalization
        // Extract 1D vector of linear continuum normalized fluxes of length numPoints-1
        double[] flux1D = new double[numPoints - 1];

        for (int i = 0; i < numPoints - 1; i++) {
            flux1D[i] = lineFlux[0][i] / lineFlux[0][numPoints - 1];
            //System.out.println("HalfPower: i, fluxiD: " + i + " " + flux1D[i] );
        }

    // To be on the safe side, let's "rediscover" the line centre of wavelength of minimum brightness:
        int[] minmax = MinMax.minMax(flux1D);
        keyLambdas[0] = minmax[0]; // line centre index
        //System.out.println("HalfPower: numPoints, line centre minmax[0], line centre flux: " 
        //        + numPoints + " " + keyLambdas[0] + " " + flux1D[keyLambdas[0]] );

        double[] help = new double[numPoints - keyLambdas[0] - 1];
        double half;
        half = flux1D[keyLambdas[0]] + ((1.0 - flux1D[keyLambdas[0]]) / 2.0);
    //System.out.println("HalfPower: half power flux: " + half);

        for (int i = keyLambdas[0]; i < numPoints - 1; i++) {

            // The last minimum of help should be the red half-depth point
            help[i - keyLambdas[0]] = Math.abs(flux1D[i] - half);
        //System.out.println("HalfPower: i, i - keyLambdas[0], fluxiD, help: " 
            //        + i + " " + (i - keyLambdas[0]) + " " + flux1D[i] + " " + help[i - keyLambdas[0]] );
        }

        minmax = MinMax.minMax(help);
        keyLambdas[1] = minmax[0] + keyLambdas[0]; // red half-power index
        //System.out.println("HalfPower: minmax[0]: " + keyLambdas[1]);

        return keyLambdas;

    }

}
