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
 * Computes the Gray kinetic temperature structure, on the Rosseland optical
 * depth scale T_kin(Tau_Ros) = Teff * (0.75Tau_Ros + Hopf)^0.25
 */
public class Temperature {

    /**
     * 
     * @param numDeps
     * @param teff
     * @param tauRos
     * @return 
     */
    public static double[][] temperature(int numDeps, double teff, double[][] tauRos) {

        //Gray kinetic temeprature structure:
        double[][] temp = new double[2][numDeps];

        double hopf, deltaLogTau, ii;

        for (int i = 0; i < numDeps; i++) {

            // Interpolate approximate Hopf function:
            deltaLogTau = (tauRos[1][i] - tauRos[1][0]) / (tauRos[1][numDeps - 1] - tauRos[1][0]);
            hopf = 0.55 + deltaLogTau * (0.710 - 0.55);

            //temp[1][i] = Math.log(teff) + 
            //             0.25 * Math.log(0.75*tauRos[0][i] + 0.5);
            temp[1][i] = Math.log(teff)
                    + 0.25 * Math.log(0.75 * (tauRos[0][i] + hopf));
            temp[0][i] = Math.exp(temp[1][i]);

        }

        return temp;

    }

}
