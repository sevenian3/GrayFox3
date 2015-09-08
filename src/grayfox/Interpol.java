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

/*
 Linear interpolation to a new abscissa - mainly for interpolating flux to a specific lambda
 */
public class Interpol {

    /**
     *
     * @param x
     * @param y
     * @param newX
     * @return
     */
    public static double interpol(double[] x, double[] y, double newX) {

        double newY;

        // Bracket newX:
        double x1, x2;
        int p1, p2;
        p1 = 0;
        p2 = 1;
        x1 = x[p1];
        x2 = x[p2];

        for (int i = 1; i < x.length; i++) {
            if (x[i] >= newX) {
                // Found upper bracket
                p2 = i;
                p1 = i - 1;
                x2 = x[p2];
                x1 = x[p1];
                break;
            }
        }

        double step = x2 - x1;

    //Interpolate
        //First order Lagrange formula
        //   newY = y[1][p2] * (newX - x1) / step
        //           + y[1][p1] * (x2 - newX) / step;
        newY = y[p2] * (newX - x1) / step
                + y[p1] * (x2 - newX) / step;

        //System.out.println("Interpol: p1, p2, x1, x2, y1, y2, newX, newY: " + 
        //        p1 + " " + p2 + " " + x1 + " " + x2 + " " + y[1][p1] + " " + y[1][p2] + " " + newX + " " + newY + " ");
        return newY;

    }

}
