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
 * Return the minimum and maximum values of an input 1D array CAUTION; Will
 * return the *first* occurence if min and/or max values occur in multiple
 * places iMinMax[0] = first occurence of minimum iMinMax[1] = first occurence
 * of maximum
 */
public class MinMax {

    /**
     * 
     * @param x
     * @return 
     */
    public static int[] minMax(double[] x) {

        int[] iMinMax = new int[2];

        int num = x.length;
        //System.out.println("MinMax: num: " + num);

        int iMin = 0;
        int iMax = 0;
        double min = x[iMin];
        double max = x[iMax];

        for (int i = 1; i < num; i++) {

            //System.out.println("MinMax: i , current min, x : " + i + " " + min + " " + x[i]);
            if (x[i] < min) {
                //System.out.println("MinMax: new min: if branch triggered" );
                min = x[i];
                iMin = i;
            }
            //System.out.println("MinMax: new min: " + min);

            if (x[i] > max) {
                max = x[i];
                iMax = i;
            }

        }
        //System.out.println("MinMax: " + iMin + " " + iMax);

        iMinMax[0] = iMin;
        iMinMax[1] = iMax;

        return iMinMax;

    }

}
