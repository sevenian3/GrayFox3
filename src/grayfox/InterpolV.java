/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package grayfox;

/**
 *
 * @author Ian vectorized version of simple linear 1st order interpolation
 * // Caution: Assumes new abscissae to which we're interpolating are entirey within the 
 *  range of the old abscissae
 */
public class InterpolV {

    public static double[] interpolV(double[] y, double[] x, double[] newX) {

        int num = newX.length;
        double[] newY = new double[num];

        int j = 0; //initialize old abscissae index
        //outer loop over new acscissae
        for (int i = 0; i < num; i++) {

            if (x[j] < newX[i]) {
                j++;
            }
            j--; //Passed the first newX - back up one

            //1st order Lagrange method:
            newY[i] = y[j+1]*(newX[i]-x[j]) + y[j]*(x[j+1]-newX[i]);
        }

        return newY;
    }

}
