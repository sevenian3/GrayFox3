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
 * Solve hydrostatic eq for P scale on the tau scale - need to pick a depth
 * dependent kappa value! - dP/dTau = g/kappa --> dP/dlogTau = Tau*g/kappa press
 * is a 4 x numDeps array: rows 0 & 1 are linear and log *gas* pressure,
 * respectively rows 2 & 3 are linear and log *radiation* pressure Split
 * pressure into gas and radiation contributions as we calculate it:
 */
public class Hydrostat {

    /**
     *
     * @param numDeps
     * @param grav
     * @param tauRos
     * @param kappa
     * @param temp
     * @return
     */
    public static double[][] hydrostatic(int numDeps, double grav, double[][] tauRos, double[][] kappa, double[][] temp) {

        //double c = Useful.c;
        double logC = Useful.logC();
        //double sigma = Useful.sigma;
        double logSigma = Useful.logSigma();

        //double c = 9.9989E+10; // light speed in vaccuum in cm/s
        //double sigma = 5.670373E-5;   //Stefan-Boltzmann constant ergs/s/cm^2/K^4   
        double radFac = Math.log(4.0) + logSigma - Math.log(3.0) - logC;

        double logRadiusSun = 0.0;  //solar units

        double press[][] = new double[4][numDeps];

        //double ln10 = Math.log(10.0); //handy wee quantity       
        //Upper boundary condition: total pressure at top of atmosphere
        double p1 = 1.0E-4;

        // System.out.println("HYDROSTAT: ln10= " + ln10 + " p1 " + p1 + "\r\n");
        //Finite differences in log(Tau) space - deltaX should be uniform, 
        //   but compute it on the fly anyway in case we play with other tau scales
        press[0][0] = p1;
        press[1][0] = Math.log(p1);
        press[2][0] = p1;
        press[3][0] = Math.log(p1);

        // Decalare scratch variables:
        double deltaX, deltaP, help, p2, p3, thisKap, logThisKap;
        double logPrad, pRad, helpSub, h, k1, k2, k3, k4;

        // Calculate P at the 2nd depth point in using Euler's method:
        //deltaX = tauRos[1][1] - tauRos[1][0];
        //help = (tauRos[0][0] / kappa[0][0]) * grav;
        //deltaP = help * (deltaX);
        //p2 = p1 + deltaP;
        // Compute LTE bolometric radiation contribution to total HSE pressure
        //logPrad = radFac + 4.0 * temp[1][1];
        //pRad = Math.exp(logPrad);
        //// Avoid zero or negative Pgas values in subtraction below:
        //if (pRad >= 0.99 * p2) {
        //    pRad = 0.99 * p2;
        //}
        // Avoid a direct subtraction in case Prad is close to Pgas for deeper 
        // layers of hotter stars, and both values are large:
        //pGas = p2 - pRad;
        //helpSub = 1.0E0 - (pRad / p2);
        //press[0][1] = helpSub * p2;
        //press[1][1] = Math.log(press[0][1]);
        //press[2][1] = pRad;
        //press[3][1] = Math.log(pRad);
        //System.out.println("HYDROSTAT: i " + i + " Pgas " + press[0][i] + " Prad " + pRad);
        //Set lower boundary of next step:
        //p1 = p2;
//RK4      for (int i = 2; i < numDeps; i++) {   //RK4
        for (int i = 1; i < numDeps; i++) {

            // Euler's method:
            // Delta log(tau):
            deltaX = tauRos[1][i] - tauRos[1][i - 1];
            //// I have no idea why this is necessary, bu tit seems to be needed to get the right pressure scaling with logg:
            //logThisKap = kappa[1][i] - (logRadiusSun - logRadius);             
            //thisKap = Math.exp(logThisKap);
            thisKap = kappa[0][i];
            help = (tauRos[0][i] / thisKap) * grav;
            //// deltaP = help * ( deltaX ); //log10
            deltaP = help * (deltaX);
            p2 = p1 + deltaP;
            //// 4th order Runge-Kutte (mid-point), p. 705, Numerical Recipes in F77, 2nd Ed.
            //h = tauRos[1][i] - tauRos[1][i - 2];
            //k1 = h * (tauRos[0][i - 2] / kappa[0][i - 2]) * grav;
            //k2 = h * (tauRos[0][i - 1] / kappa[0][i - 1]) * grav;
            //k3 = k2;
            //k4 = h * (tauRos[0][i] / kappa[0][i]) * grav;
///
            //p3 = p1 + (k1 / 6.0) + (k2 / 3.0) + (k3 / 3.0) + (k4 / 6.0);

            //System.out.println("HYDROSTAT: i " + i + " deltaX " + deltaX + 
            //                   " help " + help + " deltaP " + deltaP + " p1 " + p1 + " p2 " + p2);
            // Compute LTE bolometric radiation contribution to total HSE pressure
            logPrad = radFac + 4.0 * temp[1][i];
            pRad = Math.exp(logPrad);

            // Avoid zero or negative Pgas values in subtraction below:
            if (pRad >= 0.99 * p2) {
                pRad = 0.99 * p2;
            }  //Euler
            //if (pRad >= 0.99 * p3) {
            //    pRad = 0.99 * p3;
            //} // 2nd O R-K

            // Avoid a direct subtraction in case Prad is close to Pgas for deeper 
            // layers of hotter stars, and both values are large:
            //pGas = p2 - pRad;
            helpSub = 1.0E0 - (pRad / p2); //Euler
            //helpSub = 1.0E0 - (pRad / p3);  // 2nd O R-K

            press[0][i] = helpSub * p2;  //Euler
            //press[0][i] = helpSub * p3;  // 2nd O R-K
            press[1][i] = Math.log(press[0][i]);
            press[2][i] = pRad;
            press[3][i] = Math.log(pRad);
            //press[2][i] = 0.0;  //test
            //press[3][i] = -99.0;  //test

            //System.out.println("HYDROSTAT: temp " + temp[0][i] + " Pgas " + press[0][i] + " Prad " + pRad);
            //Set lower boundary of next step:
            //p1 = p2; //Euler
            p1 = p2;  // 2nd O R-K
            //p2 = p3;  // 2nd O R-K

        }

        return press;

    }

}
