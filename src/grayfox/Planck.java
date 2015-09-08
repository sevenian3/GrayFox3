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
 * Inputs: lambda: a single scalar wavelength in cm temp: a single scalar
 * temperature in K Returns log of Plank function in logBBlam - B_lambda
 * distribution in pure cgs units: ergs/s/cm^2/ster/cm
 */
public class Planck {

    //public static double[][][] planck(int numDeps, double[] lamSetup, double[] lambda, double[][] temp){
    /**
     *
     * @param temp
     * @param lambda
     * @return
     */
    public static double planck(double temp, double lambda) {

        //int numLams = (int) (( lamSetup[1] - lamSetup[0] ) / lamSetup[2]) + 1; 
        double logBBlam; //, BBlam;

        //double c = Useful.c; //linear
        double logC = Useful.logC();  //log
        //double k = Useful.k;  //linear
        double logK = Useful.logK();  //log
        //double h = Useful.h;  //linear
        double logH = Useful.logH();  //log

        double logPreFac = Math.log(2.0) + logH + 2.0 * logC;  //log
        double logExpFac = logH + logC - logK;      //log
        //double preFac = 2.0 * h * ( c * c );  //linear
        //double expFac = ( h / k ) * c;      //linear

        //System.out.println("logC " + logC + " logK " + logK + " logH " + logH);
        //System.out.println("logPreFac " + logPreFac + " logExpFac " + logExpFac);
        //Declare scratch variables:
        double logLam, logPreLamFac, logExpLamFac, expon, logExpon, eTerm, denom, logDenom;  //log
        //double preLamFac, expLamFac, expon, denom; //linear

        //for (int il = 0; il < numLams; il++){
        //lambda = lambda[il] * 1.0E-7;  // convert nm to cm
        //lambda = lambda * 1.0E-7;  // convert nm to cm
        logLam = Math.log(lambda); // Do the call to log for lambda once //log
        //System.out.println("lambda " + lambda + " logLam " + logLam);

        logPreLamFac = logPreFac - 5.0 * logLam;  //log
        logExpLamFac = logExpFac - logLam;    //log
        //System.out.println("logPreLamFac " + logPreLamFac + " logExpLamFac " + logExpLamFac);
        // Be VERY careful about how we divide by lambda^5:
        //preLamFac = preFac / ( lambda * lambda ); //linear
        //preLamFac = preLamFac / ( lambda * lambda );  //linear
        //preLamFac = preLamFac / lambda;   //linear
        //expLamFac = expFac / lambda;

        //for (int id = 0; id < numDeps; id++){
        //logExpon = logExpLamFac - temp[1][id];
        
        //This is very subtle and dangerous!
        logExpon = logExpLamFac - Math.log(temp);  // log of hc/kTlambda
        //System.out.println("temp " + temp + " logTemp " + Math.log(temp));
        expon = Math.exp(logExpon);  // hc/kTlambda
        //System.out.println("logExpon " + logExpon + " expon " + expon + " denom " + denom);
        // expon = expLamFac / temp;  //linear
        eTerm = Math.exp(expon); // e^hc/ktlambda
        denom = eTerm - 1.0; // e^hc/ktlambda - 1
        logDenom = Math.log(denom); // log(e^hc/ktlambda - 1)

        //BBlam[1][id][il] = logPreLamFac - logDenom;
        //BBlam[0][id][il] = Math.exp(BBlam[1][id][il]);
        logBBlam = logPreLamFac - logDenom;  //log
        // Not needed? BBlam = Math.exp(logBBlam);  //log
        //BBlam = preLamFac / denom;  //linear

        // } //id loop - depths
        // } //il loop - lambdas
        return logBBlam;
    }

    // Computes the first partial derivative of B(T) wrt T, dB/dT:
        public static double dBdT(double temp, double lambda) {

        double logdBdTlam; 

        //double c = Useful.c; //linear
        double logC = Useful.logC();  //log
        //double k = Useful.k;  //linear
        double logK = Useful.logK();  //log
        //double h = Useful.h;  //linear
        double logH = Useful.logH();  //log

        double logPreFac = Math.log(2.0) + logH + 2.0 * logC;  //log
        double logExpFac = logH + logC - logK;      //log

        //Declare scratch variables:
        double logLam, logTemp, logPreLamFac, logExpLamFac, expon, logExpon, eTerm, denom, logDenom;  //log

        //lambda = lambda * 1.0E-7;  // convert nm to cm
        logLam = Math.log(lambda); // Do the call to log for lambda once //log
        logTemp = Math.log(temp);

        logPreLamFac = logPreFac + logExpFac - 6.0 * logLam - 2.0 * logTemp;  //log
        
        logExpLamFac = logExpFac - logLam;    //log
        
        //This is very subtle and dangerous!
        logExpon = logExpLamFac - logTemp;  // log of hc/kTlambda
        expon = Math.exp(logExpon);  // hc/kTlambda
       
        eTerm = Math.exp(expon); // e^hc/ktlambda
        denom = eTerm - 1.0; // e^hc/ktlambda - 1
        logDenom = Math.log(denom); // log(e^hc/ktlambda - 1)

        
        logdBdTlam = logPreLamFac + expon - 2.0 * logDenom;  //log
        
        return logdBdTlam;
    }

    
}
