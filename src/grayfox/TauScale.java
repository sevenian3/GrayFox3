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

public class TauScale {
    
    // CAUTION: Here tau[1][] is log_10!
    /**
     * 
     * @param numDeps
     * @return 
     */
    public static double[][] tauScale(int numDeps, double log10MinDepth, double log10MaxDepth){
        
    //log_10 Rosseland optical depth scale  
        double tauRos[][] = new double[2][numDeps]; 
        
        // Construct the log ROsseland optical depth scale:
        // Try equal spacing in log depth
        
        double ln10 = Math.log(10.0);
        
//        double log10MinDepth = -4.5;
//        double log10MaxDepth = 1.5;
        
        double logMinDepth = log10MinDepth * ln10;
        double logMaxDepth = log10MaxDepth * ln10;
        
        double deltaLogTau = (logMaxDepth - logMinDepth)/(numDeps - 1.0);
        
        double ii;
        for (int i = 0; i < numDeps; i++){
            
            ii = (double)i;
            tauRos[1][i] = logMinDepth + ii*deltaLogTau;
            tauRos[0][i] = Math.exp(tauRos[1][i]);
            //System.out.println("i: " + i + " absTauDiff[1][i] " + tauRos[1][i] + " tauRos[0][i] " + tauRos[0][i]);
        }
        
        return tauRos;
        
    }
    
}
