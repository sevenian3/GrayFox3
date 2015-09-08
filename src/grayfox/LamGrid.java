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

public class LamGrid {
    
    /**
     * 
     * @param numLams
     * @param lamSetup
     * @return 
     */
    public static double[] lamgrid(int numLams, double[] lamSetup){
    
     
    double[] lambdaScale = new double[numLams];
    double logLambda;
    
    // Space lambdas logarithmically:
    double logLam1 = Math.log10(lamSetup[0]);
    double logLam2 = Math.log10(lamSetup[1]);
    double delta = ( logLam2 - logLam1 ) / numLams;
    
    double ii;
    for ( int i = 0; i < numLams; i++){
        
        ii = (double) i;
        logLambda = logLam1 + ( ii * delta );
        lambdaScale[i] = Math.pow(10.0, logLambda);
        
        //System.out.println("il " + i + " lambda: " + lambdaScale[i]); //debug
        
    }
        
    return lambdaScale;
    
    }
    
}
