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

public class Thetas {

    /**
     *
     * @return
     */
    public static double[][] thetas() {

        //int numThetas = 10; // guess
        //double[] cosTheta = new double[numThetas];
        // Try equal distribution in cos(theta) space (rather than Gaussian quadrature)
       /* 
         double ii;
         for (int i = 0; i < numThetas; i++){
            
         ii = (double) i;
         cosTheta[i] = 1.0 - ii/numThetas;
            
         }
         */
        //  cosTheta is a 2xnumThetas array:
        // row 0 is used for Gaussian quadrature weights
        // row 1 is used for cos(theta) values
        // Gaussian quadrature:
        /*
         "n = 21" Gaussian quadrature weights, w_i, and abscissae from 
         http://pomax.github.io/bezierinfo/legendre-gauss.html
         - ie. 11 point among 0 and positive abcissae
        
         This 11/21 of a 21-point formula: 0 plus the positive abscissae ,
         so I *think* it represents *half* the integral on the interval [-1,1],
         ie., on the interval[0,1].   SO: Divide the first weight (x=0) by 2 so 
         that the quadrature sum is exactly half the integral on [-1,1]. 
         */
        //int nGauss = 11;
        //int nGauss = 7;
        int nGauss = 9;
        double[] theta = new double[nGauss];
        double[] weight = new double[nGauss];
        double[][] cosTheta = new double[2][nGauss];

        // I *think* the "thetas" being assigned here (abcissae) are fractional
        // angles, theta/(pi/2).
        // 11 points on [0,1] from 21 point Gaussian quadrature on [-1,1]
        //weight[0] = 0.5 * 0.1460811336496904;  // Divide the weight of the x=0 point by 2!  
//        weight[0] = 0.1460811336496904; 
//        theta[0] = 0.0000000000000000;  //disk centre
//        weight[1] = 0.1445244039899700;
//        theta[1] = 0.1455618541608951;
//        weight[2] = 0.1398873947910731;
//        theta[2] = 0.2880213168024011;
//        weight[3] = 0.1322689386333375;
//        theta[3] = 0.4243421202074388;
//        weight[4] = 0.1218314160537285;
//        theta[4] = 0.5516188358872198;
//        weight[5] = 0.1087972991671484;
//        theta[5] = 0.6671388041974123;
//        weight[6] = 0.0934444234560339;
//        theta[6] = 0.7684399634756779;
//        weight[7] = 0.0761001136283793;
//        theta[7] = 0.8533633645833173;
//        weight[8] = 0.0571344254268572;
//        theta[8] = 0.9200993341504008;
//        weight[9] = 0.0369537897708525;
//        theta[9] = 0.9672268385663063;
//        weight[10] = 0.0160172282577743;
//        theta[10] = 0.9937521706203895;  //limb
        
               // 7 points on [0,1] from 13 point Gaussian quadrature on [-1,1]
 //       weight[0] = 0.2325515532308739;  // disk center
        //       theta[0] = 0.0000000000000000;
        //       weight[1] = 0.2262831802628972; 
//        theta[1] = 0.2304583159551348;
        //       weight[2] = 0.2078160475368885; 
        //       theta[2] = 0.4484927510364469;
        //      weight[3] = 0.1781459807619457; 
        //       theta[3] = 0.6423493394403402;
        //       weight[4] = 0.1388735102197872; 
        //       theta[4] = 0.8015780907333099;
        //       weight[5] = 0.0921214998377285; 
        //       theta[5] = 0.9175983992229779;
        //       weight[6] = 0.0404840047653159; 
        //       theta[6] = 0.9841830547185881;  //limb
        
              // 9 points on [0,1] from 17 point Gaussian quadrature on [-1,1]    
        weight[0] = 0.1794464703562065;  //disk center
        theta[0] = 0.0000000000000000;
        weight[1] = 0.1765627053669926;
        theta[1] = 0.1784841814958479;
        weight[2] = 0.1680041021564500;
        theta[2] = 0.3512317634538763;
        weight[3] = 0.1540457610768103;
        theta[3] = 0.5126905370864769;
        weight[4] = 0.1351363684685255;
        theta[4] = 0.6576711592166907;
        weight[5] = 0.1118838471934040;
        theta[5] = 0.7815140038968014;
        weight[6] = 0.0850361483171792;
        theta[6] = 0.8802391537269859;
        weight[7] = 0.0554595293739872;
        theta[7] = 0.9506755217687678;
        weight[8] = 0.0241483028685479;
        theta[8] = 0.9905754753144174;  //limb

        for (int it = 0; it < nGauss; it++) {
            cosTheta[0][it] = weight[it];
            theta[it] = theta[it] * Math.PI / 2.0;
            cosTheta[1][it] = Math.cos(theta[it]);
        }

        // Try equal distribution in cos(theta) space (rather than Gaussian quadrature)
       /* 
         double ii;
         for (int i = 0; i < numThetas; i++){
            
         ii = (double) i;
         cosTheta[i] = 1.0 - ii/numThetas;
            
         }
         */
        return cosTheta;
    }

}
