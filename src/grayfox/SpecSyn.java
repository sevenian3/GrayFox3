/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package grayfox;

/**
 *
 * Create master kappa_lambda(lambda) and tau_lambda(lambda) for
 * FormalSoln.formalSoln()
 *
 * @author Ian
 */
public class SpecSyn {

    //Merge comntinuum and line wavelength scales - for one line
    //This expects *pure* line opacity - no continuum opacity pre-added!
    public static double[] masterLambda(int numLams, int numMaster, int numNow, double[] masterLams, int numPoints, double[] listLineLambdas) {
        //                                 

        //int numCnt = lambdaScale.length;
        //skip the last wavelength point in the line lambda grid - it holds the line centre wavelength
        //int numLine = lineLambdas.length - 1;

        int numTot = numNow + numPoints; //current dynamic total

        //System.out.println("numCnt " + numCnt + " numLine " + numLine + " numTot " + numTot);
/*
         for (int i = 0; i < numCnt; i++) {
         System.out.println("i " + i + " lambdaScale[i] " + lambdaScale[i]);
         }
         for (int i = 0; i < numLine; i++) {
         System.out.println("i " + i + " lineLambdas[i] " + lineLambdas[i]);
         }
         */
        //Row 0 is merged lambda scale
        //Row 1 is log of *total* (line plus continuum kappa
        double[] masterLamsOut = new double[numTot];

        // Merge wavelengths into a sorted master list
        //initialize with first continuum lambda:
        double lastLam = masterLams[0];
        masterLamsOut[0] = masterLams[0];
        int nextCntPtr = 1;
        int nextLinePtr = 0;
        for (int iL = 1; iL < numTot; iL++) {
            if (nextCntPtr < numNow) {
                //System.out.println("nextCntPtr " + nextCntPtr + " lambdaScale[nextCntPtr] " + lambdaScale[nextCntPtr]);
                //System.out.println("nextLinePtr " + nextLinePtr + " lineLambdas[nextLinePtr] " + lineLambdas[nextLinePtr]);
                if ((masterLams[nextCntPtr] < listLineLambdas[nextLinePtr])
                        || (nextLinePtr >= numPoints - 1)) {
                    //Next point is a continuum point:
                    masterLamsOut[iL] = masterLams[nextCntPtr];
                    nextCntPtr++;

                } else if ((listLineLambdas[nextLinePtr] < masterLams[nextCntPtr])
                        && (nextLinePtr < numPoints - 1)) {
                    //Next point is a line point:
                    masterLamsOut[iL] = listLineLambdas[nextLinePtr];
                    nextLinePtr++;

                }
            }
            //System.out.println("iL " + iL + " masterLamsOut[iL] " + masterLamsOut[iL]);
        } //iL loop
        //Make sure final wavelength point in masterLams is secured:
        masterLamsOut[numTot-1] = masterLams[numNow-1];

        return masterLamsOut;
    }

    public static double[][] masterKappa(int numDeps, int numLams, int numMaster, int numNow, double[] masterLams, double[] masterLamsOut, double[][] logMasterKaps, int numPoints, double[] listLineLambdas, double[][] listLogKappaL) {
//                                          
        double logE = Math.log10(Math.E); // for debug output

        //int numLams = masterLams.length;
        int numTot = numNow + numPoints;
        
        double[][] logMasterKapsOut = new double[numTot][numDeps];
        //double[][] kappa2 = new double[2][numTot];
        //double[][] lineKap2 = new double[2][numTot];
        double kappa2, lineKap2, totKap;
        lineKap2 = 1.0e-99; //initialization

        //int numCnt = lambdaScale.length;
        //int numLine = lineLambdas.length - 1;
        double[] kappa1D = new double[numNow];
        double[] lineKap1D = new double[numPoints];
        //System.out.println("iL   masterLams    logMasterKappa");
        for (int iD = 0; iD < numDeps; iD++) {

            //Extract 1D *linear* opacity vectors for interpol()
            for (int k = 0; k < numNow; k++) {
                kappa1D[k] = Math.exp(logMasterKaps[k][iD]); //actually wavelength independent - for now
            }

            for (int k = 0; k < numPoints; k++) {
                lineKap1D[k] = Math.exp(listLogKappaL[k][iD]);
            }

            //Interpolate continuum and line opacity onto master lambda scale, and add them lambda-wise:
            for (int iL = 0; iL < numTot; iL++) {
                kappa2 = Interpol.interpol(masterLams, kappa1D, masterLamsOut[iL]);
                lineKap2 = 1.0e-99; //re-initialization
                if ( (masterLamsOut[iL] >= listLineLambdas[0]) && (masterLamsOut[iL] <= listLineLambdas[numPoints-1]) ) {
                    lineKap2 = Interpol.interpol(listLineLambdas, lineKap1D, masterLamsOut[iL]);
                    //lineKap2 = 1.0e-99;  //test
                }
                //test lineKap2 = 1.0e-99;  //test
                totKap = kappa2 + lineKap2;
                logMasterKapsOut[iL][iD] = Math.log(totKap);
                //if (iD == 36) {
                //    System.out.format("%02d   %12.8e   %12.8f%n", iL, masterLams[iL], logE * logMasterKappa[iL][iD]);
                //}
            }
        }

        return logMasterKapsOut;
    }

}
