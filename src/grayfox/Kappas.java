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
 * Compute Rosseland mean extinction coefficient (cm^2/g) structure by scaling
 * from Sun
 *
 */
public class Kappas {

    /**
     *
     * @param numDeps
     * @param kappaScale
     * @param teff
     * @param teffSun
     * @param logg
     * @param loggSun
     * @return
     */
    //** Input parameter mode: valid values are 0 or 1
    // If mode = 0: We don't yet know the mass density, rho & rhoSun, yet - the rhos passed in are meaningless and rhos must be faked
    // If mode = 1" We already have a previous in situ calculation of mass density, rho & rhoSun - passed in as parameters
    //public static double[][] kappas(int numDeps, double kappaScale,  double teff, double teffSun, double logg, double loggSun) {
    public static double[][] kappas(int mode, int numDeps, double[][] rho, double[][] rhoSun, double[][] kappaRosSun, double kappaScale, double logg, double loggSun, double teff, double teffSun, double radius,
            double massX, double massZ, double[][] tauRos, double[][] temp, double[][] tempSun, double[][] logNumsH3, double[][] logNumsH2) {

        double[][] kappa = new double[2][numDeps];

        double hotT = 6000.0;  //hotter than this in K and we use hot star formula

        double logRadiusSun = 0.0;  //solar units
        double massZSun = 0.02;

        double dilute, rhoStarFake, rhoSunFake;
        double logRhoStarFake = 0.0; //enforced initialization
        double logRhoSunFake = 0.0;  //enforced initialization

        if (mode == 0) {
            //We don't yet know rho & rhoSun - the rhos passed in are meaningless - fake up rhos:
            // Approximate mass density in atmosphere by scaling with logg and radius, then diluting:
            dilute = 5.0e-5; //tuned to give rho ~10^-1 g/cm^-3 in Sun's atmosphere
            logRhoStarFake = Math.log(3.0 / 4.0 / Math.PI) - Useful.logGConst() + logg - Math.log(Useful.rSun * radius);
            rhoStarFake = dilute * Math.exp(logRhoStarFake);

            // Do the same for Sun for consistency
            logRhoSunFake = Math.log(3.0 / 4.0 / Math.PI) - Useful.logGConst() + loggSun - Useful.logRSun();
            rhoSunFake = dilute * Math.exp(logRhoSunFake);
        }

        //System.out.println("rhoSunFake: " + rhoSunFake + " rhoStarFake: " + rhoStarFake);
        // System.out.println("logHelp " + logHelp + " help " + help + " kappaScale out " + kappaScale); //debug
        /*
         double[][] kappaRosSun = new double[2][numDeps];

        
         double minLog10KappaRosSun = -3.5;
         double maxLog10KappaRosSun = 2.0;

         double ln10 = Math.log(10.0);
         double minLogKappaRosSun = minLog10KappaRosSun * ln10;
         double maxLogKappaRosSun = maxLog10KappaRosSun * ln10;

         double deltaKappa = (maxLogKappaRosSun - minLogKappaRosSun) / numDeps;

         double ii;

         //Sun:
         for (int i = 0; i < numDeps; i++) {

         ii = (double) i;
         kappaRosSun[1][i] = minLogKappaRosSun + ii * deltaKappa;
         kappaRosSun[0][i] = Math.exp(kappaRosSun[1][i]);

         }
         */
        //Star:
        double numerator = 1.0;  //enforced initialization
        double denominator = 1.0;  //enforced initialization

        double logHelp, help, reScale, logNH3, logNH2;
        //double kappaOld, logKappaOld;
        reScale = 1.0 * kappaScale;
        //System.out.println("kappaScale: " + kappaScale);
        for (int i = 0; i < numDeps; i++) {
            logNH3 = logNumsH3[2][i];
            logNH2 = logNumsH2[2][i];
            System.out.println("i " + i);
            if (mode == 0) {
                numerator = kappaFac(numDeps, hotT, logRhoStarFake, temp[1][i], massX, massZ, logNH3, logNH2);
                denominator = kappaFac(numDeps, hotT, logRhoSunFake, tempSun[1][i], massX, massZSun, logNH3, logNH2);
            } else if (mode == 1) {
                numerator = kappaFac(numDeps, hotT, rho[1][i], temp[1][i], massX, massZ, logNH3, logNH2);
                //numerator = kappaFac(numDeps, hotT, rhoSun[1][i], temp[1][i], massX, massZ);
                denominator = kappaFac(numDeps, hotT, rhoSun[1][i], tempSun[1][i], massX, massZSun, logNH3, logNH2);
            }
            //System.out.println("i " + i + " kappaRosSun[0][i] " + kappaRosSun[0][i]);
            kappa[0][i] = reScale * kappaRosSun[0][i] * (numerator / denominator);
            kappa[1][i] = Math.log(kappa[0][i]);
            //System.out.println("kappa factor: " + (numerator / denominator) + " numerator: " + numerator + " denominator: " + denominator);
        }

        return kappa;

    }

    public static double kappaFac(int numDeps, double hotT, double logRho, double logTemp, double massX, double massZ, double logNH3, double logNH2) {

        double logE = Math.log10(Math.E); // for debug output

        double kapFac = 0.0;

        // These values tuned to produce total kappas of right order of magnitude for Sun
        double constbf = 2.34e19; // b-f pre-factor cm^2/g
        double constff = 3.68e15; // f-f pre-factor cm^2/g
        double constes = 0.2;  // Thomson scattering from free electron pre-factor cm^2/g
        double constHm = 3.9e-31 / 0.02; // H^- b-f pre-factor with 1/0.02 factor from Z term cm^2/g
        //should b-b opacity rho-and T- scaling track b-f oapcity?
        double sigmabf = 1.31e-13;  // Hydrogen b-f x-section, cm^-2
        double refLambda = 500.0; //reference lambda in nm for HI bf opacity formula

        // Paschen continuum H I opacity from n=3:
        double n3 = 3.0;
        double lamJump3 = 820.4; //Paschen jump in nm
        double logHbfFac3 = Math.log(sigmabf) - 5.0 * n3 + 3.0 * (Math.log(lamJump3) - Math.log(refLambda));
        //double hbfFac = Math.pow(lamJump / refLambda, 3.0) / Math.pow(n, 5);
        // Paschen continuum H I opacity from n=3:
        double n2 = 2.0;
        double lamJump2 = 364.0; //Paschen jump in nm
        double logHbfFac2 = Math.log(sigmabf) - 5.0 * n2 + 3.0 * (Math.log(lamJump2) - Math.log(refLambda));

        double logRhoT35, rhoT35;
        double logHmTerm, HmTerm, HmTermHot, HmHotFac;
        double logHIbfTerm3, logHIbfTerm2, HIbfTerm;
        double logHotT = Math.log(hotT);
        double thisTemp = Math.exp(logTemp);

        logRhoT35 = logRho - 3.5 * logTemp;
        rhoT35 = Math.exp(logRhoT35);

        logHmTerm = Math.log(constHm) + Math.log(massZ) + 0.5 * logRho + 9.0 * logTemp; // H^- b-f term
        HmTerm = Math.exp(logHmTerm);
        double midRange = 1500.0;  //H^- opacity ramp-down T range

        if (thisTemp < hotT) {
            // Caroll & Ostlie 2nd Ed. Ch. 9 - (1+X) factors do NOT cancel out when we divide kappa_Star/kappa_Sun
//            // Cool stars: kappa_bf + kappa_ff + kappa_H^- + kappa_es
            kapFac = rhoT35 * (1.0 + massX) * (constbf * massZ + constff * (1.0 - massZ)) + HmTerm + (1.0 + massX) * constes;
            // Cool stars: kappa_ff + kappa_H^- + kappa_es
            //kapFac = rhoT35 * (1.0 + massX) * (constff * (1.0 - massZ)) + HmTerm + (1.0 + massX) * constes;
            //kapFac =  HmTerm + (1.0 + massX) * constes;
            //System.out.println("Cool T: " + Math.exp(logTemp)
            //        + " b-f: " + logE * Math.log(rhoT35 * (1.0 + massX) * (constbf * massZ))
            //        + " f-f: " + logE * Math.log(rhoT35 * (1.0 + massX) * (constff * (1.0 - massZ)))
            //        + " H^-: " + logE * logHmTerm + " es: " + logE * Math.log((1.0 + massX) * constes)
            //        + " kapFac " + kapFac);
        }

        logHIbfTerm3 = logHbfFac3 + logNH3 - logRho;  // cm^2/g //neglects stimualted emission (for now);
        logHIbfTerm2 = logHbfFac2 + logNH2 - logRho;  // cm^2/g //neglects stimualted emission (for now)
        HIbfTerm = Math.exp(logHIbfTerm3) + Math.exp(logHIbfTerm2);

        if ( (thisTemp >= hotT) && (thisTemp < (hotT + midRange)) ) {
            HmHotFac = 1.0 - ((thisTemp - hotT) / midRange);
            HmTermHot = HmTerm * Math.sqrt(HmHotFac);
            //System.out.println("HmHotFac: " + HmHotFac);
            kapFac = rhoT35 * (constbf * massZ + constff * (1.0 - massZ)) + constes + HIbfTerm + HmTermHot;
            //System.out.println("Middle T: " + Math.exp(logTemp) + " b-f: " + rhoT35 * (constbf * massZ)
            //        + " f-f: " + rhoT35 * (constff * (1.0 - massZ))
            //        + " es: " + constes + " HIbf: " + HIbfTerm + " HmTermHot: " + HmTermHot + " kapFac " + kapFac);
        }

        if ( thisTemp >= (hotT + midRange) ) {
            // Caroll & Ostlie 2nd Ed. Ch. 9 - (1+X) factors in every term will cancel out when we divide kappa_Star/kappa_Sun
            // Hot stars: kappa_bf + kappa_ff + kappa_es
            kapFac = rhoT35 * (constbf * massZ + constff * (1.0 - massZ)) + constes + HIbfTerm;
            //System.out.println("Hot T: " + Math.exp(logTemp) + " b-f: " + rhoT35 * (constbf * massZ)
            //        + " f-f: " + rhoT35 * (constff * (1.0 - massZ))
            //       + " es: " + constes + " HIbf: " + HIbfTerm + " kapFac " + kapFac);
        }

        return kapFac;

    }

}
