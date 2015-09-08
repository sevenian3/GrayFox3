/*

 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package grayfox;

import static grayfox.TCorr.expOne;

/**
 *
 * @author Ian
 */
public class MulGrayTCorr {

    public static double[][] mgTCorr(int numDeps, double teff, double[][] tauRos, double[][] temp, double[][] rho, double[][] kappa) {

        // updated temperature structure
        double[][] newTemp = new double[2][numDeps];

        //Teff boundary between early and late-type stars:
        double isCool = 6500.0;

        //Set up multi-gray opacity:
        // lambda break-points and gray levels:
        // No. multi-gray bins = num lambda breakpoints +1
        //double[] grayLams = {30.0, 1.0e6};  //nm //test
        //double[] grayLevel = {1.0};  //test
        // ***  Late type stars, Teff < 9500 K (???):
        //
        double minLambda = 30.0;  //nm
        double maxLambda = 1.0e6;  //nm
        int maxNumBins = 11;
        double[] grayLams = new double[maxNumBins + 1];
        double[] grayLevel = new double[maxNumBins];
        double[] epsilon = new double[maxNumBins];
        //initialize everything first:
        for (int iB = 0; iB < maxNumBins; iB++) {
            grayLams[iB] = maxLambda;
            grayLevel[iB] = 1.0;
            epsilon[iB] = 0.99;
        }
        grayLams[maxNumBins] = maxLambda; //Set final wavelength

        double[][] grayLevelsEpsilons = grayLevEps(maxNumBins, minLambda, maxLambda, teff, isCool);

        //Find actual number of multi-gray bins:
        int numBins = 0; //initialization
        for (int i = 0; i < maxNumBins; i++) {
            if (grayLevelsEpsilons[0][i] < maxLambda) {
                numBins++;
            }
        }
        //Add one more for final lambda:
        //numBins++;
        //System.out.println("numBins: " + numBins);

        /*
         if (teff < isCool) {
         // physically based wavelength break-points and gray levels for Sun from Rutten Fig. 8.6
         // H I Balmer and Lyman jumps for lambda <=3640 A, H^- b-f opacity hump in visible & hole at 1.6 microns, increasing f-f beyond that
         double[] lamSet = {minLambda, 91.1, 158.5, 364.0, 794.3, 1600.0, 3.0e3, 1.0e4, 3.3e4, 1.0e5, 3.3e5, maxLambda}; //nm
         double[] levelSet = {1000.0, 100.0, 5.0, 1.0, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 1000.0};
         //photon *thermal* destruction and creation probability (as opposed to scattering)
         //WARNING:  THese cannot be set exactly = 1.0 or a Math.log() will blow up!!
         double[] epsilonSet = {0.50, 0.50, 0.50, 0.50, 0.50, 0.9, 0.99, 0.99, 0.99, 0.99, 0.99};
         int numBins = levelSet.length;

         for (int iB = 0; iB < numBins; iB++) {
         grayLams[iB] = lamSet[iB] * 1.0e-7;
         grayLevel[iB] = levelSet[iB];
         epsilon[iB] = epsilonSet[iB];
         }
         grayLams[numBins] = lamSet[numBins] * 1.0e-7; //Get final wavelength
         } else {
         // *** Early type stars, Teff > 9500 K (???)
         // It's all about H I b-f (??) What about Thomson scattering (gray)?
         // Lyman, Balmer, Paschen, Brackett jumps
         //What about He I features?
         double[] lamSet = {minLambda, 91.1, 364.0, 820.4, 1458.0, maxLambda};  //nm
         double[] levelSet = {100.0, 10.0, 2.0, 1.0, 1.0};  //???
         double[] epsilonSet = {0.5, 0.6, 0.7, 0.8, 0.5};
         int numBins = levelSet.length;
         for (int iB = 0; iB < numBins; iB++) {
         grayLams[iB] = lamSet[iB] * 1.0e-7;;
         grayLevel[iB] = levelSet[iB];
         epsilon[iB] = epsilonSet[iB];
         }
         grayLams[numBins] = lamSet[numBins] * 1.0e-7; //Get final wavelength
         }

         //Find out how many bins we really have:
         int numBins = 0;  //initialize
         for (int iB = 0; iB < maxNumBins; iB++) {
         if (grayLams[iB] < maxLambda) {
         numBins++;
         }
         }
         */
        for (int iB = 0; iB < numBins; iB++) {
            grayLams[iB] = grayLevelsEpsilons[0][iB];
            grayLevel[iB] = grayLevelsEpsilons[1][iB];
            epsilon[iB] = grayLevelsEpsilons[2][iB];
        }
        grayLams[numBins] = grayLevelsEpsilons[0][numBins]; //Get final wavelength

        //Set overall gray-level - how emissive and absorptive the gas is overall
        // a necessary "fudge" because our kappa values are arbitrary rather than "in situ"
        double graySet = 1.0;

        //double tcDamp = 0.5; // damp the temperature corrections, Delta T, by this *multiplicative* factor
        double tcDamp = 1.0;   // no damping - Lambda iteration is slow rather than oscillatory

        double logE = Math.log10(Math.E); // for debug output

        //double[][] planckBol = MulGrayTCorr.planckBin(numDeps, temp, lamStart, lamStop);
        double[] planckBol = new double[numDeps]; //just for reference - not really needed - ??
        double[] jayBol = new double[numDeps]; //just for reference - not really needed - ??
        double[] dBdTBol = new double[numDeps]; //just for reference - not really needed - ??
        double[] cool = new double[numDeps];  // cooling term in Stromgren equation
        double[] heat = new double[numDeps];  // heating term in Stromgren equation
        double[] corrDenom = new double[numDeps]; //denominator in 1st order temp correction
        //double[] accumB = new double[numDeps]; //accumulator

        //CAUTION: planckBin[2][]: Row 0 is bin-integrated B_lambda; row 1 is bin integrated dB/dT_lambda
        double[][] planckBin = new double[2][numDeps];
        double[] jayBin = new double[numDeps];
        double[] dBdTBin = new double[numDeps];
        double logCool, logHeat, logCorrDenom, logCoolTherm, logCoolScat;

        // initialize accumulators & set overell gray kappa level:
        for (int iTau = 0; iTau < numDeps; iTau++) {

            planckBol[iTau] = 0.0;  //just for reference - not really needed - ??
            jayBol[iTau] = 0.0;  //just for reference - not really needed - ??
            dBdTBol[iTau] = 0.0;  //just for reference - not really needed - ??
            cool[iTau] = 0.0;
            heat[iTau] = 0.0;
            corrDenom[iTau] = 0.0;

            kappa[1][iTau] = kappa[1][iTau] + Math.log(graySet);
            kappa[0][iTau] = Math.exp(kappa[1][iTau]);

        }

        for (int iB = 0; iB < numBins; iB++) {
            //System.out.println("iB: " + iB + " grayLams[iB] " + grayLams[iB]);
            planckBin = planckBinner(numDeps, temp, grayLams[iB], grayLams[iB + 1]);

            // We are lambda-operating on a wavelength integrated B_lambda for each multi-gray bin
            jayBin = jayBinner(numDeps, tauRos, temp, planckBin, grayLevel[iB]);
            //System.out.println("tauRos[1][iTau]   planckBin[0]   planckBin[1]   jayBin");
            for (int iTau = 0; iTau < numDeps; iTau++) {
                //System.out.format("%12.8f   %12.8f   %12.8f   %12.8f%n",
                //        logE * tauRos[1][iTau], logE * Math.log(planckBin[0][iTau]), logE * Math.log(planckBin[1][iTau]), logE * Math.log(jayBin[iTau]));
                //CAUTION: planckBin[2][]: Row 0 is bin-integrated B_lambda; row 1 is bin integrated dB/dT_lambda
                //Net LTE volume cooling rate deltaE = Integ_lam=0^infnty(4*pi*kappa*rho*B_lam)dlam - Integ_lam=0^infnty(4*pi*kappa*rho*J_lam)dlam
                // where Jlam = LambdaO[B_lam] - Rutten Eq. 7.32, 7.33 
                // CAUTION: the 4pi and rho factors cancel out when diving B-J term by dB/dT term 
                planckBol[iTau] = planckBol[iTau] + planckBin[0][iTau];  //just for reference - not really needed - ??
                //logCool = Math.log(grayLevel[iB]) + kappa[1][iTau] + Math.log(planckBin[0][iTau]);  //no scatering
                //cool[iTau] = cool[iTau] + Math.exp(logCool);   //no scattering
                logCoolTherm = Math.log(grayLevel[iB]) + Math.log(epsilon[iB]) + kappa[1][iTau] + Math.log(planckBin[0][iTau]);
                logCoolScat = Math.log(grayLevel[iB]) + Math.log((1.0 - epsilon[iB])) + kappa[1][iTau] + Math.log(jayBin[iTau]);
                cool[iTau] = cool[iTau] + Math.exp(logCoolTherm) + Math.exp(logCoolScat);
                jayBol[iTau] = jayBol[iTau] + jayBin[iTau];  //just for reference - not really needed - ??
                logHeat = Math.log(grayLevel[iB]) + kappa[1][iTau] + Math.log(jayBin[iTau]);
                heat[iTau] = heat[iTau] + Math.exp(logHeat);
                dBdTBol[iTau] = dBdTBol[iTau] + planckBin[1][iTau];  //just for reference - not really needed - ??
                logCorrDenom = Math.log(grayLevel[iB]) + kappa[1][iTau] + Math.log(planckBin[1][iTau]);
                corrDenom[iTau] = corrDenom[iTau] + Math.exp(logCorrDenom);
                // if (iTau == 10){
                //    System.out.format("iB: %02d   %12.8f   %12.8f   %12.8f   %12.8f%n", iB, logE*Math.log(planckBin[0][iTau]), logE*Math.log(cool[iTau]), logE*Math.log(heat[iTau]), logE*Math.log(corrDenom[iTau]));
                //}
            } // iTau loop
        } //iB loop

        // System.out.println("i   tauRos[1][iTau]   planckBol[0]   planckBol[1]   jayBol      cool      heat      corrDenom");
        // for (int iTau = 0; iTau < numDeps; iTau++) {
        //System.out.format("%02d   %12.8f   %12.8f   %12.8f   %12.8f   %12.8f   %12.8f   %12.8f%n", iTau,
        //        logE * tauRos[1][iTau], logE * Math.log(planckBol[iTau]), logE * Math.log(dBdTBol[iTau]), logE * Math.log(jayBol[iTau]),
        //        logE * Math.log(cool[iTau]), logE * Math.log(heat[iTau]), logE * Math.log(corrDenom[iTau]));
        // }
        double logRatio, ratio, deltaTemp, logDeltaTemp;
        double sign = 1.0; //initialize for positive JminusB

        //System.out.println("tauRos[1][iTau]   deltaTemp[iTau]");
        for (int iTau = 0; iTau < numDeps; iTau++) {
            // Compute a 1st order T correction:  Compute J-B so that DeltaT < 0 if J < B:
// avoid direct subtraction of two large almost equal numbers, J & B:

            /* 
             //Gray method:
   
             double JminusB
             logRatio = Math.log(planckBol[iTau]) - Math.log(jayBol[iTau]);
             ratio = Math.exp(logRatio);
             JminusB = jayBol[iTau] * (1.0 - ratio);
             if (JminusB < 0.0) {
             sign = -1.0;
             }

             // DeltaB/DeltaT ~ dB/dT & dB/dT = (4/pi)sigma*T^3
             logDeltaTemp = Math.log(Math.abs(JminusB)) + Math.log(Math.PI) - Math.log(4.0) - Useful.logSigma() - 3.0 * temp[1][iTau];
             deltaTemp[iTau] = sign * Math.exp(logDeltaTemp) * tcDamp;
             //System.out.format("%12.8f   %12.8f%n", tauRos[1][iTau], deltaTemp[iTau]);

             sign = 1.0; //reset sign
             */
            //Multi-Gray method:
            double deltaE;
            //double logHeatNorm, heatNorm, logCoolNorm, deltaENorm;

            ////Normalize numbers by dividing heat and cool terms each by common denominator derivative term first:
            //logHeatNorm = Math.log(heat[iTau]) - Math.log(corrDenom[iTau]);
            //heatNorm = Math.exp(logHeatNorm);
            //logCoolNorm = Math.log(cool[iTau]) - Math.log(corrDenom[iTau]);
            logRatio = Math.log(cool[iTau]) - Math.log(heat[iTau]);
            //logRatio = logCoolNorm - logHeatNorm;

            ratio = Math.exp(logRatio);
            deltaE = heat[iTau] * (1.0 - ratio);
            //deltaENorm = heatNorm * (1.0 - ratio);
            if (deltaE < 0.0) {
                sign = -1.0;
            }
            //CHEAT: Try a Tau-dependent deltaE damping here - things are flaky at tdepth where t(Tau) steepens
            deltaE = deltaE * Math.exp(1.0 * (tauRos[0][0] - tauRos[0][iTau]));

            // DeltaE/DeltaT ~ dB/dT_Bol
            logDeltaTemp = Math.log(Math.abs(deltaE)) - Math.log(corrDenom[iTau]);
            deltaTemp = sign * Math.exp(logDeltaTemp) * tcDamp;
            //deltaTemp = sign * deltaENorm * tcDamp;

            newTemp[0][iTau] = temp[0][iTau] + deltaTemp;
            newTemp[1][iTau] = Math.log(newTemp[0][iTau]);

        } //iTau loop

        return newTemp;

    }

    // method jayBolom computes bolometric angle-averaged mean intensity, J
    // This is a Lambda operation, ie. the Schwartzschild equation
    public static double[] jayBinner(int numDeps, double[][] tauRos, double[][] temp, double[][] planckBin, double grayLevel) {

        // For bolometric J on a Gray Tau scale in LTE: 
        // J(Tau) = 1/2 * Sigma_Tau=0^Infty { E_1(|t-Tau|)*Planck_Bol(Tau) }
        double logE = Math.log10(Math.E); // for debug output

        double E1; //E_1(x)

        //Set up local optical depth scale:
        double[][] tauBin = new double[2][numDeps];
        double deltaTauRos;
        tauBin[0][0] = tauRos[0][0] * grayLevel; // Is this a good idea??
        tauBin[1][0] = Math.log(tauBin[0][0]);
        for (int iTau = 1; iTau < numDeps; iTau++) {
            deltaTauRos = tauRos[0][iTau] - tauRos[0][iTau - 1];
            //grayLevel *is*, by definition, the ratio kappa_Bin/kappaRos that we need here!
            tauBin[0][iTau] = tauBin[0][iTau - 1] + grayLevel * deltaTauRos;
            tauBin[1][iTau] = Math.log(tauBin[0][iTau]);
        }

        double logInteg, integ, integ1, integ2, logInteg1, logInteg2, meanInteg, logMeanInteg, term, logTerm;
        double deltaTau, logDeltaTau; //accumulator
        double accum = 0.0; //accumulator

        double[] jayBin = new double[numDeps];

        // if E_1(t-Tau) evaluated at Tau=bottom of atmosphere, then just set Jay=B at that Tau - we're deep enough to be thermalized
        // and we don't want to let lambda operation implicitly include depths below bottom of model where B=0 implicitly 
        double tiny = 1.0e-14;  //tuned to around maxTauDIff at Tau_Ros ~ 3
        double maxTauDiff;

        //stipulate the {|t-Tau|} grid at which E_1(x)B will be evaluated - necessary to sample the 
        // sharply peaked integrand properly
        // ** CAUTION: minLog10TmTau and maxLog10TmTau are tau offsets from the centre of J integration, 
        //  NOT the optical depth scale of the atmosphere!
        //stipulate the {|t-Tau|} grid at which E_1(x)B will be evaluated - necessary to sample the 
        // sharply peaked integrand properly
        double fineFac = 3.0; // integrate E1 on a grid fineFac x finer in logTau space
        double E1Range = 36.0; // number of master tauBin intervals within which to integrate J 
        double numInteg = E1Range * fineFac;  //
        double deltaLogTauE1 = (tauBin[1][numDeps - 1] - tauBin[1][0]) / numDeps;
        deltaLogTauE1 = deltaLogTauE1 / fineFac;

        double thisTau1, logThisTau1, thisTau2, logThisTau2, logE1, deltaTauE1, logThisPlanck, iFloat;

        //Prepare 1D vectors for Interpol.interpol:
        double[] logTauBin = new double[numDeps];
        double[] logPlanck = new double[numDeps];
        //System.out.println("logTauBin  logB");
        for (int k = 0; k < numDeps; k++) {
            logTauBin[k] = tauBin[1][k];
            logPlanck[k] = Math.log(planckBin[0][k]);
            //System.out.format("%12.8f   %12.8f%n", logE*logTauBin[k], logE*logPlanck[k]);
        }

        //Outer loop over Taus where Jay(Tau) being computed:
        // Start from top and work down to around tau=1 - below that assume we're thermalized with J=B
        //System.out.println("For logTauRos = " + logE*tauRos[1][40] + ": thisTau  E1xB  E1  B");
        //System.out.println("tauRos[1][iTau]   Math.log(planckBin[iTau])   jayBin[1][iTau]");
        for (int iTau = 0; iTau < numDeps; iTau++) {
            //System.out.println("jayBinner: iTau: " + iTau + " tauRos[0] " + tauRos[0][iTau] + " tauRos[1] " + logE * tauRos[1][iTau]);
            jayBin[iTau] = planckBin[0][iTau]; //default initialization J_bin = B_bin

            if (tauRos[0][iTau] < 66.67) {
                //System.out.println("tauRos[0] < limit condition passed");
                // initial test - don't let E_1(x) factor in integrand run off bottom of atmosphere
                // - we have no emissivity down there and J will drop below B again, like at surface!
                maxTauDiff = Math.abs(tauBin[0][numDeps - 1] - tauBin[0][iTau]);
                //System.out.println("tauBin[0][numDeps - 1]: " + tauBin[0][numDeps - 1] + " tauBin[0][iTau] " + tauBin[0][iTau] + " maxTauDiff " + maxTauDiff);
                //System.out.println("maxTauDiff= " + maxTauDiff + " expOne(maxTauDiff)= " + expOne(maxTauDiff));
                if (expOne(maxTauDiff) < tiny) {

                    //System.out.println("maxTauDiff < tiny condition passed, expOne(maxTauDiff): " + expOne(maxTauDiff));
// We're above thermalization depth: J may not = B:
                    //Inner loop over depths contributing to each Jay(iTau):
                    // work outward from t=Tau (ie. i=iTau) piece-wise  
                    accum = 0.0;
                    // conribution from depths above Tau:
                    //initial integrand:
                    // start at i=1 instead of i=0 - cuts out troublesome central cusp of E_1(x) - but now we underestimate J!
                    logThisTau1 = tauBin[1][iTau] - deltaLogTauE1;
                    thisTau1 = Math.exp(logThisTau1);
                    deltaTauE1 = tauBin[0][iTau] - thisTau1;
                    E1 = expOne(deltaTauE1);
                    logE1 = Math.log(E1);
                    logThisPlanck = Interpol.interpol(logTauBin, logPlanck, logThisTau1);
                    logInteg1 = logE1 + logThisPlanck;
                    integ1 = Math.exp(logInteg1);

                    for (int i = 2; i < numInteg - 1; i++) {

                        iFloat = (double) i;
                        // Evaluate E_1(x) and log(E_1(x)) one and for all here

                        //System.out.format("%02d %12.8f %12.8f%n", j, tmTau[j], E1);
                        // LTE bolometric source function is Bolometric Planck function
                        // Extended trapezoidal rule for non-uniform abscissae - this is an exponential integrand!             
                        // We cannot evaluate E_1(x) at x=0 - singular:
                        logThisTau2 = tauBin[1][iTau] - iFloat * deltaLogTauE1;
                        thisTau2 = Math.exp(logThisTau2);

                        //if (i == numInteg - 2) {
                        //    System.out.println("i " + i + " logThisTau1 " + logE * logThisTau1 + " logThisTau2 " + logE * logThisTau2);
                        //}
                        // Make sure we're still in the atmosphere!
                        if (logThisTau2 > tauBin[1][0]) {
                            //if (i == numInteg - 2) {
                            //    System.out.println("thisTau2 > tauBin[0][0] condition passed");
                            //}
                            //if (iTau == 37) {
                            //    System.out.println("i " + i + " logThisTau1 " + logE * logThisTau1 + " logThisTau2 " + logE * logThisTau2);
                            //}

                            deltaTauE1 = tauBin[0][iTau] - thisTau2;
                            E1 = expOne(deltaTauE1);
                            logE1 = Math.log(E1);
                            // interpolate log(B(log(Tau)) to the integration abscissa
                            logThisPlanck = Interpol.interpol(logTauBin, logPlanck, logThisTau2);
                            logInteg2 = logE1 + logThisPlanck;
                            integ2 = Math.exp(logInteg2);

                            logDeltaTau = Math.log(thisTau1 - thisTau2); // logDeltaTau *NOT* the same as deltaLogTau!!

                            meanInteg = 0.5 * (integ1 + integ2); //Trapezoid rule
                            logMeanInteg = Math.log(meanInteg);
                            //if (iTau == 40) {
                            //    System.out.format("%15.8f    %15.8f    %15.8f   %15.8f%n", logE*Math.log(thisTau1), logE*logMeanInteg, logE*logE1, logE*logThisPlanck);
                            //}

                            logTerm = logMeanInteg + logDeltaTau;
                            term = Math.exp(logTerm);
                            accum = accum + term;

                            integ1 = integ2;
                            thisTau1 = thisTau2;
                            //if (iTau == 41){
                            //    System.out.println("term " + term + " accum " + accum);
                            //}
                        } // thisTau > 0
                    } // i ("t") loop, above iTau 

                    jayBin[iTau] = 0.5 * accum;  //store what we have.
                    //test jayBin[iTau] = 0.5 * planckBin[0][iTau]; // fake upper half with isotropic result
                    //test jayBin[iTau] = jayBin[iTau] + 0.5 * planckBin[0][iTau]; // test upper atmosphere part of J integration by fixing lower part with isotropic result

                    // conribution from depths below Tau:
                    // include iTau itself so we don't miss the area under the central peak of E_1(x) - the expOne function
                    // will protect itself from the x=0 singularity using variable 'tiny'
                    accum = 0.0;
                    //initial integrand:
                    // start at i=1 instead of i=0 - cuts out troublesome central cusp of E_1(x) - but now we underestimate J!
                    logThisTau1 = tauBin[1][iTau] + deltaLogTauE1;
                    thisTau1 = Math.exp(logThisTau1);
                    deltaTauE1 = thisTau1 - tauBin[0][iTau];
                    E1 = expOne(deltaTauE1);
                    logE1 = Math.log(E1);
                    logThisPlanck = Interpol.interpol(logTauBin, logPlanck, logThisTau1);
                    logInteg1 = logE1 + logThisPlanck;
                    integ1 = Math.exp(logInteg1);

                    for (int i = 2; i < numInteg - 1; i++) {

                        iFloat = (double) i;

                        logThisTau2 = tauBin[1][iTau] + iFloat * deltaLogTauE1;
                        thisTau2 = Math.exp(logThisTau2);
                     // We cannot evaluate E_1(x) at x=0 - singular:
                        // Extended trapezoidal rule for non-uniform abscissae - the is an exponential integrand! 

                        // make sure we're still in the atmosphere!
                        if (logThisTau2 < tauBin[1][numDeps - 1]) {

                            deltaTauE1 = thisTau2 - tauBin[0][iTau];
                            E1 = expOne(deltaTauE1);
                            logE1 = Math.log(E1);

                            logThisPlanck = Interpol.interpol(logTauBin, logPlanck, logThisTau2);
                            logInteg2 = logE1 + logThisPlanck;
                            integ2 = Math.exp(logInteg2);

                            logDeltaTau = Math.log(thisTau2 - thisTau1); // logDeltaTau *NOT* the same as deltaLogTau!!

                            meanInteg = 0.5 * (integ1 + integ2); //Trapezoid rule
                            logMeanInteg = Math.log(meanInteg);
                     //if (iTau == 40) {
                            //    System.out.format("%15.8f    %15.8f    %15.8f    %15.8f%n", logE*Math.log(thisTau1), logE*logMeanInteg, logE*logE1, logE*logThisPlanck);
                            //}

                            // LTE bolometric source function is Bolometric Plnack function
                            logTerm = logMeanInteg + logDeltaTau;
                            term = Math.exp(logTerm);
                            accum = accum + term;

                            integ1 = integ2;
                            thisTau1 = thisTau2;

                        }// if thisTau < tauBin[0][numDeps-1]
                    } // i ("t") loop, below iTau

                    jayBin[iTau] = jayBin[iTau] + 0.5 * accum;

                } //if branch for E_1(x) safely dwindling away before reaching bottom of atmosphere
            } // if branch for above thermalization depth of Tau=10? 

            //System.out.format("%12.8f   %12.8f  %12.8f%n",
            //       logE * tauRos[1][iTau], Math.log10(planckBin[iTau]), Math.log10(jayBin[iTau]));
        } //iTau loop

        return jayBin;

    }

    // Compute linear wave-bin-specific lambda-integrated Planck fn AND it's T derivative at all depths:
    // Row 0: B_bin(tau);   Row 1: dB/dT_bin(tau);
    public static double[][] planckBinner(int numDeps, double[][] temp, double lamStart, double lamStop) {

        double[][] planckBin = new double[2][numDeps];

        double logE = Math.log10(Math.E); // for debug output

        //MultiGray-ready:
        // Parameters of overall lambda grid (nm):
        // Planck.planck() will convert nm to cm - not any more!
        //double log10LamStart = 1.5 * 1.0e-7;  //must be < first Gray lambda break point
        //double log10LamStop = 5.0 * 1.0e-7;   //must be > last Gray lambda break point 
        double log10LamStart = Math.log10(lamStart);
        double log10LamStop = Math.log10(lamStop);
        double deltaLog10Lam = 0.1;

        int numLamAll;
        numLamAll = (int) ((log10LamStop - log10LamStart) / deltaLog10Lam);
        //System.out.println("lamStart " + lamStart + " log10LamStart " + log10LamStart + " lamStop " + lamStop + " log10LamStop " + log10LamStop + " numLamAll " + numLamAll);
        double[] lambda = new double[numLamAll];

        //Generate lambda grid separately to avoid duplicate lambda generation
        double iFloat, thisLogLam;
        //System.out.println("lambdas");
        for (int i = 0; i < numLamAll; i++) {

            iFloat = (double) i;
            thisLogLam = log10LamStart + iFloat * deltaLog10Lam;
            lambda[i] = Math.pow(10.0, thisLogLam);
            //System.out.format("%02d   %12.8f%n", i, lambda[i]);

        }

        double thisLam1, thisLam2, deltaLam, planck1, planck2, logPlanck1, logPlanck2;
        double term, integ, accum;
        double dBdT1, dBdT2, logdBdT1, logdBdT2, accum2;

        //trapezoid rule integration
        //System.out.println("Trapezoid: ");
        for (int iTau = 0; iTau < numDeps; iTau++) {
            //reset accumulators for new depth
            accum = 0.0;
            accum2 = 0.0;

            //initial integrands:
            logPlanck1 = Planck.planck(temp[0][iTau], lambda[0]);
            planck1 = Math.exp(logPlanck1);
            logdBdT1 = Planck.dBdT(temp[0][iTau], lambda[0]);
            dBdT1 = Math.exp(logdBdT1);
            for (int i = 1; i < numLamAll - 1; i++) {

                deltaLam = lambda[i + 1] - lambda[i];
                //deltaLam = deltaLam * 1.0e-7;  //nm to cm

                //Planck.planck returns log(B_lambda)
                logPlanck2 = Planck.planck(temp[0][iTau], lambda[i]);

                planck2 = Math.exp(logPlanck2);

                //if (i == 20) {
                //    System.out.println("lambda " + thisLam1 + " temp[0][iTau] " + temp[0][iTau] + " logPlanck1 " + logE*logPlanck1);
                //}
                //trapezoid rule integration
                integ = 0.5 * (planck1 + planck2) * deltaLam;
                accum = accum + integ;

                planck1 = planck2;

                //Now do the same for dB/dT:
                //Planck.dBdT returns log(dB/dT_lambda)
                logdBdT2 = Planck.dBdT(temp[0][iTau], lambda[i]);

                dBdT2 = Math.exp(logdBdT2);

                //trapezoid rule integration
                integ = 0.5 * (dBdT1 + dBdT2) * deltaLam;
                accum2 = accum2 + integ;

                dBdT1 = dBdT2;

            } // lambda i loop
            planckBin[0][iTau] = accum;
            planckBin[1][iTau] = accum2;
            //System.out.format("%02d   %12.8f%n", iTau, planckBin[iTau]);

        } //iTau loop

        //// Gray only:
        ////if (lamStart == 1000.0) {  //Could be for any gray wavelength
        //double[][] planckBol = new double[2][numDeps];
        //double[][] dBdTBol = new double[2][numDeps];
        //System.out.println("Stefan-Boltzmann:  tauRos[1]  B_Bol   dBdT_Bol");
        //for (int i = 0; i < numDeps; i++) {
        //    planckBol[1][i] = Useful.logSigma() + 4.0 * temp[1][i] - Math.log(Math.PI);
        //    planckBol[0][i] = Math.exp(planckBol[1][i]);
        //    dBdTBol[1][i] = Math.log(4.0) + Useful.logSigma() + 3.0 * temp[1][i] - Math.log(Math.PI);
        //    dBdTBol[0][i] = Math.exp(dBdTBol[1][i]);
        //    System.out.format("%02d   %12.8f   %12.8f%n", i, logE * planckBol[1][i], logE * dBdTBol[1][i]);
        //}
        //}
        return planckBin;
    }

    public static double[][] grayLevEps(int maxNumBins, double minLambda, double maxLambda, double teff, double isCool) {

        //double minLambda = 30.0;  //nm
        //double maxLambda = 1.0e6;  //nm
        //int maxNumBins = 11;
        double[][] grayLevelsEpsilons = new double[3][maxNumBins + 1];
                // The returned structure:
        //Row 0 is wavelength breakpoints
        //Row 1 is relative opacity gray levels
        //Row 2 is absolute thermal photon creation fractions, epsilon

        //initialize everything first:
        for (int iB = 0; iB < maxNumBins; iB++) {
            grayLevelsEpsilons[0][iB] = maxLambda;
            grayLevelsEpsilons[1][iB] = 1.0;
            grayLevelsEpsilons[2][iB] = 0.99;
        }
        grayLevelsEpsilons[0][maxNumBins] = maxLambda; //Set final wavelength

        if (teff < isCool) {
            // physically based wavelength break-points and gray levels for Sun from Rutten Fig. 8.6
            // H I Balmer, Lyman, and Paschen jumps for lambda <=3640 A, H^- b-f opacity hump in visible & hole at 1.6 microns, increasing f-f beyond that
            double[] lamSet = {minLambda, 91.1, 158.5, 364.0, 820.4, 1600.0, 3.0e3, 1.0e4, 3.3e4, 1.0e5, 3.3e5, maxLambda}; //nm
            //double[] levelSet =       {1000.0,100.0, 5.0,   0.5,   0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 1000.0};
            double[] levelSet =       {1000.0,100.0, 5.0,   1.0,   0.5, 0.1, 3.0, 10.0, 30.0, 100.0, 1000.0};
            //photon *thermal* destruction and creation probability (as opposed to scattering)
            //WARNING:  THese cannot be set exactly = 1.0 or a Math.log() will blow up!!
            //double[] epsilonSet =      {0.50, 0.50,  0.50,  0.50,  0.50, 0.9, 0.99, 0.99, 0.99, 0.99, 0.99};
            double[] epsilonSet =      {0.50, 0.50,  0.90,  0.99,  0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99};
            int numBins = levelSet.length;

            for (int iB = 0; iB < numBins; iB++) {
                grayLevelsEpsilons[0][iB] = lamSet[iB] * 1.0e-7;
                grayLevelsEpsilons[1][iB] = levelSet[iB];
                grayLevelsEpsilons[2][iB] = epsilonSet[iB];
            }
            grayLevelsEpsilons[0][numBins] = lamSet[numBins] * 1.0e-7; //Get final wavelength
        } else {
            // *** Early type stars, Teff > 9500 K (???)
            // It's all about H I b-f (??) What about Thomson scattering (gray)?
            // Lyman, Balmer, Paschen, Brackett jumps
            //What about He I features?
            double[] lamSet = {minLambda, 91.1, 364.0, 820.4, 1458.0, maxLambda};  //nm
            double[] levelSet = {100.0, 10.0, 2.0, 1.0, 1.0};  //???
            double[] epsilonSet = {0.5, 0.6, 0.7, 0.8, 0.5};
            int numBins = levelSet.length;
            for (int iB = 0; iB < numBins; iB++) {
                grayLevelsEpsilons[0][iB] = lamSet[iB] * 1.0e-7;;  //cm
                grayLevelsEpsilons[1][iB] = levelSet[iB];
                grayLevelsEpsilons[2][iB] = epsilonSet[iB];
            }
            grayLevelsEpsilons[0][numBins] = lamSet[numBins] * 1.0e-7; //Get final wavelength
        }

        return grayLevelsEpsilons;

    }

}
