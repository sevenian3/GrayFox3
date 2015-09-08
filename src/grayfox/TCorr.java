/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package grayfox;

/**
 * To Do: June 5: Interpolate B(Tau) onto Tau grid that's *linear* in tau (not
 * logTau) and more finely spaced (500 depths?) and use higher order integration
 * method for uniformly spaced abscissae (Simpson's 3/8 is cubic)
 *
 * // Try a temperature correction step by re-calculating radiative equilibrium
 * (RE) with angle-averaged mean // intensity, J, computed from B with the Milne
 * equation (a Lambda iteration):
 *
 * @author Ian
 */
public class TCorr {

    public static double[][] tCorr(int numDeps, double[][] tauRos, double[][] temp) {

        // updated temperature structure
        double[][] newTemp = new double[2][numDeps];

        //double tcDamp = 0.5; // damp the temperature corrections, Delta T, by this *multiplicative* factor
        double tcDamp = 1.0;   // no damping - Lambda iteration is slow rather than oscillatory

        double logE = Math.log10(Math.E); // for debug output

        double[][] planckBol = planckBolom(numDeps, temp);
        double[][] jayBol = jayBolom(numDeps, tauRos, temp);

        double logRatio, ratio, JminusB, deltaTemp, logDeltaTemp;
        double sign = 1.0; //initialize for positive JminusB

        //System.out.println("tauRos[1][iTau]   deltaTemp[iTau]");
        for (int iTau = 0; iTau < numDeps; iTau++) {
            // Compute a 1st order T correction:  Compute J-B so that DeltaT < 0 if J < B:
// avoid direct subtraction of two large almost equal numbers, J & B:

            logRatio = planckBol[1][iTau] - jayBol[1][iTau];
            ratio = Math.exp(logRatio);
            JminusB = jayBol[0][iTau] * (1.0 - ratio);
            if (JminusB < 0.0) {
                sign = -1.0;
            }

            // DeltaB/DeltaT ~ dB/dT & dB/dT = (4/pi)sigma*T^3
            logDeltaTemp = Math.log(Math.abs(JminusB)) + Math.log(Math.PI) - Math.log(4.0) - Useful.logSigma() - 3.0 * temp[1][iTau];
            deltaTemp = sign * Math.exp(logDeltaTemp) * tcDamp;
            //System.out.format("%12.8f   %12.8f%n", tauRos[1][iTau], deltaTemp);

            sign = 1.0; //reset sign           

            newTemp[0][iTau] = temp[0][iTau] + deltaTemp;
            newTemp[1][iTau] = Math.log(newTemp[0][iTau]);

        } //iTau loop

        return newTemp;

    }

    // method jayBolom computes bolometric angle-averaged mean intensity, J, for Gray LTE case:
    public static double[][] jayBolom(int numDeps, double[][] tauRos, double[][] temp) {

        // For bolometric J on a Gray Tau scale in LTE: 
        // J(Tau) = 1/2 * Sigma_Tau=0^Infty { E_1(|t-Tau|)*Planck_Bol(Tau) }
        double logE = Math.log10(Math.E); // for debug output

        double E1; //E_1(x)

        double logInteg, integ, integ1, integ2, logInteg1, logInteg2, meanInteg, logMeanInteg, term, logTerm;
        double deltaTau, logDeltaTau; //accumulator
        double accum = 0.0; //accumulator

        double[][] jayBol = new double[2][numDeps];

        double[][] planckBol = planckBolom(numDeps, temp);

        // // this is currently getting clobbered
        // // Top of atmosphere: Jay = 0.5S = 0.5B - - LTE Eddington-Barbier result for J
        // jayBol[0][0] = 0.5 * planckBol[0][0];
        // jayBol[1][0] = Math.log(jayBol[0][0]);
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
        double E1Range = 36.0; // number of master tauRos intervals within which to integrate J 
        double numInteg = E1Range * fineFac;  //
        double deltaLogTauE1 = (tauRos[1][numDeps - 1] - tauRos[1][0]) / numDeps;
        deltaLogTauE1 = deltaLogTauE1 / fineFac;

        double thisTau1, logThisTau1, thisTau2, logThisTau2, logE1, deltaTauE1, logThisPlanck, iFloat;

        //Prepare 1D vectors for Interpol.interpol:
        double[] logTauRos = new double[numDeps];
        double[] logPlanck = new double[numDeps];
        //System.out.println("logTauRos  logB");
        for (int k = 0; k < numDeps; k++) {
            logTauRos[k] = tauRos[1][k];
            logPlanck[k] = planckBol[1][k];
            //System.out.format("%12.8f   %12.8f%n", logE*logTauRos[k], logE*logPlanck[k]);
        }

        //Outer loop over Taus where Jay(Tau) being computed:
        // Start from top and work down to around tau=1 - below that assume we're thermalized with J=B
        //System.out.println("For logTauRos = " + logE*tauRos[1][40] + ": thisTau  E1xB  E1  B");
        //System.out.println("tauRos[1][iTau]   planckBol[1][iTau]   jayBol[1][iTau]");
        for (int iTau = 0; iTau < numDeps; iTau++) {
            jayBol[0][iTau] = planckBol[0][iTau]; //default initialization J=B
            jayBol[1][iTau] = planckBol[1][iTau]; //default initialization logJ = logB

            if (tauRos[0][iTau] < 0.333) {

                // initial test - don't let E_1(x) factor in integrand run off bottom of atmosphere
                // - we have no emissivity down there and J will drop below B again, like at surface!
                maxTauDiff = Math.abs(tauRos[0][numDeps - 1] - tauRos[0][iTau]);
                //System.out.println("maxTauDiff= " + maxTauDiff + " expOne(maxTauDiff)= " + expOne(maxTauDiff));
                if (expOne(maxTauDiff) < tiny) {

                    //System.out.println("Lambda operation triggered, iTau = " + iTau + " tauRos[0][iTau] " + tauRos[0][iTau]);
// We're above thermalization depth: J may not = B:
                    //Inner loop over depths contributing to each Jay(iTau):
                    // work outward from t=Tau (ie. i=iTau) piece-wise  
                    accum = 0.0;
                    // conribution from depths above Tau:

                    // start at i=1 instead of i=0 - cuts out troublesome central cusp of E_1(x) - but now we underestimate J!
                    //initial integrand:
                    logThisTau1 = tauRos[1][iTau] - deltaLogTauE1;
                    thisTau1 = Math.exp(logThisTau1);
                    deltaTauE1 = tauRos[0][iTau] - thisTau1;
                    E1 = expOne(deltaTauE1);
                    logE1 = Math.log(E1);
                    logThisPlanck = Interpol.interpol(logTauRos, logPlanck, logThisTau1);
                    logInteg1 = logE1 + logThisPlanck;
                    integ1 = Math.exp(logInteg1);

                    for (int i = 2; i < numInteg - 1; i++) {

                        iFloat = (double) i;
                        // Evaluate E_1(x) and log(E_1(x)) one and for all here

                        //System.out.format("%02d %12.8f %12.8f%n", j, tmTau[j], E1);
                        // LTE bolometric source function is Bolometric Planck function
                        // Extended trapezoidal rule for non-uniform abscissae - this is an exponential integrand!             
                        // We cannot evaluate E_1(x) at x=0 - singular:
                        logThisTau2 = tauRos[1][iTau] - iFloat * deltaLogTauE1;
                        thisTau2 = Math.exp(logThisTau2);

                        // Make sure we're still in the atmosphere!
                        if (logThisTau2 > tauRos[1][0]) {
                            //if (iTau == 20) {
                            //    System.out.println("We're in the atmosphere!");
                            //}

                            deltaTauE1 = tauRos[0][iTau] - thisTau2;
                            E1 = expOne(deltaTauE1);
                            logE1 = Math.log(E1);
                            // interpolate log(B(log(Tau)) to the integration abscissa
                            logThisPlanck = Interpol.interpol(logTauRos, logPlanck, logThisTau2);
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
                        } // thisTau > 0

                    } // i ("t") loop, above iTau 
                    jayBol[0][iTau] = 0.5 * accum;  //store what we have.

                    // conribution from depths below Tau:
                    // include iTau itself so we don't miss the area under the central peak of E_1(x) - the expOne function
                    // will protect itself from the x=0 singularity using variable 'tiny'
                    accum = 0.0;

                    //initial integrand:
                    // start at i=1 instead of i=0 - cuts out troublesome central cusp of E_1(x) - but now we underestimate J!
                    logThisTau1 = tauRos[1][iTau] + deltaLogTauE1;
                    thisTau1 = Math.exp(logThisTau1);
                    deltaTauE1 = thisTau1 - tauRos[0][iTau];
                    E1 = expOne(deltaTauE1);
                    logE1 = Math.log(E1);
                    logThisPlanck = Interpol.interpol(logTauRos, logPlanck, logThisTau1);
                    logInteg1 = logE1 + logThisPlanck;
                    integ1 = Math.exp(logInteg1);
                    
                    for (int i = 2; i < numInteg - 1; i++) {

                        iFloat = (double) i;

                        logThisTau2 = tauRos[1][iTau] + iFloat * deltaLogTauE1;
                        thisTau2 = Math.exp(logThisTau2);
                        // We cannot evaluate E_1(x) at x=0 - singular:
                        // Extended trapezoidal rule for non-uniform abscissae - the is an exponential integrand! 

                        // make sure we're still in the atmosphere!
                        if (logThisTau2 < tauRos[1][numDeps - 1]) {

                            deltaTauE1 = thisTau2 - tauRos[0][iTau];
                            E1 = expOne(deltaTauE1);
                            logE1 = Math.log(E1);

                            logThisPlanck = Interpol.interpol(logTauRos, logPlanck, logThisTau2);
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

                        }// if thisTau < tauRos[0][numDeps-1]
                    } // i ("t") loop, below iTau

                    jayBol[0][iTau] = jayBol[0][iTau] + 0.5 * accum;
                    jayBol[1][iTau] = Math.log(jayBol[0][iTau]);

                } //if branch for E_1(x) safely dwindling away before reaching bottom of atmosphere
            } // if branch for above thermalization depth of Tau=10? 

            //System.out.format("%12.8f   %12.8f  %12.8f%n",
            //       logE * tauRos[1][iTau], logE * planckBol[1][iTau], logE * jayBol[1][iTau]);
        } //iTau loop

        return jayBol;

    }

    // Compute linear and logarithmic bolometric Planck fn at all depths one and for all:
    public static double[][] planckBolom(int numDeps, double[][] temp) {

        double[][] planckBol = new double[2][numDeps];

        for (int i = 0; i < numDeps; i++) {
            planckBol[1][i] = Useful.logSigma() + 4.0 * temp[1][i] - Math.log(Math.PI);
            planckBol[0][i] = Math.exp(planckBol[1][i]);
        }

        return planckBol;
    }

    // Approximate first exponential integral function E_1(x) = -Ei(-x)
    public static double expOne(double x) {

        // From http://en.wikipedia.org/wiki/Exponential_integral 
        // Series expansion for first exponential integral function, E_1(x) = -Ei(-x)
        // Ee_one(x) = -gamma - ln(abs(x)) - Sigma_k=1^infnty{(-x)^k)/(k*k!)}
        // where: gamma =  Euler–Mascheroni constant = 0.577215665...
        double E1;

        x = Math.abs(x); // x must be positive
        // E1(x) undefined at x=0 - singular:
        //double tiny = 1.25;  //tuned to give J ~ 0.5B @ tau=0
        double tiny = 1.0e-6;
        if (x < tiny) {
            x = tiny;
        }

        // Caution: even at 11th order acuracy (k=11), approximation starts to diverge for x . 3.0:
        if (x > 3.0) {

            E1 = Math.exp(-1.0 * x) / x; // large x approx

        } else {
            double gamma = 0.577215665; //Euler–Mascheroni constant
            double kTerm = 0.0;
            double order = 11; //order of approximation
            double kFloat;
            double accum = 0.0;  //accumulator
            double kFac = 1.0; // initialize k! (k factorial)

            for (int k = 1; k <= order; k++) {
                kFloat = (double) k;
                kFac = kFac * kFloat;
                accum = accum + Math.pow((-1.0 * x), kFloat) / (k * kFac);
                //System.out.println("k: " + k + " kFac: " + kFac);
                //System.out.println("k: " + k + " Math.pow(x, kFloat): " + Math.pow(x, kFloat));
            }
            kTerm = accum;

            E1 = -1.0 * gamma - Math.log(Math.abs(x)) - kTerm;
        }

        //System.out.println("x: " + x + " exp1(x): " + E1);
        return E1;

    }

}
