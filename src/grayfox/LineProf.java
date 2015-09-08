
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
 * Line profile, phi_lambda(lambda): Assume Voigt function profile - need H(a,v)
 * Assumes CRD, LTE, ??? Input parameters: lam0 - line center wavelength in nm
 * mass - mass of absorbing particle (amu) logGammaCol - log_10(gamma) - base 10
 * logarithmic collisional (pressure) damping co-efficient (s^-1) epsilon -
 * convective microturbulence- non-thermal broadening parameter (km/s) Also
 * needs atmospheric structure information: numDeps WON'T WORK - need observer's
 * frame fixed lambda at all depths: temp structure for depth-dependent thermal
 * line broadening Teff as typical temp instead of above pressure structure,
 * pGas, if scaling gamma
 */
public class LineProf {

    public static double[][] voigt(double[][] linePoints, double lam0In, double logGammaCol,
            int numDeps, double teff, double[][] tauRos, double[][] temp, double[][] press,
            double[][] tempSun, double[][] pressSun) {

        double c = Useful.c;
        double logC = Useful.logC();
        //double k = Useful.k;
        double logK = Useful.logK();
        //double e = Useful.e;
        //double mE = Useful.mE;

        double lam0 = lam0In; // * 1.0E-7; //nm to cm
        double logLam0 = Math.log(lam0);

        double ln10 = Math.log(10.0);
        double ln2 = Math.log(2.0);
        double ln4pi = Math.log(4.0 * Math.PI);
        double lnSqRtPi = 0.5 * Math.log(Math.PI);
        double sqPi = Math.sqrt(Math.PI);
        //double ln100 = 2.0*Math.log(10.0);

        double logE = Math.log10(Math.E); // for debug output

        double doppler = linePoints[0][1] / linePoints[1][1];
        double logDopp = Math.log(doppler);
        //System.out.println("LineProf: doppler, logDopp: " + doppler + " " + logE*logDopp);

        //Put input parameters into linear cgs units:
        //double gammaCol = Math.pow(10.0, logGammaCol);
        // Lorentzian broadening:
        // Assumes Van der Waals dominates radiative damping
        // log_10 Gamma_6 for van der Waals damping around Tau_Cont = 1 in Sun 
        //  - p. 57 of Radiative Transfer in Stellar Atmospheres (Rutten)
        double logGammaSun = 9.0 * ln10; // Convert to base e 
        //double logFudge = Math.log(2.5);  // Van der Waals enhancement factor

        int tau1 = TauPoint.tauPoint(numDeps, tauRos, 1.0);

        //System.out.println("LINEGRID: Tau1: " + tau1);
        //logA = 2.0 * logLam0 + logGamma - ln4pi - logC - logDopp;
        //a = Math.exp(logA);
        //System.out.println("LINEGRID: logA: " + logE * logA);
        //Set up a half-profile Delta_lambda grid in Doppler width units 
        //from line centre to wing
        int numPoints = linePoints[0].length;
        //System.out.println("LineProf: numPoints: " + numPoints);

        // Return a 2D numPoints X numDeps array of normalized line profile points (phi)
        double[][] lineProf = new double[numPoints][numDeps];

        // Line profiel points in Doppler widths - needed for Voigt function, H(a,v):
        double[] v = new double[numPoints];
        double logV, ii;

//        lineProf[0][0] = 0.0; v[0] = 0.0; //Line centre - cannot do logaritmically!
        double gamma, logGamma, a, logA, voigt, core, wing, logWing, logVoigt;
        int il0 = 36;
        //System.out.println("il0 " + il0 + " temp[il] " + temp[0][il0] + " press[il] " + logE*press[1][il0]);
        for (int id = 0; id < numDeps; id++) {

            //Formula from p. 56 of Radiative Transfer in Stellar Atmospheres (Rutten),
            // logarithmically with respect to solar value:
            logGamma = press[1][id] - pressSun[1][tau1] + 0.7 * (tempSun[1][tau1] - temp[1][id]) + logGammaSun;
            //logGamma = logGamma + logFudge + logGammaCol;
            logGamma = logGamma + logGammaCol;
            //System.out.println("LineGrid: logGamma: " + id + " " + logE * logGamma);

            //Voigt "a" parameter with line centre wavelength:
            logA = 2.0 * logLam0 + logGamma - ln4pi - logC - logDopp;
            a = Math.exp(logA);
            //System.out.println("LineGrid: logGamma: " + id + " " + logE * logGamma + " " + logE * logA);
            //if (id == 30) {
            //    //System.out.println("il   v[il]   iy   y   logNumerator   logDenominator   logInteg ");
            //    System.out.println("voigt:   v   logVoigt: ");
            //}
            for (int il = 0; il < numPoints; il++) {

                v[il] = linePoints[1][il];
                //System.out.println("LineProf: il, v[il]: " + il + " " + v[il]);

                //if (il <= numCore) {
                if (v[il] <= 2.0 && v[il] >= -2.0) {

                    // - Gaussian ONLY - at line centre Lorentzian will diverge!
                    core = Math.exp(-1.0 * (v[il] * v[il]));
                    voigt = core;
                    //System.out.println("LINEGRID- CORE: core: " + core);

                } else {

                    logV = Math.log(Math.abs(v[il]));

                    //Gaussian core:
                    core = Math.exp(-1.0 * (v[il] * v[il]));
                    //System.out.println("LINEGRID- WING: core: " + core);

                    //Lorentzian wing:
                    logWing = logA - lnSqRtPi - (2.0 * logV);
                    wing = Math.exp(logWing);

                    voigt = core + wing;
                    //System.out.println("LINEGRID- WING: wing: " + wing + " logV " + logV);

                } // end else

                //System.out.println("LINEGRID: il, v[il]: " + il + " " + v[il] + " lineProf[0][il]: " + lineProf[0][il]);
                //System.out.println("LINEGRID: il, Voigt, H(): " + il + " " + voigt);
                //Convert from H(a,v) in dimensionless Voigt units to physical phi((Delta lambda) profile:
                logVoigt = Math.log(voigt) + 2.0 * logLam0 - lnSqRtPi - logDopp - logC;

                lineProf[il][id] = Math.exp(logVoigt);
                //if (id == 36) {
                //    System.out.println("il " + il + " linePoints " + 1.0e7 * linePoints[0][il] + " id " + id + " lineProf[il][id] " + lineProf[il][id]);
                //}

                //System.out.println("LineProf: il, id, lineProf[il][id]: " + il + " " + id + " " + lineProf[il][id]);
            } // il lambda loop

            // if (id == 20) {
            //     for (int il = 0; il < numPoints; il++) {
            //        System.out.format("Voigt: %20.16f   %20.16f%n", linePoints[1][il], logE * Math.log(lineProf[il][id]));
            //    }
            // }
        } //id loop

        /*  Line doubling now done in LineGrid: No longer done here:
         // Add the negative DeltaLambda half of the line:
         int numPoints2 = (2 * numPoints) - 1;
         //System.out.println("LineGrid: numpoints2: " + numPoints2);

         // Return a 2XnumPoints vector of Delta Lambdas and area-normalized line profile points (phi)
         // Row 0 : Delta lambdas in cm - will need to be in ? for Planck and Rad Trans?
         // Row 1 : Area-normalized line profile ("phi")
         double[][] lineProf2 = new double[numPoints2][numDeps];

         //wavelengths are depth-independent - just put them in the 0th depth slot:
         for (int il2 = 0; il2 < numPoints2; il2++) {

         if (il2 < numPoints - 1) {

         int il = (numPoints - 1) - il2;
         lineProf2[0][il2][0] = -1.0 * lineProf[0][il][0];

         } else {

         //Positive DelataLambda half:   
         int il = il2 - (numPoints - 1);
         lineProf2[0][il2][0] = lineProf[0][il][0];

         }

         for (int id = 0; id < numDeps; id++) {

         //Negative DelataLambda half:
         if (il2 < numPoints - 1) {

         int il = (numPoints - 1) - il2;
         lineProf2[1][il2][id] = lineProf[1][il][id];

         } else {

         //Positive DelataLambda half:   
         int il = il2 - (numPoints - 1);
         lineProf2[1][il2][id] = lineProf[1][il][id];

         }

         } //id loop

         } //il2 loop
         */

        /* Debug
         // Check that line profile is area-normalized (it is NOT, area = 1.4845992503443734E-19!, but IS constant with depth - !?:
         double delta;
         for (int id = 0; id < numDeps; id++) {
         double sum = 0.0;
         for (int il = 1; il < numPoints2; il++) {
         delta = lineProf2[0][il][id] - lineProf2[0][il - 1][id];
         sum = sum + (lineProf2[1][il][id] * delta);
         }
         System.out.println("LineGrid: id, Profile area = " + id + " " + sum );
         }
         */
        return lineProf;

    }

    public static double[][] voigt2(double[][] linePoints, double lam0In, double logGammaCol,
            int numDeps, double teff, double[][] tauRos, double[][] temp, double[][] press,
            double[][] tempSun, double[][] pressSun) {

        double c = Useful.c;
        double logC = Useful.logC();
        //double k = Useful.k;
        double logK = Useful.logK();
        //double e = Useful.e;
        //double mE = Useful.mE;

        double lam0 = lam0In; // * 1.0E-7; //nm to cm
        double logLam0 = Math.log(lam0);

        double ln10 = Math.log(10.0);
        double ln2 = Math.log(2.0);
        double ln4pi = Math.log(4.0 * Math.PI);
        double lnSqRtPi = 0.5 * Math.log(Math.PI);
        double sqPi = Math.sqrt(Math.PI);
        //double ln100 = 2.0*Math.log(10.0);

        double logE = Math.log10(Math.E); // for debug output

        //Recover Doppler width in cm:
        double doppler = linePoints[0][1] / linePoints[1][1];
        double logDopp = Math.log(doppler);
        //System.out.println("LineProf: doppler, logDopp: " + doppler + " " + logE*logDopp);

        //Put input parameters into linear cgs units:
        //double gammaCol = Math.pow(10.0, logGammaCol);
        // Lorentzian broadening:
        // Assumes Van der Waals dominates radiative damping
        // log_10 Gamma_6 for van der Waals damping around Tau_Cont = 1 in Sun 
        //  - p. 57 of Radiative Transfer in Stellar Atmospheres (Rutten)
        double logGammaSun = 9.0 * ln10; // Convert to base e 
        //double logFudge = Math.log(2.5);  // Van der Waals enhancement factor

        //Get index of depth nearest to tauRos=1:
        int tau1 = TauPoint.tauPoint(numDeps, tauRos, 1.0);

        //System.out.println("LINEGRID: Tau1: " + tau1);
        //logA = 2.0 * logLam0 + logGamma - ln4pi - logC - logDopp;
        //a = Math.exp(logA);
        //System.out.println("LINEGRID: logA: " + logE * logA);
        //Set up a half-profile Delta_lambda grid in Doppler width units 
        //from line centre to wing
        int numPoints = linePoints[0].length;
        //double maxV = linePoints[1][numPoints - 1]; //maximum value of v in Doppler widths
        //System.out.println("LineProf: numPoints: " + numPoints);

        // Return a 2D numPoints X numDeps array of normalized line profile points (phi)
        double[][] lineProf = new double[numPoints][numDeps];

        // Line profiel points in Doppler widths - needed for Voigt function, H(a,v):
        double[] v = new double[numPoints];
        double logV, ii;

//        lineProf[0][0] = 0.0; v[0] = 0.0; //Line centre - cannot do logaritmically!
        double gamma, logGamma, a, logA, voigt, core, wing, logWing, logVoigt;

        //Variables for managing honest Voigt profile convolution:
        //Find value of damping a parameter around tauRos=1.0 for guidance:
        logGamma = press[1][tau1] - pressSun[1][tau1] + 0.7 * (tempSun[1][tau1] - temp[1][tau1]) + logGammaSun;
        logGamma = logGamma + logGammaCol;
        //Voigt "a" parameter in Doppler widths with line centre wavelength:
        logA = 2.0 * logLam0 + logGamma - ln4pi - logC - logDopp;
        a = Math.exp(logA);
        double yLim0, yLim2, nApprox, deltaY, yLimLower;
        double deltaYApprox = 0.2;
        int numY0, numY1, numY2, numY;
        //yLim0 = Math.max(3.5 * a, 3.5); //"a" is characteristic width of Lorentizian in Doppler widths
        yLim0 = Math.max(a, 1.0); //"a" is characteristic width of Lorentizian in Doppler widths
        //CAUTION: integration variable, y, is in Doppler widths
        numY0 = (int) Math.round(yLim0 / deltaYApprox);
        //numY0 = Math.round(yLim0 / deltaYApprox);
        //System.out.println("a " + a + " yLim0 " + yLim0 + " numY0 " + numY0);
        yLim2 = Math.max(3.0 * a, 3.0); //"a" is characteristic width of Lorentizian
        numY2 = (int) Math.round(yLim2 / deltaYApprox);
        //numY2 = Math.round(yLim2 / deltaYApprox);
        double y, logNumerator, logDenominator, logDenomTerm1, logDenomTerm2, denomTerm1, denomTerm2;
        double integ, logInteg, logNextInteg, nextInteg, term;

        for (int id = 0; id < numDeps; id++) {

            //Formula from p. 56 of Radiative Transfer in Stellar Atmospheres (Rutten),
            // logarithmically with respect to solar value:
            logGamma = press[1][id] - pressSun[1][tau1] + 0.7 * (tempSun[1][tau1] - temp[1][id]) + logGammaSun;
            //logGamma = logGamma + logFudge + logGammaCol;
            logGamma = logGamma + logGammaCol;
            //System.out.println("LineGrid: logGamma: " + id + " " + logE * logGamma);

            //Voigt "a" parameter in Doppler widths with line centre wavelength:
            logA = 2.0 * logLam0 + logGamma - ln4pi - logC - logDopp;
            a = Math.exp(logA);

            //System.out.println("LineGrid: logGamma: " + id + " " + logE * logGamma + " " + logE * logA);
            //limits and sampling of integration variable y=(xi/c)*(lambda0/Doppler) (Rutten p. 59) - is this Doppler shift in Doppler widths?
            //Notes: the e^-y^2 numerator is even about y=0, BUT the ((v-y)^2 + a^2 denominator is even about y=v!
            //In the sum over y, the leading term according the numerator alone is always at y=0, 
            // BUT the leading term according to the denominator alone is always at y=v
            // --> Therefore the y integration range should always include both y=0 and y=v (for both +ve and -ve v)
            // --> Adapt y range to each v??
            // Can the y sampling be uniform???
            //int numY = 2 * (int) (yLim / deltaY);
            //System.out.println("numY " + numY);
//            int numY = 20; //I dunno - this is how many i want :-)
//            double deltaY = 2.0 * yLim / ((double) numY);
            //double yMin = -1.0 * yLim;
            //if (id == 30) {
            //    System.out.println("il   v[il]  iy   y   logNumerator   logDenominator   logInteg ");
            //    //System.out.println("voigt2:   v   logVoigt: ");
            //}
            //for (int il = 0; il < numPoints; il++) {
            //Negative half of profiel only - numPoints is always odd and symmetrical about lambda=lambda0
            for (int il = 0; il < (numPoints / 2) + 1; il++) {

                v[il] = linePoints[1][il];
                //System.out.println("LineProf: il, v[il]: " + il + " " + v[il]);
                if (Math.abs(v[il]) < deltaYApprox) {
                    numY1 = 0;
                    deltaY = deltaYApprox;
                } else {
                    nApprox = -1.0 * v[il] / deltaYApprox; // -1.0* becasue we are only treating v<=0 half of profile
                    numY1 = (int) Math.round(nApprox);
                    //numY1 = Math.round(nApprox);
                    deltaY = -1.0 * v[il] / numY1;
                }
                yLimLower = v[il] - deltaY * numY0;
                //yLimUpper = deltaY * numY2; //Do we need this??

                //Total integration range is from yLimLower to v[il] to 0 to yLim:
                numY = numY0 + numY1 + numY2;
                double yLimUpper = yLimLower + numY * deltaY;
                // if (id == 2) {
                //     if (v[il] == -4.0) {
                //         System.out.println("il " + il + " v " + v[il] + " numY0 " + numY0 + " numY1 " + numY1 + " numY2 " + numY2 + " numY " + numY);
                //         System.out.println("il " + il + " v " + v[il] + " yLimLower " + yLimLower + " deltaY " + deltaY + " yLimUpper "
                //                 + yLimUpper);
                //    }
                // }
                voigt = 0.0;// re-initialize H(a,v):
                // loop over 
                //if (id == 35) {
                //    System.out.println("il " + il + " v " + v[il]);
                //}

                // Trapezoid integration: Compute integrand at first y point:
                y = yLimLower;
                logNumerator = -1.0 * y * y;
                logDenomTerm1 = 2.0 * Math.log(Math.abs(v[il] - y));
                logDenomTerm2 = 2.0 * logA;
                denomTerm1 = Math.exp(logDenomTerm1);
                denomTerm2 = Math.exp(logDenomTerm2);
                logDenominator = Math.log(denomTerm1 + denomTerm2);

                logInteg = logNumerator - logDenominator;
                integ = Math.exp(logInteg);

                for (int iy = 1; iy < numY; iy++) {

                    y = ((double) iy) * deltaY + yLimLower;

                    logNumerator = -1.0 * y * y;
                    logDenomTerm1 = 2.0 * Math.log(Math.abs(v[il] - y));
                    logDenomTerm2 = 2.0 * logA;
                    denomTerm1 = Math.exp(logDenomTerm1);
                    denomTerm2 = Math.exp(logDenomTerm2);
                    logDenominator = Math.log(denomTerm1 + denomTerm2);

                    logNextInteg = logNumerator - logDenominator;
                    nextInteg = Math.exp(logNextInteg);
                    //Trapezoid integration:
                    term = 0.5 * (integ + nextInteg) * deltaY;

                    // Rectangular pick integration for now:
                    //Skip Lorentzian peak - it may be giving us too large a line opacity:
                    if (Math.abs(v[il] - y) > 2.0 * deltaYApprox) {
                        voigt = voigt + term;
                        //voigt = voigt + integ * deltaY;
                    }

                    // if (id == 2) {
                    //     if (v[il] == -4.0) {
                    //         System.out.println("il, v[il], iy, y, logNumerator, logDenominator, logInteg");
                    //         System.out.format("%02d   %12.8f  %02d   %12.8f   %12.8f   %12.8f   %12.8f%n", il, v[il], iy, y, logNumerator, logDenominator, logInteg);
                    //         System.out.println("iy, y, logE*logInteg");
                    //        System.out.format("%03d   %12.8f   %12.8f%n", iy, y, logE * logInteg);
                    //    }
                    // }
                    //Only treating negative v half of profile, so Lorentzian denominator is always minimal at negative or zero y
                    // --> If y > 0 be prepared to stop when contribution becomes negligible:
                    //if (y > 0 && (term / voigt < 0.01)) {
                    //     break;  //break out of iy loop
                    //}
                    integ = nextInteg;
                }  //iy loop

                //Pre-factor for H(a,v) integral
                logVoigt = Math.log(voigt) + logA - Math.log(Math.PI);

                //System.out.println("LINEGRID: il, v[il]: " + il + " " + v[il] + " lineProf[0][il]: " + lineProf[0][il]);
                //System.out.println("LINEGRID: il, Voigt, H(): " + il + " " + voigt);
                //Convert from H(a,v) in dimensionless Voigt units to physical phi((Delta lambda) profile:
                logVoigt = logVoigt + 2.0 * logLam0 - lnSqRtPi - logDopp - logC;
                //if (id == 35) {
                //    System.out.format("%12.8f   %12.8f%n", v[il], logE * logVoigt);
                //}
                lineProf[il][id] = Math.exp(logVoigt);

                //Copy over to positive half of profile - avoid doubling the central wavelength!
                if (il != ((numPoints - 1) - il)) {
                    lineProf[(numPoints - 1) - il][id] = Math.exp(logVoigt);
                }
                //System.out.println("LineProf: il, id, lineProf[il][id]: " + il + " " + id + " " + lineProf[il][id]);
            } // il lambda loop

            //if (id == 20) {
            //    for (int il = 0; il < numPoints; il++) {
            //        System.out.format("Voigt2: %20.16f   %20.16f%n", linePoints[1][il], logE * Math.log(lineProf[il][id]));
            //    }
            //}
        } //id loop

        return lineProf;

    }

    // Make line source function:
    // Equivalenth two-level atom (ETLA) approx
    //CAUTION: input lambda in nm
    public static double[] lineSource(int numDeps, double[][] tau, double[][] temp, double lambda) {

        double[] lineSource = new double[numDeps];

        //thermal photon creation/destruction probability
        double epsilon = 0.01; //should decrease with depth??

        //This is an artifact of jayBinner's original purpose:
        double grayLevel = 1.0;

        //int iLam0 = numLams / 2; //+/- 1 deltaLambda
        //double lam0 = linePoints[0][iLam0];  //line centre lambda in cm - not needed:
        //double lamStart = lambda - 0.1; // nm
        //double lamStop = lambda + 0.1; // nm
        //double lamRange = (lamStop - lamStart); // * 1.0e-7; // line width in cm
        //System.out.println("lamStart " + lamStart + " lamStop " + lamStop + " lamRange " + lamRange);
        double[] jayLambda = new double[numDeps];
        double[][] BLambda = new double[2][numDeps];
        double linSrc;

        // Dress up Blambda to look like what jayBinner expects:
        for (int i = 0; i < numDeps; i++) {
            //Planck.planck return log(B_lambda):
            BLambda[0][i] = Math.exp(Planck.planck(temp[0][i], lambda));
            BLambda[1][i] = 1.0;  //supposed to be dB/dT, but not needed. 
        }

        //CAUTION: planckBin Row 0 is linear lambda-integrated B_lambda; Row 1 is same for dB_lambda/dT
        //planckBin = MulGrayTCorr.planckBinner(numDeps, temp, lamStart, lamStop);
        jayLambda = MulGrayTCorr.jayBinner(numDeps, tau, temp, BLambda, grayLevel);
        //To begin with, coherent scattering - we're not computing line profile-weighted average Js and Bs
        for (int i = 0; i < numDeps; i++) {

            //planckBin[0][i] = planckBin[0][i] / lamRange;  //line average
            //jayBin[i] = jayBin[i];  
            linSrc = (1.0 - epsilon) * jayLambda[i] + epsilon * BLambda[0][i];
            lineSource[i] = Math.log(linSrc);
        }

        return lineSource;
    }

}
