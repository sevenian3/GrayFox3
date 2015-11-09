/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package grayfox;

/**
 *
 * @author Ian
 */
public class ScaleSolar {

    public static double phxSunTeff = 5777.0;
    public static double phxSunLogEg = Math.log(10.0) * 4.44;  //base e!

    //Corresponding Tau_1200 grid (ie. lambda_0 = 1200 nm):    
    public static double[] phxSunTau64 = {
        0.00000000000000000e+00, 9.99999999999999955e-07, 1.34596032415536424e-06,
        1.81160919420041334e-06, 2.43835409826882661e-06, 3.28192787251147086e-06,
        4.41734470314007309e-06, 5.94557070854439435e-06, 8.00250227816105150e-06,
        1.07710505603676912e-05, 1.44974067037263168e-05, 1.95129342263596224e-05,
        2.62636352765333536e-05, 3.53498110503010944e-05, 4.75794431400941376e-05,
        6.40400427119728256e-05, 8.61953566475303296e-05, 1.16015530173997152e-04,
        1.56152300600049664e-04, 2.10174801133248704e-04, 2.82886943462596928e-04,
        3.80754602122237184e-04, 5.12480587696093120e-04, 6.89778537938765824e-04,
        9.28414544519474432e-04, 1.24960914129198688e-03, 1.68192432488086880e-03,
        2.26380340952144672e-03, 3.04698957090350784e-03, 4.10112707055130048e-03,
        5.51995432128156800e-03, 7.42963950759494912e-03, 1.00000000000000002e-02,
        1.34596032415536416e-02, 1.81160919420041312e-02, 2.43835409826882656e-02,
        3.28192787251147072e-02, 4.41734470314006464e-02, 5.94557070854439424e-02,
        8.00250227816105216e-02, 1.07710505603676912e-01, 1.44974067037263136e-01,
        1.95129342263596224e-01, 2.62636352765332992e-01, 3.53498110503010240e-01,
        4.75794431400941440e-01, 6.40400427119728384e-01, 8.61953566475303168e-01,
        1.16015530173997152e+00, 1.56152300600049664e+00, 2.10174801133248704e+00,
        2.82886943462596640e+00, 3.80754602122236800e+00, 5.12480587696092608e+00,
        6.89778537938765824e+00, 9.28414544519474432e+00, 1.24960914129198672e+01,
        1.68192432488086880e+01, 2.26380340952144640e+01, 3.04698957090350528e+01,
        4.10112707055129536e+01, 5.51995432128157312e+01, 7.42963950759495040e+01,
        1.00000000000000000e+02
    };

    public static double[] logPhxSunTau64() {

        double logE = Math.log10(Math.E);

        int numPhxDep = ScaleSolar.phxSunTau64.length;
        double logPhxSunTau64[] = new double[numPhxDep];
        for (int i = 1; i < numPhxDep; i++) {
            logPhxSunTau64[i] = Math.log(ScaleSolar.phxSunTau64[i]);
        }
        logPhxSunTau64[0] = logPhxSunTau64[1] - (logPhxSunTau64[numPhxDep - 1] - logPhxSunTau64[1]) / numPhxDep;
        return logPhxSunTau64;
    }

    public static double[][] phxSunTemp(double teff, int numDeps, double[][] tauRos) {

        double logE = Math.log10(Math.E);

        //Theoretical radiative/convective model from Phoenix V15:
        double[] phxSunTemp64 = {
            3.75778887392339840e+03, 3.75778887392339840e+03, 3.78480175327941504e+03,
            3.81385432525541760e+03, 3.84360130602512768e+03, 3.87340585446516608e+03,
            3.90300184305606656e+03, 3.93231689265254528e+03, 3.96137919852984000e+03,
            3.99027119028325824e+03, 4.01910484194699648e+03, 4.04798292490651008e+03,
            4.07699548886169152e+03, 4.10623218035810816e+03, 4.13574364539801920e+03,
            4.16548101060783104e+03, 4.19541371831173824e+03, 4.22551121760088000e+03,
            4.25571229065970624e+03, 4.28594188575783232e+03, 4.31613168919769152e+03,
            4.34620698440244928e+03, 4.37603327507328960e+03, 4.40564394765877952e+03,
            4.43507740841559296e+03, 4.46439148496796224e+03, 4.49375530130093952e+03,
            4.52341166116436480e+03, 4.55357281866347264e+03, 4.58446079852491520e+03,
            4.61663974201107520e+03, 4.65052341797810624e+03, 4.68623381803595456e+03,
            4.72408924142126144e+03, 4.76494152329308416e+03, 4.80984310271200128e+03,
            4.85897778977827584e+03, 4.91315894280032960e+03, 4.97390461818851328e+03,
            5.04531167969494336e+03, 5.12680296183560704e+03, 5.22061204180252480e+03,
            5.32918534350649152e+03, 5.46202432323604352e+03, 5.61966782651567040e+03,
            5.80986721241013376e+03, 6.03911828822760320e+03, 6.23433005487621120e+03,
            6.53458311644527488e+03, 6.87429103746811904e+03, 7.29999981509928192e+03,
            7.66682942009826304e+03, 7.94223816217841024e+03, 8.16133659245977728e+03,
            8.35020013757955200e+03, 8.52047273964030720e+03, 8.67812135633704064e+03,
            8.82687568743616768e+03, 8.96926538519515648e+03, 9.10706359999037824e+03,
            9.24154121553023488e+03, 9.37363000902155008e+03, 9.50427569030960000e+03,
            9.63219702937432192e+03
        };

        // interpolate onto gS3 tauRos grid and re-scale with Teff:
        double[] phxSunTemp = new double[numDeps];
        double[][] scaleTemp = new double[2][numDeps];
        for (int i = 0; i < numDeps; i++) {
            phxSunTemp[i] = Interpol.interpol(ScaleSolar.logPhxSunTau64(), phxSunTemp64, tauRos[1][i]);
            scaleTemp[0][i] = teff * phxSunTemp[i] / ScaleSolar.phxSunTeff;
            scaleTemp[1][i] = Math.log(scaleTemp[0][i]);
            //System.out.println("tauRos[1][i] " + logE * tauRos[1][i] + " scaleTemp[1][i] " + logE * scaleTemp[1][i]);
        }

        return scaleTemp;

    }

    public static double[][] phxSunPGas(double grav, int numDeps, double[][] tauRos) {

        double logE = Math.log10(Math.E);
        double logEg = Math.log(grav); //base e!

        //Theoretical radiative/convective model from Phoenix V15:
        double[] phxSunPGas64 = {
            1.00000000000000005e-04, 7.28828683006412544e+01, 8.61732126528505984e+01,
            1.01843641855932976e+02, 1.20317369304629504e+02, 1.42093296011949696e+02,
            1.67758727999644384e+02, 1.98004769223716256e+02, 2.33644726494082176e+02,
            2.75635953319757664e+02, 3.25104809120938880e+02, 3.83378880706399168e+02,
            4.52022443862726592e+02, 5.32877321364649344e+02, 6.28113128741022208e+02,
            7.40284569989930496e+02, 8.72399144145001216e+02, 1.02799724165148560e+03,
            1.21124571517496000e+03, 1.42704756928025427e+03, 1.68117132827309248e+03,
            1.98040330055171296e+03, 2.33272402439094176e+03, 2.74752260927171264e+03,
            3.23584954067544384e+03, 3.81071167175796544e+03, 4.48742481283848128e+03,
            5.28403135449994368e+03, 6.22178478543013120e+03, 7.32571484052561408e+03,
            8.62531818498740864e+03, 1.01553497268350960e+04, 1.19567104697253520e+04,
            1.40775115384306991e+04, 1.65743702896828832e+04, 1.95139034178162464e+04,
            2.29742653550211872e+04, 2.70468752817440448e+04, 3.18381441391788864e+04,
            3.74704748898233472e+04, 4.40799661582952512e+04, 5.18080650391892096e+04,
            6.07793647633492224e+04, 7.10351288049853440e+04, 8.24259773567987968e+04,
            9.44866985169806080e+04, 1.06329924298695632e+05, 1.17862219382348656e+05,
            1.28295128203359424e+05, 1.36933948396180352e+05, 1.43493910023715958e+05,
            1.48487688700034048e+05, 1.52795575243316608e+05, 1.56932489940248512e+05,
            1.61140965195830048e+05, 1.65564070780028256e+05, 1.70312554701480352e+05,
            1.75486986284790656e+05, 1.81187218697219744e+05, 1.87518050413513344e+05,
            1.94593473563783808e+05, 2.02540389901047584e+05, 2.11500759107428064e+05,
            2.21643078023966592e+05

        };

        int numPhxDeps = phxSunPGas64.length;  //yeah, I know, 64, but that could change!
        double[] logPhxSunPGas64 = new double[numPhxDeps];
        for (int i = 0; i < phxSunPGas64.length; i++) {
            logPhxSunPGas64[i] = Math.log(phxSunPGas64[i]);
        }

        // interpolate onto gS3 tauRos grid and re-scale with Teff:
        double[] phxSunPGas = new double[numDeps];
        double[] logPhxSunPGas = new double[numDeps];
        double[][] scalePGas = new double[2][numDeps];
        for (int i = 0; i < numDeps; i++) {
            logPhxSunPGas[i] = Interpol.interpol(ScaleSolar.logPhxSunTau64(), logPhxSunPGas64, tauRos[1][i]);
            scalePGas[1][i] = logEg + logPhxSunPGas[i] - ScaleSolar.phxSunLogEg;
            scalePGas[0][i] = Math.exp(scalePGas[1][i]);
            //System.out.println("scalePGas[1][i] " + logE * scalePGas[1][i]);
        }

        return scalePGas;

    }

    public static double[][] phxSunNe(double grav, int numDeps, double[][] tauRos, double[][] scaleTemp, double kappaScale) {

        double logE = Math.log10(Math.E);
        double logEg = Math.log(grav); //base e!
        double logEkappaScale = Math.log(kappaScale);

        //Theoretical radiative/convective model from Phoenix V15:
        double[] phxSunPe64 = {
            1.53086468021591745e-07, 5.66518458165471424e-03, 6.72808433760886656e-03,
            8.00271552708326656e-03, 9.51809762875982208e-03, 1.13117438884935648e-02,
            1.34299756939525680e-02, 1.59287848014678144e-02, 1.88751877391284448e-02,
            2.23491173128862976e-02, 2.64457686695698400e-02, 3.12779350532322240e-02,
            3.69791374171045888e-02, 4.37078139287801024e-02, 5.16503829681397248e-02,
            6.10221573903118336e-02, 7.20768505868849536e-02, 8.51123959415642752e-02,
            1.00475763241309840e-01, 1.18571138726675232e-01, 1.39870552376136714e-01,
            1.64923053015554560e-01, 1.94357063774820192e-01, 2.28928720249475840e-01,
            2.69525262128246720e-01, 3.17192228891198592e-01, 3.73192988074577856e-01,
            4.39058414038311360e-01, 5.16615873984964544e-01, 6.08066526878471680e-01,
            7.16264581324812416e-01, 8.44657163125294336e-01, 9.97267452897639808e-01,
            1.17915717019238848e+00, 1.39715732004723136e+00, 1.66026825646718432e+00,
            1.97886823850223904e+00, 2.36716912384854112e+00, 2.84540915928013805e+00,
            3.44853013665125120e+00, 4.21529199485384704e+00, 5.21488490421314560e+00,
            6.56660005867586432e+00, 8.55643059606379776e+00, 1.16931723772200080e+01,
            1.71629079266534368e+01, 2.75152019254691616e+01, 4.18720694941323264e+01,
            7.66283674228108288e+01, 1.45995186997127872e+02, 3.04766672331673792e+02,
            5.44151864837275328e+02, 8.17181982032739072e+02, 1.11216222784450608e+03,
            1.43633935534913856e+03, 1.79603721463325728e+03, 2.19692608617747040e+03,
            2.64548745663525184e+03, 3.14931730610757952e+03, 3.71721361233669376e+03,
            4.35932065708395904e+03, 5.08736399892079296e+03, 5.91634943413070720e+03,
            6.85104524590000384e+03
        };

        int numPhxDeps = phxSunPe64.length;  //yeah, I know, 64, but that could change!
        double[] logPhxSunPe64 = new double[numPhxDeps];
        for (int i = 0; i < phxSunPe64.length; i++) {
            logPhxSunPe64[i] = Math.log(phxSunPe64[i]);
        }

        // interpolate onto gS3 tauRos grid and re-scale with Teff:
        double[] phxSunPe = new double[numDeps];
        double[] logPhxSunPe = new double[numDeps];
        double[] logScalePe = new double[numDeps];
        double[][] scaleNe = new double[2][numDeps];

        for (int i = 0; i < numDeps; i++) {
            logPhxSunPe[i] = Interpol.interpol(ScaleSolar.logPhxSunTau64(), logPhxSunPe64, tauRos[1][i]);
            logScalePe[i] = logEg + logPhxSunPe[i] - ScaleSolar.phxSunLogEg - logEkappaScale;
            scaleNe[1][i] = logScalePe[i] - scaleTemp[1][i] - Useful.logK();
            scaleNe[0][i] = Math.exp(scaleNe[1][i]);
            //System.out.println("scaleNe[1][i] " + logE * scaleNe[1][i]);
        }

        return scaleNe;

    }

    //Try to recover the opacity as lambda_0 = 1200 nm:
    public static double[][] phxSunKappa(int numDeps, double[][] tauRos, double kappaScale) {

        double logEkappaScale = Math.log(kappaScale);

        //Theoretical radiative/convective model from Phoenix V15:
        double[] phxSunRho64 = {
            4.13782346832222649e-16, 3.02095569469690462e-10, 3.54633225055968270e-10,
            4.15928280610231993e-10, 4.87569895799879155e-10, 5.71381142733345291e-10,
            6.69468927495419999e-10, 7.84278468388299388e-10, 9.18654436245877140e-10,
            1.07590983297567878e-09, 1.25990158939278389e-09, 1.47513757382262481e-09,
            1.72688539188771193e-09, 2.02128936476074103e-09, 2.36554000030610158e-09,
            2.76809615861929229e-09, 3.23884396019102352e-09, 3.78934920783997866e-09,
            4.43317360103421215e-09, 5.18621173362546736e-09, 6.06707380164391496e-09,
            7.09757215466433105e-09, 8.30337600953291647e-09, 9.71426731449415417e-09,
            1.13650770268615465e-08, 1.32964932176367733e-08, 1.55557163673284530e-08,
            1.81974840999693492e-08, 2.12855768344032029e-08, 2.48940684847852482e-08,
            2.91068454381155637e-08, 3.40213170202104799e-08, 3.97519122004400661e-08,
            4.64290866159173997e-08, 5.41967343519845744e-08, 6.32144869975830899e-08,
            7.36729431582295057e-08, 8.57774421976652924e-08, 9.97399445761737017e-08,
            1.15721981027072251e-07, 1.33967659681056212e-07, 1.54620178670780798e-07,
            1.77690495649821781e-07, 2.02608223525831620e-07, 2.28481547026651195e-07,
            2.53309018291389784e-07, 2.74195019891415717e-07, 2.94373976046088894e-07,
            3.05614181338722779e-07, 3.09912387277346887e-07, 3.05484245799381785e-07,
            3.00519445088246902e-07, 2.98007120264342719e-07, 2.97336159154754909e-07,
            2.97854109132361140e-07, 2.99327766949861546e-07, 3.01691329467384893e-07,
            3.04944348605014908e-07, 3.09125225055924192e-07, 3.14302162196028050e-07,
            3.20569231575000568e-07, 3.28044919674719785e-07, 3.36858977566225440e-07,
            3.47271781807407172e-07
        };

        double[] phxSunRadius64 = {
            9.98760000000000000e+10, 9.98660572490945152e+10, 9.98645871807186304e+10,
            9.98631098643980160e+10, 9.98616245003269760e+10, 9.98601306458076032e+10,
            9.98586280682994048e+10, 9.98571166428681216e+10, 9.98555962828737792e+10,
            9.98540668955362944e+10, 9.98525283799022080e+10, 9.98509805586940416e+10,
            9.98494232096872704e+10, 9.98478561022866944e+10, 9.98462789875034752e+10,
            9.98446916403608064e+10, 9.98430938763377024e+10, 9.98414855440511616e+10,
            9.98398665329129600e+10, 9.98382367818977152e+10, 9.98365962762478464e+10,
            9.98349450434777856e+10, 9.98332831693342848e+10, 9.98316107506358144e+10,
            9.98299278514395904e+10, 9.98282344996977408e+10, 9.98265306530218624e+10,
            9.98248161610832000e+10, 9.98230907740896512e+10, 9.98213541602841472e+10,
            9.98196058171550848e+10, 9.98178450629064448e+10, 9.98160711682936960e+10,
            9.98142833645708160e+10, 9.98124806655167488e+10, 9.98106617425436544e+10,
            9.98088251712680448e+10, 9.98069695137641728e+10, 9.98050931816160256e+10,
            9.98031941655960192e+10, 9.98012715884401664e+10, 9.97993275763000704e+10,
            9.97973698841808128e+10, 9.97954186071983616e+10, 9.97935139683624704e+10,
            9.97917197632915456e+10, 9.97901090800493440e+10, 9.97886636105590528e+10,
            9.97874361487011200e+10, 9.97864521731274880e+10, 9.97857037165240960e+10,
            9.97851119909964288e+10, 9.97845890321633152e+10, 9.97840832513957504e+10,
            9.97835682848038784e+10, 9.97830286360697344e+10, 9.97824528501662336e+10,
            9.97818311493811456e+10, 9.97811545315540224e+10, 9.97804143368400768e+10,
            9.97796020159347072e+10, 9.97787090120484864e+10, 9.97777268638511104e+10,
            9.97766460582020224e+10
        };

        int numPhxDeps = phxSunRadius64.length;
        double[] phxSunKappa64 = new double[numPhxDeps];
        double[] logPhxSunKappa64 = new double[numPhxDeps];
        //double[] logPhxSunRho64 = new double[numPhxDeps];
        //double[] logPhxSunRadius64 = new double[numPhxDeps];

//Fix to get right depth scale and right line strengths:
// Yeah - everywhere ya go - opacity fudge
        double fudge = 0.25;
        double logFudge = Math.log(fudge);
        double deltaRho, deltaRadius, deltaTau, logDeltaRho, logDeltaRadius, logDeltaTau;
        double logE = Math.log10(Math.E);
        for (int i = 1; i < numPhxDeps; i++) {

            //Renormalize radii before taking difference
            //Caution: Radius *decreases* with increasing i (inward) and we'll be taking the log:
            deltaRadius = (1.0e-11 * phxSunRadius64[i - 1]) - (1.0e-11 * phxSunRadius64[i]);
            deltaRadius = Math.abs(deltaRadius);
            //restore to cm:
            deltaRadius = 1.0e11 * deltaRadius;
            //Renormalize before taking rho difference
            deltaRho = (1.0e9 * phxSunRho64[i]) - (1.0e9 * phxSunRho64[i - 1]);
            deltaRho = Math.abs(deltaRho);
            //Restore g/cm^3:
            deltaRho = 1.0e-9 * deltaRho;
            //Renormalize before taking rho difference
            deltaTau = (1.0e2 * ScaleSolar.phxSunTau64[i]) - (1.0e2 * ScaleSolar.phxSunTau64[i - 1]);
            deltaTau = Math.abs(deltaTau);
            deltaTau = 1.0e-2 * deltaTau;

            logDeltaRadius = Math.log(deltaRadius);
            logDeltaRho = Math.log(deltaRho);
            logDeltaTau = Math.log(deltaTau);

            logPhxSunKappa64[i] = logDeltaTau - logDeltaRho - logDeltaRadius - logEkappaScale + logFudge;
            phxSunKappa64[i] = Math.exp(logPhxSunKappa64[i]);
            //System.out.println("logPhxSunKappa64[i] " + logE*logPhxSunKappa64[i]);

        }

        logPhxSunKappa64[0] = logPhxSunKappa64[1];
        phxSunKappa64[0] = phxSunKappa64[1];

        // interpolate onto gS3 tauRos grid and re-scale with Teff:
        double[][] phxSunKappa = new double[2][numDeps];
        for (int i = 0; i < numDeps; i++) {
            phxSunKappa[1][i] = Interpol.interpol(ScaleSolar.logPhxSunTau64(), logPhxSunKappa64, tauRos[1][i]);
            phxSunKappa[0][i] = Math.exp(phxSunKappa[1][i]);
            //System.out.println("phxSunKappa[1][i], i= " + i + " " + logE * phxSunKappa[1][i]);
        }

        return phxSunKappa;

    }

}
