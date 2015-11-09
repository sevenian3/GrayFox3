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
public class ScaleVega {

    //Castelli & Kurucz (ALMOST!)
    public static double phxVegaTeff = 9950.0;  //actual INCORRECT Teff used in Phoenix model (should have been 9550 -sigh!)
    public static double phxVegaLogEg = Math.log(10.0) * 3.95; //base e
    public static double phxVegaLogEkappaScale = Math.log(10.0) * (-0.5);  //base!

    //Corresponding Tau_500 grid (ie. lambda_0 = 500 nm):    
    public static double[] phxVegaTau64 = {
        0.00000000000000000e+00, 1.00000000000000004e-10, 1.57297315730079583e-10,
        2.47424455358884461e-10, 3.89192026739294322e-10, 6.12188611096406130e-10,
        9.62956252459902995e-10, 1.51470433677439602e-09, 2.38258926299323964e-09,
        3.74774895556145216e-09, 5.89510850740025871e-09, 9.27284744151620558e-09,
        1.45859401172503535e-08, 2.29432922784316770e-08, 3.60891828940797229e-08,
        5.67673159613064525e-08, 8.92934642191482617e-08, 1.40456222339119575e-07,
        2.20933867515307984e-07, 3.47523043140230439e-07, 5.46644418403068955e-07,
        8.59856996736334509e-07, 1.35253197498353518e-06, 2.12749649104013246e-06,
        3.34649487265776861e-06, 5.26394660573542543e-06, 8.28004671228647794e-06,
        1.30242912196233633e-05, 2.04868604813359955e-05, 3.22252816145080504e-05,
        5.06895029660801231e-05, 7.97332275225631174e-05, 1.25418226637949074e-04,
        1.97279503937761955e-04, 3.10315364179716846e-04, 4.88117738152716565e-04,
        7.67796099716601789e-04, 1.20772265513446235e-03, 1.89971531799055940e-03,
        2.98820120171229553e-03, 4.70036027890743165e-03, 7.39354054836428853e-03,
        1.16298408199920333e-02, 1.82934274335286202e-02, 2.87750703079705100e-02,
        4.52624131938807600e-02, 7.11965609886321127e-02, 1.11990279327247338e-01,
        1.76157703260379023e-01, 2.77091338680335086e-01, 4.35857237864710867e-01,
        6.85591735576461137e-01, 1.00000000000000000e+00, 1.07841739692903849e+00,
        1.69632161773558221e+00, 2.66826837084713286e+00, 4.19711452381726513e+00,
        6.60194848408189738e+00, 1.03846877513435061e+01, 1.63348350798136970e+01,
        2.56942571094824572e+01, 4.04163767300010406e+01, 6.35738757116481565e+01,
        1.00000000000000000e+02
    };

    public static double[] logPhxVegaTau64() {

        double logE = Math.log10(Math.E);

        int numPhxDep = ScaleVega.phxVegaTau64.length;
        double logPhxVegaTau64[] = new double[numPhxDep];
        for (int i = 1; i < numPhxDep; i++) {
            logPhxVegaTau64[i] = Math.log(ScaleVega.phxVegaTau64[i]);
        }
        logPhxVegaTau64[0] = logPhxVegaTau64[1] - (logPhxVegaTau64[numPhxDep - 1] - logPhxVegaTau64[1]) / numPhxDep;
        return logPhxVegaTau64;
    }

    public static double[][] phxVegaTemp(double teff, int numDeps, double[][] tauRos) {

        double logE = Math.log10(Math.E);

        //Theoretical radiative/convective model from Phoenix V15:
        double[] phxVegaTemp64 = {
            5.92666192154309101e+03, 5.92666192154309101e+03, 5.92669171597276090e+03,
            5.92673998930622656e+03, 5.92681897305491384e+03, 5.92694950529302332e+03,
            5.92716641457396599e+03, 5.92752552761684638e+03, 5.92810993487369979e+03,
            5.92903227707427686e+03, 5.93043220345250484e+03, 5.93247494804068447e+03,
            5.93535705326493917e+03, 5.93931461731802392e+03, 5.94462395405969710e+03,
            5.95160056847340456e+03, 5.96062955115910609e+03, 5.97234442427941121e+03,
            5.98799376822690192e+03, 6.00995395792729505e+03, 6.04217083863954394e+03,
            6.08927266539518769e+03, 6.15482922005066484e+03, 6.24095845993735929e+03,
            6.34702439197141939e+03, 6.46752230210132257e+03, 6.59445408379155560e+03,
            6.72080130115506290e+03, 6.84109768032372904e+03, 6.95113817061313694e+03,
            7.04860283942552360e+03, 7.13460514913545467e+03, 7.21354674932631588e+03,
            7.28978629931232808e+03, 7.36559554177949485e+03, 7.44189226318755846e+03,
            7.51906108650742408e+03, 7.59685046298776069e+03, 7.67490002230248228e+03,
            7.75361847195948849e+03, 7.83483766810037559e+03, 7.92234631761687160e+03,
            8.02067697937434696e+03, 8.13521360575207927e+03, 8.27574328647004768e+03,
            8.45212850359107870e+03, 8.67381075095960477e+03, 8.95207036703875929e+03,
            9.30689214747616825e+03, 9.74795486690789039e+03, 1.02809293652303932e+04,
            1.09262232212095878e+04, 1.12593239932788256e+04, 1.15410572182275682e+04,
            1.20102451746046791e+04, 1.31012465922441861e+04, 1.39755996506773790e+04,
            1.49601225189046345e+04, 1.60017842217904235e+04, 1.71032050644639567e+04,
            1.83299248272271325e+04, 1.96750591610957927e+04, 2.11303208351190988e+04,
            2.28754881624451700e+04
        };

        // interpolate onto gS3 tauRos grid and re-scale with Teff:
        double[] phxVegaTemp = new double[numDeps];
        double[][] scaleTemp = new double[2][numDeps];
        for (int i = 0; i < numDeps; i++) {
            phxVegaTemp[i] = Interpol.interpol(ScaleVega.logPhxVegaTau64(), phxVegaTemp64, tauRos[1][i]);
            scaleTemp[0][i] = teff * phxVegaTemp[i] / ScaleVega.phxVegaTeff;
            scaleTemp[1][i] = Math.log(scaleTemp[0][i]);
            //System.out.println("tauRos[1][i] " + logE * tauRos[1][i] + " scaleTemp[1][i] " + logE * scaleTemp[1][i]);
        }

        return scaleTemp;

    }

    public static double[][] phxVegaPGas(double grav, int numDeps, double[][] tauRos) {

        double logE = Math.log10(Math.E);
        double logEg = Math.log(grav); //base e!

        //Theoretical radiative/convective model from Phoenix V15:
        double[] phxVegaPGas64 = {
            1.00000000000000005e-04, 1.03180061021383636e-04, 1.05002968345422655e-04,
            1.07871552081389134e-04, 1.12386725052260572e-04, 1.19496313343916537e-04,
            1.30697666356185981e-04, 1.48361992972487070e-04, 1.76258386690715085e-04,
            2.20411975250140330e-04, 2.90535998351765538e-04, 4.02482306184926818e-04,
            5.82567140052077013e-04, 8.75483847795477193e-04, 1.35932010786618583e-03,
            2.17504753883526704e-03, 3.58604002686978935e-03, 6.09998771718679896e-03,
            1.07169005213047595e-02, 1.94122344665836263e-02, 3.59841794245942953e-02,
            6.72855684478238375e-02, 1.24554375514099314e-01, 2.24101975995654207e-01,
            3.86646959394992218e-01, 6.36529574106352247e-01, 1.00343700649801648e+00,
            1.52697716962899532e+00, 2.26349996201742520e+00, 3.29629279099446704e+00,
            4.75030546972967382e+00, 6.80809320973782395e+00, 9.71920175580820889e+00,
            1.38039294970946660e+01, 1.94633310169022842e+01, 2.72009206319020187e+01,
            3.76537976616197625e+01, 5.16414528869924041e+01, 7.02366160217307538e+01,
            9.48405015720874474e+01, 1.27226394963505882e+02, 1.69478952948159872e+02,
            2.23823762506721806e+02, 2.92345015096225097e+02, 3.76222048844978588e+02,
            4.74870843243441925e+02, 5.85652914503416241e+02, 7.04296132837323171e+02,
            8.25563282958302466e+02, 9.46928267959175855e+02, 1.07335181478674599e+03,
            1.22086252441237139e+03, 1.38209820348047128e+03, 1.42038354608713939e+03,
            1.70862293207878497e+03, 2.14866236364070392e+03, 2.81854335114876267e+03,
            3.79724473416333012e+03, 5.19324230258075386e+03, 7.13948492727422672e+03,
            9.83704209840320618e+03, 1.35943082419561761e+04, 1.88366045459964553e+04,
            2.62524123256841995e+04
        };

        int numPhxDeps = phxVegaPGas64.length;  //yeah, I know, 64, but that could change!
        double[] logPhxVegaPGas64 = new double[numPhxDeps];
        for (int i = 0; i < phxVegaPGas64.length; i++) {
            logPhxVegaPGas64[i] = Math.log(phxVegaPGas64[i]);
        }

        // interpolate onto gS3 tauRos grid and re-scale with Teff:
        double[] phxVegaPGas = new double[numDeps];
        double[] logPhxVegaPGas = new double[numDeps];
        double[][] scalePGas = new double[2][numDeps];
        for (int i = 0; i < numDeps; i++) {
            logPhxVegaPGas[i] = Interpol.interpol(ScaleVega.logPhxVegaTau64(), logPhxVegaPGas64, tauRos[1][i]);
            scalePGas[1][i] = logEg + logPhxVegaPGas[i] - ScaleVega.phxVegaLogEg;
            scalePGas[0][i] = Math.exp(scalePGas[1][i]);
            //System.out.println("scalePGas[1][i] " + logE * scalePGas[1][i]);
        }

        return scalePGas;

    }

    public static double[][] phxVegaNe(double grav, int numDeps, double[][] tauRos, double[][] scaleTemp, double kappaScale) {

        double logE = Math.log10(Math.E);

        double logEg = Math.log(grav); //base e!
        double logEkappaScale = Math.log(kappaScale);

        //Theoretical radiative/convective model from Phoenix V15:
        double[] phxVegaPe64 = {
            4.72002072485995232e-05, 4.86860599791532979e-05, 4.95374579675847434e-05,
            5.08766509824373813e-05, 5.29830852602723504e-05, 5.62962571928286436e-05,
            6.15073327430559335e-05, 6.97031622460308396e-05, 8.25926107249302883e-05,
            1.02862199286890015e-04, 1.34734927270846358e-04, 1.84848712459978671e-04,
            2.63639186140655435e-04, 3.87523342400134094e-04, 5.82367141737887679e-04,
            8.89037494764373859e-04, 1.37240124998838403e-03, 2.13631397713436293e-03,
            3.34991041307789979e-03, 5.29691302460575959e-03, 8.47193682504659984e-03,
            1.37503269072523880e-02, 2.26494500083962054e-02, 3.77313887778134710e-02,
            6.30573524863389245e-02, 1.04294141537741802e-01, 1.68555447892395988e-01,
            2.64385476659441288e-01, 4.01678991133522623e-01, 5.91536883461206586e-01,
            8.46971246490695551e-01, 1.18620359311333967e+00, 1.63751721597473354e+00,
            2.24004471648700143e+00, 3.04197315912985289e+00, 4.10229233843691699e+00,
            5.49478489526292702e+00, 7.30946159460564449e+00, 9.65851834249412100e+00,
            1.26929321198114806e+01, 1.66335694784246328e+01, 2.18272692034725146e+01,
            2.87966183424124722e+01, 3.83435065340652415e+01, 5.18840423796927794e+01,
            7.16147842418639584e+01, 1.00791921801285326e+02, 1.44113729824036938e+02,
            2.08146410010596099e+02, 2.96078444680182940e+02, 4.02722746098501375e+02,
            5.18794407567871303e+02, 6.06231844915437705e+02, 6.38241492986035837e+02,
            7.83101211887763725e+02, 1.01011386935192570e+03, 1.33500379236050821e+03,
            1.81044511775468641e+03, 2.49791613029077826e+03, 3.47361550203666138e+03,
            4.84059855334337863e+03, 6.73990325606538136e+03, 9.37502940878883783e+03,
            1.30918950254820193e+04
        };

        int numPhxDeps = phxVegaPe64.length;  //yeah, I know, 64, but that could change!
        double[] logPhxVegaPe64 = new double[numPhxDeps];
        for (int i = 0; i < phxVegaPe64.length; i++) {
            logPhxVegaPe64[i] = Math.log(phxVegaPe64[i]);
        }

        // interpolate onto gS3 tauRos grid and re-scale with Teff:
        double[] phxVegaPe = new double[numDeps];
        double[] logPhxVegaPe = new double[numDeps];
        double[] logScalePe = new double[numDeps];
        double[][] scaleNe = new double[2][numDeps];

        for (int i = 0; i < numDeps; i++) {
            logPhxVegaPe[i] = Interpol.interpol(ScaleVega.logPhxVegaTau64(), logPhxVegaPe64, tauRos[1][i]);
            logScalePe[i] = logEg + ScaleVega.phxVegaLogEkappaScale + logPhxVegaPe[i] - ScaleVega.phxVegaLogEg - logEkappaScale;
            scaleNe[1][i] = logScalePe[i] - scaleTemp[1][i] - Useful.logK();
            scaleNe[0][i] = Math.exp(scaleNe[1][i]);
            //System.out.println("scaleNe[1][i] " + logE * scaleNe[1][i]);
        }

        return scaleNe;

    }

    //Try to recover the opacity as lambda_0 = 1200 nm:
    public static double[][] phxVegaKappa(int numDeps, double[][] tauRos, double kappaScale) {

        double logEkappaScale = Math.log(kappaScale);

        //Theoretical radiative/convective model from Phoenix V15:
        double[] phxVegaRho64 = {
            1.37223023525555033e-16, 1.41626154144729084e-16, 1.44150314678269172e-16,
            1.48123868450381874e-16, 1.54381872696320795e-16, 1.64244615143250253e-16,
            1.79805591659528720e-16, 2.04399072972825006e-16, 2.43371389381634857e-16,
            3.05381535266035120e-16, 4.04658625218062506e-16, 5.65060079088544924e-16,
            8.27656644735441331e-16, 1.26547228507548020e-15, 2.01314207075828015e-15,
            3.32823827759883533e-15, 5.72029150223031827e-15, 1.02224805472660837e-14,
            1.89501093450353969e-14, 3.61761582879661581e-14, 7.01351717400695016e-14,
            1.35418361851359771e-13, 2.55025155983697073e-13, 4.59971412683176497e-13,
            7.85290223156380020e-13, 1.26757212771981250e-12, 1.95008670419967345e-12,
            2.89368402824842727e-12, 4.19202236851314695e-12, 5.99359837525162254e-12,
            8.53008593628702862e-12, 1.21376785175710081e-11, 1.72572973426857710e-11,
            2.44343113771720162e-11, 3.43400746676555601e-11, 4.78069147968408779e-11,
            6.58740603564234330e-11, 8.98765392638669218e-11, 1.21560637642731269e-10,
            1.63165474617735588e-10, 2.17382402795025265e-10, 2.87014403511926534e-10,
            3.74453897536819904e-10, 4.80819443523762100e-10, 6.03543249187692372e-10,
            7.34746346150234889e-10, 8.60873424351652219e-10, 9.63722272464669253e-10,
            1.02170752697610886e-09, 1.02830142548489420e-09, 1.00463534341353045e-09,
            9.89790772818917655e-10, 1.06177003365404541e-09, 1.04384077274744122e-09,
            1.18730551703367904e-09, 1.33892557911181482e-09, 1.63509907922796515e-09,
            2.04571150388433141e-09, 2.59448905764241340e-09, 3.30151864821760360e-09,
            4.19852739489000737e-09, 5.36614157540176651e-09, 6.89691211725147016e-09,
            8.86141299252599472e-09
        };

        double[] phxVegaRadius64 = {
            1.67000000000000000e+11, 1.66997434824341736e+11, 1.66996000021148224e+11,
            1.66993792340530334e+11, 1.66990434901249573e+11, 1.66985415593276306e+11,
            1.66978091504010590e+11, 1.66967747703703705e+11, 1.66953729005060303e+11,
            1.66935618940256409e+11, 1.66913380644500641e+11, 1.66887368804995331e+11,
            1.66858204504870911e+11, 1.66826596519038635e+11, 1.66793201315707397e+11,
            1.66758557698015289e+11, 1.66723082439006622e+11, 1.66687100184744232e+11,
            1.66650891192972473e+11, 1.66614752631880219e+11, 1.66579066154884308e+11,
            1.66544321171909760e+11, 1.66511040302863800e+11, 1.66479661006683319e+11,
            1.66450420200505585e+11, 1.66423258626350128e+11, 1.66397856575875244e+11,
            1.66373783121292786e+11, 1.66350617201162903e+11, 1.66327999415688385e+11,
            1.66305658414097168e+11, 1.66283448699324646e+11, 1.66261378401786652e+11,
            1.66239566586469635e+11, 1.66218148050891907e+11, 1.66197208457392120e+11,
            1.66176772400093903e+11, 1.66156811150862274e+11, 1.66137260364144562e+11,
            1.66118049813053711e+11, 1.66099136305695648e+11, 1.66080531728826263e+11,
            1.66062302365836456e+11, 1.66044553510948181e+11, 1.66027449353869537e+11,
            1.66011174236419525e+11, 1.65995846949225647e+11, 1.65981457156168640e+11,
            1.65967813697925201e+11, 1.65954317129624390e+11, 1.65939737291530853e+11,
            1.65922026262855011e+11, 1.65903466968833923e+11, 1.65899161977920959e+11,
            1.65868787846637787e+11, 1.65828275959904175e+11, 1.65776089155915100e+11,
            1.65714855838088715e+11, 1.65645832748185791e+11, 1.65570171656197845e+11,
            1.65487724367720032e+11, 1.65397726452453491e+11, 1.65299665551203156e+11,
            1.65191869291644409e+11
        };

        int numPhxDeps = phxVegaRadius64.length;
        double[] phxVegaKappa64 = new double[numPhxDeps];
        double[] logPhxVegaKappa64 = new double[numPhxDeps];
        //double[] logPhxSunRho64 = new double[numPhxDeps];
        //double[] logPhxSunRadius64 = new double[numPhxDeps];

//Fix to get right depth scale and right line strengths:
// Yeah - everywhere ya go - opacity fudge
        double fudge = 0.02;
        double logFudge = Math.log(fudge);
        double deltaRho, deltaRadius, deltaTau, logDeltaRho, logDeltaRadius, logDeltaTau;
        double logE = Math.log10(Math.E);
        for (int i = 1; i < numPhxDeps; i++) {

            //Renormalize radii before taking difference
            //Caution: Radius *decreases* with increasing i (inward) and we'll be taking the log:
            deltaRadius = (1.0e-12 * phxVegaRadius64[i - 1]) - (1.0e-12 * phxVegaRadius64[i]);
            deltaRadius = Math.abs(deltaRadius);
            //restore to cm:
            deltaRadius = 1.0e12 * deltaRadius;
            //Renormalize before taking rho difference
            deltaRho = (1.0e12 * phxVegaRho64[i]) - (1.0e12 * phxVegaRho64[i - 1]);
            deltaRho = Math.abs(deltaRho);
            //Restore g/cm^3:
            deltaRho = 1.0e-12 * deltaRho;
            //Renormalize before taking rho difference
            deltaTau = (1.0e2 * ScaleVega.phxVegaTau64[i]) - (1.0e2 * ScaleVega.phxVegaTau64[i - 1]);
            deltaTau = Math.abs(deltaTau);
            deltaTau = 1.0e-2 * deltaTau;

            logDeltaRadius = Math.log(deltaRadius);
            logDeltaRho = Math.log(deltaRho);
            logDeltaTau = Math.log(deltaTau);

            logPhxVegaKappa64[i] = logDeltaTau - logDeltaRho - logDeltaRadius - logEkappaScale + logFudge;
            phxVegaKappa64[i] = Math.exp(logPhxVegaKappa64[i]);
            //System.out.println("logPhxSunKappa64[i] " + logE*logPhxSunKappa64[i]);

        }

        logPhxVegaKappa64[0] = logPhxVegaKappa64[1];
        phxVegaKappa64[0] = phxVegaKappa64[1];

        // interpolate onto gS3 tauRos grid and re-scale with Teff:
        double[][] phxVegaKappa = new double[2][numDeps];
        for (int i = 0; i < numDeps; i++) {
            phxVegaKappa[1][i] = Interpol.interpol(ScaleVega.logPhxVegaTau64(), logPhxVegaKappa64, tauRos[1][i]);
            phxVegaKappa[0][i] = Math.exp(phxVegaKappa[1][i]);
            //System.out.println("phxVegaKappa[1][i] " + logE * phxVegaKappa[1][i]);
        }

        return phxVegaKappa;

    }

}
