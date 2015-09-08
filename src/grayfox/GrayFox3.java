/*
 * The openStar project: stellar atmospheres and spectra
 *
 * grayFox3
 * V3.0, May 2015
 * 
 * C. Ian Short
 * Saint Mary's University
 * Department of Astronomy and Physics
 * Institute for Computational Astrophysics (ICA)
 * Halifax, NS, Canada
 *
 * Open source pedagogical computational stellar astrophysics
 *
 * 1D, static, plane-parallel, LTE, multi-gray stellar atmospheric model
 * Voigt spectral line profile
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

import java.text.DecimalFormat;
import javafx.application.Application;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Scene;
import javafx.scene.chart.LineChart;
import javafx.scene.control.Button;
import javafx.scene.control.Label;
import javafx.scene.control.PasswordField;
import javafx.scene.control.TextField;
import javafx.scene.control.Tooltip;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.StackPane;
import javafx.scene.paint.Color;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import javafx.scene.text.Text;
import javafx.stage.Stage;

public class GrayFox3 extends Application {

    @Override
    public void start(Stage primaryStage) {

        //double teff;
        primaryStage.setTitle("GrayFox 3.0 ");

        GridPane grid = new GridPane();
        //grid.setAlignment(Pos.CENTER);
        grid.setHgap(6);
        grid.setVgap(6);
        grid.setPadding(new Insets(15, 15, 15, 15));
        //grid.setAlignment(Pos.TOP_RIGHT);

        Text scenetitle = new Text("grayFox 3.0 (Hover for tool tips)");
        scenetitle.setFont(Font.font("Tahoma", FontWeight.NORMAL, 16));
        grid.add(scenetitle, 0, 0, 2, 1);

        Label mainLbl = new Label("openStar Static 1D plane-parallel LTE multi-gray atmosphere + Voigt profile spectral line)");
        mainLbl.setFont(Font.font("Tahoma", FontWeight.NORMAL, 11));
        grid.add(mainLbl, 2, 0);

        ImageView imageLogo;
        imageLogo = new ImageView(
                new Image(GrayFox3.class.getResourceAsStream("graphics/SMULOGO2.png")));
        grid.add(imageLogo, 4, 0, 2, 1);

        // Model atmosphere parameters:
        Label atmosLbl = new Label("Model atmosphere parameters:");
        atmosLbl.setFont(Font.font("Tahoma", FontWeight.NORMAL, 12));
        grid.add(atmosLbl, 0, 1);
        Tooltip atmosTip = new Tooltip();
        atmosTip.setText(
                "Default parameters are for the Sun, a G2 V dwarf"
        );
        atmosLbl.setTooltip(atmosTip);

        // Teff input:
        Label teffLbl = new Label("Effective temperature (3000 - 40,000 K)");
        //teffLbl.setAlignment(Pos.CENTER_RIGHT);
        grid.add(teffLbl, 0, 2);
        Tooltip teffTip = new Tooltip();
        teffTip.setText(
                "Equivalent solid spherical blackbody surface temperature, in Kelvins"
        );
        teffLbl.setTooltip(teffTip);

        TextField teffIn = new TextField("5800");
        //teffIn.maxWidth(1.0);
        grid.add(teffIn, 1, 2);

        //logg input
        Label loggLbl = new Label("log_10 Surface gravity (1.0 - 5.5 cm/s/s)");
        grid.add(loggLbl, 0, 3);
        Tooltip loggTip = new Tooltip();
        loggTip.setText(
                "log(GM/R^2); M=mass, R=radius \r\n"
                + "1.0: Giant, 5.5: Dwarf"
        );
        loggLbl.setTooltip(loggTip);

        TextField loggIn = new TextField("4.5");
        grid.add(loggIn, 1, 3);

        //kappaScale input
        Label kappaLbl = new Label("Mean opacity (0.1 - 3.0 x solar)");
        grid.add(kappaLbl, 0, 4);
        Tooltip kappaTip = new Tooltip();
        kappaTip.setText(
                "Rosseland mean extinction multiplier - simulates varying metallicity"
        );
        kappaLbl.setTooltip(kappaTip);

        TextField kappaIn = new TextField("1.0");
        grid.add(kappaIn, 1, 4);

        //massStar input
        Label massStarLbl = new Label("Mass (0.25 - 5.0 x M_Sun)");
        grid.add(massStarLbl, 0, 5);
        Tooltip massStarTip = new Tooltip();
        massStarTip.setText(
                "Star mass in solar masses"
        );
        massStarLbl.setTooltip(massStarTip);

        TextField massStarIn = new TextField("1.0");
        grid.add(massStarIn, 1, 5);

        //
        //
        // Spectral line parameters:
        //
        //
        Label lineLbl = new Label("Spectral line parameters:");
        lineLbl.setFont(Font.font("Tahoma", FontWeight.NORMAL, 12));
        grid.add(lineLbl, 2, 1);

        // wavelength (lambda):
        Label lamLbl = new Label("Wavelength (200 - 2000 nm)");
        grid.add(lamLbl, 2, 2);

        TextField lamIn = new TextField("500.0");
        grid.add(lamIn, 3, 2);

        // xiT:
        Label xitLbl = new Label("Microturbulence (1.0 - 8.0 km/s)");
        grid.add(xitLbl, 2, 3);
        Tooltip xitTip = new Tooltip();
        xitTip.setText(
                "RMS velocity of optically thin gas cells\r\n"
                + " - affects amount of Doppler broadening in line core caused by convective motion"
        );
        xitLbl.setTooltip(xitTip);

        TextField xitIn = new TextField("1.0");
        grid.add(xitIn, 3, 3);

        // mass:
        Label massLbl = new Label("Particle mass (1 - 238 amu)");
        grid.add(massLbl, 2, 4);
        Tooltip massTip = new Tooltip();
        massTip.setText(
                "Mass of absorbing particles in atomic mass units (amu)\r\n"
                + "- affects amount of Doppler broadening in line core caused by kinetic thermal motion"
        );
        massLbl.setTooltip(massTip);

        TextField massIn = new TextField("12.0");
        grid.add(massIn, 3, 4);

        // Van der Waals broadening
        Label gamLbl = new Label("log_10 Broadening factor (0.0 - 1.0)");
        grid.add(gamLbl, 2, 5);
        Tooltip gamTip = new Tooltip();
        gamTip.setText(
                "Enhances the amount of pressure (collisional) line broadening \r\n"
                + " - affects the Lorentzian line wings"
        );
        gamLbl.setTooltip(gamTip);

        TextField gamIn = new TextField("0.0");
        grid.add(gamIn, 3, 5);

        // Abundance log n_l:
        Label a12Lbl = new Label("A12 Abundance (3.0 - 9.0)");
        grid.add(a12Lbl, 4, 2);
        Tooltip a12Tip = new Tooltip();
        a12Tip.setText(
                "Abundance absorbing element on A12 scale"
        );
        a12Lbl.setTooltip(a12Tip);

        TextField a12In = new TextField("3.5");
        grid.add(a12In, 5, 2);

        // Oscillator strength -  log f:
        Label logfLbl = new Label("log_10 Oscillator strength (-6.0 - 1.0)");
        grid.add(logfLbl, 4, 3);
        Tooltip logfTip = new Tooltip();
        logfTip.setText(
                "The quantum mechanical probablity for the line transition to absorb a photon (unitless parameter)"
        );
        logfLbl.setTooltip(logfTip);

        TextField logfIn = new TextField("-1.0");
        grid.add(logfIn, 5, 3);

        // Ground state ionization E (eV) - Stage I:
        Label ion1Lbl = new Label("Stage I Ground ion. E (ev) (5.0 - 25.0)");
        grid.add(ion1Lbl, 4, 4);

        TextField ion1In = new TextField("8.0");
        grid.add(ion1In, 5, 4);

        // Lower level excitation E (ev):
        Label excLbl = new Label("Excitation E (ev) (0  - Ion.E - 10+Ion.E)");
        grid.add(excLbl, 4, 5);
        Tooltip excTip = new Tooltip();
        excTip.setText(
                "The excitation energy with respect to the ground state of the lower atomic E-level\r\n"
                + " of the line transition in electron Volts (eV)\r\n"
                + " - affects line strength in a Teff-dependent way\r\n"
                + " ** Values < IonE: Neutral stage, Values > Ion.E.: Ionized stage "
        );
        excLbl.setTooltip(excTip);

        TextField excIn = new TextField("0.0");
        grid.add(excIn, 5, 5);

        //Ground state ionization E (eV) - Stage II: :
        Label ion2Lbl = new Label("Stage II Ground ion. E (ev) (5.0 - 50.0)");
        grid.add(ion2Lbl, 6, 2);

        TextField ion2In = new TextField("8.0");
        grid.add(ion2In, 7, 2);

        //Ground state statistical weight (or partition fn) - Stage I: :
        Label gw1Lbl = new Label("Stage I Ground stat. wght");
        grid.add(gw1Lbl, 6, 3);
        Tooltip gw1Tip = new Tooltip();
        gw1Tip.setText(
                "Statistical weight, g, of ground state of ion stage I\r\n"
                + " OR partition fn of ion stage I (unitless)"
        );
        gw1Lbl.setTooltip(gw1Tip);
        TextField gw1In = new TextField("1");
        grid.add(gw1In, 7, 3);

        //Ground state statistical weight (or partition fn) - Stage II: :
        Label gw2Lbl = new Label("Stage II Ground stat. wght");
        grid.add(gw2Lbl, 6, 4);
        Tooltip gw2Tip = new Tooltip();
        gw2Tip.setText(
                "Statistical weight, g, of ground state of ion stage I\r\n"
                + " OR partition fn of ion stage I (unitless)"
        );
        gw2Lbl.setTooltip(gw2Tip);
        TextField gw2In = new TextField("1");
        grid.add(gw2In, 7, 4);

        //Lower E-level statistical weight:
        Label gwLLbl = new Label("Lower E level stat. wght");
        grid.add(gwLLbl, 6, 5);
        Tooltip gwLTip = new Tooltip();
        gwLTip.setText(
                "Statistical weight, g, of lower E level of b-b transition"
        );
        gwLLbl.setTooltip(gwLTip);
        TextField gwLIn = new TextField("1");
        grid.add(gwLIn, 7, 5);

        //
        //
        //Label kappaRng = new Label("Range: 0.1 - 3.0");
        //grid.add(kappaRng, 2, 4);
        Button btn = new Button("Model");
        HBox hbBtn = new HBox(10);
        hbBtn.setAlignment(Pos.BOTTOM_RIGHT);
        hbBtn.getChildren().add(btn);
        grid.add(hbBtn, 5, 7);

        final Text actiontarget = new Text();
        grid.add(actiontarget, 0, 8);

        ////Rectangle r = new Rectangle(10,600,1590,890);
        //Rectangle r = new Rectangle(0, 0, 50, 50);
        //r.setFill(Color.WHITE);
        //grid.getChildren().add(r);
        // Need to add logic to reset fields?? 
        btn.setOnAction(new EventHandler<ActionEvent>() {

            @Override
            public void handle(ActionEvent e) {
                actiontarget.setFill(Color.FIREBRICK);
                //actiontarget.setText("Output here");
                actiontarget.setFont(Font.font("Tahoma", FontWeight.NORMAL, 14));

                String teffStr = teffIn.getText();
                String loggStr = loggIn.getText();
                String kappaStr = kappaIn.getText();
                String massStarStr = massStarIn.getText();

                String lamStr = lamIn.getText();
                String xitStr = xitIn.getText();
                String massStr = massIn.getText();
                String gamStr = gamIn.getText();
                String a12Str = a12In.getText();
                String fStr = logfIn.getText();
                String ion1Str = ion1In.getText();
                String ion2Str = ion2In.getText();
                String gw1Str = gw1In.getText();
                String gw2Str = gw2In.getText();
                String gwLStr = gwLIn.getText();
                String excStr = excIn.getText();

                // If *any* model atmosphere input string is empty, set teff to empty to simplify 
                // conditional logic below.
                if (loggStr == null || loggStr.isEmpty()) {
                    teffStr = null;
                }
                if (kappaStr == null || kappaStr.isEmpty()) {
                    teffStr = null;
                }
                if (massStarStr == null || massStarStr.isEmpty()) {
                    teffStr = null;
                }
                // If *any* spectral line input string is empty, set teff to empty to simplify 
                // conditional logic below.  
                if (lamStr == null || lamStr.isEmpty()) {
                    teffStr = null;
                }
                if (xitStr == null || xitStr.isEmpty()) {
                    teffStr = null;
                }
                if (massStr == null || massStr.isEmpty()) {
                    teffStr = null;
                }
                if (gamStr == null || gamStr.isEmpty()) {
                    teffStr = null;
                }
                if (a12Str == null || a12Str.isEmpty()) {
                    teffStr = null;
                }
                if (fStr == null || fStr.isEmpty()) {
                    teffStr = null;
                }
                if (ion1Str == null || ion1Str.isEmpty()) {
                    teffStr = null;
                }
                if (ion2Str == null || ion2Str.isEmpty()) {
                    teffStr = null;
                }
                if (gw1Str == null || gw1Str.isEmpty()) {
                    teffStr = null;
                }
                if (gw2Str == null || gw2Str.isEmpty()) {
                    teffStr = null;
                }
                if (gwLStr == null || gwLStr.isEmpty()) {
                    teffStr = null;
                }
                if (excStr == null || excStr.isEmpty()) {
                    teffStr = null;
                }

                if ((teffStr != null && !teffStr.isEmpty())) {

                    // Argument 1: Effective temperature, Teff, in K:
                    double teff = (Double.valueOf(teffStr)).doubleValue();

                    // Argument 2: Logarithmic surface gravity, g, in cm/s/s:
                    double logg = (Double.valueOf(loggStr)).doubleValue();

                    //Argument 3: Linear sclae factor for solar Rosseland oapcity distribution
                    // mimics "metallicity" parameter - ??  (unitless)
                    double kappaScale = (Double.valueOf(kappaStr)).doubleValue();

                    //Argument 4: Stellar mass, M, in solar masses
                    double massStar = (Double.valueOf(massStarStr)).doubleValue();

                    // Sanity check:
                    if (teff < 3000.0) {
                        teff = 3000.0;
                        teffStr = "3000";
                    }
                    if (teff > 50000.0) {
                        teff = 50000.0;
                        teffStr = "50000";
                    }
                    //logg limit is strongly Teff-dependent:
                    double minLogg = 3.5; //safe initialization
                    String minLoggStr = "3.5";
                    if (teff <= 4000.0) {
                        minLogg = 1.5;
                        minLoggStr = "1.5";
                    } else if ((teff > 4000.0) && (teff <= 5500.0)) {
                        minLogg = 2.0;
                        minLoggStr = "2.0";
                    } else if ((teff > 5000.0) && (teff <= 7000.0)) {
                        minLogg = 2.5;
                        minLoggStr = "2.5";
                    } else if ((teff > 7000.0) && (teff < 9000.0)) {
                        minLogg = 3.0;
                        minLoggStr = "3.0";
                    } else if (teff >= 9000.0) {
                        minLogg = 3.5;
                        minLoggStr = "3.5";
                    }
                    if (logg < minLogg) {
                        logg = minLogg;
                        loggStr = minLoggStr;
                    }
                    if (logg > 8.5) {
                        logg = 8.5;
                        loggStr = "8.5";
                    }
                    if (kappaScale < 0.1) {
                        kappaScale = 0.1;
                        kappaStr = "0.1";
                    }
                    if (kappaScale > 3.0) {
                        kappaScale = 3.0;
                        kappaStr = "3.0";
                    }
                    if (massStar < 0.25) {
                        massStar = 0.25;
                        massStarStr = "0.25";
                    }
                    if (massStar > 5.0) {
                        massStar = 5.0;
                        massStarStr = "5.0";
                    }

                    double grav = Math.pow(10.0, logg);

                    // Argument 4: Wavelength, lambda, in nm:
                    double lam0 = (Double.valueOf(lamStr)).doubleValue();

                    // Argument 5: microturbulence, xi_T, in km/s:
                    double xiT = (Double.valueOf(xitStr)).doubleValue();

                    // Argument 6: Mass, m, in atomic mass units (amu):
                    double mass = (Double.valueOf(massStr)).doubleValue();

                    // Argument 7: linear Van de Waals broadening enhancement factor, "gamma" (unitless):
                    double logGammaCol = (Double.valueOf(gamStr)).doubleValue();

                    // Argument 8:  lower E-level num density  (log(N_l): 
                    //N_l in cm^-3
                    double A12 = (Double.valueOf(a12Str)).doubleValue();

                    // Argument 9:  oscillator strength (log f): 
                    //unitless
                    double logf = (Double.valueOf(fStr)).doubleValue();

                    //// Argument 10: ground state ionization E (eV) - Stage I
                    double chiI1 = (Double.valueOf(ion1Str)).doubleValue();

                    // Argument 11: lower level excitation E (eV)
                    double chiL = (Double.valueOf(excStr)).doubleValue();

                    //// Argument 12: ground state ionization E (eV) - Stage II
                    double chiI2 = (Double.valueOf(ion2Str)).doubleValue();

                    //// Argument 13: ground state statistical weight (or partition fn) - Stage I
                    double gw1 = (Double.valueOf(gw1Str)).doubleValue();

                    //// Argument 14: ground state statistical weight (or partition fn) - Stage II
                    double gw2 = (Double.valueOf(gw2Str)).doubleValue();

                    //// Argument 15: ground state statistical weight of lower E-level of b-b transition
                    double gwL = (Double.valueOf(gwLStr)).doubleValue();

                    if (lam0 < 200.0) {
                        lam0 = 200.0;
                        lamStr = "200";
                    }
                    if (lam0 > 20000.0) {
                        lam0 = 20000.0;
                        lamStr = "20000";
                    }
                    if (xiT < 0.0) {
                        xiT = 0.0;
                        xitStr = "0.0";
                    }
                    if (xiT > 8.0) {
                        xiT = 8.0;
                        xitStr = "8.0";
                    }
                    if (mass < 1.0) {
                        mass = 1.0;
                        massStr = "1.0";
                    }
                    if (mass > 238.0) {
                        mass = 238.0;
                        massStr = "238";
                    }
                    if (logGammaCol < 0.0) {
                        logGammaCol = 0.0;
                        gamStr = "0.0";
                    }
                    if (logGammaCol > 1.0) {
                        logGammaCol = 1.0;
                        gamStr = "1.0";
                    }
                    if (A12 < 3.0) {
                        A12 = 3.0;
                        a12Str = "3.0";
                    }
                    if (A12 > 9.0) {
                        A12 = 9.0;
                        a12Str = "9.0";
                    }
                    if (logf < -6.0) {
                        logf = -6.0;
                        fStr = "-6.0";
                    }
                    if (logf > 1.0) {
                        logf = 1.0;
                        fStr = "1.0";
                    }
                    if (chiI1 < 5.0) {
                        chiI1 = 5.0;
                        ion1Str = "5.0";
                    }
                    if (chiI1 > 25.0) {
                        chiI1 = 25.0;
                        ion1Str = "25.0";
                    }
                    if (chiI2 < 5.0) {
                        chiI2 = 5.0;
                        ion2Str = "5.0";
                    }
                    if (chiI2 > 50.0) {
                        chiI2 = 50.0;
                        ion2Str = "50.0";
                    }
                    if (gw1 < 1.0) {
                        gw1 = 1.0;
                        gw1Str = "1";
                    }
                    if (gw1 > 18.0) {
                        gw1 = 18.0;
                        gw1Str = "18";
                    }
                    if (gw2 < 1.0) {
                        gw2 = 1.0;
                        gw2Str = "1";
                    }
                    if (gw2 > 100.0) {
                        gw2 = 100.0;
                        gw2Str = "100";
                    }
                    if (gwL < 1.0) {
                        gwL = 1.0;
                        gwLStr = "1";
                    }
                    if (gwL > 100.0) {
                        gwL = 100.0;
                        gwLStr = "100";
                    }

                    boolean ionized = false;
                    // Note: Upper limit of chiL depends on value of chiI above!
                    if (chiL < 0.0) {
                        chiL = 0.0;   // Ground state case!
                        excStr = "0.0";
                    }
                    // choice of neutral or singly ionized stage
                    if (chiL < chiI1) {
                        ionized = false;
                    }
                    if (chiL == chiI1) {
                        ionized = false;
                        chiL = 0.99 * chiI1;
                        excStr = ion1Str;
                    }
                    if (chiL > chiI1) {
                        ionized = true;
                        // chiL = chiL - chiI; Not here!
                        if (chiL > chiI1 + 50.0) {
                            chiL = 1 + 50.0;
                            excStr = " " + chiL;
                        }
                    }

                    //   if (chiL > 0.99 * chiI) {
                    //       chiL = 0.99 * chiI;
                    //      excStr = ionStr;
                    //}
                    //if (chiL > 10.0) {
                    //    chiL = 10.0;
                    //   excStr = "10.0";
                    //}
                    // debug
                    /*actiontarget.setText(
                     "Stellar atmosphere parameters: \r\n"
                     + "T_eff:  " + teffStr + " K,  "
                     + "log(g):  " + loggStr + " cm/s/s,  "
                     + "Kappa_Ros:  " + kappaStr + " x solar\r\n\r\n"
                     + "Spectral line parameters: \r\n"
                     + "lambda:  " + lamStr + " nm,  "
                     + "xi_T:  " + xitStr + " km/s,  "
                     + "mass:  " + massStr + " amu\r\n"
                     + "Gamma:  " + gamStr + ", log(N): " + nStr + " cm^-3\r\n"
                     + "log(f): " + fStr + "excit E_l: " + excStr + " (eV)");
                     //+ "grnd ioniz E: " + ionStr + " eV " + "excit E_l: " + excStr + " (eV)");
                     //debug
                     */
                    //Gray structure and Voigt line code code begins here:
// Initial set-up:
                    // optical depth grid
                    int numDeps = 48;
                    double log10MinDepth = -6.0;
                    double log10MaxDepth = 2.0;
                    //int numThetas = 10;  // Guess

                    //wavelength grid (cm):
                    double[] lamSetup = new double[3];
                    /*
                     lamSetup[0] = 100.0 * 1.0e-7;  // Start wavelength, cm
                     lamSetup[1] = 1000.0 * 1.0e-7; // End wavelength, cm
                     lamSetup[2] = 40;  // number of lambda
                     */
                    lamSetup[0] = 100.0 * 1.0e-7;  // test Start wavelength, cm
                    lamSetup[1] = 2000.0 * 1.0e-7; // test End wavelength, cm
                    lamSetup[2] = 40;  // test number of lambda
                    lam0 = lam0 * 1.0e-7;  // line centre lambda from nm to cm
                    //int numLams = (int) (( lamSetup[1] - lamSetup[0] ) / lamSetup[2]) + 1;  
                    int numLams = (int) lamSetup[2];

// Solar parameters:
                    double teffSun = 5778.0;
                    double loggSun = 4.44;
                    double gravSun = Math.pow(10.0, loggSun);
                    double kappaScaleSun = 1.0;
//Solar units:
                    double massSun = 1.0;
                    double radiusSun = 1.0;
                    //double massStar = 1.0; //solar masses // test
                    double logRadius = 0.5 * (Math.log(massStar) + Math.log(gravSun) - Math.log(grav));
                    double radius = Math.exp(logRadius); //solar radii
                    //double radius = Math.sqrt(massStar * gravSun / grav); // solar radii
                    double logLum = 2.0 * Math.log(radius) + 4.0 * Math.log(teff / teffSun);
                    double bolLum = Math.exp(logLum); // L_Bol in solar luminosities 

                    //Composition by mass fraction - needed for opacity approximations
                    //   and interior structure
                    double massX = 0.70; //Hydrogen
                    double massY = 0.28; //Helium
                    double massZSun = 0.02; // "metals"
                    double massZ = massZSun * kappaScale; //approximation

                    double logNH = 17.0;
                    double logN = (A12 - 12.0) + logNH;

                    //Output files:
                    String outfile = "gray_structure."
                            + teffStr + "-" + loggStr + "-" + kappaStr + ".out";
                    String specFile = "gray_spectrum."
                            + teffStr + "-" + loggStr + "-" + kappaStr + ".out";
                    String lineFile = "voigt_line."
                            + teffStr + "-" + loggStr + "-" + kappaStr + "-" + xitStr + ".out";

                    double logE = Math.log10(Math.E); // for debug output

                    //log_10 Rosseland optical depth scale  
                    double tauRos[][] = TauScale.tauScale(numDeps, log10MinDepth, log10MaxDepth);

                    //Gray kinetic temeprature structure:
                    double temp[][] = Temperature.temperature(numDeps, teff, tauRos);

                    //Now do the same for the Sun, for reference:
                    double tempSun[][] = Temperature.temperature(numDeps, teffSun, tauRos);

                    //
                    // BEGIN Initial guess for Sun section:
                    //
                    //Data amalgamated for several stars from Table 9.2, Observation and Analysis of Stellar Photospheres, 3rd Ed.,
                    // David F. Gray ("Dfg")
                    //** CAUTION: last two values in list are for logg=4.0, first 4 are for logg=4.6, and rest are for solar logg
                    double[] tempDfg = {3017.0, 3111.0, 3262.0, 3592.0, 4310.0, 4325.0, 4345.0, 4370.0, 4405.0, 4445.0, 4488.0, 4524.0, 4561.0, 4608.0, 4660.0, 4720.0, 4800.0, 4878.0, 4995.0, 5132.0, 5294.0, 5490.0, 5733.0, 6043.0, 6429.0, 6904.0, 7467.0, 7962.0, 8358.0, 8630.0, 8811.0, 9643.0, 12945.0};
                    double[] log10PgDfg = {3.22, 3.89, 4.45, 5.00, 2.87, 3.03, 3.17, 3.29, 3.41, 3.52, 3.64, 3.75, 3.86, 3.97, 4.08, 4.19, 4.30, 4.41, 4.52, 4.63, 4.74, 4.85, 4.95, 5.03, 5.10, 5.15, 5.18, 5.21, 5.23, 5.26, 5.29, 3.76, 3.88};
                    double[] log10PeDfg = {-2.12, -1.51, -0.95, -0.26, -1.16, -1.02, -0.89, -0.78, -0.66, -0.55, -0.44, -0.33, -0.23, -0.12, -0.01, 0.10, 0.22, 0.34, 0.47, 0.61, 0.76, 0.93, 1.15, 1.43, 1.78, 2.18, 2.59, 2.92, 3.16, 3.32, 3.42, 2.96, 3.43};
                    double[] log10KapOverPeDfg = {-0.46, -0.53, -0.64, -0.85, -1.22, -1.23, -1.24, -1.25, -1.26, -1.28, -1.30, -1.32, -1.33, -1.35, -1.37, -1.40, -1.43, -1.46, -1.50, -1.55, -1.60, -1.66, -1.73, -1.81, -1.91, -2.01, -2.11, -2.18, -2.23, -2.25, -2.27, -1.82, -1.73};

                    int numDfg = tempDfg.length;
                    double[] log10KapDfg = new double[numDfg];
                    double[] logKapDfg = new double[numDfg];
                    double[] logPgDfg = new double[numDfg];
                    double[] logPeDfg = new double[numDfg];
                    double pgDfg;

                    for (int i = 0; i < numDfg; i++) {
                        //Rescale pressures to logg=4.44; assume logP scales with logg through HSE
                        if (i <= 3) {
                            log10PgDfg[i] = log10PgDfg[i] - 0.2;
                            log10PeDfg[i] = log10PeDfg[i] - 0.2;
                        }
                        if (i >= numDfg - 2) {
                            log10PgDfg[i] = log10PgDfg[i] + 0.44;
                            log10PeDfg[i] = log10PeDfg[i] + 0.44;
                        }
                        log10KapDfg[i] = log10KapOverPeDfg[i] + log10PeDfg[i];
                        logKapDfg[i] = Math.log(Math.pow(10.0, log10KapDfg[i]));

                        //Dress up DFG temp and pressure to look like what State.massDensity expects...
                        logPgDfg[i] = Math.log(Math.pow(10.0, log10PgDfg[i]));
                        logPeDfg[i] = Math.log(Math.pow(10.0, log10PeDfg[i]));
                    }

                    //Interpolate DFG data onto our Gray Teff structure:
                    double[][] kapDfg2 = new double[2][numDeps];
                    double[][] rhoDfg2 = new double[2][numDeps];
                    //double[][] tempDfg2 = new double[2][numDeps];
                    double[][] pressDfg2 = new double[4][numDeps];
                    double[] logPeDfg2 = new double[numDeps];
                    double[][] NeDfg2 = new double[2][numDeps];

                    //Prepare simple temperature vector for input to interpol():
                    double[] tempSimp = new double[numDeps];
                    for (int i = 0; i < numDeps; i++) {
                        tempSimp[i] = tempSun[0][i];
                    }

                    //System.out.println("tauRos[1][i]   temp[0][i]   logE*logKapDfg2[i]   logE*pressDfg2[1][i]   logE*logPeDfg2[i]");
                    for (int i = 0; i < numDeps; i++) {

                        if (tempSimp[i] <= tempDfg[0]) {
                            kapDfg2[1][i] = logKapDfg[0];
                            pressDfg2[1][i] = logPgDfg[0];
                            logPeDfg2[i] = logPeDfg[0];
                        } else if (tempSimp[i] >= tempDfg[numDfg - 1]) {
                            kapDfg2[1][i] = logKapDfg[numDfg - 1];
                            pressDfg2[1][i] = logPgDfg[numDfg - 1];
                            logPeDfg2[i] = logPeDfg[numDfg - 1];
                        } else {
                            kapDfg2[1][i] = Interpol.interpol(tempDfg, logKapDfg, tempSimp[i]);
                            pressDfg2[1][i] = Interpol.interpol(tempDfg, logPgDfg, tempSimp[i]);
                            logPeDfg2[i] = Interpol.interpol(tempDfg, logPeDfg, tempSimp[i]);
                        }

                        kapDfg2[0][i] = Math.exp(kapDfg2[1][i]);
                        pressDfg2[0][i] = Math.exp(pressDfg2[1][i]);
                        pressDfg2[2][i] = 0.0;
                        pressDfg2[3][i] = 0.0;

                        //Electron number density, Ne:
                        NeDfg2[1][i] = logPeDfg2[i] - tempSun[1][i] - Useful.logK();
                        NeDfg2[0][i] = Math.exp(NeDfg2[1][i]);

                        //System.out.format("%12.8f   %12.8f   %12.8f   %12.8f   %12.8f%n", logE * tauRos[1][i], tempSun[0][i], logE * kapDfg2[1][i], logE * pressDfg2[1][i], logE * logPeDfg2[i]);
                    }

                    //
                    // END initial guess for Sun section
                    //
                    // *********************
                    //
                    // mean molecular weight and Ne for Star & Sun
                    double[] mmw = State.mmwFn(numDeps, temp, kappaScale);
                    double[] mmwSun = State.mmwFn(numDeps, tempSun, kappaScaleSun);
                    double[][] Ne = State.NeFn(numDeps, temp, NeDfg2, kappaScale);

                    rhoDfg2 = State.massDensity(numDeps, tempSun, pressDfg2, mmw, kappaScale);

                    // Create kappa structure here: Initialize solar kappa_Ross structure and
                    // scale it by logg, radius, and kappaScale:
                    // We do not know rho yet, so cannot compute kappa "in situ" - also in situ does not account for b-b oapcity scaling
                    //double kappa[][] = Kappas.kappas(numDeps, kappaScale, teff, teffSun, logg, loggSun);
                    /* Needed now that we have DFG guess?
                     double[][] rho = new double[2][numDeps];
                     double[][] rhoSun = new double[2][numDeps];
                     //initialize rho arrays
                     for (int i = 0; i < numDeps; i++) {
                     rho[0][i] = 0.0;
                     rho[1][i] = 0.0;
                     rhoSun[0][i] = 0.0;
                     rhoSun[1][i] = 0.0;
                     }
                     */
                    //Get H I n=2 & n=3 number densities for Balmer and Pashen continuum  for kappa calculation
                    // Paschen:
                    boolean ionizedHI = false;
                    double chiI1H = 13.6; //eV
                    double chiI2H = 1.0e6;  //eV //H has no third ionization stage!
                    double gw1H = 2.0;
                    double gw2H = 1.0;  // umm... doesn't exist - no "HIII"
                    // n=3 level - Paschen jump
                    double lamJump3 = 820.4 * 1.0e-7; //Paschen jump - cm
                    double chiLH3 = 12.1; //eV
                    double gwLH3 = 2 * 3 * 3; // 2n^2
                    double logNumsH3[][];
                    // n=2 level - Balmer jump
                    double lamJump2 = 364.0 * 1.0e-7; //Paschen jump - cm
                    double chiLH2 = 10.2; //eV
                    double gwLH2 = 2 * 2 * 2; // 2n^2   
                    double logNumsH2[][];
                    int mode;
                    /* Needed now that we have DFG guess?
                     mode = 0;  //call kappas without knowledge of rho
                     double kappa[][] = Kappas.kappas(mode, numDeps, rho, rhoSun, kappaScale, logg, loggSun, teff, teffSun, radius, massX, massZ, tauRos, temp, tempSun);
                     //double kappaSun[][] = Kappas.kappas(numDeps, kappaScaleSun, teffSun, teffSun, loggSun, loggSun);
                     double kappaSun[][] = Kappas.kappas(mode, numDeps, rho, rhoSun, kappaScaleSun, loggSun, loggSun, teffSun, teffSun, radiusSun, massX, massZSun, tauRos, tempSun, tempSun);
                     */
                    mode = 1;  //call kappas without knowledge of rho

                    logNumsH3 = LevelPops.levelPops(lamJump3, logNH, Ne, ionizedHI, chiI1H, chiI2H, chiLH3, gw1H, gw2H, gwLH3,
                            numDeps, kappaScale, tauRos, temp, rhoDfg2);
                    logNumsH2 = LevelPops.levelPops(lamJump2, logNH, Ne, ionizedHI, chiI1H, chiI2H, chiLH2, gw1H, gw2H, gwLH2,
                            numDeps, kappaScale, tauRos, temp, rhoDfg2);
                    System.out.println("logNumsH2 " + logE*logNumsH2[2][36] + " logNumsH3 " + logE*logNumsH3[2][36]);
                    double kappa[][] = Kappas.kappas(mode, numDeps, rhoDfg2, rhoDfg2, kapDfg2, kappaScale, logg, loggSun, teff, teffSun, radius, massX, massZ, tauRos, temp, tempSun, logNumsH3, logNumsH2);
                    //double kappaSun[][] = Kappas.kappas(numDeps, kappaScaleSun, teffSun, teffSun, loggSun, loggSun);
                    logNumsH3 = LevelPops.levelPops(lamJump3, logNH, Ne, ionizedHI, chiI1H, chiI2H, chiLH3, gw1H, gw2H, gwLH3,
                            numDeps, kappaScale, tauRos, tempSun, rhoDfg2);
                    logNumsH2 = LevelPops.levelPops(lamJump2, logNH, Ne, ionizedHI, chiI1H, chiI2H, chiLH2, gw1H, gw2H, gwLH2,
                            numDeps, kappaScale, tauRos, tempSun, rhoDfg2);                    
                    double kappaSun[][] = Kappas.kappas(mode, numDeps, rhoDfg2, rhoDfg2, kapDfg2, kappaScaleSun, loggSun, loggSun, teffSun, teffSun, radiusSun, massX, massZSun, tauRos, tempSun, tempSun, logNumsH3, logNumsH2);

                    //Next solve hydrostatic eq for P scale on the tau scale - need to pick a depth dependent kappa value!
                    //   - scale kapp_Ross with log(g) from solar value? - Kramers opacity law?
                    //   - dP/dTau = g/kappa
                    //press is a 4 x numDeps array:
                    // rows 0 & 1 are linear and log *gas* pressure, respectively
                    // rows 2 & 3 are linear and log *radiation* pressure
                    double press[][] = Hydrostat.hydrostatic(numDeps, grav, tauRos, kappa, temp);

                    // Then solve eos for the rho scale - need to pick a mean molecular weight, mu
                    double[][] rho = State.massDensity(numDeps, temp, press, mmw, kappaScale);

                    //Now do the same for the Sun, for reference:
                    double pressSun[][] = Hydrostat.hydrostatic(numDeps, gravSun, tauRos, kappaSun, tempSun);
                    double[][] rhoSun = State.massDensity(numDeps, tempSun, pressSun, mmwSun, kappaScaleSun);
                    //double depthsSun[] = DepthScale.depthScale(numDeps, tauRos, kappaSun, rhoSun);
                    // Special one-time print-out of Sun's structure:
                    // System.out.println("Sun: " + " i " + " temp " + " kappa " + " pressGas " + " pressRad " + " rho " + " depths ");
                    //for (int k = 0; k < numDeps; k++) {
                    //    //System.out.println(" " + k + " " + temp[0][k] + " " + kappa[0][k] + " " + press[0][k] + " " + press[2][k] + " " + rho[0][k] + " " + depths[k]);
                    //    System.out.println(" " + k + " " + temp[1][k] + " " + kappa[1][k] + " " + press[1][k] + " " + press[3][k] + " " + rho[1][k]);
                    // }
                    // Limb darkening:
                    // Establish wavelength grid:

                    //compute kappas again with in situ densities thsi time:
                    mode = 1;  //call kappas ** with ** knowledge of rho
                    logNumsH3 = LevelPops.levelPops(lamJump3, logNH, Ne, ionizedHI, chiI1H, chiI2H, chiLH3, gw1H, gw2H, gwLH3,
                            numDeps, kappaScale, tauRos, temp, rho);
                    logNumsH2 = LevelPops.levelPops(lamJump2, logNH, Ne, ionizedHI, chiI1H, chiI2H, chiLH2, gw1H, gw2H, gwLH2,
                            numDeps, kappaScale, tauRos, temp, rho);                    
                    kappa = Kappas.kappas(mode, numDeps, rho, rhoSun, kapDfg2, kappaScale, logg, loggSun, teff, teffSun, radius, massX, massZ, tauRos, temp, tempSun, logNumsH3, logNumsH2);
                    //double kappaSun[][] = Kappas.kappas(numDeps, kappaScaleSun, teffSun, teffSun, loggSun, loggSun);
                    //kappaSun = Kappas.kappas(mode, numDeps, rho, rhoSun, kappaScaleSun, loggSun, loggSun, teffSun, teffSun, radiusSun, massX, massZSun, tauRos, tempSun, tempSun);
                    // Then construct geometric depth scale from tau, kappa and rho
                    double depths[] = DepthScale.depthScale(numDeps, tauRos, kappa, rho);

                    double newTemp[][] = new double[2][numDeps];
                    int numTCorr = 10;  //test
                    //int numTCorr = 0;
                    for (int i = 0; i < numTCorr; i++) {
                        //newTemp = TCorr.tCorr(numDeps, tauRos, temp);
                        newTemp = MulGrayTCorr.mgTCorr(numDeps, teff, tauRos, temp, rho, kappa);
                        for (int iTau = 0; iTau < numDeps; iTau++) {
                            //System.out.format("%12.8f   %12.8f%n", logE*tauRos[1][iTau], temp[0][iTau]);
                            temp[1][iTau] = newTemp[1][iTau];
                            temp[0][iTau] = newTemp[0][iTau];
                            //System.out.format("%12.8f   %12.8f%n", logE*tauRos[1][iTau], temp[0][iTau]);
                            //System.out.println("TCorr 1: iTau: " + iTau + " deltaT " + (newTemp[0][iTau] - temp[0][iTau]));
                        }
                    }
                    //  System.out.println("tauRos[1][iTau]   temp[0][iTau]   temp[0][iTau]-grayTemp[0][iTau]");
                    //  for (int iTau = 0; iTau < numDeps; iTau++) {
                    //      System.out.format("%12.8f   %12.8f   %12.8f%n", logE*tauRos[1][iTau], temp[0][iTau], temp[0][iTau]-grayTemp[0][iTau]);
                    //  }

                    /*
                     //Convection:
                     // Teff below which stars are convective.  
                     //  - has to be finessed because Convec.convec() does not work well :-(
                     double convTeff = 6500.0;
                     double[][] convTemp = new double[2][numDeps];
                     if (teff < convTeff) {
                     convTemp = Convec.convec(numDeps, tauRos, depths, temp, press, rho, kappa, kappaSun, kappaScale, teff, logg);

                     for (int iTau = 0; iTau < numDeps; iTau++) {
                     temp[1][iTau] = convTemp[1][iTau];
                     temp[0][iTau] = convTemp[0][iTau];
                     }
                     }
                     */
                    /*
                     //Recall hydrostat with updates temps            
                     //Recall state withupdated Press                    
                     //recall kappas withupdates rhos
                     //Recall depths with re-updated kappas
                     press = Hydrostat.hydrostatic(numDeps, grav, tauRos, kappa, temp);
                     rho = State.massDensity(numDeps, temp, press, mmw, kappaScale);
                     //pressSun = Hydrostat.hydrostatic(numDeps, gravSun, tauRos, kappaSun, tempSun, logRadius);
                     //rhoSun = State.massDensity(numDeps, tempSun, pressSun, kappaScaleSun);
                     mode = 1;  //call kappas ** with ** knowledge of rho
                     Ne = State.NeFn(numDeps, temp, NeDfg2, kappaScale);
                     logNumsH3 = LevelPops.levelPops(lamJump3, logNH, Ne, ionizedHI, chiI1H, chiI2H, chiLH3, gw1H, gw2H, gwLH3,
                     numDeps, kappaScale, tauRos, temp, rho);
                     kappa = Kappas.kappas(mode, numDeps, rho, rhoSun, kapDfg2, kappaScale, logg, loggSun, teff, teffSun, radius, massX, massZ, tauRos, temp, tempSun);
                     //double kappaSun[][] = Kappas.kappas(numDeps, kappaScaleSun, teffSun, teffSun, loggSun, loggSun);
                     //kappaSun = Kappas.kappas(mode, numDeps, rho, rhoSun, kappaScaleSun, loggSun, loggSun, teffSun, teffSun, radiusSun, massX, massZSun, tauRos, tempSun, tempSun);
                     */
                    depths = DepthScale.depthScale(numDeps, tauRos, kappa, rho);

                    //Okay - Now all the emergent radiation stuff:
                    // Set up theta grid
                    //  cosTheta is a 2xnumThetas array:
                    // row 0 is used for Gaussian quadrature weights
                    // row 1 is used for cos(theta) values
                    // Gaussian quadrature:
                    // Number of angles, numThetas, will have to be determined after the fact
                    double cosTheta[][] = Thetas.thetas();
                    int numThetas = cosTheta[0].length;

                    boolean lineMode;

                    //
                    // ************
                    //
                    //  Spectrum synthesis section:
                    // Set up multi-Gray continuum info:
                    double isCool = 7300.0;  //Class A0

                    //Set up multi-gray opacity:
                    // lambda break-points and gray levels:
                    // No. multi-gray bins = num lambda breakpoints +1
                    double minLambda = 30.0;  //nm
                    double maxLambda = 1.0e6;  //nm
                    int maxNumBins = 11;
                    double[][] grayLevelsEpsilons = MulGrayTCorr.grayLevEps(maxNumBins, minLambda, maxLambda, teff, isCool);
                    //Find actual number of multi-gray bins:
                    int numBins = 0; //initialization
                    for (int i = 0; i < maxNumBins; i++) {
                        if (grayLevelsEpsilons[0][i] < maxLambda) {
                            numBins++;
                        }
                    }

                    //
                    //Line list:
                    int numLines = 14;
                    //int numLines = 1;
                    double[] listLam0 = new double[numLines];  // nm
                    double[] listMass = new double[numLines]; // amu
                    double[] listLogGammaCol = new double[numLines];
                    //abundance in logarithmic A12 sysytem
                    double[] listA12 = new double[numLines];
                    //Log of unitless oscillator strength, f 
                    double[] listLogf = new double[numLines];
                    //Ground state ionization E - Stage I (eV) 
                    double[] listChiI1 = new double[numLines];
                    //Ground state ionization E - Stage II (eV)
                    double[] listChiI2 = new double[numLines];
                    //Excitation E of lower E-level of b-b transition (eV)
                    double[] listChiL = new double[numLines];
                    //Unitless statisital weight, Ground state ionization E - Stage I
                    double[] listGw1 = new double[numLines];
                    //Unitless statisital weight, Ground state ionization E - Stage II
                    double[] listGw2 = new double[numLines];
                    //Unitless statisital weight, lower E-level of b-b transition                 
                    double[] listGwL = new double[numLines];
                    //double[] listGwU For now we'll just set GwU to 1.0
                    // Is stage II?
                    boolean[] listIonized = new boolean[numLines];

                    //Atomic Data sources:
                    //http://www.nist.gov/pml/data/asd.cfm
                    // Masses: http://www.chemicalelements.com/show/mass.html
                    //Solar abundances: http://arxiv.org/pdf/0909.0948.pdf
                    //   - Asplund, M., Grevesse, N., Sauval, A., Scott, P., 2009, arXiv:0909.0948v1
                    //
                    //
                    //
                    // Line list
                    //    *** CAUTION: THese should be in order of increasing wavelength (lam0) for labeling purposes in graphical output
                    //    
                    //    
                    //    
                    //    
                    //CaII K
                    //listName[0] = "Ca II K";
                    listLam0[0] = 393.366;
                    listA12[0] = 6.34;
                    listLogf[0] = -0.166;
                    listChiI1[0] = 6.113;
                    listChiI2[0] = 11.872;
                    //This is necessary for consistency with Stage II treatment of user-defined spectral line:
                    listChiL[0] = 0.01 + listChiI1[0];
                    listMass[0] = 40.078;
                    listLogGammaCol[0] = 1.0;
                    listGw1[0] = 1.0;
                    listGw2[0] = 2.0;
                    listGwL[0] = 2.0;
                    listIonized[0] = true;

                    //CaII H
                    //listName[1] = "Ca II H";
                    listLam0[1] = 396.847;
                    listA12[1] = 6.34;
                    listLogf[1] = -0.482;
                    listChiI1[1] = 6.113;
                    listChiI2[1] = 11.872;
                    //This is necessary for consistency with Stage II treatment of user-defined spectral line:
                    listChiL[1] = 0.01 + listChiI1[1];
                    listMass[1] = 40.078;
                    listLogGammaCol[1] = 1.0;
                    listGw1[1] = 1.0;
                    listGw2[1] = 2.0;
                    listGwL[1] = 2.0;
                    listIonized[1] = true;

                    //Fe I 4045
                    //listName[2] = "Fe I";
                    listLam0[2] = 404.581;
                    listA12[2] = 7.50; //??????
                    listLogf[2] = -0.674;
                    listChiI1[2] = 7.902;
                    listChiI2[2] = 16.199;
                    listChiL[2] = 1.485;
                    listMass[2] = 55.845;
                    listLogGammaCol[2] = 0.0;
                    listGw1[2] = 1.0;
                    listGw2[2] = 1.0;
                    listGwL[2] = 9.0;
                    listIonized[2] = false;

                    //Hdelta
                    //listName[3] = "H I &#948";
                    listLam0[3] = 410.174;
                    listA12[3] = 12.0;    //By definition - it's Hydrogen
                    listLogf[3] = -1.655;
                    listChiI1[3] = 13.6;
                    listChiI2[3] = 1.0e6;   //Set very high arbitrary value - there is no "H III"!
                    listChiL[3] = 10.2;
                    listMass[3] = 1.0;
                    listLogGammaCol[3] = 1.0;
                    listGw1[3] = 2.0; // 2n^2
                    listGw2[3] = 1.0;
                    listGwL[3] = 8.0; // 2n^2
                    listIonized[3] = false;

                    //CaI 4227
                    //listName[4] = "Ca I";
                    listLam0[4] = 422.673;
                    listA12[4] = 6.34;
                    listLogf[4] = 0.243;
                    listChiI1[4] = 6.113;
                    listChiI2[4] = 11.872;
                    listChiL[4] = 0.00;
                    listMass[4] = 40.078;
                    listLogGammaCol[4] = 1.0;
                    listGw1[4] = 1.0;
                    listGw2[4] = 1.0;
                    listGwL[4] = 1.0;
                    listIonized[4] = false;

                    //Fe I 4271
                    //listName[5] = "Fe I";
                    listLam0[5] = 427.176;
                    listA12[5] = 7.50; //??????
                    listLogf[5] = -1.118;
                    listChiI1[5] = 7.902;
                    listChiI2[5] = 16.199;
                    listChiL[5] = 1.485;
                    listMass[5] = 55.845;
                    listLogGammaCol[5] = 0.0;
                    listGw1[5] = 1.0;
                    listGw2[5] = 1.0;
                    listGwL[5] = 9.0;
                    listIonized[5] = false;

                    //Hgamma
                    //[6] = "H I &#947";
                    listLam0[6] = 434.047;
                    listA12[6] = 12.0;    //By definition - it's Hydrogen
                    listLogf[6] = -1.350;
                    listChiI1[6] = 13.6;
                    listChiI2[6] = 1.0e6;   //Set very high arbitrary value - there is no "H III"!
                    listChiL[6] = 10.2;
                    listMass[6] = 1.0;
                    listLogGammaCol[6] = 1.0;
                    listGw1[6] = 2.0; // 2n^2
                    listGw2[6] = 1.0;
                    listGwL[6] = 8.0; // 2n^2
                    listIonized[6] = false;

                    //He I 4387
                    //listName[7] = "He I";
                    listLam0[7] = 438.793;
                    listA12[7] = 10.93; //??????
                    listLogf[7] = -1.364;
                    listChiI1[7] = 24.587;
                    listChiI2[7] = 54.418;
                    listChiL[7] = 21.218;
                    listMass[7] = 4.003;
                    listLogGammaCol[7] = 0.0;
                    listGw1[7] = 1.0;
                    listGw2[7] = 1.0;
                    listGwL[7] = 3.0;
                    listIonized[7] = false;

                    //He I 4471
                    //listName[8] = "He I";
                    listLam0[8] = 447.147;
                    listA12[8] = 10.93; //??????
                    listLogf[8] = -0.986;
                    listChiI1[8] = 24.587;
                    listChiI2[8] = 54.418;
                    listChiL[8] = 20.964;
                    listMass[8] = 4.003;
                    listLogGammaCol[8] = 0.0;
                    listGw1[8] = 1.0;
                    listGw2[8] = 1.0;
                    listGwL[8] = 5.0;
                    listIonized[8] = false;

                    //Hbeta
                    //listName[9] = "H I &#946";
                    listLam0[9] = 486.128;
                    listA12[9] = 12.0;    //By definition - it's Hydrogen
                    listLogf[9] = -0.914;
                    listChiI1[9] = 13.6;
                    listChiI2[9] = 1.0e6;   //Set very high arbitrary value - there is no "H III"!
                    listChiL[9] = 10.2;
                    listMass[9] = 1.0;
                    listLogGammaCol[9] = 1.0;
                    listGw1[9] = 2.0;  // 2n^2
                    listGw2[9] = 1.0;
                    listGwL[9] = 8.0;  // 2n^2
                    listIonized[9] = false;

                    //MgIb1
                    //listName[10] = "Mg I <em>b</em><sub>1</sub>";
                    listLam0[10] = 518.360;  //nm
                    listA12[10] = 7.60;    // Grevesse & Sauval 98
                    listLogf[10] = -0.867;
                    listChiI1[10] = 7.646;
                    listChiI2[10] = 15.035;
                    listChiL[10] = 2.717;
                    listMass[10] = 24.305;
                    listLogGammaCol[10] = 1.0;
                    listGw1[10] = 1.0;
                    listGw2[10] = 1.0;
                    listGwL[10] = 5.0;
                    listIonized[10] = false;

                    //NaID2
                    //listName[11] = "Na I <em>D</em><sub>2</sub>";
                    listLam0[11] = 588.995;
                    listA12[11] = 6.24;    // Grevesse & Sauval 98
                    listLogf[11] = -0.193;
                    listChiI1[11] = 5.139;
                    listChiI2[11] = 47.286;
                    listChiL[11] = 0.0;
                    listMass[11] = 22.990;
                    listLogGammaCol[11] = 1.0;
                    listGw1[11] = 2.0;
                    listGw2[11] = 1.0;
                    listGwL[11] = 2.0;
                    listIonized[11] = false;

                    //NaID1
                    //listName[12] = "Na I <em>D</em><sub>1</sub>";
                    listLam0[12] = 589.592;  //nm
                    listA12[12] = 6.24;    // Grevesse & Sauval 98    
                    listLogf[12] = -0.495;
                    listChiI1[12] = 5.139;
                    listChiI2[12] = 47.286;
                    listChiL[12] = 0.0;
                    listMass[12] = 22.990;
                    listLogGammaCol[12] = 1.0;
                    listGw1[12] = 2.0;
                    listGw2[12] = 1.0;
                    listGwL[12] = 2.0;
                    listIonized[12] = false;

                    //Halpha
                    //listName[13] = "H I &#945";
                    listLam0[13] = 656.282;
                    listA12[13] = 12.0;    //By definition - it's Hydrogen
                    listLogf[13] = -0.193;
                    listChiI1[13] = 13.6;
                    listChiI2[13] = 1.0e6;   //Set very high arbitrary value - there is no "H III"!
                    listChiL[13] = 10.2;
                    listMass[13] = 1.0;
                    listLogGammaCol[13] = 1.0;
                    listGw1[13] = 2.0; // 2n^2
                    listGw2[13] = 1.0;
                    listGwL[13] = 8.0; // 2n^2
                    listIonized[13] = false;

                    //Notes
                    //
                    //CAUTION: This treatment expects numPoints (number of wavelengths, lambda) to be the same for *all* spectral lines!
                    //CONTINUUM lambda scale (nm)
                    double lambdaScale[] = LamGrid.lamgrid(numLams, lamSetup); //cm

                    int listNumCore = 5;  //per wing
                    int listNumWing = 5;  // half-core
                    //int numWing = 0;  //debug
                    int listNumPoints = 2 * (listNumCore + listNumWing) - 1; // + 1;  //Extra wavelength point at end for monochromatic continuum tau scale
                    int numMaster = numLams + (numLines * listNumPoints); //total size (number of wavelengths) of master lambda & total kappa arrays 
                    int numNow = numLams;  //initialize dynamic coutner of how many array elements are in use
                    double[] masterLams = new double[numMaster];
                    double[][] logMasterKaps = new double[numMaster][numDeps];
                    //System.out.println("numLams " + numLams + " numLines " + numLines + " listNumPoints " + listNumPoints);
                    //System.out.println("numMaster " + numMaster + " numNow " + numNow);
                    //seed masterLams and logMasterKaps with continuum SED lambdas and kapaps:
                    //This just initializes the first numLams of the numMaster elements
                    //Also - put in multi-Gray opacities here:
                    //Find which gray level bin the spectrum synthesis region starts in - assume that the first gray-level bin
                    // is always at a shorter wavelength than the start of the synthesis region:
                    int whichBin = 0;  //initialization
                    for (int iB = 0; iB < numBins; iB++) {
                        if (grayLevelsEpsilons[0][iB] >= lambdaScale[0]) {
                            System.out.println("grayLevelsEpsilons[0][iB] " + grayLevelsEpsilons[0][iB] + " lambdaScale[0] " + lambdaScale[0]);
                            whichBin = iB;  //found it!
                            break;
                        }
                    }
                    //System.out.println("starting whichBin " + whichBin);

                    //First wavelength definitely falls in first found bin:
                    masterLams[0] = lambdaScale[0];
                    for (int iD = 0; iD < numDeps; iD++) {
                        logMasterKaps[0][iD] = kappa[1][iD] + Math.log(grayLevelsEpsilons[1][0]);
                    }
                    for (int iL = 1; iL < numLams; iL++) {
                        masterLams[iL] = lambdaScale[iL];
                        //System.out.println("iL " + iL + " lambdaScale[iL] " + lambdaScale[iL] + " whichBin+1 " + (whichBin + 1) + " grayLevelsEpsilons[0][whichBin + 1] " + grayLevelsEpsilons[0][whichBin + 1]);
                        if ((lambdaScale[iL] >= grayLevelsEpsilons[0][whichBin + 1])
                                && (lambdaScale[iL - 1] < grayLevelsEpsilons[0][whichBin + 1])
                                && (whichBin < numBins - 1)) {
                            whichBin++;
                            //System.out.println("whichBin " + whichBin);
                        }
                        for (int iD = 0; iD < numDeps; iD++) {
                            logMasterKaps[iL][iD] = kappa[1][iD] + Math.log(grayLevelsEpsilons[1][whichBin]);
                        }
                    }
                    //initialize the rest with dummy values
                    for (int iL = numLams; iL < numMaster; iL++) {
                        masterLams[iL] = lambdaScale[numLams - 1];
                        for (int iD = 0; iD < numDeps; iD++) {
                            logMasterKaps[iL][iD] = kappa[1][iD];
                        }
                    }

                    //Stuff for the the Teff recovery test:
                    double lambda1, lambda2, fluxSurfBol, logFluxSurfBol;
                    fluxSurfBol = 0;

                    for (int iLine = 0; iLine < numLines; iLine++) {
                        //System.out.println("iLine " + iLine + " numNow " + numNow);
                        double listLogN = (listA12[iLine] - 12.0) + logNH;
                        listLam0[iLine] = listLam0[iLine] * 1.0e-7;  // nm to cm
                        double[][] listLinePoints = LineGrid.lineGrid(listLam0[iLine], listMass[iLine], xiT, numDeps, teff, listNumCore, listNumWing);
                        // // Real Voigt fn profile (voigt2()):        
                        //double[][] listLineProf = LineProf.voigt2(listLinePoints, listLam0[iLine], listLogGammaCol[iLine],
                        //        numDeps, teff, tauRos, temp, press, tempSun, pressSun);
                        // Gaussian + Lorentzian approximation to profile (voigt()):
                        double[][] listLineProf = LineProf.voigt(listLinePoints, listLam0[iLine], listLogGammaCol[iLine],
                                numDeps, teff, tauRos, temp, press, tempSun, pressSun);
                        //double[][] listLogNums = LevelPops.levelPops(listLam0[iLine], listLogN, Ne, ionized, listChiI[iLine], listChiL[iLine], 
                        //        numDeps, kappaScale, tauRos, temp, rho);
                        double[][] listLogNums = LevelPops.levelPops(listLam0[iLine], listLogN, Ne, listIonized[iLine], listChiI1[iLine],
                                listChiI2[iLine], listChiL[iLine], listGw1[iLine], listGw2[iLine], listGwL[iLine],
                                numDeps, kappaScale, tauRos, temp, rho);
                        double[][] listLogKappaL = LineKappa.lineKap(listLam0[iLine], listLogNums, listLogf[iLine], listLinePoints, listLineProf,
                                numDeps, kappaScale, tauRos, temp, rhoSun);
                        //int listNumPoints = listLinePoints[0].length; // + 1;  //Extra wavelength point at end for monochromatic continuum tau scale
                        //double logTauL[][] = LineTau2.tauLambda(numDeps, listNumPoints, logKappaL,
                        //        kappa, tauRos, rho, logg);
                        double[] listLineLambdas = new double[listNumPoints];
                        for (int il = 0; il < listNumPoints; il++) {
                            // // lineProf[iLine][*] is DeltaLambda from line centre in cm
                            // if (il == listNumPoints - 1) {
                            //    listLineLambdas[il] = listLam0[iLine]; // Extra row for line centre continuum taus scale
                            // } else {
                            //lineLambdas[il] = (1.0E7 * linePoints[iLine][il]) + lam0; //convert to nm
                            listLineLambdas[il] = listLinePoints[0][il] + listLam0[iLine];
                            // }
                        }

                        double[] masterLamsOut = SpecSyn.masterLambda(numLams, numMaster, numNow, masterLams, listNumPoints, listLineLambdas);
                        double[][] logMasterKapsOut = SpecSyn.masterKappa(numDeps, numLams, numMaster, numNow, masterLams, masterLamsOut, logMasterKaps, listNumPoints, listLineLambdas, listLogKappaL);
                        numNow = numNow + listNumPoints;

                        //update masterLams and logMasterKaps:
                        for (int iL = 0; iL < numNow; iL++) {
                            masterLams[iL] = masterLamsOut[iL];
                            for (int iD = 0; iD < numDeps; iD++) {
                                //Still need to put in multi-Gray levels here:
                                logMasterKaps[iL][iD] = logMasterKapsOut[iL][iD];
                                //if (iD == 36) {
                                //    System.out.println("iL " + iL + " masterLams[iL] " + masterLams[iL] + " logMasterKaps[iL][iD] " + logMasterKaps[iL][iD]);
                                //}
                            }
                        }
                    } //numLines loop

                    //int numMaster = masterLams.length;
                    double logTauMaster[][] = LineTau2.tauLambda(numDeps, numMaster, logMasterKaps,
                            kappa, tauRos);

                    //Evaluate formal solution of rad trans eq at each lambda throughout line profile
                    // Initial set to put lambda and tau arrays into form that formalsoln expects
                    //double[] masterLambdas = new double[numMaster];
                    double[][] masterIntens = new double[numMaster][numThetas];
                    double[] masterIntensLam = new double[numThetas];

                    double[][] masterFlux = new double[2][numMaster];
                    double[] masterFluxLam = new double[2];

                    double[][] thisTau = new double[2][numDeps];

                    //                   double[][] lineJay = new double[2][listNumPoints];
//                    double[] lineJayLam = new double[2];
                    //double[][] thisTau = new double[2][numDeps];
                    lineMode = false;  //no scattering for overall SED

                    for (int il = 0; il < numMaster; il++) {

//                        // lineProf[0][*] is DeltaLambda from line centre in cm
//                        if (il == listNumPoints - 1) {
//                            lineLambdas[il] = lam0; // Extra row for line centre continuum taus scale
//                        } else {
                        //lineLambdas[il] = (1.0E7 * linePoints[0][il]) + lam0; //convert to nm
                        //masterLambdas[il] = masterLams[il];
//                        }
                        for (int id = 0; id < numDeps; id++) {
                            thisTau[1][id] = logTauMaster[il][id];
                            thisTau[0][id] = Math.exp(logTauMaster[il][id]);
                        } // id loop

                        masterIntensLam = FormalSoln.formalSoln(numDeps,
                                cosTheta, masterLams[il], thisTau, temp, lineMode);

                        masterFluxLam = Flux.flux(masterIntensLam, cosTheta);

                        for (int it = 0; it < numThetas; it++) {
                            masterIntens[il][it] = masterIntensLam[it];
                            //System.out.println(" il " + il + " it " + it + " logIntens " + logE*Math.log(lineIntensLam[it]) );
                        } //it loop - thetas

                        masterFlux[0][il] = masterFluxLam[0];
                        masterFlux[1][il] = masterFluxLam[1];

                        //System.out.println("il " + il + " masterLams[il] " + masterLams[il] + " masterFlux[1][il] " + logE * masterFlux[1][il]);
                        //// Teff test - Also needed for convection module!:
                        if (il > 1) {
                            lambda2 = masterLams[il]; // * 1.0E-7;  // convert nm to cm
                            lambda1 = masterLams[il - 1]; // * 1.0E-7;  // convert nm to cm
                            fluxSurfBol = fluxSurfBol
                                    + masterFluxLam[0] * (lambda2 - lambda1);
                        }
                    } //il loop

                    logFluxSurfBol = Math.log(fluxSurfBol);
                    double logTeffFlux = (logFluxSurfBol - Useful.logSigma()) / 4.0;
                    double teffFlux = Math.exp(logTeffFlux);
                    String pattern = "0000.00";
                    ////String pattern = "#####.##";
                    DecimalFormat myFormatter = new DecimalFormat(pattern);

                    ////Teff test
                    System.out.println("FLUX: Recovered Teff = " + myFormatter.format(teffFlux));
                    //Compute JOhnson-Cousins photometric color indices:

                    double colors[] = Photometry.UBVRI(masterLams, masterFlux, numDeps, tauRos, temp);
                    // test double icolors[][] = Photometry.iColors(lambdaScale, intens, numDeps, numThetas, numLams, tauRos, temp); //testing here - really for disk image rendering
//
                    //
                    //
                    // Line profile section:
                    //
                    //
                    // Set up line grid of lambda points sampling entire profile:
                    int numCore = 5;  //per wing
                    int numWing = 15;  // half-core
                    //int numWing = 0;  //debug
                    int numPoints = 2 * (numCore + numWing) - 1; // + 1;  //Extra wavelength point at end for monochromatic continuum tau scale
                    //linePoints: Row 0 in cm (will need to be in nm for Plack.planck), Row 1 in Doppler widths
                    double linePoints[][] = LineGrid.lineGrid(lam0, mass, xiT, numDeps, teff, numCore, numWing);

                    //Compute area-normalized depth-independent line profile "phi_lambda(lambda)"
                    // In class LineProf, method voigt() is the Gaussian + Lorentzian approximation
                    //                  , method voigt2() is the actual Voigt profile
                    double lineProf[][] = LineProf.voigt(linePoints, lam0, logGammaCol,
                            numDeps, teff, tauRos, temp, press, tempSun, pressSun);
                    //double lineProf[][] = LineProf.voigt2(linePoints, lam0, logGammaCol,
                    //        numDeps, teff, tauRos, temp, press, tempSun, pressSun);

                    // Level population now computed in LevelPops.levelPops()
                    //double logNums[][] = LevelPops.levelPops(lam0, logN, Ne, ionized, chiI1, chiL, 
                    //        numDeps, kappaScale, tauRos, temp, rho);
                    double logNums[][] = LevelPops.levelPops(lam0, logN, Ne, ionized, chiI1, chiI2, chiL, gw1, gw2, gwL,
                            numDeps, kappaScale, tauRos, temp, rho);

                    //Compute depth-dependent logarithmic monochromatic extinction co-efficient, kappa_lambda(lambda, tauRos):
                    //Handing in rhoSun instead of rho here is a *weird* fake to get line broadening to scale with logg 
                    //approximately okay for saturated lines:   There's something wrong!  
                    double logKappaL[][] = LineKappa.lineKap(lam0, logNums, logf, linePoints, lineProf,
                            numDeps, kappaScale, tauRos, temp, rhoSun);
                    double[][] logTotKappa = LineKappa.lineTotalKap(linePoints, logKappaL, numDeps, teff, lam0, kappaScale, kappa);

                    //For spectrum synthesis:
                    //public static double[][] tauLambda(int numDeps, double[][] lineProf, double[][] logKappaL, 
                    //                             double[][] kappa, double[][] tauRos, double[][] rho, double[] depths) {
                    //Compute monochromatic optical depth scale, Tau_lambda throughout line profile
                    //CAUTION: Returns numPoints+1 x numDeps array: the numPoints+1st row holds the line centre continuum tau scale
                    // Method 1: double logTauL[][] = LineTau.tauLambda(numDeps, lineProf, logKappaL,
                    // Method 1:        kappa, tauRos, rho, depths);
                    // Method 2:
                    double logTauL[][] = LineTau2.tauLambda(numDeps, numPoints, logTotKappa,
                            kappa, tauRos);

                    //Evaluate formal solution of rad trans eq at each lambda throughout line profile
                    // Initial set to put lambda and tau arrays into form that formalsoln expects
                    double[] lineLambdas = new double[numPoints];

                    double[][] lineIntens = new double[numPoints][numThetas];
                    double[] lineIntensLam = new double[numThetas];

                    double[][] lineFlux = new double[2][numPoints];
                    double[] lineFluxLam = new double[2];

                    lineMode = true;  //scattering not working correctly?

                    for (int il = 0; il < numPoints; il++) {

                        // lineProf[0][*] is DeltaLambda from line centre in cm
                        // if (il == numPoints - 1) {
                        //     lineLambdas[il] = lam0; // Extra row for line centre continuum taus scale
                        //} else {
                        //lineLambdas[il] = (1.0E7 * linePoints[0][il]) + lam0; //convert to nm
                        lineLambdas[il] = linePoints[0][il] + lam0;
                        //}

                        for (int id = 0; id < numDeps; id++) {
                            thisTau[1][id] = logTauL[il][id];
                            thisTau[0][id] = Math.exp(logTauL[il][id]);
                            //if (il == 9) {
                            //System.out.println("lineLambdas[il]" + lineLambdas[il]);
                            //System.out.println("id " + id + " thisTau[1] " + logE*thisTau[1][id] + " tauRos[1] " + logE*tauRos[1][id]);
                            //}
                        } // id loop

                        lineIntensLam = FormalSoln.formalSoln(numDeps,
                                cosTheta, lineLambdas[il], thisTau, temp, lineMode);

                        lineFluxLam = Flux.flux(lineIntensLam, cosTheta);

                        for (int it = 0; it < numThetas; it++) {
                            lineIntens[il][it] = lineIntensLam[it];
                            //System.out.println(" il " + il + " it " + it + " logIntens " + logE*Math.log(lineIntensLam[it]) );
                        } //it loop - thetas

                        lineFlux[0][il] = lineFluxLam[0];
                        lineFlux[1][il] = lineFluxLam[1];

                    } //il loop

                    //Continuum flux at line centre for Eq width calculation:
                    //int ilLam0 = LamPoint.lamPoint(numLams, lambdaScale, lam0);
                    // Solve formal sol of rad trans eq for outgoing surfaace I(0, theta)
                    //double intens[][] = FormalSoln.formalsoln(numLams, numDeps,
                    //        cosTheta, lambdaScale, tauRos, temp);
                    double[] intensCont = new double[numThetas];
                    double[] fluxCont = new double[2];

                    //double[][] intens = new double[numLams][numThetas];
                    //double[][] flux = new double[2][numLams];
                    //  double[][] intens = new double[3][numThetas];
                    //double[][] flux = new double[2][3];
                    //double lambda1, lambda2, fluxSurfBol, logFluxSurfBol;
                    //fluxSurfBol = 0;
                    lineMode = false;

                    //for (int il = ilLam0-1; il <= ilLam0+1; il++) {
                    //System.out.println("ilLam0 " + ilLam0 + " lambdaScale[ilLam0] " + lambdaScale[ilLam0] + " lam0 " + lam0);
                    //intensCont = FormalSoln.formalSoln(numDeps,
                    //        cosTheta, lambdaScale[ilLam0], tauRos, temp, lineMode);
                    //intensCont = FormalSoln.formalSoln(numDeps,
                    //        cosTheta, lam0, tauRos, temp, lineMode);
//
                    //fluxCont = Flux.flux(intensCont, cosTheta);
                    //System.out.println("fluxCont[1] " + logE*fluxCont[1]);
                    //for (int it = 0; it < numThetas; it++) {
                    //intens[il][it] = intensLam[it];
                    //} //it loop - thetas
                    //flux[0][il] = fluxLam[0];
                    //flux[1][il] = fluxLam[1];
                    //System.out.println("il " + il + " lambdaScale[il] " + lambdaScale[il] + " flux[1][il] " + logE*flux[1][il]);
                    //// Teff test - Also needed for convection module!:
                    //if (il > 1) {
                    //lambda2 = lambdaScale[il]; // * 1.0E-7;  // convert nm to cm
                    //lambda1 = lambdaScale[il - 1]; // * 1.0E-7;  // convert nm to cm
                    //fluxSurfBol = fluxSurfBol
                    // + fluxLam[0] * (lambda2 - lambda1);
                    //}
                    //} //il loop - lambdas
                    //Compute flux at each lambda throughout line profile
                    //lineFlux = Flux.flux(lineIntens, cosTheta, numPoints, lineLambdas);
                    //Get equivalent width, W_lambda, in pm - picometers:
                    double Wlambda = EqWidth.eqWidth(lineFlux, linePoints, lam0);

                    String patternCol = "0.00";
                    //String pattern = "#####.##";
                    DecimalFormat colFormatter = new DecimalFormat(patternCol);
                    String patternWl = "0.00";
                    //String pattern = "#####.##";
                    DecimalFormat WlFormatter = new DecimalFormat(patternWl);

                    actiontarget.setText("Photometric color indices: "
                            + "U-B: " + colFormatter.format(colors[0])
                            + " B-V: " + colFormatter.format(colors[1])
                            + " V-R: " + colFormatter.format(colors[2])
                            + " V-I: " + colFormatter.format(colors[3])
                            + " R-I: " + colFormatter.format(colors[4]) + "\r\n"
                            + "Spectral line equivalent width: " + WlFormatter.format(Wlambda) + " pm");

                    // No! //grid.getChildren().add(r);
                    // Graphical output section:
                    // Plot 1: T_Kin(depth):
                    LineChart<Number, Number> lineChartT2 = LineCharts.t2Plot(numDeps, depths, temp, tauRos);
                    grid.add(lineChartT2, 0, 9);

                    // Plot 2: T_Kin(log(Tau_Ros)):
                    //grid.add(r, 0, 9);
                    LineChart<Number, Number> lineChartT = LineCharts.tPlot(numDeps, tauRos, temp);
                    grid.add(lineChartT, 2, 9);

                    // Plot 3: log(P(log(Tau_Ros))):
                    LineChart<Number, Number> lineChartP = LineCharts.pPlot(numDeps, tauRos, press);
                    grid.add(lineChartP, 4, 9);

                    // Plot 4: log(I_lambda(lambda, cos(theta))):
                    LineChart<Number, Number> lineChartLimb = LineCharts.limbPlot(numLams, masterLams, cosTheta, masterIntens, masterFlux, numPoints, lineLambdas, lineIntens, lineFlux);
                    grid.add(lineChartLimb, 0, 10);

                    // Plot 5: Continuum log(F_lambda(log(lambda), log(I_lambda(log(lambda))):
                    //LineChart<Number, Number> lineChartSpec = LineCharts.specPlot(numLams, lambdaScale, cosTheta, intens, flux);
                    //grid.add(lineChartSpec, 2, 10);
                    LineChart<Number, Number> lineChartSpec = LineCharts.specPlot(numMaster, masterLams, cosTheta, masterIntens, masterFlux);
                    grid.add(lineChartSpec, 2, 10);

                    // Plot 6: Line log(F_lambda(log(lambda), log(I_lambda(log(lambda))):
                    LineChart<Number, Number> lineChartLine = LineCharts.linePlot(numPoints, lineLambdas, cosTheta, lineIntens, lineFlux);
                    //LineChart<Number, Number> lineChartLine = LineCharts.linePlot(numPoints, lineLambdas, cosTheta, lineProf, lam0);
                    grid.add(lineChartLine, 4, 10);

                    //Debug versions of the function with parameters for plotting up quqntities used in line profile
                    // calculation:
                    //LineChart<Number, Number> lineChartLine = LineCharts.linePlot(numPoints, lineLambdas, cosTheta, logKappaL, lam0, kappa);
                    //grid.add(lineChartLine, 4, 10);
                    //LineChart<Number, Number> lineChartLine = LineCharts.linePlot(numPoints, lineLambdas, cosTheta, logTauL, lam0);
                    //grid.add(lineChartLine, 4, 10);
                    // Final scene / stage stuff:
                    //Scene scene = new Scene(grid, 1600, 900);
                    //scene.getStylesheets().add("../../GrayCascadeStyleSheet.css");
                    //
                    //primaryStage.setScene(scene);
                    //
                    //primaryStage.show();
                } else {
                    actiontarget.setText("All fields must have values");
                }

            }
        }
        );

        // Final scene / stage stuff:
        Scene scene = new Scene(grid, 1600, 900);

        //Stylesheet.css must go in /src/ directory
        scene.getStylesheets()
                .add("GrayCascadeStyleSheet.css");

        primaryStage.setScene(scene);

        primaryStage.show();

    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        launch(args);
    }

}
