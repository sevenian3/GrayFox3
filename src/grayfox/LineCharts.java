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

import java.text.DecimalFormat;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
//import javafx.scene.chart.XYChart.Data;
//import javafx.scene.chart.XYChart.getXValue;

public class LineCharts {

    /**
     * Plot 1: T_Kin(log(depth)):
     *
     * @param numDeps
     * @param depths
     * @param temp
     * @param tauRos
     * @return
     */
    public static LineChart<Number, Number> t2Plot(int numDeps, double[] depths, double[][] temp, double[][] tauRos) {

        double logE = Math.log10(Math.E);

        double minX = 1.0E-5 * depths[0];
        double maxX = 1.0E-5 * depths[numDeps - 1];
        double deltaX = 100.0;

        int minYint = (int) ((0.9 * temp[0][0]) / 1000.0);
        double minY = (double) (minYint) * 1000.0;
        int maxYint = (int) ((1.1 * temp[0][numDeps - 1]) / 1000.0);
        double maxY = (double) (maxYint) * 1000.0;
        double deltaY = 1000.0;

        //final NumberAxis xAxis = new NumberAxis(minX, maxX,deltaX);
        //final NumberAxis yAxis = new NumberAxis(minY, maxY,deltaY);
        final NumberAxis xAxis;
        xAxis = new NumberAxis();
        final NumberAxis yAxis;
        yAxis = new NumberAxis();

        String xTitle;
        xTitle = "depth (km)";
        xAxis.setLabel(xTitle);

        String yTitle;
        yTitle = "T_Kin (K)";
        yAxis.setLabel(yTitle);

        //String teffLbl = "4500";
        final LineChart<Number, Number> lineChartT2;
        lineChartT2 = new LineChart<>(xAxis, yAxis);

        lineChartT2.setTitle("Kinetic temperature vs Depth"); //\r\n" + teffLbl + "K");
        lineChartT2.setId("tempDpth");

        //javafx.scene.chart.ValueAxis.lowerBound = 0.0;
        //lineChartT.lowerBound = 0.0;
        XYChart.Series seriesT;
        seriesT = new XYChart.Series();
        seriesT.setName("T_Kin");

        double depths2[] = new double[numDeps];
        double conv = 1.0E-5;
        for (int id = 0; id < numDeps; id++) {

            depths2[id] = conv * depths[id];
            XYChart.Data<?, ?> xyPair;
            xyPair = new XYChart.Data(depths2[id], temp[0][id]);
            boolean add;
            //add = seriesT.getData().add(new XYChart.Data(conv * depths[id], temp[0][id]));
            add = seriesT.getData().add(xyPair);
        }

        int tTau1;
        tTau1 = TauPoint.tauPoint(numDeps, tauRos, 1.0);

        // vertical marker for tau=1
        XYChart.Series seriesdTau1;
        seriesdTau1 = new XYChart.Series();
        seriesdTau1.setName("Tau=1");
        boolean add;
        add = seriesdTau1.getData().add(new XYChart.Data(1.0E-5 * depths[tTau1], minY));
        boolean add1;
        add1 = seriesdTau1.getData().add(new XYChart.Data(1.0E-5 * depths[tTau1], maxY));

        // horizontal marker for tau=1
        XYChart.Series seriestTau1;
        seriestTau1 = new XYChart.Series();
        seriestTau1.setName(" ");
        boolean add2;
        add2 = seriestTau1.getData().add(new XYChart.Data(minX, temp[0][tTau1]));
        boolean add3;
        add3 = seriestTau1.getData().add(new XYChart.Data(maxX, temp[0][tTau1]));

        //lineChartT2.getData().add(seriesT);
        boolean addAll; //, seriesTau23);
        addAll = lineChartT2.getData().addAll(seriesT, seriesdTau1, seriestTau1);

        lineChartT2.setCreateSymbols(false);

        return lineChartT2;

    }

    /**
     * Plot 2: T_Kin(log(Tau_Ros)):
     *
     * @param numDeps
     * @param tauRos
     * @param temp
     * @return
     */
    public static LineChart<Number, Number> tPlot(int numDeps, double[][] tauRos, double[][] temp) {

        double logE = Math.log10(Math.E);

        double minX = logE * tauRos[1][0];
        double maxX = logE * tauRos[1][numDeps - 1];
        double deltaX = 1.0;

        int minYint = (int) ((0.9 * temp[0][0]) / 1000.0);
        double minY = (double) (minYint) * 1000.0;
        int maxYint = (int) ((1.1 * temp[0][numDeps - 1]) / 1000.0);
        double maxY = (double) (maxYint) * 1000.0;
        double deltaY = 1000.0;

        final NumberAxis xAxis = new NumberAxis(minX, maxX, deltaX);
        final NumberAxis yAxis = new NumberAxis(minY, maxY, deltaY);

        String xTitle = "Log Tau_Ros";
        xAxis.setLabel(xTitle);

        String yTitle = "T_Kin (K)";
        yAxis.setLabel(yTitle);

        //String teffLbl = "4500";
        final LineChart<Number, Number> lineChartT
                = new LineChart<Number, Number>(xAxis, yAxis);

        lineChartT.setTitle("Kinetic temperature vs log(Tau_Ros)"); //\r\n" + teffLbl + "K");
        lineChartT.setId("tempTau");

        XYChart.Series seriesT = new XYChart.Series();
        seriesT.setName("T_Kin");

        for (int id = 0; id < numDeps; id++) {

            seriesT.getData().add(new XYChart.Data(logE * tauRos[1][id], temp[0][id]));

        }

        // vertical marker for tau=1
        XYChart.Series seriesTau1 = new XYChart.Series();
        seriesTau1.setName("Tau=1");
        seriesTau1.getData().add(new XYChart.Data(0.0, minY));
        seriesTau1.getData().add(new XYChart.Data(0.0, maxY));

        // horizontal marker for tau=1
        int tTau1 = TauPoint.tauPoint(numDeps, tauRos, 1.0);
        XYChart.Series seriestTau1 = new XYChart.Series();
        seriestTau1.setName(" ");
        //System.out.println("logTauRos range: " + logE*tauRos[1][0] + " " + logE*tauRos[1][numDeps-1]);
        seriestTau1.getData().add(new XYChart.Data(minX, temp[0][tTau1]));
        seriestTau1.getData().add(new XYChart.Data(maxX, temp[0][tTau1]));

        //seriesTau23.getData().add(new XYChart.Data(Math.log10(2.0/3.0), 3000.0));
        //seriesTau23.getData().add(new XYChart.Data(Math.log10(2.0/3.0), 9000.0));
        //lineChartT.getData().add(seriesT);
        lineChartT.getData().addAll(seriesT, seriesTau1, seriestTau1); //, seriesTau23);

        lineChartT.setCreateSymbols(false);

        return lineChartT;

    }

    /**
     * Plot 3: log(P(log(Tau_Ros))):
     *
     * @param numDeps
     * @param tauRos
     * @param press
     * @return
     */
    public static LineChart<Number, Number> pPlot(int numDeps, double[][] tauRos, double[][] press) {

        double logE = Math.log10(Math.E);

        double minX = logE * tauRos[1][0];
        double maxX = logE * tauRos[1][numDeps - 1];
        double deltaX = 1.0;

        // From Hydrostat.hydrostat:
        //press is a 4 x numDeps array:
        // rows 0 & 1 are linear and log *gas* pressure, respectively
        // rows 2 & 3 are linear and log *radiation* pressure
        // Don't use upper boundary condition as lower y-limit - use a couple of points below surface:
        //int minYint = (int) ( Math.log10(0.9 * (press[0][2] + press[2][2])) / 10.0 );
        //double minY = (double) (minYint) * 10.0;
        //int maxYint = (int) ( Math.log10(1.1 * (press[0][numDeps-1] + press[2][numDeps-1])) / 10.0 );
        //double maxY = (double) (maxYint) * 10.0;
        int minYint = (int) (Math.log10(0.01 * (press[0][1] + press[2][1])));
        double minY = (double) (minYint);
        int maxYint = (int) (Math.log10(10.0 * (press[0][numDeps - 1] + press[2][numDeps - 1])));
        double maxY = (double) (maxYint);
        double deltaY = 1.0;

        final NumberAxis xAxis = new NumberAxis(minX, maxX, deltaX);
        final NumberAxis yAxis = new NumberAxis(minY, maxY, deltaY);

        xAxis.setLabel("Log Tau_Ros");

        String yTitle = "log P (dynes/cm/cm)";
        yAxis.setLabel(yTitle);

        //String teffLbl = "4500";
        final LineChart<Number, Number> lineChartP
                = new LineChart<Number, Number>(xAxis, yAxis);

        lineChartP.setTitle("Log Pressure vs log(Tau_Ros)"); // \r\n" + teffLbl + "K");
        lineChartP.setId("pressTau");

        XYChart.Series seriesP = new XYChart.Series();
        seriesP.setName("log P_Tot");

        XYChart.Series seriesPg = new XYChart.Series();
        seriesPg.setName("log P_Gas");

        XYChart.Series seriesPr = new XYChart.Series();
        seriesPr.setName("log P_Rad");

        // From Hydrostat.hydrostat:
        //press is a 4 x numDeps array:
        // rows 0 & 1 are linear and log *gas* pressure, respectively
        // rows 2 & 3 are linear and log *radiation* pressure
        for (int id = 0; id < numDeps; id++) {

            seriesP.getData().add(new XYChart.Data(logE * tauRos[1][id], Math.log10(press[0][id] + press[2][id]))); // Total pressure
            seriesPg.getData().add(new XYChart.Data(logE * tauRos[1][id], logE * press[1][id])); // gas pressure
            seriesPr.getData().add(new XYChart.Data(logE * tauRos[1][id], logE * press[3][id])); // radiation pressure

        }

        // Vertical marker for tau=1
        XYChart.Series seriesTau1 = new XYChart.Series();
        seriesTau1.setName("Tau=1");
        seriesTau1.getData().add(new XYChart.Data(0.0, logE * minY));
        seriesTau1.getData().add(new XYChart.Data(0.0, logE * maxY));

        // horizontal marker for tau=1
        int pTau1 = TauPoint.tauPoint(numDeps, tauRos, 1.0);
        XYChart.Series seriespTau1 = new XYChart.Series();
        seriespTau1.setName(" ");
        seriespTau1.getData().add(new XYChart.Data(minX, logE * press[1][pTau1]));
        seriespTau1.getData().add(new XYChart.Data(maxX, logE * press[1][pTau1]));

        //lineChartP.getData().add(seriesP);
        lineChartP.getData().addAll(seriesPr, seriesPg, seriesP, seriesTau1, seriespTau1);

        lineChartP.setCreateSymbols(false);

        return lineChartP;

    }

    /**
     * Plot 4: Limb darkening: log(I_lambda(lambda, cos(theta))):
     *
     * @param numLams
     * @param lambdaScale
     * @param cosTheta
     * @param intens
     * @param flux
     * @param numPoints
     * @param lineLambdas
     * @param lineIntens
     * @param lineFlux
     * @return
     */
    public static LineChart<Number, Number> limbPlot(int numLams,
            double[] lambdaScale, double[][] cosTheta, double[][] intens, double[][] flux,
            int numPoints, double[] lineLambdas, double[][] lineIntens, double[][] lineFlux) {
        //limbPlot(numLams, 
        //lambdaScale, cosTheta, intens, 
        //numPoints, lineIntens)    

        double logE = Math.log10(Math.E);

        final NumberAxis xAxis = new NumberAxis();
        final NumberAxis yAxis = new NumberAxis();

        xAxis.setLabel("theta (deg)");

        String yTitle = "I_lambda(theta)/I_lambda(0)";
        yAxis.setLabel(yTitle);

        int numMus = cosTheta[0].length;

        final LineChart<Number, Number> lineChartLimb
                = new LineChart<Number, Number>(xAxis, yAxis);

        lineChartLimb.setTitle("Limb darkening profiles"); //\r\n" + teffLbl + "K");
        lineChartLimb.setId("limb");

        String pattern = "0.0";
        //String pattern = "###.####";
        DecimalFormat myFormatter = new DecimalFormat(pattern);

        String lamStr = myFormatter.format(1.0e7 * lambdaScale[0]);
        XYChart.Series seriesL0 = new XYChart.Series();
        seriesL0.setName(lamStr + " nm");

        // Intermediate line series should be near lambda_max:
        int[] iLamMinMax = MinMax2.minMax(flux);
        int iLamMax = iLamMinMax[1];
        lamStr = myFormatter.format(1.0e7 * lambdaScale[iLamMax]);
        XYChart.Series seriesLn2 = new XYChart.Series();
        seriesLn2.setName("Max flux: " + lamStr + " nm");

        lamStr = myFormatter.format(1.0e7 * lambdaScale[numLams - 1]);
        XYChart.Series seriesLn = new XYChart.Series();
        seriesLn.setName(lamStr + " nm");

        // Line center limb darkening:
        int[] keyLambds = HalfPower.halfPower(numPoints, lineFlux);
        lamStr = myFormatter.format(lineLambdas[keyLambds[0]]);
        //System.out.println("index lam_0, lam_0: " + keyLambds[0] + " " + lineLambdas[keyLambds[0]]);
        XYChart.Series seriesLine = new XYChart.Series();
        seriesLine.setName("Line: " + lamStr + " nm");

        double theta;
        for (int it = 0; it < numMus; it++) {

            theta = (180.0 / Math.PI) * Math.acos(cosTheta[1][it]);
            //System.out.println("it, cos(theta): " + it + " " + cosTheta[0][it] + " theta: " + theta);
            seriesL0.getData().add(new XYChart.Data(theta, intens[0][it] / intens[0][0]));
            seriesLn2.getData().add(new XYChart.Data(theta, intens[iLamMax][it] / intens[iLamMax][0]));
            seriesLn.getData().add(new XYChart.Data(theta, intens[numLams - 1][it] / intens[numLams - 1][0]));

            // Line centre limb darkening:
            seriesLine.getData().add(new XYChart.Data(theta, lineIntens[keyLambds[0]][it] / lineIntens[keyLambds[0]][0]));

        }
        //lineChartLimb.getData().add(seriesL0);
        lineChartLimb.getData().addAll(seriesL0, seriesLn2, seriesLn, seriesLine);

        return lineChartLimb;

    }

    /**
     * Plot 5: SED: log(F_lambda(log(lambda)), I_lambda(log(lambda))):
     *
     * @param numLams
     * @param lambdaScale
     * @param cosTheta
     * @param intens
     * @param flux
     * @param lamUBVRI
     * @return
     */
    public static LineChart<Number, Number> specPlot(int numLams,
            double[] lambdaScale, double[][] cosTheta, double[][] intens,
            double[][] flux) {

        double logE = Math.log10(Math.E);
        int numMus = cosTheta[0].length;

        final NumberAxis xAxis = new NumberAxis();
        final NumberAxis yAxis = new NumberAxis();

        xAxis.setLabel("lambda (nm)");

        String yTitle = "F_lam, I_lam (/F_max)";
        yAxis.setLabel(yTitle);

        final LineChart<Number, Number> lineChartSpec
                = new LineChart<Number, Number>(xAxis, yAxis);

        lineChartSpec.setId("spec");

        String pattern = "0.00";
        //String pattern = "###.####";
        DecimalFormat myFormatter = new DecimalFormat(pattern);

        XYChart.Series seriesF = new XYChart.Series();
        seriesF.setName("Flux");

        XYChart.Series seriesI0 = new XYChart.Series();
        seriesI0.setName("I(theta=0)");
        String thetaStr = myFormatter.format((180.0 / Math.PI) * Math.acos(cosTheta[1][numMus - 2]));
        XYChart.Series seriesIn = new XYChart.Series();
        seriesIn.setName("I(theta=" + thetaStr + ")");

        // Only linear units will work (!) - divide flux by lambda_max:
        int[] iLamMinMax = MinMax2.minMax(flux);
        int iLamMax = iLamMinMax[1];
        double fMax = flux[0][iLamMax];
        String lamStr = myFormatter.format(1.0e7 * lambdaScale[iLamMax]);

        lineChartSpec.setTitle("Spectral energy distribution (SED)\r\n" + "lambda_Max: " + lamStr + " nm");

        for (int il = 0; il < numLams - 1; il++) {

            lambdaScale[il] = lambdaScale[il] * 1.0e7;
            //System.out.println("il, lambdaScale[il], flux[0][il], logE*flux[1][il]: " + il + " " + lambdaScale[il] + " " + flux[0][il] + " " + logE * flux[1][il]);
            //seriesF.getData().add( new XYChart.Data( Math.log10(lambdaScale[il]), logE*flux[1][il] ) );
            seriesF.getData().add(new XYChart.Data(lambdaScale[il], flux[0][il] / fMax));
            //seriesI0.getData().add( new XYChart.Data( Math.log10(lambdaScale[il]), Math.log10(intens[il][0]) ) );
            //seriesIn.getData().add( new XYChart.Data( Math.log10(lambdaScale[il]), Math.log10(intens[il][numMus-2]) ) );
            seriesI0.getData().add(new XYChart.Data(lambdaScale[il], intens[il][0] / fMax));
            seriesIn.getData().add(new XYChart.Data(lambdaScale[il], intens[il][numMus - 2] / fMax));
        }

        //Add UBVRI band centres
        //Photometric bands - Johnson UBVRI
        double[][][] filters = FilterSet.filterSet();
        int lam0_ptr = 11; // approximate band centre
        int numBands = filters.length;
        double[] lamUBVRI = new double[numBands];
        for (int ib = 0; ib < numBands; ib++) {
            lamUBVRI[ib] = filters[ib][0][lam0_ptr] * 1.0e7;
        }
        //System.out.println("lamUBVRI[0]: " + lamUBVRI[0][1]);

        XYChart.Series seriesU = new XYChart.Series();
        seriesU.setName("U");
        seriesU.getData().add(new XYChart.Data(lamUBVRI[0], 0.0));
        seriesU.getData().add(new XYChart.Data(lamUBVRI[0], 1.1));
        XYChart.Series seriesB = new XYChart.Series();
        seriesB.setName("B");
        seriesB.getData().add(new XYChart.Data(lamUBVRI[2], 0.0));
        seriesB.getData().add(new XYChart.Data(lamUBVRI[2], 1.1));
        XYChart.Series seriesV = new XYChart.Series();
        seriesV.setName("V");
        seriesV.getData().add(new XYChart.Data(lamUBVRI[3], 0.0));
        seriesV.getData().add(new XYChart.Data(lamUBVRI[3], 1.1));
        XYChart.Series seriesR = new XYChart.Series();
        seriesR.setName("R");
        seriesR.getData().add(new XYChart.Data(lamUBVRI[4], 0.0));
        seriesR.getData().add(new XYChart.Data(lamUBVRI[4], 1.1));
        XYChart.Series seriesI = new XYChart.Series();
        seriesI.setName("I");
        seriesI.getData().add(new XYChart.Data(lamUBVRI[5], 0.0));
        seriesI.getData().add(new XYChart.Data(lamUBVRI[5], 1.1));

        //lineChartSpec.getData().add(seriesSpec);
        lineChartSpec.getData().addAll(seriesF, seriesI0, seriesIn,
                seriesU, seriesB, seriesV, seriesR, seriesI);

        lineChartSpec.setCreateSymbols(false);

        return lineChartSpec;

    }

    /**
     * Plot 6: Normalized line profile (log(F_lambda), log(I_lambda(theta)):
     *
     * @param numPoints
     * @param lineLambdas
     * @param cosTheta
     * @param lineIntens
     * @param lineFlux
     * @return
     */
    public static LineChart<Number, Number> linePlot(int numPoints,
            double[] lineLambdas, double[][] cosTheta,
            double[][] lineIntens, double[][] lineFlux) {

        //Debug versions of the function with parameters for plotting up quqntities used in line profile
        // calculation:
        //public static LineChart<Number, Number> linePlot(int numPoints,
        //        double[] lineLambdas, double[][] cosTheta,
        //       double[][][] lineProf, double lam0) {
        //public static LineChart<Number, Number> linePlot(int numPoints,
        //        double[] lineLambdas, double[][] cosTheta,
        //        double[][] logKappaL, double lam0, double[][] kappa) {    
        //    public static LineChart<Number, Number> linePlot(int numPoints,
        //        double[] lineLambdas, double[][] cosTheta,
        //       double[][] logTauL, double lam0) {  
        double logE = Math.log10(Math.E);
        int numMus = cosTheta[0].length;

        final NumberAxis xAxis = new NumberAxis();
        final NumberAxis yAxis = new NumberAxis(0.0, 1.1, 0.2);

        xAxis.setLabel("Delta lambda (nm)");

        String yTitle = "Normalized F_lambda, I_lambda";
        yAxis.setLabel(yTitle);

        String pattern = "0.0";
        //String pattern = "###.####";
        DecimalFormat myFormatter = new DecimalFormat(pattern);

        final LineChart<Number, Number> lineChartLine
                = new LineChart<Number, Number>(xAxis, yAxis);
        lineChartLine.setId("line");

        int[] keyLambds = HalfPower.halfPower(numPoints, lineFlux);
        double lam0 = lineLambdas[keyLambds[0]];
        lam0 = lam0 * 1.0e7;
        String lamStr = myFormatter.format(lineLambdas[keyLambds[0]]);
        lineChartLine.setTitle("Spectral line\r\n" + "lambda_0: " + lamStr + " nm");
        XYChart.Series seriesFlux = new XYChart.Series();
        seriesFlux.setName("flux");

        //XYChart.Series seriesI0 = new XYChart.Series();
        //String thetaStr = myFormatter.format((180.0 / Math.PI) * Math.acos(cosTheta[1][5]));
        //seriesI0.setName("I(theta=" + thetaStr + " deg)");
        //XYChart.Series seriesIn = new XYChart.Series();
        //thetaStr = myFormatter.format((180.0 / Math.PI) * Math.acos(cosTheta[1][numMus - 2]));
        //seriesIn.setName("I(theta=" + thetaStr + " deg)");
        // Debug series:
        //XYChart.Series seriesPhi = new XYChart.Series();
        //seriesPhi.setName("phi");
        //XYChart.Series seriesKap = new XYChart.Series();
        //seriesKap.setName("kap");
        XYChart.Series seriesTau = new XYChart.Series();
        seriesTau.setName("tauL");

        int iCount = ((int) (numPoints / 2)) - 1; //initialize
        //System.out.println("numPoints " + numPoints + " iCount " + iCount);
        for (int il = ((int) (numPoints / 2)); il < numPoints; il++) {
//console.log("il " + il + " lineFlux[0][il]/lineFlux[0][numPoints - 1] " + lineFlux[0][il]/lineFlux[0][numPoints - 1]);
            if (lineFlux[0][il] < 0.95 * lineFlux[0][numPoints - 1]) {
//console.log("Condition triggered");
                iCount++;
            }
        }
        //System.out.println("iCount " + iCount);
//One to three more if we can accomodate them:
        //console.log("iCount " + iCount + " numPoints " + numPoints);
        if (iCount < numPoints - 1) {
            //console.log("First count++");
            iCount++;
        }
        //console.log("iCount " + iCount + " numPoints " + numPoints);
        if (iCount < numPoints - 1) {
            //console.log("Second count++");
            iCount++;
        }
        //console.log("iCount " + iCount + " numPoints " + numPoints);
        if (iCount < numPoints - 1) {
            //console.log("Third count++");
            iCount++;
        }
        //System.out.println("iCount " + iCount + " numPoints " + numPoints);

        int iStart = numPoints - iCount - 1;
        int iStop = iCount;
        //Set minimum range of x-axis to 0.1 nm:
        while ((lineLambdas[iStop] - lineLambdas[iStart]) < 1.0e-7 * 0.1) {
            iStart--;
            iStop++;
        }
        //System.out.println("iStart: " + iStart + " iStop: " + iStop);

        // CAUTION; numPoints-1st value holds the line centre monochromatic *continuum* flux for normalization
        for (int il = iStart; il <= iStop; il++) {

            lineLambdas[il] = lineLambdas[il] * 1.0e7;
            seriesFlux.getData().add(new XYChart.Data(lineLambdas[il] - lam0, lineFlux[0][il] / lineFlux[0][numPoints - 1]));
            //seriesIn.getData().add(new XYChart.Data(lineLambdas[il] - lam0, lineIntens[il][numMus - 2] / lineIntens[numPoints - 1][numMus - 2]));
            //seriesI0.getData().add(new XYChart.Data(lineLambdas[il] - lam0, lineIntens[il][5] / lineIntens[numPoints - 1][5]));
            //seriesI0.getData().add(new XYChart.Data(lineLambdas[il]-lam0, lineIntens[il][0]/1.0E+15 ));

            //Debug versions of the function with parameters for plotting up quqntities used in line profile
            // calculation:
            //seriesPhi.getData().add( new XYChart.Data(lineProf[0][il][0], Math.log10(lineProf[1][il][0])) );
            //seriesPhi.getData().add( new XYChart.Data(lineLambdas[il]-lam0, Math.log10(lineProf[1][il][0])) );
            //seriesKap.getData().add( new XYChart.Data(lineLambdas[il]-lam0, logE*logKappaL[il][7]) );
            //seriesKap.getData().add( new XYChart.Data(lineLambdas[il]-lam0, Math.log10(Math.exp(logKappaL[il][27]) + kappa[0][27]) ) );
            //seriesKap.getData().add( new XYChart.Data(lineLambdas[il]-lam0, Math.exp(logKappaL[il][27]) ) );
            //seriesTau.getData().add( new XYChart.Data(lineLambdas[il]-lam0, logE*logTauL[il][47]) );
            //System.out.println("il " + il + " lineLambdas[il]-lam0 " + (lineLambdas[il] - lam0)
            //        + " lineFlux[0][il] / lineFlux[0][numPoints - 1] " + lineFlux[0][il] / lineFlux[0][numPoints - 1]);
        }

        // Add line centre marker:
        //XYChart.Series seriesLam0 = new XYChart.Series();
        //seriesLam0.setName("lambda_0");
        //seriesLam0.getData().add(new XYChart.Data(lam0, 0.0));
        //seriesLam0.getData().add(new XYChart.Data(lam0, 1.05));
        //seriesLam0.getData().add(new XYChart.Data(0.0, 0.0));
        //seriesLam0.getData().add(new XYChart.Data(0.0, 1.0));
        // Add blue half-depth width vertical marker:
        //XYChart.Series seriesBlue2 = new XYChart.Series();
        //seriesBlue2.setName("lam_1/2");
        //seriesBlue2.getData().add(new XYChart.Data(lineLambdas[keyLambds[1]] - lam0, 0.0));
        //seriesBlue2.getData().add(new XYChart.Data(lineLambdas[keyLambds[1]] - lam0, 0.5));
        // Add red half-depth width vertical marker:
        //XYChart.Series seriesRed2 = new XYChart.Series();
        //seriesRed2.setName("lambda_1/2");
        //int index = numPoints - keyLambds[1];
        //int index = keyLambds[1];
        //System.out.println("LineCharts: numPoints, keyLambds[1], index: " + numPoints + " " + keyLambds[1] + " " + index);
        //seriesRed2.getData().add(new XYChart.Data(lineLambdas[index] - lam0, 0.0));
        //seriesRed2.getData().add(new XYChart.Data(lineLambdas[index] - lam0, lineFlux[0][index] / lineFlux[0][numPoints - 1]));
        //// Add half-power horizonal marker
        //XYChart.Series seriesHalf = new XYChart.Series();
        ////seriesRed2.setName("lam_1/2");
        //seriesHalf.getData().add(new XYChart.Data(lineLambdas[0] - lam0, lineFlux[0][index] / lineFlux[0][numPoints - 1]));
        //seriesHalf.getData().add(new XYChart.Data(lineLambdas[numPoints - 2] - lam0, lineFlux[0][index] / lineFlux[0][numPoints - 1]));
        // Add rectified continuum horizonal line
        XYChart.Series seriesCont = new XYChart.Series();
        seriesCont.setName(" ");
        seriesCont.getData().add(new XYChart.Data(lineLambdas[0] - lam0, 1.0));
        seriesCont.getData().add(new XYChart.Data(lineLambdas[numPoints - 2] - lam0, 1.0));

        //lineChartLine.getData().add(seriesFlux);
        //lineChartLine.getData().add(seriesI0);
        //lineChartLine.getData().addAll(seriesI0, seriesIn, seriesFlux, seriesLam0, seriesRed2, seriesCont);//, seriesHalf );
        lineChartLine.getData().addAll(seriesFlux); //, seriesCont);
        //        seriesLam0, seriesBlue2, seriesRed2, seriesHalf);

        //Debug versions of the function with parameters for plotting up quqntities used in line profile
        // calculation:
        //lineChartLine.getData().add(seriesPhi);
        //lineChartLine.getData().add(seriesKap);
        //lineChartLine.getData().add(seriesTau);
        //lineChartLine.setCreateSymbols(false);
        return lineChartLine;

    }

}
