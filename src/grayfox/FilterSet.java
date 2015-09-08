/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package grayfox;

/**
 *
 * @author Ian Put definition of filter transparency curves in separate function
 * so it can be called by plotting part of program, as well as by photometry
 * routine
 *
 */
public class FilterSet {

    public static double[][][] filterSet() {

        int numBands = 6; // Bessell-Johnson UxBxBVRI
        int numLambdaFilt = 25;  //test for now

        double[][][] filterCurves = new double[numBands][2][numLambdaFilt];

        //Initialize all filterCurves - the real data below won't fill in all the array elements:
        for (int ib = 0; ib < numBands; ib++) {
            for (int il = 0; il < numLambdaFilt; il++) {
                filterCurves[ib][0][il] = 1000.0; //placeholder wavelength (nm)
                filterCurves[ib][1][il] = 0.0e0;  // initialize filter transparency to 0.0
            }
        }

        /*
         //Test: central lambda and half power hald-widths:
        
         //U: filter 0: //lambda0=360, Delta_lambda_1/2=35
         filterCurves[0][0][0] = 325.0;
         filterCurves[0][1][0] = 0.5;
         filterCurves[0][0][1] = 360.0;
         filterCurves[0][1][1] = 1.0;
         filterCurves[0][0][2] = 395.0;
         filterCurves[0][1][2] = 0.5;

         //B: filter 1: //lambda0=440, Delta_lambda_1/2=50
         filterCurves[1][0][0] = 390.0;
         filterCurves[1][1][0] = 0.5;
         filterCurves[1][0][1] = 440.0;
         filterCurves[1][1][1] = 1.0;
         filterCurves[1][0][2] = 490.0;
         filterCurves[1][1][2] = 0.5;

         //V: filter 2: //lambda0=550, Delta_lambda_1/2=45
         filterCurves[2][0][0] = 505.0;
         filterCurves[2][1][0] = 0.5;
         filterCurves[2][0][1] = 550.0;
         filterCurves[2][1][1] = 1.0;
         filterCurves[2][0][2] = 595.0;
         filterCurves[2][1][2] = 0.5;

         //R: filter 3: //lambda0=640, Delta_lambda_1/2=87.5
         filterCurves[3][0][0] = 552.5;
         filterCurves[3][1][0] = 0.5;
         filterCurves[3][0][1] = 640.0;
         filterCurves[3][1][1] = 1.0;
         filterCurves[3][0][2] = 727.5;
         filterCurves[3][1][2] = 0.5;

         //I: filter 4: //lambda0=790, Delta_lambda_1/2=70
         filterCurves[4][0][0] = 720.0;
         filterCurves[4][1][0] = 0.5;
         filterCurves[4][0][1] = 790.0;
         filterCurves[4][1][1] = 1.0;
         filterCurves[4][0][2] = 860.0;
         filterCurves[4][1][2] = 0.5;
         */
//http://ulisse.pd.astro.it/Astro/ADPS/Systems/Sys_136/index_136.html
//Bessell, M. S., 1990, PASP, 102, 1181
//photometric filter data for Bessell UxBxBVRI system from Asiago database in Java & JavaScript syntax
//Individual bands are below master table
//        UX
        filterCurves[0][0][0] = 300.0;
        filterCurves[0][1][0] = 0.000;
        filterCurves[0][0][1] = 305.0;
        filterCurves[0][1][1] = 0.016;
        filterCurves[0][0][2] = 310.0;
        filterCurves[0][1][2] = 0.068;
        filterCurves[0][0][3] = 315.0;
        filterCurves[0][1][3] = 0.167;
        filterCurves[0][0][4] = 320.0;
        filterCurves[0][1][4] = 0.287;
        filterCurves[0][0][5] = 325.0;
        filterCurves[0][1][5] = 0.423;
        filterCurves[0][0][6] = 330.0;
        filterCurves[0][1][6] = 0.560;
        filterCurves[0][0][7] = 335.0;
        filterCurves[0][1][7] = 0.673;
        filterCurves[0][0][8] = 340.0;
        filterCurves[0][1][8] = 0.772;
        filterCurves[0][0][9] = 345.0;
        filterCurves[0][1][9] = 0.841;
        filterCurves[0][0][10] = 350.0;
        filterCurves[0][1][10] = 0.905;
        filterCurves[0][0][11] = 355.0;
        filterCurves[0][1][11] = 0.943;
        filterCurves[0][0][12] = 360.0;
        filterCurves[0][1][12] = 0.981;
        filterCurves[0][0][13] = 365.0;
        filterCurves[0][1][13] = 0.993;
        filterCurves[0][0][14] = 370.0;
        filterCurves[0][1][14] = 1.000;
        filterCurves[0][0][15] = 375.0;
        filterCurves[0][1][15] = 0.989;
        filterCurves[0][0][16] = 380.0;
        filterCurves[0][1][16] = 0.916;
        filterCurves[0][0][17] = 385.0;
        filterCurves[0][1][17] = 0.804;
        filterCurves[0][0][18] = 390.0;
        filterCurves[0][1][18] = 0.625;
        filterCurves[0][0][19] = 395.0;
        filterCurves[0][1][19] = 0.423;
        filterCurves[0][0][20] = 400.0;
        filterCurves[0][1][20] = 0.238;
        filterCurves[0][0][21] = 405.0;
        filterCurves[0][1][21] = 0.114;
        filterCurves[0][0][22] = 410.0;
        filterCurves[0][1][22] = 0.051;
        filterCurves[0][0][23] = 415.0;
        filterCurves[0][1][23] = 0.019;
        filterCurves[0][0][24] = 420.0;
        filterCurves[0][1][24] = 0.000;

//BX
        filterCurves[1][0][0] = 360.0;
        filterCurves[1][1][0] = 0.000;
        filterCurves[1][0][1] = 370.0;
        filterCurves[1][1][1] = 0.026;
        filterCurves[1][0][2] = 380.0;
        filterCurves[1][1][2] = 0.120;
        filterCurves[1][0][3] = 390.0;
        filterCurves[1][1][3] = 0.523;
        filterCurves[1][0][4] = 400.0;
        filterCurves[1][1][4] = 0.875;
        filterCurves[1][0][5] = 410.0;
        filterCurves[1][1][5] = 0.956;
        filterCurves[1][0][6] = 420.0;
        filterCurves[1][1][6] = 1.000;
        filterCurves[1][0][7] = 430.0;
        filterCurves[1][1][7] = 0.998;
        filterCurves[1][0][8] = 440.0;
        filterCurves[1][1][8] = 0.972;
        filterCurves[1][0][9] = 450.0;
        filterCurves[1][1][9] = 0.901;
        filterCurves[1][0][10] = 460.0;
        filterCurves[1][1][10] = 0.793;
        filterCurves[1][0][11] = 470.0;
        filterCurves[1][1][11] = 0.694;
        filterCurves[1][0][12] = 480.0;
        filterCurves[1][1][12] = 0.587;
        filterCurves[1][0][13] = 490.0;
        filterCurves[1][1][13] = 0.470;
        filterCurves[1][0][14] = 500.0;
        filterCurves[1][1][14] = 0.362;
        filterCurves[1][0][15] = 510.0;
        filterCurves[1][1][15] = 0.263;
        filterCurves[1][0][16] = 520.0;
        filterCurves[1][1][16] = 0.169;
        filterCurves[1][0][17] = 530.0;
        filterCurves[1][1][17] = 0.107;
        filterCurves[1][0][18] = 540.0;
        filterCurves[1][1][18] = 0.049;
        filterCurves[1][0][19] = 550.0;
        filterCurves[1][1][19] = 0.010;
        filterCurves[1][0][20] = 560.0;
        filterCurves[1][1][20] = 0.000;
        filterCurves[1][0][21] = 560.0;
        filterCurves[1][1][21] = 0.000;
        filterCurves[1][0][22] = 560.0;
        filterCurves[1][1][22] = 0.000;
        filterCurves[1][0][23] = 560.0;
        filterCurves[1][1][23] = 0.000;
        filterCurves[1][0][24] = 560.0;
        filterCurves[1][1][24] = 0.000;

//B
        filterCurves[2][0][0] = 360.0;
        filterCurves[2][1][0] = 0.000;
        filterCurves[2][0][1] = 370.0;
        filterCurves[2][1][1] = 0.030;
        filterCurves[2][0][2] = 380.0;
        filterCurves[2][1][2] = 0.134;
        filterCurves[2][0][3] = 390.0;
        filterCurves[2][1][3] = 0.567;
        filterCurves[2][0][4] = 400.0;
        filterCurves[2][1][4] = 0.920;
        filterCurves[2][0][5] = 410.0;
        filterCurves[2][1][5] = 0.978;
        filterCurves[2][0][6] = 420.0;
        filterCurves[2][1][6] = 1.000;
        filterCurves[2][0][7] = 430.0;
        filterCurves[2][1][7] = 0.978;
        filterCurves[2][0][8] = 440.0;
        filterCurves[2][1][8] = 0.935;
        filterCurves[2][0][9] = 450.0;
        filterCurves[2][1][9] = 0.853;
        filterCurves[2][0][10] = 460.0;
        filterCurves[2][1][10] = 0.740;
        filterCurves[2][0][11] = 470.0;
        filterCurves[2][1][11] = 0.640;
        filterCurves[2][0][12] = 480.0;
        filterCurves[2][1][12] = 0.536;
        filterCurves[2][0][13] = 490.0;
        filterCurves[2][1][13] = 0.424;
        filterCurves[2][0][14] = 500.0;
        filterCurves[2][1][14] = 0.325;
        filterCurves[2][0][15] = 510.0;
        filterCurves[2][1][15] = 0.235;
        filterCurves[2][0][16] = 520.0;
        filterCurves[2][1][16] = 0.150;
        filterCurves[2][0][17] = 530.0;
        filterCurves[2][1][17] = 0.095;
        filterCurves[2][0][18] = 540.0;
        filterCurves[2][1][18] = 0.043;
        filterCurves[2][0][19] = 550.0;
        filterCurves[2][1][19] = 0.009;
        filterCurves[2][0][20] = 560.0;
        filterCurves[2][1][20] = 0.000;
        filterCurves[2][0][21] = 560.0;
        filterCurves[2][1][21] = 0.000;
        filterCurves[2][0][22] = 560.0;
        filterCurves[2][1][22] = 0.000;
        filterCurves[2][0][23] = 560.0;
        filterCurves[2][1][23] = 0.000;
        filterCurves[2][0][24] = 560.0;
        filterCurves[2][1][24] = 0.000;

//V
        filterCurves[3][0][0] = 470.0;
        filterCurves[3][1][0] = 0.000;
        filterCurves[3][0][1] = 480.0;
        filterCurves[3][1][1] = 0.030;
        filterCurves[3][0][2] = 490.0;
        filterCurves[3][1][2] = 0.163;
        filterCurves[3][0][3] = 500.0;
        filterCurves[3][1][3] = 0.458;
        filterCurves[3][0][4] = 510.0;
        filterCurves[3][1][4] = 0.780;
        filterCurves[3][0][5] = 520.0;
        filterCurves[3][1][5] = 0.967;
        filterCurves[3][0][6] = 530.0;
        filterCurves[3][1][6] = 1.000;
        filterCurves[3][0][7] = 540.0;
        filterCurves[3][1][7] = 0.973;
        filterCurves[3][0][8] = 550.0;
        filterCurves[3][1][8] = 0.898;
        filterCurves[3][0][9] = 560.0;
        filterCurves[3][1][9] = 0.792;
        filterCurves[3][0][10] = 570.0;
        filterCurves[3][1][10] = 0.684;
        filterCurves[3][0][11] = 580.0;
        filterCurves[3][1][11] = 0.574;
        filterCurves[3][0][12] = 590.0;
        filterCurves[3][1][12] = 0.461;
        filterCurves[3][0][13] = 600.0;
        filterCurves[3][1][13] = 0.359;
        filterCurves[3][0][14] = 610.0;
        filterCurves[3][1][14] = 0.270;
        filterCurves[3][0][15] = 620.0;
        filterCurves[3][1][15] = 0.197;
        filterCurves[3][0][16] = 630.0;
        filterCurves[3][1][16] = 0.135;
        filterCurves[3][0][17] = 640.0;
        filterCurves[3][1][17] = 0.081;
        filterCurves[3][0][18] = 650.0;
        filterCurves[3][1][18] = 0.045;
        filterCurves[3][0][19] = 660.0;
        filterCurves[3][1][19] = 0.025;
        filterCurves[3][0][20] = 670.0;
        filterCurves[3][1][20] = 0.017;
        filterCurves[3][0][21] = 680.0;
        filterCurves[3][1][21] = 0.013;
        filterCurves[3][0][22] = 690.0;
        filterCurves[3][1][22] = 0.009;
        filterCurves[3][0][23] = 700.0;
        filterCurves[3][1][23] = 0.000;
        filterCurves[3][0][24] = 700.0;
        filterCurves[3][1][24] = 0.000;

//R
        filterCurves[4][0][0] = 550.0;
        filterCurves[4][1][0] = 0.00;
        filterCurves[4][0][1] = 560.0;
        filterCurves[4][1][1] = 0.23;
        filterCurves[4][0][2] = 570.0;
        filterCurves[4][1][2] = 0.74;
        filterCurves[4][0][3] = 580.0;
        filterCurves[4][1][3] = 0.91;
        filterCurves[4][0][4] = 590.0;
        filterCurves[4][1][4] = 0.98;
        filterCurves[4][0][5] = 600.0;
        filterCurves[4][1][5] = 1.00;
        filterCurves[4][0][6] = 610.0;
        filterCurves[4][1][6] = 0.98;
        filterCurves[4][0][7] = 620.0;
        filterCurves[4][1][7] = 0.96;
        filterCurves[4][0][8] = 630.0;
        filterCurves[4][1][8] = 0.93;
        filterCurves[4][0][9] = 640.0;
        filterCurves[4][1][9] = 0.90;
        filterCurves[4][0][10] = 650.0;
        filterCurves[4][1][10] = 0.86;
        filterCurves[4][0][11] = 660.0;
        filterCurves[4][1][11] = 0.81;
        filterCurves[4][0][12] = 670.0;
        filterCurves[4][1][12] = 0.78;
        filterCurves[4][0][13] = 680.0;
        filterCurves[4][1][13] = 0.72;
        filterCurves[4][0][14] = 690.0;
        filterCurves[4][1][14] = 0.67;
        filterCurves[4][0][15] = 700.0;
        filterCurves[4][1][15] = 0.61;
        filterCurves[4][0][16] = 710.0;
        filterCurves[4][1][16] = 0.56;
        filterCurves[4][0][17] = 720.0;
        filterCurves[4][1][17] = 0.51;
        filterCurves[4][0][18] = 730.0;
        filterCurves[4][1][18] = 0.46;
        filterCurves[4][0][19] = 740.0;
        filterCurves[4][1][19] = 0.40;
        filterCurves[4][0][20] = 750.0;
        filterCurves[4][1][20] = 0.35;
        filterCurves[4][0][21] = 800.0;
        filterCurves[4][1][21] = 0.14;
        filterCurves[4][0][22] = 850.0;
        filterCurves[4][1][22] = 0.03;
        filterCurves[4][0][23] = 900.0;
        filterCurves[4][1][23] = 0.00;
        filterCurves[4][0][24] = 900.0;
        filterCurves[4][1][24] = 0.000;

//I
        filterCurves[5][0][0] = 700.0;
        filterCurves[5][1][0] = 0.000;
        filterCurves[5][0][1] = 710.0;
        filterCurves[5][1][1] = 0.024;
        filterCurves[5][0][2] = 720.0;
        filterCurves[5][1][2] = 0.232;
        filterCurves[5][0][3] = 730.0;
        filterCurves[5][1][3] = 0.555;
        filterCurves[5][0][4] = 740.0;
        filterCurves[5][1][4] = 0.785;
        filterCurves[5][0][5] = 750.0;
        filterCurves[5][1][5] = 0.910;
        filterCurves[5][0][6] = 760.0;
        filterCurves[5][1][6] = 0.965;
        filterCurves[5][0][7] = 770.0;
        filterCurves[5][1][7] = 0.985;
        filterCurves[5][0][8] = 780.0;
        filterCurves[5][1][8] = 0.990;
        filterCurves[5][0][9] = 790.0;
        filterCurves[5][1][9] = 0.995;
        filterCurves[5][0][10] = 800.0;
        filterCurves[5][1][10] = 1.000;
        filterCurves[5][0][11] = 810.0;
        filterCurves[5][1][11] = 1.000;
        filterCurves[5][0][12] = 820.0;
        filterCurves[5][1][12] = 0.990;
        filterCurves[5][0][13] = 830.0;
        filterCurves[5][1][13] = 0.980;
        filterCurves[5][0][14] = 840.0;
        filterCurves[5][1][14] = 0.950;
        filterCurves[5][0][15] = 850.0;
        filterCurves[5][1][15] = 0.910;
        filterCurves[5][0][16] = 860.0;
        filterCurves[5][1][16] = 0.860;
        filterCurves[5][0][17] = 870.0;
        filterCurves[5][1][17] = 0.750;
        filterCurves[5][0][18] = 880.0;
        filterCurves[5][1][18] = 0.560;
        filterCurves[5][0][19] = 890.0;
        filterCurves[5][1][19] = 0.330;
        filterCurves[5][0][20] = 900.0;
        filterCurves[5][1][20] = 0.150;
        filterCurves[5][0][21] = 910.0;
        filterCurves[5][1][21] = 0.030;
        filterCurves[5][0][22] = 920.0;
        filterCurves[5][1][22] = 0.000;
        filterCurves[5][0][23] = 920.0;
        filterCurves[5][1][23] = 0.000;
        filterCurves[5][0][24] = 920.0;
        filterCurves[5][1][24] = 0.000;

        for (int ib = 0; ib < numBands; ib++) {
//wavelength loop is over photometric filter data wavelengths
            for (int il = 0; il < numLambdaFilt; il++) {
                filterCurves[ib][0][il] = filterCurves[ib][0][il] * 1.0e-7; // nm to cm
            }
        }

        return filterCurves;
    }

}
