// #include <string>
// #include <iostream>
// #include <fstream>
// #include <FileCtrl.hpp>
// #include <math.h>
// #include <vcl.h>
// #include <System.hpp>
// #pragma hdrstop
//
// #include "HistoPrepMain.h"
// #include "HistoAllThread.h"
// #pragma package(smart_init)
// using namespace std;
//
//
//
// Rcpp::NumericMatrix clean()
// {
//     ifstream InFileStream;
//     ifstream UHDRFileStream;
//     ofstream OutFileStream;
//     ofstream RHDRFileStream;
//     string s_buffer;
//     string s_path;
//     int n_rows, n_cols;
//     int i_nodata;
//     int i_row, i_col;
//     int i_rad;
//     int j_row, j_col;
//     double f_sum, f_med_sum, f_mod_max;
//     float f_cellsize;
//     float f_xllcorner, f_yllcorner; // coords of the bottom left corner of the lower left hand cell
//     float f_nodata;
//     float f_uns, f_res;
//     double d_res;
//     TCursor Save_Cursor;
//     int i_lock;
//     float *f_InRow_ptr;
//     float *f_OutRow_ptr;
//     int i_prop;
//     int i_scenario, n_scenarios;
//     AnsiString s_scenario;
//     AnsiString s_in_root;
//     AnsiString s_year, s_shortyear;
//     // AnsiString s_path,s_name;
//     int iAttributes = 0;
//     int i_flag;
//     TSearchRec sr;
//     int i_in, i_count;
//     double d_gaussian_wt[5][5] = {{1, 4, 7, 4, 1}, {4, 16, 26, 16, 4}, {7, 26, 41, 26, 7}, {4, 16, 26, 16, 4}, {1, 4, 7, 4, 1}};
//     double d_gauss_wt[3][3] = {{1, 2, 1}, {2, 4, 2}, {1, 2, 1}};
//     // TBitmap BitPlot;
//     int i_red, i_green, i_blue;
//     double d_frac, d_middle, d_max;
//     double d_sum_val, d_sum_wt;
//     double d_new_val;
//     double d_up_val, d_down_val;
//     int i_up_row, i_down_row;
//     double d_grid_max;
//     double d_bin_val;
//
//     // Load initial grid to 650 input grid
//     // Raw grids are 650 by 650
//     const int n_rows = 650;
//     const int n_cols = 650;
//
//
//     for (int i_row = 0; i_row < n_rows; i_row++)
//     {
//         for (i_col = 0; i_col < n_cols; i_col++)
//         {
//             d_in_grid[i_row][i_col] = 0;
//         }
//     }
//     // Read in grid
//     //  This is saved with 0,0 in bottom left [649][0], so observed (row) is inverted.
//     //  Put it back so [0][0] is 0,0
//     d_max = 0;
//     for (i_row = (n_rows - 1); i_row >= 0; i_row--)
//     {
//         for (i_col = 0; i_col < n_cols; i_col++)
//         {
//             InFileStream >> i_in;
//             d_in_grid[i_row][i_col] = i_in;
//             if (i_in > 0)
//             {
//                 if (i_in > d_max)
//                     d_max = i_in;
//             }
//         }
//     }
//     InFileStream.close();
//     // Special treatments!
//     // Zero zero has 20000 automatic identity comparisons where cells are
//     // compared against themselves
//     d_in_grid[0][0] = 0;
//     // add back in afterwards?
//
//     // Row two (near the top, has all the accumulated totals for stuff going
//     //  off the top. This bit is all cut off, but clean anyway, in case of
//     // future change of mind assume tail off
//     for (i_col = 0; i_col < n_cols; i_col++)
//     {
//         if (d_in_grid[647][i_col] > 0)
//         {
//             d_new_val = d_in_grid[2][i_col] / 2;
//             d_in_grid[648][i_col] = d_new_val;
//         }
//         else
//             d_in_grid[648][i_col] = 0;
//     }
//
//     // Effectively Trim to 400 * 400  (although the j indices can look beyond,
//     // so the gauss averaging can take place across  the whole 650 by 650 grid!
//     int n_rows = 400;
//     int n_cols = 400;
//
//     // Gaussian Kernel
//     d_max = 0;
//     d_grid_max = 0;
//     for (i_row = 0; i_row < n_rows; i_row++)
//     {
//         for (i_col = 0; i_col < n_cols; i_col++)
//         {
//             d_sum_val = 0;
//             d_sum_wt = 0;
//             for (j_row = -2; j_row < 3; j_row++)
//             {
//                 for (j_col = -2; j_col < 3; j_col++)
//                 {
//                     // Sum weights
//                     if ((i_row + j_row) < 0 || (i_col + j_col) < 0 || (i_row + j_row) > 649 || (i_col + j_col) > 649)
//                     {
//                         d_sum_wt += 0;
//                     } // end if kernel off grid
//                     else
//                     {
//                         d_sum_val += d_in_grid[i_row + j_row][i_col + j_col];
//                         d_sum_wt += d_gaussian_wt[j_row][j_col];
//                     } // end if kernel in grid
//                 } // end for j_col
//             } // enf for j_row
//             // Half zero for edge..full wt =273.. so only has any effect when d_sum_wt!=273
//             d_sum_wt = (273 + d_sum_wt) / 2;
//             d_new_grid[i_row][i_col] = d_sum_val / d_sum_wt;
//             if (d_new_grid[i_row][i_col] > d_max)
//                 d_max = d_new_grid[i_row][i_col];
//
//             if (i_row > 395)
//                 d_new_grid[i_row][i_col] = 0; // catch noise at the bottom
//         } // end for i_col
//     } // end for i_row
//
//     d_max = 0;
//     for (i_col = 0; i_col < n_cols; i_col++)
//     {
//         // Get sum in column (to 400)
//         f_sum = 0;
//         for (i_row = 0; i_row < n_rows; i_row++)
//         {
//             // Remove values below a threshold of 0.1        (Max 0.003% ish)
//             if (d_new_grid[i_row][i_col] <= 3)
//                 d_new_grid[i_row][i_col] = 0;
//             f_sum += d_new_grid[i_row][i_col];
//         }
//         // Now we have the sum (f_sum) of all rows in column i_col
//         // Loop through and normalise , and find median row  (sum normalised frequencies==0.5)
//         f_med_sum = 0; // zero at start of each column
//         f_mod_max = 0;
//         for (i_row = 0; i_row < n_rows; i_row++)
//         {
//             if (f_sum > 0)
//             {
//                 d_new_val = d_new_grid[i_row][i_col] / f_sum;
//                 d_new_grid[i_row][i_col] = d_new_val;
//                 if (d_new_val > d_grid_max)
//                     d_grid_max = d_new_val; // this will remain as this value unless we switch to the cumulative approaches (where it will become 1)
//                 if (HPMainWin->MedianRB->Checked)
//                 { // Calculate Median
//                     if (f_med_sum <= 0.5)
//                     {
//                         f_med_sum += d_new_grid[i_row][i_col];
//                         i_med_row[i_col] = i_row; // row in which median falls (beyond this
//                                                   //	f_med_gap[i_col]=0.5-f_med_sum; // how far from current row to the true median ( this is eating into i_med_row[]
//                     }
//                 }
//                 else
//                 { // Calculate Mode
//                     if (d_new_grid[i_row][i_col] > f_mod_max)
//                     {
//                         f_mod_max = d_new_grid[i_row][i_col]; // mode is tallest
//                         i_med_row[i_col] = i_row;             // row in which median falls (beyond this
//                     }
//                 }
//             }
//             if (d_new_grid[i_row][i_col] > d_max)
//                 d_max = d_new_grid[i_row][i_col];
//         }
//         f_sum = 0;
//         for (i_row = 0; i_row < n_rows; i_row++)
//         {
//             f_sum += d_new_grid[i_row][i_col];
//         }
//         if (f_sum > 0)
//             f_sum = 0;
//     } // end for i_col
//
//
// }
