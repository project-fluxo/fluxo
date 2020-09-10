#ifndef CONNECTIVITY_H

#define CONNECTIVITY_H
//!mapping from HOPEST node of local sides to P4EST nodes of local sides
// int const H2P_FaceNodeMap[1 : 4, 1 : 6) = &RESHAPE((/ 0, 2, 3, 1, &0, 1, 3, 2, &0, 1, 3, 2, &1, 0, 2, 3, &0, 2, 3, 1, &0, 1, 3, 2 /), (/ 4, 6 /))
//int const H2P_FaceNodeMap[1 : 4, 1 : 6) = &RESHAPE((/ 0, 2, 3, 1, &0, 1, 3, 2, &0, 1, 3, 2, &1, 0, 2, 3, &0, 2, 3, 1, &0, 1, 3, 2 /), (/ 4, 6 /))

extern const int H2P_FaceNodeMap[7][5];
//! mapping from P4EST node of local sides to HOPEST node of local sides
extern const int P2H_FaceNodeMap[6][4];

extern const int Rmatrix[6][6];

extern const int Qmatrix[3][4];

extern const int Pmatrix[8][4];
extern const int P2H_side[];
extern const int H2P_side[]; //"-1" - There is no H2P_side[0] in HOPR

extern const int H_MortarCase[5][5]; //              ! (4,1)->7, (4,3)->8


extern const int P2H_MortarMap[9][4];

int GetHFlip(int PSide0, int PSide1, int PFlip);

int GetHMortar(int PMortar, int PSide, int PnbSide, int PFlip);


#endif