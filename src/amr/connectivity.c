//!==================================================================================================================================
//! Copyright (c) 2018 - 2020 Alexander Astanin
//!
//! This file is part of FLUXO (github.com/project-fluxo/fluxo). FLUXO is free software: you can redistribute it and/or modify
//! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
//! of the License, or (at your option) any later version.
//!
//! FLUXO is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
//! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
//!
//! You should have received a copy of the GNU General Public License along with FLUXO. If not, see <http://www.gnu.org/licenses/>.
//!==================================================================================================================================
#include "connectivity.h"

const int H2P_FaceNodeMap[7][5] =
    {{-1, -1, -1, -1, -1},
     {-1, 0, 2, 3, 1},
     {-1, 0, 1, 3, 2},
     {-1, 0, 1, 3, 2},
     {-1, 1, 0, 2, 3},
     {-1, 0, 2, 3, 1},
     {-1, 0, 1, 3, 2}};
//! mapping from P4EST node of local sides to HOPEST node of local sides
const int P2H_FaceNodeMap[6][4] =
    {{1, 4, 2, 3},
     {1, 2, 4, 3},
     {1, 2, 4, 3},
     {2, 1, 3, 4},
     {1, 4, 2, 3},
     {1, 2, 4, 3}};

const int Rmatrix[6][6] =
    {{0, 1, 1, 0, 0, 1},
     {2, 0, 0, 1, 1, 0},
     {2, 0, 0, 1, 1, 0},
     {0, 2, 2, 0, 0, 1},
     {0, 2, 2, 0, 0, 1},
     {2, 0, 0, 2, 2, 0}};

const int Qmatrix[3][4] =
    {{1, 2, 5, 6},
     {0, 3, 4, 7},
     {0, 4, 3, 7}};

const int Pmatrix[8][4] =
    {{0, 1, 2, 3},
     {0, 2, 1, 3},
     {1, 0, 3, 2},
     {1, 3, 0, 2},
     {2, 0, 3, 1},
     {2, 3, 0, 1},
     {3, 1, 2, 0},
     {3, 2, 1, 0}};
const int P2H_side[] = {5, 3, 2, 4, 1, 6};
const int H2P_side[] = {-1, 4, 2, 1, 3, 0, 5}; //"-1" - There is no H2P_side[0] in HOPR

const int H_MortarCase[5][5] = //  first CGNS node and second CGNS node->Mortar Case [1:8]
    {{-1, -1, -1, -1, -1},
     {-1, 0, 1, 0, 2},  //           ! (1,2)->1, (1,4)->2
     {-1, 3, 0, 4, 0},  //                         ! (2,1)->3, (2,3)->4
     {-1, 0, 5, 0, 6},  //                       ! (3,2)->5, (3,4)->6
     {-1, 7, 0, 8, 0}}; //              ! (4,1)->7, (4,3)->8


const int P2H_MortarMap[9][4] = //   !p4est mortar ID, MortarCase -> iMortar CGNS
    {{-1, -1, -1, -1},          //! iMortar = P2H_MortarMap(iPMortar, H_MortarCase( node1, node2) )
     {1, 2, 3, 4},
     {1, 3, 2, 4},
     {2, 1, 4, 3},
     {2, 4, 1, 3},
     {4, 2, 3, 1},
     {4, 3, 2, 1},
     {3, 1, 4, 2},
     {3, 4, 1, 2}};

int GetHFlip(int PSide0, int PSide1, int PFlip)
{
       // !1. Get CGNS side from P4 side
    
    int HSide0 = P2H_side[PSide0];
    // !2. First node CGNS -> p4est
    int PNode0 = H2P_FaceNodeMap[HSide0][1];
    // !3. Get oriented node on neighbour side, formula and matrices see paper Burstedde p4est, 2011
    int PNode1 = Pmatrix[ Qmatrix[ Rmatrix[PSide0][PSide1] ][PFlip] ][PNode0];
    // !4. P4EST node -> CGNS
    int GetHFlip = P2H_FaceNodeMap[PSide1][PNode1];
    // return PNode1;
    return GetHFlip;
}

int GetHMortar(int PMortar, int PSide, int PnbSide, int PFlip)
{

    int PNodeA = Pmatrix[Qmatrix[Rmatrix[PnbSide][PSide]][PFlip]][0];
    int PNodeB = Pmatrix[Qmatrix[Rmatrix[PnbSide][PSide]][PFlip]][1];
    int HNode1 = P2H_FaceNodeMap[PSide][PNodeA];
    int HNode2 = P2H_FaceNodeMap[PSide][PNodeB];
    int GetHMortar = P2H_MortarMap[H_MortarCase[HNode1][HNode2]][PMortar];
    return GetHMortar;
}
