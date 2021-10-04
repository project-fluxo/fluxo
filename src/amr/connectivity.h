//!=================================================================================================================================
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
//!=================================================================================================================================
#ifndef CONNECTIVITY_H

#define CONNECTIVITY_H

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
