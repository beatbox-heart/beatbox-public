/**
 * Copyright (C) (2010-2016) Vadim Biktashev, Irina Biktasheva et al. 
 * (see ../AUTHORS for the full list of contributors)
 *
 * This file is part of Beatbox.
 *
 * Beatbox is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Beatbox is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Beatbox.  If not, see <http://www.gnu.org/licenses/>.
 */

/* Disclaimer: This file is not properly tested */


//---------------------------------------------------------------------------
// Constants. If really required, you can type in the units of the constants
// in place of the , 1 ) as indicated in the comments.
//---------------------------------------------------------------------------

_(   CaMKo  ,  0.05,  1  )  // dimensionless (in CaMK)
_(   KmCaM  , 0.0015,  1  )  // millimolar (in CaMK)
_(   KmCaMK  , 0.15,  1  )  // millimolar (in CaMK)
_(   aCaMK  , 0.05,  1  )  // per_millimolar_per_millisecond (in CaMK)
_(   bCaMK  , 0.00068,  1  )  // per_millisecond (in CaMK)
_(   Kmn  , 0.002,  1  )  // millimolar (in ICaL)
_(   PCa_b  , 0.0001,  1  )  // dimensionless (in ICaL)
_(   k2n  , 1000.0,  1  )  // per_millisecond (in ICaL)
_(   PCab  , 2.5e-8,  1  )  // milliS_per_microF (in ICab)
_(   GK1_b  , 0.1908,  1  )  // milliS_per_microF (in IK1)
_(   GKb_b  , 0.003,  1  )  // milliS_per_microF (in IKb)
_(   GKr_b  , 0.046,  1  )  // milliS_per_microF (in IKr)
_(   GKs_b  , 0.0034,  1  )  // milliS_per_microF (in IKs)
_(   Gncx_b  , 0.0008,  1  )  // milliS_per_microF (in INaCa_i)
_(   KmCaAct  , 150.0e-6,  1  )  // millimolar (in INaCa_i)
_(   kasymm  , 12.5,  1  )  // dimensionless (in INaCa_i)
_(   kcaoff  , 5.0e3,  1  )  // per_millisecond (in INaCa_i)
_(   kcaon  , 1.5e6,  1  )  // per_millisecond (in INaCa_i)
_(   kna1  , 15.0,  1  )  // per_millisecond (in INaCa_i)
_(   kna2  , 5.0,  1  )  // per_millisecond (in INaCa_i)
_(   kna3  , 88.12,  1  )  // per_millisecond (in INaCa_i)
_(   qca  , 0.167,  1  )  // dimensionless (in INaCa_i)
_(   qna  , 0.5224,  1  )  // dimensionless (in INaCa_i)
_(   wca  , 6.0e4,  1  )  // dimensionless (in INaCa_i)
_(   wna  , 6.0e4,  1  )  // dimensionless (in INaCa_i)
_(   wnaca  , 5.0e3,  1  )  // dimensionless (in INaCa_i)
_(   H  , 1.0e-7,  1  )  // millimolar (in INaK)
_(   Khp  , 1.698e-7,  1  )  // millimolar (in INaK)
_(   Kki  , 0.5,  1  )  // per_millisecond (in INaK)
_(   Kko  , 0.3582,  1  )  // per_millisecond (in INaK)
_(   Kmgatp  , 1.698e-7,  1  )  // millimolar (in INaK)
_(   Knai0  , 9.073,  1  )  // millimolar (in INaK)
_(   Knao0  , 27.78,  1  )  // millimolar (in INaK)
_(   Knap  , 224.0,  1  )  // millimolar (in INaK)
_(   Kxkur  , 292.0,  1  )  // millimolar (in INaK)
_(   MgADP  , 0.05,  1  )  // millimolar (in INaK)
_(   MgATP  , 9.8,  1  )  // millimolar (in INaK)
_(   Pnak_b  , 30.0,  1  )  // milliS_per_microF (in INaK)
_(   delta  , -0.155,  1  )  // millivolt (in INaK)
_(   eP  , 4.2,  1  )  // dimensionless (in INaK)
_(   k1m  , 182.4,  1  )  // per_millisecond (in INaK)
_(   k1p  , 949.5,  1  )  // per_millisecond (in INaK)
_(   k2m  , 39.4,  1  )  // per_millisecond (in INaK)
_(   k2p  , 687.2,  1  )  // per_millisecond (in INaK)
_(   k3m  , 79300.0,  1  )  // per_millisecond (in INaK)
_(   k3p  , 1899.0,  1  )  // per_millisecond (in INaK)
_(   k4m  , 40.0,  1  )  // per_millisecond (in INaK)
_(   k4p  , 639.0,  1  )  // per_millisecond (in INaK)
_(   GNaL_b  , 0.0075,  1  )  // milliS_per_microF (in INaL)
_(   thL  , 200.0,  1  )  // millisecond (in INaL)
_(   PNab  , 3.75e-10,  1  )  // milliS_per_microF (in INab)
_(   Ahf  , 0.99,  1  )  // dimensionless (in INa)
_(   GNa  , 75.0,  1  )  // milliS_per_microF (in INa)
_(   hssV1  , 82.9,  1  )  // millivolt (in INa)
_(   hssV2  , 6.086,  1  )  // millivolt (in INa)
_(   mssV1  , 39.57,  1  )  // millivolt (in INa)
_(   mssV2  , 9.871,  1  )  // millivolt (in INa)
_(   mtD1  , 6.765,  1  )  // dimensionless (in INa)
_(   mtD2  , 8.552,  1  )  // dimensionless (in INa)
_(   mtV1  , 11.64,  1  )  // millivolt (in INa)
_(   mtV2  , 34.77,  1  )  // millivolt (in INa)
_(   mtV3  , 77.42,  1  )  // millivolt (in INa)
_(   mtV4  , 5.955,  1  )  // millivolt (in INa)
_(   GpCa  , 0.0005,  1  )  // milliS_per_microF (in IpCa)
_(   KmCap  , 0.0005,  1  )  // millimolar (in IpCa)
_(   Gto_b  , 0.02,  1  )  // milliS_per_microF (in Ito)
_(   L  , 0.01,  1  )  // centimeter (in cell_geometry)
_(   rad  , 0.0011,  1  )  // centimeter (in cell_geometry)
_(   celltype  , 0,  1  )  // dimensionless (in environment)
_(   cao  , 1.8,  1  )  // millimolar (in extracellular)
_(   ko  , 5.4,  1  )  // millimolar (in extracellular)
_(   nao  , 140.0,  1  )  // millimolar (in extracellular)
_(   BSLmax  , 1.124,  1  )  // millimolar (in intracellular_ions)
_(   BSRmax  , 0.047,  1  )  // millimolar (in intracellular_ions)
_(   KmBSL  , 0.0087,  1  )  // millimolar (in intracellular_ions)
_(   KmBSR  , 0.00087,  1  )  // millimolar (in intracellular_ions)
_(   cm  , 1.0,  1  )  // microF_per_centimeter_squared (in intracellular_ions)
_(   cmdnmax_b  , 0.05,  1  )  // millimolar (in intracellular_ions)
_(   csqnmax  , 10.0,  1  )  // millimolar (in intracellular_ions)
_(   kmcmdn  , 0.00238,  1  )  // millimolar (in intracellular_ions)
_(   kmcsqn  , 0.8,  1  )  // millimolar (in intracellular_ions)
_(   kmtrpn  , 0.0005,  1  )  // millimolar (in intracellular_ions)
_(   trpnmax  , 0.07,  1  )  // millimolar (in intracellular_ions)
_(   amp  , -80.0,  1  )  // microA_per_microF (in membrane)
_(   duration  , 0.5,  1  )  // millisecond (in membrane)
_(   F  , 96485.0,  1  )  // coulomb_per_mole (in physical_constants)
_(   R  , 8314.0,  1  )  // joule_per_kilomole_kelvin (in physical_constants)
_(   T  , 310.0,  1  )  // kelvin (in physical_constants)
_(   zca  , 2.0,  1  )  // dimensionless (in physical_constants)
_(   zk  , 1.0,  1  )  // dimensionless (in physical_constants)
_(   zna  , 1.0,  1  )  // dimensionless (in physical_constants)
_(   PKNa  , 0.01833,  1  )  // dimensionless (in reversal_potentials)
_(   bt  , 4.75,  1  )  // millisecond (in ryr)
