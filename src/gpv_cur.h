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


  _(ENa ,     RT_F*log(Nao/Nai))
  _(EK  ,     RT_F*log(Ko/Ki))
  _(ECa , 0.5*RT_F*log(Cao/Cai))
  _(Emh ,     RT_F*log((Nao+0.12*Ko)/(Nai+0.12*Ki)))
  _(Ik    , IKx / 140. * (Ki - Ko * exp(-V/RT_F)) * x )
  _(Iki   , GKi * (V - EK) * Ko/(Ko+kmKi) / (1 + exp( 2.*(V-EK+10.-Vsh)/RT_F )))
  _(Ito   , Gto * (V - EK) * q*r)
  _(IbK   , GbK * (V - EK))
  _(INa   , GNa * (V - Emh) * cub(m)*h)
  _(IbNa  , GbNa* (V - ENa))
  _(IbCa  , GbCa* (V - ECa))
  _(IsiCa ,     4*PCa*(V-50.)/RT_F / (1-exp(-2*(V-50.)/RT_F)) * (Cai*exp(2*50./RT_F) - Cao*exp(-2*(V-50)/RT_F))*d*f)
  _(IsiK  ,  PCaK*PCa*(V-50.)/RT_F / (1-exp(-  (V-50.)/RT_F)) * (Ki *exp(  50./RT_F) - Ko *exp(-  (V-50)/RT_F))*d*f)
  _(IsiNa , PCaNa*PCa*(V-50.)/RT_F / (1-exp(-  (V-50.)/RT_F)) * (Nai*exp(  50./RT_F) - Nao*exp(-  (V-50)/RT_F))*d*f)
  _(INaK  , INaKx*Ko/(Ko+kmK) * Nai/(Nai+kmNa))
  _(INaCa , kNaCa*(exp(gamma*V/RT_F)*cub(Nai)*Cao-exp(-(1-gamma)*V/RT_F)*cub(Nao)*Cai)/(1.+dNaCa*(Cai*cub(Nao)+Cao*cub(Nai))))
  _(Iup   , (3.*Cai - 0.23*Caup*kcyca*kxcs/ksrca) / (Cai + Caup*kcyca*kxcs/ksrca + kcyca*kxcs + kcyca))
  _(Itr   , 50.*(Caup - Carel))
  _(Irel  , sqr(fact/(fact+0.25)) * kmCad * Carel) 
