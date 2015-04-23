/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
inline const Foam::discreteVelocity& 
Foam::fvDVM::DVi(const label I) const
{
    return DV_[I];
}

inline const Foam::discreteVelocity& Foam::fvDVM::DVxyz
(
    const label ix,
    const label iy,
    const label iz
 ) const
{
    const label I = iz*nXiY_*nXiX_ + iy*nXiX_ + ix;
    return DV_[I];
}

inline Foam::volScalarField& Foam::fvDVM::rhoVol()
{
    return rhoVol_;
}

inline const Foam::volVectorField& Foam::fvDVM::Uvol() const
{
    return Uvol_;
}

inline Foam::label Foam::fvDVM::nXi() const
{
    return nXi_;
}

inline Foam::label Foam::fvDVM::nXiPerDim() const
{
    return nXiPerDim_;
}

inline Foam::dimensionedScalar Foam::fvDVM::xiMax() const
{
    return xiMax_;
}

inline Foam::dimensionedScalar Foam::fvDVM::xiMin() const
{
    return xiMin_;
}

//inline Foam::scalar Foam::fvDVM::res() const
//{
    //return res_;
//}

//inline Foam::label Foam::fvDVM::checkSteps() const
//{
    //return checkSteps_;
//}

inline Foam::dimensionedScalar Foam::fvDVM::nu() const
{
    return nu_;
}

inline Foam::dimensionedScalar Foam::fvDVM::CsSqr() const
{
    return CsSqr_;
}

inline Foam::dimensionedScalar Foam::fvDVM::Cs() const
{
    return Cs_;
}

inline Foam::dimensionedScalar Foam::fvDVM::tau() const
{
    return tau_;
}

inline Foam::scalar Foam::fvDVM::omega() const
{
    return omega_;
}
