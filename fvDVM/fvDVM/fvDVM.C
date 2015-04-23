/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "fvDVM.H"
#include "constants.H"
#include "fvm.H"
//#include "calculatedMaxwellFvPatchField.H"
#include "scalarIOList.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvDVM, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvDVM::setDVgrid
(
    scalarField& weights,
    scalarField& Xis,
    scalar xiMin,
    scalar xiMax,
    label nXi
)
{
    // Read from file ./constant/Xis and ./constant/weights
    Info << "Reading xi list" << endl;
    scalarIOList xiList
    (
        IOobject
        (
             "Xis",
             time_.time().constant(),
             mesh_,
             IOobject::MUST_READ,
             IOobject::NO_WRITE
        )
    );

    Info << "Reading weight list" << endl;
    scalarIOList weightList
    (
        IOobject
        (
             "weights",
             time_.time().constant(),
             mesh_,
             IOobject::MUST_READ,
             IOobject::NO_WRITE
        )
    );

    // Check
    if (
            xiList.size() != nXi 
         || weightList.size() != nXi
       )
    {
        WarningIn("setDVgrid")
            << "Num of discreteVelocity not consistent" << endl;
        std::exit(-1);
    }


    if ( 
           xiList[0]  != xiMin
        || xiList[nXi-1] != xiMax)
    {
        WarningIn("setDVgrid")
            << "xi not consistant" <<endl;
        std::exit(-1);
    }


    for (label i = 0; i < nXi ; i++)
    {
        weights[i] = weightList[i];
        Xis[i] = xiList[i];
    }
}

void Foam::fvDVM::initialiseDV()
{
    scalarField weights1D(nXiPerDim_);
    scalarField Xis(nXiPerDim_);

    //get discrete velocity points and weights
    setDVgrid
    (
         weights1D,
         Xis, 
         xiMin_.value(), 
         xiMax_.value(), 
         nXiPerDim_
    );

    if (mesh_.nSolutionD() == 3)    //3D(X & Y & Z)
    {
        nXiX_ = nXiY_ = nXiZ_ = nXiPerDim_;
        nXi_ = nXiX_*nXiY_*nXiZ_;
        DV_.setSize(nXi_);
        label i = 0;
        for (label iz = 0; iz < nXiZ_; iz++)
        {
            for (label iy = 0; iy < nXiY_; iy++)
            {
                for (label ix = 0; ix < nXiZ_; ix++)
                {
                    scalar weight = weights1D[iz]*weights1D[iy]*weights1D[ix];
                    dimensionedVector xi
                    (
                        "xi", 
                        dimLength/dimTime, 
                        vector(Xis[ix], Xis[iy], Xis[iz])
                    );

                    label symXtargetDVid = iz*nXiY_*nXiX_ + iy*nXiX_ + (nXiX_ - ix -1);
                    label symYtargetDVid = iz*nXiY_*nXiX_ + (nXiY_ - iy - 1)*nXiX_ + ix;
                    label symZtargetDVid = (nXiZ_ - iz -1)*nXiY_*nXiX_ + iy*nXiX_ + ix;
                    DV_.set
                    (
                        i,
                        new discreteVelocity
                        (
                            *this,
                            mesh_,
                            time_,
                            weight,
                            xi,
                            i,
                            symXtargetDVid,
                            symYtargetDVid,
                            symZtargetDVid
                        )
                    );
                    i++;
                }
            }
        }
    }
    else
    {
        if (mesh_.nSolutionD() == 2)    //2D (X & Y)
        {
            nXiX_ = nXiY_ = nXiPerDim_;
            nXiZ_ = 1;
            nXi_ = nXiX_*nXiY_*nXiZ_;
            DV_.setSize(nXi_);
            label i = 0;
            for (label iy = 0; iy < nXiY_; iy++)
            {
                for (label ix = 0; ix < nXiX_; ix++)
                {
                    scalar weight = weights1D[iy]*weights1D[ix]*1;
                    dimensionedVector xi
                    (
                        "xi",
                        dimLength/dimTime, 
                        vector(Xis[ix], Xis[iy], 0.0)
                    );
                    label symXtargetDVid = iy*nXiX_ + (nXiX_ - ix -1);
                    label symYtargetDVid = (nXiY_ - iy - 1)*nXiX_ + ix;
                    label symZtargetDVid = 0;
                    DV_.set
                    (
                        i,
                        new discreteVelocity
                        (
                            *this,
                            mesh_,
                            time_,
                            weight,
                            xi,
                            i,
                            symXtargetDVid,
                            symYtargetDVid,
                            symZtargetDVid
                        )
                    );
                    i++;
                }
            }
        }
        else    //1D (X)
        {
            nXiX_ = nXiPerDim_;
            nXiY_ = nXiZ_ = 1;
            nXi_ = nXiX_*nXiY_*nXiZ_;
            DV_.setSize(nXi_);
            label i = 0;
            for (label ix = 0; ix < nXiX_; ix++)
            {
                scalar weight = weights1D[ix]*1*1;
                dimensionedVector xi
                (
                    "xi",
                    dimLength/dimTime,
                    vector(Xis[ix], 0.0, 0.0)
                );
                label symXtargetDVid = (nXiX_ - ix -1);
                label symYtargetDVid = 0;
                label symZtargetDVid = 0;
                DV_.set
                (
                    i,
                    new discreteVelocity
                    (
                        *this,
                        mesh_,
                        time_,
                        weight,
                        xi,
                        i,
                        symXtargetDVid,
                        symYtargetDVid,
                        symZtargetDVid
                    )
                );
                i++;
            }
        }
    }

    Info<< "fvDVM : Allocated " << DV_.size()
        << " discrete velocities" << endl;
}

void Foam::fvDVM::updateGvol()
{
    forAll(DV_, DVid)
        DV_[DVid].updateGvol();
}

void Foam::fvDVM::updateGsurf()
{
    forAll(DV_, DVid)
        DV_[DVid].updateGsurf();
}

void Foam::fvDVM::updateGnewVol()
{
    forAll(DV_, DVid)
        DV_[DVid].updateGnewVol();
}

void Foam::fvDVM::updateMacroVol()
{
    volVectorField rhoUvol = rhoVol_*Uvol_;

    rhoVol_ =  dimensionedScalar("0", rhoVol_.dimensions(), 0);
    rhoUvol = dimensionedVector("0", rhoUvol.dimensions(), vector(0, 0, 0));

    forAll(DV_, dvi)
    {
        discreteVelocity& dv = DV_[dvi];
        rhoVol_ += dv.gVol();
        rhoUvol += dv.xi()*dv.gVol();
    }

    //- get Prim. from Consv.
    Uvol_ = rhoUvol/rhoVol_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvDVM::fvDVM
(
 volScalarField& rho,
 volVectorField& U
 )
    :
        IOdictionary
        (
         IOobject
         (
          "DVMProperties",
          U.time().constant(),
          U.mesh(),
          IOobject::MUST_READ,
          IOobject::NO_WRITE
         )
        ),
        mesh_(rho.mesh()),
        time_(rho.time()),
        rhoVol_(rho),
        Uvol_(U),
        fvDVMparas_(subOrEmptyDict("fvDVMparas")),
        gasProperties_(subOrEmptyDict("gasProperties")),
        nXiPerDim_(readLabel(fvDVMparas_.lookup("nDV"))),
        xiMax_(fvDVMparas_.lookup("xiMax")),
        xiMin_(fvDVMparas_.lookup("xiMin")),
        nu_(gasProperties_.lookup("nu")),
        CsSqr_(dimensionedScalar("CsSqr",pow(dimLength/dimTime,2), 1.0/3.0)),
        Cs_(dimensionedScalar("Cs",dimLength/dimTime, 1.0/sqrt(3.0))),
        tau_(nu_/CsSqr_),
        omega_(scalar(1.0/(tau_.value()/time_.deltaTValue()+0.5))),
        DV_(0)
{
    initialiseDV();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvDVM::~fvDVM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvDVM::evolution()
{
    //DEBUG
    const labelUList& faceCells = mesh_.boundary()[1].faceCells();
    const fvsPatchField<vector>& CfPatch = mesh_.Cf().boundaryField()[1];
    updateGvol();
    updateGsurf();
    //Check gSurf at farField
    forAll(DV_, dvi)
    {
        Info << "Dv[" << dvi <<"] " << DV_[dvi].gSurf().boundaryField()[1][10] << endl;
    }
    Info << "Probe point x, y, z" << CfPatch[faceCells[10]] << endl;
    updateGnewVol();
    updateMacroVol();
}

// ************************************************************************* //
