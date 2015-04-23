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

#include "fvDVM.H"
#include "discreteVelocity.H"
#include "constants.H"
#include "fixedGradientFvPatchField.H"

#include "farFieldFvPatchField.H"
#include "farFieldFvsPatchField.H"

#include "bounceBackFvPatchField.H"
#include "bounceBackFvsPatchField.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::discreteVelocity::discreteVelocity
(
    fvDVM& dvm,
    const fvMesh& mesh,
    const Time& time,
    const scalar weight,
    const dimensionedVector xi,
    const label DVid,
    const label symXtargetDVid,
    const label symYtargetDVid,
    const label symZtargetDVid
)
:
    dvm_(dvm),
    mesh_(mesh),
    time_(time),
    weight_(weight),
    xi_(xi),
    myDVid_(DVid),
    bDVid_(dvm_.nXi() - DVid - 1),
    symXtargetDVid_(symXtargetDVid),
    symYtargetDVid_(symYtargetDVid),
    symZtargetDVid_(symZtargetDVid),
    gVol_
    (
        IOobject
        (
            "gVol" + name(DVid),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "0", dimMass/pow3(dimLength), 0.0
        ),
        "fixedGradient"
    ),
    gSurf_
    (
        IOobject
        (
            "gSurf" + name(DVid),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "0", dimMass/pow3(dimLength), 0.0
        )
    ),
    gGrad_
    (
        IOobject
        (
            "gGrad" + name(DVid),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "0", gVol_.dimensions()/dimLength, vector(0,0,0)
        ),
        "zeroGradient"
    )
{
    initDFtoEq();
    setBCtype();
    initBoundaryField();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::discreteVelocity::~discreteVelocity()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::discreteVelocity::initDFtoEq()
{
    equilibrium
    (
        gVol_,
        dvm_.rhoVol(), 
        dvm_.Uvol() 
    );
}

void Foam::discreteVelocity::setBCtype()
{
    // Only from rho's B.C can we determine surfaceScalarField f's B.C.
    // for all patchi fo h/g barP, set to zeroGradient
    // May improve for inlet ingoing DF

    const GeometricField<vector, fvPatchField, volMesh>::GeometricBoundaryField& 
        uBCs = dvm_.Uvol().boundaryField();

    forAll(uBCs, patchi)
    {
        if (uBCs[patchi].type() == "farField") //inlet
        {
            gSurf_.boundaryField().set
            (
                patchi, 
                fvsPatchField<scalar>::New
                (
                    "farField", mesh_.boundary()[patchi], gSurf_
                )
            );
        }
        else if (uBCs[patchi].type() == "bounceBack") //maxwellWall
        {
            gSurf_.boundaryField().set
            (
                patchi, 
                fvsPatchField<scalar>::New
                (
                    "bounceBack", mesh_.boundary()[patchi], gSurf_
                )
            );
        }
    }
}

void Foam::discreteVelocity::initBoundaryField()
{
    GeometricField<scalar, fvsPatchField, surfaceMesh>::GeometricBoundaryField& 
        gBCs = gSurf_.boundaryField();

    forAll(gBCs, patchi)
    {
        if(gBCs[patchi].type() == "farField" )
        {
            equilibrium
            (
                gBCs[patchi],
                dvm_.rhoVol().boundaryField()[patchi],
                dvm_.Uvol().boundaryField()[patchi]
            );
        }
    }
}

void Foam::discreteVelocity::updateGvol()
{
    //- get delta t
    volScalarField gEq
    (
        IOobject
        (
            "gEq",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", gVol_.dimensions(), 0)
    );


    //- get gEq and hEq
    equilibrium
    (
        gEq,
        dvm_.rhoVol(), 
        dvm_.Uvol() 
    );
    gVol_ = (1.0 - dvm_.omega())*gVol_ + dvm_.omega()*gEq;
    gVol_.correctBoundaryConditions(); // NOTE: check if the newly defined zeroGradientFvsPatchField 
}

void Foam::discreteVelocity::updateGsurf()
{
  // Pre-setup
    // 1. correct the boundary value of gBarPvol_
    //    the gradient at boundary is known
    // 2. get the gradient
    //
    gGrad_ = fvc::grad(gVol_); 
    // The DVMsymmetry is  rocessed automatically in fvc::grad operator

    //
    // 3. correct the boundary value of the grad field
    //    to be used at next time
    gGrad_.correctBoundaryConditions();

    // 4. patch the normal component of boundary value of the grad 
    //    to the gradient field of the fixed 
    //    gradient feild of the gBarPvol_ and 
    //    hBarPvol_ ...
    //    NOTE: we need the surfaceNormal Gradient

    forAll(gGrad_.boundaryField(), patchi)
    {
        const vectorField n
        (
            mesh_.Sf().boundaryField()[patchi]
           /mesh_.magSf().boundaryField()[patchi]
        );

        if ( 
               gVol_.boundaryField()[patchi].type() != "empty" 
            && gVol_.boundaryField()[patchi].type() != "processor"
            && gVol_.boundaryField()[patchi].type() != "symmetryPlane"
           ) // only for fixed gradient g/hBarPvol
        {
            // normal component of the grad field
            fixedGradientFvPatchField<scalar>& gVolPatch = 
                refCast<fixedGradientFvPatchField<scalar> >
                (gVol_.boundaryField()[patchi]);

            forAll(gVolPatch, pFacei)
            {
                gVolPatch.gradient()[pFacei] =
                    gGrad_.boundaryField()[patchi][pFacei]&n[pFacei];
            }
        }
    }

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const vectorField& Cf = mesh_.Cf();
    const vectorField& Sf = mesh_.Sf();
    const vectorField& C = mesh_.C();

    const Field<scalar>& iGvol = gVol_;
    const Field<vector>& iGgrad = gGrad_;
    const Field<scalar>& ibGvol = dvm_.DVi(bDVid_).gVol();
    const Field<vector>& ibGgrad = dvm_.DVi(bDVid_).gGrad();

    // This is what we want to update in this function
    Field<scalar>& iGsurf = gSurf_;

    scalar dt = time_.deltaTValue();
    vector xii = xi_.value();

    // internal faces first
    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];
        if ((xii&Sf[facei]) >=  0) // comming from own
        {

            iGsurf[facei] = iGvol[own] 
              + (iGgrad[own]&(Cf[facei] - C[own] - 0.5*xii*dt));
        }
        // Debug, no = 0, =0 put to > 0
        else// if ((xii&Sf[facei]) < 0) // comming form nei
        {
            iGsurf[facei] = iGvol[nei]
              + (iGgrad[nei]&(Cf[facei] - C[nei] - 0.5*xii*dt));
        }
    }

  // boundary faces
    forAll(gSurf_.boundaryField(), patchi)
    {
        word type = gSurf_.boundaryField()[patchi].type();
        fvsPatchField<scalar>& gSurfPatch = gSurf_.boundaryField()[patchi];
        //const fvPatchField<vector>& Upatch = dvm_.Uvol().boundaryField()[patchi];
        const fvsPatchField<vector>& SfPatch =
            mesh_.Sf().boundaryField()[patchi];
        const fvsPatchField<vector>& CfPatch =
            mesh_.Cf().boundaryField()[patchi];
        const labelUList& faceCells = mesh_.boundary()[patchi].faceCells();
        
        //- NOTE: outging DF can be treate unifily for all BCs, including processor BC
        if (type == "zeroGradient")
        {
            gSurfPatch == gSurf_.boundaryField()[patchi];//.patchInternalField();
        }
        else if (type == "farField")
        {
            //check each boundary face in the patch
            forAll(gSurfPatch, facei)
            {
                //out or in ?
                if ((xii&SfPatch[facei]) > 0 ) // outgoing
                {
                    gSurfPatch[facei] = iGvol[faceCells[facei]] 
                      + ((iGgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - xii*dt));
                //incoming and parallel to face, not changed.
                }
            }
            Info << "farField's patch id: " << patchi << endl;
        }
        else if (type == "bounceBack")
        {
            //check each boundary face in the patch
            forAll(gSurfPatch, facei)
            {
                if ((xii&SfPatch[facei]) > 0 ) // outgoing
                {
                    gSurfPatch[facei] = iGvol[faceCells[facei]] 
                      + ((iGgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - xii*dt));
                //incoming and parallel to face, not changed.
                }
                //out or in ?
                else //((xii&SfPatch[facei]) < 0 ) // incomming
                {
                    gSurfPatch[facei] = ibGvol[faceCells[facei]] 
                      + ((ibGgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] + xii*dt));
                      //+ 2.0*weight_*(xii&Upatch[facei]);

                //incoming and parallel to face, not changed.
                }
            }
        }

        else if (type == "processor") // parallel
        {
            forAll(gSurfPatch, facei)
            {
                vector faceSf= SfPatch[facei];
                if ((xii&faceSf) >  0 ) // outgoing
                {
                    gSurfPatch[facei] = iGvol[faceCells[facei]] 
                      + ((iGgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - xii*dt));
                } 
                else //incomming from processor boundaryField
                {
                    gSurfPatch[facei] = gVol_.boundaryField()[patchi][facei]
                      + ((gGrad_.boundaryField()[patchi][facei])
                       &(CfPatch[facei] - mesh_.C().boundaryField()[patchi][facei] - xii*dt));
                }
            }
        }
    }
}

void Foam::discreteVelocity::updateGnewVol()
{
    // store the gTildePlus in gTilde
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const vectorField Sf = mesh_.Sf();
    const scalarField V = mesh_.V();
    const scalar dt = time_.deltaTValue();
    const vector xii = xi_.value();

    // internal faces
    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];
        gVol_[own] -= ((xii&Sf[facei])*gSurf_[facei]*dt/V[own]);
        gVol_[nei] += ((xii&Sf[facei])*gSurf_[facei]*dt/V[nei]);
    }
    // boundary faces
    forAll(gSurf_.boundaryField(), patchi)
    {
        const fvsPatchField<scalar>& gSurfPatch =
            gSurf_.boundaryField()[patchi];
        const fvsPatchField<vector>& SfPatch =
            mesh_.Sf().boundaryField()[patchi];
        const labelUList& pOwner = mesh_.boundary()[patchi].faceCells();
        forAll(gSurfPatch, pFacei)
        {
            const label own = pOwner[pFacei];
            gVol_[own] -= (xii&SfPatch[pFacei]) 
               *gSurfPatch[pFacei]*dt/V[own];
        }
    }

}

void Foam::discreteVelocity::equilibrium
(
    volScalarField& geq,
    const volScalarField&  rho,
    const volVectorField&  U
)
{
    dimensionedVector xii = xi_;//.value();
    dimensionedScalar CsSqr = dvm_.CsSqr();//.value();
    geq =  weight_*rho*(1.0 + (xii&U)/CsSqr + 0.5*((xii&U)*(xii&U))/CsSqr/CsSqr - 0.50*magSqr(U)/CsSqr);
    //geq =  weight_*rho*(1.0 + (xii&U)/CsSqr);// + 0.5*((xii&U)*(xii&U))/CsSqr/CsSqr);// - 0.50*magSqr(U)/CsSqr);
}

void Foam::discreteVelocity::equilibrium
(
    fvsPatchField<scalar>& geq,
    const fvPatchField<scalar>&  rho,
    const fvPatchField<vector>&  U
)
{
    vector xii = xi_.value();
    scalar CsSqr = dvm_.CsSqr().value();
    geq =  weight_*rho*(1.0 + (xii&U)/CsSqr + 0.5*((xii&U)*(xii&U))/CsSqr/CsSqr - 0.50*magSqr(U)/CsSqr);
}
// ************************************************************************* //
