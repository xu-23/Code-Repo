#include "GenFvMatrix.H"
#include "multivariateGaussConvectionScheme.H"
#include "gaussConvectionScheme.H"
#include "snGradScheme.H"
#include "linear.H"
#include "orthogonalSnGrad.H"

#include "pimpleControl.H"
#include "fvOptions.H"
#include "turbulentFluidThermoModel.H"

#include "incompressibleTwoPhaseMixture.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"

#include "turbulenceModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

#include "gaussGrad.H"
#include "extrapolatedCalculatedFvPatchField.H"
#include <mpi.h>

namespace Foam{

tmp<fvScalarMatrix>
GenMatrix_TT(
    const volScalarField& rhoCp,
    volScalarField& TT,
    const surfaceScalarField& rhoPhiCP,   // phi  @dfLowMachFoam
    const volScalarField& kk
){
    Info << "TTEqn" << endl;
    
    const fvMesh& mesh = TT.mesh();

    tmp<fvScalarMatrix> tfvm_res
    (
        new fvScalarMatrix
        (
            TT,
            rhoCp.dimensions()*TT.dimensions()*dimVol/dimTime
        )
    );
    fvScalarMatrix& fvm = tfvm_res.ref();

    // -------------------------------------------------------





    // -------------------------------------------------------
    
    Info << "TTEqn ConvectionSchemeFvmDiv" << endl;
    word name("div(phi,T)");
    // ---------
    tmp<fv::convectionScheme<scalar>> cs = fv::convectionScheme<scalar>::New(mesh,rhoPhiCP,mesh.divScheme(name));
    fv::gaussConvectionScheme<scalar>& gcs = dynamic_cast<fv::gaussConvectionScheme<scalar>&>(cs.ref());
    tmp<surfaceScalarField> tweights = gcs.interpScheme().weights(TT);
    const surfaceScalarField& weights = tweights();
    // ---------
    fvm.lower() = -weights.primitiveField()*rhoPhiCP.primitiveField();
    fvm.upper() = fvm.lower() + rhoPhiCP.primitiveField();
    fvm.negSumDiag();
    forAll(TT.boundaryField(), patchi)
    {
        const fvPatchField<scalar>& psf = TT.boundaryField()[patchi];
        const fvsPatchScalarField& patchFlux = rhoPhiCP.boundaryField()[patchi];
        const fvsPatchScalarField& pw = weights.boundaryField()[patchi];

        fvm.internalCoeffs()[patchi] = patchFlux*psf.valueInternalCoeffs(pw);
        fvm.boundaryCoeffs()[patchi] = -patchFlux*psf.valueBoundaryCoeffs(pw);
    }
    if (gcs.interpScheme().corrected())
    {
        fvm += fvcSurfaceIntegrate(rhoPhiCP*gcs.interpScheme().correction(TT));
    }

    // -------------------------------------------------------
    
    Info << "TTEqn EulerDdtSchemeFvmDdt" << endl;

    scalar rDeltaT = 1.0/mesh.time().deltaTValue();

    fvm.diag() += rDeltaT*rhoCp.primitiveField()*mesh.Vsc();

    if (mesh.moving())
    {
        fvm.source() += rDeltaT
            *rhoCp.oldTime().primitiveField()
            *TT.oldTime().primitiveField()*mesh.Vsc0();
    }
    else
    {
        fvm.source() += rDeltaT
            *rhoCp.oldTime().primitiveField()
            *TT.oldTime().primitiveField()*mesh.Vsc();
    }

    // -------------------------------------------------------
    
    Info << "TEqn fvmSP_fvcDDT_fvcDIV" << endl;

    // ddt(rhoCp)
    tmp<volScalarField> tddtRhoCp
    (
        new volScalarField
        (
            IOobject
            (
                "ddtRhoCp",
                rhoCp.time().timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar(rhoCp.dimensions()/dimTime, 0.0)
        )
    );
    volScalarField& ddtRhoCp = tddtRhoCp.ref();
    
    if (mesh.moving())
    {
        ddtRhoCp.primitiveFieldRef() = rDeltaT *
            (
                rhoCp.primitiveField()
              - rhoCp.oldTime().primitiveField() * 
                mesh.Vsc0() / mesh.Vsc()
            );
        
        forAll(rhoCp.boundaryField(), patchi)
        {
            ddtRhoCp.boundaryFieldRef()[patchi] = rDeltaT *
                (
                    rhoCp.boundaryField()[patchi]
                  - rhoCp.oldTime().boundaryField()[patchi]
                );
        }
    }
    else
    {
        ddtRhoCp = rDeltaT * (rhoCp - rhoCp.oldTime());
    }


 


    // -------------------------------------------------------

    


    // -------------------------------------------------------




    // -------------------------------------------------------
    
    tmp<fvScalarMatrix> tfvm
    (
        fvm
        // fvm::ddt(rhoCp,TT)
        // + fvm::div(rhoPhiCP,TT, "div(phi,T)")
        // - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoPhiCP), TT)
        - fvm::Sp(ddtRhoCp + fvc::div(rhoPhiCP), TT) 
        - fvm::laplacian(kk, TT,  "laplacian(kk,T)")
    );

    return tfvm;
}
























































}