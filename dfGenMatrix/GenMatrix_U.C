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

double UEqn_build_pre_time = 0., UEqn_build_item1_intern_time = 0., UEqn_build_item1_bound_time = 0., UEqn_build_item2_intern_time = 0., UEqn_build_item2_bound_time = 0., UEqn_build_item3_time = 0., UEqn_build_item3_last_time = 0;
double startU, endU;

namespace Foam{

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    >
>
gaussGradCalcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
);

template<class Type>
void gaussGradCorrectBoundaryConditions
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >& gGrad
);

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    >
>
gaussGradGradf
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf,
    const word& name
);

// // 专门用于计算 grad(U) 的函数
// tmp<volTensorField> calcGradU(const volVectorField& U)
// {
//     const fvMesh& mesh = U.mesh();
    
//     // 创建线性插值方案
//     tmp<surfaceInterpolationScheme<vector>> tinterpScheme =
//         tmp<surfaceInterpolationScheme<vector>>
//         (
//             new linear<vector>(mesh)
//         );

//     // 插值到场表面
//     tmp<surfaceVectorField> tinterpolate = tinterpScheme().interpolate(U);

//     // 计算梯度
//     tmp<volTensorField> tgGrad = gaussGradGradf(tinterpolate(), "grad(" + U.name() + ')');
//     volTensorField& gGrad = tgGrad.ref();

//     // 修正边界条件
//     gaussGradCorrectBoundaryConditions(U, gGrad);

//     return tgGrad;
// }

// 梯度计算核心函数实现
template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    >
>
gaussGradCalcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
)
{
    const fvMesh& mesh = vsf.mesh();

    tmp<surfaceInterpolationScheme<Type>> tinterpScheme_ =
    tmp<surfaceInterpolationScheme<Type>>
    (
        new linear<Type>(mesh)
    );

    typedef typename outerProduct<vector, Type>::type GradType;

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tinterpolate = tinterpScheme_().interpolate(vsf);

    tmp<GeometricField<GradType, fvPatchField, volMesh>> tgGrad
    (
        gaussGradGradf(tinterpolate.ref(), name)
    );
    GeometricField<GradType, fvPatchField, volMesh>& gGrad = tgGrad.ref();

    gaussGradCorrectBoundaryConditions(vsf, gGrad);

    return tgGrad;
}

template<class Type>
void gaussGradCorrectBoundaryConditions
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >& gGrad
)
{
    typename GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >::Boundary& gGradbf = gGrad.boundaryFieldRef();

    forAll(vsf.boundaryField(), patchi)
    {
        if (!vsf.boundaryField()[patchi].coupled())
        {
            const vectorField n
            (
                vsf.mesh().Sf().boundaryField()[patchi]
              / vsf.mesh().magSf().boundaryField()[patchi]
            );

            gGradbf[patchi] += n *
            (
                vsf.boundaryField()[patchi].snGrad()
              - (n & gGradbf[patchi])
            );
        }
     }
}

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    >
>
gaussGradGradf
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf,
    const word& name
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = ssf.mesh();

    tmp<GeometricField<GradType, fvPatchField, volMesh>> tgGrad
    (
        new GeometricField<GradType, fvPatchField, volMesh>
        (
            IOobject
            (
                name,
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<GradType>
            (
                "0",
                ssf.dimensions()/dimLength,
                Zero
            ),
            extrapolatedCalculatedFvPatchField<GradType>::typeName
        )
    );
    GeometricField<GradType, fvPatchField, volMesh>& gGrad = tgGrad.ref();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const vectorField& Sf = mesh.Sf();

    Field<GradType>& igGrad = gGrad;
    const Field<Type>& issf = ssf;

    forAll(owner, facei)
    {
        GradType Sfssf = Sf[facei]*issf[facei];

        igGrad[owner[facei]] += Sfssf;
        igGrad[neighbour[facei]] -= Sfssf;
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const vectorField& pSf = mesh.Sf().boundaryField()[patchi];

        const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            igGrad[pFaceCells[facei]] += pSf[facei]*pssf[facei];
        }
    }

    igGrad /= mesh.V();

    gGrad.correctBoundaryConditions();

    return tgGrad;
}

tmp<fvVectorMatrix>
GenMatrix_U(
    const volScalarField& rho,
    volVectorField& U,
    const surfaceScalarField& rhoPhi,   // phi  @dfLowMachFoam
    const volScalarField& p, 
    incompressible::turbulenceModel& turbulence
){
    UEqn_build_pre_time = 0., UEqn_build_item1_intern_time = 0., UEqn_build_item1_bound_time = 0., UEqn_build_item2_intern_time = 0., UEqn_build_item2_bound_time = 0., UEqn_build_item3_time = 0., UEqn_build_item3_last_time = 0;
    startU = MPI_Wtime();

    word name("div("+rhoPhi.name()+','+U.name()+')');

    const fvMesh& mesh = U.mesh();
    // div
    tmp<fv::convectionScheme<vector>> cs = fv::convectionScheme<vector>::New(mesh,rhoPhi,mesh.divScheme(name));
    fv::gaussConvectionScheme<vector>& gcs = dynamic_cast<fv::gaussConvectionScheme<vector>&>(cs.ref());
    tmp<surfaceScalarField> tweights = gcs.interpScheme().weights(U);
    const surfaceScalarField& weights = tweights();

    // -------------------------------------------------------

    Info << "UEqn EulerDdtSchemeFvmDdt" << endl;
    Info << "UEqn gaussConvectionSchemeFvmDiv" << endl;

    tmp<fvVectorMatrix> tfvm_DDT
    (
        new fvVectorMatrix
        (
            U,
            rho.dimensions()*U.dimensions()*dimVol/dimTime
        )
    );
    fvVectorMatrix& fvm = tfvm_DDT.ref();

    // -------------
    scalar rDeltaT = 1.0/mesh.time().deltaTValue();

    Info << "UEqn gaussLaplacianSchemeFvmLaplacian" << endl;
    
    const volScalarField alphaEff = turbulence.alpha()*rho*turbulence.nuEff();
    // tmp<fv::snGradScheme<vector>> tsnGradScheme_(new fv::orthogonalSnGrad<vector>(mesh));  // 正交修正方案
    tmp<fv::snGradScheme<vector>> tsnGradScheme_(new fv::correctedSnGrad<vector>(mesh));  // 对应snGradSchemes{default         corrected;}

    surfaceScalarField alphaEfff = fvc::interpolate(alphaEff);    ///  
    
    // gammaMagSf = alphaEfff * magSf
    GeometricField<scalar, fvsPatchField, surfaceMesh> gammaMagSf
    (
        alphaEfff * mesh.magSf()
    );
    
    const surfaceScalarField& deltaCoeffs = tsnGradScheme_().deltaCoeffs(U);
    // -------------

    endU = MPI_Wtime();
    UEqn_build_pre_time += endU - startU;

    // interField
    // fvm.lower() = -weights.primitiveField()*rhoPhi.primitiveField();
    // // fvm.upper() = fvm.lower() + rhoPhi.primitiveField();
    // fvm.upper() = -weights.primitiveField()*rhoPhi.primitiveField() + rhoPhi.primitiveField();
    // fvm.negSumDiag();   
    // fvm.diag() += rDeltaT*rho.primitiveField()*mesh.Vsc(); 
    // fvm.source() = rDeltaT
    //     *rho.oldTime().primitiveField()
    //     *U.oldTime().primitiveField()*mesh.Vsc();


    // interField
    // fvm_laplace.upper() = deltaCoeffs.primitiveField()*gammaMagSf.primitiveField();
    // fvm_laplace.negSumDiag();
    // fvm -=fvm_laplace;

    startU = MPI_Wtime();

    vector* __restrict__ sourcePtr = fvm.source().begin();
    scalar* __restrict__ diagPtr = fvm.diag().begin();  
    scalar* __restrict__ lowerPtr = fvm.lower().begin();
    scalar* __restrict__ upperPtr = fvm.upper().begin();

    const labelUList& l = fvm.lduAddr().lowerAddr();
    const labelUList& u = fvm.lduAddr().upperAddr();

    const scalar* const __restrict__ gammaMagSfPtr = gammaMagSf.primitiveField().begin();
    const scalar* const __restrict__ deltaCoeffsPtr = deltaCoeffs.primitiveField().begin();

    const scalar* const __restrict__ weightsPtr = weights.primitiveField().begin();
    const scalar* const __restrict__ rhoPhiPtr = rhoPhi.primitiveField().begin();
    const scalar* const __restrict__ rhoPtr = rho.primitiveField().begin();
    const scalar* const __restrict__ meshVscPtr = mesh.Vsc()().begin();
    const scalar* const __restrict__ rhoOldTimePtr = rho.oldTime().primitiveField().begin();

    const vector* const __restrict__ UOldTimePtr = U.oldTime().primitiveField().begin();

    const label nFaces = fvm.lower().size();
    const label nCells = fvm.diag().size();

    for (label celli = 0; celli < nCells; celli++)
    {
        diagPtr[celli] = 0.0;
    }

    for (label facei = 0; facei < nFaces; facei++)
    {
        scalar flux = weightsPtr[facei] * rhoPhiPtr[facei];
        lowerPtr[facei] = -flux;
        upperPtr[facei] = -flux + rhoPhiPtr[facei];
        
        diagPtr[l[facei]] += flux;         
        diagPtr[u[facei]] += flux - rhoPhiPtr[facei];  
    }

    for (label celli = 0; celli < nCells; celli++)
    {
        diagPtr[celli] += rDeltaT * rhoPtr[celli] * meshVscPtr[celli];
        sourcePtr[celli] = rDeltaT * rhoOldTimePtr[celli] * UOldTimePtr[celli] * meshVscPtr[celli];
    }

    for (label facei = 0; facei < nFaces; facei++)
    {
        scalar coeffs = deltaCoeffsPtr[facei] * gammaMagSfPtr[facei];

        upperPtr[facei] -= coeffs;
        lowerPtr[facei] -= coeffs;
        
        diagPtr[l[facei]] += coeffs; 
        diagPtr[u[facei]] += coeffs;  
    }

    endU = MPI_Wtime();
    UEqn_build_item1_intern_time += endU - startU;

    startU = MPI_Wtime();

    // boundaryField
    // forAll(U.boundaryField(), patchi)
    // {
    //     const fvPatchField<vector>& psf = U.boundaryField()[patchi];
    //     const fvsPatchScalarField& patchFlux = rhoPhi.boundaryField()[patchi];
    //     const fvsPatchScalarField& pw = weights.boundaryField()[patchi];

    //     fvm.internalCoeffs()[patchi] = patchFlux*psf.valueInternalCoeffs(pw);
    //     fvm.boundaryCoeffs()[patchi] = -patchFlux*psf.valueBoundaryCoeffs(pw);
    // }

    const label nPatches = U.boundaryField().size();
    for (label patchi = 0; patchi < nPatches; patchi++)
    {
        const fvPatchField<vector>& psf = U.boundaryField()[patchi];
        const fvsPatchScalarField& patchFlux = rhoPhi.boundaryField()[patchi];
        const fvsPatchScalarField& pw = weights.boundaryField()[patchi];
        
        // 获取边界系数数组指针
        Field<vector>* internalCoeffsPtr = &fvm.internalCoeffs()[patchi];
        Field<vector>* boundaryCoeffsPtr = &fvm.boundaryCoeffs()[patchi];
        
        // 计算边界系数
        Field<vector> valueInternalCoeffs = psf.valueInternalCoeffs(pw);
        Field<vector> valueBoundaryCoeffs = psf.valueBoundaryCoeffs(pw);
        
        const label patchSize = patchFlux.size();
        const scalar* const __restrict__ patchFluxPtr = patchFlux.begin();
        const vector* const __restrict__ valueInternalCoeffsPtr = valueInternalCoeffs.begin();
        const vector* const __restrict__ valueBoundaryCoeffsPtr = valueBoundaryCoeffs.begin();
        vector* const __restrict__ internalCoeffsPatchPtr = internalCoeffsPtr->begin();
        vector* const __restrict__ boundaryCoeffsPatchPtr = boundaryCoeffsPtr->begin();
        
        for (label i = 0; i < patchSize; i++)
        {
            internalCoeffsPatchPtr[i] = patchFluxPtr[i] * valueInternalCoeffsPtr[i];
            boundaryCoeffsPatchPtr[i] = -patchFluxPtr[i] * valueBoundaryCoeffsPtr[i];
        }
    }

    // correct
    if (gcs.interpScheme().corrected())
    {
        fvm += fvc::surfaceIntegrate(rhoPhi*gcs.interpScheme().correction(U));
    }

    endU = MPI_Wtime();
    UEqn_build_item1_bound_time += endU - startU;
    // --------------------------------------------------------------------------------------------------------------
    startU = MPI_Wtime();

    // forAll(U.boundaryField(), patchi)
    // {
    //     const fvPatchVectorField& pvf = U.boundaryField()[patchi];
    //     const fvsPatchScalarField& pGamma = gammaMagSf.boundaryField()[patchi];
    //     const fvsPatchScalarField& pDeltaCoeffs = deltaCoeffs.boundaryField()[patchi];
        
    //     if (pvf.coupled())
    //     {
    //         fvm.internalCoeffs()[patchi] -=
    //             pGamma * pvf.gradientInternalCoeffs(pDeltaCoeffs);
    //         fvm.boundaryCoeffs()[patchi] -=
    //            -pGamma * pvf.gradientBoundaryCoeffs(pDeltaCoeffs);
    //     }
    //     else
    //     {
    //         fvm.internalCoeffs()[patchi] -= pGamma * pvf.gradientInternalCoeffs();
    //         fvm.boundaryCoeffs()[patchi] -= -pGamma * pvf.gradientBoundaryCoeffs();
    //     }
    // }

    for (label patchi = 0; patchi < nPatches; patchi++)
    {
        const fvPatchVectorField& pvf = U.boundaryField()[patchi];
        const fvsPatchScalarField& pGamma = gammaMagSf.boundaryField()[patchi];
        const fvsPatchScalarField& pDeltaCoeffs = deltaCoeffs.boundaryField()[patchi];
        
        // 获取边界系数数组指针
        Field<vector>* internalCoeffsPtr = &fvm.internalCoeffs()[patchi];
        Field<vector>* boundaryCoeffsPtr = &fvm.boundaryCoeffs()[patchi];
        
        const label patchSize = pGamma.size();
        const scalar* const __restrict__ pGammaPtr = pGamma.begin();
        const scalar* const __restrict__ pDeltaCoeffsPtr = pDeltaCoeffs.begin();
        vector* const __restrict__ internalCoeffsPatchPtr = internalCoeffsPtr->begin();
        vector* const __restrict__ boundaryCoeffsPatchPtr = boundaryCoeffsPtr->begin();
        
        if (pvf.coupled())
        {
            Field<vector> gradientInternalCoeffs = pvf.gradientInternalCoeffs(pDeltaCoeffs);
            Field<vector> gradientBoundaryCoeffs = pvf.gradientBoundaryCoeffs(pDeltaCoeffs);
            
            const vector* const __restrict__ gradientInternalCoeffsPtr = gradientInternalCoeffs.begin();
            const vector* const __restrict__ gradientBoundaryCoeffsPtr = gradientBoundaryCoeffs.begin();
            
            for (label i = 0; i < patchSize; i++)
            {
                internalCoeffsPatchPtr[i] -= pGammaPtr[i] * gradientInternalCoeffsPtr[i];
                boundaryCoeffsPatchPtr[i] -= -pGammaPtr[i] * gradientBoundaryCoeffsPtr[i];
            }
        }
        else
        {
            Field<vector> gradientInternalCoeffs = pvf.gradientInternalCoeffs();
            Field<vector> gradientBoundaryCoeffs = pvf.gradientBoundaryCoeffs();
            
            const vector* const __restrict__ gradientInternalCoeffsPtr = gradientInternalCoeffs.begin();
            const vector* const __restrict__ gradientBoundaryCoeffsPtr = gradientBoundaryCoeffs.begin();
            
            for (label i = 0; i < patchSize; i++)
            {
                internalCoeffsPatchPtr[i] -= pGammaPtr[i] * gradientInternalCoeffsPtr[i];
                boundaryCoeffsPatchPtr[i] -= -pGammaPtr[i] * gradientBoundaryCoeffsPtr[i];
            }
        }
    }

    if (mesh.fluxRequired(U.name()))
    {
        tmp<GeometricField<vector, fvsPatchField, surfaceMesh>> tCorrection
        (
            gammaMagSf * tsnGradScheme_().correction(U)
        );
        
        tmp<GeometricField<vector, fvPatchField, volMesh>> tDivCorrection
        (
            fvc::div(tCorrection())
        );
        
        const vector* const __restrict__ divCorrectionPtr = 
            tDivCorrection().primitiveField().begin();
        
        for (label celli = 0; celli < nCells; celli++)
        {
            sourcePtr[celli] += meshVscPtr[celli] * divCorrectionPtr[celli];
        }
        
        fvm.faceFluxCorrectionPtr() = new
        GeometricField<vector, fvsPatchField, surfaceMesh>
        (
            -tCorrection
        );
    }
    else
    {
        tmp<GeometricField<vector, fvsPatchField, surfaceMesh>> tCorrection
        (
            gammaMagSf * tsnGradScheme_().correction(U)
        );
        
        tmp<GeometricField<vector, fvPatchField, volMesh>> tDivCorrection
        (
            fvc::div(tCorrection())
        );
        
        const vector* const __restrict__ divCorrectionPtr = 
            tDivCorrection().primitiveField().begin();
        
        for (label celli = 0; celli < nCells; celli++)
        {
            sourcePtr[celli] += meshVscPtr[celli] * divCorrectionPtr[celli];
        }
    }

    endU = MPI_Wtime();
    UEqn_build_item2_bound_time += endU - startU;

    // --------------------------------------------------------------------------------------------------------------
    startU = MPI_Wtime();

    Info << "UEqn Calculating grad(U) using extracted functions" << endl;
    tmp<surfaceInterpolationScheme<vector>> tinterpScheme =
    tmp<surfaceInterpolationScheme<vector>>
    (
        new linear<vector>(mesh)
    );
    tmp<surfaceVectorField> tinterpolate = tinterpScheme().interpolate(U);
    tmp<volTensorField> tgGradU = gaussGradGradf(tinterpolate(), "grad(" + U.name() + ')');
    volTensorField& gGrad = tgGradU.ref();
    gaussGradCorrectBoundaryConditions(U, gGrad);

    const volTensorField& gradU = tgGradU();
    // const volScalarField alphaEff = turbulence.alpha()*rho*turbulence.nuEff();

    endU = MPI_Wtime();
    UEqn_build_item3_time += endU - startU;

    startU = MPI_Wtime();

    // -------------------------------------------------------

    // -------------------------------------------------------  
    // interFoam    
    // const alphaField& alpha;

    tmp<fvVectorMatrix> tfvm
    (
            tfvm_DDT
        // + fvm::ddt(rho, U) 
        // + fvm::div(rhoPhi, U)
        // + MRF.DDt(rho, U) 
        // + turbulence.divDevRhoReff(rho, U)  
        // - fvc::div((turbulence.alpha()*rho*turbulence.nuEff())*dev2(T(fvc::grad(U))))
        - fvc::div((turbulence.alpha()*rho*turbulence.nuEff())*dev2(T(gradU)))
        // - fvc::div((turbulence.alpha()*rho*turbulence.nuEff())*dev2T_gradU)
        
        // - fvm::laplacian(turbulence.alpha()*rho*turbulence.nuEff(), U)
        // - fvOptions(rho, U)
    );

    // @dfLowMachFoam    
    // tmp<fvVectorMatrix> tUEqn
    // (
    //     fvm::ddt(rho, U) + fvm::div(phi, U)
    //     + turbulence->divDevRhoReff(U)
    //     == 
    //     -fvc::grad(p)
    // );     

    endU = MPI_Wtime();
    UEqn_build_item3_last_time += endU - startU;

    Info<< "========Time Spent in UEqn build once=================="<< endl;
    Info << "U pre Time :   " << UEqn_build_pre_time << endl;
    Info << "U item1 intern Time : " << UEqn_build_item1_intern_time << endl;
    Info << "U item1 bound Time : " << UEqn_build_item1_bound_time << endl;
    // Info << "U item2 intern Time : " << UEqn_build_item2_intern_time << endl;
    Info << "U item2 bound Time : " << UEqn_build_item2_bound_time << endl;
    Info << "U item3 Time : " << UEqn_build_item3_time << endl;
    Info << "U item3 last Time : " << UEqn_build_item3_last_time << endl;
    Info<< "============================================"<<nl<< endl;

    // return tfvm_DDT;
    return tfvm;
}

}

// tmp<fvVectorMatrix> tUEqn = GenMatrix_U(rho, U, phi, p, turbulence());
// fvVectorMatrix& UEqn = tUEqn.ref();

// time_monitor_UEqn_pre += UEqnClock.timeIncrement();

// tmp<fvVectorMatrix> tUEqn_answer
// (
//         fvm::ddt(rho, U) 
//     + fvm::div(rhoPhi, U)
//     + MRF.DDt(rho, U)
//     + turbulence->divDevRhoReff(rho, U)
//     ==
//         fvOptions(rho, U)
// );
// fvVectorMatrix& UEqn_answer = tUEqn_answer.ref();
// check_fvmatrix_equal(UEqn_answer, UEqn, "UEqn");




