#!/usr/bin/env python

#--------------------------------------------------------------------------------------
## pythonFlu - Python wrapping for OpenFOAM C++ API
## Copyright (C) 2010- Alexey Petrov
## Copyright (C) 2009-2010 Pebble Bed Modular Reactor (Pty) Limited (PBMR)

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## See http://sourceforge.net/projects/pythonflu
##
## Author : Alexey Petrov, Andrey Simurzin
##

#---------------------------------------------------------------------------
from Foam import ref, man


#---------------------------------------------------------------------------
def createFields( runTime, mesh, pimple ):
    ref.ext_Info()<< "Reading thermophysical properties\n" << ref.nl
    
    pThermo = man.basicPsiThermo.New( mesh )
    
    p = man.volScalarField( pThermo.p(), man.Deps( pThermo ) )
    h = man.volScalarField( pThermo.h(), man.Deps( pThermo ) )
    psi = man.volScalarField( pThermo.psi(), man.Deps( pThermo ) )
    
    rho = man.volScalarField( man.IOobject( ref.word( "rho" ),
                                            ref.fileName( runTime.timeName() ),
                                            mesh,
                                            ref.IOobject.READ_IF_PRESENT,
                                            ref.IOobject.AUTO_WRITE ),
                              man.volScalarField( pThermo.rho(), man.Deps( pThermo ) ) )
    
    ref.ext_Info()<< "Reading field U\n" << ref.nl
    U = man.volVectorField( man.IOobject( ref.word( "U" ),
                                          ref.fileName( runTime.timeName() ),
                                          mesh,
                                          ref.IOobject.MUST_READ,
                                          ref.IOobject.AUTO_WRITE ),
                            mesh )

    phi = man.compressibleCreatePhi( runTime, mesh, U, rho )

    rhoMax = ref.dimensionedScalar( pimple.dict().lookup( ref.word( "rhoMax" ) ) )
    rhoMin = ref.dimensionedScalar( pimple.dict().lookup( ref.word( "rhoMin" ) ) )

    ref.ext_Info()<< "Creating turbulence model\n" << ref.nl
    turbulence = man.compressible.turbulenceModel.New( rho, U, phi, pThermo );
  
    ref.ext_Info()<< "Creating field dpdt\n" << ref.nl
    dpdt = man.volScalarField( ref.word( "dpdt" ), man.fvc.ddt(p) )
    
    ref.ext_Info() << "Creating field kinetic energy K\n" << ref.nl
    K = man.volScalarField( ref.word( "K" ), man.volScalarField( 0.5 * U.magSqr(), man.Deps( U ) ) )
  
    return pThermo, p, h, psi, rho, U, phi, rhoMax, rhoMin, turbulence, dpdt, K


#---------------------------------------------------------------------------
def fun_Ueqn( pimple, rho, p, U, phi, K, turbulence ):
    # The initial C++ expression does not work properly, because of
    #  1. turbulence.divDevRhoReff( U ) - changes values for the U boundaries
    #  2. the order of expression arguments computation differs with C++
    # UEqn = man.fvm.ddt( rho, U ) + man.fvm.div( phi, U ) + man.fvVectorMatrix( turbulence.divDevRhoReff( U ), man.Deps( turbulence, U ) )
    
    UEqn = man.fvVectorMatrix( turbulence.divDevRhoReff( U ), man.Deps( turbulence, U ) ) + ( man.fvm.div( phi, U ) + man.fvm.ddt( rho, U ) )
  
    UEqn.relax()
  
    rAU = 1.0 / UEqn.A()
  
    if pimple.momentumPredictor():
        ref.solve( UEqn == -man.fvc.grad( p ) )
        K << 0.5 * U.magSqr()
        pass
    
    return UEqn, rAU


#---------------------------------------------------------------------------
def fun_hEqn( thermo, rho, p, h, phi, turbulence, dpdt, K ):
    hEqn = ( ref.fvm.ddt( rho, h ) + ref.fvm.div( phi, h ) - ref.fvm.laplacian( turbulence.alphaEff(), h ) 
             == dpdt() - ( ref.fvc.ddt( rho, K ) + ref.fvc.div( phi, K ) ) ) # mixed calculations

    hEqn.relax()
    hEqn.solve()

    thermo.correct()
    pass


#---------------------------------------------------------------------------
def fun_pEqn( mesh, runTime, pimple, thermo, rho, p, h, psi, U, phi, turbulence, UEqn, rAU, dpdt, K, cumulativeContErr, rhoMax, rhoMin ):

    rho << thermo.rho()
    rho << rho().ext_max( rhoMin )
    rho << rho().ext_min( rhoMax )
    rho.relax()

    U << rAU() * UEqn.H()

    if pimple.transonic():
        phid = ref.surfaceScalarField( ref.word( "phid" ), ref.fvc.interpolate( psi ) * ( ( ref.fvc.interpolate( U ) & mesh.Sf() ) \
                                                                                          + ref.fvc.ddtPhiCorr( rAU, rho, U, phi ) ) );
        while pimple.correctNonOrthogonal():
            pEqn = ref.fvm.ddt( psi, p) + ref.fvm.div( phid, p ) - ref.fvm.laplacian( rho() * rAU, p ) # mixed calculations
            
            pEqn.solve( mesh.solver( p.select( pimple.finalInnerIter() ) ) )

            if pimple.finalNonOrthogonalIter():
                phi == pEqn.flux()
                pass
            pass
        pass
    else:
        phi << ref.fvc.interpolate( rho ) * ( ( ref.fvc.interpolate( U ) & mesh.Sf() ) + ref.fvc.ddtPhiCorr( rAU, rho, U, phi ) )
        
        while pimple.correctNonOrthogonal():
            pEqn = ref.fvm.ddt( psi, p ) + ref.fvc.div( phi ) - ref.fvm.laplacian( rho()*rAU, p )  # mixed calculations
            
            pEqn.solve( mesh.solver( p.select( pimple.finalInnerIter() ) ) )
            
            if pimple.finalNonOrthogonalIter():
                phi += pEqn.flux()
                pass
            pass    
        pass

    ref.rhoEqn( rho, phi )
    cumulativeContErr = ref.compressibleContinuityErrs( rho(), thermo, cumulativeContErr ) # mixed calculations
    # Explicitly relax pressure for momentum corrector
    p.relax()
    
    # Recalculate density from the relaxed pressure
    rho << thermo.rho()
    rho << rho().ext_max( rhoMin )
    rho << rho().ext_min( rhoMax )
    rho.relax()
    ref.ext_Info()<< "rho max/min : " << rho.ext_max().value()  << " " << rho.ext_min().value() << ref.nl

    U -= rAU * ref.fvc.grad( p )
    U.correctBoundaryConditions()

    K << 0.5 * U.magSqr()
    dpdt << ref.fvc.ddt( p )
    
    return cumulativeContErr



#---------------------------------------------------------------------------
def main_standalone( argc, argv ):
 
    args = ref.setRootCase( argc, argv )
    
    runTime = man.createTime( args )

    mesh = man.createMesh( runTime )
    
    pimple = man.pimpleControl(mesh)
    
    pThermo, p, h, psi, rho, U, phi, rhoMax, rhoMin, turbulence, dpdt, K = createFields( runTime, mesh, pimple )
  
    cumulativeContErr = ref.initContinuityErrs()

    ref.ext_Info()<< "\nStarting time loop\n" << ref.nl;

    while runTime.run():
    
        adjustTimeStep, maxCo, maxDeltaT = ref.readTimeControls( runTime )
        CoNum, meanCoNum = ref.compressibleCourantNo( mesh, phi, rho, runTime )

        runTime = ref.setDeltaT( runTime, adjustTimeStep, maxCo, maxDeltaT, CoNum )
        runTime.increment()

        ref.ext_Info()<< "Time = " << runTime.timeName() << ref.nl << ref.nl;

        ref.rhoEqn( rho, phi )

        # --- Pressure-velocity PIMPLE corrector loop
        while pimple.loop():

            UEqn, rAU = fun_Ueqn( pimple, rho, p, U, phi, K, turbulence )
      
            fun_hEqn( pThermo, rho, p, h, phi, turbulence, dpdt, K )

            # --- PISO loop
            while pimple.correct():
                cumulativeContErr = fun_pEqn( mesh, runTime, pimple, pThermo, rho, p, h, psi, U, phi, 
                                              turbulence, UEqn, rAU, dpdt, K, cumulativeContErr, rhoMax, rhoMin )
                pass
            
            if pimple.turbCorr():
                turbulence.correct()
                pass
            
            pass
        
        runTime.write()
            
        ref.ext_Info()<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" \
                      << "  ClockTime = " << runTime.elapsedClockTime() << " s" \
                      << ref.nl << ref.nl
        pass

    ref.ext_Info()<< "End\n" << ref.nl
        
    import os
    return os.EX_OK


#--------------------------------------------------------------------------------------
import sys, os
from Foam import FOAM_REF_VERSION
if FOAM_REF_VERSION( ">=", "020100" ):
   if __name__ == "__main__" :
      argv = sys.argv
      os._exit( main_standalone( len( argv ), argv ) )
      pass
   pass
else:
   from Foam.OpenFOAM import ext_Info
   ext_Info()<< "\nTo use this solver, It is necessary to SWIG OpenFoam2.1.0 or higher \n "


#--------------------------------------------------------------------------------------
