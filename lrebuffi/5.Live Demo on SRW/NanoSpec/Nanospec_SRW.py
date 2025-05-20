# -*- coding: utf-8 -*-
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
from srwlib import *
import sys, numpy
import time


#**********************Input Parameters:
strIntSE_OutFileName     = 'nanospec_res_int_se.dat' #file name for output initial single-electron SR intensity data
strIntPropSE_OutFileName = 'nanospec_res_int_prop_se.dat' #file name for output propagated single-electron SR intensity data
strIntPropME_OutFileName = 'nanospec_res_int_prop_me.dat' #file name for output propagated multi-electron SR intensity data

numPer = 40           # Number of ID Periods (without counting for terminations
undPer = 0.1          # Period Length [m]


class Orientation:
    UP = 0
    DOWN = 1
    LEFT = 2
    RIGHT = 3

import scipy.constants as codata
def magnetic_field_from_K(K, period_length):
    return K * 2 * pi * codata.m_e * codata.c / (codata.e * period_length)


def createUndulator(Kv, undPer, numPer):
    #***********Undulator
    By = magnetic_field_from_K(Kv, undPer) #Peak Vertical field [T]
    print("By calculated: " + str(By) + " T")
    Bx = 0.0 #Peak Vertical field [T]

    phBy = 0 #Initial Phase of the Vertical field component
    sBy = -1 #Symmetry of the Vertical field component vs Longitudinal position
    xcID = 0 #Transverse Coordinates of Undulator Center [m]
    ycID = 0
    zcID = 0 #Longitudinal Coordinate of Undulator Center [m]

    und = SRWLMagFldU([SRWLMagFldH(1, 'h', Bx, phBy, sBy, 1), SRWLMagFldH(1, 'v', By, phBy, sBy, 1)], undPer, numPer) #Planar Undulator
    magFldCnt = SRWLMagFldC([und], array('d', [xcID]), array('d', [ycID]), array('d', [zcID])) #Container of all Field Elements

    return magFldCnt

def createElectronBeam(undPer, numPer, elettra="E1", electron_energy=2.0):
    #***********Electron Beam
    elecBeam = SRWLPartBeam()

    elecBeam.partStatMom1.x = 0. #Initial Transverse Coordinates (initial Longitudinal Coordinate will be defined later on) [m]
    elecBeam.partStatMom1.y = 0. #-0.00025
    # Roughly ! check!
    elecBeam.partStatMom1.z = -0.5*undPer*(numPer + 8) #Initial Longitudinal Coordinate (set before the ID)
    elecBeam.partStatMom1.xp = 0. #Initial Relative Transverse Velocities
    elecBeam.partStatMom1.yp = 0.
    elecBeam.partStatMom1.gamma = electron_energy/0.51099890221e-03 #Relative Energy

    if electron_energy == 2.0:
        if elettra == "E1" or elettra is None:
            print ("Using phase space of E1 2.0")
            elecBeam.Iavg = 0.32 #Average Current [A]

            #2nd order statistical moments
            elecBeam.arStatMom2[0] = (0.2529e-3)**2 #<(x-x0)^2>
            elecBeam.arStatMom2[1] = 0
            elecBeam.arStatMom2[2] = (0.02881e-3)**2 #<(x'-x'0)^2>
            elecBeam.arStatMom2[3] = (0.01844e-3)**2 #<(y-y0)^2>
            elecBeam.arStatMom2[4] = 0
            elecBeam.arStatMom2[5] = (5.235e-6)**2 #<(y'-y'0)^2>
            # energy spread
            elecBeam.arStatMom2[10] = (0.80e-03)**2 #<(E-E0)^2>/E0^2
        elif elettra == "E24B":
            print ("Using phase space of E24B 2.0")
            elecBeam.Iavg = 0.4 #Average Current [A]

            #2nd order statistical moments
            elecBeam.arStatMom2[0] = (8.071e-05)**2 #<(x-x0)^2>
            elecBeam.arStatMom2[1] = 0
            elecBeam.arStatMom2[2] = (8.968e-06)**2 #<(x'-x'0)^2>
            elecBeam.arStatMom2[3] = (4.66e-06)**2 #<(y-y0)^2>
            elecBeam.arStatMom2[4] = 0
            elecBeam.arStatMom2[5] = (1.553e-06)**2 #<(y'-y'0)^2>
            # energy spread
            elecBeam.arStatMom2[10] = (0.80e-03)**2 #<(E-E0)^2>/E0^2
        elif elettra == "E26B":
            print ("Using phase space of E26B 2.0")
            elecBeam.Iavg = 0.4 #Average Current [A]

            #2nd order statistical moments
            elecBeam.arStatMom2[0] = (5.545e-05)**2 #<(x-x0)^2>
            elecBeam.arStatMom2[1] = 0
            elecBeam.arStatMom2[2] = (4.508e-06)**2 #<(x'-x'0)^2>
            elecBeam.arStatMom2[3] = (2.784e-06)**2 #<(y-y0)^2>
            elecBeam.arStatMom2[4] = 0
            elecBeam.arStatMom2[5] = (8.98e-07)**2 #<(y'-y'0)^2>
            # energy spread
            elecBeam.arStatMom2[10] = (0.80e-03)**2 #<(E-E0)^2>/E0^2
    elif electron_energy == 2.4:
        if elettra == "E1" or elettra is None:
            print ("Using phase space of E1 2.4")
            elecBeam.Iavg = 0.15 #Average Current [A]

            #2nd order statistical moments
            elecBeam.arStatMom2[0] = (0.2842e-3)**2 #<(x-x0)^2>
            elecBeam.arStatMom2[1] = 0
            elecBeam.arStatMom2[2] = (0.03483e-3)**2 #<(x'-x'0)^2>
            elecBeam.arStatMom2[3] = (0.01604e-3)**2 #<(y-y0)^2>
            elecBeam.arStatMom2[4] = 0
            elecBeam.arStatMom2[5] = (6.171e-6)**2 #<(y'-y'0)^2>
            # energy spread
            elecBeam.arStatMom2[10] = (0.95e-03)**2 #<(E-E0)^2>/E0^2


        elif elettra == "E24B":
            print ("Using phase space of E24B 2.4")
            elecBeam.Iavg = 0.4 #Average Current [A]

            #2nd order statistical moments
            elecBeam.arStatMom2[0] = (0.09387e-03)**2 #<(x-x0)^2>
            elecBeam.arStatMom2[1] = 0
            elecBeam.arStatMom2[2] = (0.01055e-03)**2 #<(x'-x'0)^2>
            elecBeam.arStatMom2[3] = (5.513e-06)**2 #<(y-y0)^2>
            elecBeam.arStatMom2[4] = 0
            elecBeam.arStatMom2[5] = (1.796e-06)**2 #<(y'-y'0)^2>
            # energy spread
            elecBeam.arStatMom2[10] = (0.90e-03)**2 #<(E-E0)^2>/E0^2
        elif elettra == "E26B":
            print ("Using phase space of E26B 2.4")
            elecBeam.Iavg = 0.4 #Average Current [A]

            #2nd order statistical moments
            elecBeam.arStatMom2[0] = (0.05601e-3)**2 #<(x-x0)^2>
            elecBeam.arStatMom2[1] = 0
            elecBeam.arStatMom2[2] = (6.364e-06)**2 #<(x'-x'0)^2>
            elecBeam.arStatMom2[3] = (2.814e-06)**2 #<(y-y0)^2>
            elecBeam.arStatMom2[4] = 0
            elecBeam.arStatMom2[5] = (1.267e-06)**2 #<(y'-y'0)^2>
            # energy spread
            elecBeam.arStatMom2[10] = (0.81e-03)**2 #<(E-E0)^2>/E0^2


    return elecBeam

def createInitialWavefrontMesh(elecBeam, photon_energy):
    #****************** Initial Wavefront
    wfr = SRWLWfr() #For intensity distribution at fixed photon energy
    wfr.allocate(1, 101, 101) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
    wfr.mesh.zStart = 11.517 #Longitudinal Position [m] from Center of Straight Section at which SR has to be calculated
    wfr.mesh.eStart = photon_energy #Initial Photon Energy [eV]
    wfr.mesh.eFin = wfr.mesh.eStart #Final Photon Energy [eV]

    #OC: it makes sense to choose the initial ranges equal or a bit larger that the first aperture of the beamline (not much larger)

    half_window = 0.5

    wfr.mesh.xStart = -half_window*1e-3 #-0.00015 #Initial Horizontal Position [m]
    wfr.mesh.xFin = -1 * wfr.mesh.xStart #0.00015 #Final Horizontal Position [m]
    wfr.mesh.yStart = -half_window*1e-3 #-0.00015 #Initial Vertical Position [m]
    wfr.mesh.yFin = -1 * wfr.mesh.yStart#0.00015 #Final Vertical Position [m]

    wfr.partBeam = elecBeam

    return wfr

def createCalculationPrecisionSettings():
    #***********Precision Parameters for SR calculation
    meth = 1 #SR calculation method: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
    relPrec = 0.01 #relative precision
    zStartInteg = 0 #longitudinal position to start integration (effective if < zEndInteg)
    zEndInteg = 0 #longitudinal position to finish integration (effective if > zStartInteg)
    npTraj = 100000 #Number of points for trajectory calculation
    useTermin = 1 #Use "terminating terms" (i.e. asymptotic expansions at zStartInteg and zEndInteg) or not (1 or 0 respectively)
    arPrecParSpec = [meth, relPrec, zStartInteg, zEndInteg, npTraj, useTermin, 0]

    return arPrecParSpec

def get_orientation_vectors(grazing_angle=0.003, orientation_of_reflection_plane=Orientation.UP):
    if orientation_of_reflection_plane == Orientation.LEFT:
        nvx = -numpy.cos(grazing_angle)
        nvy = 0
        nvz = -numpy.sin(grazing_angle)
        tvx = numpy.sin(grazing_angle)
        tvy = 0
    elif orientation_of_reflection_plane == Orientation.RIGHT:
        nvx = numpy.cos(grazing_angle)
        nvy = 0
        nvz = -numpy.sin(grazing_angle)
        tvx = -numpy.sin(grazing_angle)
        tvy = 0
    elif orientation_of_reflection_plane == Orientation.UP:
        nvx = 0
        nvy = numpy.cos(grazing_angle)
        nvz = -numpy.sin(grazing_angle)
        tvx = 0
        tvy = numpy.sin(grazing_angle)
    elif orientation_of_reflection_plane == Orientation.DOWN:
        nvx = 0
        nvy = -numpy.cos(grazing_angle)
        nvz = -numpy.sin(grazing_angle)
        tvx = 0
        tvy = -numpy.sin(grazing_angle)

    return nvx, nvy, nvz, tvx, tvy

def createBeamline(n_elements, photon_energy, figure_errors=0):

    wavelength = 1.239842e-6/photon_energy

    #***************** Optical Elements and Propagation Parameters

    # FRONT END MASK --------------------------------------------------
    frontend_mask = SRWLOptA('r', 'a', 650e-06, 730e-06)

    # DRIFT TO TORUS --------------------------------------------------
    drift_mask_torus = SRWLOptD(4.483)

    # TORUS -----------------------------------------------------------
    incident_angle = 88.25
    grazing_angle = numpy.radians(90.-incident_angle)

    nvx, nvy, nvz, tvx, tvy = get_orientation_vectors(grazing_angle=grazing_angle,
                                                      orientation_of_reflection_plane=Orientation.LEFT)

    torus = SRWLOptMirTor(_rt=116.50,
                          _rs=0.1543,
                          _size_tang=0.25,
                          _size_sag=0.01,
                          _nvx = nvx,
                          _nvy = nvy,
                          _nvz = nvz,
                          _tvx = tvx,
                          _tvy = tvy)

    #height_profile_data = srwl_uti_read_data_cols("../torus_height_error_profile_1D.dat",
    #                                              _str_sep='\t',
    #                                              _i_col_start=0,
    #                                              _i_col_end=1)

    height_profile_data = srwl_uti_read_data_cols("../torus_height_error_profile_2D.dat",
                                                  _str_sep='\t')

    optTrEr_torus = srwl_opt_setup_surf_height_1d(_height_prof_data=height_profile_data,
                                                  _ang=grazing_angle,
                                                  _dim='x',
                                                  _amp_coef=1.0)

    # DRIFT TO HORIZONTAL FOCUS ----------------------------------------
    drift_torus_hf = SRWLOptD(2.0)

    # DRIFT TO VERTICAL FOCUS ------------------------------------------
    drift_torus_vf = SRWLOptD(1.0)

    # ENTRANCE SLIT ----------------------------------------------------
    entrance_slit = SRWLOptA('r', 'a', 0.1, 20e-06)

    # DRIFT TO VLS GRATING ---------------------------------------------
    drift_slit_vls = SRWLOptD(2.0)

    # VLS GRATING ------------------------------------------------------
    alpha_shadow = 87.98094647
    vls_grooving = [400., 0.8, -0.000299, 0.0, 0.0] #Polynomial coefficients for VLS Grating Groove Density
    alpha_grating = numpy.radians(90.-alpha_shadow)

    beta_grating = asin(wavelength*vls_grooving[0]*1.e+03 - cos(alpha_grating)) # Grating Output Angle
    deflection_angle = alpha_grating + beta_grating + 1.57079632679 # Grating Deflection Angle

    print("BETA GRATING", -numpy.degrees(beta_grating))

    nvx, nvy, nvz, tvx, tvy = get_orientation_vectors(grazing_angle=alpha_grating,
                                                      orientation_of_reflection_plane=Orientation.DOWN)

    grating_substrate = SRWLOptMirPl(_size_tang=0.08,
                                     _size_sag=0.01,
                                     _ap_shape='r',
                                     _nvx = nvx,
                                     _nvy = nvy,
                                     _nvz = nvz,
                                     _tvx = tvx,
                                     _tvy = tvy)

    #The VLS Grating itself:
    vls_grating = SRWLOptG(_mirSub=grating_substrate,
                           _m=1,
                           _grDen=vls_grooving[0],
                           _grDen1=vls_grooving[1],
                           _grDen2=vls_grooving[2],
                           _grDen3=vls_grooving[3],
                           _grDen4=vls_grooving[4],
                           _grAng=0.0)

    height_profile_data = srwl_uti_read_data_cols("../VLS_height_error_profile_2D.dat",
                                                  _str_sep='\t')

    optTrEr_grating = srwl_opt_setup_surf_height_2d(_height_prof_data=height_profile_data,
                                                    _ang=alpha_grating,
                                                    _ang_r=deflection_angle,
                                                    _dim='y',
                                                    _amp_coef=1.0)

    # DRIFT TO EXIT SLIT ---------------------------------------------
    drift_vls_slit = SRWLOptD(2.0)

    # EXIT SLIT (GRATING FOCUS) --------------------------------------
    exit_slit = SRWLOptA('r', 'a', 0.02, 10e-06)

    # DRIFT TO KB-VRM ------------------------------------------------
    drift_slit_kbvrm = SRWLOptD(8.3)

    # KB-VRM ---------------------------------------------------------
    incident_angle = 88.5
    grazing_angle = numpy.radians(90.-incident_angle)

    nvx, nvy, nvz, tvx, tvy = get_orientation_vectors(grazing_angle=grazing_angle,
                                                      orientation_of_reflection_plane=Orientation.UP)
    kbvrm = SRWLOptMirEl(_p=8.3,
                         _q=1.65,
                         _ang_graz=grazing_angle,
                         _size_tang=0.3,
                         _size_sag=0.01,
                         _nvx = nvx,
                         _nvy = nvy,
                         _nvz = nvz,
                         _tvx = tvx,
                         _tvy = tvy)

    height_profile_data = srwl_uti_read_data_cols("../KBVRM_height_error_profile_2D.dat",
                                                  _str_sep='\t')

    optTrEr_kbvrm = srwl_opt_setup_surf_height_2d(_height_prof_data=height_profile_data,
                                                  _ang=grazing_angle,
                                                  _dim='y',
                                                  _amp_coef=1.0)

    # DRIFT TO KB-HRM ------------------------------------------------

    drift_kbvrm_kbhrm = SRWLOptD(0.45)

    # KB-HRM ---------------------------------------------------------
    incident_angle = 88.0
    grazing_angle = numpy.radians(90.-incident_angle)

    nvx, nvy, nvz, tvx, tvy = get_orientation_vectors(grazing_angle=grazing_angle,
                                                      orientation_of_reflection_plane=Orientation.LEFT)

    kbhrm = SRWLOptMirEl(_p=13.75,
                         _q=1.2,
                         _ang_graz=grazing_angle,
                         _size_tang=0.3,
                         _size_sag=0.01,
                         _nvx = nvx,
                         _nvy = nvy,
                         _nvz = nvz,
                         _tvx = tvx,
                         _tvy = tvy)

    height_profile_data = srwl_uti_read_data_cols("../KBHRM_height_error_profile_2D.dat",
                                                  _str_sep='\t')

    optTrEr_kbhrm = srwl_opt_setup_surf_height_2d(_height_prof_data=height_profile_data,
                                                  _ang=grazing_angle,
                                                  _dim='x',
                                                  _amp_coef=1.0)

    # DRIFT TO SAMPLE ------------------------------------------------
    drift_kbhrm_sample = SRWLOptD(1.2)

    #Wavefront Propagation Parameters:
    #[ 0]: Auto-Resize (1) or not (0) Before propagation
    #[ 1]: Auto-Resize (1) or not (0) After propagation
    #[ 2]: Relative Precision for propagation with Auto-Resizing (1. is nominal)
    #[ 3]: Allow (1) or not (0) for semi-analytical treatment of the quadratic (leading) phase terms at the propagation
    #[ 4]: Do any Resizing on Fourier side, using FFT, (1) or not (0)
    #[ 5]: Horizontal Range modification factor at Resizing (1. means no modification)
    #[ 6]: Horizontal Resolution modification factor at Resizing
    #[ 7]: Vertical Range modification factor at Resizing
    #[ 8]: Vertical Resolution modification factor at Resizing
    #[ 9]: Type of wavefront Shift before Resizing (not yet implemented)
    #[10]: New Horizontal wavefront Center position after Shift (not yet implemented)
    #[11]: New Vertical wavefront Center position after Shift (not yet implemented)
    #[12]: Optional: Orientation of the Output Optical Axis vector in the Incident Beam Frame: Horizontal Coordinate
    #[13]: Optional: Orientation of the Output Optical Axis vector in the Incident Beam Frame: Vertical Coordinate
    #[14]: Optional: Orientation of the Output Optical Axis vector in the Incident Beam Frame: Longitudinal Coordinate
    #[15]: Optional: Orientation of the Horizontal Base vector of the Output Frame in the Incident Beam Frame: Horizontal Coordinate
    #[16]: Optional: Orientation of the Horizontal Base vector of the Output Frame in the Incident Beam Frame: Vertical Coordinate

    #OC: propagation through drifts is often better done with "semi-analytical treatment of the quadratic (leading) phase terms"

    #                       [ 0] [1] [2]  [3] [4] [5]  [6]    [7]  [8]    [9] [10] [11]
    pp_frontend_mask      = [ 0,  0, 1.0,  0,  0, 1.5,  2.0,  1.5,  5.0,  0,  0,  0] #OC_cor
    pp_drift_mask_torus   = [ 0,  0, 1.0,  1,  0, 1.0,  1.0,  1.0,  1.0,  0,  0,  0] #OC_cor
    pp_torus              = [ 0,  0, 1.0,  0,  0, 1.0,  1.0,  1.0,  1.0,  0,  0,  0]
    pp_error_torus        = [ 0,  0, 1.0,  0,  0, 1.0,  1.0,  1.0,  1.0,  0,  0,  0]
    pp_drift_torus_hf     = [ 0,  0, 1.0,  1,  0, 1.0,  4.0,  1.0,  8.0,  0,  0,  0]
    pp_drift_torus_vf     = [ 0,  0, 1.0,  1,  0, 2.0,  1.0,  1.0,  1.0,  0,  0,  0]
    #OC_cor: increasing resolution would be better done before the slit (edges will be "sharper"):
    pp_entrance_slit      = [ 0,  0, 1.0,  0,  0, 1.0,  1.0,  0.02, 5.0,  0,  0,  0]
    pp_drift_slit_vls     = [ 0,  0, 1.0,  1,  0, 1.0,  1.0,  1.0,  1.0,  0,  0,  0]
    #OC_cor: deflection of the optical axis should be DOWN as well!
    pp_vls_grating        = [ 0,  0, 1.0,  0,  0, 1.0,  1.0,  1.0,  1.0,  0,  0,  0,  0,  -sin(deflection_angle),  cos(deflection_angle),  1,  0]
    pp_error_grating      = [ 0,  0, 1.0,  0,  0, 1.0,  1.0,  1.0,  1.0,  0,  0,  0]
    pp_drift_vls_slit     = [ 0,  0, 1.0,  1,  0, 1.0,  1.0,  1.0,  1.0,  0,  0,  0]
    pp_exit_slit          = [ 0,  0, 1.0,  0,  0, 1.0,  1.0,  0.5,  10.0, 0,  0,  0] #OC_cor
    pp_drift_slit_kbvrm   = [ 0,  0, 1.0,  1,  0, 1.0,  1.0,  1.0,  1.0,  0,  0,  0] #OC_cor
    pp_kbvrm              = [ 0,  0, 1.0,  0,  0, 1.0,  1.0,  0.1,  1.0,  0,  0,  0] #OC_cor
    pp_error_kbvrm        = [ 0,  0, 1.0,  0,  0, 1.0,  1.0,  1.0,  1.0,  0,  0,  0]
    pp_drift_kbvrm_kbhrm  = [ 0,  0, 1.0,  1,  0, 1.0,  1.0,  1.0,  1.0,  0,  0,  0]
    pp_kbhrm              = [ 0,  0, 1.0,  0,  0, 1.0,  1.0,  1.0,  1.0,  0,  0,  0]
    pp_error_kbhrm        = [ 0,  0, 1.0,  0,  0, 1.0,  1.0,  1.0,  1.0,  0,  0,  0]
    #OC_cor (special propagator - "to waist", impact of range and resolution resizing params are reversed for it!):
    pp_drift_kbhrm_sample = [ 0,  0, 1.0,  4,  0, 1.0,  1.0,  3.0,  1.0,  0,  0,  0]


    #"Beamline" - Container of Optical Elements (together with the corresponding wavefront propagation instructions)

    if figure_errors == 0:
        print("BEAMLINE WITH NO HEIGHT ERROR PROFILES")

        optBL = SRWLOptC([frontend_mask,        #1
                          drift_mask_torus,     #2
                          torus,                #3
                          drift_torus_hf,       #4
                          drift_torus_vf,       #5
                          entrance_slit,        #6
                          drift_slit_vls,       #7
                          vls_grating,          #8
                          drift_vls_slit,       #9
                          exit_slit,            #10
                          drift_slit_kbvrm,     #11
                          kbvrm,                #12
                          drift_kbvrm_kbhrm,    #13
                          kbhrm,                #14
                          drift_kbhrm_sample    #15
                          ][:n_elements],
                         [pp_frontend_mask,     #1
                          pp_drift_mask_torus,  #2
                          pp_torus,             #3
                          pp_drift_torus_hf,    #4
                          pp_drift_torus_vf,    #5
                          pp_entrance_slit,     #6
                          pp_drift_slit_vls,    #7
                          pp_vls_grating,       #8
                          pp_drift_vls_slit,    #9
                          pp_exit_slit,         #10
                          pp_drift_slit_kbvrm,  #11
                          pp_kbvrm,             #12
                          pp_drift_kbvrm_kbhrm, #13
                          pp_kbhrm,             #14
                          pp_drift_kbhrm_sample #15
                          ][:n_elements])
    elif figure_errors == 1:
        print("BEAMLINE WITH HEIGHT ERROR PROFILES - NO VLS")

        optBL = SRWLOptC([frontend_mask,        #1
                          drift_mask_torus,     #2
                          torus,                #3
                          optTrEr_torus,        #4
                          drift_torus_hf,       #5
                          drift_torus_vf,       #6
                          entrance_slit,        #7
                          drift_slit_vls,       #8
                          vls_grating,          #9
                          drift_vls_slit,       #10
                          exit_slit,            #11
                          drift_slit_kbvrm,     #12
                          kbvrm,                #13
                          optTrEr_kbvrm,        #14
                          drift_kbvrm_kbhrm,    #15
                          kbhrm,                #16
                          optTrEr_kbhrm,        #17
                          drift_kbhrm_sample    #18
                          ][:n_elements],
                         [pp_frontend_mask,     #1
                          pp_drift_mask_torus,  #2
                          pp_torus,             #3
                          pp_error_torus,       #4
                          pp_drift_torus_hf,    #5
                          pp_drift_torus_vf,    #6
                          pp_entrance_slit,     #7
                          pp_drift_slit_vls,    #8
                          pp_vls_grating,       #9
                          pp_drift_vls_slit,    #10
                          pp_exit_slit,         #11
                          pp_drift_slit_kbvrm,  #12
                          pp_kbvrm,             #13
                          pp_error_kbvrm,       #14
                          pp_drift_kbvrm_kbhrm, #15
                          pp_kbhrm,             #16
                          pp_error_kbhrm,       #17
                          pp_drift_kbhrm_sample #18
                          ][:n_elements])
    elif figure_errors == 2:
        print("BEAMLINE WITH HEIGHT ERROR PROFILES + VLS")

        optBL = SRWLOptC([frontend_mask,        #1
                          drift_mask_torus,     #2
                          torus,                #3
                          optTrEr_torus,        #4
                          drift_torus_hf,       #5
                          drift_torus_vf,       #6
                          entrance_slit,        #7
                          drift_slit_vls,       #8
                          vls_grating,          #9
                          optTrEr_grating,      #10
                          drift_vls_slit,       #11
                          exit_slit,            #12
                          drift_slit_kbvrm,     #13
                          kbvrm,                #14
                          optTrEr_kbvrm,        #15
                          drift_kbvrm_kbhrm,    #16
                          kbhrm,                #17
                          optTrEr_kbhrm,        #18
                          drift_kbhrm_sample    #19
                          ][:n_elements],
                         [pp_frontend_mask,     #1
                          pp_drift_mask_torus,  #2
                          pp_torus,             #3
                          pp_error_torus,       #4
                          pp_drift_torus_hf,    #5
                          pp_drift_torus_vf,    #6
                          pp_entrance_slit,     #7
                          pp_drift_slit_vls,    #8
                          pp_vls_grating,       #9
                          pp_error_grating,     #10
                          pp_drift_vls_slit,    #11
                          pp_exit_slit,         #12
                          pp_drift_slit_kbvrm,  #13
                          pp_kbvrm,             #14
                          pp_error_kbvrm,       #15
                          pp_drift_kbvrm_kbhrm, #16
                          pp_kbhrm,             #17
                          pp_error_kbhrm,       #18
                          pp_drift_kbhrm_sample #19
                          ][:n_elements])
    return optBL


if __name__== "__main__":

    ##############################
    # 1 front-end Mask
    # 2 to torus
    # 5 to HF
    # 6 to VF
    # 18 to final focus
    n_propagation_elements = int(sys.argv[1])

    ##############################
    # E1   = Elettra 1.0
    # E24B = Elettra 2.0 4-bend achromat
    # E26B = Elettra 2.0 6-bend achromat
    elettra = sys.argv[2]

    ##############################
    # 20 = 2.0 GeV
    # 24 = 2.4 GeV
    ee = sys.argv[3]
    electron_energy = int(ee)/10

    ##############################
    # 4  = 400 eV
    # 10 = 1000 eV
    photon_energy = int(sys.argv[4])

    pe = str(photon_energy)

    if photon_energy == 400:
        photon_energy = 399.8 # 3rd Harmonic
        if electron_energy == 2.0:
            Kv = 1.92
        elif electron_energy == 2.4:
            Kv = 2.4864
    elif photon_energy == 1000:
        photon_energy = 1000 # 5th Harmonic
        if electron_energy == 2.0:
            Kv = 0.0 #TODO: missing value
        elif electron_energy == 2.4:
            Kv = 1.8569

    print("Undulator", photon_energy, Kv)

    ##############################
    # 0   = Total Intensity
    # 4   = Mutual intensity cuts x/y
    plots = int(sys.argv[5])

    ##############################
    # 0   = No Errors
    # 1   = Errors, No VLS
    # 2   = Errors ALL
    figure_errors = int(sys.argv[6])

    print("input parameters")
    print("NE", n_propagation_elements)
    print("Elettra", elettra)
    print("EE", electron_energy)
    print("PE", photon_energy)
    print("Plots", plots)
    print("Figure Error", figure_errors)

    do_ME = (plots >= 0)

    outdir = "./" + elettra + "_" + ("ME" if do_ME else "SE")

    outdir += "_EE" + ee
    outdir += "_PE" + pe

    if n_propagation_elements == 2:
        outdir += "_FE"
        if plots == 4: outdir += "_degcoh"
    elif n_propagation_elements == 5:
        outdir += "_HF"
    elif n_propagation_elements == 6:
        outdir += "_VF"
    else:
        outdir += "_NE" + str(n_propagation_elements)

    print("Switching working directory to: " + outdir)

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    os.chdir(outdir)

    magFldCnt = createUndulator(Kv, undPer, numPer)
    elecBeam = createElectronBeam(undPer, numPer, elettra, electron_energy)
    wfr = createInitialWavefrontMesh(elecBeam, photon_energy)

    meshInitPartCoh = deepcopy(wfr.mesh)

    optBL = createBeamline(n_propagation_elements, photon_energy, figure_errors)

    arPrecParSpec = createCalculationPrecisionSettings()

    # This is the convergence parameter. Higher is more accurate but slower!!
    # 0.2 is used in the original example. But I think it should be higher. The calculation may then however need too much memory.
    sampFactNxNyForProp = 1.0 #0.6 #sampling factor for adjusting nx, ny (effective if > 0)

    if(not do_ME and srwl_uti_proc_is_master()):
        print('   Performing Initial Single-E Electric Field calculation ... ', end='')
        arPrecParSpec[6] = sampFactNxNyForProp #sampling factor for adjusting nx, ny (effective if > 0)
        srwl.CalcElecFieldSR(wfr, 0, magFldCnt, arPrecParSpec)
        print('done')
        print('   Extracting Intensity from the Calculated Initial Electric Field ... ', end='')
        arI = array('f', [0]*wfr.mesh.nx*wfr.mesh.ny) #"flat" array to take 2D intensity data
        srwl.CalcIntFromElecField(arI, wfr, 6, 0, 3, wfr.mesh.eStart, 0, 0)
        print('done')
        print('   Saving the Initial Wavefront Intensity into a file ... ', end='')
        srwl_uti_save_intens_ascii(arI, wfr.mesh, strIntSE_OutFileName)
        print('done')

        print('   Simulating Electric Field Wavefront Propagation ... ', end='')
        t0 = time.time();
        srwl.PropagElecField(wfr, optBL)
        print('done; lasted', round(time.time() - t0), 's')

        print('   Extracting Intensity from the Propagated Electric Field  ... ', end='')
        arI1 = array('f', [0]*wfr.mesh.nx*wfr.mesh.ny) #"flat" 2D array to take intensity data
        srwl.CalcIntFromElecField(arI1, wfr, 6, 0, 3, wfr.mesh.eStart, 0, 0)

        print('done')
        print('   Saving the Propagated Wavefront Intensity data to a file ... ', end='')
        srwl_uti_save_intens_ascii(arI1, wfr.mesh, strIntPropSE_OutFileName)
        print('done')


    if do_ME:
        # Remove comments to start multi electron.
        # print('   Starting simulation of Partially-Coherent Wavefront Propagation (takes a lot of time)... ')
        nMacroElec = 500000 #T otal number of Macro-Electrons (Wavefronts)
        nMacroElecAvgOneProc = 5 # Number of Macro-Electrons (Wavefronts) to average on each node (for MPI calculations)
        nMacroElecSavePer = 20 # Saving periodicity (in terms of Macro-Electrons) for the Resulting Intensity
        srCalcMeth = 1 # SR calculation method (1- undulator)
        srCalcPrec = 0.01 # SR calculation rel. accuracy

        radStokesProp = srwl_wfr_emit_prop_multi_e(elecBeam,
                                                   magFldCnt,
                                                   meshInitPartCoh,
                                                   srCalcMeth,
                                                   srCalcPrec,
                                                   nMacroElec,
                                                   nMacroElecAvgOneProc,
                                                   nMacroElecSavePer,
                                                   strIntPropME_OutFileName,
                                                   sampFactNxNyForProp,
                                                   optBL,
                                                   _char=plots)
        print('done')
