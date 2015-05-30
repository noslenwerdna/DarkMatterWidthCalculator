#!/usr/bin/env python

import sys
import argparse as ap

def quit(p):
    p.print_help()
    sys.exit()

def wCalc(type, gQ, gDM, mMed, mDM):
    ''' 
        Calculates the total width for s-channel mediated dark matter:
        q              dark matter
         \\   med.   //
           ========== 
         //          \\
        q              dark matter
        Gamma_min = Gamma_DM + sum_q Gamma_qq

        type: indicates whether the mediator has vector or axial-vector symmetry
        gQ: coupling of mediator to quarks
        gDM: coupling of mediator to dark matter
        mMed: mass of mediator in GeV
        mDM: mass of dark matter in GeV
    '''

    import sys
    import math as m

    if not type == 'vector' and not type == 'axial':
        print 'type must be axial or vector'
        sys.exit()

    # Partial width has the form: normFactor * ratioFunc * beta**exp

    # Quark Partial Width Calculation
    normFactorQQ = (3. * gQ**2 * mMed)/(12. * m.pi)
    quarkName = ['u', 'd', 'c', 's', 't', 'b']
    quarkMass = [2.55e-3, 5.04e-3, 1.42, 1.01e-1, 172., 4.70] # in GeV
    ratioFuncQQ, betaQQ = list(), list()
    betaExp = 0.
    if type == 'vector':
        betaExp = 0.5
        ratioFuncQQ = [(1. + (2. * mass**2)/mMed**2) for mass in quarkMass]
    else: # type is axial
        betaExp = 1.5
        ratioFuncQQ = [1. for mass in quarkMass]
    for mass in quarkMass:
        if mMed > 2.*mass:
            betaQQ.append((1. - (2. * mass)**2/mMed**2)**betaExp)
        else: 
            betaQQ.append(0.)
    Gamma_qq = [normFactorQQ * r * b for r, b in zip(ratioFuncQQ, betaQQ)]
    sumGamma_qq = 0.
    for w in Gamma_qq:
        sumGamma_qq += w

    # Dark Matter Partial Width Calculation
    normFactorDM = (gDM**2 * mMed)/(12. * m.pi)
    ratioFuncDM = 0.
    if type == 'vector':
        ratioFuncDM = (1. + (2. * mDM**2)/mMed**2)
    else: # type is axial 
        ratioFuncDM = 1.
    if mMed > 2.*mDM:
        betaDM = (1. - (2. * mDM)**2/mMed**2)**betaExp
    else:
        betaDM = 0.
    Gamma_DM = normFactorDM * ratioFuncDM * betaDM

    return Gamma_DM + sumGamma_qq

if __name__ == '__main__':

    parser = ap.ArgumentParser(description='Calculates the total width '
    ' for s-channel mediated dark matter: \n'
    '   q              dark matter \n'
    '    \\\\   med.   // \n'
    '      ========== \n'
    '    //          \\\\ \n'
    '   q              dark matter \n'
    '   Gamma_min = Gamma_DM + sum_q Gamma_qq \n\n'
    '   type: indicates whether the mediator has vector or '
    ' axial-vector symmetry \n'
    '   gQ: coupling of mediator to quarks \n'
    '   gDM: coupling of mediator to dark matter \n'
    '   mMed: mass of mediator in GeV \n'
    '   mDM: mass of dark matter in GeV \n', 
    formatter_class=ap.RawTextHelpFormatter)
    parser.add_argument(
        '--type', type=str,  
        help='Symmetry of mediator. Must be vector or axial')
    parser.add_argument(
        '--gQ', type=float, default=1., help='Coupling of mediator to quarks')
    parser.add_argument(
        '--gDM', type=float, default=1.,
        help='Coupling of mediator to dark matter')
    parser.add_argument(
        '--mMed', type=float, help='Mass of mediator')
    parser.add_argument(
        '--mDM', type=float, help='Mass of mediator')
    args = parser.parse_args()

    if len(sys.argv) == 1:
        quit(parser)
    if not args.type == 'vector' and not args.type =='axial':
        print '\ntype must be axial or vector\n\n'
        quit(parser)

    print 'full mediator width={0:0.2f} GeV'.format(
        wCalc(args.type, args.gQ, args.gDM, args.mMed, args.mDM))
