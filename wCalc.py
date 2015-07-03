#!/usr/bin/env python

import sys
import argparse as ap

def quit(p):
    p.print_help()
    sys.exit()

def swCalc(type, gQ, gDM, mMed, mDM):
    ''' 
        Calculates the total width for s-channel mediated dark matter:
        q              dark matter
         \\   med.   //
           ========== 
         //          \\
        q              dark matter
        Gamma_min = Gamma_DM + sum_q Gamma_qq

        type: indicates type of the mediator (scalar, pseudo, vector, axial)
        gQ: coupling of mediator to quarks
        gDM: coupling of mediator to dark matter
        mMed: mass of mediator in GeV
        mDM: mass of dark matter in GeV
    '''

    import sys
    import math as m

    # Partial width has the form: normFactor * ratioFunc * betaFact

    # Quark Partial Width Calculation. Calculate for each quark in list.
    normFactorQQ = 1.
    vev = 246.
    betaExp = 0.
    quarkName = ['u', 'd', 'c', 's', 't', 'b']
    quarkMass = [2.3e-3, 4.8e-3, 1.275, 9.5e-2, 173.5, 4.65] # in GeV
    ratioFuncQQ, betaQQ = list(), list() # lists containing entry for each quark
    if type == 'vector' or type == 'axial':
        normFactorQQ = (3. * gQ**2 * mMed)/(12. * m.pi)
        if type == 'vector':
            betaExp = 0.5
            ratioFuncQQ = [(1. + (2. * mass**2)/mMed**2) for mass in quarkMass]
        elif type == 'axial': 
            betaExp = 1.5
            ratioFuncQQ = [1. for mass in quarkMass]
        for mass in quarkMass:
            if mMed > 2.*mass:
                betaQQ.append((1. - (2. * mass)**2/mMed**2)**betaExp)
            else: 
                betaQQ.append(0.)
    elif type == 'scalar' or type == 'pseudo':
        normFactorQQ = (3. * gQ**2 * mMed)/(8. * m.pi * vev**2)
        if type == 'scalar':
            betaExp = 1.5
        elif type == 'pseudo':
            betaExp = 0.5
        for mass in quarkMass:
            if mMed > 2.*mass:
                betaQQ.append(1.)
                ratioFuncQQ.append(mass**2 * 
                                   (1. - (2. * mass)**2/mMed**2)**betaExp)
            else: 
                betaQQ.append(0.)
                ratioFuncQQ.append(0.)
    else:
        print "Invalid type"
        return 0
    # sum over all quarks
    Gamma_qq = [normFactorQQ * b * r for r, b in zip(ratioFuncQQ, betaQQ)]
    sumGamma_qq = 0.
    for w in Gamma_qq:
        sumGamma_qq += w

    # Dark Matter Partial Width Calculation
    ratioFuncDM = 0.
    if type == 'vector' or type == 'axial':
        normFactorDM = (gDM**2 * mMed)/(12. * m.pi)
        if type == 'vector':
            ratioFuncDM = (1. + (2. * mDM**2)/mMed**2)
        elif type == 'axial': 
            ratioFuncDM = 1.
        if mMed > 2.*mDM:
            betaDM = (1. - (2. * mDM)**2/mMed**2)**betaExp
        else:
            betaDM = 0.
    if type == 'scalar' or type == 'pseudo':
        normFactorDM = (gDM**2 * mMed)/(8. * m.pi)
        if type == 'scalar':
            betaExp = 1.5
        elif type == 'pseudo':
            betaExp = 0.5
        ratioFuncDM = 0.
        if mMed > 2.*mDM:
            betaDM = 1.
            ratioFuncDM = (1. - (2. * mDM)**2/mMed**2)**betaExp
        else:
            betaDM = 0.
            ratioFuncDM = 0.
    Gamma_DM = normFactorDM * ratioFuncDM * betaDM

    return Gamma_DM + sumGamma_qq

def twCalc(type, gQ, gDM, mMed, mDM):
    ''' 
        Calculates the total width for t-channel mediated dark matter:
        q ============ dark matter
               ||
               || med.
               ||
        q ============ dark matter
        Mediator is different for each type of quark. 
        Simplify and assume that the mediators have all the same mass.

        type: indicates type of the mediator (tchan)
        gQ: coupling of mediator to quarks
        gDM: coupling of mediator to dark matter
        mMed: mass of mediator in GeV
        mDM: mass of dark matter in GeV
    '''

    import sys
    import math as m

    if gQ != gDM:
        print "gQ must equal gDM for t-channel"
        return 0

    # Partial width has the form: normFactor * massFactor * rootFactor

    # Width Calculation. One mediator, one width for each quark in list.
    normFactor = 1.
    quarkMass = { 
        'u': 2.3e-3,
        'd': 4.8e-3,
        'c': 1.275,
        's': 9.5e-2,
        't': 173.5,
        'b': 4.65,
    }
    width, massFactor, rootFactor = dict(), dict(), dict() 
    if type == 'tchan':
        normFactor = (gDM**2 * mMed)/(16. * m.pi)
        mDMRatio = mDM/mMed
        for quark in quarkMass:
            mQRatio = quarkMass[quark]/mMed
            if mMed > 2.*mDM and mMed > 2.*quarkMass[quark]:
                massFactor[quark] = 1 - mDMRatio**2 - mQRatio**2
                massDiff = mQRatio - mDMRatio
                massSum = mQRatio + mDMRatio
                rootFactor[quark] = (1. - massSum**2) * (1. - massDiff**2)
            elif mMed > 2.*mDM: # implying mMed <= 2.*quarkMass[quark]
                massFactor[quark] = 1 - mDMRatio**2
                rootFactor[quark] = (1. - mDMRatio**2) * (1. - mDMRatio**2)
            elif mMed > 2.*quarkMass[quark]: # implying mMed <= 2.*mDM
                massFactor[quark] = 1 - mQRatio**2
                rootFactor[quark] = (1. - mQRatio**2) * (1. - mQRatio**2)
            else:
                massFactor[quark] = 0.
                rootFactor[quark] = 0.
            width[quark] = (normFactor * massFactor[quark] * 
                            m.sqrt(rootFactor[quark]))

    return width

if __name__ == '__main__':

    parser = ap.ArgumentParser(description='Calculates the total width '
    ' for s-channel mediated dark matter: \n'
    '   q              dark matter \n'
    '    \\\\   med.   // \n'
    '      ========== \n'
    '    //          \\\\ \n'
    '   q              dark matter \n'
    '   Gamma_min = Gamma_DM + sum_q Gamma_qq \n\n'

    ' for t-channel mediated dark matter: \n'
    '   q ============ dark matter \n'
    '          || \n'
    '          || med. \n'
    '          || \n'
    '   q ============ dark matter \n'
    '   Mediator is different for each type of quark. \n'
    '   Simplify and assume that t-chan mediators have all the same mass. \n\n'

    '   type: indicates type of mediator: vector, axial, scalar, pseudo, '
    ' tchan \n'
    '   gQ: coupling of mediator to quarks \n'
    '   gDM: coupling of mediator to dark matter \n'
    '   mMed: mass of mediator in GeV \n'
    '   mDM: mass of dark matter in GeV \n', 
    formatter_class=ap.RawTextHelpFormatter)
    parser.add_argument(
        '--type', type=str,  
        help=('Symmetry of mediator. '
              'Must be vector, axial, scalar, pseudo, or tchan'))
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
    if not (args.type == 'vector' or args.type == 'axial' or
            args.type == 'scalar' or args.type == 'pseudo' or
            args.type == 'tchan'):
        print '\ntype must be axial, vector, scalar, pseudo, tchan\n\n'
        quit(parser)

    if args.type != 'tchan':
        print 'full mediator width={0:0.3e} GeV'.format(
            swCalc(args.type, args.gQ, args.gDM, args.mMed, args.mDM))
    else: # tchannel 
        width = twCalc(args.type, args.gQ, args.gDM, args.mMed, args.mDM)
        print 'eta_{0} width={1:0.8e} GeV'.format('u', width['u'])
        print 'eta_{0} width={1:0.8e} GeV'.format('d', width['d'])
        print 'eta_{0} width={1:0.8e} GeV'.format('c', width['c'])
        print 'eta_{0} width={1:0.8e} GeV'.format('s', width['s'])
        print 'eta_{0} width={1:0.8e} GeV'.format('t', width['t'])
        print 'eta_{0} width={1:0.8e} GeV'.format('b', width['b'])
