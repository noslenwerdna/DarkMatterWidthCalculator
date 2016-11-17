# DarkMatterWidthCalculator

        Calculates the total width for s-channel mediated dark matter:
        q              dark matter
         \\   med.   //
           ========== 
         //          \\
        q              dark matter
        Gamma_min = Gamma_DM + sum_q Gamma_qq + sum_l Gamma_ll

        type: indicates whether the mediator has vector or axial-vector symmetry
        gQ: coupling of mediator to quarks
        gDM: coupling of mediator to dark matter
        gL: coupling of mediator to charged leptons (neutrino coupling fixed by gauge invariance; only implemented for V/A models)
        mMed: mass of mediator in GeV
        mDM: mass of dark matter in GeV

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
